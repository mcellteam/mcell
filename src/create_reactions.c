/***********************************************************************************
 *                                                                                 *
 * Copyright (C) 2006-2013 by                                                      *
 * The Salk Institute for Biological Studies and                                   *
 * Pittsburgh Supercomputing Center, Carnegie Mellon University                    *
 *                                                                                 *
 * This program is free software; you can redistribute it and/or                   *
 * modify it under the terms of the GNU General Public License                     *
 * as published by the Free Software Foundation; either version 2                  *
 * of the License, or (at your option) any later version.                          *
 *                                                                                 *
 * This program is distributed in the hope that it will be useful,                 *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of                  *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                   *
 * GNU General Public License for more details.                                    *
 *                                                                                 *
 * You should have received a copy of the GNU General Public License               *
 * along with this program; if not, write to the Free Software                     *
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA. *
 *                                                                                 *
 ***********************************************************************************/

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "config.h"
#include "create_reactions.h"
#include "strfunc.h"



/* static functions */
static char* create_prod_signature(struct product **product_head);
static int sort_product_list_compare(struct product *list_item, 
    struct product *new_item);
static struct product* sort_product_list(struct product *product_head);


/*************************************************************************
 *
 * extract_reactants extracts the reactant info into a pathway structure
 *
 *************************************************************************/
MCELL_STATUS 
extract_reactants(struct pathway *path, struct species_opt_orient *reactants,
  int *num_reactants, int *num_vol_mols, int *num_grid_mols, int *all_3d) 
{
  for (struct species_opt_orient *current_reactant = reactants; 
       current_reactant != NULL; current_reactant = current_reactant->next)
  {
    short orient = current_reactant->orient_set ? current_reactant->orient : 0;
    struct species *reactant_species = (struct species *)current_reactant->mol_type->value;
  
    if ((reactant_species->flags & NOT_FREE) == 0)
    {
      (*num_vol_mols)++;
    } 
    else 
    {
      *all_3d = 0;
    }

    if (reactant_species->flags & ON_GRID)
    {
      (*num_grid_mols)++;
    }

    switch (*num_reactants)
    {
      case 0:
        path->reactant1 = reactant_species;
        path->orientation1 = orient;
        break;

      case 1:
        path->reactant2 = reactant_species;
        path->orientation2 = orient;
        break;

      case 2:
        path->reactant3 = reactant_species;
        path->orientation3 = orient;
        break;

      /* too many reactants */
      default: return MCELL_FAIL;
    }

    (*num_reactants)++;
  }

  return MCELL_SUCCESS;
}


/*************************************************************************
 *
 * extract_catalytic_arrow extracts the info for a catalytic arrow
 * into a pathway structure
 *
 *************************************************************************/
MCELL_STATUS 
extract_catalytic_arrow(struct pathway *path, 
    struct reaction_arrow *react_arrow, int *num_reactants, 
    int *num_vol_mols, int *num_grid_mols, int *all_3d) 
{
  if (*num_reactants >= 3) 
  {
    return MCELL_FAIL;
  }

  struct species *catalyst_species = (struct species *)react_arrow->catalyst.mol_type->value;
  short orient = react_arrow->catalyst.orient_set ? react_arrow->catalyst.orient : 0;

  /* XXX: Should surface class be allowed inside a catalytic arrow? */
  if (catalyst_species->flags & IS_SURFACE)
  {
    //mdlerror(parse_state, "A surface classes may not appear inside a catalytic arrow");
    return MCELL_FAIL;
  }

  /* Count the type of this reactant */
  if ((catalyst_species->flags & NOT_FREE) == 0)
  {
    (*num_vol_mols)++;
  }
  else 
  {
    *all_3d = 0;
  }

  if (catalyst_species->flags & ON_GRID) 
  {
    (*num_grid_mols)++;
  }

  /* Copy in catalytic reactant */
  switch (*num_reactants)
  {
    case 1:
      path->reactant2 = (struct species*)react_arrow->catalyst.mol_type->value;
      path->orientation2 = orient;
      break;

    case 2:
      path->reactant3 = (struct species*)react_arrow->catalyst.mol_type->value;
      path->orientation3 = orient;
      break;

    case 0:
    default:
      //mcell_internal_error("Catalytic reagent ended up in an invalid slot (%d).", reactant_idx);
      return MCELL_FAIL;
  }
  (*num_reactants)++;

  return MCELL_SUCCESS;
}



/*************************************************************************
 *
 * extract_surface extracts the info for a surface included in the
 * reaction specification
 *
 *************************************************************************/
MCELL_STATUS 
extract_surface(struct pathway *path, struct species_opt_orient *surf_class,
    int *num_reactants, int *num_surfaces, int *oriented_count) 
{
  short orient = surf_class->orient_set ? surf_class->orient : 0;
  if (surf_class->orient_set) {
    oriented_count++;
  }

  /* Copy reactant into next available slot */
  switch (*num_reactants)
  {
    case 0:
      //mdlerror(parse_state, "Before defining reaction surface class at least one reactant should be defined.");
      return MCELL_FAIL;

    case 1:
      path->reactant2 = (struct species*)surf_class->mol_type->value;
      path->orientation2 = orient;
      break;

    case 2:
      path->reactant3 = (struct species*)surf_class->mol_type->value;
      path->orientation3 = orient;
      break;

    default:
      //mdlerror(parse_state, "Too many reactants--maximum number is two plus reaction surface class.");
      return MCELL_FAIL;
  }

  num_reactants++;
  num_surfaces++;

  return MCELL_SUCCESS;
}



/*************************************************************************
 *
 * add_catalytic_species_to_products adds all species that are part of a
 * catalytic reaction to the list of products.
 *
 *************************************************************************/
MCELL_STATUS
add_catalytic_species_to_products(struct pathway *path, int catalytic,
    int bidirectional, int all_3d)
{
  struct species *catalyst;
  short catalyst_orient;
  switch (catalytic)
  {
    case 0: 
      catalyst = path->reactant1; 
      catalyst_orient = path->orientation1; 
      break;
    case 1: 
      catalyst = path->reactant2; 
      catalyst_orient = path->orientation2; 
      break;
    case 2: 
      catalyst = path->reactant3; 
      catalyst_orient = path->orientation3; 
      break;
    default:
      //mcell_internal_error("Catalytic reagent index is invalid.");
      return MCELL_FAIL;
  }

  if (bidirectional || !(catalyst->flags & IS_SURFACE))
  {
    struct product *prodp = (struct product*)CHECKED_MALLOC_STRUCT(
        struct product, "reaction product");
    if (prodp == NULL)
    {
      return MCELL_FAIL;
    }

    prodp->is_complex = 0;
    prodp->prod = catalyst;
    if (all_3d) 
    {
      prodp->orientation = 0;
    }
    else 
    {
      prodp->orientation = catalyst_orient;
    }
    prodp->next = path->product_head;
    path->product_head = prodp;
  }

  return MCELL_SUCCESS;
}



/*************************************************************************
 *
 * extract_products extracts the product info into a pathway structure
 *
 *************************************************************************/
MCELL_STATUS 
extract_products(struct pathway *path, struct species_opt_orient *products,
  int *num_surf_products, int *bidirectional, int *all_3d) 
{
  struct species_opt_orient *current_product;
  for (current_product = products;
       current_product != NULL;
       current_product = current_product->next)
  {
    /* Nothing to do for NO_SPECIES */
    if (current_product->mol_type == NULL)
    {
      continue;
    }

    /* Create new product */
    struct product *prodp = (struct product*)CHECKED_MALLOC_STRUCT(
       struct product, "reaction product"); 
    if (prodp == NULL)
    {
      return MCELL_FAIL;
    }

    /* Set product species and orientation */
    prodp->prod = (struct species *)current_product->mol_type->value;
    if (*all_3d) 
    {
      prodp->orientation = 0;
    }
    else 
    {
      prodp->orientation = current_product->orient;
    }

    /* Disallow surface as product unless reaction is bidirectional */
    if (!*bidirectional && (prodp->prod->flags & IS_SURFACE))
    {
      return MCELL_FAIL;
    }

    /* Append product to list */
    prodp->next = path->product_head;
    path->product_head = prodp;

    if (prodp->prod->flags & ON_GRID)
    {
      num_surf_products++;
    }

    /* Add product if it isn't a surface */
    if (!(prodp->prod->flags&IS_SURFACE))
    {
      if (all_3d == 0 && (!current_product->orient_set))
      {
        return MCELL_FAIL;  // product orientation not specified
      }
      else
      {
        if ((prodp->prod->flags&NOT_FREE)!=0)
        {
          return MCELL_FAIL; // trying to create surface product in presence 
                             // of only volume reactants
        }
        if (current_product->orient_set)
        {
            return MCELL_FAIL; // orientation specified for only volume reactants
        }
      }
    }
  }

  return MCELL_SUCCESS;
}


/*************************************************************************
 *
 * extract_pathname extracts the pathname (if one was given into 
 * a pathway structure
 *
 *************************************************************************/
MCELL_STATUS 
extract_pathname(struct pathway *path, struct rxn *rxnp, 
    struct sym_table *pathname)
{
  struct rxn_pathname *rxpnp = (struct rxn_pathname *)pathname->value;
  rxpnp->rx = rxnp;
  path->pathname = rxpnp;

  return MCELL_FAIL;
}



/*************************************************************************
 *
 * extract_rath extracts the forward rate of the reaction 
 *
 *************************************************************************/
MCELL_STATUS
extract_forward_rate(struct pathway *path, struct reaction_rates* rate,
    const char *rate_filename)
{
  switch (rate->forward_rate.rate_type)
  {
    case RATE_UNSET:
      return MCELL_FAIL;  // no rate set

    case RATE_CONSTANT:
      path->km = rate->forward_rate.v.rate_constant;
      path->km_filename = NULL;
      path->km_complex = NULL;
      break;

    case RATE_FILE:
      path->km = 0.0;
      path->km_filename = (char*)rate_filename;
      free(rate->forward_rate.v.rate_file);
      path->km_complex = NULL;
      break;

    case RATE_COMPLEX:
      path->km = 0.0;
      path->km_filename = NULL;
      path->km_complex = rate->forward_rate.v.rate_complex;
      break;

    default: 
      //UNHANDLED_CASE(rate->forward_rate.rate_type);
      return MCELL_FAIL;
  }

  return MCELL_SUCCESS;
}



/*************************************************************************
 *
 * create_product_signature for the pathway 
 *
 *************************************************************************/
MCELL_STATUS
create_product_signature(struct pathway *path) 
{
  if (path->product_head != NULL)
  {
    path->prod_signature = create_prod_signature(&path->product_head);
    if (path->prod_signature == NULL)
    {
      return MCELL_FAIL;  // creation of field failed
    }
  }
  else
  {
    path->prod_signature = NULL;
  }

  return MCELL_SUCCESS;
}



/*************************************************************************
 *
 * grid_space_available_for_surface_products checks for enough available
 * grid space for surface products. 
 * If the vacancy search distance is zero and this reaction produces more
 * grid molecules than it comsumes, it can never succeed, except if it is a
 * volume molecule hitting the surface and producing a single grid molecule.
 * Fail with an error message.
 *
 *************************************************************************/
MCELL_STATUS
grid_space_available_for_surface_products(double vacancy_search_dist2,
    int num_grid_mols, int num_vol_mols, int num_surf_products) 
{
  if ((vacancy_search_dist2 == 0) && (num_surf_products > num_grid_mols))
  {
    /* The case with one volume molecule reacting with the surface and
     * producing one grid molecule is okay.
     */
    if (num_grid_mols == 0 && num_vol_mols == 1 && num_surf_products == 1)
    {
      /* do nothing */
    }
    else
    {
      return MCELL_FAIL;  // number of surface products exceeds number of 
                          // surface reactants but VACANCY_SEARCH_DISTANCE
                          // is not specified
    }
  }
  return MCELL_SUCCESS;
}





/*************************************************************************
 create_rx_name:
    Assemble reactants alphabetically into a reaction name string.

 In:  p: reaction pathway whose reaction name we are to create
 Out: a string to be used as a symbol name for the reaction
*************************************************************************/
char*
create_rx_name( struct pathway *p)
{

  struct species *reagents[3];
  int n_reagents = 0;
  int is_complex = 0;

  /* Store reagents in an array. */
  reagents[0] = p->reactant1;
  reagents[1] = p->reactant2;
  reagents[2] = p->reactant3;

  /* Count non-null reagents. */
  for (n_reagents = 0; n_reagents < 3; ++ n_reagents)
    if (reagents[n_reagents] == NULL)
      break;
    else if (p->is_complex[n_reagents])
      is_complex = 1;

  /* Sort reagents. */
  for (int i = 0; i<n_reagents; ++i)
  {
    for (int j = i+1; j<n_reagents; ++ j)
    {
      /* If 'i' is a subunit, 'i' wins. */
      if (p->is_complex[i])
        break;

      /* If 'j' is a subunit, 'j' wins. */
      else if (p->is_complex[j])
      {
        struct species *tmp = reagents[j];
        reagents[j] = reagents[i];
        reagents[i] = tmp;
      }

      /* If 'j' precedes 'i', 'j' wins. */
      else if (strcmp(reagents[j]->sym->name, reagents[i]->sym->name) < 0)
      {
        struct species *tmp = reagents[j];
        reagents[j] = reagents[i];
        reagents[i] = tmp;
      }
    }
  }

  /* Now, produce a name! */
  if (is_complex)
  {
    switch (n_reagents)
    {
      case 1: return alloc_sprintf("(%s)", reagents[0]->sym->name);
      case 2: return alloc_sprintf("(%s)+%s", reagents[0]->sym->name, reagents[1]->sym->name);
      case 3: return alloc_sprintf("(%s)+%s+%s", reagents[0]->sym->name, reagents[1]->sym->name, reagents[2]->sym->name);
      default:
        //mcell_internal_error("Invalid number of reagents in reaction pathway (%d).", n_reagents);
        return NULL;
    }
  }
  else
  {
    switch (n_reagents)
    {
      case 1: return alloc_sprintf("%s", reagents[0]->sym->name);
      case 2: return alloc_sprintf("%s+%s", reagents[0]->sym->name, reagents[1]->sym->name);
      case 3: return alloc_sprintf("%s+%s+%s", reagents[0]->sym->name, reagents[1]->sym->name, reagents[2]->sym->name);
      default:
        //mcell_internal_error("Invalid number of reagents in reaction pathway (%d).", n_reagents);
        return NULL;
    }
  }
}



/************************************************************************
 * static helper functions
 ************************************************************************/

/*************************************************************************
 sort_product_list_compare:
    Comparison function for products to be sorted when generating the product
    signature.

 In:  list_item: first item to compare
      new_item:  second item to compare
 Out: -1 if list_item < new_item, 1 if list_item > new_item, 0 if they are
      equal

  XXX Currently this function also appears in mdlparse_util.c. It should
      eventually be removed from there and only appear in this file.
*************************************************************************/
static int 
sort_product_list_compare(struct product *list_item, struct product *new_item)
{
  int cmp = list_item->is_complex - new_item->is_complex;
  if (cmp != 0)
    return cmp;

  cmp = strcmp(list_item->prod->sym->name, new_item->prod->sym->name);
  if (cmp == 0)
  {
    if (list_item->orientation > new_item->orientation)
      cmp = -1;
    else if (list_item->orientation < new_item->orientation)
      cmp = 1;
    else
      cmp = 0;
  }
  return cmp;
}

/*************************************************************************
 sort_product_list:
    Sorts product_head in alphabetical order, and descending orientation order.
    Current algorithm uses insertion sort.

 In:  product_head: list to sort
 Out: the new list head

  XXX Currently this function also appears in mdlparse_util.c. It should
      eventually be removed from there and only appear in this file.
*************************************************************************/
static struct product*
sort_product_list(struct product *product_head)
{
  struct product *next;             /* Saved next item (next field in product is overwritten) */
  struct product *iter;             /* List iterator */
  struct product *result = NULL;    /* Sorted list */
  int cmp;

  /* Use insertion sort to sort the list of products */
  for (struct product *current = product_head;
       current != NULL;
       current = next)
  {
    next = current->next;

    /* First item added always goes at the head */
    if (result == NULL)
    {
      current->next = result;
      result = current;
      continue;
    }

    /* Check if the item belongs at the head */
    cmp = sort_product_list_compare(result, current);
    if (cmp >= 0)
    {
      current->next = result;
      result = current;
      continue;
    }

    /* Otherwise, if it goes after the current entry, scan forward to find the insert point */
    else
    {
      /* locate the node before the point of insertion */
      iter = result;
      while (iter->next != NULL  &&  sort_product_list_compare(iter, current) < 0)
        iter = iter->next;

      current->next = iter->next;
      iter->next = current;
    }
  }

  return result;
}


/*************************************************************************
 create_prod_signature:
    Returns a string containing all products in the product_head list,
    separated by '+', and sorted in alphabetical order by name and descending
    orientation order.

 In:  product_head: list of products
 Out: product signature as a string.  *product_head list is sorted in
      alphabetical order by name, and descending order by orientation.  Returns
      NULL on failure.

  XXX Currently this function also appears in mdlparse_util.c. It should
      eventually be removed from there and only appear in this file.
*************************************************************************/
static char*
create_prod_signature( struct product **product_head)
{
  /* points to the head of the sorted alphabetically list of products */
  char *prod_signature = NULL;

  *product_head = sort_product_list(*product_head);

  /* create prod_signature string */
  struct product *current = *product_head;
  prod_signature = CHECKED_STRDUP(current->prod->sym->name, "product name");

  /* Concatenate to create product signature */
  char *temp_str = NULL;
  while (current->next != NULL)
  {
    temp_str = prod_signature;
    prod_signature = CHECKED_SPRINTF("%s+%s",
                                     prod_signature,
                                     current->next->prod->sym->name);

    if (prod_signature == NULL)
    {
      if (temp_str != NULL) free(temp_str);
      return NULL;
    }
    if (temp_str != NULL) free(temp_str);

    current = current->next;
  }

  return prod_signature;
}


