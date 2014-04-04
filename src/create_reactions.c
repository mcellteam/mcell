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

#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "config.h"
#include "create_reactions.h"
#include "logging.h"
#include "macromolecule.h"
#include "react_util.h"
#include "strfunc.h"
#include "sym_table.h"



/* static functions */
//static char* create_prod_signature(struct product **product_head);
static int sort_product_list_compare(struct product *list_item, 
    struct product *new_item);
static struct product* sort_product_list(struct product *product_head);


/*************************************************************************
 *
 * extract_reactants extracts the reactant info into a pathway structure
 *
 *************************************************************************/
MCELL_STATUS 
extract_reactants(struct pathway *pathp, struct species_opt_orient *reactants,
  int *num_reactants, int *num_vol_mols, int *num_grid_mols, 
  int *num_complex_reactants, int *all_3d, int *oriented_count, int *complex_type) 
{
  int reactant_idx = 0;
  struct species_opt_orient *current_reactant;
  for (current_reactant = reactants;
       reactant_idx < 3 && current_reactant != NULL;
       ++ reactant_idx, current_reactant = current_reactant->next)
  {
    /* Extract orientation and species */
    short orient = current_reactant->orient_set ? current_reactant->orient : 0;
    struct species *reactant_species = (struct species *)current_reactant->mol_type->value;

    /* Count the type of this reactant */
    if (current_reactant->orient_set)
    {
      ++(*oriented_count);
    }
    
    if (reactant_species->flags & NOT_FREE)
    {
      *all_3d = 0;
      if (reactant_species->flags & ON_GRID)
      {
        ++(*num_grid_mols);
      }
    }
    else
    {
      ++(*num_vol_mols);
    }

    /* Sanity check this reactant */
    if (current_reactant->is_subunit)
    {
      if ((reactant_species->flags & NOT_FREE) == 0)
      {
        *complex_type = TYPE_3D;
      }
      else if (reactant_species->flags & ON_GRID)
      {
        *complex_type = TYPE_GRID;
      }
      else
      {
        //mdlerror(parse_state, "Only a molecule may be used as a macromolecule subunit in a reaction.");
        return MCELL_FAIL;
      }
      ++(*num_complex_reactants);
    }
    else if (reactant_species->flags & IS_SURFACE)
    {
      //mdlerror(parse_state, "Surface class can be listed only as the last reactant on the left-hand side of the reaction with the preceding '@' sign.");
      return MCELL_FAIL;
    }

    /* Copy in reactant info */
    pathp->is_complex[reactant_idx] = current_reactant->is_subunit;
    switch (reactant_idx)
    {
      case 0:
        pathp->reactant1 = reactant_species;
        pathp->orientation1 = orient;
        break;

      case 1:
        pathp->reactant2 = reactant_species;
        pathp->orientation2 = orient;
        break;

      case 2:
        pathp->reactant3 = reactant_species;
        pathp->orientation3 = orient;
        break;

      default: UNHANDLED_CASE(reactant_idx);
    }
  }
  *num_reactants = reactant_idx;

  /* we had more than 3 reactants */
  if (current_reactant != NULL) 
  { 
    return MCELL_FAIL;
  }
  else
  {
    return MCELL_SUCCESS;
  }
}


/*************************************************************************
 *
 * extract_catalytic_arrow extracts the info for a catalytic arrow
 * into a pathway structure
 *
 *************************************************************************/
MCELL_STATUS 
extract_catalytic_arrow(struct pathway *pathp, 
    struct reaction_arrow *react_arrow, int *reactant_idx, 
    int *num_vol_mols, int *num_grid_mols, int *all_3d,
    int *oriented_count) 
{
  struct species *catalyst_species = (struct species *) react_arrow->catalyst.mol_type->value;
  short orient = react_arrow->catalyst.orient_set ? react_arrow->catalyst.orient : 0;

  /* XXX: Should surface class be allowed inside a catalytic arrow? */
  if (catalyst_species->flags & IS_SURFACE)
  {
    //mdlerror(parse_state, "A surface classes may not appear inside a catalytic arrow");
    return MCELL_FAIL;
  }

  /* Count the type of this reactant */
  if (react_arrow->catalyst.orient_set)
  {
    ++(*oriented_count);
  }

  if (catalyst_species->flags & NOT_FREE)
  {
    *all_3d = 0;
    if (catalyst_species->flags & ON_GRID)
    {
      ++(*num_grid_mols);
    }
  }
  else
  {
    ++(*num_vol_mols);
  }

    /* Copy in catalytic reactant */
  switch (*reactant_idx)
  {
    case 1:
      pathp->reactant2 = (struct species*) react_arrow->catalyst.mol_type->value;
      pathp->orientation2 = orient;
      break;

    case 2:
      pathp->reactant3 = (struct species*) react_arrow->catalyst.mol_type->value;
      pathp->orientation3 = orient;
      break;

    case 0:
    default:
      //mcell_internal_error("Catalytic reagent ended up in an invalid slot (%d).", reactant_idx);
      return MCELL_FAIL;
  }
  ++(*reactant_idx);

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
    int *num_reactants, unsigned int *num_surfaces, int *oriented_count)
{
  short orient = surf_class->orient_set ? surf_class->orient : 0;
  if (surf_class->orient_set) {
    (*oriented_count)++;
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

  (*num_reactants)++;
  (*num_surfaces)++;

  return MCELL_SUCCESS;
}


/*************************************************************************
 *
 * check_surface_specs performs a number of sanity checks to make sure
 * the surface specifications are sane
 *
 *************************************************************************/
MCELL_STATUS 
check_surface_specs(MCELL_STATE *state, int num_reactants, int num_surfaces, 
    int num_vol_mols, int all_3d, int oriented_count)
{
  if (num_surfaces > 1)
  {
    /* Shouldn't happen */
    mcell_internal_error("Too many surfaces--reactions can take place on at most one surface.");
    return MCELL_FAIL;
  }

  if (num_surfaces == num_reactants)
  {
    mcell_error("Reactants cannot consist entirely of surfaces.  Use a surface release site instead!");
    return MCELL_FAIL;
  }

  if ((num_vol_mols == 2) && (num_surfaces == 1))
  {
    mcell_error("Reaction between two volume molecules and a surface is not defined.");
    return MCELL_FAIL;
  }

  if (all_3d)
  {
    if (oriented_count != 0)
    {
      if (state->notify->useless_vol_orient==WARN_ERROR)
      {
        mcell_error("Error: orientation specified for molecule in reaction in volume");
        return MCELL_FAIL;
      }
      else if (state->notify->useless_vol_orient==WARN_WARN)
      {
        mcell_error("Warning: orientation specified for molecule in reaction in volume");
      }
    }
  }
  else
  {
    if (num_reactants != oriented_count)
    {
      if (state->notify->missed_surf_orient==WARN_ERROR)
      {
        mcell_error("Error: orientation not specified for molecule in reaction at surface\n  (use ; or ', or ,' for random orientation)");
        return MCELL_FAIL;
      }
      else if (state->notify->missed_surf_orient==WARN_WARN)
      {
        mcell_error("Warning: orientation not specified for molecule in reaction at surface\n  (use ; or ', or ,' for random orientation)");
      }
    }
  }

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
      mcell_internal_error("Catalytic reagent index is invalid.");
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
extract_products(MCELL_STATE *state, struct pathway *pathp, 
  struct species_opt_orient *products, int *num_surf_products, 
  int *num_complex_products, int bidirectional, int complex_type, 
  int all_3d) 
{
  struct species_opt_orient *current_product;
  for (current_product = products;
       current_product != NULL;
       current_product = current_product->next)
  {
    /* Nothing to do for NO_SPECIES */
    if (current_product->mol_type == NULL)
      continue;

    /* Create new product */
    struct product *prodp = (struct product*)CHECKED_MALLOC_STRUCT(struct product, 
        "reaction product");

    if (prodp == NULL)
    {
      //mcell_error_raw("Out of memory while creating reaction: %s -> ... ",
      //                rxnp->sym->name);
      return MCELL_FAIL;
    }

    /* Set product species and orientation */
    prodp->prod = (struct species *) current_product->mol_type->value;
    if (all_3d) 
    {
      prodp->orientation = 0;
    }
    else 
    {
      prodp->orientation = current_product->orient;
    }

    /* Disallow surface as product unless reaction is bidirectional */
    if (!bidirectional)
    {
      if (prodp->prod->flags & IS_SURFACE)
      {
        mcell_error_raw("Surface_class '%s' is not allowed to be on the product side of the reaction.", prodp->prod->sym->name);
        return MCELL_FAIL;
      }
    }

    /* Copy over complex-related state for product */
    prodp->is_complex = current_product->is_subunit;
    if (current_product->is_subunit)
    {
      ++(*num_complex_products);
      if ((prodp->prod->flags & NOT_FREE) != 0)
      {
        if (complex_type == TYPE_3D)
        {
          mcell_error_raw("Volume subunit cannot become a surface subunit '%s' in a macromolecular reaction.", prodp->prod->sym->name);
          return MCELL_FAIL;
        }
      }
      else if ((prodp->prod->flags & ON_GRID) == 0)
      {
        if (complex_type == TYPE_GRID)
        {
          mcell_error_raw("Surface subunit cannot become a volume subunit '%s' in a macromolecular reaction.", prodp->prod->sym->name);
          return MCELL_FAIL;
        }
      }
      else
      {
        mcell_error_raw("Only a molecule may be used as a macromolecule subunit in a reaction.");
        return MCELL_FAIL;
      }
    }

    /* Append product to list */
    prodp->next = pathp->product_head;
    pathp->product_head = prodp;

    if (prodp->prod->flags & ON_GRID)
    {
      ++(*num_surf_products);
    }

    /* Add product if it isn't a surface */
    if (! (prodp->prod->flags&IS_SURFACE))
    {
      if (all_3d == 0)
      {
        if (!current_product->orient_set)
        {
          if (state->notify->missed_surf_orient==WARN_ERROR)
          {
            mcell_error_raw("Error: product orientation not specified in reaction with orientation\n  (use ; or ', or ,' for random orientation)");
            return MCELL_FAIL;
          }
          else if (state->notify->missed_surf_orient==WARN_WARN)
          {
            mcell_error_raw("Warning: product orientation not specified for molecule in reaction at surface\n  (use ; or ', or ,' for random orientation)");
          }
        }
      }
      else
      {
        if ((prodp->prod->flags&NOT_FREE)!=0)
        {
          mcell_error("Reaction has only volume reactants but is trying to create a surface product");
          return MCELL_FAIL;
        }
        if (current_product->orient_set)
        {
          if (state->notify->useless_vol_orient==WARN_ERROR)
          {
            mcell_error("Error: orientation specified for molecule in reaction in volume");
            return MCELL_FAIL;
          }
          else if (state->notify->useless_vol_orient==WARN_WARN)
          {
            mcell_error("Warning: orientation specified for molecule in reaction at surface");
          }
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
create_rx_name(struct pathway *p)
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


/*************************************************************************
 concat_rx_name:
    Concatenates reactants onto a reaction name.  Reactants which are subunits
    in macromolecular complexes will have their names parenthesized.

 In:  parse_state: parser state
      name1: name of first reactant (or first part of reaction name)
      is_complex1: 0 unless the first reactant is a subunit in a complex
      name2: name of second reactant (or second part of reaction name)
      is_complex2: 0 unless the second reactant is a subunit in a complex
 Out: reaction name as a string, or NULL if an error occurred
*************************************************************************/
static char*
concat_rx_name(char *name1, int is_complex1, char *name2, int is_complex2)
{
  char *rx_name;

  /* Make sure they aren't both subunits  */
  if (is_complex1  &&  is_complex2)
  {
    //mdlerror_fmt(parse_state, "File '%s', Line %ld: Internal error -- a reaction cannot have two reactants which are subunits of a macromolecule.", __FILE__, (long)__LINE__);
    return NULL;
  }

  /* Sort them */
  if (is_complex2  ||  strcmp(name2, name1) <= 0)
  {
    char *nametmp = name1;
    int is_complextmp = is_complex1;
    name1 = name2;
    is_complex1 = is_complex2;
    name2 = nametmp;
    is_complex2 = is_complextmp;
    assert(is_complex2 == 0);
  }

  /* Build the name */
  if (is_complex1)
    rx_name = CHECKED_SPRINTF("(%s)+%s", name1, name2);
  else
    rx_name = CHECKED_SPRINTF("%s+%s", name1, name2);

  /* Die if we failed to allocate memory */
  if (rx_name == NULL)
    return NULL;

  return rx_name;
}


/***********************************************************************
 invert_current_reaction_pathway:
    Creates a new reversed pathway, where the reactants of new pathway are the
    products of the current pathway, and the products of new pathway are the
    reactants of the current pathway.

 In:  parse_state: parser state
      pathp: pathway to invert
      reverse_rate: the reverse reaction rate
 Out: Returns 1 on error and 0 - on success.  The new pathway is added to the
      linked list of the pathways for the current reaction.
***********************************************************************/
MCELL_STATUS invert_current_reaction_pathway(MCELL_STATE *state,
    struct pathway *pathp, struct reaction_rate *reverse_rate,
    const char *rate_filename)
{
  struct rxn *rx;
  struct pathway *path;
  struct product *prodp;
  struct sym_table *sym;
  char *inverse_name;
  int nprods;  /* number of products */
  int all_3d;  /* flag that tells whether all products are volume_molecules */
  int num_surf_products = 0;
  int num_grid_mols = 0;
  int num_vol_mols = 0;


  /* flag that tells whether there is a surface_class
     among products in the direct reaction */
  int is_surf_class = 0;

  all_3d=1;
  for (nprods=0,prodp=pathp->product_head ; prodp!=NULL ; prodp=prodp->next)
  {
    nprods++;
    if ((prodp->prod->flags&NOT_FREE)!=0) all_3d=0;
    if ((prodp->prod->flags & IS_SURFACE) != 0)
    {
           is_surf_class = 1;
    }
  }

  if (nprods==0)
  {
    //mdlerror(parse_state, "Can't create a reverse reaction with no products");
    return MCELL_FAIL;
  }
  if (nprods==1 && (pathp->product_head->prod->flags&IS_SURFACE))
  {
    //mdlerror(parse_state, "Can't create a reverse reaction starting from only a surface");
    return MCELL_FAIL;
  }
  if (nprods>3)
  {
    //mdlerror(parse_state, "Can't create a reverse reaction involving more than three products. Please note that surface_class from the reaction reactant side also counts as a product.");
    return MCELL_FAIL;
  }

  if (pathp->pathname != NULL)
  {
    //mdlerror(parse_state, "Can't name bidirectional reactions--write each reaction and name them separately");
    return MCELL_FAIL;
  }
  if (all_3d)
  {
    if ((pathp->reactant1->flags&NOT_FREE)!=0) all_3d = 0;
    if (pathp->reactant2!=NULL && (pathp->reactant2->flags&NOT_FREE)!=0) all_3d = 0;
    if (pathp->reactant3!=NULL && (pathp->reactant3->flags&NOT_FREE)!=0) all_3d = 0;

    if (!all_3d)
    {
      //mdlerror(parse_state, "Cannot reverse orientable reaction with only volume products");
      return MCELL_FAIL;
    }
  }

  prodp = pathp->product_head;
  if (nprods==1)
  {
    if (prodp->is_complex)
    {
      inverse_name = CHECKED_SPRINTF("(%s)",
                                     prodp->prod->sym->name);
    }
    else
      inverse_name = strdup(prodp->prod->sym->name);

    if (inverse_name == NULL)
      return MCELL_FAIL;
  }
  else if (nprods == 2)
  {
    inverse_name = concat_rx_name(prodp->prod->sym->name, prodp->is_complex, prodp->next->prod->sym->name, prodp->next->is_complex);
  }
  else
  {
    if (prodp->is_complex || prodp->next->is_complex || prodp->next->next->is_complex)
    {
      //mdlerror(parse_state, "MCell does not currently support trimolecular reactions for macromolecules");
      return MCELL_FAIL;
    }
    inverse_name = concat_rx_name(prodp->prod->sym->name, 0, prodp->next->prod->sym->name, 0);
    inverse_name = concat_rx_name(inverse_name, 0, prodp->next->next->prod->sym->name, 0);
  }
  if (inverse_name==NULL)
  {
    //mdlerror(parse_state, "Out of memory forming reaction name");
    return MCELL_FAIL;
  }

  sym = retrieve_sym(inverse_name, state->rxn_sym_table);
  if (sym==NULL)
  {
    sym = store_sym(inverse_name,RX,state->rxn_sym_table, NULL);
    if (sym==NULL)
    {
      //mdlerror_fmt(parse_state, "File '%s', Line %ld: Out of memory while storing reaction pathway.", __FILE__, (long)__LINE__);
      return MCELL_FAIL;
    }
  }
  free(inverse_name);
  rx = (struct rxn*)sym->value;
  rx->n_reactants = nprods;
  rx->n_pathways++;

  path = (struct pathway*)CHECKED_MALLOC_STRUCT(struct pathway, "reaction pathway");
  if (path == NULL)
  {
    return MCELL_FAIL;
  }
  path->pathname=NULL;
  path->flags = 0;
  path->reactant1=prodp->prod;
  if ((path->reactant1->flags & NOT_FREE) == 0)
  {
         ++ num_vol_mols;
  }
  else 
  {
     if (path->reactant1->flags & ON_GRID)
     {
         ++ num_grid_mols;
     }
  }
  path->is_complex[0] = prodp->is_complex;
  path->is_complex[1] = 0;
  path->is_complex[2] = 0;
  path->orientation1 = prodp->orientation;
  path->reactant2=NULL;
  path->reactant3=NULL;
  path->prod_signature = NULL;
  if (nprods > 1)
  {
      path->reactant2 = prodp->next->prod;
      if ((path->reactant2->flags & NOT_FREE) == 0)
      {
         ++ num_vol_mols;
      }
      else 
      {
         if (path->reactant2->flags & ON_GRID)
         {
           ++ num_grid_mols;
         }
      }
      path->orientation2 = prodp->next->orientation;
      path->is_complex[1] = prodp->next->is_complex;
  }
  if (nprods > 2)
  {
      path->reactant3 = prodp->next->next->prod;
      if ((path->reactant3->flags & NOT_FREE) == 0)
      {
         ++ num_vol_mols;
      }
      else 
      {
         if (path->reactant3->flags & ON_GRID)
         {
           ++ num_grid_mols;
         }
      }
      path->orientation3 = prodp->next->next->orientation;
  }

  switch (reverse_rate->rate_type)
  {
    case RATE_UNSET:
      //mdlerror_fmt(parse_state, "File %s, Line %d: Internal error: Reverse rate is not set", __FILE__, __LINE__);
      return MCELL_FAIL;

    case RATE_CONSTANT:
      path->km = reverse_rate->v.rate_constant;
      path->km_filename = NULL;
      path->km_complex = NULL;
      break;

    case RATE_FILE:
      path->km = 0.0;
      path->km_filename = (char*)rate_filename;
      free(reverse_rate->v.rate_file);
      path->km_complex = NULL;
      break;

    case RATE_COMPLEX:
      path->km = 0.0;
      path->km_filename = NULL;
      path->km_complex = reverse_rate->v.rate_complex;
      break;

    default: UNHANDLED_CASE(reverse_rate->rate_type);
  }

  path->product_head = (struct product*)CHECKED_MALLOC_STRUCT(struct product, "reaction product");
  if (path->product_head == NULL)
    return 1;

  path->product_head->orientation = pathp->orientation1;
  path->product_head->prod = pathp->reactant1;
  path->product_head->is_complex = pathp->is_complex[0];
  path->product_head->next = NULL;
  if (path->product_head->prod->flags & ON_GRID)
    ++ num_surf_products;

  if ((pathp->reactant2!=NULL) && ((pathp->reactant2->flags & IS_SURFACE) == 0))
  {
    path->product_head->next = (struct product*)CHECKED_MALLOC_STRUCT(struct product, "reaction product");
    if (path->product_head->next == NULL)
      return 1;
    path->product_head->next->orientation = pathp->orientation2;
    path->product_head->next->prod = pathp->reactant2;
    path->product_head->next->is_complex = pathp->is_complex[1];
    path->product_head->next->next = NULL;
    if (path->product_head->next->prod->flags & ON_GRID)
      ++ num_surf_products;
  }

  if ((pathp->reactant3!=NULL) && ((pathp->reactant3->flags & IS_SURFACE) == 0))
  {
    path->product_head->next->next = (struct product*)CHECKED_MALLOC_STRUCT(struct product, "reaction product");
    if (path->product_head->next->next == NULL)
      return 1;
    path->product_head->next->next->orientation = pathp->orientation3;
    path->product_head->next->next->prod = pathp->reactant3;
    path->product_head->next->next->is_complex = pathp->is_complex[2];
    path->product_head->next->next->next = NULL;
    if (path->product_head->next->next->prod->flags & ON_GRID)
      ++ num_surf_products;
  }

  path->prod_signature = create_prod_signature(&path->product_head);
  if (path->prod_signature == NULL)
  {
      //mdlerror(parse_state, "Error creating 'prod_signature' field for reaction pathway.");
      return MCELL_FAIL;
  }


  if ((state->vacancy_search_dist2 == 0) &&
     (num_surf_products > num_grid_mols))
  {
      /* the case with one volume molecule reacting with the surface
         and producing one grid molecule is excluded */
      if (!((num_grid_mols == 0) && (num_vol_mols == 1)))
      {
          //mdlerror(parse_state, "Error: number of surface products exceeds number of surface reactants, but VACANCY_SEARCH_DISTANCE is not specified or set to zero.");
           return MCELL_FAIL;
      }
  }

  /* Now go back to the original reaction and if there is a "surface_class"
     among products - remove it.  We do not need it now on the product side
     of the reaction */
   if (is_surf_class)
   {
     prodp = pathp->product_head;
     if (prodp->prod->flags & IS_SURFACE)
     {
       pathp->product_head = prodp->next;
       prodp->next = NULL;
       //mem_put(parse_state->prod_mem, (void *)prodp);
     }
     else if (prodp->next->prod->flags & IS_SURFACE)
     {
       //struct product *temp = prodp->next;
       prodp->next = prodp->next->next;
       //mem_put(parse_state->prod_mem, temp);
     }
     else
     {
       //struct product *temp = prodp->next->next;
       prodp->next->next = prodp->next->next->next;
       //mem_put(parse_state->prod_mem, temp);
     }
   }

  path->next = rx->pathway_head;
  rx->pathway_head = path;
  return 0;
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
char*
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


/*************************************************************************
 * init_reactions and related machinery
 *************************************************************************/

/*************************************************************************
 check_duplicate_special_reactions:
   Check for duplicate special reaction pathways (e.g. TRANSPARENT = molecule).

 In: path: Parse-time structure for reaction pathways
 Out: Nothing. 
 Note: I'm not sure if this code is ever actually called.
*************************************************************************/
static void 
check_duplicate_special_reactions(struct pathway *path)
{
  /* if it is a special reaction - check for the duplicates pathways */
  if (path->next != NULL)
  {
    if ((path->flags & PATHW_TRANSP) && (path->next->flags & PATHW_TRANSP))
    {
      if ((path->orientation2 == path->next->orientation2) ||
         (path->orientation2 == 0) || (path->next->orientation2 == 0))
      {
         mcell_error("Exact duplicates of special reaction TRANSPARENT = %s are not allowed.  Please verify the contents of DEFINE_SURFACE_CLASS statement.", path->reactant2->sym->name);
      }
    }

    if ((path->flags & PATHW_REFLEC) && (path->next->flags & PATHW_REFLEC))
    {
      if ((path->orientation2 == path->next->orientation2) ||
         (path->orientation2 == 0) || (path->next->orientation2 == 0))
      {
         mcell_error("Exact duplicates of special reaction REFLECTIVE = %s are not allowed.  Please verify the contents of DEFINE_SURFACE_CLASS statement.", path->reactant2->sym->name);
      }
    }
    if ((path->flags & PATHW_ABSORP) && (path->next->flags & PATHW_ABSORP))
    {
      if ((path->orientation2 == path->next->orientation2) ||
         (path->orientation2 == 0) || (path->next->orientation2 == 0))
      {
        mcell_error("Exact duplicates of special reaction ABSORPTIVE = %s are not allowed.  Please verify the contents of DEFINE_SURFACE_CLASS statement.", path->reactant2->sym->name);
      }
    }
  }
}


/*************************************************************************
 set_product_geometries:
 
  Walk through the list, setting the geometries of each of the products. We do
  this by looking for an earlier geometric match and pointing there or we just
  point to 0 if there is no match.

 In: path: Parse-time structure for reaction pathways
     rx: Pathways leading away from a given intermediate
     prod: Parse-time structure for products of reaction pathways
 Out: max_num_surf_products: Maximum number of surface products 
*************************************************************************/
static int 
set_product_geometries(struct pathway *path, struct rxn *rx, struct product *prod)
{
  int recycled1, recycled2, recycled3;
  int k, kk, k2;
  short geom;
  struct product *prod2;
  int max_num_surf_products;         /* maximum number of surface products */
  int num_surf_products_per_pathway;

  max_num_surf_products = 0;
  for (int n_pathway=0; path!=NULL ; n_pathway++ , path = path->next)
  {
    recycled1 = 0;
    recycled2 = 0;
    recycled3 = 0;
    k = rx->product_idx[n_pathway] + rx->n_reactants;
    num_surf_products_per_pathway = 0;
    for (prod=path->product_head ; prod != NULL ; prod = prod->next)
    {
      if (recycled1==0 && prod->prod == path->reactant1)
      {
        recycled1 = 1;
        kk = rx->product_idx[n_pathway] + 0;
      }
      else if (recycled2==0 && prod->prod == path->reactant2)
      {
        recycled2 = 1;
        kk = rx->product_idx[n_pathway] + 1;
      }
      else if (recycled3==0 && prod->prod == path->reactant3)
      {
        recycled3 = 1;
        kk = rx->product_idx[n_pathway] + 2;
      }
      else
      {
        kk = k;
        k++;
      }

      if (prod->prod->flags & ON_GRID) num_surf_products_per_pathway++;

      rx->players[kk] = prod->prod;
      if (rx->is_complex) rx->is_complex[kk] = prod->is_complex;

      if ((prod->orientation+path->orientation1)*(prod->orientation-path->orientation1)==0 && prod->orientation*path->orientation1!=0)
      {
        if (prod->orientation == path->orientation1) rx->geometries[kk] = 1;
        else rx->geometries[kk] = -1;
      }
      else if (rx->n_reactants > 1 &&
                (prod->orientation+path->orientation2)*(prod->orientation-path->orientation2)==0 && prod->orientation*path->orientation2!=0
             )
      {
        if (prod->orientation == path->orientation2) rx->geometries[kk] = 2;
        else rx->geometries[kk] = -2;
      }
      else if (rx->n_reactants > 2 &&
                (prod->orientation+path->orientation3)*(prod->orientation-path->orientation3)==0 && prod->orientation*path->orientation3!=0
             )
      {
        if (prod->orientation == path->orientation3) rx->geometries[kk] = 3;
        else rx->geometries[kk] = -3;
      }
      else
      {
        k2 = 2*rx->n_reactants + 1;  /* Geometry index of first non-reactant product, counting from 1. */
        geom = 0;
        for (prod2=path->product_head ; prod2!=prod && prod2!=NULL && geom==0 ; prod2 = prod2->next)
        {
          if ((prod2->orientation+prod->orientation)*(prod2->orientation-prod->orientation)==0 && prod->orientation*prod2->orientation!=0)
          {
            if (prod2->orientation == prod->orientation) geom = 1;
            else geom = -1;
          }
          else geom = 0;

          if (recycled1 == 1)
          {
            if (prod2->prod == path->reactant1)
            {
              recycled1 = 2;
              geom *= rx->n_reactants+1;
            }
          }
          else if (recycled2==1)
          {
            if (prod2->prod == path->reactant2)
            {
              recycled2 = 2;
              geom *= rx->n_reactants+2;
            }
          }
          else if (recycled3==1)
          {
            if (prod2->prod == path->reactant3)
            {
              recycled3 = 2;
              geom *= rx->n_reactants+3;
            }
          }
          else
          {
            geom *= k2;
            k2++;
          }
        }
        rx->geometries[kk] = geom;
      }
      if (num_surf_products_per_pathway > max_num_surf_products) max_num_surf_products = num_surf_products_per_pathway;
    }

    k = rx->product_idx[n_pathway];
    if (recycled1==0) rx->players[k] = NULL;
    if (recycled2==0 && rx->n_reactants>1) rx->players[k+1] = NULL;
    if (recycled3==0 && rx->n_reactants>2) rx->players[k+2] = NULL;
  } /* end for (n_pathway = 0, ...) */
  return max_num_surf_products;
}



/*************************************************************************
 alphabetize_pathway:
    The reaction pathway (path) is alphabetized.

 In: path: Parse-time structure for reaction pathways
     reaction: Reaction pathways leading away from a given intermediate
 Out: Nothing. 
*************************************************************************/
static void 
alphabetize_pathway(struct pathway *path, struct rxn *reaction)
{
  unsigned char temp_is_complex;
  short geom, geom2;
  struct species *temp_sp, *temp_sp2;

  /* Alphabetize if we have two molecules */
  if ((path->reactant2->flags&IS_SURFACE)==0)
  {
    if (strcmp(path->reactant1->sym->name, path->reactant2->sym->name) > 0)
    {
      temp_sp = path->reactant1;
      path->reactant1 = path->reactant2;
      path->reactant2 = temp_sp;
      geom = path->orientation1;
      path->orientation1 = path->orientation2;
      path->orientation2 = geom;
      temp_is_complex = path->is_complex[0];
      path->is_complex[0] = path->is_complex[1];
      path->is_complex[1] = temp_is_complex;
    }
    else if (strcmp(path->reactant1->sym->name, path->reactant2->sym->name) == 0)
    {
      if (path->orientation1 < path->orientation2)
      {
        geom = path->orientation1;
        path->orientation1 = path->orientation2;
        path->orientation2 = geom;
        temp_is_complex = path->is_complex[0];
        path->is_complex[0] = path->is_complex[1];
        path->is_complex[1] = temp_is_complex;
      }
    }
  }

  /* Alphabetize if we have three molecules */
  if (reaction->n_reactants == 3)
  {
    if ((path->reactant3->flags&IS_SURFACE)==0)
    {
      if (strcmp(path->reactant1->sym->name, path->reactant3->sym->name) > 0)
      {
         /* Put reactant3 at the beginning */
         temp_sp = path->reactant1;
         geom = path->orientation1;
         path->reactant1 = path->reactant3;
         path->orientation1 = path->orientation3;

         /* Put former reactant1 in place of reactant2 */
         temp_sp2 = path->reactant2;
         geom2 = path->orientation2;
         path->reactant2 = temp_sp;
         path->orientation2 = geom;

         /* Put former reactant2 in place of reactant3 */
         path->reactant3 = temp_sp2;
         path->orientation3 = geom2;
         /* XXX: Update to deal with macromolecules? */

      }
      else if (strcmp(path->reactant2->sym->name, path->reactant3->sym->name) > 0)
      {

         /* Put reactant3 after reactant1 */
         temp_sp = path->reactant2;
         path->reactant2 = path->reactant3;
         path->reactant3 = temp_sp;
         geom = path->orientation2;
         path->orientation2 = path->orientation3;
         path->orientation3 = geom;

      }
    } /*end */
  }
}



/*************************************************************************
 warn_about_high_rates:
    If HIGH_REACTION_PROBABILITY is set to WARNING or ERROR, and the reaction
    probability is high, give the user a warning or error respectively.

 In: parse_state: parser state
     warn_file: The log/error file. Can be stdout/stderr 
     rate_warn: If 1, warn the user about high reaction rates (or give error)
     print_once: If the warning has been printed once, don't repeat it
 Out: print_once. Also print out reaction probabilities (with warning/error)
*************************************************************************/
static int 
warn_about_high_rates(MCELL_STATE *state, FILE *warn_file, int rate_warn, int print_once)
{
  if (rate_warn)
  {
    if (state->notify->high_reaction_prob==WARN_ERROR)
    {
      warn_file = mcell_get_error_file();
      if (!print_once)
      {
        fprintf(warn_file, "\n");
        fprintf(warn_file, "Reaction probabilities generated for the following reactions:\n");
        print_once = 1;
      }
      fprintf(warn_file,"\tError: High ");
    }
    else
    {
      if (!print_once)
      {
        fprintf(warn_file, "\n");
        fprintf(warn_file, "Reaction probabilities generated for the following reactions:\n");
        print_once = 1;
      }
      if (state->notify->high_reaction_prob==WARN_WARN) fprintf(warn_file,"\tWarning: High ");
      else fprintf(warn_file,"\t");
    }
  }
  else 
  {
      if (!print_once)
      {
        fprintf(warn_file, "\n");
        fprintf(warn_file, "Reaction probabilities generated for the following reactions:\n");
        print_once = 1;
      }
      fprintf(warn_file,"\t");
  }
  return print_once;
}



/*************************************************************************
 add_surface_reaction_flags:
 
 In: parse_state: parser state
 Out: Nothing
*************************************************************************/
static void 
add_surface_reaction_flags(MCELL_STATE *state)
{
  struct species *temp_sp;

  /* Add flags for surface reactions with ALL_MOLECULES */
  if (state->all_mols->flags & (CAN_MOLWALL|CAN_GRIDWALL))
  {
    for (int n_mol_bin=0; n_mol_bin<state->mol_sym_table->n_bins; n_mol_bin++)
    {
      for (struct sym_table *symp = state->mol_sym_table->entries[n_mol_bin];
           symp != NULL;
           symp = symp->next)
      {
        temp_sp = (struct species*) symp->value;
        if (temp_sp == state->all_mols) continue;
        if (temp_sp == state->all_volume_mols) continue;
        if (temp_sp == state->all_surface_mols) continue;

        if (((temp_sp->flags & NOT_FREE) == 0) && ((temp_sp->flags & CAN_MOLWALL) == 0))
        {
          temp_sp->flags |= CAN_MOLWALL;
        }
        else if ((temp_sp->flags & ON_GRID) && ((temp_sp->flags & CAN_REGION_BORDER) == 0))
        {
          temp_sp->flags |= CAN_REGION_BORDER;
        }
      }
    }
  }

  /* Add flags for surface reactions with ALL_VOLUME_MOLECULES */
  if (state->all_volume_mols->flags & CAN_MOLWALL)
  {
    for (int n_mol_bin=0; n_mol_bin<state->mol_sym_table->n_bins; n_mol_bin++)
    {
      for (struct sym_table *symp = state->mol_sym_table->entries[n_mol_bin];
           symp != NULL;
           symp = symp->next)
      {
        temp_sp = (struct species*) symp->value;
        if (temp_sp == state->all_mols) continue;
        if (temp_sp == state->all_volume_mols) continue;
        if (temp_sp == state->all_surface_mols) continue;
        if (((temp_sp->flags & NOT_FREE) == 0) && ((temp_sp->flags & CAN_MOLWALL) == 0))
        {
          temp_sp->flags |= CAN_MOLWALL;
        }
      }
    }
  }

  /* Add flags for surface reactions with ALL_SURFACE_MOLECULES */
  if (state->all_surface_mols->flags & CAN_GRIDWALL)
  {
    for (int n_mol_bin=0; n_mol_bin<state->mol_sym_table->n_bins; n_mol_bin++)
    {
      for (struct sym_table *symp = state->mol_sym_table->entries[n_mol_bin];
           symp != NULL;
           symp = symp->next)
      {
        temp_sp = (struct species*) symp->value;
        if (temp_sp == state->all_mols) continue;
        if (temp_sp == state->all_volume_mols) continue;
        if (temp_sp == state->all_surface_mols) continue;
        if (((temp_sp->flags & ON_GRID) && ((temp_sp->flags & CAN_REGION_BORDER) == 0)))
        {
          temp_sp->flags |= CAN_REGION_BORDER;
        }
      }
    }
  }
}

/*************************************************************************
 scale_probabilities:
 
  Scale probabilities, notifying and warning as appropriate.

 In: path: Parse-time structure for reaction pathways
     rx: Pathways leading away from a given intermediate
     parse_state: parser state
     pb_factor:
 Out: Return 1 if rates are high and HIGH_REACTION_PROBABILITY is set to ERROR
 Note: This does not work properly right now. Even if rates are high and
       HIGH_REACTION_PROBABILITY is set to ERROR, the error is ignored
*************************************************************************/
static int 
scale_probabilities(MCELL_STATE *state, struct pathway *path, struct rxn *rx, 
    double pb_factor)
{
  int print_once = 0;  /* flag */
  FILE *warn_file;
  int is_gigantic;
  double rate;

  for (int n_pathway=0;path != NULL;n_pathway++, path = path->next)
  {
    int rate_notify=0, rate_warn=0;
    if (rx->cum_probs[n_pathway]==GIGANTIC) is_gigantic=1;
    else is_gigantic=0;

    /* automatic surface reactions will be printed out from 'init_sim()'. */
    if (is_gigantic) continue;

    if (! rx->rates  ||  ! rx->rates[n_pathway])
      rate = pb_factor*rx->cum_probs[n_pathway];
    else
      rate = 0.0;
    rx->cum_probs[n_pathway] = rate;

    if ((state->notify->reaction_probabilities==NOTIFY_FULL && ((rate>=state->notify->reaction_prob_notify) || (state->notify->reaction_prob_notify==0.0))))
      rate_notify = 1;
    if ((state->notify->high_reaction_prob != WARN_COPE && ((rate>=state->notify->reaction_prob_warn) || ((state->notify->reaction_prob_warn==0.0)))))
      rate_warn = 1;

    if ((rate > 1.0) && (!state->reaction_prob_limit_flag))
    {
      state->reaction_prob_limit_flag = 1;
    }


    if (rate_warn || rate_notify)
    {

      warn_file = mcell_get_log_file();

      print_once = warn_about_high_rates(state, warn_file, rate_warn, print_once);

      if (rx->rates  &&  rx->rates[n_pathway])
        fprintf(warn_file,"Varying probability \"%s\" set for ", rx->rates[n_pathway]->name);
      else
        fprintf(warn_file,"Probability %.4e set for ",rate);
      if (rx->n_reactants==1) fprintf(warn_file,"%s{%d} -> ",rx->players[0]->sym->name,rx->geometries[0]);
      else if (rx->n_reactants == 2)
      {
        if (rx->players[1]->flags & IS_SURFACE)
        {
          fprintf(warn_file,"%s{%d} @ %s{%d} -> ",
                  rx->players[0]->sym->name,rx->geometries[0],
                  rx->players[1]->sym->name,rx->geometries[1]);
         }
         else
         {
           fprintf(warn_file,"%s{%d} + %s{%d} -> ",
                   rx->players[0]->sym->name,rx->geometries[0],
                   rx->players[1]->sym->name,rx->geometries[1]);
         }
      }
      else
      {
        if (rx->players[2]->flags & IS_SURFACE)
        {
          fprintf(warn_file,"%s{%d} + %s{%d}  @ %s{%d} -> ",
                  rx->players[0]->sym->name,rx->geometries[0],
                  rx->players[1]->sym->name,rx->geometries[1],
                  rx->players[2]->sym->name,rx->geometries[2]);
        }
        else
        {
          fprintf(warn_file,"%s{%d} + %s{%d}  + %s{%d} -> ",
                  rx->players[0]->sym->name,rx->geometries[0],
                  rx->players[1]->sym->name,rx->geometries[1],
                  rx->players[2]->sym->name,rx->geometries[2]);
        }
      }
      if (path->product_head == NULL)
      {
        fprintf(warn_file,"NULL ");
      }
      else
      {
        for (struct product *prod = path->product_head ; prod != NULL ; prod = prod->next)
        {
         fprintf(warn_file,"%s{%d} ",prod->prod->sym->name, prod->orientation);
        }
      }

      fprintf(warn_file,"\n");

      if (rate_warn && state->notify->high_reaction_prob==WARN_ERROR)
        return 1;
    }
  }
  return 0;
}



/*************************************************************************
 equivalent_geometry_for_two_reactants:

 In: o1a: orientation of the first reactant from first reaction
     o1b: orientation of the second reactant from first reaction
     o2a: orientation of the first reactant from second reaction
     o2b: orientation of the second reactant from second reaction
 Out: Returns 1 if the two pathways (defined by pairs o1a-o1b and o2a-o2b)
      have equivalent geometry, 0 otherwise.
*************************************************************************/
static int equivalent_geometry_for_two_reactants(int o1a, int o1b, int o2a, int o2b)
{

    /* both reactants for each pathway are in the same
       orientation class and parallel one another */
    if ((o1a == o1b) && (o2a == o2b))
    {
       return 1;
    /* both reactants for each pathway are in the same
       orientation class and opposite one another */
    }
    else if ((o1a == -o1b) && (o2a == -o2b))
    {
       return 1;
    }
    /* reactants are not in the same orientation class */
    if (abs(o1a) != abs(o1b))
    {
       if ((abs(o2a) != abs(o2b)) || ((o2a == 0) && (o2b == 0)))
       {
          return 1;
       }
    }
    if (abs(o2a) != abs(o2b))
    {
       if ((abs(o1a) != abs(o1b)) || ((o1a == 0) && (o1b == 0)))
       {
          return 1;
       }
    }

    return 0;
}

/*************************************************************************
 equivalent_geometry:

 In: p1, p2: pathways to compare
     n: The number of reactants for the pathways
 Out: Returns 1 if the two pathways are the same (i.e. have equivalent
      geometry), 0 otherwise.
*************************************************************************/
static int equivalent_geometry(struct pathway *p1, struct pathway *p2, int n)
{

  short o11,o12,o13,o21,o22,o23; /* orientations of individual reactants */
  /* flags for 3-reactant reactions signaling whether molecules orientations
   * are parallel one another and molecule and surface orientaions are parallel
   * one another
   */
  int mols_parallel_1 = SHRT_MIN + 1; /* for first pathway */
  int mols_parallel_2 = SHRT_MIN + 2; /* for second pathway */
  int mol_surf_parallel_1 = SHRT_MIN + 3; /* for first pathway */
  int mol_surf_parallel_2 = SHRT_MIN + 4; /* for second pathway */

  if (memcmp(p1->is_complex, p2->is_complex, 3))
    return 0;

  if (n < 2)
  {
     /* one reactant case */
     /* RULE: all one_reactant pathway geometries are equivalent */

      return 1;

  }
  else if (n < 3)
  {
    /* two reactants case */

    /* RULE - Two pathways have equivalent geometry when:
       1) Both pathways have exactly the same number of reactants;
       2) There exists an identity mapping between reactants from Pathway 1 and
          Pathway 2 such that for each pair of reactants, r1a and r1b from Pathway
          1, and r2a, and r2b from Pathway 2:
         - r1a is the same species as r2a (likewise for r1b and r2b);
         - r1a and r1b have the same orientation in the same orientation class
           if and only if r2a and r2b do;
         - r1a and r1b have the opposite orientation in the same orientation
           class if and only if r2a and r2b do;
         - r1a and r1b are not in the same orientation class, either because
           they have different orientation classes or both are in the zero
           orientation class, if and only if r2a and r2b are not in the same
           orientation class or both are in the zero orientation class
     */

    o11 = p1->orientation1;
    o12 = p1->orientation2;
    o21 = p2->orientation1;
    o22 = p2->orientation2;

    o13 = o23 = 0;


    return equivalent_geometry_for_two_reactants(o11, o12, o21, o22);

  }
  else if (n < 4)
  {
     /* three reactants case */

    o11 = p1->orientation1;
    o12 = p1->orientation2;
    o13 = p1->orientation3;
    o21 = p2->orientation1;
    o22 = p2->orientation2;
    o23 = p2->orientation3;

    /* special case: two identical reactants */
      if ((p1->reactant1 == p1->reactant2)
          && (p2->reactant1 == p2->reactant2))
      {

       /* Case 1: two molecules and surface are in the same orientation class */
        if ((abs(o11) == abs(o12)) && (abs(o11) == abs(o13)))
        {
          if (o11 == o12) mols_parallel_1 = 1;
          else mols_parallel_1 = 0;

          if (mols_parallel_1)
          {
            if ((o11 == -o13) || (o12 == -o13))
            {
               mol_surf_parallel_1 = 0;
            }
            else 
            {
               mol_surf_parallel_1 = 1;
            }
          }
          else 
          {
               mol_surf_parallel_1 = 0;
          }

          if ((abs(o21) == abs(o22)) && (abs(o21) == abs(o23)))
          {
             if (o21 == o22) mols_parallel_2 = 1;
             else mols_parallel_2 = 0;

             if (mols_parallel_2)
             {
               if ((o21 == -o23) || (o22 == -o23))
               {
                  mol_surf_parallel_2 = 0;
               }
               else 
               {
                  mol_surf_parallel_2 = 1;
               }
             }
             else 
             {
                  mol_surf_parallel_2 = 0;
             }

          }

          if ((mols_parallel_1 == mols_parallel_2) &&
              (mol_surf_parallel_1 == mol_surf_parallel_2))
              {
                 return 1;
          }

       } /* end case 1 */

       /* Case 2: one molecule and surface are in the same orientation class */
       else if ((abs(o11) == abs(o13)) || (abs(o12) == abs(o13)))
       {
          if ((o11 == o13) || (o12 == o13)) mol_surf_parallel_1 = 1;
          else mol_surf_parallel_1 = 0;

          /* check that pathway2 is also in the case2 */

          if ((abs(o21) != abs(o23)) || (abs(o22) != abs(o23)))
          {
             if ((abs(o21) == abs(o23)) || (abs(o22) == abs(o23)))
             {
                if ((o21 == o23) || (o22 == o23)) mol_surf_parallel_2 = 1;
                else mol_surf_parallel_2 = 0;

             }
          }
          if (mol_surf_parallel_1 == mol_surf_parallel_2)
          {
             return 1;
          }

       } /* end case 2 */

       /* Case 3: two molecules but not surface are in the same
                  orientation class */
       else if ((abs(o11) == abs(o12)) && (abs(o11) != abs(o13)))
       {
          if (o11 == o12) mols_parallel_1 = 1;
          else mols_parallel_1 = 0;

          if ((abs(o21) == abs(o22)) && (abs(o21) != abs(o23)))
          {
             if (o21 == o22) mols_parallel_2 = 1;
             else mols_parallel_2 = 0;
          }
          if (mols_parallel_1 == mols_parallel_2)
          {
                 return 1;
          }

       }
       /* Case 4: all molecules and surface are in different orientation classes */
       else if ((abs(o11) != abs(o13)) && (abs(o12) != abs(o13)) &&
                 (abs(o11) != abs(o12)))
                 {

          if ((abs(o21) != abs(o23)) && (abs(o22) != abs(o23)) &&
                 (abs(o21) != abs(o22)))
                 {
               return 1;
          }
       } /* end all cases */

    }
    else { /* no identical reactants */

       if ((equivalent_geometry_for_two_reactants(o11, o12, o21, o22))
           && (equivalent_geometry_for_two_reactants(o12, o13, o22, o23))
           && (equivalent_geometry_for_two_reactants(o11, o13, o21, o23)))
           {
                return 1;
       }

    }

  } // end if (n < 4)


  return 0;
}

/*************************************************************************
 create_sibling_reaction:
    Create a sibling reaction to the given reaction -- a reaction into which
    some of the pathways may be split by split_reaction.

 In:  rx:   reaction for whom to create sibling
 Out: sibling reaction, or NULL on error
*************************************************************************/
static struct rxn *create_sibling_reaction(struct rxn *rx)
{

  struct rxn *reaction = CHECKED_MALLOC_STRUCT(struct rxn, "reaction");
  if (reaction == NULL)
    return NULL;
  reaction->next = NULL;
  reaction->sym = rx->sym;
  reaction->n_reactants = rx->n_reactants;
  reaction->n_pathways = 0;
  reaction->cum_probs = NULL;
  reaction->product_idx = NULL;
  reaction->rates = NULL;
  reaction->max_fixed_p = 0.0;
  reaction->min_noreaction_p = 0.0;
  reaction->pb_factor = 0.0;
  reaction->players = NULL;
  reaction->geometries = NULL;
  reaction->is_complex = NULL;
  reaction->n_occurred = 0;
  reaction->n_skipped = 0.0;
  reaction->prob_t = NULL;
  reaction->pathway_head = NULL;
  reaction->info = NULL;
  return reaction;
}

/*************************************************************************
 split_reaction:
 In:  parse_state: parser state
      rx: reaction to split
 Out: Returns head of the linked list of reactions where each reaction
      contains only geometrically equivalent pathways
*************************************************************************/
static struct rxn *split_reaction(struct rxn *rx)
{
  struct rxn  *curr_rxn_ptr = NULL,  *head = NULL, *end = NULL;
  struct rxn *reaction;
  struct pathway *to_place, *temp;

  /* keep reference to the head of the future linked_list */
  head = end = rx;
  to_place = head->pathway_head->next;
  head->pathway_head->next = NULL;
  head->n_pathways = 1;
  while (to_place != NULL)
  {
    if (to_place->flags & (PATHW_TRANSP | PATHW_REFLEC | PATHW_ABSORP | PATHW_CLAMP_CONC))
    {
      reaction = create_sibling_reaction(rx);
      if (reaction == NULL)
        return NULL;

      reaction->pathway_head = to_place;
      to_place = to_place->next;
      reaction->pathway_head->next = NULL;
      ++ reaction->n_pathways;

      end->next = reaction;
      end = reaction;
    }
    else
    {
      for (curr_rxn_ptr = head; curr_rxn_ptr != NULL; curr_rxn_ptr = curr_rxn_ptr->next)
      {
        if (curr_rxn_ptr->pathway_head->flags & (PATHW_TRANSP | PATHW_REFLEC | PATHW_ABSORP))
          continue;
        if (equivalent_geometry(to_place, curr_rxn_ptr->pathway_head, curr_rxn_ptr->n_reactants))
          break;
      }

      if (! curr_rxn_ptr)
      {
        reaction = create_sibling_reaction(rx);
        if (reaction == NULL)
          return NULL;

        end->next = reaction;
        end = reaction;

        curr_rxn_ptr = end;
      }

      temp = to_place;
      to_place = to_place->next;

      temp->next = curr_rxn_ptr->pathway_head;
      curr_rxn_ptr->pathway_head = temp;
      ++ curr_rxn_ptr->n_pathways;
    }
  }

  return head;
}

/*************************************************************************
 check_reaction_for_duplicate_pathways:
 In:  head: head of linked list of pathways
 Out: Sorts linked list of pathways in alphabetical order according to the
      "prod_signature" field.  Checks for the duplicate pathways.  Prints error
      message and exits simulation if duplicates found.
 Note: This function is called after 'split_reaction()' function so all
       pathways have equivalent geometry from the reactant side.  Here we check
       whether relative orientation of all players (both reactants and
       products) is the same for the two seemingly identical pathways.
 RULE: Two reactions pathways are duplicates if and only if
        (a) they both have the same number and species of reactants;
        (b) they both have the same number and species of products;
        (c) there exists a bijective mapping between the reactants and products
            of the two pathways such that reactants map to reactants, products
            map to products, and the two pathways have equivalent geometry
            under mapping.
            Two pathways R1 and R2 have an equivalent geometry under a mapping
            M if and only if for every pair of players "i" and "j" in R1, the
            corresponding players M(i) and M(j) in R2 have the same orientation
            relation as do "i" and "j" in R1.
            Two players "i" and "j" in a reaction pathway have the following
            orientation:
              parallel - if both "i" and "j" are in the same nonzero orientation
              class with the same sign;
              antiparallel (opposite) - if they are both in the same nonzero
              orientation class but have opposite sign;
              independent - if they are in different orientation classes or both
              in the zero orientation class.

 PostNote: In this function we check only the validity of Rule (c) since
           conditions of Rule (a) and (b) are already satisfied when the
           function is called.
*************************************************************************/
static void check_reaction_for_duplicate_pathways(struct pathway **head)
{

   struct pathway *result = NULL; /* build the sorted list here */
   struct pathway *null_result = NULL; /* put pathways with NULL
                                          prod_signature field here */
   struct pathway *current, *next, **pprev;
   struct product *iter1, *iter2;
   int pathways_equivalent;  /* flag */
   int i, j;
   int num_reactants; /* number of reactants in the pathway */
   int num_products; /* number of products in the pathway */
   int num_players; /* total number of reactants and products in the pathway */
   int *orient_players_1, *orient_players_2; /* array of orientations of players */
   int o1a, o1b, o2a,o2b;

   /* extract  pathways with "prod_signature" field equal to NULL
     into "null_result" list */
   current = *head;
   pprev = head;
   while (current != NULL)
   {
     if (current->prod_signature == NULL)
     {
       *pprev = current->next;
       current->next = null_result;
       null_result = current;
       current = *pprev;
     }
     else
     {
       pprev = &current->next;
       current = current->next;
     }
   }

   /* check for duplicate pathways in null_result */
     current = null_result;
     if ((current != NULL) && (current->next != NULL))
     {
       /* From the previously called function "split_reaction()"
          we know that reactant-reactant pairs in two pathways
          are equivalent. Because there are no products the pathways
          are duplicates.
          RULE: There may be no more than one pathway with zero (--->NULL)
                products in the reaction->pathway_head
                after calling the function "split_reaction()"
       */
       if (current->reactant2 == NULL)
         mcell_error("Exact duplicates of reaction %s  ----> NULL are not allowed.  Please verify that orientations of reactants are not equivalent.",
                     current->reactant1->sym->name);
       else if (current->reactant3 == NULL)
         mcell_error("Exact duplicates of reaction %s + %s  ----> NULL are not allowed.  Please verify that orientations of reactants are not equivalent.",
                     current->reactant1->sym->name,
                     current->reactant2->sym->name);
       else
         mcell_error("Exact duplicates of reaction %s + %s + %s  ----> NULL are not allowed.  Please verify that orientations of reactants are not equivalent.",
                     current->reactant1->sym->name,
                     current->reactant2->sym->name,
                     current->reactant3->sym->name);
     }

    /* now sort the remaining pathway list by "prod_signature" field
       and check for the duplicates */
     current = *head;

  while(current != NULL)
  {
     next = current->next;

     /* insert in sorted order into the "result" */
     if (result == NULL || (strcmp(result->prod_signature, current->prod_signature) >= 0))
     {
        current->next = result;
        result = current;
     }
     else 
     {
        struct pathway *iter = result;
        while(iter->next != NULL && (strcmp(iter->next->prod_signature, current->prod_signature) < 0))
        {
             iter = iter->next;
        }
        current->next = iter->next;
        iter->next = current;
     }

     /* move along the original list */
     current = next;
  }

   /* Now check for the duplicate pathways */
   /* Since the list is sorted we can proceed down the list
      and compare the adjacent nodes */

   current = result;

   if (current != NULL)
   {
     while(current->next != NULL) 
     {
       if (strcmp(current->prod_signature, current->next->prod_signature) == 0)
       {

         pathways_equivalent  = 1;
         /* find total number of players in the pathways */
         num_reactants = 0;
         num_products = 0;
         if (current->reactant1 != NULL) num_reactants++;
         if (current->reactant2 != NULL) num_reactants++;
         if (current->reactant3 != NULL) num_reactants++;

         iter1 = current->product_head;
         while(iter1 != NULL)
         {
           num_products++;
           iter1 = iter1->next;
         }

         num_players = num_reactants + num_products;

         /* create arrays of players orientations */
         orient_players_1 = CHECKED_MALLOC_ARRAY(int, num_players, "reaction player orientations");
         if (orient_players_1 == NULL)
           mcell_die();
         orient_players_2 = CHECKED_MALLOC_ARRAY(int, num_players, "reaction player orientations");
         if (orient_players_2 == NULL)
           mcell_die();

         if (current->reactant1!=NULL) orient_players_1[0]=current->orientation1;
         if (current->reactant2!=NULL) orient_players_1[1]=current->orientation2;
         if (current->reactant3!=NULL) orient_players_1[2]=current->orientation3;
         if (current->next->reactant1!=NULL) orient_players_2[0]=current->next->orientation1;
         if (current->next->reactant2!=NULL) orient_players_2[1]=current->next->orientation2;
         if (current->next->reactant3!=NULL) orient_players_2[2]=current->next->orientation3;


         iter1 = current->product_head;
         iter2 = current->next->product_head;

         for (i = num_reactants; i < num_players; i++)
         {
           orient_players_1[i] = iter1->orientation;
           orient_players_2[i] = iter2->orientation;
           iter1 = iter1->next;
           iter2 = iter2->next;
         }


         /* below we will compare only reactant-product
            and product-product combinations
            because reactant-reactant combinations
            were compared previously in the function
            "equivalent_geometry()"
            */

         /* Initial assumption - pathways are equivalent.
            We check whether this assumption is
            valid by  comparing pairs as described
            above */

         i = 0;
         while((i < num_players) && (pathways_equivalent))
         {
           if (i < num_reactants)
           {
             j = num_reactants;
           }
           else 
           {
             j = i + 1;
           }
           for (; j < num_players; j++)
           {
             o1a = orient_players_1[i];
             o1b = orient_players_1[j];
             o2a = orient_players_2[i];
             o2b = orient_players_2[j];
             if (!equivalent_geometry_for_two_reactants(o1a, o1b, o2a, o2b))
             {
               pathways_equivalent = 0;
               break;
             }
           }
           i++;
         }

         if (pathways_equivalent)
         {
           if (current->reactant2 == NULL)
             mcell_error("Exact duplicates of reaction %s  ----> %s are not allowed.  Please verify that orientations of reactants are not equivalent.",
                         current->reactant1->sym->name,
                         current->prod_signature);
           else if (current->reactant3 == NULL)
             mcell_error("Exact duplicates of reaction %s + %s  ----> %s are not allowed.  Please verify that orientations of reactants are not equivalent.",
                         current->reactant1->sym->name,
                         current->reactant2->sym->name,
                         current->prod_signature);
           else
             mcell_error("Exact duplicates of reaction %s + %s + %s  ----> %s are not allowed.  Please verify that orientations of reactants are not equivalent.",
                         current->reactant1->sym->name,
                         current->reactant2->sym->name,
                         current->reactant3->sym->name,
                         current->prod_signature);
         }
       }

       current = current->next;
     }
   }

    if (null_result == NULL)
    {
       *head = result;
    }
    else if (result == NULL)
    {
       *head = null_result;
    }
    else
    {
       current = result;
       while(current->next != NULL)
       {
          current = current->next;
       }
       current->next = null_result;
       null_result->next = NULL;

       *head = result;
    }

}

/*************************************************************************
 reaction_has_complex_rates:
    Check if a reaction has any complex rates.


    In:  struct rxn *rx - the reaction to check
    Out: 1 if the reaction has complex pathways, 0 otherwise
*************************************************************************/
static int reaction_has_complex_rates(struct rxn *rx)
{
  struct pathway *path;
  for (path = rx->pathway_head;
       path != NULL;
       path = path->next)
  {
    if (path->km_complex)
      return 1;
  }

  return 0;
}

/*************************************************************************
 reorder_varying_pathways:
    Sort pathways so that all complex rates come at the end.  This allows us to
    quickly determine whether a reaction definitely occurs, definitely does not
    occur, or may occur depending on the states of the subunits in the complex.

    In:  struct rxn *rx - the reaction whose pathways to sort
    Out: 1 if the reaction has complex pathways, 0 otherwise

    XXX: Worthwhile sorting pathways by probability?
*************************************************************************/
static int reorder_varying_pathways(struct rxn *rx)
{

  int num_fixed = 0, num_varying = 0;
  int num_fixed_players = 0, num_varying_players = 0;
  int pathway_idx;
  int already_sorted = 1;

  /* If we have no rates, we're done */
  if (! rx->rates)
    return 0;

  /* Count fixed and varying pathways and players */
  for (pathway_idx = 0; pathway_idx < rx->n_pathways; ++ pathway_idx)
  {
    int player_count = rx->product_idx[pathway_idx+1] - rx->product_idx[pathway_idx];
    if (! rx->rates[pathway_idx])
    {
      ++ num_fixed;
      num_fixed_players += player_count;
      if (num_varying) already_sorted = 0;
    }
    else
    {
      ++ num_varying;
      num_varying_players += player_count;
    }
  }

  /* If all are fixed or all are varying, we're done */
  if (! num_fixed  || ! num_varying)
    return 0;

  /* If all fixed pathways already precede all varying pathways, we're done
   */
  if (already_sorted)
    return 0;

  /* Allocate space for sorted info */
  int pathway_mapping[rx->n_pathways];
  struct species **newplayers = NULL;
  short *newgeometries = NULL;
  unsigned char *new_is_complex = NULL;
  u_int *new_product_index = NULL;
  double *new_cum_probs = NULL;
  struct complex_rate **new_complex_rates = NULL;
  struct pathway_info *new_pathway_info = NULL;

  if ((newplayers = CHECKED_MALLOC_ARRAY(struct species*, rx->product_idx[rx->n_pathways], "reaction players array")) == NULL)
    goto failure;
  if ((newgeometries = CHECKED_MALLOC_ARRAY(short, rx->product_idx[rx->n_pathways], "reaction geometries array")) == NULL)
    goto failure;
  if (rx->is_complex)
    if ((new_is_complex = CHECKED_MALLOC_ARRAY(unsigned char, rx->product_idx[rx->n_pathways], "reaction 'is complex' flag array")) == NULL)
      goto failure;
  if ((new_product_index = CHECKED_MALLOC_ARRAY(u_int, rx->product_idx[rx->n_pathways] + 1, "reaction product index array")) == NULL)
    goto failure;
  if ((new_cum_probs = CHECKED_MALLOC_ARRAY(double, rx->n_pathways, "reaction cumulative probabilities array")) == NULL)
    goto failure;
  if ((new_complex_rates = CHECKED_MALLOC_ARRAY(struct complex_rate *, rx->n_pathways, "reaction complex rates array")) == NULL)
    goto failure;
  if ((new_pathway_info = CHECKED_MALLOC_ARRAY(struct pathway_info, rx->n_pathways, "reaction pathway info")) == NULL)
    goto failure;

  memcpy(newplayers, rx->players, rx->n_reactants * sizeof(struct species *));

  /* Now, step through the array until all fixed rates are at the beginning
   */
  int placed_fixed = 0, placed_varying = 0;
  int idx=0;
  int next_player_fixed = rx->n_reactants;
  int next_player_varying = rx->n_reactants + num_fixed_players;
  for (idx = 0; idx < rx->n_pathways; ++ idx)
  {
    int dest_player_idx;
    int dest_pathway;
    int num_players_to_copy = rx->product_idx[idx+1] - rx->product_idx[idx];

    /* Figure out where to put this pathway */
    if (! rx->rates[idx])
    {
      dest_player_idx = next_player_fixed;
      dest_pathway = placed_fixed;

      ++ placed_fixed;
      next_player_fixed += rx->product_idx[idx+1] - rx->product_idx[idx];
    }
    else
    {
      dest_player_idx = next_player_varying;
      dest_pathway = num_fixed + placed_varying;

      ++ placed_varying;
      next_player_varying += rx->product_idx[idx+1] - rx->product_idx[idx];
    }
    pathway_mapping[idx] = next_player_fixed;

    /* Copy everything in */
    memcpy(newplayers + dest_player_idx,
           rx->players + rx->product_idx[idx],
           sizeof(struct species *) * num_players_to_copy);
    memcpy(newgeometries + dest_player_idx,
           rx->geometries + rx->product_idx[idx],
           sizeof(short) * num_players_to_copy);
    if (new_is_complex)
      memcpy(new_is_complex + dest_player_idx,
             rx->is_complex + rx->product_idx[idx],
             sizeof(unsigned char) * num_players_to_copy);
    new_product_index[dest_pathway] = dest_player_idx;
    new_cum_probs[dest_pathway] = rx->cum_probs[idx];
    new_complex_rates[dest_pathway] = rx->rates[idx];
    new_pathway_info[dest_pathway].count = 0.0;
    new_pathway_info[dest_pathway].pathname = rx->info[idx].pathname;
    if (rx->info[idx].pathname) rx->info[idx].pathname->path_num = dest_pathway;
  }
  new_product_index[rx->n_pathways] = rx->product_idx[rx->n_pathways];

  /* Now, fix up varying rates */
  struct t_func *tf;
  for (tf = rx->prob_t; tf != NULL; tf = tf->next)
    tf->path = pathway_mapping[tf->path];

  /* Swap in newly ordered items */
  free(rx->players);
  free(rx->geometries);
  if (rx->is_complex)
    free(rx->is_complex);
  free(rx->product_idx);
  free(rx->cum_probs);
  free(rx->rates);
  free(rx->info);

  rx->players = newplayers;
  rx->geometries = newgeometries;
  rx->is_complex = new_is_complex;
  rx->product_idx = new_product_index;
  rx->cum_probs = new_cum_probs;
  rx->rates = new_complex_rates;
  rx->info = new_pathway_info;

  return 0;

failure:
  if (newplayers) free(newplayers);
  if (newgeometries) free(newgeometries);
  if (new_is_complex) free(new_is_complex);
  if (new_product_index) free(new_product_index);
  if (new_cum_probs) free(new_cum_probs);
  if (new_complex_rates) free(new_complex_rates);
  if (new_pathway_info) free(new_pathway_info);
  return 1;
}

/*************************************************************************
 set_reaction_player_flags:
    Set the reaction player flags for all participating species in this
    reaction.

 In:  rx: the reaction
 Out: Nothing.  Species flags may be updated.
*************************************************************************/
static void set_reaction_player_flags(struct rxn *rx)
{
  switch (rx->n_reactants)
  {
    case 1:
      /* do nothing */
      return;

    case 2:
      if (strcmp(rx->players[0]->sym->name, "ALL_MOLECULES")==0)
      {
          rx->players[0]->flags |= (CAN_MOLWALL|CAN_GRIDWALL);
      }
      else if (strcmp(rx->players[0]->sym->name, "ALL_VOLUME_MOLECULES")==0)
      {
          rx->players[0]->flags |= CAN_MOLWALL;
      }
      else if (strcmp(rx->players[0]->sym->name, "ALL_SURFACE_MOLECULES")==0)
      {
          rx->players[0]->flags |= CAN_GRIDWALL;
      }
      else if ((rx->players[0]->flags & NOT_FREE)==0)
      {
        /* two volume molecules */
        if ((rx->players[1]->flags & NOT_FREE)==0)
        {
          rx->players[0]->flags |= CAN_MOLMOL;
          rx->players[1]->flags |= CAN_MOLMOL;
        }
        /* one volume molecules and one wall */
        else if ((rx->players[1]->flags & IS_SURFACE)!=0)
        {
          rx->players[0]->flags |= CAN_MOLWALL;
        }
        /* one volume molecule and one grid molecule */
        else if ((rx->players[1]->flags & ON_GRID)!= 0)
        {
          rx->players[0]->flags |= CAN_MOLGRID;
        }
      }
      else if ((rx->players[0]->flags & IS_SURFACE)!=0)
      {
        /* one volume molecule and one wall */
        if ((rx->players[1]->flags & NOT_FREE)==0)
        {
          rx->players[1]->flags |= CAN_MOLWALL;
        }
        /* one grid molecule and one wall */
        else if ((rx->players[1]->flags & ON_GRID) != 0)
        {
          rx->players[1]->flags |= CAN_GRIDWALL;
        }
      }
      else if ((rx->players[0]->flags & ON_GRID)!= 0)
      {
        /* one volume molecule and one grid molecule */
        if ((rx->players[1]->flags & NOT_FREE)==0)
        {
          rx->players[1]->flags |= CAN_MOLGRID;
        }
        /* two grid molecules */
        else if ((rx->players[1]->flags & ON_GRID) != 0)
        {
          rx->players[0]->flags |= CAN_GRIDGRID;
          rx->players[1]->flags |= CAN_GRIDGRID;
        }
        /* one grid molecule and one wall */
        else if ((rx->players[1]->flags & IS_SURFACE) != 0)
        {
          rx->players[0]->flags |= CAN_GRIDWALL;
        }
      }
      break;

    case 3:
      if ((rx->players[2]->flags & IS_SURFACE) != 0)
      {
        /* two molecules and surface  */
        if ((rx->players[0]->flags & NOT_FREE)==0)
        {
          /* one volume molecule, one grid molecule, one surface */
          if ((rx->players[1]->flags & ON_GRID)!= 0)
          {
            rx->players[0]->flags |= CAN_MOLGRID;
          }
        }
        else if ((rx->players[0]->flags & ON_GRID)!= 0)
        {
          /* one volume molecule, one grid molecule, one surface */
          if ((rx->players[1]->flags & NOT_FREE)==0)
          {
            rx->players[1]->flags |= CAN_MOLGRID;
          }
          /* two grid molecules, one surface */
          else if ((rx->players[1]->flags & ON_GRID) != 0)
          {
            rx->players[0]->flags |= CAN_GRIDGRID;
            rx->players[1]->flags |= CAN_GRIDGRID;
          }
        }
      }
      else
      {
        if ((rx->players[0]->flags & NOT_FREE)==0)
        {
          if ((rx->players[1]->flags & NOT_FREE)==0)
          {
            /* three volume molecules */
            if ((rx->players[2]->flags & NOT_FREE)==0)
            {
              rx->players[0]->flags |= CAN_MOLMOLMOL;
              rx->players[1]->flags |= CAN_MOLMOLMOL;
              rx->players[2]->flags |= CAN_MOLMOLMOL;
            }
            /* two volume molecules and one grid molecule */
            else if ((rx->players[2]->flags & ON_GRID) !=0)
            {
              rx->players[0]->flags |= CAN_MOLMOLGRID;
              rx->players[1]->flags |= CAN_MOLMOLGRID;
            }
          }
          else if ((rx->players[1]->flags & ON_GRID) !=0)
          {
            /* two volume molecules and one grid molecule */
            if ((rx->players[2]->flags & NOT_FREE)==0)
            {
              rx->players[0]->flags |= CAN_MOLMOLGRID;
              rx->players[2]->flags |= CAN_MOLMOLGRID;
            }
            /* one volume molecules and two grid molecules */
            else if ((rx->players[2]->flags & ON_GRID) !=0)
            {
              rx->players[0]->flags |= CAN_MOLGRIDGRID;
            }
          }
        }
        else if ((rx->players[0]->flags & ON_GRID) != 0)
        {
          if ((rx->players[1]->flags & NOT_FREE)==0)
          {
            /* two volume molecules and one grid molecule */
            if ((rx->players[2]->flags & NOT_FREE)==0)
            {
              rx->players[1]->flags |= CAN_MOLMOLGRID;
              rx->players[2]->flags |= CAN_MOLMOLGRID;
            }
            /* one volume molecule and two grid molecules */
            else if ((rx->players[2]->flags & ON_GRID) !=0)
            {
              rx->players[1]->flags |= CAN_MOLGRIDGRID;
            }
          }
          else if ((rx->players[1]->flags & ON_GRID) !=0)
          {
            /* one volume molecule and two grid molecules */
            if ((rx->players[2]->flags & NOT_FREE)==0)
            {
              rx->players[2]->flags |= CAN_MOLGRIDGRID;
            }
            /* three grid molecules */
            else if ((rx->players[2]->flags & ON_GRID) !=0)
            {
              rx->players[0]->flags |= CAN_GRIDGRIDGRID;
              rx->players[1]->flags |= CAN_GRIDGRIDGRID;
              rx->players[2]->flags |= CAN_GRIDGRIDGRID;
            }
          }
        }
      }
      break;

    default:
      //assert(0);
      break;
  }
}

/*************************************************************************
 build_reaction_hash_table:
    Scan the symbol table, copying all reactions found into the reaction hash.

 In:  parse_state: parser state
      num_rx: num reactions expected
 Out: 0 on success, 1 if we fail to allocate the table
*************************************************************************/
static int build_reaction_hash_table(MCELL_STATE *state, int num_rx)
{
  struct rxn **rx_tbl = NULL;
  int rx_hash;
  for (rx_hash=2; rx_hash<=num_rx && rx_hash != 0; rx_hash <<= 1)
    ;
  rx_hash <<= 1;

  if (rx_hash == 0) rx_hash = MAX_RX_HASH;
  if (rx_hash > MAX_RX_HASH) rx_hash = MAX_RX_HASH;
#ifdef REPORT_RXN_HASH_STATS
  mcell_log("Num rxns: %d", num_rx);
  mcell_log("Size of hash: %d", rx_hash);
#endif

  /* Create the reaction hash table */
  state->rx_hashsize = rx_hash;
  rx_hash -= 1;
  rx_tbl = CHECKED_MALLOC_ARRAY(struct rxn*, state->rx_hashsize, "reaction hash table");
  if (rx_tbl==NULL)
     return 1;
  state->reaction_hash = rx_tbl;
  for (int i=0;i<=rx_hash;i++) rx_tbl[i] = NULL;

#ifdef REPORT_RXN_HASH_STATS
  int numcoll = 0;
#endif
  for (int i=0;i<state->rxn_sym_table->n_bins;i++)
  {
    for (struct sym_table *sym = state->rxn_sym_table->entries[i]; sym != NULL; sym = sym->next)
    {
      if (sym == NULL) continue;

      struct rxn *rx = (struct rxn*) sym->value;
      int table_slot;
      if (rx->n_reactants == 1)
      {
        table_slot = rx->players[0]->hashval & rx_hash;
      }
      else
      {
        table_slot = (rx->players[0]->hashval + rx->players[1]->hashval) & rx_hash;
      }


#ifdef REPORT_RXN_HASH_STATS
      if (rx_tbl[table_slot] != NULL)
      {
        mcell_log("Collision: %s and %s", rx_tbl[table_slot]->sym->name, sym->name);
        ++ numcoll;
      }
#endif
      state->n_reactions++;
      while (rx->next != NULL) rx = rx->next;
      rx->next = rx_tbl[table_slot];
      rx_tbl[table_slot] = (struct rxn*)sym->value;
    }
  }
#ifdef REPORT_RXN_HASH_STATS
  mcell_log("Num collisions: %d", numcoll);
#endif

  return 0;
}

/*************************************************************************
 load_rate_file:
    Read in a time-varying reaction rates file.

 In:  parse_state:  parser state
      rx:    Reaction structure that we'll load the rates into.
      fname: Filename to read the rates from.
      path:  Index of the pathway that these rates apply to.
 Out: Returns 1 on error, 0 on success.
      Rates are added to the prob_t linked list.  If there is a rate given for
      time <= 0, then this rate is stuck into cum_probs and the (time <= 0)
      entries are not added to the list.  If no initial rate is given in the
      file, it is assumed to be zero.
 Note: The file format is assumed to be two columns of numbers; the first
      column is time (in seconds) and the other is rate (in appropriate
      units) that starts at that time.  Lines that are not numbers are
      ignored.
*************************************************************************/
#define RATE_SEPARATORS "\f\n\r\t\v ,;"
#define FIRST_DIGIT "+-0123456789"
static int
load_rate_file(MCELL_STATE *state, struct rxn *rx, char *fname, int path)
{
  int i;
  FILE *f = fopen(fname,"r");

  if (!f) return 1;
  else
  {
    struct t_func *tp,*tp2;
    double t,rate;
    char buf[2048];
    char *cp;
    int linecount = 0;
#ifdef DEBUG
    int valid_linecount = 0;
#endif

    tp2 = NULL;
    while (fgets(buf,2048,f))
    {
      linecount++;
      for (i=0;i<2048;i++) { if (!strchr(RATE_SEPARATORS,buf[i])) break; }

      if (i<2048 && strchr(FIRST_DIGIT,buf[i]))
      {
        t = strtod((buf+i) , &cp);
        if (cp == (buf+i)) continue;  /* Conversion error. */

        for (i=cp-buf ; i<2048 ; i++) { if (!strchr(RATE_SEPARATORS,buf[i])) break; }
        rate = strtod((buf+i) , &cp);
        if (cp == (buf+i)) continue;  /* Conversion error */

        /// XXX: MARKUS - adapt the below warnings
#if 0
        /* at this point we need to handle negative reaction rates */
        if (rate < 0.0)
        {
          if (parse_state->vol->notify->neg_reaction==WARN_ERROR)
          {
            mdlerror(parse_state, "Error: reaction rates should be zero or positive.");
            return 1;
          }
          else if (parse_state->vol->notify->neg_reaction == WARN_WARN) {
            mcell_warn("Warning: negative reaction rate %f; setting to zero and continuing.", rate);
            rate = 0.0;
          }
        }
#endif

        tp = CHECKED_MEM_GET(state->tv_rxn_mem, "time-varying reaction rate");
        if (tp == NULL)
          return 1;
        tp->next = NULL;
        tp->path = path;
        tp->time = t / state->time_unit;
        tp->value = rate;
#ifdef DEBUG
        valid_linecount++;
#endif

        if (rx->prob_t == NULL)
        {
          rx->prob_t = tp;
          tp2 = tp;
        }
        else
        {
          if (tp2==NULL)
          {
            tp2 = tp;
            tp->next = rx->prob_t;
            rx->prob_t = tp;
          }
          else
          {
            if (tp->time < tp2->time)
              mcell_warn("In rate file '%s', line %d is out of sequence.  Resorting.", fname, linecount);
            tp->next = tp2->next;
            tp2->next = tp;
            tp2 = tp;
          }
        }
      }
    }

#ifdef DEBUG
    mcell_log("Read %d rates from file %s.", valid_linecount, fname);
#endif

    fclose(f);
  }
  return 0;
}
#undef FIRST_DIGIT
#undef RATE_SEPARATORS



/*************************************************************************
 prepare_reactions:
    Postprocess the parsed reactions, moving them to the reaction hash table,
    and transferring information from the pathway structures to a more compact,
    runtime-optimized form.

 In: parse_state: parser state
 Out: Returns 1 on error, 0 on success.
      Reaction hash table is built and geometries are set properly.  Unlike in
      the parser, reactions with different reactant geometries are _different
      reactions_, and are stored as separate struct rxns.

 Note: The user inputs _geometric equivalence classes_, but here we convert
       from that to _output construction geometry_.  A geometry of 0 means to
       choose a random orientation.  A geometry of k means to adopt the
       geometry of the k'th species in the list (reactants start at #1,
       products are in order after reactants).  A geometry of -k means to adopt
       the opposite of the geometry of the k'th species.  The first n_reactants
       products determine the fate of the reactants (NULL = destroyed), and the
       rest are real products.
 PostNote: The reactants are used for triggering, and those have equivalence
       class geometry even in here.
*************************************************************************/
int init_reactions(MCELL_STATE *state)
{
  struct pathway *path;
  struct product *prod;
  struct rxn *rx;
  struct t_func *tp;
  //double D_tot, t_step;
  short geom;
  int k, kk;
  /* flags that tell whether reactant_1 is also on the product list,
     same for reactant_2 and reactant_3 */
  int recycled1, recycled2, recycled3;
  int num_rx, num_players;
  struct species *temp_sp;
  int n_prob_t_rxns; /* # of pathways with time-varying rates */
  struct rxn *reaction;

  num_rx = 0;

  state->vacancy_search_dist2 *= state->r_length_unit;         /* Convert units */
  state->vacancy_search_dist2 *= state->vacancy_search_dist2;  /* Take square */

  if (state->rx_radius_3d <= 0.0)
  {
    state->rx_radius_3d = 1.0/sqrt(MY_PI*state->grid_density);
  }
  state->tv_rxn_mem = create_mem(sizeof(struct t_func) , 100);
  if (state->tv_rxn_mem == NULL) return 1;

  for (int n_rxn_bin=0; n_rxn_bin<state->rxn_sym_table->n_bins; n_rxn_bin++)
  {
    for (struct sym_table *sym = state->rxn_sym_table->entries[n_rxn_bin];
         sym != NULL;
         sym = sym->next)
    {
      reaction = (struct rxn*)sym->value;
      reaction->next = NULL;

      for (path=reaction->pathway_head ; path != NULL ; path = path->next)
      {
        check_duplicate_special_reactions(path);

        /* if one of the reactants is a surface, move it to the last reactant.
         * Also arrange reactant1 and reactant2 in alphabetical order */
        if (reaction->n_reactants>1)
        {
          /* Put surface last */
          if ((path->reactant1->flags & IS_SURFACE) != 0)
          {
            temp_sp = path->reactant1;
            path->reactant1 = path->reactant2;
            path->reactant2 = temp_sp;
            geom = path->orientation1;
            path->orientation1 = path->orientation2;
            path->orientation2 = geom;
          }
          if (reaction->n_reactants>2)
          {
            if ((path->reactant2->flags & IS_SURFACE) != 0)
            {
              temp_sp = path->reactant3;
              path->reactant3 = path->reactant2;
              path->reactant2 = temp_sp;
              geom = path->orientation3;
              path->orientation3 = path->orientation2;
              path->orientation2 = geom;
            }
          }
          alphabetize_pathway(path, reaction);
        } /* end if (n_reactants > 1) */

      }  /* end for (path = reaction->pathway_head; ...) */


      /* if reaction contains equivalent pathways, split this reaction into a
       * linked list of reactions each containing only equivalent pathways.
       */

      rx = split_reaction(reaction);

      /* set the symbol value to the head of the linked list of reactions */
      sym->value = (void *)rx;

      while (rx != NULL)
      {
        double pb_factor = 0.0;
        /* Check whether reaction contains pathways with equivalent product
         * lists.  Also sort pathways in alphabetical order according to the
         * "prod_signature" field.
         */
        check_reaction_for_duplicate_pathways(&rx->pathway_head);

        num_rx++;

        /* At this point we have reactions of the same geometry and can collapse them
         * and count how many non-reactant products are in each pathway. */

        /* Search for reactants that appear as products */
        /* Any reactants that don't appear are set to be destroyed. */
        rx->product_idx = CHECKED_MALLOC_ARRAY(u_int, rx->n_pathways+1, "reaction product index array");
        rx->cum_probs = CHECKED_MALLOC_ARRAY(double, rx->n_pathways, "reaction cumulative probabilities array");

        /* Note, that the last member of the array "rx->product_idx"
         * contains size of the array "rx->players" */

        if (rx->product_idx == NULL  || rx->cum_probs == NULL)
          return 1;

        if (reaction_has_complex_rates(rx))
        {
          int pathway_idx;
          rx->rates = CHECKED_MALLOC_ARRAY(struct complex_rate *, rx->n_pathways, "reaction complex rates array");
          if (rx->rates == NULL)
            return 1;
          for (pathway_idx = 0; pathway_idx < rx->n_pathways; ++ pathway_idx)
            rx->rates[pathway_idx] = NULL;
        }

        n_prob_t_rxns = 0;
        path = rx->pathway_head;

        for (int n_pathway=0; path!=NULL ; n_pathway++ , path = path->next)
        {

          rx->product_idx[n_pathway] = 0;
          if (rx->rates)
            rx->rates[n_pathway] = path->km_complex;

          /* Look for concentration clamp */
          if (path->reactant2!=NULL && (path->reactant2->flags&IS_SURFACE)!=0 &&
              path->km >= 0.0 && path->product_head==NULL && ((path->flags & PATHW_CLAMP_CONC) != 0))
          {
            struct ccn_clamp_data *ccd;

            if (n_pathway!=0 || path->next!=NULL)
              mcell_warn("Mixing surface modes with other surface reactions.  Please don't.");

            if (path->km>0)
            {
              ccd = CHECKED_MALLOC_STRUCT(struct ccn_clamp_data, "concentration clamp data");
              if (ccd==NULL)
                return 1;

              ccd->surf_class = path->reactant2;
              ccd->mol = path->reactant1;
              ccd->concentration = path->km;
              if (path->orientation1*path->orientation2==0)
              {
                ccd->orient = 0;
              }
              else
              {
                ccd->orient = (path->orientation1==path->orientation2) ? 1 : -1;
              }
              ccd->sides = NULL;
              ccd->next_mol = NULL;
              ccd->next_obj = NULL;
              ccd->objp = NULL;
              ccd->n_sides = 0;
              ccd->side_idx = NULL;
              ccd->cum_area = NULL;
              ccd->scaling_factor = 0.0;
              ccd->next = state->clamp_list;
              state->clamp_list = ccd;
            }
            path->km = GIGANTIC;
          }
          else if ((path->flags & PATHW_TRANSP) != 0)
          {
            rx->n_pathways = RX_TRANSP;
            if (path->reactant2!=NULL && (path->reactant2->flags&IS_SURFACE) &&
               (path->reactant1->flags & ON_GRID))
            {
                    path->reactant1->flags |= CAN_REGION_BORDER;
            }
          }
          else if ((path->flags & PATHW_REFLEC) != 0)
          {
            rx->n_pathways = RX_REFLEC;
            if (path->reactant2!=NULL && (path->reactant2->flags&IS_SURFACE) &&
               (path->reactant1->flags & ON_GRID))
            {
                    path->reactant1->flags |= CAN_REGION_BORDER;
            }
          }
          else if (path->reactant2!=NULL && (path->reactant2->flags&IS_SURFACE) && (path->reactant1->flags & ON_GRID) && (path->product_head==NULL) && (path->flags & PATHW_ABSORP))
          {
             rx->n_pathways = RX_ABSORB_REGION_BORDER;
             path->reactant1->flags |= CAN_REGION_BORDER;
          }
          else if ((strcmp(path->reactant1->sym->name, "ALL_SURFACE_MOLECULES") == 0))
          {
            if (path->reactant2!=NULL && (path->reactant2->flags&IS_SURFACE)  && (path->product_head==NULL) && (path->flags & PATHW_ABSORP))
            {
              rx->n_pathways = RX_ABSORB_REGION_BORDER;
              path->reactant1->flags |= CAN_REGION_BORDER;
            }
          }
          if (path->km_filename == NULL) rx->cum_probs[n_pathway] = path->km;
          else
          {
            rx->cum_probs[n_pathway]=0;
            n_prob_t_rxns++;
          }

          recycled1 = 0;
          recycled2 = 0;
          recycled3 = 0;

          for (prod=path->product_head ; prod != NULL ; prod = prod->next)
          {
            if (recycled1 == 0 && prod->prod == path->reactant1) recycled1 = 1;
            else if (recycled2 == 0 && prod->prod == path->reactant2) recycled2 = 1;
            else if (recycled3 == 0 && prod->prod == path->reactant3) recycled3 = 1;
            else rx->product_idx[n_pathway]++;
          }


        } /* end for (n_pathway=0,path=rx->pathway_head; ...) */

        /* Now that we know how many products there really are, set the index array */
        /* and malloc space for the products and geometries. */
        num_players = rx->n_reactants;
        kk = rx->n_pathways;
        if (kk<=RX_SPECIAL) kk = 1;
        for (int n_pathway=0;n_pathway<kk;n_pathway++)
        {
          k = rx->product_idx[n_pathway] + rx->n_reactants;
          rx->product_idx[n_pathway] = num_players;
          num_players += k;
        }
        rx->product_idx[kk] = num_players;

        rx->players = CHECKED_MALLOC_ARRAY(struct species*, num_players, "reaction players array");
        rx->geometries = CHECKED_MALLOC_ARRAY(short, num_players, "reaction geometries array");
        if (rx->pathway_head->is_complex[0] ||
            rx->pathway_head->is_complex[1] ||
            rx->pathway_head->is_complex[2])
        {
          rx->is_complex = CHECKED_MALLOC_ARRAY(unsigned char, num_players, "reaction 'is complex' flag");
          if (rx->is_complex == NULL)
            return 1;
          memset(rx->is_complex, 0, sizeof(unsigned char) * num_players);
        }
        else
          rx->is_complex = NULL;

        if (rx->players==NULL || rx->geometries==NULL)
          return 1;

        /* Load all the time-varying rates from disk (if any), merge them into */
        /* a single sorted list, and pull off any updates for time zero. */
        if (n_prob_t_rxns > 0)
        {
          path = rx->pathway_head;
          for (int n_pathway = 0; path!=NULL ; n_pathway++, path=path->next)
          {
            if (path->km_filename != NULL)
            {
              if (load_rate_file(state, rx, path->km_filename, n_pathway))
                mcell_error("Failed to load rates from file '%s'.", path->km_filename);
            }
          }
          rx->prob_t = (struct t_func*) ae_list_sort((struct abstract_element*)rx->prob_t);

          while (rx->prob_t != NULL && rx->prob_t->time <= 0.0)
          {
            rx->cum_probs[ rx->prob_t->path ] = rx->prob_t->value;
            rx->prob_t = rx->prob_t->next;
          }
        } /* end if (n_prob_t_rxns > 0) */


        /* Set the geometry of the reactants.  These are used for triggering.                 */
        /* Since we use flags to control orientation changes, just tell everyone to stay put. */
        path = rx->pathway_head;
        rx->players[0] = path->reactant1;
        rx->geometries[0] = path->orientation1;
        if (rx->is_complex) rx->is_complex[0] = path->is_complex[0];
        if (rx->n_reactants > 1)
        {
          rx->players[1] = path->reactant2;
          rx->geometries[1] = path->orientation2;
          if (rx->is_complex) rx->is_complex[1] = path->is_complex[1];
          if (rx->n_reactants > 2)
          {
            rx->players[2] = path->reactant3;
            rx->geometries[2] = path->orientation3;
            if (rx->is_complex) rx->is_complex[2] = path->is_complex[2];
          }
        }

        /* maximum number of surface products */
        path = rx->pathway_head;
        int max_num_surf_products = set_product_geometries(path, rx, prod);

        pb_factor = compute_pb_factor(state, rx, max_num_surf_products);
        rx->pb_factor = pb_factor;
        path = rx->pathway_head;

        if (scale_probabilities(state, path, rx, pb_factor))
          return 1;

        if (n_prob_t_rxns > 0)
        {
          for (tp = rx->prob_t ; tp != NULL ; tp = tp->next)
            tp->value *= pb_factor;
        }

        /* Move counts from list into array */
        if (rx->n_pathways > 0)
        {
          rx->info = CHECKED_MALLOC_ARRAY(struct pathway_info, rx->n_pathways, "reaction pathway info");
          if (rx->info == NULL)
            return 1;

          path = rx->pathway_head;
          for (int n_pathway=0; path!=NULL ; n_pathway++,path=path->next)
          {
            rx->info[n_pathway].count = 0;
            rx->info[n_pathway].pathname = path->pathname;    /* Keep track of named rxns */
            if (path->pathname!=NULL)
            {
              rx->info[n_pathway].pathname->path_num = n_pathway;
              rx->info[n_pathway].pathname->rx = rx;
            }
          }
        }
        else /* Special reaction, only one exit pathway */
        {
          rx->info = CHECKED_MALLOC_STRUCT(struct pathway_info,
                                               "reaction pathway info");
          if (rx->info == NULL)
            return 1;
          rx->info[0].count = 0;
          rx->info[0].pathname = rx->pathway_head->pathname;
          if (rx->pathway_head->pathname!=NULL)
          {
            rx->info[0].pathname->path_num = 0;
            rx->info[0].pathname->rx = rx;
          }
        }

        /* Sort pathways so all fixed pathways precede all varying pathways */
        if (rx->rates  &&  rx->n_pathways > 0)
          reorder_varying_pathways(rx);

        /* Compute cumulative properties */
        for (int n_pathway=1; n_pathway<rx->n_pathways; ++n_pathway)
          rx->cum_probs[n_pathway] += rx->cum_probs[n_pathway-1];
        if (rx->n_pathways > 0)
          rx->min_noreaction_p = rx->max_fixed_p = rx->cum_probs[rx->n_pathways - 1];
        else
          rx->min_noreaction_p = rx->max_fixed_p = 1.0;
        if (rx->rates)
          for (int n_pathway=0; n_pathway<rx->n_pathways; ++n_pathway)
            if (rx->rates[n_pathway])
              rx->min_noreaction_p += macro_max_rate(rx->rates[n_pathway], pb_factor);

        rx = rx->next;
      }
    }
  }

  if (state->grid_grid_reaction_flag  || state->grid_grid_grid_reaction_flag)
  {
    if (state->notify->reaction_probabilities==NOTIFY_FULL)
      mcell_log("For reaction between two (or three) surface molecules the upper probability limit is given. The effective reaction probability will be recalculated dynamically during simulation.");
  }

  if (build_reaction_hash_table(state, num_rx))
    return 1;

  state->rx_radius_3d *= state->r_length_unit; /* Convert into length units */

  for (int n_rxn_bin=0;n_rxn_bin<state->rx_hashsize;n_rxn_bin++)
  {
    for (struct rxn *this_rx = state->reaction_hash[n_rxn_bin];
         this_rx != NULL;
         this_rx = this_rx->next)
    {
      /* Here we deallocate all memory used for creating pathways. */
      path = this_rx->pathway_head;
      struct pathway *next_path = path;
      while (next_path != NULL)
      {
        next_path = path->next;
        if (path->prod_signature != NULL) 
        {
          free(path->prod_signature);
        }

        struct product *dead_prod = path->product_head;
        struct product *nxt = dead_prod;
        while (nxt != NULL) 
        {
          nxt = dead_prod->next;
          free(dead_prod);
          dead_prod = nxt;
        }

        free(path);
        path = next_path;
      }

      set_reaction_player_flags(this_rx);
      this_rx->pathway_head = NULL;
    }
  }

  add_surface_reaction_flags(state);

  if (state->notify->reaction_probabilities==NOTIFY_FULL)
    mcell_log_raw("\n");

  return 0;
}

