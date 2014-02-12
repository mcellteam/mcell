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
#include <string.h>

//#include "logging.h"
//#include "mem_util.h"
//#include "sym_table.h"
#include "create_reactions.h"
#include "strfunc.h"


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


