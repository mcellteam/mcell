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

//#include "libmcell.h"
//#include "logging.h"
#include "sym_table.h"
#include "create_reactions.h"



/*************************************************************************
 *
 * extract_reactants extracts the reactant info into a pathway structure
 *
 *************************************************************************/
MCELL_STATUS 
extract_reactants(struct pathway *path, struct species_opt_orient *reactants,
  int *num_reactants, int *num_vol_mols, int *num_grid_mols) 
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
    int *num_vol_mols, int *num_grid_mols) 
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

