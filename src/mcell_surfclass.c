/******************************************************************************
 *
 * Copyright (C) 2006-2017 by
 * The Salk Institute for Biological Studies and
 * Pittsburgh Supercomputing Center, Carnegie Mellon University
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301,
 * USA.
 *
******************************************************************************/

#include <string.h>
#include <stdlib.h>

#include "config.h"

#include "logging.h"
#include "sym_table.h"
#include "mcell_species.h"
#include "init.h"
#include "mcell_init.h"
#include "mcell_structs.h"
#include "mcell_reactions.h"
#include "mcell_surfclass.h"

static int check_valid_molecule_release(
  MCELL_STATE *state, struct mcell_species *mol_type);

/**************************************************************************
 mcell_add_surf_class_properties:
    Add a property to a surface class

 In: state:         the simulation state
     reaction_type: RFLCT, TRANSP, or SINK
     sc_sym:        symbol for surface class
     reactant_sym:  symbol for reactant molecule
     orient:        orientation for molecule
 Out: 0 on success, 1 on failure.
**************************************************************************/
MCELL_STATUS mcell_add_surf_class_properties(
    MCELL_STATE *state, int reaction_type, mcell_symbol *sc_sym,
    mcell_symbol *reactant_sym, short orient) {

  struct species *surf_class = (struct species *)sc_sym->value;
  if (mcell_add_surface_reaction(state->rxn_sym_table, reaction_type,
                                 surf_class, reactant_sym, orient)) {
    return MCELL_FAIL;
  }
  return MCELL_SUCCESS;

}

/**************************************************************************
 mcell_create_surf_class:
    Start a surface class declaration.

 In: state:           the simulation state
     surf_class_name: the surface class name
     sc_sym:          surface class symbol
 Out: 0 on success, 1 on failure. The surface class is created
**************************************************************************/
MCELL_STATUS mcell_create_surf_class(
    MCELL_STATE *state, char *surf_class_name, mcell_symbol **sc_sym) {

  struct sym_entry *sym = NULL;
  int error_code = new_mol_species(state, surf_class_name, &sym);
  if (error_code) {
    return error_code; 
  }

  struct species *surf_class = (struct species *)sym->value;
  surf_class->flags = IS_SURFACE;
  surf_class->refl_mols = NULL;
  surf_class->transp_mols = NULL;
  surf_class->absorb_mols = NULL;
  surf_class->clamp_conc_mols = NULL;

  if (sc_sym != NULL) {
    *sc_sym = sym;
  }

  return MCELL_SUCCESS;

}

/**************************************************************************
 mcell_add_mol_release_to_surf_class:
    Create a surface class "property" to release molecules. The surface class
    still needs to be assigned to a region. Multiple properties can be added to
    a single surface class. 

 In: state:          the simulation state
     sc_sym:         surface class symbol
     sm_info:        info about the surface molecule to be released
     quantity:       the amount of surface molecules to release
     density_or_num: 0 if density, 1 if number
     smd_list:
 Out: The surface molecule data
**************************************************************************/
struct sm_dat *mcell_add_mol_release_to_surf_class(
    MCELL_STATE *state, struct sym_entry *sc_sym, struct mcell_species *sm_info,
    double quantity, int density_or_num, struct sm_dat *smd_list) {

  struct species *species_ptr = (struct species *)sm_info->mol_type->value;
  if (!(species_ptr->flags & ON_GRID)) {
    return NULL;
  } else if (check_valid_molecule_release(state, sm_info)) {
    return NULL;
  }

  struct sm_dat *sm_dat_ptr = CHECKED_MALLOC_STRUCT(
      struct sm_dat, "surface molecule data");
  if (sm_dat_ptr == NULL) {
    return NULL; 
  }

  sm_dat_ptr->next = smd_list;
  sm_dat_ptr->sm = species_ptr;
  sm_dat_ptr->quantity_type = density_or_num;
  sm_dat_ptr->quantity = quantity;
  sm_dat_ptr->orientation = sm_info->orient;

  struct species *sc = (struct species *)sc_sym->value;
  sc->sm_dat_head = sm_dat_ptr;

  return sm_dat_ptr;
}

/**************************************************************************
 check_valid_molecule_release:
    Check that a particular molecule type is valid for inclusion in a release
    site.  Checks that orientations are present if required, and absent if
    forbidden, and that we aren't trying to release a surface class.

 In: state:    the simulation state
     mol_type: molecule species and (optional) orientation for release
 Out: 0 on success, 1 on failure
**************************************************************************/
static int check_valid_molecule_release(MCELL_STATE *state, 
                                        struct mcell_species *mol_type) {

  struct species *mol = (struct species *)mol_type->mol_type->value;
  if (mol->flags & ON_GRID) {
    if (!mol_type->orient_set) {
      if (state->notify->missed_surf_orient == WARN_ERROR) {
        return 1;
      }
    }
  } else if ((mol->flags & NOT_FREE) == 0) {
    if (mol_type->orient_set) {
      if (state->notify->useless_vol_orient == WARN_ERROR) {
        return 1;
      }
    }
  } else {
    return 1;
  }

  return 0;
}

/**************************************************************************
 mcell_assign_surf_class_to_region:
    Assign a surface class to a region

 In: sg_name:  the name of the surface class to be assigned
     rgn:      the region which will have the surface class assigned to it
 Out: 0 on success, 1 on failure. Surface class is assigned
**************************************************************************/
MCELL_STATUS mcell_assign_surf_class_to_region(
    struct sym_entry *sc_sym, struct region *rgn) {

  if (rgn->surf_class != NULL)
    return MCELL_FAIL;
  rgn->surf_class = (struct species *)sc_sym->value;
  return MCELL_SUCCESS;
}
