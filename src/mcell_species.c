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

#include "config.h"
#include "logging.h"
#include <math.h>
#include <stdlib.h>

#include "diffuse_util.h"
#include "sym_table.h"
#include "mcell_species.h"
#include "logging.h"

/* static helper functions */
static struct species *assemble_mol_species(MCELL_STATE *state,
                                            struct sym_entry *sym_ptr,
                                            struct mcell_species_spec *species);

static int ensure_rdstep_tables_built(MCELL_STATE *state);

/*************************************************************************
 mcell_create_species:
    Create a new species. This uses the same helper functions as the parser,
    but is meant to be used independent of the parser.

 In: state: the simulation state
     name:  molecule name
     D:     diffusion constant
     is_2d: 1 if the species is a 2D molecule, 0 if 3D
     custom_time_step: time_step for the molecule (< 0.0 for a custom space
                       step, >0.0 for custom timestep, 0.0 for default
                       timestep)
     target_only: 1 if the molecule cannot initiate reactions
     max_step_length:
 Out: Returns 0 on sucess and 1 on error
*************************************************************************/
MCELL_STATUS
mcell_create_species(MCELL_STATE *state, struct mcell_species_spec *species,
                     mcell_symbol **species_ptr) {
  struct sym_entry *sym = NULL;
  int error_code = new_mol_species(state, species->name, &sym);
  if (error_code) {
    return error_code;
  }

  assemble_mol_species(state, sym, species);

  error_code = ensure_rdstep_tables_built(state);
  if (error_code) {
    return error_code;
  }

  if (species_ptr != NULL) {
    *species_ptr = sym;
  }

  return MCELL_SUCCESS;
}

/*****************************************************************************
 *
 * mcell_add_to_species_list creates a linked list of mcell_species from
 * mcell_symbols.
 *
 * The list of mcell_species is for example used to provide the list
 * of reactants, products and surface classes needed for creating
 * reactions.
 *
 * During the first invocation of this function, NULL should be provided for
 * the species_list to initialize a new mcell_species list with mcell_symbol.
 * On subsecquent invocations the current mcell_species list should
 * be provided as species_list to which the new mcell_symbol will be appended
 * with the appropriate flags for orientation status.
 *
 *****************************************************************************/
struct mcell_species *
mcell_add_to_species_list(mcell_symbol *species_ptr, bool is_oriented,
                          int orientation, struct mcell_species *species_list) {
  struct mcell_species *species = (struct mcell_species *)CHECKED_MALLOC_STRUCT(
      struct mcell_species, "species list");
  if (species == NULL) {
    return NULL;
  }

  species->next = NULL;
  species->mol_type = species_ptr;
  species->orient_set = is_oriented ? 1 : 0;
  species->orient = orientation;

  if (species_list != NULL) {
    species->next = species_list;
  }

  return species;
}

/*****************************************************************************
 *
 * mcell_delete_species_list frees all memory associated with a list of
 * mcell_species
 *
 *****************************************************************************/
void mcell_delete_species_list(struct mcell_species *species) {
  struct mcell_species *tmp = species;
  while (species) {
    tmp = species->next;
    free(species);
    species = tmp;
  }
}

/**************************************************************************
 new_mol_species:
    Create a new species. There must not yet be a molecule or named reaction
    pathway with the supplied name.

 In: state: the simulation state
     name:  name for the new species
     sym_ptr:   symbol for the species
 Out: 0 on success, positive integer on failure
**************************************************************************/
int new_mol_species(MCELL_STATE *state, char *name, struct sym_entry **sym_ptr) {
  // Molecule already defined
  if (retrieve_sym(name, state->mol_sym_table) != NULL) {
    return 2;
  }
  // Molecule already defined as a named reaction pathway
  else if (retrieve_sym(name, state->rxpn_sym_table) != NULL) {
    return 3;
  }
  *sym_ptr = store_sym(name, MOL, state->mol_sym_table, NULL);
  // Out of memory while creating molecule
  if (*sym_ptr == NULL) {
    return 4;
  }

  return 0;
}

/*****************************************************************************
 *
 * static helper functions
 *
 *****************************************************************************/

/**************************************************************************
assemble_mol_species:
   Helper function to assemble a molecule species from its component pieces.

   NOTE: A couple of comments regarding the unit conversions below:
   Internally, mcell works with with the per species length
   normalization factor

      new_spec->space_step = sqrt(4*D*t), D = diffusion constant (1)

   If the user supplies a CUSTOM_SPACE_STEP or SPACE_STEP then
   it is assumed to correspond to the average diffusion step and
   is hence equivalent to lr_bar in 2 or 3 dimensions for surface and
   volume molecules, respectively:

   lr_bar_2D = sqrt(pi*D*t)       (2)
   lr_bar_3D = 2*sqrt(4*D*t/pi)   (3)

   Hence, given a CUSTOM_SPACE_STEP/SPACE_STEP we need to
   solve eqs (2) and (3) for t and obtain new_spec->space_step
   via equation (1)

   2D:
    lr_bar_2D = sqrt(pi*D*t) => t = (lr_bar_2D^2)/(pi*D)

   3D:
    lr_bar_3D = 2*sqrt(4*D*t/pi) => t = pi*(lr_bar_3D^2)/(16*D)

   The remaining coefficients are:

    - 1.0e8 : needed to convert D from cm^2/s to um^2/s
    - global_time_unit, length_unit, r_length_unit: mcell
      internal time/length conversions.

In: state: the simulation state
    sym_ptr:   symbol for the species
    D:     diffusion constant
    is_2d: 1 if the species is a 2D molecule, 0 if 3D
    custom_time_step: time_step for the molecule (<0.0 for a custom space
                      step, >0.0 for custom timestep, 0.0 for default
                      timestep)
    target_only: 1 if the molecule cannot initiate reactions
Out: the species, or NULL if an error occurred
**************************************************************************/
struct species *assemble_mol_species(MCELL_STATE *state,
                                     struct sym_entry *sym_ptr,
                                     struct mcell_species_spec *species) {
  // Fill in species info

  // The global time step must be defined before creating any species since it
  // is used in calculations involving custom time and space steps
  double global_time_unit = state->time_unit;
  struct species *new_spec = (struct species *)sym_ptr->value;

  if (species->is_2d) {
    new_spec->flags |= ON_GRID;
  } else {
    new_spec->flags &= ~ON_GRID;
  }

  new_spec->D = species->D;
  new_spec->time_step = species->custom_time_step;

  if (species->target_only) {
    new_spec->flags |= CANT_INITIATE;
  }
  if (species->max_step_length > 0) {
    new_spec->flags |= SET_MAX_STEP_LENGTH;
  }

  // Determine the actual space step and time step

  // Immobile (boring)
  if (!distinguishable(new_spec->D, 0, EPS_C)) {
    new_spec->space_step = 0.0;
    new_spec->time_step = 1.0;
  }
  // Custom timestep or spacestep
  else if (new_spec->time_step != 0.0) {
    // Hack--negative value means custom space step
    if (new_spec->time_step < 0) {
      double lr_bar = -new_spec->time_step;
      if (species->is_2d) {
        new_spec->time_step =
            lr_bar * lr_bar / (MY_PI * 1.0e8 * new_spec->D * global_time_unit);
      } else {
        new_spec->time_step =
            lr_bar * lr_bar * MY_PI /
            (16.0 * 1.0e8 * new_spec->D * global_time_unit);
      }
      new_spec->space_step =
          sqrt(4.0 * 1.0e8 * new_spec->D * new_spec->time_step *
               global_time_unit) *
          state->r_length_unit;
    }
    else {
      new_spec->space_step =
          sqrt(4.0 * 1.0e8 * new_spec->D * new_spec->time_step) *
          state->r_length_unit;
      new_spec->time_step /= global_time_unit;
    }
  }
  // Global timestep (this is the typical case)
  else if (!distinguishable(state->space_step, 0, EPS_C)) { 
    new_spec->space_step =
        sqrt(4.0 * 1.0e8 * new_spec->D * global_time_unit) *
        state->r_length_unit;
    new_spec->time_step = 1.0;
  }
  // Global spacestep
  else { 
    double space_step = state->space_step * state->length_unit;
    if (species->is_2d) {
      new_spec->time_step =
          space_step * space_step /
          (MY_PI * 1.0e8 * new_spec->D * global_time_unit);
    }
    else {
      new_spec->time_step =
          space_step * space_step * MY_PI /
          (16.0 * 1.0e8 * new_spec->D * global_time_unit);
    }
    new_spec->space_step = sqrt(4.0 * 1.0e8 * new_spec->D *
                                new_spec->time_step * global_time_unit) *
                                state->r_length_unit;
  }

  new_spec->refl_mols = NULL;
  new_spec->transp_mols = NULL;
  new_spec->absorb_mols = NULL;
  new_spec->clamp_conc_mols = NULL;

  species->custom_time_step = new_spec->time_step;
  species->space_step = new_spec->space_step;

  return new_spec;
}

/**************************************************************************
 ensure_rdstep_tables_built:
    Build the r_step/d_step tables if they haven't been built yet.

 In: state: the simulation state
 Out: 0 on success, positive integer on failure
**************************************************************************/
int ensure_rdstep_tables_built(MCELL_STATE *state) {
  if (state->r_step != NULL && state->r_step_surface != NULL &&
      state->d_step != NULL) {
    return 0;
  }

  if (state->r_step == NULL) {
    // Out of memory while creating r_step data for molecule
    if ((state->r_step = init_r_step(state->radial_subdivisions)) == NULL) {
      return 5;
    }
  }

  if (state->r_step_surface == NULL) {
    state->r_step_surface = init_r_step_surface(state->radial_subdivisions);
    // Cannot store r_step_surface data.
    if (state->r_step_surface == NULL) {
      return 6;
    }
  }

  if (state->d_step == NULL) {
    // Out of memory while creating d_step data for molecule
    if ((state->d_step = init_d_step(state->radial_directions,
                                     &state->num_directions)) == NULL) {
      return 7;
    }

    // Num directions, rounded up to the nearest 2^n - 1
    state->directions_mask = state->num_directions;
    state->directions_mask |= (state->directions_mask >> 1);
    state->directions_mask |= (state->directions_mask >> 2);
    state->directions_mask |= (state->directions_mask >> 4);
    state->directions_mask |= (state->directions_mask >> 8);
    state->directions_mask |= (state->directions_mask >> 16);
    // Internal error: bad number of default RADIAL_DIRECTIONS (max 131072).
    if (state->directions_mask > (1 << 18)) {
      return 8;
    }
  }

  return 0;
}
