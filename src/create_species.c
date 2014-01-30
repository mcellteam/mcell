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

#include "libmcell.h"
#include "logging.h"
#include "sym_table.h"
#include "create_species.h"



/**************************************************************************
 assemble_mol_species:
    Helper function to assemble a molecule species from its component pieces.
   
    NOTE: A couple of comments regarding the unit conversions below:
    Internally, mcell works with with the per species length
    normalization factor
   
       new_species->space_step = sqrt(4*D*t), D = diffusion constant (1)
   
    If the user supplies a CUSTOM_SPACE_STEP or SPACE_STEP then
    it is assumed to correspond to the average diffusion step and
    is hence equivalent to lr_bar in 2 or 3 dimensions for surface and
    volume molecules, respectively:
   
    lr_bar_2D = sqrt(pi*D*t)       (2)
    lr_bar_3D = 2*sqrt(4*D*t/pi)   (3)
   
    Hence, given a CUSTOM_SPACE_STEP/SPACE_STEP we need to
    solve eqs (2) and (3) for t and obtain new_species->space_step
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
     sym:   symbol for the species
     D_ref: reference diffusion constant
     D:     diffusion constant
     is_2d: 1 if the species is a 2D molecule, 0 if 3D
     custom_time_step: time_step for the molecule (<0.0 for a custom space
                       step, >0.0 for custom timestep, 0.0 for default
                       timestep)
     target_only: 1 if the molecule cannot initiate reactions
 Out: the species, or NULL if an error occurred
**************************************************************************/
struct species*
assemble_mol_species(MCELL_STATE* state,
                     struct sym_table *sym,
                     double D_ref,
                     double D,
                     int is_2d,
                     double custom_time_step,
                     int target_only,
                     double max_step_length)
{
  // Fill in species info 
  
  // The global time step must be defined before creating any species since it
  // is used in calculations involving custom time and space steps
  double global_time_unit = state->time_unit;
  struct species *new_species = (struct species *) sym->value;

  if (is_2d) {
    new_species->flags |= ON_GRID;
  }
  else {
    new_species->flags &= ~ON_GRID;
  }

  new_species->D = D;
  new_species->D_ref = D_ref;
  new_species->time_step = custom_time_step;

  if (new_species->D_ref == 0) {
    new_species->D_ref = new_species->D;
  }
  if (target_only) {
    new_species->flags |= CANT_INITIATE;
  }
  if (max_step_length > 0) {
    new_species->flags |= SET_MAX_STEP_LENGTH;
  }

  // Determine the actual space step and time step

  if (new_species->D == 0) /* Immobile (boring) */
  {
    new_species->space_step = 0.0;
    new_species->time_step = 1.0;
  }
  else if (new_species->time_step != 0.0) /* Custom timestep */
  {
    if (new_species->time_step < 0) /* Hack--negative value means custom space step */
    {
      double lr_bar = -new_species->time_step;
      if (is_2d)
      {
        new_species->time_step = 
          lr_bar * lr_bar / (MY_PI * 1.0e8 * new_species->D * global_time_unit);
      }
      else
      {
        new_species->time_step =
          lr_bar * lr_bar * MY_PI / 
          (16.0 * 1.0e8 * new_species->D * global_time_unit);
      }
      new_species->space_step =
        sqrt(4.0 * 1.0e8 * new_species->D * new_species->time_step * global_time_unit)
        * state->r_length_unit;
    }
    else
    {
      new_species->space_step =
        sqrt(4.0 * 1.0e8 * new_species->D * new_species->time_step)
        * state->r_length_unit;
      new_species->time_step /= global_time_unit;
    }
  }
  else if (state->space_step == 0) /* Global timestep */
  {
    new_species->space_step =
      sqrt(4.0 * 1.0e8 * new_species->D * global_time_unit)
      * state->r_length_unit;
    new_species->time_step=1.0;
  }
  else /* Global spacestep */
  {
    double space_step = state->space_step * state->length_unit;
    if (is_2d)
    {
      new_species->time_step = space_step * space_step /
        (MY_PI* 1.0e8 * new_species->D * global_time_unit);
    }
    else
    {
      new_species->time_step =
        space_step * space_step * MY_PI /
        (16.0 * 1.0e8 * new_species->D * global_time_unit);
    }
    new_species->space_step =
      sqrt(4.0 * 1.0e8 * new_species->D * new_species->time_step * global_time_unit)
      * state->r_length_unit;
  }

  new_species->refl_mols = NULL;
  new_species->transp_mols = NULL;
  new_species->absorb_mols = NULL;
  new_species->clamp_conc_mols = NULL;

  return new_species;

}



/*************************************************************************
 add_to_species_list:
    Helper function to add a species to a species list.

 In:  species_list_mem:
      list: the list of species
      spec: the species to be added
 Out: 0 on success, 1 on failure
*************************************************************************/
int
add_to_species_list(struct mem_helper *species_list_mem,
                    struct species_list *list,
                    struct species *spec)
{
  struct species_list_item *ptrl = (struct species_list_item *) \
    CHECKED_MEM_GET(species_list_mem, "species list");
  if (ptrl == NULL) {
    return 1;
  }

  ptrl->spec = spec;
  ptrl->next = NULL;
  if (list->species_count == 0)
  {
    list->species_tail = list->species_head = ptrl;
    list->species_count = 1;
  }
  else
  {
    list->species_tail = list->species_tail->next = ptrl;
    ++ list->species_count;
  }
  return 0;
}



/*************************************************************************
 report_diffusion_distances:
    Helper function to print average diffusion distances per species.

 In:  spec: the species
      time_unit:
      length_unit:
      lvl:
 Out: 0 on success, 1 on failure
*************************************************************************/
static void
report_diffusion_distances(struct species *spec,
                           double time_unit,
                           double length_unit,
                           int lvl)
{

  double l_perp_bar = 0;
  double l_perp_rms = 0;
  double l_r_bar = 0;
  double l_r_rms = 0;

  if (spec->time_step == 1.0)
  {
    /* Theoretical average diffusion distances for the molecule need to
     * distinguish between 2D and 3D molecules for computing l_r_bar and
     * firiends */
    if ((spec->flags & NOT_FREE) == 0)
    {
      l_perp_bar = sqrt(4 * 1.0e8 * spec->D * time_unit / MY_PI);
      l_perp_rms = sqrt(2 * 1.0e8 * spec->D * time_unit);
      l_r_bar = 2 * l_perp_bar;
      l_r_rms = sqrt(6 * 1.0e8 * spec->D * time_unit);
    }
    else {
      l_r_bar = sqrt(MY_PI * 1.0e8 * spec->D * time_unit);
    }


    if (lvl == NOTIFY_FULL)
    {

      mcell_log("MCell: Theoretical average diffusion distances for molecule %s:\n"
                "\tl_r_bar = %.9g microns\n"
                "\tl_r_rms = %.9g microns\n"
                "\tl_perp_bar = %.9g microns\n"
                "\tl_perp_rms = %.9g microns",
                spec->sym->name,
                l_r_bar,
                l_r_rms,
                l_perp_bar,
                l_perp_rms);
    }
    else if (lvl == NOTIFY_BRIEF)
      mcell_log("  l_r_bar=%.9g um for %s", l_r_bar, spec->sym->name);
  }
  else
  {
    if (lvl == NOTIFY_FULL)
    {
      /* the size of the length unit depends on if the molecule is
       * 2D or 3D; the values for step length simply follow from
       * converting space_step = sqrt(4Dt) into the 2D/3D expression
       * for l_r_bar */
      double step_length = 0.0;
      if ((spec->flags & NOT_FREE) == 0) {
        step_length = length_unit * spec->space_step * 2.0 / sqrt(MY_PI);
      }
      else {
        step_length = length_unit * spec->space_step * sqrt(MY_PI) / 2.0;
      }

      mcell_log("MCell: Theoretical average diffusion time for molecule %s:\n"
                "\tl_r_bar fixed at %.9g microns\n"
                "\tPosition update every %.3e seconds (%.3g timesteps)",
                spec->sym->name,
                step_length,
                spec->time_step * time_unit, spec->time_step);
    }
    else if (lvl == NOTIFY_BRIEF)
    {
      mcell_log("  delta t=%.3g timesteps for %s", spec->time_step, spec->sym->name);
    }
  }
}



/*************************************************************************
 print_species_summary:
    Helper function to finish the creation of a single molecule, undoing any
    state changes we made during the creation of the molecule. Presently, this
    just means "print the diffusion distances report".

 In:  state: the simulation state
      mol:   species finished
 Out: A report is printed to the file handle.
*************************************************************************/
void
print_species_summary(MCELL_STATE* state, struct species *mol)
{
  if (state->procnum == 0)
  {
    if (state->notify->diffusion_constants == NOTIFY_BRIEF) {
      mcell_log("Defining molecule with the following diffusion constant:");
    }
    report_diffusion_distances(
      mol, state->time_unit, state->length_unit,
      state->notify->diffusion_constants);
    no_printf("Molecule %s defined with D = %g\n", mol->sym->name, mol->D);
  }
}



/*************************************************************************
 print_species_summaries:
    Helper function to finish the creation of multiple molecules, undoing any
    state changes we made during the creation of the molecule. Presently, this
    just means "print the diffusion distances report".

 In: state: the simulation state
     mols: species finished
 Out: A report is printed to the file handle.
*************************************************************************/
void
print_species_summaries(MCELL_STATE* state,
                        struct species_list_item *mols)
{
  if (state->procnum == 0)
  {
    struct species_list_item *ptrl;
    if (state->notify->diffusion_constants == NOTIFY_BRIEF) {
      mcell_log("Defining molecules with the following theoretical average diffusion distances:");
    }
    for (ptrl = mols; ptrl != NULL; ptrl = ptrl->next)
    {
      struct species *spec = (struct species *) ptrl->spec;
      report_diffusion_distances(
        spec, state->time_unit, state->length_unit,
        state->notify->diffusion_constants);
      no_printf("Molecule %s defined with D = %g\n", spec->sym->name, spec->D);
    }
    if (state->notify->diffusion_constants == NOTIFY_BRIEF) {
      mcell_log_raw("\n");
    }
  }
}
