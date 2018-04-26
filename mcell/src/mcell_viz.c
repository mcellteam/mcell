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

#include <assert.h>
#include <string.h>
#include <stdlib.h>

#include "sym_table.h"
#include "mcell_species.h"
#include "mcell_viz.h"
#include "mcell_misc.h"
#include "logging.h"

static int select_viz_molecules(struct mcell_species *mol_viz_list,
                                struct viz_output_block *vizblk);

static struct frame_data_list *create_viz_frame(long long iterations,
                                                long long start, long long end,
                                                long long step);

/*************************************************************************
 mcell_create_viz_output:
    Create a new set of viz output.

 In:  state: MCell state
      filename: the path and filename prefix where the viz data will be
                (e.g. "./viz_data/my_viz")
      mol_viz_list: a list of the molecules to be visualized
      start: the first frame of the viz data
      end: the last frame of the viz data
      step: the delta between iterations
 Out: Returns 1 on error and 0 on success
 Note: Right now, only iterations (not time points) can be specified.
*************************************************************************/
MCELL_STATUS
mcell_create_viz_output(MCELL_STATE *state, char *filename,
                        struct mcell_species *mol_viz_list, long long start,
                        long long end, long long step) {

  struct viz_output_block *vizblk = CHECKED_MALLOC_STRUCT(
      struct viz_output_block, "visualization data output parameters");
  if (vizblk == NULL)
    return MCELL_FAIL;

  mcell_new_viz_output_block(vizblk);
  // In principal, it's possible to have multiple viz blocks, but this isn't
  // supported in the API yet.
  vizblk->next = state->viz_blocks;
  state->viz_blocks = vizblk;

  // Only CELLBLENDER mode is supported right now
  vizblk->viz_mode = CELLBLENDER_MODE;

  // Set the viz output path and filename prefix
  vizblk->file_prefix_name = CHECKED_STRDUP(filename, "file_prefix_name");

  // Select which molecules will be visualized
  if (select_viz_molecules(mol_viz_list, vizblk))
    return MCELL_FAIL;

  // Select which iterations will be visualized
  struct frame_data_list *new_frame =
      create_viz_frame(state->iterations, start, end, step);
  if (new_frame == NULL)
    return MCELL_FAIL;

  new_frame->next = NULL;
  state->viz_blocks->frame_data_head = new_frame;

  return MCELL_SUCCESS;
}

/**************************************************************************
 mcell_new_viz_output_block:
    Build a new VIZ output block, containing parameters for an output set for
    visualization.
**************************************************************************/
void mcell_new_viz_output_block(struct viz_output_block *vizblk) {
  vizblk->frame_data_head = NULL;
  vizblk->viz_mode = -1;
  vizblk->file_prefix_name = NULL;
  vizblk->viz_output_flag = 0;
  vizblk->species_viz_states = NULL;

  if (pointer_hash_init(&vizblk->parser_species_viz_states, 32))
    mcell_allocfailed("Failed to initialize viz species states table.");
}

/**************************************************************************
 mcell_mcell_create_viz_frame:
    Create a frame for output in the visualization.

 In: time_type: either OUTPUT_BY_TIME_LIST or OUTPUT_BY_ITERATION_LIST
     type: the type (MOL_POS, etc.)
     iteration_list: list of iterations/times at which to output
 Out: the frame_data_list object, if successful, or NULL if we ran out of
      memory
**************************************************************************/
struct frame_data_list *
mcell_create_viz_frame(int time_type, int type,
                       struct num_expr_list *iteration_list) {

  struct frame_data_list *fdlp;
  fdlp = CHECKED_MALLOC_STRUCT(struct frame_data_list, "VIZ_OUTPUT frame data");
  if (fdlp == NULL)
    return NULL;

  fdlp->list_type = time_type;
  fdlp->type = type;
  fdlp->viz_iteration = -1;
  fdlp->n_viz_iterations = 0;
  fdlp->iteration_list = iteration_list;
  fdlp->curr_viz_iteration = iteration_list;
  return fdlp;
}

/**************************************************************************
 mcell_set_molecule_viz_state:
    Sets a flag on a viz block, requesting that a molecule is visualized.

 In: vizblk: the viz block to check
     specp: the molecule species
     viz_state: the desired viz state
 Out: 0 on success, 1 on failure
**************************************************************************/
MCELL_STATUS
mcell_set_molecule_viz_state(struct viz_output_block *vizblk,
                             struct species *specp, int viz_state) {

  /* Make sure not to override a specific state with a generic state. */
  if (viz_state == INCLUDE_OBJ) {
    void *const exclude = (void *)(intptr_t)EXCLUDE_OBJ;

    void *oldval = pointer_hash_lookup_ext(&vizblk->parser_species_viz_states,
                                           specp, specp->hashval, exclude);
    if (oldval != exclude)
      return 0;
  } else
    vizblk->viz_output_flag |= VIZ_MOLECULES_STATES;

  /* Store new value in the hashtable or die trying. */
  void *val = (void *)(intptr_t)viz_state;
  assert(viz_state == (int)(intptr_t)val);
  if (pointer_hash_add(&vizblk->parser_species_viz_states, specp,
                       specp->hashval, val)) {
    mcell_allocfailed(
        "Failed to store viz state for molecules of species '%s'.",
        specp->sym->name);
    return MCELL_FAIL;
  }
  return MCELL_SUCCESS;
}

/******************************************************************************
 *
 * static helper functions
 *
 ******************************************************************************/
int select_viz_molecules(struct mcell_species *mol_viz_list,
                         struct viz_output_block *vizblk) {

  // Select individual molecules to visualize
  struct mcell_species *current_molecule;
  for (current_molecule = mol_viz_list; current_molecule != NULL;
       current_molecule = current_molecule->next) {
    struct species *specp = (struct species *)current_molecule->mol_type->value;
    if (mcell_set_molecule_viz_state(vizblk, specp, INCLUDE_OBJ))
      return 1;
  }
  return 0;
}

struct frame_data_list *create_viz_frame(long long iterations, long long start,
                                         long long end, long long step) {
  struct num_expr_list_head *list =
      CHECKED_MALLOC_STRUCT(struct num_expr_list_head, "numeric list head");
  if (list == NULL)
    return NULL;
  list->value_head = NULL;
  list->value_tail = NULL;
  list->value_count = 0;
  list->shared = 0;

  // Build a list of iterations for VIZ output
  if (end > iterations)
    end = iterations;
  for (long long current = start; current <= end; current += step) {
    struct num_expr_list *nel =
        CHECKED_MALLOC_STRUCT(struct num_expr_list, "VIZ_OUTPUT iteration");
    if (nel == NULL) {
      mcell_free_numeric_list(list->value_head);
      free(list);
      return NULL;
    }

    ++list->value_count;
    if (list->value_tail)
      list->value_tail = list->value_tail->next = nel;
    else
      list->value_head = list->value_tail = nel;
    list->value_tail->value = current;
    list->value_tail->next = NULL;
  }

  // Sorting is unnecessary at the moment
  // mcell_sort_numeric_list(list->value_head);
  // struct num_expr_list *times_sorted = list->value_head;

  // Create the viz frames using the list of sorted times
  struct frame_data_list *new_frame;
  if ((new_frame = mcell_create_viz_frame(
           OUTPUT_BY_ITERATION_LIST, ALL_MOL_DATA, list->value_head)) == NULL) {
    free(list);
    return NULL;
  }
  free(list);
  return new_frame;
}
