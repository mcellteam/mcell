/***********************************************************************************
 *                                                                                 *
 * Copyright (C) 2006-2014 by *
 * The Salk Institute for Biological Studies and *
 * Pittsburgh Supercomputing Center, Carnegie Mellon University *
 *                                                                                 *
 * This program is free software; you can redistribute it and/or *
 * modify it under the terms of the GNU General Public License *
 * as published by the Free Software Foundation; either version 2 *
 * of the License, or (at your option) any later version. *
 *                                                                                 *
 * This program is distributed in the hope that it will be useful, *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the *
 * GNU General Public License for more details. *
 *                                                                                 *
 * You should have received a copy of the GNU General Public License *
 * along with this program; if not, write to the Free Software *
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301,
 *USA. *
 *                                                                                 *
 ***********************************************************************************/

#include "create_viz_output.h"

#include <stdlib.h>

int select_viz_molecules(struct mcell_species *mol_viz_list,
                         struct viz_output_block *vizblk) {
  // Set flags to visualize all molecules
  //vizblk->viz_output_flag = VIZ_ALL_MOLECULES;
  //vizblk->default_mol_state = INT_MAX;

  // Select individual molecules to visualize
  struct mcell_species *current_molecule;
  for (current_molecule = mol_viz_list;
       current_molecule != NULL;
       current_molecule = current_molecule->next) {
    struct species *specp = (struct species *)current_molecule->mol_type->value;
    if (mcell_set_molecule_viz_state(vizblk, specp, INCLUDE_OBJ))
      return 1;
  }
  return 0;
}
  
struct frame_data_list *create_viz_frame(long long iterations,
                                         long long start,
                                         long long end,
                                         long long step) {
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
  for (long long current = start; current <= end; current+=step) {
    struct num_expr_list *nel = CHECKED_MALLOC_STRUCT(
      struct num_expr_list, "VIZ_OUTPUT iteration");
    if (nel == NULL)
      return NULL;

    ++list->value_count;
    if (list->value_tail)
      list->value_tail = list->value_tail->next = nel;
    else
      list->value_head = list->value_tail = nel;
    list->value_tail->value = current;
    list->value_tail->next = NULL;
  }

  // Sorting is unnecessary at the moment
  //mcell_sort_numeric_list(list->value_head);
  //struct num_expr_list *times_sorted = list->value_head;
  
  // Create the viz frames using the list of sorted times
  struct frame_data_list *new_frame;
  if ((new_frame = mcell_create_viz_frame(
      OUTPUT_BY_ITERATION_LIST, ALL_MOL_DATA, list->value_head)) == NULL) {
    return NULL;
  }
  free(list);
  return new_frame;
}
