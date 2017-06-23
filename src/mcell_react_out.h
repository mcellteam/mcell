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

#pragma once

#include "mcell_init.h"

struct output_column_list {
  struct output_column *column_head;
  struct output_column *column_tail;
};

struct output_set_list {
  struct output_set *set_head;
  struct output_set *set_tail;
};

struct output_times_inlist {
  enum output_timer_type_t type;
  double step;
  struct num_expr_list_head values;
};

int mcell_get_count(char *mol_name, char *reg_name, struct volume *world);

struct output_request *mcell_new_output_request(MCELL_STATE *state,
                                                struct sym_entry *target,
                                                short orientation,
                                                struct sym_entry *location,
                                                struct periodic_image *img,
                                                int report_flags);

struct output_set *mcell_create_new_output_set(char *comment, int exact_time,
                                               struct output_column *col_head,
                                               int file_flags,
                                               char *outfile_name);

MCELL_STATUS mcell_prepare_single_count_expr(struct output_column_list *list,
                                             struct output_expression *expr,
                                             char *custom_header);

MCELL_STATUS
mcell_add_reaction_output_block(MCELL_STATE *state,
                                struct output_set_list *osets, int buffer_size,
                                struct output_times_inlist *otimes);

MCELL_STATUS mcell_create_count(MCELL_STATE *state, struct sym_entry *target,
                                short orientation, struct sym_entry *location,
                                int report_flags, char *custom_header,
                                struct output_column_list *count_list);

MCELL_STATUS mcell_get_counter_value(MCELL_STATE *state,
                                     const char *counter_name, int column,
                                     double *count_data,
                                     enum count_type_t *count_data_type);
