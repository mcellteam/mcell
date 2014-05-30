/******************************************************************************
 *
 * Copyright (C) 2006-2014 by
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

#ifndef CREATE_REACTION_OUTPUT_H
#define CREATE_REACTION_OUTPUT_H

#include "libmcell.h"

struct output_column *new_output_column();

struct output_block *new_output_block(int buffersize);

void set_reaction_output_timer_step(MCELL_STATE *state, struct output_block *obp,
  double step);

int set_reaction_output_timer_iterations(MCELL_STATE *state,
  struct output_block *obp, struct num_expr_list_head *step_values);

int set_reaction_output_timer_times(MCELL_STATE *state, struct output_block *obp,
  struct num_expr_list_head *step_values);

int output_block_finalize(MCELL_STATE *state, struct output_block *obp);

#endif
