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

#include "mcell_structs.h"

/* status of libMCell API calls */
typedef int MCELL_STATUS;

#define MCELL_SUCCESS 0
#define MCELL_FAIL 1

/* state of mcell simulation */
typedef struct volume MCELL_STATE;

void mcell_set_seed(MCELL_STATE *state, int seed);

MCELL_STATE *mcell_create(void);

MCELL_STATUS mcell_init_state(MCELL_STATE *state);

MCELL_STATUS mcell_init_simulation(MCELL_STATE *state);

MCELL_STATUS mcell_redo_geom(MCELL_STATE *state);

MCELL_STATUS mcell_init_read_checkpoint(MCELL_STATE *state);

MCELL_STATUS mcell_init_output(MCELL_STATE *state);

MCELL_STATUS mcell_set_partition(MCELL_STATE *state, int dim,
                                 struct num_expr_list_head *head);

MCELL_STATUS mcell_set_time_step(MCELL_STATE *state, double step);

MCELL_STATUS mcell_set_iterations(MCELL_STATE *state, long long iterations);

MCELL_STATUS mcell_silence_notifications(MCELL_STATE *state);
MCELL_STATUS mcell_enable_notifications(MCELL_STATE *state);
MCELL_STATUS mcell_silence_warnings(MCELL_STATE *state);
MCELL_STATUS mcell_enable_warnings(MCELL_STATE *state);
