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

#ifndef LIBMCELL_H
#define LIBMCELL_H

#include <stdbool.h>

#include "config.h"

#include "mcell_engine.h"
#include "mcell_init.h"
#include "mcell_structs.h"

/**********************************************************************
 * type declarations
 **********************************************************************/

/****************************************************************
 * routines for running simulations
 ****************************************************************/

/* this function runs the whole simulations */
MCELL_STATUS mcell_run_simulation(MCELL_STATE *state);

/* returns the recommended output frequence either based on
 * a user request in the MDL or via some heuristics */
//long long mcell_determine_output_frequency(MCELL_STATE *state);

/* this function runs a single iteration of simulations */
MCELL_STATUS mcell_run_iteration(MCELL_STATE *state, long long output_frequency,
                                 int *restarted_from_checkpoint);

/* flush all output buffers to disk to disk after the simulation
 * run is complete */
MCELL_STATUS mcell_flush_data(MCELL_STATE *state);

/* print any warnings that were gererated during the simulation
 * run */
MCELL_STATUS mcell_print_final_warnings(MCELL_STATE *state);

/* print the final simulation statistics */
MCELL_STATUS mcell_print_final_statistics(MCELL_STATE *state);

/*****************************************************************
 * non API helper functions
 *
 * NOTE: These functions should *not* be called directly and are
 *       not part of the API. These functions are currently used by
 *       both libmcell and the parser and need to eventually be
 *       internalized by libmcell once the parser is fully API
 *       compliant.
 *****************************************************************/
void mcell_print_version();
void mcell_print_usage(const char *executable_name);
void mcell_print_stats();
int mcell_argparse(int argc, char **argv, MCELL_STATE *state);


/* helper functions for dealing with expression lists - mostly during parsing */
struct num_expr_list * mcell_copysort_numeric_list(struct num_expr_list *head);

void mcell_sort_numeric_list(struct num_expr_list *head);

void mcell_free_numeric_list(struct num_expr_list *nel);

MCELL_STATUS mcell_generate_range(struct num_expr_list_head *list,
                                  double start, double end, double step);

// Maybe move this somewhere else
int advance_range(struct num_expr_list_head *list, double tmp_dbl);

int mcell_generate_range_singleton(struct num_expr_list_head *lh, double value);


/* helper functions for IO */

// XXX this is a temporary hack to be able to print in mcell.c
// since mcell disables regular printf
void mcell_print(const char *message);

#endif
