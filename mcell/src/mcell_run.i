/******************************************************************************
 *
 * Copyright (C) 2006-2015 by
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
/* this function runs the whole simulations */
MCELL_STATUS mcell_run_simulation(MCELL_STATE *state);

/* this function runs n iterations */
MCELL_STATUS mcell_run_n_iterations(MCELL_STATE *state, long long output_frequency,
                    int *INPUT, int n_iter);

/* this function runs a single iteration of simulations */
MCELL_STATUS mcell_run_iteration(MCELL_STATE *state, long long output_frequency,
                                 int *INPUT);

/* flush all output buffers to disk to disk after the simulation
 * run is complete */
MCELL_STATUS mcell_flush_data(MCELL_STATE *state);

/* print any warnings that were gererated during the simulation
 * run */
MCELL_STATUS mcell_print_final_warnings(MCELL_STATE *state);

/* print the final simulation statistics */
MCELL_STATUS mcell_print_final_statistics(MCELL_STATE *state);


