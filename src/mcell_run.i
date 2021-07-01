/******************************************************************************
 *
 * Copyright (C) 2006-2015 by
 * The Salk Institute for Biological Studies and
 * Pittsburgh Supercomputing Center, Carnegie Mellon University
 *
 * Use of this source code is governed by an MIT-style
 * license that can be found in the LICENSE file or at
 * https://opensource.org/licenses/MIT.
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


