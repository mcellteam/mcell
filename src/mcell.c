/******************************************************************************
 *
 * Copyright (C) 2006-2016 by
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

#include <stdlib.h>

#include "mcell_init.h"
#include "mcell_misc.h"
#include "mcell_run.h"
#include "logging.h"
#include "unistd.h"

//#include "api_test.h"

#define CHECKED_CALL_EXIT(function, error_message)                             \
  {                                                                            \
    if (function) {                                                            \
      mcell_print(error_message);                                              \
      exit(1);                                                                 \
    }                                                                          \
  }

int main(int argc, char **argv) {
  u_int procnum = 0;

  // Print the command line parameters.
  mcell_log ( "" );
  int arg_index;
  char cwd_buf[10001];
  mcell_log ( "MCell called from \"%s\" with:", getcwd(cwd_buf,10000) );
  for (arg_index=0; arg_index<argc; arg_index++) {
    mcell_log ( "  Arg %d = %s", arg_index, argv[arg_index] );
  }
  mcell_log ( "" );

  // initialize the mcell simulation
  MCELL_STATE *state = mcell_create();
  CHECKED_CALL_EXIT(!state, "Failed to initialize MCell simulation.");

  // Parse the command line arguments and print out errors if necessary.
  if (mcell_argparse(argc, argv, state)) {
    if (procnum == 0) {
      mcell_print_version();
      mcell_print_usage(argv[0]);
    }
    exit(1);
  }

  CHECKED_CALL_EXIT(
      mcell_init_state(state),
      "An error occured during set up of the initial simulation state");

  if (state->notify->progress_report != NOTIFY_NONE) {
    mcell_print_version();
  }

  // test_api(state);

  // Comment out MDL parsing when testing the API
  CHECKED_CALL_EXIT(mcell_parse_mdl(state),
                    "An error occured during parsing of the mdl file.");

  CHECKED_CALL_EXIT(mcell_init_simulation(state),
                    "An error occured during simulation creation.");

  CHECKED_CALL_EXIT(
      mcell_init_read_checkpoint(state),
      "An error occured during initialization and reading of checkpoint.");

  CHECKED_CALL_EXIT(mcell_init_output(state),
                    "An error occured during setting up of output.");

  CHECKED_CALL_EXIT(mcell_run_simulation(state),
                    "Error running mcell simulation.");

  if (state->notify->progress_report != NOTIFY_NONE) {
    mcell_print("Done running.");
  }

  mcell_print_stats();

  exit(0);
}
