/******************************************************************************
 *
 * Copyright (C) 2006-2017 by
 * The Salk Institute for Biological Studies and
 * Pittsburgh Supercomputing Center, Carnegie Mellon University
 *
 * Use of this source code is governed by an MIT-style
 * license that can be found in the LICENSE file or at
 * https://opensource.org/licenses/MIT.
 *
******************************************************************************/

#include "config.h"

#include <stdlib.h>

#include "../src4/mcell4_iface_for_mcell3.h"
#include "mcell_init.h"
#include "mcell_misc.h"
#include "mcell_run.h"
#include "init.h"

#include "dump_state.h"

#define CHECKED_CALL_EXIT(function, error_message)                             \
  {                                                                            \
    if (function) {                                                            \
      mcell_print(error_message);                                              \
      exit(1);                                                                 \
    }                                                                          \
  }

int main(int argc, char **argv) {
  u_int procnum = 0;

  // initialize the mcell simulation
  MCELL_STATE *state = mcell_create();
  CHECKED_CALL_EXIT(!state, "Failed to initialize MCell simulation.");

  // Parse the command line arguments and print out errors if necessary.
  if (mcell_argparse(argc, argv, state)) {
    if (procnum == 0) {
      mcell_print("\n\n***************\nArgument error.\n***************\n\n");
      mcell_print_version();
      mcell_print_usage(argv[0]);
    }
    exit(1);
  }

  // Somehow these variables are correct here, but become 0 in mcell_init_state, so save and restore.
  long saved_dump_level = state->dump_level;
  long saved_viz_options = state->viz_options;
  double saved_bond_angle = state->bond_angle; // Might as well do the same for the bond angle!!

  CHECKED_CALL_EXIT(
      mcell_init_state(state),
      "An error occured during set up of the initial simulation state");

  // Restore variables that were set to 0 by mcell_init_state.
  state->dump_level = saved_dump_level;
  state->viz_options = saved_viz_options;
  state->bond_angle = saved_bond_angle;

  if (state->notify->progress_report != NOTIFY_NONE) {
    mcell_print_version();
  }

  CHECKED_CALL_EXIT(parse_input(state),
                    "An error occured during parsing of the mdl file.");

  // full checkpoint read must be done after full initialization,
  // however some values from it are already needed earlier
  CHECKED_CALL_EXIT(
      mcell_init_read_checkpoint_time_and_iteration(state),
      "An error occured during initialization and reading of checkpoint.");

  CHECKED_CALL_EXIT(mcell_init_simulation(state),
                    "An error occured during simulation creation.");

  // read all data from the checkpoint now
  CHECKED_CALL_EXIT(
      mcell_init_read_checkpoint(state),
      "An error occured during initialization and reading of checkpoint.");

  CHECKED_CALL_EXIT(mcell_init_output(state),
                    "An error occured during setting up of output.");

  if (state->dump_mcell3) {
    dump_volume(state, "initial", DUMP_EVERYTHING);
  }
  
  if (state->use_mcell4 || state->mdl2datamodel4) {
    if (!mcell4_convert_mcell3_volume(state)) {
      exit(EXIT_FAILURE);
    }

		if (state->mdl2datamodel4) {
		  mcell4_convert_to_data_model(state->mdl2datamodel4_only_viz);
      mcell4_delete_world();
      return 0;
		}

    mcell4_run_simulation(
        state->dump_mcell4 || state->dump_mcell4_with_geometry, state->dump_mcell4_with_geometry);

    mcell4_delete_world();
  }
  else {
    CHECKED_CALL_EXIT(mcell_run_simulation(state),
                      "Error running mcell simulation.");

    if (state->notify->progress_report != NOTIFY_NONE) {
      mcell_print("Done running.");
    }

    mcell_print_stats();
  }
  exit(0);
}
