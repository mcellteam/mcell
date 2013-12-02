/***********************************************************************************
 *                                                                                 *
 * Copyright (C) 2006-2013 by                                                      *
 * The Salk Institute for Biological Studies and                                   *
 * Pittsburgh Supercomputing Center, Carnegie Mellon University                    *
 *                                                                                 *
 * This program is free software; you can redistribute it and/or                   *
 * modify it under the terms of the GNU General Public License                     *
 * as published by the Free Software Foundation; either version 2                  *
 * of the License, or (at your option) any later version.                          *
 *                                                                                 *
 * This program is distributed in the hope that it will be useful,                 *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of                  *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                   *
 * GNU General Public License for more details.                                    *
 *                                                                                 *
 * You should have received a copy of the GNU General Public License               *
 * along with this program; if not, write to the Free Software                     *
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA. *
 *                                                                                 *
 ***********************************************************************************/

#if defined(__linux__)
#define _GNU_SOURCE 1
#endif

#ifndef _WIN32
#include <sys/resource.h>
#endif
#include <stdlib.h>
#if defined(__linux__)
#include <fenv.h>
#endif

#include <signal.h>
#include <string.h>
#include <time.h>

#include "argparse.h"
#include "chkpt.h"
#include "count_util.h"
#include "init.h"
#include "libmcell.h"
#include "logging.h"
#include "mem_util.h"
#include "react_output.h"
#include "react_util.h"
#include "sym_table.h"
#include "version_info.h"
#include "species.h"


/* declaration of static functions */
static int install_usr_signal_handlers(void);

struct output_column* get_counter_trigger_column(MCELL_STATE* state, 
    const char *counter_name, int column_id);


/************************************************************************
 * 
 * function for initializing the main mcell simulator. MCELL_STATE 
 * keeps track of the state of the simulation.
 *
 * Returns NULL on error and a pointer to MCELL_STATE otherwis
 *
 ************************************************************************/
MCELL_STATE* 
mcell_create() 
{
  // signal handlers
  if (install_usr_signal_handlers()) 
  {
    return NULL;
  }

  // logging
  mcell_set_log_file(stdout);
  mcell_set_error_file(stderr);

  MCELL_STATE *state = CHECKED_MALLOC_STRUCT_NODIE(struct volume, "world");
  if (state == NULL) {
    return NULL;
  }
  memset(state, 0, sizeof(struct volume));

#if defined(__linux__)
  feenableexcept(FE_DIVBYZERO);
#endif

  state->procnum=0;
  state->iterations=INT_MIN; /* indicates iterations not set */
  state->chkpt_infile = NULL;
  state->chkpt_outfile = NULL;
  state->chkpt_init = 1;
  state->log_freq = ULONG_MAX; /* Indicates that this value has not been set by user */
  state->seed_seq = 1;
  state->with_checks_flag = 1;

  time_t begin_time_of_day;
  time(&begin_time_of_day);
  state->begin_timestamp = begin_time_of_day;
  state->initialization_state = "initializing";

  if (!(state->var_sym_table = init_symtab(1024)))
  {
    mcell_log("Failed to initialize MDL variable symbol table.");
    return NULL;
  }

  return state;
}


/************************************************************************
 * 
 * function for initializing the intial simulation state (variables,
 * notifications, data structures)
 *
 * Returns 1 on error and 0 on success 
 *
 ************************************************************************/
MCELL_STATUS
mcell_init_state(MCELL_STATE* state)
{
  if (init_notifications(state))
  {
    mcell_log("Unknown error while initializing user-notification data "
              "structures.");
    return MCELL_FAIL;
  }

  
  if (init_variables(state))
  {
    mcell_log("Unknown error while initializing system variables.");
    return MCELL_FAIL;
  }


  if (init_data_structures(state))
  {
    mcell_log("Unknown error while initializing system data structures.");
    return MCELL_FAIL;
  }

  return MCELL_SUCCESS;
}



/************************************************************************
 * 
 * function for parsing the models underlying mdl file. The function
 * updates the state accordingly.
 *
 * Returns 0 on sucess and 1 on error 
 *
 * NOTE: This is currently just a very thin wrapper around parse_input()
 *
 ************************************************************************/
MCELL_STATUS 
mcell_parse_mdl(MCELL_STATE* state)
{
  return parse_input(state);
}



/************************************************************************
 * 
 * function for setting up all the internal data structure to get the
 * simulation into a runnable state. 
 *
 * NOTE: Before this function can be called the engine user code
 *       either needs to call
 *       - mcell_parse_mdl() to parse a valid MDL file or
 *       - the individiual API functions for adding model elements
 *         (molecules, geometry, ...) 
 *         XXX: These functions don't exist yet!
 *
 * Returns 0 on sucess and 1 on error 
 *
 * NOTE: This is currently just a very thin wrapper around parse_input()
 *
 ************************************************************************/
MCELL_STATUS 
mcell_init_simulation(MCELL_STATE* state) 
{
  if (init_species(state)) 
  {
    mcell_error_nodie("Error initializing species.");
    return MCELL_FAIL;
  }


  if (state->notify->progress_report != NOTIFY_NONE)
    mcell_log("Creating geometry (this may take some time)");
  if (init_geom(state)) 
  {
    mcell_error_nodie("Error initializing geometry.");
    return MCELL_FAIL;
  }

  
  if (init_partitions(state))
  {
    mcell_error_nodie("Error initializing partitions.");
    return MCELL_FAIL;
  }


  if (init_vertices_walls(state))
  {
    mcell_error_nodie("Error initializing vertices and walls.");
    return MCELL_FAIL;
  }


  if (init_regions(state))
  {
    mcell_error_nodie("Error initializing regions.");
    return MCELL_FAIL;
  }


  if (state->place_waypoints_flag)
  {
    if (place_waypoints(state))
    {
      mcell_error_nodie("Error while placing waypoints.");
      return MCELL_FAIL;
    }
  }


  if (state->with_checks_flag)
  {
    if(check_for_overlapped_walls(state->n_subvols, state->subvol))
    {
      mcell_error_nodie("Error while checking for overlapped walls.");
      return MCELL_FAIL;
    }
  }


  if (init_effectors(state))
  {
    mcell_error_nodie("Error while placing effectors on regions.");
    return MCELL_FAIL;
  }


  if (init_releases(state))
  {
    mcell_error_nodie("Error while initializing release sites.");
    return MCELL_FAIL;
  }


  if (init_counter_name_hash(state))
  {
    mcell_error_nodie("Error while initializing counter name hash.");
    return MCELL_FAIL;
  }

  return MCELL_SUCCESS;
}



/************************************************************************
 * 
 * function for reading and initializing the checkpoint if requested
 *
 * Returns 1 on error and 0 on success 
 *
 ************************************************************************/
MCELL_STATUS 
mcell_read_checkpoint(MCELL_STATE* state)
{
  if (state->chkpt_infile)
  {
    if (load_checkpoint(state)) 
    {
      mcell_error_nodie("Error while loading previous checkpoint.");
      return MCELL_FAIL;
    }

    long long exec_iterations;
    if (init_checkpoint_state(state, &exec_iterations))
    {
      mcell_error_nodie("Error while initializing checkpoint.");
      return MCELL_FAIL;
    }

    /* XXX This is a hack to be backward compatible with the previous
     * MCell behaviour. Basically, as soon as exec_iterations <= 0 
     * MCell will stop and we emulate this by returning 1 even though
     * this is not an error (as implied by returning 1). */
    if (exec_iterations <= 0) 
    {
      mem_dump_stats(mcell_get_log_file());
      return MCELL_FAIL;
    }
  }
  else 
  {
    state->chkpt_seq_num=1;
  }

  // set the iteration time to the start time of the checkpoint 
  state->it_time = state->start_time;

  return MCELL_SUCCESS;
}



/************************************************************************
 * 
 * function for initializing the viz and reaction data output
 *
 * XXX: This function has to be called last, i.e. after the 
 *      simulation has been initialized and checkpoint information
 *      been read.
 *
 * Returns 1 on error and 0 on success 
 *
 ************************************************************************/
MCELL_STATUS 
mcell_init_output(MCELL_STATE* state)
{
  if (init_viz_data(state)) 
  {
    mcell_error_nodie("Error while initializing viz data.");
    return MCELL_FAIL;
  }


  if (init_reaction_data(state)) 
  {
    mcell_error_nodie("Error while initializing reaction data.");
    return MCELL_FAIL;
  }


  if (init_timers(state)) 
  {
    mcell_error_nodie("Error initializing the simulation timers.");
    return MCELL_FAIL;
  }

  // signal successful end of simulation
  state->initialization_state = NULL;

  return MCELL_SUCCESS;
}



/************************************************************************
 * 
 * function for retrieving the current value of a given count
 * expression
 *
 * The call expects:
 *
 * - MCELL_STATE
 * - counter_name: a string containing the name of the count statement to 
 *   be retrieved. Currently, the name is identical to the full path to which 
 *   the corresponding reaction output will be written but this may change
 *   in the future
 * - column: int describing the column to be retrieved
 * - count_data: a *double which will receive the actual value
 * - count_data_type: a *count_type_t which will receive the type of the 
 *   data (for casting of count_data)
 *
 * NOTE: This function can be called anytime after the 
 *       REACTION_DATA_OUTPUT has been either parsed or
 *       set up with API calls.
 *
 * Returns 1 on error and 0 on success 
 *
 ************************************************************************/
MCELL_STATUS
mcell_get_counter_value(MCELL_STATE* state, const char *counter_name,
    int column_id, double *count_data, enum count_type_t *count_data_type)
{
  struct output_column *column = NULL;
  if ((column = get_counter_trigger_column(state, counter_name, column_id))
        == NULL)
  {
    return MCELL_FAIL;
  }
     
  // if we happen to encounter trigger data we bail
  if (column->data_type == COUNT_TRIG_STRUCT) 
  {
    return MCELL_FAIL;
  }

  // evaluate the expression and retrieve it
  eval_oexpr_tree(column->expr,1);
  *count_data = (double)column->expr->value;
  *count_data_type = column->data_type;

  return MCELL_SUCCESS;
}



/************************************************************************
 * 
 * function for changing the reaction rate constant of a given named
 * reaction.
 *
 * The call expects:
 *
 * - MCELL_STATE
 * - reaction name: const char* containing the name of reaction
 * - new rate: a double with the new reaction rate constant
 *
 * NOTE: This function can be called anytime after the 
 *       REACTION_DATA_OUTPUT has been either parsed or
 *       set up with API calls.
 *
 * Returns 1 on error and 0 on success 
 *
 ************************************************************************/
MCELL_STATUS
mcell_change_reaction_rate(MCELL_STATE* state, const char *reaction_name,
    double new_rate)
{
  // sanity check
  if (new_rate < 0.0) 
  {
    return MCELL_FAIL;
  }

  // retrive reaction corresponding to name if it exists
  struct rxn *rx = NULL;
  int path_id = 0;
  if (get_rxn_by_name(state->reaction_hash, state->rx_hashsize,
        reaction_name, &rx, &path_id)) {
    return MCELL_FAIL;
  }

  // now change the rate
  if (change_reaction_probability(state, rx, path_id, new_rate))
  {
    return MCELL_FAIL;
  }

  return MCELL_SUCCESS;
}



/**************************************************************************
 *
 * what follows are helper functions *not* part of the actual API.
 *
 * XXX: These functions should absolutely not be called from client
 *      code and will be removed eventually.
 *
 **************************************************************************/


/***********************************************************************
 * install_usr_signal_handlers:
 *
 *   Set signal handlers for checkpointing on SIGUSR signals.
 *
 *   In:  None
 *   Out: 0 on success, 1 on failure.
 ***********************************************************************/
static int 
install_usr_signal_handlers(void)
{
#ifndef _WIN32 /* fixme: Windows does not support USR signals */
  struct sigaction sa, saPrev;
  sa.sa_sigaction = NULL;
  sa.sa_handler = &chkpt_signal_handler;
  sa.sa_flags = SA_RESTART;
  sigfillset(&sa.sa_mask);

  if (sigaction(SIGUSR1, &sa, &saPrev) != 0)
  {
    mcell_error("Failed to install USR1 signal handler.");
    return 1;
  }
  if (sigaction(SIGUSR2, &sa, &saPrev) != 0)
  {
    mcell_error("Failed to install USR2 signal handler.");
    return 1;
  }
#endif

  return 0;
}



/************************************************************************
 * 
 * helper function printing the version string
 *
 ************************************************************************/
void 
mcell_print_version()
{
  print_version(mcell_get_log_file());
}



/************************************************************************
 * 
 * helper function printing the usage information
 *
 ************************************************************************/
void 
mcell_print_usage(const char *executable_name)
{
  print_usage(mcell_get_log_file(), executable_name);
}



/************************************************************************
 * 
 * helper function printing the simulation stats
 *
 ************************************************************************/
void 
mcell_print_stats()
{
  mem_dump_stats(mcell_get_log_file());
}



/************************************************************************
 * 
 * function for printing a string
 *
 * XXX: This is a temporary hack to be able to print in mcell.c
 *      since mcell disables regular printf
 *
 ************************************************************************/
void 
mcell_print(const char *message)
{
  mcell_log("%s", message);
}


/************************************************************************
 * 
 * helper function for parsing the commandline and setting the
 * relevant parts of the state (seed #, logging, ...)
 *
 ************************************************************************/
int
mcell_argparse(int argc, char **argv, MCELL_STATE* state)
{
  return argparse_init(argc, argv, state);
}



/************************************************************************
 * 
 * helper function for retrieving the output_column corresponding
 * to a given count or trigger statement.
 *
 ************************************************************************/
struct output_column*
get_counter_trigger_column(MCELL_STATE* state, const char *counter_name,
    int column_id)
{
  // retrieve the counter for the requested counter_name
  struct sym_table *counter_sym = retrieve_sym(counter_name, 
      state->counter_by_name);
  if (counter_sym == NULL) {
    mcell_log("Failed to retrieve symbol for counter %s.", counter_name);
    return NULL;
  }
  struct output_set *counter = (struct output_set*)(counter_sym->value);
 
  // retrieve the requested column
  struct output_column *column = counter->column_head;
  int count = 0;
  while (count < column_id && column != NULL)
  {
    count++;
    column = column->next;
  }
  if (count != column_id || column == NULL)
  {
    return NULL;
  }

  return column;
}



/*************************************************************************
 mcell_create_species:
    Create a new species. This uses the same helper functions as the parser,
    but is meant to be used independent of the parser.

 In: state: the simulation state
     name:  molecule name
     D:     diffusion constant
     is_2d: 1 if the species is a 2D molecule, 0 if 3D
     custom_time_step: time_step for the molecule (< 0.0 for a custom space
                       step, >0.0 for custom timestep, 0.0 for default
                       timestep)
     target_only: 1 if the molecule cannot initiate reactions
 Out: Returns 0 on sucess and 1 on error 
*************************************************************************/
MCELL_STATUS
mcell_create_species(MCELL_STATE* state,
                     char *name,
                     double D,
                     int is_2d,
                     double custom_time_step,
                     int target_only,
                     double max_step_length)

{
  // Store the new molecule in the symbol table.
  struct sym_table *sym = new_mol_species(state, name);
  if (!sym) {
    return MCELL_FAIL;
  }
  // Perhaps we should consider getting rid of D_ref. It doesn't seem to be
  // used for anything.
  int D_ref = D; 
  struct species* spec = assemble_mol_species(
    state, sym, D_ref, D, is_2d, custom_time_step, target_only,
    max_step_length);
  // Print out information about the diffusion distances
  finish_molecule(state, spec);
  return MCELL_SUCCESS;
}
