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

#include <chkpt.h>
#include <signal.h>
#include <string.h>
#include <time.h>

#include "argparse.h"
#include "count_util.h"
#include "init.h"
#include "libmcell.h"
#include "logging.h"
#include "mem_util.h"
#include "version_info.h"
#include "sym_table.h"


/* declaration of static functions */
static int install_usr_signal_handlers(void);


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
    return 1;
  }

  
  if (init_variables(state))
  {
    mcell_log("Unknown error while initializing system variables.");
    return 1;
  }


  if (init_data_structures(state))
  {
    mcell_log("Unknown error while initializing system data structures.");
    return 1;
  }

  return 0;
}



/************************************************************************
 * 
 * function for running the mcell simulation engine given MCELL_STATE 
 *
 * Returns 0 on sucess and 1 on error 
 *
 * NOTE: This is currently just a very thin wrapper around run_sim()
 *
 ************************************************************************/
/*
MCELL_STATUS 
mcell_run_simulation(MCELL_STATE* state)
{
  run_sim(state);
  return 0;
}
*/


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
    return 1;
  }


  if (state->notify->progress_report != NOTIFY_NONE)
    mcell_log("Creating geometry (this may take some time)");
  if (init_geom(state)) 
  {
    mcell_error_nodie("Error initializing geometry.");
    return 1;
  }

  
  if (init_partitions(state))
  {
    mcell_error_nodie("Error initializing partitions.");
    return 1;
  }


  if (init_vertices_walls(state))
  {
    mcell_error_nodie("Error initializing vertices and walls.");
    return 1;
  }


  if (init_regions(state))
  {
    mcell_error_nodie("Error initializing regions.");
    return 1;
  }


  if (state->place_waypoints_flag)
  {
    if (place_waypoints(state))
    {
      mcell_error_nodie("Error while placing waypoints.");
      return 1;
    }
  }


  if (state->with_checks_flag)
  {
    if(check_for_overlapped_walls(state->n_subvols, state->subvol))
    {
      mcell_error_nodie("Error while checking for overlapped walls.");
      return 1;
    }
  }


  if (init_effectors(state))
  {
    mcell_error_nodie("Error while placing effectors on regions.");
    return 1;
  }


  if (init_releases(state))
  {
    mcell_error_nodie("Error while initializing release sites.");
    return 1;
  }

  return 0;
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
      return 1;
    }

    long long exec_iterations;
    if (init_checkpoint_state(state, &exec_iterations))
    {
      mcell_error_nodie("Error while initializing checkpoint.");
      return 1;
    }

    /* XXX This is a hack to be backward compatible with the previous
     * MCell behaviour. Basically, as soon as exec_iterations <= 0 
     * MCell will stop and we emulate this by returning 1 even though
     * this is not an error (as implied by returning 1). */
    if (exec_iterations <= 0) 
    {
      mem_dump_stats(mcell_get_log_file());
      return 1;
    }
  }
  else 
  {
    state->chkpt_seq_num=1;
  }

  // set the iteration time to the start time of the checkpoint 
  state->it_time = state->start_time;

  return 0;
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
    return 1;
  }


  if (init_reaction_data(state)) 
  {
    mcell_error_nodie("Error while initializing reaction data.");
    return 1;
  }


  if (init_timers(state)) 
  {
    mcell_error_nodie("Error initializing the simulation timers.");
    return 1;
  }

  // signal successful end of simulation
  state->initialization_state = NULL;

  return 0;
}



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
