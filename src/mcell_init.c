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
#include <stdlib.h>
#include <string.h>

#include "config.h"
#include "mcell_structs.h"
#include "chkpt.h"
#include "count_util.h"
#include "init.h"
#include "logging.h"
#include "sym_table.h"
#include "mcell_init.h"
#include "mcell_misc.h"
#include "mcell_reactions.h"
#include "dyngeom.h"
#include "chkpt.h"

/* simple wrapper for executing the supplied function call. In case
 * of an error returns with MCELL_FAIL and prints out error_message */
#define CHECKED_CALL(function, error_message)                                  \
  {                                                                            \
    if (function) {                                                            \
      mcell_log(error_message);                                                \
      return MCELL_FAIL;                                                       \
    }                                                                          \
  }

// static helper functions
static int install_usr_signal_handlers(void);

static bool has_micro_rev_and_trimol_rxns(struct species **species_list,
  int n_species, byte has_vol_rev, byte has_surf_rev);

/************************************************************************
 *
 * change the seed
 *
 ************************************************************************/

void mcell_set_seed(MCELL_STATE *state, int seed) {
  u_int signed_seed = (u_int) seed;
  state->seed_seq = signed_seed;
}

/************************************************************************
 *
 * function for initializing the main mcell simulator. MCELL_STATE
 * keeps track of the state of the simulation.
 *
 * Returns NULL on error and a pointer to MCELL_STATE otherwise
 *
 ************************************************************************/
MCELL_STATE *mcell_create() {
  // signal handlers
  if (install_usr_signal_handlers()) {
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

  state->procnum = 0;
  state->rx_hashsize = 0;
  state->iterations = INT_MIN; /* indicates iterations not set */
  state->chkpt_infile = NULL;
  state->chkpt_outfile = NULL;
  state->chkpt_init = 1;
  state->log_freq =
      ULONG_MAX; /* Indicates that this value has not been set by user */
  state->seed_seq = 1;
  state->with_checks_flag = 1;

  time_t begin_time_of_day;
  time(&begin_time_of_day);
  state->begin_timestamp = begin_time_of_day;
  state->initialization_state = "initializing";

  if (!(state->var_sym_table = init_symtab(1024))) {
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
mcell_init_state(MCELL_STATE *state) {
  CHECKED_CALL(
      init_notifications(state),
      "Unknown error while initializing user-notification data structures.");

  CHECKED_CALL(init_variables(state),
               "Unknown error while initializing system variables.");

  CHECKED_CALL(init_data_structures(state),
               "Unknown error while initializing system data structures.");

  return MCELL_SUCCESS;
}

/************************************************************************
 *
 * function for setting up all the internal data structure to get the
 * simulation into a runnable state.
 *
 * NOTE: Before this function can be called the engine user code
 *       either needs to call
 *       - parse_input() to parse a valid MDL file or
 *       - the individiual API functions for adding model elements
 *         (molecules, geometry, ...)
 *
 * Returns 0 on sucess and 1 on error
 *
 ************************************************************************/
MCELL_STATUS
mcell_init_simulation(MCELL_STATE *state) {
  CHECKED_CALL(init_reactions(state), "Error initializing reactions.");

  CHECKED_CALL(init_species(state), "Error initializing species.");

  if (has_micro_rev_and_trimol_rxns(state->species_list, state->n_species,
    state->volume_reversibility, state->surface_reversibility)) {
    mcell_error("Tri-molecular reactions can not be combined with microscopic "
      "reversibility turned on. Please set MICROSCOPIC_REVESIBILITY = NO");
  }

  if (state->notify->progress_report != NOTIFY_NONE)
    mcell_log("Creating geometry (this may take some time)");

  CHECKED_CALL(init_bounding_box(state), "Error initializing bounding box.");
  CHECKED_CALL(init_partitions(state), "Error initializing partitions.");
  CHECKED_CALL(init_vertices_walls(state),
               "Error initializing vertices and walls.");
  CHECKED_CALL(init_regions(state), "Error initializing regions.");

  if (state->place_waypoints_flag) {
    CHECKED_CALL(place_waypoints(state), "Error while placing waypoints.");
  }

  if (state->with_checks_flag) {
    CHECKED_CALL(check_for_overlapped_walls(
        state->rng, state->n_subvols, state->subvol),
        "Error while checking for overlapped walls.");
  }

  CHECKED_CALL(init_surf_mols(state),
               "Error while placing surface molecules on regions.");

  CHECKED_CALL(init_releases(state->releaser), "Error while initializing release sites.");

  CHECKED_CALL(init_species_mesh_transp(state),
               "Error while initializing species-mesh transparency list.");

  CHECKED_CALL(init_counter_name_hash(
      &state->counter_by_name, state->output_block_head),
      "Error while initializing counter name hash.");
  
  /*CHECKED_CALL(init_dynamic_geometry(state),*/
  /*             "Error while initializing scheduled changes in geometry.");*/

  return MCELL_SUCCESS;
}

/************************************************************************
 *
 * Function for recreating the geometry when using dynamic meshes.
 *
 * NOTE: This entails destroying the existing geometry (and a number of related
 *       dependencies), reparsing the appropriate MDLs, and re-initializing the
 *       partitions, geometry, etc.
 *
 * Returns 1 on error and 0 on success
 *
 * NOTE: This is doing a little more than just destroying and reinitializing
 * geometry, so it probably makes sense to rename this and/or split it up into
 * multiple functions.
 *
 ************************************************************************/
MCELL_STATUS
mcell_redo_geom(MCELL_STATE *state) {
  // We set this mainly to take care of some issues with counting, triggers,
  // memory cleanup.
  state->dynamic_geometry_flag = 1;

  CHECKED_CALL(reset_current_counts(
    state->mol_sym_table,
    state->count_hashmask,
    state->count_hash),
    "Error when reseting counters.");

  struct vector3 llf;
  struct vector3 urb;
  if (state->periodic_box_obj) {
    struct polygon_object* p = (struct polygon_object*)(state->periodic_box_obj->contents);
    struct subdivided_box* sb = p->sb;
    llf.x = sb->x[0];
    llf.y = sb->y[0];
    llf.z = sb->z[0];
       
    urb.x = sb->x[1];
    urb.y = sb->y[1];
    urb.z = sb->z[1];
  }

  CHECKED_CALL(destroy_everything(state), "Error when freeing memory.");
  // We need to reenable the ability to parse geometry
  state->disable_polygon_objects = 0;
  // Reparse the geometry and instantiations. Nothing else should be included
  // in these other MDLs.
#ifdef NOSWIG
  CHECKED_CALL(parse_input(state),
               "An error occured during parsing of the mdl file.");
#endif
  CHECKED_CALL(init_bounding_box(state), "Error initializing bounding box.");
  // This should ideally be in destroy_everything
  free(state->subvol);
  if (state->periodic_box_obj) {
    mcell_create_periodic_box(state, "PERIODIC_BOX_INST", &llf, &urb);
  }
  CHECKED_CALL(init_partitions(state), "Error initializing partitions.");
  CHECKED_CALL(init_vertices_walls(state),
               "Error initializing vertices and walls.");
  CHECKED_CALL(init_regions(state), "Error initializing regions.");


  if (state->place_waypoints_flag) {
    CHECKED_CALL(place_waypoints(state), "Error while placing waypoints.");
  }

  if (state->with_checks_flag) {
    CHECKED_CALL(check_for_overlapped_walls(
        state->rng, state->n_subvols, state->subvol),
        "Error while checking for overlapped walls.");
  }
  CHECKED_CALL(init_species_mesh_transp(state),
               "Error while initializing species-mesh transparency list.");
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
mcell_init_read_checkpoint(MCELL_STATE *state) {

  // set up global state in chkpt.c. This is needed to provided
  // the state for the signal triggered checkpointing
  CHECKED_CALL(set_checkpoint_state(state),
    "An error occured during setting the state of the checkpointing routine.");

  if (state->chkpt_infile) { //state->chkpt_flag == 1) {
    CHECKED_CALL(load_checkpoint(state),
      "Error while loading previous checkpoint.");

    long long exec_iterations;
    CHECKED_CALL(init_checkpoint_state(state, &exec_iterations),
      "Error while initializing checkpoint.");

    CHECKED_CALL(reschedule_release_events(state),
      "Error while rescheduling release events");

    /* XXX This is a hack to be backward compatible with the previous
     * MCell behaviour. Basically, as soon as exec_iterations <= 0
     * MCell will stop and we emulate this by returning 1 even though
     * this is not an error (as implied by returning 1). */
    if (exec_iterations <= 0) {
      mem_dump_stats(mcell_get_log_file());
      return MCELL_FAIL;
    }
  } else {
    state->chkpt_seq_num = 1;
  }

  if (state->chkpt_infile) {
  }

  // set the iteration time to the start time of the checkpoint
  state->current_iterations = state->start_iterations;

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
mcell_init_output(MCELL_STATE *state) {
  CHECKED_CALL(init_viz_data(state), "Error while initializing viz data.");
  CHECKED_CALL(init_reaction_data(state),
               "Error while initializing reaction data.");
  CHECKED_CALL(init_timers(state), "Error initializing the simulation timers.");

  // signal successful end of simulation
  state->initialization_state = NULL;

  return MCELL_SUCCESS;
}

/*************************************************************************
 mcell_set_partition:
    Set the partitioning in a particular dimension.

 In:  state: the simulation state
      dim: the dimension whose partitions we'll set
      head: the partitioning
 Out: 0 on success, 1 on failure
*************************************************************************/
MCELL_STATUS mcell_set_partition(MCELL_STATE *state, int dim,
                                 struct num_expr_list_head *head) {
  /* Allocate array for partitions */
  double *dblp = CHECKED_MALLOC_ARRAY(double, (head->value_count + 2),
                                      "volume partitions");
  if (dblp == NULL)
    return MCELL_FAIL;

  /* Copy partitions in sorted order to the array */
  unsigned int num_values = 0;
  dblp[num_values++] = -GIGANTIC;
  struct num_expr_list *nel;
  for (nel = head->value_head; nel != NULL; nel = nel->next)
    dblp[num_values++] = nel->value * state->r_length_unit;
  dblp[num_values++] = GIGANTIC;
  qsort(dblp, num_values, sizeof(double), &double_cmp);

  /* Copy the partitions into the model */
  switch (dim) {
  case X_PARTS:
    if (state->x_partitions != NULL)
      free(state->x_partitions);
    state->nx_parts = num_values;
    state->x_partitions = dblp;
    break;

  case Y_PARTS:
    if (state->y_partitions != NULL)
      free(state->y_partitions);
    state->ny_parts = num_values;
    state->y_partitions = dblp;
    break;

  case Z_PARTS:
    if (state->z_partitions != NULL)
      free(state->z_partitions);
    state->nz_parts = num_values;
    state->z_partitions = dblp;
    break;

  default:
    UNHANDLED_CASE(dim);
  }

  if (!head->shared)
    mcell_free_numeric_list(head->value_head);

  return MCELL_SUCCESS;
}

/*************************************************************************
 mcell_set_iterations:
    Set the number of iterations for the simulation.

 In: state: the simulation state
     iterations: number of iterations to run
 Out: 0 on success; 1 on failure.
      number of iterations is set.
*************************************************************************/
MCELL_STATUS
mcell_set_iterations(MCELL_STATE *state, long long iterations) {
  if (iterations < 0) {
    return MCELL_FAIL;
  }
  state->iterations = iterations;
  return MCELL_SUCCESS;
}

/*************************************************************************
 mcell_set_time_step:
    Set the global timestep for the simulation.

 In: state: the simulation state
      step: timestep to set
 Out: 0 on success; any other integer value is a failure.
      global timestep is updated.
*************************************************************************/
MCELL_STATUS
mcell_set_time_step(MCELL_STATE *state, double step) {
  if (step <= 0) {
    return 2;
  }
  // Timestep was already set. Could introduce subtle problems if we let it
  // change after defining the species, since it is used in calculations there.
  if (distinguishable(state->time_unit, 0, EPS_C)) {
    return 3;
  }
  state->time_unit = step;
  return MCELL_SUCCESS;
}

/*************************************************************************
 mcell_silence_notifications:

    Silence notifications

 In: state: the simulation state
 Out: 0 on success; any other integer value is a failure.
*************************************************************************/
MCELL_STATUS
mcell_silence_notifications(MCELL_STATE *state) {
  /*state->quiet_flag = 1;*/
  state->notify->progress_report = NOTIFY_NONE;
  state->notify->diffusion_constants = NOTIFY_NONE;
  state->notify->reaction_probabilities = NOTIFY_NONE;
  state->notify->time_varying_reactions = NOTIFY_NONE;
  state->notify->reaction_prob_notify = 0.0;
  state->notify->partition_location = NOTIFY_NONE;
  state->notify->box_triangulation = NOTIFY_NONE;
  state->notify->iteration_report = NOTIFY_NONE;
  state->notify->custom_iteration_value = 0;
  state->notify->release_events = NOTIFY_NONE;
  state->notify->file_writes = NOTIFY_NONE;
  state->notify->final_summary = NOTIFY_NONE;
  state->notify->throughput_report = NOTIFY_NONE;
  state->notify->checkpoint_report = NOTIFY_NONE;
  state->notify->reaction_output_report = NOTIFY_NONE;
  state->notify->volume_output_report = NOTIFY_NONE;
  state->notify->viz_output_report = NOTIFY_NONE;
  state->notify->molecule_collision_report = NOTIFY_NONE;
  return MCELL_SUCCESS;
}

/*************************************************************************
 mcell_enable_notifications:

    Enable notifications

 In: state: the simulation state
 Out: 0 on success; any other integer value is a failure.
*************************************************************************/
MCELL_STATUS
mcell_enable_notifications(MCELL_STATE *state) {
  state->notify->progress_report = NOTIFY_FULL;
  state->notify->diffusion_constants = NOTIFY_BRIEF;
  state->notify->reaction_probabilities = NOTIFY_FULL;
  state->notify->time_varying_reactions = NOTIFY_FULL;
  state->notify->reaction_prob_notify = 0.0;
  state->notify->partition_location = NOTIFY_NONE;
  state->notify->box_triangulation = NOTIFY_NONE;
  state->notify->iteration_report = NOTIFY_FULL;
  state->notify->custom_iteration_value = 0;
  state->notify->release_events = NOTIFY_FULL;
  state->notify->file_writes = NOTIFY_NONE;
  state->notify->final_summary = NOTIFY_FULL;
  state->notify->throughput_report = NOTIFY_FULL;
  state->notify->checkpoint_report = NOTIFY_FULL;
  state->notify->reaction_output_report = NOTIFY_NONE;
  state->notify->volume_output_report = NOTIFY_NONE;
  state->notify->viz_output_report = NOTIFY_NONE;
  state->notify->molecule_collision_report = NOTIFY_NONE;
  return MCELL_SUCCESS;
}


/*************************************************************************
 mcell_enable_warnings:

    Enable warnings

 In: state: the simulation state
 Out: 0 on success; any other integer value is a failure.
*************************************************************************/
MCELL_STATUS
mcell_enable_warnings(MCELL_STATE *state) {
  state->notify->neg_diffusion = WARN_WARN;
  state->notify->neg_reaction = WARN_WARN;
  state->notify->high_reaction_prob = WARN_COPE;
  state->notify->reaction_prob_warn = 1.0;
  state->notify->close_partitions = WARN_WARN;
  state->notify->degenerate_polys = WARN_WARN;
  state->notify->overwritten_file = WARN_COPE;
  state->notify->short_lifetime = WARN_WARN;
  state->notify->short_lifetime_value = 50;
  state->notify->missed_reactions = WARN_WARN;
  state->notify->missed_reaction_value = 0.001;
  state->notify->missed_surf_orient = WARN_ERROR;
  state->notify->useless_vol_orient = WARN_WARN;
  state->notify->mol_placement_failure = WARN_WARN;
  state->notify->invalid_output_step_time = WARN_WARN;
  state->notify->large_molecular_displacement = WARN_WARN;
  state->notify->add_remove_mesh_warning = WARN_WARN;
  return MCELL_SUCCESS;
}

/*************************************************************************
 mcell_silence_warnings:

    Silence warnings

 In: state: the simulation state
 Out: 0 on success; any other integer value is a failure.
*************************************************************************/
MCELL_STATUS
mcell_silence_warnings(MCELL_STATE *state) {
  state->notify->neg_diffusion = WARN_COPE;
  state->notify->neg_reaction = WARN_COPE;
  state->notify->high_reaction_prob = WARN_COPE;
  state->notify->close_partitions = WARN_COPE;
  state->notify->degenerate_polys = WARN_COPE;
  state->notify->overwritten_file = WARN_COPE;
  state->notify->short_lifetime = WARN_COPE;
  state->notify->missed_reactions = WARN_COPE;
  state->notify->missed_surf_orient = WARN_ERROR;
  state->notify->useless_vol_orient = WARN_COPE;
  state->notify->mol_placement_failure = WARN_COPE;
  state->notify->invalid_output_step_time = WARN_COPE;
  state->notify->large_molecular_displacement = WARN_COPE;
  state->notify->add_remove_mesh_warning = WARN_COPE;
  return MCELL_SUCCESS;
}


/*****************************************************************************
 *
 * static helper functions
 *
 *****************************************************************************/

/***********************************************************************
 * install_usr_signal_handlers:
 *
 *   Set signal handlers for checkpointing on SIGUSR signals.
 *
 *   In:  None
 *   Out: 0 on success, 1 on failure.
 ***********************************************************************/
int install_usr_signal_handlers(void) {
#ifndef _WIN32 /* fixme: Windows does not support USR signals */
  struct sigaction sa, saPrev;
  sa.sa_sigaction = NULL;
  sa.sa_handler = &chkpt_signal_handler;
  sa.sa_flags = SA_RESTART;
  sigfillset(&sa.sa_mask);

  if (sigaction(SIGUSR1, &sa, &saPrev) != 0) {
    mcell_error("Failed to install USR1 signal handler.");
    /*return 1;*/
  }
  if (sigaction(SIGUSR2, &sa, &saPrev) != 0) {
    mcell_error("Failed to install USR2 signal handler.");
    /*return 1;*/
  }
#endif

  return 0;
}


/*
 * has_micro_rev_and_trimol_rxns tests if the model has surface or
 * volume reversibility selected and also contains trimolecular reactions.
 * In that case it returns true and false otherwise.
 */
bool has_micro_rev_and_trimol_rxns(struct species **species_list,
  int n_species, byte has_vol_rev, byte has_surf_rev) {

  if (!has_vol_rev && !has_surf_rev) {
    return false;
  }

  for (int i = 0; i < n_species; i++) {
    struct species *sp = species_list[i];

    if (sp->flags & (CAN_VOLVOLVOL | CAN_VOLVOLSURF | CAN_VOLSURFSURF |
      CAN_SURFSURFSURF)) {
      return true;
    }
  }
  return false;
}
