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

#include <assert.h>
#include <float.h>
#include <math.h>
#include <unistd.h>

#ifndef _WIN32
#include <sys/resource.h>
#endif
#if defined(__linux__)
#include <fenv.h>
#endif

#include "sym_table.h"
#include "logging.h"
#include "vol_util.h"
#include "react_output.h"
#include "viz_output.h"
#include "volume_output.h"
#include "diffuse.h"
#include "init.h"
#include "chkpt.h"
#include "argparse.h"
#include "dyngeom.h"

#include "mcell_run.h"
#include "mcell_reactions.h"
#include "mcell_react_out.h"

// static helper functions
static long long mcell_determine_output_frequency(MCELL_STATE *state);

/***********************************************************************
 process_volume_output:

    Produce this round's volume output, if any.

    In:  struct volume *wrld - the world
         double not_yet - earliest time which should not yet be output
    Out: none.  volume output files are updated as appropriate.
 ***********************************************************************/
static void process_volume_output(struct volume *wrld, double not_yet) {
  struct volume_output_item *vo;
  for (vo = (struct volume_output_item *)schedule_next(
           wrld->volume_output_scheduler);
       vo != NULL || not_yet >= wrld->volume_output_scheduler->now;
       vo = (struct volume_output_item *)schedule_next(
           wrld->volume_output_scheduler)) {
    if (vo == NULL)
      continue;
    if (update_volume_output(wrld, vo))
      mcell_error("Failed to update volume output.");
  }
}

/***********************************************************************
 process_reaction_output:

    Produce this round's reaction output, if any.

 In: wrld: the world
     not_yet: earliest time which should not yet be output
 Out: none.  reaction output data structs/files are updated as appropriate.
 ***********************************************************************/
static void process_reaction_output(struct volume *wrld, double not_yet) {
  struct output_block *obp;
  for (obp = schedule_next(wrld->count_scheduler);
       obp != NULL || not_yet >= wrld->count_scheduler->now;
       obp = schedule_next(wrld->count_scheduler)) {
    if (obp == NULL)
      continue;
    if (update_reaction_output(wrld, obp))
      mcell_error("Failed to update reaction output.");
  }
  if (wrld->count_scheduler->error)
    mcell_internal_error("Scheduler reported an out-of-memory error while "
                         "retrieving next scheduled reaction output, but this "
                         "should never happen.");
}

/***********************************************************************
 process_molecule_releases:

    Produce this round's release events, if any.

 In: wrld: the world
     not_yet: earliest time which should not yet be processed
 Out: none.  molecules are released into the world.
 ***********************************************************************/
void process_molecule_releases(struct volume *wrld, double not_yet) {
  for (struct release_event_queue *req = schedule_next(wrld->releaser);
       req != NULL || not_yet >= wrld->releaser->now;
       req = schedule_next(wrld->releaser)) {
    if (req == NULL ||
        !distinguishable(req->release_site->release_prob, MAGIC_PATTERN_PROBABILITY, EPS_C))
      continue;
    if (release_molecules(wrld, req))
      mcell_error("Failed to release molecules of type '%s'.",
                  req->release_site->mol_type->sym->name);
  }
  if (wrld->releaser->error)
    mcell_internal_error("Scheduler reported an out-of-memory error while "
                         "retrieving next scheduled release event, but this "
                         "should never happen.");
}

/***********************************************************************
 process_geometry_changes:

    Produce this round's geometry changes, if any.

 In: state: MCell state 
     not_yet: earliest time which should not yet be processed
 Out: Nothing. Save molecules (with their positions, etc), trash the old
      geometry, initialize new geometry, place all the molecules (moving them
      if necessary to keep them in/out of the appropriate compartments).
 ***********************************************************************/
void process_geometry_changes(struct volume *state, double not_yet) {
  for (struct dg_time_filename *dg_time_fname = schedule_next(
       state->dynamic_geometry_scheduler);
       dg_time_fname != NULL || not_yet >= state->dynamic_geometry_scheduler->now;
       dg_time_fname = schedule_next(state->dynamic_geometry_scheduler)) {

    if (dg_time_fname == NULL)
      continue;
    update_geometry(state, dg_time_fname);
  }
  if (state->dynamic_geometry_scheduler->error)
    mcell_internal_error("Scheduler reported an out-of-memory error while "
                         "retrieving next scheduled geometry change, but this "
                         "should never happen.");
}

/***********************************************************************
 make_checkpoint:

    Produce a checkpoint file.

    In:  struct volume *wrld - the world
    Out: 0 on success, 1 on failure.
         On success, checkpoint file is created.
         On failure, old checkpoint file, if any, is left intact.
 ***********************************************************************/
static int make_checkpoint(struct volume *wrld) {
  /* Make sure we have a filename */
  if (wrld->chkpt_outfile == NULL)
    wrld->chkpt_outfile = CHECKED_SPRINTF("checkpt.%d", getpid());

  /* Print a useful status message */
  switch (wrld->checkpoint_requested) {
  case CHKPT_ITERATIONS_CONT:
  case CHKPT_ALARM_CONT:
    if (wrld->notify->checkpoint_report != NOTIFY_NONE)
      mcell_log("MCell: time = %lld, writing to checkpoint file %s (periodic).",
                wrld->current_iterations, wrld->chkpt_outfile);
    break;

  case CHKPT_ALARM_EXIT:
    if (wrld->notify->checkpoint_report != NOTIFY_NONE)
      mcell_log("MCell: time = %lld, writing to checkpoint file %s (time limit "
                "elapsed).",
                wrld->current_iterations, wrld->chkpt_outfile);
    break;

  case CHKPT_SIGNAL_CONT:
  case CHKPT_SIGNAL_EXIT:
    if (wrld->notify->checkpoint_report != NOTIFY_NONE)
      mcell_log("MCell: time = %lld, writing to checkpoint file %s (user "
                "signal detected).",
                wrld->current_iterations, wrld->chkpt_outfile);
    break;

  case CHKPT_ITERATIONS_EXIT:
    if (wrld->notify->checkpoint_report != NOTIFY_NONE)
      mcell_log("MCell: time = %lld, writing to checkpoint file %s.",
                wrld->current_iterations, wrld->chkpt_outfile);
    break;

  default:
    wrld->checkpoint_requested = CHKPT_NOT_REQUESTED;
    return 0;
  }

  /* Make the checkpoint */
  create_chkpt(wrld, wrld->chkpt_outfile);
  wrld->last_checkpoint_iteration = wrld->current_iterations;

  /* Break out of the loop, if appropriate */
  if (wrld->checkpoint_requested == CHKPT_ALARM_EXIT ||
      wrld->checkpoint_requested == CHKPT_SIGNAL_EXIT ||
      wrld->checkpoint_requested == CHKPT_ITERATIONS_EXIT)
    return 1;

  /* Schedule the next checkpoint, if appropriate */
  if (wrld->checkpoint_requested == CHKPT_ALARM_CONT) {
    if (wrld->continue_after_checkpoint)
      alarm(wrld->checkpoint_alarm_time);
    else
      return 1;
  }

  wrld->checkpoint_requested = CHKPT_NOT_REQUESTED;
  return 0;
}

static double find_next_viz_output_frame(struct frame_data_list *fdl) {
  double next_time = DBL_MAX;
  for (; fdl != NULL; fdl = fdl->next) {
    if (fdl->curr_viz_iteration == NULL)
      continue;

    if (fdl->viz_iteration < next_time)
      next_time = fdl->viz_iteration;
  }

  return next_time;
}

static double find_next_viz_output(struct viz_output_block *vizblk) {
  double next_time = DBL_MAX;
  while (vizblk != NULL) {
    double this_time = find_next_viz_output_frame(vizblk->frame_data_head);
    if (this_time < next_time)
      next_time = this_time;
    vizblk = vizblk->next;
  }

  return next_time;
}

static int print_molecule_collision_report(
    enum notify_level_t molecule_collision_report,
    long long vol_vol_colls,
    long long vol_surf_colls,
    long long surf_surf_colls,
    long long vol_wall_colls,
    long long vol_vol_vol_colls,
    long long vol_vol_surf_colls,
    long long vol_surf_surf_colls,
    long long surf_surf_surf_colls,
    struct reaction_flags *rxn_flags) {
  if (molecule_collision_report == NOTIFY_FULL) {
    mcell_log_raw("\n");
    mcell_log("\tCounts of Reaction Triggered Molecule Collisions");
    mcell_log("(VM = volume molecule, SM = surface molecule, W = wall)");
    if (rxn_flags->vol_vol_reaction_flag) {
      mcell_log("Total number of VM-VM collisions: %lld", vol_vol_colls);
    }
    if (rxn_flags->vol_surf_reaction_flag) {
      mcell_log("Total number of VM-SM collisions: %lld", vol_surf_colls);
    }
    if (rxn_flags->surf_surf_reaction_flag) {
      mcell_log("Total number of SM-SM collisions: %lld", surf_surf_colls);
    }
    if (rxn_flags->vol_wall_reaction_flag) {
      mcell_log("Total number of VM-W collisions: %lld", vol_wall_colls);
    }
    if (rxn_flags->vol_vol_vol_reaction_flag) {
      mcell_log("Total number of VM-VM-VM collisions: %lld", vol_vol_vol_colls);
    }
    if (rxn_flags->vol_vol_surf_reaction_flag) {
      mcell_log("Total number of VM-VM-SM collisions: %lld",
                vol_vol_surf_colls);
    }
    if (rxn_flags->vol_surf_surf_reaction_flag) {
      mcell_log("Total number of VM-SM-SM collisions: %lld",
                vol_surf_surf_colls);
    }
    if (rxn_flags->surf_surf_surf_reaction_flag) {
      mcell_log("Total number of SM-SM-SM collisions: %lld",
                surf_surf_surf_colls);
    }
    mcell_log_raw("\n");
  }

  return 0;
}

/***********************************************************************
 run_sim:

    Simulation main loop.

    In:  pointer to simulation state
    Out: None
 ***********************************************************************/
MCELL_STATUS
mcell_run_simulation(MCELL_STATE *world) {
  if (world->notify->progress_report != NOTIFY_NONE)
    mcell_log("Running simulation.");

  /* If we're reloading a checkpoint, we want to skip all of the processing
  * which happened on the last iteration before checkpointing.  To do this, we
  * skip the first part of run_timestep.
  */
  int restarted_from_checkpoint = 0;
  if (world->start_iterations != 0) {
    assert(world->current_iterations - world->start_iterations == 0.0);
    restarted_from_checkpoint = 1;
  }

  long long frequency = mcell_determine_output_frequency(world);
  int status = 0;
  while (world->current_iterations <= world->iterations) {
    // XXX: A return status of 1 from mcell_run_iterations does not
    // indicate an error but is used to break out of the loop.
    // This behavior is non-conformant and should be changed.
    if (mcell_run_iteration(world, frequency, &restarted_from_checkpoint) ==
        1) {
      break;
    }
  }

  if (mcell_flush_data(world)) {
    mcell_error_nodie("Failed to flush reaction and visualization data.");
    status = 1;
  }

  if (mcell_print_final_warnings(world)) {
    mcell_error_nodie("Failed to print final warnings.");
    status = 1;
  }

  if (mcell_print_final_statistics(world)) {
    mcell_error_nodie("Failed to print final statistics.");
    status = 1;
  }

  return status;
}

/**************************************************************************
*
* Run multiple iterations at once
*
 *************************************************************************/
MCELL_STATUS
mcell_run_n_iterations(MCELL_STATE *world, long long frequency,
                    int *restarted_from_checkpoint, int n_iter) {

  int i_iter = 0;
  while (i_iter < n_iter) {
    // XXX: A return status of 1 from mcell_run_iterations does not
    // indicate an error but is used to break out of the loop.
    // This behavior is non-conformant and should be changed.
    if (mcell_run_iteration(world, frequency, restarted_from_checkpoint) == 1) {
      break;
    }
    i_iter += 1;
  }

  return 0;
}

/**************************************************************************
 *
 * this function runs a single mcell iteration
 *
 * parameters:
 *    frequency                : the requested output frequency
 *    restart_from_checkpoint  : does this iteration follow a checkpoint
 *
 *************************************************************************/
MCELL_STATUS
mcell_run_iteration(MCELL_STATE *world, long long frequency,
                    int *restarted_from_checkpoint) {
  emergency_output_hook_enabled = 1;

  long long iter_report_phase = world->current_iterations % frequency;
  double not_yet = world->current_iterations + 1.0;

  if (world->current_iterations != 0)
    world->elapsed_time = world->current_iterations;
  else
    world->elapsed_time = 1.0;

  if (!*restarted_from_checkpoint) {

    /* Change geometry if needed */
    process_geometry_changes(world, not_yet);

    /* Release molecules */
    process_molecule_releases(world, not_yet);

    /* Produce output */
    process_reaction_output(world, not_yet);
    process_volume_output(world, not_yet);
    for (struct viz_output_block *vizblk = world->viz_blocks; vizblk != NULL;
         vizblk = vizblk->next) {
      if (vizblk->frame_data_head && update_frame_data_list(world, vizblk))
        mcell_error("Unknown error while updating frame data list.");
    }

    /* Produce iteration report */
    if (iter_report_phase == 0 &&
        world->notify->iteration_report != NOTIFY_NONE) {
      mcell_log_raw("Iterations: %lld of %lld ", world->current_iterations,
                    world->iterations);

      if (world->notify->throughput_report != NOTIFY_NONE) {
        struct timeval cur_time;
        gettimeofday(&cur_time, NULL);
        if (world->last_timing_time.tv_sec > 0) {
          double time_diff =
              (double)(cur_time.tv_sec - world->last_timing_time.tv_sec) *
                  1000000.0 +
              (double)(cur_time.tv_usec - world->last_timing_time.tv_usec);
          time_diff /= (double)(world->current_iterations - world->last_timing_iteration);
          mcell_log_raw(" (%.6lg iter/sec)", 1000000.0 / time_diff);
          world->last_timing_iteration = world->current_iterations;
          world->last_timing_time = cur_time;
        } else {
          world->last_timing_iteration = world->current_iterations;
          world->last_timing_time = cur_time;
        }
      }

      mcell_log_raw("\n");
    }

    /* Check for a checkpoint on this iteration */
    if (world->chkpt_iterations && world->current_iterations != world->start_iterations &&
        ((world->current_iterations - world->start_iterations) % world->chkpt_iterations == 0)) {
      if (world->continue_after_checkpoint) {
        world->checkpoint_requested = CHKPT_ITERATIONS_CONT;
      } else {
        world->checkpoint_requested = CHKPT_ITERATIONS_EXIT;
      }
    }

    /* No checkpoint signalled.  Keep going. */
    if (world->checkpoint_requested != CHKPT_NOT_REQUESTED) {
      // This won't work with (non-trad) PBCs until we start saving the
      // molecules' periodic box in the checkpoint file. In principle, it
      // should probably work with the traditional form, but I'm disabling it
      // here to be safe.
      if (world->periodic_box_obj) {
        mcell_error(
          "periodic boundary conditions do not currently work with "
          "checkpointing.");
      }
      /* Make a checkpoint, exiting the loop if necessary */
      if (make_checkpoint(world))
        return 1;
    }

    /* Even if no checkpoint, the last iteration is a half-iteration. */
    if (world->current_iterations >= world->iterations)
      return 1;
  }

  // reset this flag to zero
  *restarted_from_checkpoint = 0;

  run_concentration_clamp(world, world->current_iterations);

  double next_release_time;
  if (!schedule_anticipate(world->releaser, &next_release_time))
    next_release_time = world->iterations + 1;

  if (next_release_time < world->current_iterations + 1)
    next_release_time = world->current_iterations + 1;

  double next_vol_output;
  if (!schedule_anticipate(world->volume_output_scheduler, &next_vol_output))
    next_vol_output = world->iterations + 1;
  double next_viz_output = find_next_viz_output(world->viz_blocks);
  double next_barrier =
      min3d(next_release_time, next_vol_output, next_viz_output);

  while (world->storage_head != NULL &&
         world->storage_head->store->current_time <= not_yet) {
    int done = 0;
    while (!done) {
      done = 1;
      for (struct storage_list *local = world->storage_head; local != NULL;
           local = local->next) {
        if (local->store->timer->current != NULL) {
          run_timestep(world, local->store, next_barrier,
                       (double)world->iterations + 1.0);
          done = 0;
        }
      }
    }

    for (struct storage_list *local = world->storage_head; local != NULL;
         local = local->next) {
      /* Not using the return value -- just trying to advance the scheduler */
      void *o = schedule_next(local->store->timer);
      if (o != NULL)
        mcell_internal_error("Scheduler dropped a molecule on the floor!");
      local->store->current_time += 1.0;
    }
  }

  world->current_iterations++;

  return 0;
}

/**************************************************************************
 *
 * This function returns the output frequency.
 *
 * The returned value is either the one requested in the mdl file or
 * is determined via a heuristic.
 *
 **************************************************************************/
long long mcell_determine_output_frequency(MCELL_STATE *world) {
  long long frequency;
  if (world->notify->custom_iteration_value != 0) {
    frequency = world->notify->custom_iteration_value;
  } else {
    if (world->iterations < 10)
      frequency = 1;
    else if (world->iterations < 1000)
      frequency = 10;
    else if (world->iterations < 100000)
      frequency = 100;
    else if (world->iterations < 10000000)
      frequency = 1000;
    else if (world->iterations < 1000000000)
      frequency = 10000;
    else
      frequency = 100000;
  }

  return frequency;
}

/************************************************************************
 *
 * function responsible for flushing any remaining reaction
 * and viz data output to disk. Also sets a final simulation checkpoint
 * if necessary.
 *
 ***********************************************************************/
MCELL_STATUS
mcell_flush_data(MCELL_STATE *world) {
  int status = 0;
  if (world->chkpt_iterations &&
      world->current_iterations > world->last_checkpoint_iteration) {
    status = make_checkpoint(world);
  }

  emergency_output_hook_enabled = 0;
  int num_errors = flush_reaction_output(world);
  if (num_errors != 0) {
    mcell_warn("%d errors occurred while flushing buffered reaction output.\n"
               "  Simulation complete anyway--continuing as normal.",
               num_errors);
    status = 1;
  }

  if (world->notify->progress_report != NOTIFY_NONE)
    mcell_log("Exiting run loop.");

  int warned = 0;
  for (struct viz_output_block *vizblk = world->viz_blocks; vizblk != NULL;
       vizblk = vizblk->next) {
    if (finalize_viz_output(world, vizblk) && !warned) {
      mcell_warn("VIZ output was not successfully finalized.\n"
                 "  Visualization of results may not work correctly.");
      warned = 1;
      status = 1;
    }
  }

  return status;
}

/***********************************************************************
 *
 * function printing final simulation warning about missed reactions,
 * high reaction probabilities, etc.
 *
 ***********************************************************************/
MCELL_STATUS
mcell_print_final_warnings(MCELL_STATE *world) {
  int first_report = 1;
  if (world->notify->missed_reactions != WARN_COPE) {
    for (int i = 0; i < world->rx_hashsize; i++) {
      for (struct rxn *rxp = world->reaction_hash[i]; rxp != NULL;
           rxp = rxp->next) {
        int do_not_print = 0;
        /* skip printing messages if one of the reactants is a surface */
        for (unsigned int j = 0; j < rxp->n_reactants; j++) {
          if ((rxp->players[j]->flags & IS_SURFACE) != 0)
            do_not_print = 1;
        }

        if (do_not_print == 1)
          continue;
        if (rxp->n_occurred * world->notify->missed_reaction_value <
            rxp->n_skipped) {
          if (first_report) {
            mcell_log("Warning: Some reactions were missed because reaction "
                      "probability exceeded 1.");
            first_report = 0;
          }
          mcell_log_raw("  ");
          for (unsigned int j = 0; j < rxp->n_reactants; j++) {
            mcell_log_raw("%s%s[%d]", j ? " + " : "",
                          rxp->players[j]->sym->name, rxp->geometries[j]);
          }
          mcell_log_raw("  --  %g%% of reactions missed.\n",
                        0.001 * round(1000 * rxp->n_skipped * 100 /
                                      (rxp->n_skipped + rxp->n_occurred)));
        }
      }
    }
    if (!first_report)
      mcell_log_raw("\n");
  }

  first_report += 1;

  if (world->notify->short_lifetime != WARN_COPE) {
    for (int i = 0; i < world->n_species; i++) {
      if (world->species_list[i] == world->all_mols)
        continue;
      if ((world->species_list[i] == world->all_volume_mols) ||
          (world->species_list[i] == world->all_surface_mols))
        continue;
      if (world->species_list[i]->n_deceased <= 0)
        continue;

      double mean_lifetime_seconds = world->species_list[i]->cum_lifetime_seconds /
                 world->species_list[i]->n_deceased;
      double short_lifetime_seconds = convert_iterations_to_seconds(
          world->start_iterations, world->time_unit,
          world->simulation_start_seconds, world->notify->short_lifetime_value);
      double mean_lifetime_timesteps = mean_lifetime_seconds/world->time_unit;
      if (mean_lifetime_seconds < short_lifetime_seconds) {
        if (first_report) {
          if (first_report > 1)
            mcell_log_raw("\n");
          mcell_log("Warning: Some molecules had a lifetime short relative to "
                    "the timestep.");
          first_report = 0;
        }
        mcell_log("  Mean lifetime of %s was %g timesteps.",
                  world->species_list[i]->sym->name, mean_lifetime_timesteps);
      }
    }
    if (!first_report)
      mcell_log_raw("\n");
  }

  return 0;
}

/***********************************************************************
 *
 * function printing final simulation statistics and timings
 *
 ***********************************************************************/
MCELL_STATUS
mcell_print_final_statistics(MCELL_STATE *world) {
  if (world->reaction_prob_limit_flag)
    mcell_warn("during the simulation some reaction probabilities were greater "
               "than 1. You may want to rerun the simulation with the WARNINGS "
               "block enabled to get more detail.\n");

  if (world->notify->final_summary == NOTIFY_FULL) {
    mcell_log("iterations = %lld ; elapsed time = %1.15g seconds",
              world->current_iterations,
              world->chkpt_start_time_seconds +
                  ((world->current_iterations - world->start_iterations) * world->time_unit));

    if (world->diffusion_number > 0)
      mcell_log("Average diffusion jump was %.2f timesteps\n",
                world->diffusion_cumtime / (double)world->diffusion_number);
    mcell_log("Total number of random number use: %lld", rng_uses(world->rng));
    mcell_log("Total number of ray-subvolume intersection tests: %lld",
              world->ray_voxel_tests);
    mcell_log("Total number of ray-polygon intersection tests: %lld",
              world->ray_polygon_tests);
    mcell_log("Total number of ray-polygon intersections: %lld",
              world->ray_polygon_colls);
    mcell_log("Total number of dynamic geometry molecule displacements: %lld",
              world->dyngeom_molec_displacements);
    print_molecule_collision_report(
        world->notify->molecule_collision_report,
        world->vol_vol_colls,
        world->vol_surf_colls,
        world->surf_surf_colls,
        world->vol_wall_colls,
        world->vol_vol_vol_colls,
        world->vol_vol_surf_colls,
        world->vol_surf_surf_colls,
        world->surf_surf_surf_colls,
        &world->rxn_flags);

    struct rusage run_time = { .ru_utime = { 0, 0 }, .ru_stime = { 0, 0 } };
    time_t t_end; /* global end time of MCell run */
    double u_init_time,
        s_init_time;               /* initialization time (user) and (system) */
    double u_run_time, s_run_time; /* run time (user) and (system) */

    u_init_time = world->u_init_time.tv_sec +
                  (world->u_init_time.tv_usec / MAX_TARGET_TIMESTEP);
    s_init_time = world->s_init_time.tv_sec +
                  (world->s_init_time.tv_usec / MAX_TARGET_TIMESTEP);

    mcell_log("Initialization CPU time = %f (user) and %f (system)",
              u_init_time, s_init_time);

    getrusage(RUSAGE_SELF, &run_time);
    u_run_time = run_time.ru_utime.tv_sec +
                 (run_time.ru_utime.tv_usec / MAX_TARGET_TIMESTEP);
    s_run_time = run_time.ru_stime.tv_sec +
                 (run_time.ru_stime.tv_usec / MAX_TARGET_TIMESTEP);

    mcell_log("Simulation CPU time = %f (user) and %f (system)",
              u_run_time - u_init_time, s_run_time - s_init_time);
    t_end = time(NULL);
    mcell_log("Total wall clock time = %ld seconds",
              (long)difftime(t_end, world->t_start));
  }

  return 0;
}
