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

#include <assert.h>
#include <float.h>
#include <math.h>
#include <unistd.h>

#ifndef _WIN32
#include <sys/resource.h>
#endif
#include <stdlib.h>
#if defined(__linux__)
#include <fenv.h>
#endif

#include "delayed_trigger.h"
#include "sym_table.h"
#include "logging.h"
#include "rng.h"
#include "strfunc.h"
#include "vol_util.h"
#include "react_output.h"
#include "viz_output.h"
#include "volume_output.h"
#include "diffuse.h"
#include "init.h"
#include "chkpt.h"
#include "version_info.h"
#include "argparse.h"

#include "mcell_run.h"
#include "mcell_init.h"

// helper struct to pass data to worker threads
struct worker_data {
  struct volume *global_state;
  thread_state_t *thread_state;
};


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

#if 0
/*
 * Get the log frequency specified in this world.
 */
static long long get_log_frequency(struct volume *wrld) {
  if (wrld->notify->custom_iteration_value != 0) {
   return ( wrld->notify->custom_iteration_value );
  }

  if (wrld->iterations < 10)              return ( 1 );
  else if (wrld->iterations < 1000)       return ( 10 );
  else if (wrld->iterations < 100000)     return ( 100 );
  else if (wrld->iterations < 10000000)   return ( 1000 );
  else if (wrld->iterations < 1000000000) return ( 10000 );
  else                                    return ( 100000 );
}
#endif

/*
 * Produce the iteration/timing report.
 */
static void produce_iteration_report(struct volume *wrld, struct timing_info *timing,
  long long frequency) {

  if (wrld->notify->iteration_report == NOTIFY_NONE)
    return;

  if (timing->iter_report_phase == 0) {
    mcell_log_raw("Iterations: %lld of %lld ", wrld->current_iterations, wrld->iterations);

    if (wrld->notify->throughput_report != NOTIFY_NONE) {
      struct timeval cur_time;
      gettimeofday(&cur_time, NULL);
      if (timing->last_timing_time.tv_sec > 0) {
        double time_diff = (double) (cur_time.tv_sec - timing->last_timing_time.tv_sec) *
          1000000.0 + (double) (cur_time.tv_usec - timing->last_timing_time.tv_usec);
        time_diff /= (double)(wrld->current_iterations - timing->last_timing_iteration);
        mcell_log_raw(" (%.6lg iter/sec)", 1000000.0 / time_diff);
        timing->last_timing_iteration = wrld->current_iterations;
        timing->last_timing_time = cur_time;
      } else {
        timing->last_timing_iteration = wrld->current_iterations;
        timing->last_timing_time = cur_time;
      }
    }

    mcell_log_raw("\n");
  }

  if (++ timing->iter_report_phase == frequency) {
    timing->iter_report_phase = 0;
  }
}


#ifdef MCELL_COMPLETE_END_STATS
typedef struct {
  double mean_ray_voxel_tests;
  double mean_ray_polygon_tests;
  double mean_ray_polygon_colls;
  double mean_mol_mol_colls;
  double mean_mol_grid_colls;
  double mean_grid_grid_colls;
  double mean_mol_wall_colls;
  double mean_mol_mol_mol_colls;
  double mean_mol_mol_grid_colls;
  double mean_mol_grid_grid_colls;
  double mean_grid_grid_grid_colls;
  double mean_random_numbers;

  double stdev_ray_voxel_tests;
  double stdev_ray_polygon_tests;
  double stdev_ray_polygon_colls;
  double stdev_mol_mol_colls;
  double stdev_mol_grid_colls;
  double stdev_grid_grid_colls;
  double stdev_mol_wall_colls;
  double stdev_mol_mol_mol_colls;
  double stdev_mol_mol_grid_colls;
  double stdev_mol_grid_grid_colls;
  double stdev_grid_grid_grid_colls;
  double stdev_random_numbers;
} runtime_mean_variance_t;

static int print_molecule_collision_report_full(struct volume *world,
  runtime_statistics_t *overall, runtime_mean_variance_t *mean_var,
  runtime_statistics_t *mins, runtime_statistics_t *maxes,
  runtime_statistics_t *zeros) {

  mcell_log_raw("\n");
  mcell_log("\tCounts of Reaction Triggered Molecule Collisions");
  mcell_log("(VM = volume molecule, SM = surface molecule, W = wall)");

#define PRINT_REPORT(type, lbl) do {                                      \
  if (world->rxn_flags.type##_reaction_flag) {                                      \
    mcell_log("Total number of " lbl " collisions: %lld (per-subdiv: "    \
      "mean=%.2f, stdev=%.2f, min=%lld, max=%lld, zero=%d)",              \
              overall->type##_colls,                                      \
              mean_var->mean_type##_colls,                                \
              mean_var->stdev_type##_colls,                               \
              mins->type##_colls,                                         \
              maxes->type##_colls,                                        \
              (int) zeros->type##_colls);                                 \
  } } while(0)
  PRINT_REPORT(vol_vol,        "VM-VM");
  PRINT_REPORT(vol_surf,       "VM-SM");
  PRINT_REPORT(surf_surf,      "SM-SM");
  PRINT_REPORT(vol_surf,       "VM-W");
  PRINT_REPORT(vol_vol_vol,    "VM-VM-VM");
  PRINT_REPORT(vol_vol_surf,   "VM-VM-SM");
  PRINT_REPORT(vol_surf_surf,  "VM-SM-SM");
  PRINT_REPORT(surf_surf_surf, "SM-SM-SM");
#undef PRINT_REPORT

  mcell_log_raw("\n");
  return ( 0 );
}
#else


static int print_molecule_collision_report(struct volume *world,
  runtime_statistics_t *stats) {

  if (world->notify->molecule_collision_report == NOTIFY_FULL) {
    mcell_log_raw("\n");
    mcell_log("\tCounts of Reaction Triggered Molecule Collisions");
    mcell_log("(VM = volume molecule, SM = surface molecule, W = wall)");

#define PRINT_REPORT(type, lbl) do {                                      \
  if (world->rxn_flags.type##_reaction_flag) {                                      \
    mcell_log("Total number of " lbl " collisions: %lld",                 \
              stats->type##_colls);                                       \
  } } while(0)
  PRINT_REPORT(vol_vol,        "VM-VM");
  PRINT_REPORT(vol_surf,       "VM-SM");
  PRINT_REPORT(surf_surf,      "SM-SM");
  PRINT_REPORT(vol_surf,       "VM-W");
  PRINT_REPORT(vol_vol_vol,    "VM-VM-VM");
  PRINT_REPORT(vol_vol_surf,   "VM-VM-SM");
  PRINT_REPORT(vol_surf_surf,  "VM-SM-SM");
  PRINT_REPORT(surf_surf_surf, "SM-SM-SM");
#undef PRINT_REPORT
      mcell_log_raw("\n");
  }
  return 0;
}
#endif

/**
 * Display a final statistics report for this run.
 *
 * Currently, the full end-stats are compile-time selectable.
 * If they are useful to anyone else, we can make them runtime
 * selectable, but they probably need some formatting work so
 * that they aren't hideous.
 *
 * Apologies in advance for the following lengthy function,
 * which computes fairly extensive run-time statistics over
 * the storages.  Macros have been used to make the code more
 * readable, wherever possible.
 *
 * FIXME (Markus): Check if we really need these macros. Replace by
 * functions if necessary.
 */
static void display_final_summary(struct volume *world) {
  struct rusage run_time;
  double u_run_time, s_run_time; /* run time (user) and (system) */
  double u_init_time, s_init_time; /* initialization time (user) and (system) */

  mcell_log("iterations = %lld ; elapsed time = %1.15g seconds",
   world->current_iterations, world->chkpt_start_time_seconds +
   ((world->current_iterations - world->start_iterations)*world->time_unit));

  runtime_statistics_t overall_stats;
  overall_stats = world->subdivisions[0].stats;

#ifdef MCELL_COMPLETE_END_STATS
  /* Initialize aggregate statistics. */
  /* These aggregate statistics are computed, largely, by
   * initializing each statistics block to the statistics from
   * the first storage, and then processing them with each
   * subsequent storage -- for overall stats, using summation,
   * for min/max stats using min/max functions, etc.
   */
  world->subdivisions[0].stats.random_numbers = rng_uses(world->subdivisions[0].rng);
  runtime_statistics_t min_stats = world->subdivisions[0].stats;
  runtime_statistics_t max_stats = world->subdivisions[0].stats;
  runtime_statistics_t num_zeros;
  runtime_mean_variance_t mean_var;
  num_zeros.diffusion_number =
    (world->subdivisions[0].stats.diffusion_number == 0) ? 1 : 0;
  num_zeros.ray_voxel_tests =
    (world->subdivisions[0].stats.ray_voxel_tests == 0) ? 1 : 0;
  num_zeros.ray_polygon_tests =
    (world->subdivisions[0].stats.ray_polygon_tests == 0) ? 1 : 0;
  num_zeros.ray_polygon_colls =
    (world->subdivisions[0].stats.ray_polygon_colls == 0) ? 1 : 0;
  num_zeros.mol_mol_colls =
    (world->subdivisions[0].stats.mol_mol_colls == 0) ? 1 : 0;
  num_zeros.mol_grid_colls =
    (world->subdivisions[0].stats.mol_grid_colls == 0) ? 1 : 0;
  num_zeros.grid_grid_colls =
    (world->subdivisions[0].stats.mol_grid_colls == 0) ? 1 : 0;
  num_zeros.mol_wall_colls =
    (world->subdivisions[0].stats.mol_wall_colls == 0) ? 1 : 0;
  num_zeros.mol_mol_mol_colls =
    (world->subdivisions[0].stats.mol_mol_mol_colls == 0) ? 1 : 0;
  num_zeros.mol_mol_grid_colls =
    (world->subdivisions[0].stats.mol_mol_grid_colls == 0) ? 1 : 0;
  num_zeros.mol_grid_grid_colls =
    (world->subdivisions[0].stats.mol_grid_grid_colls == 0) ? 1 : 0;
  num_zeros.grid_grid_grid_colls =
    (world->subdivisions[0].stats.grid_grid_grid_colls == 0) ? 1 : 0;
  num_zeros.random_numbers =
    (world->subdivisions[0].stats.random_numbers == 0) ? 1 : 0;

  /* Compute aggregate statistics. */
#define UPDATE_OVERALL(s1, s2, fld) do { (s1).fld += (s2).fld; } while (0)

#define UPDATE_ZERO(s1, s2, fld)    do {                                    \
  if ((s2).fld == 0) ++ (s1).fld;                                           \
} while (0)

#define UPDATE_MIN(s1, s2, fld)     do {                                    \
  if ((s2).fld != 0  &&  (s1).fld > (s2).fld) (s1).fld = (s2).fld;          \
} while (0)

#define UPDATE_MAX(s1, s2, fld)     do {                                    \
  if ((s2).fld != 0  &&  (s1).fld < (s2).fld) (s1).fld = (s2).fld;          \
} while (0)

  for (int i=1; i<world->num_subdivisions; ++i) {
    struct storage *store = & world->subdivisions[i];
    runtime_statistics_t *substat = & store->stats;
    substat->random_numbers = rng_uses(store->rng);

    /* Update zeros. */
    UPDATE_ZERO(num_zeros, *substat, diffusion_number);
    UPDATE_ZERO(num_zeros, *substat, ray_voxel_tests);
    UPDATE_ZERO(num_zeros, *substat, ray_polygon_tests);
    UPDATE_ZERO(num_zeros, *substat, ray_polygon_colls);
    UPDATE_ZERO(num_zeros, *substat, vol_vol_colls);
    UPDATE_ZERO(num_zeros, *substat, vol_surf_colls);
    UPDATE_ZERO(num_zeros, *substat, surf_surf_colls);
    UPDATE_ZERO(num_zeros, *substat, vol_wall_colls);
    update_zero(num_zeros, *substat, vol_vol_vol_colls);
    UPDATE_ZERO(num_zeros, *substat, vol_vol_surf_colls);
    UPDATE_ZERO(num_zeros, *substat, vol_surf_surf_colls);
    UPDATE_ZERO(num_zeros, *substat, surf_surf_surf_colls);
    UPDATE_ZERO(num_zeros, *substat, random_numbers);

    /* Update overall stats. */
    UPDATE_OVERALL(overall_stats, *substat, diffusion_number);
    UPDATE_OVERALL(overall_stats, *substat, diffusion_cumtime);
    UPDATE_OVERALL(overall_stats, *substat, ray_voxel_tests);
    UPDATE_OVERALL(overall_stats, *substat, ray_polygon_tests);
    UPDATE_OVERALL(overall_stats, *substat, ray_polygon_colls);
    UPDATE_OVERALL(overall_stats, *substat, vol_vol_colls);
    UPDATE_OVERALL(overall_stats, *substat, vol_surf_colls);
    UPDATE_OVERALL(overall_stats, *substat, surf_surf_colls);
    UPDATE_OVERALL(overall_stats, *substat, vol_surf_colls);
    UPDATE_OVERALL(overall_stats, *substat, vol_vol_vol_colls);
    UPDATE_OVERALL(overall_stats, *substat, vol_vol_surf_colls);
    UPDATE_OVERALL(overall_stats, *substat, vol_surf_surf_colls);
    UPDATE_OVERALL(overall_stats, *substat, surf_surf_surf_colls);
    UPDATE_OVERALL(overall_stats, *substat, random_numbers);

    /* Update minimum stats. */
    if (substat->diffusion_number != 0
      && (min_stats.diffusion_cumtime / (double) min_stats.diffusion_number) >
     (substat->diffusion_cumtime  / (double) substat->diffusion_number)) {

      min_stats.diffusion_cumtime = substat->diffusion_cumtime;
      min_stats.diffusion_number  = substat->diffusion_number;
    }
    UPDATE_MIN(min_stats, *substat, ray_voxel_tests);
    UPDATE_MIN(min_stats, *substat, ray_polygon_tests);
    UPDATE_MIN(min_stats, *substat, ray_polygon_colls);
    UPDATE_MIN(min_stats, *substat, vol_vol_colls);
    UPDATE_MIN(min_stats, *substat, vol_surf_colls);
    UPDATE_MIN(min_stats, *substat, surf_surf_colls);
    UPDATE_MIN(min_stats, *substat, vol_wall_colls);
    UPDATE_MIN(min_stats, *substat, vol_vol_vol_colls);
    UPDATE_MIN(min_stats, *substat, vol_vol_surf_colls);
    UPDATE_MIN(min_stats, *substat, vol_surf_surf_colls);
    UPDATE_MIN(min_stats, *substat, surf_surf_surf_colls);
    UPDATE_MIN(min_stats, *substat, random_numbers);

    /* Update maximum stats. */
    if (substat->diffusion_number != 0
        && (max_stats.diffusion_cumtime / (double) max_stats.diffusion_number) <
           (substat->diffusion_cumtime  / (double) substat->diffusion_number)) {

      max_stats.diffusion_cumtime = substat->diffusion_cumtime;
      max_stats.diffusion_number  = substat->diffusion_number;
    }
    UPDATE_MAX(max_stats, *substat, ray_voxel_tests);
    UPDATE_MAX(max_stats, *substat, ray_polygon_tests);
    UPDATE_MAX(max_stats, *substat, ray_polygon_colls);
    UPDATE_MAX(max_stats, *substat, vol_vol_colls);
    UPDATE_MAX(max_stats, *substat, vol_surf_colls);
    UPDATE_MAX(max_stats, *substat, surf_surf_colls);
    UPDATE_MAX(max_stats, *substat, vol_wall_colls);
    UPDATE_MAX(max_stats, *substat, vol_vol_vol_colls);
    UPDATE_MAX(max_stats, *substat, vol_vol_surf_colls);
    UPDATE_MAX(max_stats, *substat, vol_surf_surf_colls);
    UPDATE_MAX(max_stats, *substat, surf_surf_surf_colls);
    UPDATE_MAX(max_stats, *substat, random_numbers);
  }
#undef UPDATE_ZERO
#undef UPDATE_OVERALL
#undef UPDATE_MIN
#undef UPDATE_MAX

  /* Compute means. */
#define COMPUTE_MEAN(st) (overall_stats.st / (double) (world->num_subdivisions - num_zeros.st))
  double mean_diffusion_jump = overall_stats.diffusion_cumtime /
    (double) overall_stats.diffusion_number;
  double min_mean_diffusion_jump = min_stats.diffusion_cumtime /
    (double) min_stats.diffusion_number;
  double max_mean_diffusion_jump = max_stats.diffusion_cumtime /
    (double) max_stats.diffusion_number;
  mean_var.mean_ray_voxel_tests = COMPUTE_MEAN(ray_voxel_tests);
  mean_var.mean_ray_polygon_tests = COMPUTE_MEAN(ray_polygon_tests);
  mean_var.mean_ray_polygon_colls = COMPUTE_MEAN(ray_polygon_colls);
  mean_var.mean_mol_mol_colls = COMPUTE_MEAN(vol_vol_colls);
  mean_var.mean_mol_grid_colls = COMPUTE_MEAN(vol_surf_colls);
  mean_var.mean_grid_grid_colls = COMPUTE_MEAN(surf_surf_colls);
  mean_var.mean_mol_wall_colls = COMPUTE_MEAN(vol_wall_colls);
  mean_var.mean_mol_mol_mol_colls = COMPUTE_MEAN(vol_vol_vol_colls);
  mean_var.mean_mol_mol_grid_colls = COMPUTE_MEAN(vol_vol_surf_colls);
  mean_var.mean_mol_grid_grid_colls = COMPUTE_MEAN(vol_surf_surf_colls);
  mean_var.mean_grid_grid_grid_colls = COMPUTE_MEAN(surf_surf_surf_colls);
  mean_var.mean_random_numbers = COMPUTE_MEAN(random_numbers);
#undef COMPUTE_MEAN

  /* Compute stdevs. */
#define COMPUTE_STDEV(st) do {                                              \
  stdev_##st = 0;                                                           \
  int num_samples = world->num_subdivisions - num_zeros.st;                 \
  if (num_samples > 1) {                                                    \
    for (int i=0; i<world->num_subdivisions; ++i) {                         \
      double diff = world->subdivisions[i].stats.st - mean_##st;            \
      mean_var.stdev_##st += diff*diff;                                     \
    }                                                                       \
    mean_var.stdev_##st /= (double) (num_samples - 1);                      \
  }                                                                         \
} while (0)
  COMPUTE_STDEV(ray_voxel_tests);
  COMPUTE_STDEV(ray_polygon_tests);
  COMPUTE_STDEV(ray_polygon_colls);
  COMPUTE_STDEV(vol_vol_colls);
  COMPUTE_STDEV(vol_surf_colls);
  COMPUTE_STDEV(surf_surf_colls);
  COMPUTE_STDEV(vol_wall_colls);
  COMPUTE_STDEV(vol_vol_vol_colls);
  COMPUTE_STDEV(vol_vol_surf_colls);
  COMPUTE_STDEV(vol_surf_surf_colls);
  COMPUTE_STDEV(surf_surf_surf_colls);
  COMPUTE_STDEV(random_numbers);
#undef COMPUTE_STDEV

  if (overall_stats.diffusion_number > 0) {
    mcell_log("Average diffusion jump was %.2f timesteps (per-subdiv: min=%.2f, "
      "max=%.2f, zero=%d)\n", mean_diffusion_jump, min_mean_diffusion_jump,
      max_mean_diffusion_jump, (int) num_zeros.diffusion_number);
  }
  mcell_log("Total number of random number use: %lld (per-subdiv: total=%lld, "
    "mean=%.2f, stdev=%.2f, min=%lld, max=%lld, zero=%d)",
    overall_stats.random_numbers + rng_uses(world->rng_global),
    overall_stats.random_numbers, mean_var.mean_random_numbers,
    mean_var.stdev_random_numbers, min_stats.random_numbers,
    max_stats.random_numbers, (int) num_zeros.random_numbers);
  mcell_log("Total number of ray-subvolume intersection tests: %lld (per-subdiv: "
    "mean=%.2f, stdev=%.2f, min=%lld, max=%lld, zero=%d)",
    overall_stats.ray_voxel_tests, mean_var.mean_ray_voxel_tests,
    mean_var.stdev_ray_voxel_tests, min_stats.ray_voxel_tests,
    max_stats.ray_voxel_tests, (int) num_zeros.ray_voxel_tests);
  mcell_log("Total number of ray-polygon intersection tests: %lld (per-subdiv: "
    "mean=%.2f, stdev=%.2f, min=%lld, max=%lld, zero=%d)", overall_stats.ray_polygon_tests,
    mean_var.mean_ray_polygon_tests, mean_var.stdev_ray_polygon_tests,
    min_stats.ray_polygon_tests, max_stats.ray_polygon_tests,
    (int) num_zeros.ray_polygon_tests);
  mcell_log("Total number of ray-polygon intersections: %lld (per-subdiv: mean=%.2f, ",
    "stdev=%.2f, min=%lld, max=%lld, zero=%d)", overall_stats.ray_polygon_colls,
    mean_var.mean_ray_polygon_colls, mean_var.stdev_ray_polygon_colls,
    min_stats.ray_polygon_colls, max_stats.ray_polygon_colls,
    (int) num_zeros.ray_polygon_colls);
  print_molecule_collision_report_full(world, & overall_stats, & mean_var, & min_stats,
    & max_stats, & num_zeros);
#else
#define UPDATE_OVERALL(s1, s2, fld) (s1).fld += (s2).fld
  for (int i=1; i<world->num_subdivisions; ++i) {

    struct storage *store = & world->subdivisions[i];
    runtime_statistics_t *substat = & store->stats;
    substat->random_numbers = rng_uses(store->rng);

    /* Update overall stats. */
    UPDATE_OVERALL(overall_stats, *substat, diffusion_number);
    UPDATE_OVERALL(overall_stats, *substat, diffusion_cumtime);
    UPDATE_OVERALL(overall_stats, *substat, ray_voxel_tests);
    UPDATE_OVERALL(overall_stats, *substat, ray_polygon_tests);
    UPDATE_OVERALL(overall_stats, *substat, ray_polygon_colls);
    UPDATE_OVERALL(overall_stats, *substat, vol_vol_colls);
    UPDATE_OVERALL(overall_stats, *substat, vol_surf_colls);
    UPDATE_OVERALL(overall_stats, *substat, surf_surf_colls);
    UPDATE_OVERALL(overall_stats, *substat, vol_wall_colls);
    UPDATE_OVERALL(overall_stats, *substat, vol_vol_vol_colls);
    UPDATE_OVERALL(overall_stats, *substat, vol_vol_surf_colls);
    UPDATE_OVERALL(overall_stats, *substat, vol_surf_surf_colls);
    UPDATE_OVERALL(overall_stats, *substat, surf_surf_surf_colls);
    UPDATE_OVERALL(overall_stats, *substat, random_numbers);
  }
#undef UPDATE_OVERALL

  print_molecule_collision_report(world, &overall_stats);
#endif

  long t_final = time(NULL);
  getrusage(RUSAGE_SELF, & run_time);

  /* Compute user/system time for initialization and simulation. */
  u_init_time = world->u_init_time.tv_sec + (world->u_init_time.tv_usec/MAX_TARGET_TIMESTEP);
  s_init_time = world->s_init_time.tv_sec + (world->s_init_time.tv_usec/MAX_TARGET_TIMESTEP);
  u_run_time = run_time.ru_utime.tv_sec + (run_time.ru_utime.tv_usec/MAX_TARGET_TIMESTEP);
  s_run_time = run_time.ru_stime.tv_sec + (run_time.ru_stime.tv_usec/MAX_TARGET_TIMESTEP);

  mcell_log("Initialization CPU time = %f (user) and %f (system)",
            u_init_time, s_init_time);
  mcell_log("Simulation CPU time = %f (user) and %f (system)",
            u_run_time - u_init_time, s_run_time - s_init_time);
  mcell_log("Total wall clock time = %ld seconds",
            (long) difftime(t_final, world->t_start) );
}


static int is_subdivision_complete(struct storage *nation) {
  if (nation->inbound != NULL && nation->inbound->fill != 0) {
    return 0;
  }

  if (nation->timer->current_count != 0) {
    return 0;
  }

  return 1;
}


/* is_iteration_complete returns 1 when there are no blocked
 * or queued tasks left to work on for this iteration. Otherwise
 * the function returns 0 */
static int is_iteration_complete(struct volume *wrld) {
  if (wrld->task_queue.ready_head != NULL) {
    return 0;
  }

  if (wrld->task_queue.blocked_head != NULL) {
    return 0;
  }

  return 1;
}


/* wake_worker_pool tells all workers that there is work to do via a
 * broadcast on the dispatch_ready condition variable */
static void wake_worker_pool(struct volume *wrld) {
  pthread_cond_broadcast(& wrld->dispatch_ready);
  pthread_mutex_unlock(& wrld->dispatch_lock);
}


/* wait_for_sequential_section waits on the dispatch_empty condition variable
 * signaling that all events from the task_queue have been processed. */
static void wait_for_sequential_section(struct volume *wrld) {
  pthread_mutex_lock(& wrld->dispatch_lock);
  while (wrld->task_queue.blocked_head != NULL ||
         wrld->task_queue.ready_head != NULL ||
         wrld->task_queue.num_pending != 0) {
    pthread_cond_wait(& wrld->dispatch_empty, & wrld->dispatch_lock);
  }
}


/* perform_sequential_section is responsible for executing the following
 * events in serial:
 * 1) molecule transfers between memory partitions
 * 2) reaction triggered molecule releases
 */
static int perform_sequential_section(struct volume *wrld) {
  for (int i=0; i<wrld->num_threads; ++i) {
    outbound_molecules_play(wrld, & wrld->threads[i].outbound);
  }

  return 1;
}


/* perform_delayed_sequential_actions executes count and trigger events in
 * a sequential fashion */
static int perform_delayed_sequential_actions(struct volume *wrld) {
  /* Update counts. */
  for (int i=0; i<wrld->num_threads; ++i) {
    delayed_count_play(& wrld->threads[i].count_updates);
  }

  /* Flush triggers. */
  for (int i=0; i<wrld->num_threads; ++i) {
    delayed_trigger_flush(wrld, & wrld->threads[i].triggers, 1);
  }

  return 1;
}


/* transfer_to_queue transfers the provided storage to the requested task
 * queue. */
static void transfer_to_queue(struct storage *store, struct storage **head) {
  /* Unlink it. */
  * (store->pprev) = store->next;
  if (store->next) {
    store->next->pprev = store->pprev;
  }

  /* Put it in the given queue. */
  store->next = *head;
  if (store->next != NULL) {
    store->next->pprev = & store->next;
  }
  store->pprev = head;
  *head = store;
}


static void unblock_neighbors(struct volume *wrld, thread_state_t *state,
  struct storage *last) {
  int xmin, xmax;
  int ymin, ymax;
  int zmin, zmax;

  /* Find the range over which to iterate (relative to 'last'. */
  xmin = - min2i(2, last->subdiv_x);
  ymin = - min2i(2, last->subdiv_y);
  zmin = - min2i(2, last->subdiv_z);
  xmax = wrld->subdivisions_nx - max2i((wrld->subdivisions_nx - 3), last->subdiv_x);
  ymax = wrld->subdivisions_ny - max2i((wrld->subdivisions_ny - 3), last->subdiv_y);
  zmax = wrld->subdivisions_nz - max2i((wrld->subdivisions_nz - 3), last->subdiv_z);

  /* Compute the y stride and z stride. */
  int ystride = wrld->subdivision_ystride - (xmax - xmin);
  int zstride = wrld->subdivision_zstride - (ymax - ymin) * wrld->subdivision_ystride;

  /* Find the starting storage. */
  struct storage *corner = last + zmin * wrld->subdivision_zstride
    + ymin * wrld->subdivision_ystride + xmin;

  /* Now, unlock and check each neighbor. */
  struct storage *store = corner;
  for (int z = zmin; z < zmax; ++ z) {
    for (int y = ymin; y < ymax; ++ y) {
      for (int x = xmin; x < xmax; ++ x) {
        /* Unlock. */
        if (-- store->lock_count == 0) {
          /* We were the last locker for this region.  Since no lockers remain,
           * toss the subdivision either into 'complete' or into 'ready',
           * according to whether work remains to be done. */
          if (is_subdivision_complete(store)) {
            // PARALLELDEBUG: */ mcell_log("Requeueing subdiv %p as COMPLETE (ubn).", store);
            transfer_to_queue(store, & wrld->task_queue.complete_head);
          } else {
            // PARALLELDEBUG: */ mcell_log("Requeueing subdiv %p as READY (ubn).", store);
            transfer_to_queue(store, & wrld->task_queue.ready_head);
          }
        }
        /* Advance to next store. */
        store += 1;
      }

      /* Advance to next row of stores. */
      store += ystride;
    }

    /* Advance to next slab of stores. */
    store += zstride;
  }
}


static void block_neighbors(struct volume *world, thread_state_t *state,
  struct storage *last) {
  int xmin, xmax;
  int ymin, ymax;
  int zmin, zmax;

  /* Find the range over which to iterate (relative to 'last'. */
  xmin = - min2i(2, last->subdiv_x);
  ymin = - min2i(2, last->subdiv_y);
  zmin = - min2i(2, last->subdiv_z);
  xmax = world->subdivisions_nx - max2i((world->subdivisions_nx - 3), last->subdiv_x);
  ymax = world->subdivisions_ny - max2i((world->subdivisions_ny - 3), last->subdiv_y);
  zmax = world->subdivisions_nz - max2i((world->subdivisions_nz - 3), last->subdiv_z);

  /* Compute the y stride and z stride. */
  int ystride = world->subdivision_ystride - (xmax - xmin);
  int zstride = world->subdivision_zstride - (ymax - ymin) * world->subdivision_ystride;

  /* Find the starting storage. */
  struct storage *corner = last
        + zmin * world->subdivision_zstride
        + ymin * world->subdivision_ystride
        + xmin;

  /* Now, lock each neighbor. */
  struct storage *store = corner;
  for (int z = zmin; z < zmax; ++ z)
  {
    for (int y = ymin; y < ymax; ++ y)
    {
      for (int x = xmin; x < xmax; ++ x)
      {
        /* Lock. */
        if (++ store->lock_count == 1)
        {
          /* We were the first locker for this region.  Toss the subdivision
           * into 'blocked'.
           */
          transfer_to_queue(store, & world->task_queue.blocked_head);
          // PARALLELDEBUG: */ mcell_log("Requeueing subdiv %p as BLOCKED (bn).", store);
        }

        /* Advance to next store. */
        store += 1;
      }

      /* Advance to next row of stores. */
      store += ystride;
    }

    /* Advance to next slab of stores. */
    store += zstride;
  }
}


static struct storage *schedule_subdivision(struct volume *wrld,
  thread_state_t *state, struct storage *last) {

  struct storage *subdiv = NULL;

  /* Acquire scheduler lock. */
  pthread_mutex_lock(& wrld->dispatch_lock);
  assert(wrld->sequential == 0);
  /* If we have a "last" subdivision from the previous schedule... */
  if (last != NULL) {
    // PARALLELDEBUG: thread_log("released task!");

    /* Remove "last" from the pending count. */
    --wrld->task_queue.num_pending;

    /* Transition newly unblocked subdivisions to 'ready' queue. */
    unblock_neighbors(wrld, state, last);

    /* If all subdivisions are in the completed queue... */
    if (wrld->task_queue.blocked_head == NULL && wrld->task_queue.ready_head == NULL &&
        wrld->task_queue.num_pending == 0) {
      // PARALLELDEBUG: thread_log("no tasks.  waking master.");

      /*   Wake the master. */
      pthread_cond_signal(&wrld->dispatch_empty);

      /*   Sleep until tasks are available. */
      pthread_cond_wait(&wrld->dispatch_ready, &wrld->dispatch_lock);
    }
  }

  /* Until we get a non-NULL subdivision... */
  while (subdiv == NULL) {

    /* subdiv <- Grab next subdivision. */
    subdiv = wrld->task_queue.ready_head;
    if (subdiv != NULL) {
      // PARALLELDEBUG: thread_log("got task!");
      wrld->task_queue.ready_head = subdiv->next;

      /* When we take the task, we transfer it to the blocked list until it
       * finishes executing. */
      transfer_to_queue(subdiv, & wrld->task_queue.blocked_head);
      // PARALLELDEBUG: */ mcell_log("Requeueing subdiv %p as BLOCKED (ssd).", subdiv);
    } else { /* if (subdiv == NULL) */
    /* else if no subdivision is available, we'll have to wait. */

      // PARALLELDEBUG: */ thread_log("no tasks.  sleeping.");
      /* wait for dispatch ready. */
      pthread_cond_wait(&wrld->dispatch_ready, &wrld->dispatch_lock);
    }
  }

  /* Move newly unavailable subdivisions to 'blocked' queue. */
  block_neighbors(wrld, state, subdiv);
  ++wrld->task_queue.num_pending;

  /* Release scheduler lock. */
  pthread_mutex_unlock(& wrld->dispatch_lock);

  return subdiv;
}


static void *worker_loop(struct worker_data *data) { //thread_state_t *state) {
  // PARALLELDEBUG: char filename[1024];
  // PARALLELDEBUG: snprintf(filename, 1024, "/tmp/thdlog.%08lx", pthread_self());
  // PARALLELDEBUG: FILE *outfile = fopen(filename, "w");
  // PARALLELDEBUG: setlinebuf(outfile);
  // PARALLELDEBUG: fprintf(outfile, "Worker beginning.\n");

  /* Stash our state. */
  struct volume *world = data->global_state;
  pthread_setspecific(world->thread_data, (void *)data->thread_state);

  // PARALLELDEBUG: fprintf(outfile, "Entering task loop.\n");

  /* Worker loop doesn't exit until the process exits. */
  struct storage *current = NULL;
  while (1) {
    // PARALLELDEBUG: fprintf(outfile, "Waiting for some work.\n");

    /* Return our current subdivision, if any, and grab the next scheduled
     * subdivision. */
    current = schedule_subdivision(world, data->thread_state, current);

    // PARALLELDEBUG: fprintf(outfile, "Running subdivision %p.\n", current);

    /* Play out the remainder of the iteration in this subdivision. */
    run_timestep(world, current, world->next_barrier, (double) (world->iterations + 1));
  }

  return NULL;
}


static void start_worker_pool(struct volume *wrld) {

  int num_workers = wrld->num_threads;

  /* Initialize sync primitives. */
  /* create or initialize thread-local variables */
  pthread_key_create(& wrld->thread_data, NULL);
  pthread_mutex_init(& wrld->trig_lock, NULL);
  pthread_cond_init(& wrld->dispatch_empty, NULL);
  pthread_cond_init(& wrld->dispatch_ready, NULL);

  /* Initialize thread states. */
  wrld->threads = CHECKED_MALLOC_ARRAY(thread_state_t, num_workers,
    "state structures for worker thread pool");
  for (int i=0; i<num_workers; ++i) {
    delayed_count_init(& wrld->threads[i].count_updates, 4096);
    delayed_trigger_init(& wrld->threads[i].triggers, 65536);
    outbound_molecules_init(& wrld->threads[i].outbound);
    struct worker_data data = {.global_state = wrld, .thread_state = &wrld->threads[i]};
    pthread_create(& wrld->threads[i].thread_id, NULL, (void *(*)(void *)) worker_loop,
                   (void *)&data);
  }
}


static void queue_subdivisions(struct volume *world) {
  world->task_queue.ready_head    = NULL;
  world->task_queue.complete_head = NULL;
  world->task_queue.blocked_head  = NULL;
  world->task_queue.num_pending   = 0;

  struct storage *subdiv = world->subdivisions;
  for (int i=0; i<world->num_subdivisions; ++i) {
    if (is_subdivision_complete(subdiv)) {
      // PARALLELDEBUG: */ mcell_log("Queueing subdiv %p as COMPLETE.", subdiv);
      if (world->task_queue.complete_head) {
        world->task_queue.complete_head->pprev = & subdiv->next;
      }
      subdiv->next = world->task_queue.complete_head;
      subdiv->pprev = & world->task_queue.complete_head;
      world->task_queue.complete_head = subdiv;
    } else {
      // PARALLELDEBUG: */ mcell_log("Queueing subdiv %p as READY.", subdiv);
      if (world->task_queue.ready_head) {
        world->task_queue.ready_head->pprev = & subdiv->next;
      }
      subdiv->next = world->task_queue.ready_head;
      subdiv->pprev = & world->task_queue.ready_head;
      world->task_queue.ready_head = subdiv;
    }
    ++ subdiv;
  }
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
  struct timing_info timing = {{ 0, 0 }, 0, world->current_iterations % frequency};

  /* Whether we are running in parallel or not, we begin in a sequential
   * section. */
  world->sequential = 1;

  /* Start up worker thread pool. */
  if (world->num_threads > 0) {
    pthread_mutex_init(&world->dispatch_lock, NULL);
    pthread_mutex_lock(&world->dispatch_lock);
    start_worker_pool(world);
    queue_subdivisions(world);
  }

  while (world->current_iterations <= world->iterations) {
    // XXX: A return status of 1 from mcell_run_iterations does not
    // indicate an error but is used to break out of the loop.
    // This behavior is non-conformant and should be changed.
    // NOTE: mcell_run_iteration requires the dispatch_lock mutex to be locked
    // upon entry.
    if (mcell_run_iteration(world, frequency, &timing, &restarted_from_checkpoint) == 1) {
      break;
    }
  }

  int status = 0;
  if (mcell_flush_data(world)) {
    mcell_error_nodie("Failed to flush reaction and visualization data.");
    status = 1;
  }

  // SMP: this may need to be updated to deal with the SMP code.
  if (mcell_print_final_warnings(world)) {
    mcell_error_nodie("Failed to print final warnings.");
    status = 1;
  }

  display_final_summary(world);

  return status;
}

/**************************************************************************
 *
 * this function runs a single mcell iteration
 *
 * parameters:
 *    frequency                : the requested output frequency
 *    restart_from_checkpoint  : does this iteration follow a checkpoint
 *
 *
 * PTHREAD NOTE: mcell_run_iteration requires the dispatch_lock mutex to be
 *               locked upon entry.
 *************************************************************************/
MCELL_STATUS
mcell_run_iteration(MCELL_STATE *world, long long frequency,
  struct timing_info *timing, int *restarted_from_checkpoint) {

  emergency_output_hook_enabled = 1;

  //long long iter_report_phase = world->current_iterations % frequency;
  long long not_yet = world->current_iterations + 1.0;

  if (world->current_iterations != 0) {
    world->elapsed_time = world->current_iterations;
  } else {
    world->elapsed_time = 1.0;
  }

  if (!*restarted_from_checkpoint) {

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

    produce_iteration_report(world, timing, frequency);

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
      /* Make a checkpoint, exiting the loop if necessary */
      if (make_checkpoint(world))
        return MCELL_FAIL;
    }

    /* Even if no checkpoint, the last iteration is a half-iteration. */
    if (world->current_iterations >= world->iterations)
      return MCELL_FAIL;
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

  /* Run in parallel mode */
  if (world->num_threads > 0) {
    /* Save next barrier for all workers. */
    world->next_barrier = next_barrier;

    /* While work remains */
    while (! is_iteration_complete(world)) {
      /* Begin non-sequential section */
      world->sequential = 0;
      wake_worker_pool(world);
      wait_for_sequential_section(world);

      /* End non-sequential section */
      world->sequential = 1;
      if (! perform_sequential_section(world)) {
        mcell_internal_error("Error while performing sequential section.");
        return MCELL_FAIL;
       }
     }

    if (! perform_delayed_sequential_actions(world)) {
      mcell_internal_error("Error while performing delayed sequential actions.");
      return MCELL_FAIL;
    }

    /* Advance the time. */
    for (int i=0; i<world->num_subdivisions; ++i) {
      struct storage *local = &world->subdivisions[i];
       /* Not using the return value -- just trying to advance the scheduler */
      void *o = schedule_next(local->timer);
      if (o != NULL) {
        mcell_internal_error("Scheduler dropped a molecule on the floor!");
      }
      local->current_time += 1.0;
      if (!is_subdivision_complete(local)) {
        transfer_to_queue(local, & world->task_queue.ready_head);
      }
    }
  } else  { /* Run in non-parallel mode */
    while (world->subdivisions[0].current_time <= not_yet) {
      int done = 0;
      while (! done) {
        done = 1;
        for (int i=0; i<world->num_subdivisions; ++i) {
          struct storage *local = & world->subdivisions[i];
          if (local->timer->current != NULL) {
            run_timestep(world, local, next_barrier, (double) (world->iterations + 1));
            done = 0;
          }
        }
      }

      for (int i=0; i<world->num_subdivisions; ++i) {
        struct storage *local = & world->subdivisions[i];
        /* Not using the return value -- just trying to advance the scheduler */
        void *o = schedule_next(local->timer);
        if (o != NULL)
          mcell_internal_error("Scheduler dropped a molecule on the floor!");
        local->current_time += 1.0;
      }
    }
  }

  world->current_iterations++;
  return MCELL_SUCCESS;
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

