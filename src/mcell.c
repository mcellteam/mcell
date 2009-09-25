#if defined(__linux__)
#define _GNU_SOURCE 1
#endif

#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <signal.h>
#if defined(__linux__)
#include <fenv.h>
#endif
#include <float.h>

#include "sym_table.h"
#include "logging.h"
#include "rng.h"
#include "strfunc.h"
#include "argparse.h"
#include "vol_util.h"
#include "react_output.h"
#include "viz_output.h"
#include "volume_output.h"
#include "diffuse.h"
#include "init.h"
#include "chkpt.h"
#include "version_info.h"
#include "argparse.h"

struct volume *world;

/***********************************************************************
 process_volume_output:

    Produce this round's volume output, if any.

    In:  struct volume *wrld - the world
         double not_yet - earliest time which should not yet be output
    Out: none.  volume output files are updated as appropriate.
 ***********************************************************************/
static void process_volume_output(struct volume *wrld, double not_yet)
{
  struct volume_output_item *vo;
  for (vo = (struct volume_output_item *) schedule_next(wrld->volume_output_scheduler);
       vo != NULL  ||  not_yet >= wrld->volume_output_scheduler->now;
       vo = (struct volume_output_item *) schedule_next(wrld->volume_output_scheduler))
  {
    if (vo == NULL) continue;
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
static void process_reaction_output(struct volume *wrld, double not_yet)
{
  struct output_block *obp;
  for ( obp=schedule_next(wrld->count_scheduler) ;
        obp!=NULL || not_yet>=wrld->count_scheduler->now ;
        obp=schedule_next(wrld->count_scheduler) )
  {
    if (obp==NULL) continue;
    if (update_reaction_output(obp))
      mcell_error("Failed to update reaction output.");
  }
  if (wrld->count_scheduler->error)
    mcell_internal_error("Scheduler reported an out-of-memory error while retrieving next scheduled reaction output, but this should never happen.");
}

/***********************************************************************
 process_molecule_releases:

    Produce this round's release events, if any.

 In: wrld: the world
     not_yet: earliest time which should not yet be processed
 Out: none.  molecules are released into the world.
 ***********************************************************************/
static void process_molecule_releases(struct volume *wrld, double not_yet)
{
  for (struct release_event_queue *req = schedule_next(wrld->releaser);
       req != NULL || not_yet >= wrld->releaser->now;
       req = schedule_next(wrld->releaser)) 
  {
    if (req == NULL || req->release_site->release_prob==MAGIC_PATTERN_PROBABILITY) continue;
    if (release_molecules(req))
      mcell_error("Failed to release molecules of type '%s'.",
                  req->release_site->mol_type->sym->name);
  }
  if (wrld->releaser->error)
    mcell_internal_error("Scheduler reported an out-of-memory error while retrieving next scheduled release event, but this should never happen.");
}

/***********************************************************************
 make_checkpoint:

    Produce a checkpoint file.

    In:  struct volume *wrld - the world
    Out: 0 on success, 1 on failure.
         On success, checkpoint file is created.
         On failure, old checkpoint file, if any, is left intact.
 ***********************************************************************/
static int make_checkpoint(struct volume *wrld)
{
  /* Make sure we have a filename */
  if (wrld->chkpt_outfile == NULL)
    wrld->chkpt_outfile = CHECKED_SPRINTF("checkpt.%d", getpid());

  /* Print a useful status message */
  switch (wrld->checkpoint_requested)
  {
    case CHKPT_ITERATIONS_CONT:
    case CHKPT_ALARM_CONT:
      if (wrld->notify->checkpoint_report != NOTIFY_NONE)
        mcell_log("MCell: time = %lld, writing to checkpoint file %s (periodic).",
                  wrld->it_time,
                  wrld->chkpt_outfile);
      break;

    case CHKPT_ALARM_EXIT:
      if (wrld->notify->checkpoint_report != NOTIFY_NONE)
        mcell_log("MCell: time = %lld, writing to checkpoint file %s (time limit elapsed).",
                  wrld->it_time,
                  wrld->chkpt_outfile);
      break;

    case CHKPT_SIGNAL_CONT:
    case CHKPT_SIGNAL_EXIT:
      if (wrld->notify->checkpoint_report != NOTIFY_NONE)
        mcell_log("MCell: time = %lld, writing to checkpoint file %s (user signal detected).",
                  wrld->it_time,
                  wrld->chkpt_outfile);
      break;

    case CHKPT_ITERATIONS_EXIT:
      if (wrld->notify->checkpoint_report != NOTIFY_NONE)
        mcell_log("MCell: time = %lld, writing to checkpoint file %s.",
                  wrld->it_time,
                  wrld->chkpt_outfile);
      break;

    default:
      wrld->checkpoint_requested = CHKPT_NOT_REQUESTED;
      return 0;
  }

  /* Make the checkpoint */
  create_chkpt(wrld->chkpt_outfile);
  wrld->last_checkpoint_iteration = wrld->it_time;

  /* Break out of the loop, if appropriate */
  if (wrld->checkpoint_requested == CHKPT_ALARM_EXIT   ||
      wrld->checkpoint_requested == CHKPT_SIGNAL_EXIT  ||
      wrld->checkpoint_requested == CHKPT_ITERATIONS_EXIT)
    return 1;

  /* Schedule the next checkpoint, if appropriate */
  if (wrld->checkpoint_requested == CHKPT_ALARM_CONT)
  {
    if (wrld->continue_after_checkpoint)
      alarm(wrld->checkpoint_alarm_time);
    else
      return 1;
  }

  wrld->checkpoint_requested = CHKPT_NOT_REQUESTED;
  return 0;
}

static double find_next_viz_output_frame(struct frame_data_list *fdl)
{
  double next_time = DBL_MAX;
  for (; fdl != NULL; fdl = fdl->next)
  {
    if (fdl->curr_viz_iteration == NULL)
      continue;

    if (fdl->viz_iteration < next_time)
      next_time = fdl->viz_iteration;
  }

  return next_time;
}

static double find_next_viz_output(struct viz_output_block *vizblk)
{
  double next_time = DBL_MAX;
  while (vizblk != NULL)
  {
    double this_time = find_next_viz_output_frame(vizblk->frame_data_head);
    if (this_time < next_time)
      next_time = this_time;
    vizblk = vizblk->next;
  }

  return next_time;
}

/*
 * Get the log frequency specified in this world.
 */
static long long get_log_frequency(struct volume *wrld)
{
  if (wrld->notify->custom_iteration_value != 0)
    return wrld->notify->custom_iteration_value;

  if (wrld->iterations < 10)              return 1;
  else if (wrld->iterations < 1000)       return 10;
  else if (wrld->iterations < 100000)     return 100;
  else if (wrld->iterations < 10000000)   return 1000;
  else if (wrld->iterations < 1000000000) return 10000;
  else                                    return 100000;
}

struct timing_info
{
  struct timeval last_timing_time;
  long long last_timing_iteration;
  long long iter_report_phase;
};

/*
 * Produce the iteration/timing report.
 */
static void produce_iteration_report(struct volume *wrld,
                                     struct timing_info *timing,
                                     long long frequency)
{
  if (wrld->notify->iteration_report == NOTIFY_NONE)
    return;

  if (timing->iter_report_phase == 0)
  {
    mcell_log_raw("Iterations: %lld of %lld ", wrld->it_time, wrld->iterations);

    if (wrld->notify->throughput_report != NOTIFY_NONE)
    {
      struct timeval cur_time;
      gettimeofday(&cur_time, NULL);
      if (timing->last_timing_time.tv_sec > 0)
      {
        double time_diff = (double) (cur_time.tv_sec - timing->last_timing_time.tv_sec) * 1000000.0 +
              (double) (cur_time.tv_usec - timing->last_timing_time.tv_usec);
        time_diff /= (double)(wrld->it_time - timing->last_timing_iteration);
        mcell_log_raw(" (%.6lg iter/sec)", 1000000.0 / time_diff);
        timing->last_timing_iteration = wrld->it_time;
        timing->last_timing_time = cur_time;
      }
      else
      {
        timing->last_timing_iteration = wrld->it_time;
        timing->last_timing_time = cur_time;
      }
    }

    mcell_log_raw("\n");
  }

  if (++ timing->iter_report_phase == frequency)
    timing->iter_report_phase = 0;
}

static void display_final_summary(double t_initial)
{
  struct rusage run_time;

  mcell_log("iterations = %lld ; elapsed time = %1.15g seconds",
            world->it_time,
            world->chkpt_elapsed_real_time_start + ((world->it_time - world->start_time)*world->time_unit));

  /* Initialize aggregate statistics. */
  world->subdivisions[0].stats.random_numbers = rng_uses(world->subdivisions[0].rng);
  runtime_statistics_t min_stats;
  runtime_statistics_t max_stats;
  runtime_statistics_t overall_stats;
  runtime_statistics_t num_zeros;
  overall_stats = max_stats = min_stats = world->subdivisions[0].stats;
  num_zeros.diffusion_number  = (world->subdivisions[0].stats.diffusion_number  == 0) ? 1 : 0;
  num_zeros.ray_voxel_tests   = (world->subdivisions[0].stats.ray_voxel_tests   == 0) ? 1 : 0;
  num_zeros.ray_polygon_tests = (world->subdivisions[0].stats.ray_polygon_tests == 0) ? 1 : 0;
  num_zeros.ray_polygon_colls = (world->subdivisions[0].stats.ray_polygon_colls == 0) ? 1 : 0;
  num_zeros.mol_mol_colls     = (world->subdivisions[0].stats.mol_mol_colls     == 0) ? 1 : 0;
  num_zeros.mol_mol_mol_colls = (world->subdivisions[0].stats.mol_mol_mol_colls == 0) ? 1 : 0;
  num_zeros.random_numbers    = (world->subdivisions[0].stats.random_numbers    == 0) ? 1 : 0;

  /* Compute aggregate statistics. */
#define UPDATE_OVERALL(s1, s2, fld) (s1).fld += (s2).fld
#define UPDATE_ZERO(s1, s2, fld)    do {                                    \
  if ((s2).fld == 0) ++ (s1).fld;                                           \
} while (0)
#define UPDATE_MIN(s1, s2, fld)     do {                                    \
  if ((s2).fld != 0  &&  (s1).fld > (s2).fld) (s1).fld = (s2).fld;          \
} while (0)
#define UPDATE_MAX(s1, s2, fld)     do {                                    \
  if ((s2).fld != 0  &&  (s1).fld < (s2).fld) (s1).fld = (s2).fld;          \
} while (0)
  for (int i=1; i<world->num_subdivisions; ++i)
  {
    struct storage *store = & world->subdivisions[i];
    runtime_statistics_t *substat = & store->stats;
    substat->random_numbers = rng_uses(store->rng);

    /* Update zeros. */
    UPDATE_ZERO(num_zeros, *substat, diffusion_number);
    UPDATE_ZERO(num_zeros, *substat, ray_voxel_tests);
    UPDATE_ZERO(num_zeros, *substat, ray_polygon_tests);
    UPDATE_ZERO(num_zeros, *substat, ray_polygon_colls);
    UPDATE_ZERO(num_zeros, *substat, mol_mol_colls);
    UPDATE_ZERO(num_zeros, *substat, mol_mol_mol_colls);
    UPDATE_ZERO(num_zeros, *substat, random_numbers);

    /* Update overall stats. */
    UPDATE_OVERALL(overall_stats, *substat, diffusion_number);
    UPDATE_OVERALL(overall_stats, *substat, diffusion_cumtime);
    UPDATE_OVERALL(overall_stats, *substat, ray_voxel_tests);
    UPDATE_OVERALL(overall_stats, *substat, ray_polygon_tests);
    UPDATE_OVERALL(overall_stats, *substat, ray_polygon_colls);
    UPDATE_OVERALL(overall_stats, *substat, mol_mol_colls);
    UPDATE_OVERALL(overall_stats, *substat, mol_mol_mol_colls);
    UPDATE_OVERALL(overall_stats, *substat, random_numbers);

    /* Update minimum stats. */
    if (substat->diffusion_number != 0 
        && (min_stats.diffusion_cumtime / (double) min_stats.diffusion_number) >
           (substat->diffusion_cumtime  / (double) substat->diffusion_number))
    {
      min_stats.diffusion_cumtime = substat->diffusion_cumtime;
      min_stats.diffusion_number  = substat->diffusion_number;
    }
    UPDATE_MIN(min_stats, *substat, ray_voxel_tests);
    UPDATE_MIN(min_stats, *substat, ray_polygon_tests);
    UPDATE_MIN(min_stats, *substat, ray_polygon_colls);
    UPDATE_MIN(min_stats, *substat, mol_mol_colls);
    UPDATE_MIN(min_stats, *substat, mol_mol_mol_colls);
    UPDATE_MIN(min_stats, *substat, random_numbers);

    /* Update maximum stats. */
    if (substat->diffusion_number != 0 
        && (max_stats.diffusion_cumtime / (double) max_stats.diffusion_number) <
           (substat->diffusion_cumtime  / (double) substat->diffusion_number))
    {
      max_stats.diffusion_cumtime = substat->diffusion_cumtime;
      max_stats.diffusion_number  = substat->diffusion_number;
    }
    UPDATE_MAX(max_stats, *substat, ray_voxel_tests);
    UPDATE_MAX(max_stats, *substat, ray_polygon_tests);
    UPDATE_MAX(max_stats, *substat, ray_polygon_colls);
    UPDATE_MAX(max_stats, *substat, mol_mol_colls);
    UPDATE_MAX(max_stats, *substat, mol_mol_mol_colls);
    UPDATE_MAX(max_stats, *substat, random_numbers);
  }
#undef UPDATE_ZERO
#undef UPDATE_OVERALL
#undef UPDATE_MIN
#undef UPDATE_MAX

  /* Compute means. */
#define COMPUTE_MEAN(st) (overall_stats.st / (double) (world->num_subdivisions - num_zeros.st))
  double mean_diffusion_jump     = overall_stats.diffusion_cumtime / (double) overall_stats.diffusion_number;
  double min_mean_diffusion_jump = min_stats.diffusion_cumtime     / (double) min_stats.diffusion_number;
  double max_mean_diffusion_jump = max_stats.diffusion_cumtime     / (double) max_stats.diffusion_number;
  double mean_ray_voxel_tests    = COMPUTE_MEAN(ray_voxel_tests);
  double mean_ray_polygon_tests  = COMPUTE_MEAN(ray_polygon_tests);
  double mean_ray_polygon_colls  = COMPUTE_MEAN(ray_polygon_colls);
  double mean_mol_mol_colls      = COMPUTE_MEAN(mol_mol_colls);
  double mean_mol_mol_mol_colls  = COMPUTE_MEAN(mol_mol_mol_colls);
  double mean_random_numbers     = COMPUTE_MEAN(random_numbers);
#undef COMPUTE_MEAN

  /* Compute stdevs. */
#define COMPUTE_STDEV(st) do {                                              \
  stdev_##st = 0;                                                           \
  int num_samples = world->num_subdivisions - num_zeros.st;                 \
  if (num_samples > 1) {                                                    \
    for (int i=0; i<world->num_subdivisions; ++i) {                         \
      double diff = world->subdivisions[i].stats.st - mean_##st;            \
      stdev_##st += diff*diff;                                              \
    }                                                                       \
    stdev_##st /= (double) (num_samples - 1);                               \
  }                                                                         \
} while (0)
  double stdev_ray_voxel_tests;
  double stdev_ray_polygon_tests;
  double stdev_ray_polygon_colls;
  double stdev_mol_mol_colls;
  double stdev_mol_mol_mol_colls;
  double stdev_random_numbers;
  COMPUTE_STDEV(ray_voxel_tests);
  COMPUTE_STDEV(ray_polygon_tests);
  COMPUTE_STDEV(ray_polygon_colls);
  COMPUTE_STDEV(mol_mol_colls);
  COMPUTE_STDEV(mol_mol_mol_colls);
  COMPUTE_STDEV(random_numbers);
#undef COMPUTE_STDEV

  if (overall_stats.diffusion_number > 0)
    mcell_log("Average diffusion jump was %.2f timesteps (per-subdiv: min=%.2f, max=%.2f, zero=%d)\n",
              mean_diffusion_jump,
              min_mean_diffusion_jump,
              max_mean_diffusion_jump,
              (int) num_zeros.diffusion_number);
  mcell_log("Total number of random number use: %lld (per-subdiv: total=%lld, mean=%.2f, stdev=%.2f, min=%lld, max=%lld, zero=%d)",
            overall_stats.random_numbers + rng_uses(world->rng_global),
            overall_stats.random_numbers,
            mean_random_numbers,
            stdev_random_numbers,
            min_stats.random_numbers,
            max_stats.random_numbers,
            (int) num_zeros.random_numbers);
  mcell_log("Total number of ray-subvolume intersection tests: %lld (per-subdiv: mean=%.2f, stdev=%.2f, min=%lld, max=%lld, zero=%d)",
            overall_stats.ray_voxel_tests,
            mean_ray_voxel_tests,
            stdev_ray_voxel_tests,
            min_stats.ray_voxel_tests,
            max_stats.ray_voxel_tests,
            (int) num_zeros.ray_voxel_tests);
  mcell_log("Total number of ray-polygon intersection tests: %lld (per-subdiv: mean=%.2f, stdev=%.2f, min=%lld, max=%lld, zero=%d)",
            overall_stats.ray_polygon_tests,
            mean_ray_polygon_tests,
            stdev_ray_polygon_tests,
            min_stats.ray_polygon_tests,
            max_stats.ray_polygon_tests,
            (int) num_zeros.ray_polygon_tests);
  mcell_log("Total number of ray-polygon intersections: %lld (per-subdiv: mean=%.2f, stdev=%.2f, min=%lld, max=%lld, zero=%d)",
            overall_stats.ray_polygon_colls,
            mean_ray_polygon_colls,
            stdev_ray_polygon_colls,
            min_stats.ray_polygon_colls,
            max_stats.ray_polygon_colls,
            (int) num_zeros.ray_polygon_colls);
  mcell_log("Total number of molecule-molecule collisions: %lld (per-subdiv: mean=%.2f, stdev=%.2f, min=%lld, max=%lld, zero=%d)",
            overall_stats.mol_mol_colls,
            mean_mol_mol_colls,
            stdev_mol_mol_colls,
            min_stats.mol_mol_colls,
            max_stats.mol_mol_colls,
            (int) num_zeros.mol_mol_colls);
  mcell_log("Total number of molecule-molecule-molecule collisions: %lld (per-subdiv: mean=%.2f, stdev=%.2f, min=%lld, max=%lld, zero=%d)",
            overall_stats.mol_mol_mol_colls,
            mean_mol_mol_mol_colls,
            stdev_mol_mol_mol_colls,
            min_stats.mol_mol_mol_colls,
            max_stats.mol_mol_mol_colls,
            (int) num_zeros.mol_mol_mol_colls);

  long t_final = time(NULL);
  getrusage(RUSAGE_SELF, & run_time);
  mcell_log("Total CPU time = %f (user) and %f (system)",
            run_time.ru_utime.tv_sec + (run_time.ru_utime.tv_usec/MAX_TARGET_TIMESTEP),
            run_time.ru_stime.tv_sec + (run_time.ru_stime.tv_usec/MAX_TARGET_TIMESTEP) );
  mcell_log("Total wall clock time = %d seconds",
            (int)(t_final - t_initial) );
}

static int is_subdivision_complete(struct storage *nation)
{
  if (nation->inbound != NULL  &&
      nation->inbound->fill != 0)
    return 0;

  if (nation->timer->current_count != 0)
    return 0;

  return 1;
}

static int is_iteration_complete(struct volume *wrld)
{
  if (task_queue.ready_head == NULL)
    return 0;

  if (task_queue.blocked_head == NULL)
    return 0;

  return 1;
}

static void wake_worker_pool(struct volume *wrld)
{
  pthread_mutex_lock(& wrld->dispatch_lock);
  pthread_cond_broadcast(& wrld->dispatch_ready);
  pthread_mutex_unlock(& wrld->dispatch_lock);
}

static void wait_for_sequential_section(struct volume *wrld)
{
  pthread_mutex_lock(& wrld->dispatch_lock);
  pthread_cond_wait(& wrld->dispatch_empty, & wrld->dispatch_lock);
  pthread_mutex_unlock(& wrld->dispatch_lock);
}

static int perform_sequential_section(struct volume *wrld)
{
  /* Perform molecule transfers. */
  /* Perform rxntrig releases. */
  for (int i=0; i<wrld->num_threads; ++i)
    outbound_molecules_play(wrld, & wrld->threads[i].outbound);

  return 1;
}

static int perform_delayed_sequential_actions(struct volume *wrld)
{
  /* Update counts. */
  for (int i=0; i<wrld->num_threads; ++i)
    delayed_count_play(& wrld->threads[i].count_updates);

  /* Flush triggers. */
  for (int i=0; i<wrld->num_threads; ++i)
    delayed_trigger_flush(& wrld->threads[i].triggers, 1);

  return 1;
}

static void transfer_to_queue(struct storage *store,
                              struct storage **head)
{
  /* Unlink it. */
  * (last->pprev) = last->next;
  if (last->next)
    last->next->pprev = last->pprev;

  /* Put it in the given queue. */
  last->next = *head;
  if (last->next != NULL)
    last->next->pprev = & last->next;
  last->pprev = head;
}

static void unblock_neighbors(struct volume *wrld,
                              thread_state_t *state,
                              struct storage *last)
{
  int xmin, xmax;
  int ymin, ymax;
  int zmin, zmax;

  /* Find the region over which to iterate. */
  xmin = - min(1, last->subdiv_x);
  ymin = - min(1, last->subdiv_y);
  zmin = - min(1, last->subdiv_z);
  xmax = wrld->subdivisions_nx - max((wrld->subdivisions_nx - 2), last->subdiv_x);
  ymax = wrld->subdivisions_ny - max((wrld->subdivisions_ny - 2), last->subdiv_y);
  zmax = wrld->subdivisions_nz - max((wrld->subdivisions_nz - 2), last->subdiv_z);

  /* Compute the y stride and z stride. */
  int ystride = wrld->subdivision_ystride - (xmax - xmin);
  int zstride = wrld->subdivision_zstride - (ymax - ymin) * wrld->subdivision_ystride;

  /* Find the starting storage. */
  struct storage *corner = world->subdivisions
        + zmin * world->subdivision_zstride
        + ymin * world->subdivision_ystride;

  /* Now, unlock and check each neighbor. */
  struct storage *store = corner;
  for (int z = zmin; z < zmax; ++ z)
  {
    for (int y = ymin; y < ymax; ++ y)
    {
      for (int x = xmin; x < xmax; ++ x)
      {
        /* Unlock. */
        store->locker = NULL;

        /* If no work remains to be done... */
        if (is_subdivision_complete(store))
          transfer_to_queue(store, & wrld->task_queue.complete_head);

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
                                            thread_state_t *state,
                                            struct storage *last)
{
  struct storage *subdiv = NULL;

  /* Acquire scheduler lock. */
  pthread_mutex_lock(& wrld->dispatch_lock);

  /* XXX */
  /* If we have a "last" subdivision from the previous schedule... */
  if (last != NULL)
  {
    /* Put "current" in the completed queue. */
    last->next = & wrld->task_queue.complete_head;
    if (last->next != NULL)
      last->next->pprev = & last->next;
    last->pprev = & wrld->task_queue.complete_head;
    -- wrld->task_queue.num_pending;

    /* Transition newly unblocked subdivisions to 'ready' queue. */
    unblock_neighbors(wrld, state, last);

    /* If all subdivisions are in the completed queue... */
    if (wrld->task_queue.blocked_head == NULL  &&
        wrld->task_queue.ready_head == NULL    &&
        wrld->task_queue.num_pending == 0)
    {
      /*   Wake the master. */
      pthread_cond_signal(& wrld->dispatch_empty);

      /*   Sleep until tasks are available. */
      pthread_cond_wait(& wrld->dispatch_ready,
                        & wrld->dispatch_lock);
    }
  }

  /* While subdiv is NULL */
  while (subdiv != NULL)
  {
    /* XXX */
    /* subdiv <- Grab next subdivision. */

    /* if no subdivision is available... */
    if (subdiv == NULL)
    {
      /* wait for dispatch ready. */
      pthread_cond_wait(& wrld->dispatch_ready,
                        & wrld->dispatch_lock);
    }
  }

  /* XXX */
  /* Move newly unavailable subdivisions to 'blocked' queue. */

  /* Release scheduler lock. */
  pthread_mutex_unlock(& wrld->dispatch_lock);

  return subdiv;
}

static void worker_loop(struct volume *wrld,
                        thread_state_t *state)
{
  /* Worker loop doesn't exit until the process exits. */
  struct storage *current = NULL;
  while (1)
  {
    /* Return our current subdivision, if any, and grab the next scheduled
     * subdivision. */
    current = schedule_subdivision(wrld, state, current);

    /* Play out the remainder of the iteration in this subdivision. */
    run_timestep(current, wrld->next_barrier, (double) (wrld->iterations + 1));
  }
}

/***********************************************************************
 run_sim:

    Simulation main loop.

    In:  None
    Out: None
 ***********************************************************************/
static void run_sim(void)
{
  long t_initial;

  double next_release_time, next_viz_output, next_vol_output;
  int first_report;
  /* used to suppress printing some warning messages when the reactant is a surface */
  int do_not_print;
  long long not_yet;
  long long frequency = get_log_frequency(world);
  
  emergency_output_hook_enabled = 1;
  if (world->notify->progress_report!=NOTIFY_NONE)
    mcell_log("Running simulation.");

  t_initial = time(NULL);

  world->it_time = world->start_time;
  world->last_checkpoint_iteration = 0;

  struct timing_info timing = {
    { 0, 0 },
    0,
    world->it_time % frequency
  };

  /* If we're reloading a checkpoint, we want to skip all of the processing
   * which happened on the last iteration before checkpointing.  To do this, we
   * skip the first part of the loop using a goto.  This is primarily relevant
   * to release events which occur on the final iteration.
   */
  if (world->start_time != 0)
  {
    not_yet = world->it_time + 1.0;
    /* TIME_STEP may have been changed between checkpoints */
    world->current_real_time = world->current_start_real_time + (world->it_time - world->start_time)*world->time_unit;
    world->elapsed_time = world->it_time;

    goto resume_after_checkpoint;
  }
 
  while (world->it_time <= world->iterations) 
  {
    not_yet = world->it_time + 1.0;

    if (world->it_time!=0) world->elapsed_time=world->it_time;
    else world->elapsed_time=1.0;
  
    /* Release molecules */
    process_molecule_releases(world, not_yet);

    /* Produce output */
    process_reaction_output(world, not_yet);
    process_volume_output(world, not_yet);
    for (struct viz_output_block *vizblk = world->viz_blocks;
         vizblk != NULL;
         vizblk = vizblk->next)
    {
      if (vizblk->frame_data_head  &&  update_frame_data_list(vizblk))
      mcell_error("Unknown error while updating frame data list.");
    }

    /* Produce iteration report */
    produce_iteration_report(world, & timing, frequency);

    /* Check for a checkpoint on this iteration */
    if (world->chkpt_iterations  &&  (world->it_time - world->start_time) == world->chkpt_iterations)
      world->checkpoint_requested = CHKPT_ITERATIONS_EXIT;
 
    /* No checkpoint signalled.  Keep going. */
    if (world->checkpoint_requested != CHKPT_NOT_REQUESTED)
    {
      /* Make a checkpoint, exiting the loop if necessary */
      if (make_checkpoint(world))
        break;
    }

    /* Even if no checkpoint, the last iteration is a half-iteration. */
    if (world->it_time >= world->iterations)
      break;

resume_after_checkpoint:    /* Resuming loop here avoids extraneous releases */
    
    run_concentration_clamp(world->it_time);

    if (! schedule_anticipate( world->releaser , &next_release_time))
      next_release_time = world->iterations + 1;
    if (next_release_time < world->it_time+1) next_release_time = world->it_time+1;
    if (! schedule_anticipate( world->volume_output_scheduler , &next_vol_output))
      next_vol_output = world->iterations + 1;
    next_viz_output = find_next_viz_output(world->viz_blocks);
    double next_barrier = min3d(next_release_time, next_vol_output, next_viz_output);

    /* Stuff happens. */
#if 1
    /* Save next barrier for all workers. */
    world->next_barrier = next_barrier;

    /* While work remains */
    while (! is_iteration_complete(world))
    {
      /* Wake Worker Pool */
      wake_worker_pool(world);

      /* Sleep until sequential section */
      wait_for_sequential_section(world);

      /* Perform sequential actions */
      if (perform_sequential_section(world))
      {
        mcell_internal_error("Error while performing sequential section.");
        return;
      }
    }

    /* Perform delayed sequential actions. */
    if (perform_delayed_sequential_actions(world))
    {
      mcell_internal_error("Error while performing delayed sequential actions.");
      return;
    }

#else
    while (world->subdivisions[0].current_time <= not_yet)
    {
      int done = 0;
      while (! done)
      {
        done = 1;
        for (int i=0; i<world->num_subdivisions; ++i)
        {
          struct storage *local = & world->subdivisions[i];
          if (local->store->timer->current != NULL)
          {
            run_timestep(local->store, next_barrier, (double) (world->iterations + 1));
            done = 0;
          }
        }
      }

      for (int i=0; i<world->num_subdivisions; ++i)
      {
        struct storage *local = & world->subdivisions[i];
        /* Not using the return value -- just trying to advance the scheduler */
        void *o = schedule_next(local->store->timer);
        if (o != NULL)
          mcell_internal_error("Scheduler dropped a molecule on the floor!");
        local->store->current_time += 1.0;
      }
    }
#endif

    world->it_time++;
  }
  
  /* If we didn't make a final iteration checkpoint, make one */
  if (world->chkpt_iterations  &&  world->it_time > world->last_checkpoint_iteration)
    make_checkpoint(world);
  
  emergency_output_hook_enabled = 0;
  int num_errors = flush_reaction_output();
  if (num_errors != 0)
  {
    mcell_warn("%d errors occurred while flushing buffered reaction output.\n"
               "  Simulation complete anyway--continuing as normal.",
               num_errors);
  }

  if (world->notify->progress_report!=NOTIFY_NONE)
    mcell_log("Exiting run loop.");
  int warned = 0;
  for (struct viz_output_block *vizblk = world->viz_blocks;
       vizblk != NULL;
       vizblk = vizblk->next)
  {
    if (finalize_viz_output(vizblk)  &&  ! warned)
    {
    mcell_warn("VIZ output was not successfully finalized.\n"
               "  Visualization of results may not work correctly.");
      warned = 1;
    }
  }
 
  first_report=1;
  
  if (world->notify->missed_reactions != WARN_COPE)
  {
    for (int i=0;i<world->rx_hashsize;i++)
    {
      for (struct rxn *rxp = world->reaction_hash[i]; rxp != NULL; rxp = rxp->next)
      {
        do_not_print = 0;
        /* skip printing messages if one of the reactants is a surface */
        for (unsigned int j=0;j<rxp->n_reactants;j++)
        {
          if ((rxp->players[j]->flags & IS_SURFACE) != 0)
            do_not_print = 1;
        }

        if (do_not_print == 1) continue;
        if (rxp->n_occurred*world->notify->missed_reaction_value < rxp->n_skipped)
        {
          if (first_report)
          {
            mcell_log("Warning: Some reactions were missed because reaction probability exceeded 1.");
            first_report=0; 
          }
          mcell_log_raw("  ");
          for (unsigned int j=0; j<rxp->n_reactants; j++)
          {
            mcell_log_raw("%s%s[%d]",
                          j ? " + " : "",
                          rxp->players[j]->sym->name,
                          rxp->geometries[j]);
          }
          mcell_log_raw("  --  %g%% of reactions missed.\n",
                        0.001*round(1000*rxp->n_skipped*100/(rxp->n_skipped+rxp->n_occurred)));
        }
      }
    }
    if (!first_report) mcell_log_raw("\n");
  }
  
  first_report+=1;
  
  if (world->notify->short_lifetime != WARN_COPE)
  {
    for (int i=0;i<world->n_species;i++)
    {
      if (world->species_list[i]->n_deceased <= 0)
        continue;

      double f = world->species_list[i]->cum_lifetime / world->species_list[i]->n_deceased;
      if (f < world->notify->short_lifetime_value)
      {
        if (first_report)
        {
          if (first_report>1)
            mcell_log_raw("\n");
          mcell_log("Warning: Some molecules had a lifetime short relative to the timestep.");
          first_report=0;
        }
        mcell_log("  Mean lifetime of %s was %g timesteps.",
                  world->species_list[i]->sym->name,
                  0.01*round(100*f));
      }
    }
    if (!first_report) mcell_log_raw("\n");
  }
    
  if (world->notify->final_summary==NOTIFY_FULL)
  {
    display_final_summary(t_initial);
  }
}

/***********************************************************************
 install_usr_signal_handlers:

    Set signal handlers for checkpointing on SIGUSR signals.

    In:  None
    Out: 0 on success, 1 on failure.
 ***********************************************************************/
static int install_usr_signal_handlers(void)
{
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

  return 0;
}

int mcell_main(int argc, char **argv)
{
  char hostname[64];
  u_int procnum;
  long long exec_iterations = 0; /* number of simulation iterations for this run */
  time_t begin_time_of_day;  /* start time of the simulation */

  mcell_set_log_file(stdout);
  mcell_set_error_file(stderr);

  /* get the process start time */
  time(&begin_time_of_day);
  if (install_usr_signal_handlers())
    mcell_die();

#if defined(__linux__)
  feenableexcept(FE_DIVBYZERO);
#endif


  world = CHECKED_MALLOC_STRUCT(struct volume, "world");
  memset(world, 0, sizeof(struct volume));

  world->procnum=0;
  procnum=world->procnum;
  gethostname(hostname,64);

  world->iterations=INT_MIN; /* indicates iterations not set */
  world->chkpt_infile = NULL;
  world->chkpt_outfile = NULL;
  world->chkpt_init = 1;
  world->log_freq = ULONG_MAX; /* Indicates that this value has not been set by user */
  world->begin_timestamp = begin_time_of_day;
  world->initialization_state = "initializing";

  if ((world->var_sym_table = init_symtab(1024)) == NULL)
    mcell_allocfailed("Failed to initialize MDL variable symbol table.");

  /*
   * Parse the command line arguments and print out errors if necessary.
   */
  if (argparse_init(argc,argv,world))
  {
    if (procnum == 0)
    {
      print_version(mcell_get_log_file());
      print_usage(mcell_get_log_file(), argv[0]);
    }

    exit(1);
  }

  if (init_notifications())
    mcell_error("Unknown error while initializing user-notification data structures.");
  
  if (world->notify->progress_report!=NOTIFY_NONE)
    print_version(mcell_get_log_file());

  if (init_sim())
    mcell_error("An unknown error occurred inside the MDL parser.\n             This was likely caused by an out-of-memory error.");
  world->initialization_state = NULL;

  if(world->chkpt_flag)
  {
    if (world->notify->checkpoint_report != NOTIFY_NONE)
      mcell_log("MCell: checkpoint sequence number %d begins at elapsed time %1.15g seconds",
                world->chkpt_seq_num,
                world->chkpt_elapsed_real_time_start);
    if (world->iterations < world->start_time)
    {
      mcell_error("Start time after checkpoint %lld is greater than total number of iterations specified %lld.",
                  world->start_time,
                  world->iterations);
    }
    if (world->chkpt_iterations)
    {
      if ((world->iterations - world->start_time) < world->chkpt_iterations)
        world->chkpt_iterations = world->iterations - world->start_time;
      else
        world->iterations = world->chkpt_iterations + world->start_time;
    }

    if (world->chkpt_iterations)
      exec_iterations = world->chkpt_iterations;
    else if (world->chkpt_infile)
      exec_iterations = world->iterations - world->start_time;
    else
      exec_iterations = world->iterations;
    if (exec_iterations < 0)
      mcell_error("Number of iterations to execute is zero or negative. Please verify ITERATIONS and/or CHECKPOINT_ITERATIONS commands.");
    if (world->notify->progress_report != NOTIFY_NONE)
      mcell_log("MCell: executing %lld iterations starting at iteration number %lld.",
                exec_iterations,
                world->start_time);
  }

  if((world->chkpt_flag) && (exec_iterations <= 0))
  {
    mem_dump_stats(mcell_get_log_file());
    exit(0);
  }

  if (world->notify->progress_report!=NOTIFY_NONE)
    mcell_log("Running...");
  run_sim();

  if (world->notify->progress_report!=NOTIFY_NONE)
    mcell_log("Done running.");
  mem_dump_stats(mcell_get_log_file());

  exit(0);
}
