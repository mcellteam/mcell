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

static int print_molecule_collision_report()
{
  if (world->notify->molecule_collision_report == NOTIFY_FULL)
  {
     mcell_log_raw("\n");
     mcell_log("\tCounts of Reaction Triggered Molecule Collisions");
     mcell_log("(VM = volume molecule, SM = surface molecule, W = wall)");
     if(world->mol_mol_reaction_flag) 
     {
       mcell_log("Total number of VM-VM collisions: %lld", world->mol_mol_colls);
     }
     if(world->mol_grid_reaction_flag)
     { 
        mcell_log("Total number of VM-SM collisions: %lld", world->mol_grid_colls);
     }
     if(world->grid_grid_reaction_flag)
     { 
       mcell_log("Total number of SM-SM collisions: %lld", world->grid_grid_colls);
     }
     if(world->mol_wall_reaction_flag)
     { 
        mcell_log("Total number of VM-W collisions: %lld", world->mol_wall_colls);
     }
     if(world->mol_mol_mol_reaction_flag)
     { 
        mcell_log("Total number of VM-VM-VM collisions: %lld", world->mol_mol_mol_colls); 
     }
     if(world->mol_mol_grid_reaction_flag)
     { 
       mcell_log("Total number of VM-VM-SM collisions: %lld", world->mol_mol_grid_colls);
     } 
     if(world->mol_grid_grid_reaction_flag)
     { 
       mcell_log("Total number of VM-SM-SM collisions: %lld", world->mol_grid_grid_colls); 
     }
     if(world->grid_grid_grid_reaction_flag)
     { 
       mcell_log("Total number of SM-SM-SM collisions: %lld", world->grid_grid_grid_colls);
     } 
     mcell_log_raw("\n");
  }

  return 0;
}


/***********************************************************************
 run_sim:

    Simulation main loop.

    In:  None
    Out: None
 ***********************************************************************/
static void run_sim(void)
{
  struct rusage run_time;
  time_t t_end;  /* global end time of MCell run */
  double u_init_time, s_init_time; /* initialization time (user) and (system) */
  double u_run_time, s_run_time; /* run time (user) and (system) */

  struct storage_list *local;
  double next_release_time, next_viz_output, next_vol_output;
  int first_report;
  /* used to suppress printing some warning messages when the reactant is a surface */
  int do_not_print;
  long long not_yet;
  long long frequency;
  
  emergency_output_hook_enabled = 1;
  if (world->notify->progress_report!=NOTIFY_NONE)
    mcell_log("Running simulation.");
  
  if (world->notify->custom_iteration_value != 0)
  {
    frequency = world->notify->custom_iteration_value;
  }
  else
  {
    if      (world->iterations < 10)         frequency = 1;
    else if (world->iterations < 1000)       frequency = 10;
    else if (world->iterations < 100000)     frequency = 100;
    else if (world->iterations < 10000000)   frequency = 1000;
    else if (world->iterations < 1000000000) frequency = 10000;
    else                                     frequency = 100000;
  }
  
  world->diffusion_number = 0;
  world->diffusion_cumtime = 0.0;
  world->it_time = world->start_time;
  world->last_checkpoint_iteration = 0;

  struct timeval last_timing_time = { 0, 0 };
  long long last_timing_iteration = 0;
  long long iter_report_phase = world->it_time % frequency;

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
    if ( iter_report_phase == 0 && world->notify->iteration_report != NOTIFY_NONE)
    {
      mcell_log_raw("Iterations: %lld of %lld ", world->it_time,world->iterations);

      if (world->notify->throughput_report != NOTIFY_NONE)
      {
        struct timeval cur_time;
        gettimeofday(&cur_time, NULL);
        if (last_timing_time.tv_sec > 0)
        {
          double time_diff = (double) (cur_time.tv_sec - last_timing_time.tv_sec) * 1000000.0 +
                (double) (cur_time.tv_usec - last_timing_time.tv_usec);
          time_diff /= (double)(world->it_time - last_timing_iteration);
          mcell_log_raw(" (%.6lg iter/sec)", 1000000.0 / time_diff);
          last_timing_iteration = world->it_time;
          last_timing_time = cur_time;
        }
        else
        {
          last_timing_iteration = world->it_time;
          last_timing_time = cur_time;
        }
      }

      mcell_log_raw("\n");
    }

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
    
    while (world->storage_head->store->current_time <= not_yet)
    {
      int done = 0;
      while (! done)
      {
        done = 1;
        for (local = world->storage_head ; local != NULL ; local = local->next)
        {
          if (local->store->timer->current != NULL)
          {
            run_timestep( local->store , next_barrier , (double)world->iterations+1.0 );
            done = 0;
          }
        }
      }

      for (local = world->storage_head ; local != NULL ; local = local->next)
      {
        /* Not using the return value -- just trying to advance the scheduler */
        void *o = schedule_next(local->store->timer);
        if (o != NULL)
          mcell_internal_error("Scheduler dropped a molecule on the floor!");
        local->store->current_time += 1.0;
      }
    }

    if (++ iter_report_phase == frequency) iter_report_phase = 0;
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
      if(world->species_list[i] == world->all_mols) continue;
      if ((world->species_list[i] == world->all_volume_mols)  ||
          (world->species_list[i] == world->all_surface_mols)) 
             continue;
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
   
  if(world->reaction_prob_limit_flag) mcell_log("Warning: During the simulation some reaction probabilities were greater than 1.  You may want to rerun the simulation with the WARNINGS block enabled to get more detail.\n");
 
  if (world->notify->final_summary==NOTIFY_FULL)
  {
    mcell_log("iterations = %lld ; elapsed time = %1.15g seconds",
              world->it_time,
              world->chkpt_elapsed_real_time_start+((world->it_time - world->start_time)*world->time_unit));

    if (world->diffusion_number > 0)
      mcell_log("Average diffusion jump was %.2f timesteps\n",
                world->diffusion_cumtime/(double)world->diffusion_number);
    mcell_log("Total number of random number use: %lld", rng_uses(world->rng));
    mcell_log("Total number of ray-subvolume intersection tests: %lld", world->ray_voxel_tests);
    mcell_log("Total number of ray-polygon intersection tests: %lld", world->ray_polygon_tests);
    mcell_log("Total number of ray-polygon intersections: %lld", world->ray_polygon_colls);
    print_molecule_collision_report();
 

    u_init_time = world->u_init_time.tv_sec + (world->u_init_time.tv_usec/MAX_TARGET_TIMESTEP);
    s_init_time = world->s_init_time.tv_sec + (world->s_init_time.tv_usec/MAX_TARGET_TIMESTEP); 

    mcell_log("Initialization CPU time = %f (user) and %f (system)", u_init_time, s_init_time);
               
    getrusage(RUSAGE_SELF,&run_time);
    u_run_time = run_time.ru_utime.tv_sec + (run_time.ru_utime.tv_usec/MAX_TARGET_TIMESTEP);
    s_run_time = run_time.ru_stime.tv_sec + (run_time.ru_stime.tv_usec/MAX_TARGET_TIMESTEP);


    mcell_log("Simulation CPU time = %f (user) and %f (system)", u_run_time - u_init_time, s_run_time - s_init_time);
    t_end = time(NULL);
    mcell_log("Total wall clock time = %ld seconds",
              (long)difftime(t_end, world->t_start) );
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

int main(int argc, char **argv)
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

  run_sim();

  if (world->notify->progress_report!=NOTIFY_NONE)
    mcell_log("Done running.");
  mem_dump_stats(mcell_get_log_file());

  exit(0);
}
