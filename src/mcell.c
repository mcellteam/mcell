#if defined(__linux__)
#define _GNU_SOURCE 1
#endif

#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <signal.h>
#if defined(__linux__)
#include <fenv.h>
#endif

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
  int i;            /* for emergency output */
  struct volume_output_item *vo;
  for (vo = (struct volume_output_item *) schedule_next(world->volume_output_scheduler);
       vo != NULL  ||  not_yet >= world->volume_output_scheduler->now;
       vo = (struct volume_output_item *) schedule_next(world->volume_output_scheduler))
  {
    if (vo == NULL) continue;
    if (update_volume_output(world, vo))
    {
      fprintf(world->err_file,"File '%s', Line %ld: Error while updating volume output. Trying to save intermediate results.\n", __FILE__, (long)__LINE__);
      i = emergency_output();
      fprintf(world->err_file, "%d error%s while saving intermediate results.\n", i, (i==1) ? "": "s");
      exit(EXIT_FAILURE);
    }
  }
}

static void process_reaction_output(struct volume *wrld, double not_yet)
{
  int i;            /* for emergency output */
  struct output_block *obp;
  for ( obp=schedule_next(wrld->count_scheduler) ;
        obp!=NULL || not_yet>=wrld->count_scheduler->now ;
        obp=schedule_next(wrld->count_scheduler) )
  {
    if (obp==NULL) continue;
    if (update_reaction_output(obp))
    {
      fprintf(wrld->err_file,"File '%s', Line %ld: Error while updating reaction output. Trying to save intermediate results.\n", __FILE__, (long)__LINE__);
      i = emergency_output();
      fprintf(wrld->err_file,"%d error%s while saving intermediate results.\n",i,(i==1)?"":"s");
      exit(EXIT_FAILURE);
    }
  }
  if (wrld->count_scheduler->error)
  {
    fprintf(wrld->err_file,"File '%s', Line %ld: Out of memory while scheduling molecule release. Trying to save intermediate results.\n", __FILE__, (long)__LINE__);
    i = emergency_output();
    fprintf(wrld->err_file,"%d error%s while saving intermediate results.\n",i,(i==1)?"":"s");
    exit(EXIT_FAILURE);
  }
}

static void process_molecule_releases(struct volume *wrld, double not_yet)
{
  int i;            /* for emergency output */
  struct release_event_queue *req;
  for ( req= schedule_next(world->releaser) ;
        req!=NULL || not_yet>=world->releaser->now ;
        req=schedule_next(world->releaser)) 
  {
    if (req==NULL || req->release_site->release_prob==MAGIC_PATTERN_PROBABILITY) continue;
    if ( release_molecules(req) )
    {
      fprintf(world->err_file,"File '%s', Line %ld: Out of memory while releasing molecules of type %s\n", __FILE__, (long)__LINE__, req->release_site->mol_type->sym->name);
      i = emergency_output();
      fprintf(world->err_file,"%d error%s while saving intermediate results.\n",i,(i==1)?"":"s");
      exit(EXIT_FAILURE);
    }
  }
  if (world->releaser->error)
  {
    fprintf(world->err_file,"File '%s', Line %ld: Out of memory while scheduling molecule release. Trying to save intermediate results.\n", __FILE__, (long)__LINE__);
    i = emergency_output();
    fprintf(world->err_file,"%d error%s while saving intermediate results.\n",i,(i==1)?"":"s");
    exit(EXIT_FAILURE);
  }
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
    wrld->chkpt_outfile = alloc_sprintf("checkpt.%d", getpid());

  /* Print a useful status message */
  switch (wrld->checkpoint_requested)
  {
    case CHKPT_ITERATIONS_CONT:
    case CHKPT_ALARM_CONT:
      if (wrld->notify->progress_report != NOTIFY_NONE)
        fprintf(wrld->log_file,
                "MCell: time = %lld, writing to checkpoint file %s (periodic)\n",
                wrld->it_time,
                wrld->chkpt_outfile);
      break;

    case CHKPT_ALARM_EXIT:
      fprintf(wrld->log_file,
              "MCell: time = %lld, writing to checkpoint file %s (time limit elapsed)\n",
              wrld->it_time,
              wrld->chkpt_outfile);
      break;

    case CHKPT_SIGNAL_CONT:
    case CHKPT_SIGNAL_EXIT:
      fprintf(wrld->log_file,
              "MCell: time = %lld, writing to checkpoint file %s (user signal detected)\n",
              wrld->it_time,
              wrld->chkpt_outfile);
      break;

      if (wrld->notify->progress_report != NOTIFY_NONE)
        fprintf(wrld->log_file,
                "MCell: time = %lld, writing to checkpoint file %s (periodic)\n",
                wrld->it_time,
                wrld->chkpt_outfile);
      break;

    case CHKPT_ITERATIONS_EXIT:
      fprintf(wrld->log_file,
              "MCell: time = %lld, writing to checkpoint file %s\n",
              wrld->it_time,
              wrld->chkpt_outfile);
      break;

    default:
      wrld->checkpoint_requested = CHKPT_NOT_REQUESTED;
      return 0;
  }

  /* Make the checkpoint */
  if (create_chkpt(wrld->chkpt_outfile))
    exit(1);
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

/***********************************************************************
 run_sim:

    Simulation main loop.

    In:  None
    Out: None
 ***********************************************************************/
void run_sim(void)
{
  struct rusage run_time;
  long t_initial,t_final;

  struct storage_list *local;
  double next_release_time;
  int i;
  int first_report;
  /* used to suppress printing some warning messages when the reactant is a surface */
  int do_not_print;
  long long not_yet;
  long long frequency;
  
  if (world->notify->progress_report!=NOTIFY_NONE) fprintf(world->log_file,"Running simulation.\n");

  t_initial = time(NULL);
  
  if (world->notify->custom_iterations==NOTIFY_CUSTOM)
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
    else                                       frequency = 100000;
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
    if (world->frame_data_head  &&  update_frame_data_list(world->frame_data_head))
    {
       fprintf(world->err_file, "Error while updating frame data list.\n");
       exit(EXIT_FAILURE);
    } 

    /* Produce iteration report */
    if ( iter_report_phase == 0 && world->notify->custom_iterations!=NOTIFY_NONE)
    {
      printf("Iterations: %lld of %lld ",world->it_time,world->iterations);

      if (world->notify->throughput_report != NOTIFY_NONE)
      {
        struct timeval cur_time;
        gettimeofday(&cur_time, NULL);
        if (last_timing_time.tv_sec > 0)
        {
          double time_diff = (double) (cur_time.tv_sec - last_timing_time.tv_sec) * 1000000.0 +
                (double) (cur_time.tv_usec - last_timing_time.tv_usec);
          time_diff /= (double)(world->it_time - last_timing_iteration);
          printf(" (%.6lg iter/sec)", 1000000.0 / time_diff);
          last_timing_iteration = world->it_time;
          last_timing_time = cur_time;
        }
        else
        {
          last_timing_iteration = world->it_time;
          last_timing_time = cur_time;
        }
      }

      printf("\n");
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

resume_after_checkpoint:    /* Resuming loop here avoids extraneous releases */
    
    run_concentration_clamp(world->it_time);

    i = schedule_anticipate( world->releaser , &next_release_time);
    if (!i) next_release_time = world->iterations + 1;
    if (next_release_time < world->it_time+1) next_release_time = world->it_time+1;
    
    
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
            run_timestep( local->store , next_release_time , (double)world->iterations+1.0 );
            done = 0;
          }
        }
      }

      for (local = world->storage_head ; local != NULL ; local = local->next)
      {
        /* Not using the return value -- just trying to advance the scheduler */
        void *o = schedule_next(local->store->timer);
        if (o != NULL)
          fprintf(stderr, "Internal error!  Scheduler dropped a molecule on the floor!\n");
        local->store->current_time += 1.0;
      }
    }

    if (++ iter_report_phase == frequency) iter_report_phase = 0;
    world->it_time++;
  }
  
  /* If we didn't make a final iteration checkpoint, make one */
  if (world->chkpt_iterations  &&  world->it_time > world->last_checkpoint_iteration)
    make_checkpoint(world);
  
  i = flush_reaction_output();
  if (i)
  {
    fprintf(world->err_file,"Error at file %s line %d\n",__FILE__,__LINE__);
    fprintf(world->err_file,"  Could not write output for triggered reactions: %d errors\n",i);
    fprintf(world->err_file,"  Simulation complete anyway--continuing as normal.\n");
  }

  if (world->notify->progress_report!=NOTIFY_NONE) fprintf(world->log_file,"Exiting run loop.\n");
  if (finalize_viz_output(world->frame_data_head))
  {
    fprintf(world->err_file, "Warning: viz output was not successfully finalized.\n");
    fprintf(world->err_file, "  Visualization of results may not work correctly.\n");
  }
 
  first_report=1;
  
  if (world->notify->missed_reactions != WARN_COPE)
  {
    for (i=0;i<world->rx_hashsize;i++)
    {
      struct rxn *rxp;
      int j;
      
      for (rxp = world->reaction_hash[i] ; rxp != NULL ; rxp = rxp->next)
      {
         do_not_print = 0;
         /* skip printing messages if one of the reactants is a surface */
          for (j=0;j<rxp->n_reactants;j++)
          {
            if((rxp->players[j]->flags & IS_SURFACE) != 0) {
                  do_not_print = 1;
            }
          }
               
          if(do_not_print == 1) continue;
        if (rxp->n_occurred*world->notify->missed_reaction_value < rxp->n_skipped)
        {
          if (first_report)
          {
            fprintf(world->log_file,"\nWARNING: some reactions were missed because reaction probability exceeded 1.\n");
            first_report=0; 
          }
          fprintf(world->log_file,"  ");
          for (j=0;j<rxp->n_reactants;j++)
          {
            fprintf(world->log_file,"%s%s[%d]",(j)?" + ":"",rxp->players[j]->sym->name,rxp->geometries[j]);
          }
          fprintf(world->log_file,"  --  %g%% of reactions missed.\n",0.001*round(1000*rxp->n_skipped*100/(rxp->n_skipped+rxp->n_occurred)));
        }
      }
    }
    if (!first_report) fprintf(world->log_file,"\n");
  }
  
  first_report+=1;
  
  if (world->notify->short_lifetime != WARN_COPE)
  {
    for (i=0;i<world->n_species;i++)
    {
      double f;
      
      if (world->species_list[i]->n_deceased > 0)
      {
        f = world->species_list[i]->cum_lifetime / world->species_list[i]->n_deceased;
        
        if (f < world->notify->short_lifetime_value)
        {
          if (first_report)
          {
            if (first_report>1) fprintf(world->log_file,"\n");
            fprintf(world->log_file,"WARNING: some molecules had a lifetime short relative to the timestep.\n");
            first_report=0;
          }
          fprintf(world->log_file,"  Mean lifetime of %s was %g timesteps.\n",world->species_list[i]->sym->name,0.01*round(100*f));
        }
      }
    }
    if (!first_report) fprintf(world->log_file,"\n");
  }
    
  if (world->notify->final_summary==NOTIFY_FULL)
  {
    fprintf(world->log_file,"iterations = %lld ; elapsed time = %1.15g seconds\n",world->it_time,world->chkpt_elapsed_real_time_start+((world->it_time - world->start_time)*world->time_unit));
    fflush(world->log_file);

    if(world->diffusion_number > 0)
    { 
      fprintf(world->log_file,"Average diffusion jump was %.2f timesteps\n",world->diffusion_cumtime/(double)world->diffusion_number);
    }
    if(world->notify->final_summary == NOTIFY_FULL){
       fprintf(world->log_file,"Total number of random number use: %lld\n",world->random_number_use);
       fprintf(world->log_file,"Total number of ray-subvolume intersection tests: %lld\n",world->ray_voxel_tests);
       fprintf(world->log_file,"Total number of ray-polygon intersection tests: %lld\n",world->ray_polygon_tests);
       fprintf(world->log_file,"Total number of ray-polygon intersections: %lld\n",world->ray_polygon_colls);
       fprintf(world->log_file,"Total number of molecule-molecule collisions: %lld\n",world->mol_mol_colls);
       fprintf(world->log_file,"Total number of molecule-molecule-molecule collisions: %lld\n",world->mol_mol_mol_colls);
    }
 
    t_final = time(NULL);
    getrusage(RUSAGE_SELF,&run_time);
    fprintf( world->log_file,"Total CPU time = %f (user) and %f (system)\n",
             run_time.ru_utime.tv_sec + (run_time.ru_utime.tv_usec/MAX_TARGET_TIMESTEP),
             run_time.ru_stime.tv_sec + (run_time.ru_stime.tv_usec/MAX_TARGET_TIMESTEP) );
    fprintf( world->log_file,"Total wall clock time = %d seconds\n",
             (int)(t_final - t_initial) );
  }

}

/***********************************************************************
 install_usr_signal_handlers:

    Set signal handlers for checkpointing on SIGUSR signals.

    In:  None
    Out: 0 on success, 1 on failure.
 ***********************************************************************/
int install_usr_signal_handlers(void)
{
  struct sigaction sa, saPrev;
  sa.sa_sigaction = NULL;
  sa.sa_handler = &chkpt_signal_handler;
  sa.sa_flags = SA_RESTART;
  sigfillset(&sa.sa_mask);

  if (sigaction(SIGUSR1, &sa, &saPrev) != 0)
  {
    fprintf(stderr, "Failed to install USR1 signal handler\n");
    return 1;
  }
  if (sigaction(SIGUSR2, &sa, &saPrev) != 0)
  {
    fprintf(stderr, "Failed to install USR2 signal handler\n");
    return 1;
  }

  return 0;
}

int main(int argc, char **argv) {

  FILE *err_file;
  FILE *log_file;
  char hostname[64];
  u_int procnum;
  long long exec_iterations = 0; /* number of simulation iterations for this run */
  time_t begin_time_of_day;  /* start time of the simulation */

  /* get the process start time */
  time(&begin_time_of_day);
  if (install_usr_signal_handlers())
    exit(EXIT_FAILURE);

#if defined(__linux__)
  feenableexcept(FE_DIVBYZERO);
#endif

  log_file=stdout;
  err_file=stderr;

  if ((world=(struct volume *)malloc(sizeof(struct volume)))==NULL) {
    fprintf(err_file,"File '%s', Line %ld: Out of memory while creating world volume data structure.\n", __FILE__, (long)__LINE__);
    exit(EXIT_FAILURE);
  }
  memset(world, 0, sizeof(struct volume));
  world->log_file=log_file;
  world->err_file=err_file;

  world->procnum=0;
  procnum=world->procnum;
  gethostname(hostname,64);

  world->iterations=INT_MIN; /* indicates iterations not set */
  world->chkpt_infile = NULL;
  world->chkpt_init = 1;
  world->log_freq = -1; /* Indicates that this value has not been set by user */
  world->begin_timestamp = begin_time_of_day;
  world->initialization_state = "initializing";

  /*
   * Parse the command line arguments and print out errors if necessary.
   */
  if (argparse_init(argc,argv,world)) {
    if (world->log_file!=NULL) {
      log_file=world->log_file;
    }

    if (procnum == 0)
    {
      print_version(log_file);
      print_usage(log_file, argv[0]);
    }

    exit(1);
  }

  log_file=world->log_file;
  
  if (init_notifications())
  {
    fprintf(world->err_file,"File '%s', Line %ld: Could not initialize user-notification data structures.\n", __FILE__, (long)__LINE__);
    exit(1);
  }
  
  if (world->notify->progress_report!=NOTIFY_NONE)
  {
    print_version(log_file);
  }

  if (init_sim()) {
    exit(EXIT_FAILURE);
  }
  world->initialization_state = NULL;

  if(world->chkpt_flag)
  {
  	fprintf(log_file,"MCell: checkpoint sequence number %d begins at elapsed time %1.15g seconds\n", world->chkpt_seq_num, world->chkpt_elapsed_real_time_start);
        if(world->iterations < world->start_time){
  	   fprintf(world->err_file,"Error: start time after checkpoint %lld is greater than total number of iterations specified %lld.\n", world->start_time, world->iterations);
           exit(EXIT_FAILURE);
          
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
        if(exec_iterations < 0) {
  	   fprintf(world->err_file,"Error: number of iterations to execute is zero or negative. Please verify ITERATIONS and/or CHECKPOINT_ITERATIONS commands.\n");
           exit(EXIT_FAILURE);
        }
  	fprintf(log_file,"MCell: executing %lld iterations starting at iteration number %lld.\n",
          exec_iterations,world->start_time);
  }

  if((world->chkpt_flag) && (exec_iterations <= 0))
  {
    mem_dump_stats(stdout);
    exit(0);
  }

  printf("Running...\n");
  run_sim();

  printf("Done running.\n");
  mem_dump_stats(stdout);

  exit(0);
}
