#if defined(__linux__)
#define _GNU_SOURCE 1
#endif

#include <time.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#if defined(__linux__)
#include <fenv.h>
#endif

#include "rng.h"
#include "strfunc.h"
#include "argparse.h"
#include "mdlparse.h"
#include "vol_util.h"
#include "react_output.h"
#include "viz_output.h"
#include "diffuse.h"
#include "init.h"

struct volume *world;

void run_sim(void)
{
  struct rusage run_time;
  long t_initial,t_final;

  struct storage_list *local;
  struct release_event_queue *req;
  struct output_block *obp;
  double next_release_time;
  int i;
  int count;
  long long total_coll_1,total_coll_2;
  double total_len;
  int frequency;
  
/*
  for (i=0;i<10000;i++)
  {
    struct vector3 v;
    pick_displacement(&v,1.0);
    printf("vector = %7.3f %7.3f %7.3f\n",v.x,v.y,v.z);
  }
  return;
*/  

  t_initial = time(NULL);
  
  if (world->iterations < 1000) frequency = 10;
  else if (world->iterations < 100000) frequency = 100;
  else if (world->iterations < 10000000) frequency = 1000;
  else frequency = 10000;
  
  world->diffusion_number = world->diffusion_cumsteps = 0.0;
  
  while (world->it_time < world->iterations
/*         && world->it_time < world->chkpt_iterations */
        )
  {
    req = schedule_next( world->releaser );
    while (req != NULL)
    {
      if ( release_molecules(req) )
      {
	printf("Out of memory while releasing molecules of type %s\n",req->release_site->mol_type->sym->name);
	return;
      }
      fprintf(world->log_file,"Releasing type = %s! \n",req->release_site->mol_type->sym->name);
      req = schedule_next( world->releaser );
    }
    if (world->releaser->error)
    {
      fprintf(world->err_file,"Out of memory while scheduling molecule release. Trying to save intermediate results.\n");
      i = emergency_output();
      fprintf(world->err_file,"%d error%s while saving intermediate results.\n",i,(i==1)?"":"s");
      exit( EXIT_FAILURE );
    }

    obp = schedule_next( world->count_scheduler );
    while (obp != NULL)
    {
      if (update_reaction_output(obp))
      {
	fprintf(world->err_file,"Error while updating reaction output. Trying to save intermediate results\n");
	i = emergency_output();
	fprintf(world->err_file,"%d error%s while saving intermediate results.\n",i,(i==1)?"":"s");
	exit( EXIT_FAILURE );
      }
      obp = schedule_next( world->count_scheduler );
    }
    if (world->count_scheduler->error)
    {
      fprintf(world->err_file,"Out of memory while scheduling molecule release. Trying to save intermediate results.\n");
      i = emergency_output();
      fprintf(world->err_file,"%d error%s while saving intermediate results.\n",i,(i==1)?"":"s");
      exit( EXIT_FAILURE );
    }

    update_frame_data_list(world->frame_data_head);

    i = schedule_anticipate( world->releaser , &next_release_time);
    if (!i) next_release_time = world->iterations + 1;
    if (next_release_time < world->it_time+1) next_release_time = world->it_time+1;
    
    for (local = world->storage_head ; local != NULL ; local = local->next)
    {
      if (local->store->current_time <= world->it_time)
      {
        run_timestep( local->store , next_release_time , (double)world->iterations );
      }
    }
    
    world->it_time++;
    
    if ( (world->it_time % frequency) == 0 )
    {
      printf("Iterations: %d of %d ",world->it_time,world->iterations);
#if 0
      printf("count ");
      for (i=0;i<world->n_species;i++)
      {
        printf("#%s=%d ",
               world->species_list[i]->sym->name,
               world->species_list[i]->population
              );
      }
#endif
      printf("\n");
    }
  }
  if(world->diffusion_number > 0){ 
  	fprintf(stderr,"Average diffusion jump was %.2f timesteps\n",world->diffusion_cumsteps/world->diffusion_number);
  }
  fprintf(stderr,"Total number of molecule-molecule collisions: %f\n",world->mol_mol_colls);

#if 0
  fprintf(stderr,"We have %d directions with a mask of %d\n",world->num_directions,world->directions_mask);
  for (i=0;i<world->num_directions;i++)
  {
    printf("xyz = %.4f %.4f %.4f\n",world->d_step[3*i],world->d_step[3*i+1],world->d_step[3*i+2]);
  }
#endif

#if 0
  count = 0;
  total_coll_1 = 0;
  total_coll_2 = 0;
  total_len = 0.0;
  for (i=0;i<world->n_subvols;i++)
  {
    struct molecule *mol;
    for (mol = world->subvol[i].mol_head ; mol != NULL ; mol = mol->next_v)
    {
      if (mol->properties != NULL)
      {
        if (strcmp(mol->properties->sym->name,"Glw")==0)
        {
          count++;
        }
        else if (strcmp(mol->properties->sym->name,"ACh")==0)
        {
          total_len += mol->path_length;
          total_coll_1 += mol->collisions;
        }
        else if (strcmp(mol->properties->sym->name,"Glu")==0)
        {
          total_coll_2 += mol->collisions;
        }
      }
    }
  }
  printf("%d %.2f %lld %lld\n",count,total_len*world->length_unit,total_coll_1,total_coll_2);
#endif
#if 1
  printf("Reaction counters:\n");
  for (i=0;i<world->rx_hashsize;i++)
  {
    struct rxn *rxp;
    int j,k;
    
    for (rxp = world->reaction_hash[i] ; rxp != NULL ; rxp = rxp->next)
    {
      for (j=0;j<rxp->n_reactants;j++)
      {
        if (j==0) fprintf(world->log_file,"Reaction %s[%d]",rxp->players[0]->sym->name,rxp->geometries[0]);
        else fprintf(world->log_file," + %s[%d]",rxp->players[j]->sym->name,rxp->geometries[j]);
      }
      if (rxp->n_pathways==RX_WINDOW) fprintf(world->log_file," WINDOW");
      else if (rxp->n_pathways==RX_GHOST) fprintf(world->log_file," GHOST");
      fprintf(world->log_file,"\n");
      for (j=0;j<rxp->n_pathways;j++)
      {
        fprintf(world->log_file,"  Count %g: ->",rxp->counter[j]);
        for (k=rxp->product_idx[j] ; k<rxp->product_idx[j+1] ; k++)
        {
          if (rxp->players[k]!=NULL) fprintf(world->log_file," %s{%d}",rxp->players[k]->sym->name,rxp->geometries[k]);
        }
        fprintf(world->log_file,"\n");
      }
    }
  }
#endif
#if 0
  for (i=0;i<world->n_subvols;i++)
  {
    struct molecule *mol;
    int j;
/*    printf("Subvolume %d (%d molecules):\n",i,world->subvol[i].mol_count); */
    
    for (mol = world->subvol[i].mol_head ; mol != NULL ; mol = mol->next_v)
    {
      if (mol->properties == NULL)
      {
/*        printf("Defunct molecule.\n"); */
        continue;
      }
      for (j=0;j<world->n_species;j++) if (mol->properties==world->species_list[j]) break;
      printf("location = %7.3f %7.3f %7.3f %d %9.3f %d\n",
             mol->pos.x*world->length_unit,mol->pos.y*world->length_unit,mol->pos.z*world->length_unit,
             j,mol->path_length*world->length_unit,mol->collisions);
    }
  }
#endif
#if 1
  t_final = time(NULL);
  getrusage(RUSAGE_SELF,&run_time);
  fprintf( world->log_file,"Total CPU time = %f (user) and %f (system)\n",
           run_time.ru_utime.tv_sec + (run_time.ru_utime.tv_usec/MAX_TARGET_TIMESTEP),
           run_time.ru_stime.tv_sec + (run_time.ru_stime.tv_usec/MAX_TARGET_TIMESTEP) );
  fprintf( world->log_file,"Total wall clock time = %d seconds\n",
           (int)(t_final - t_initial) );
#endif

  printf("Exiting run loop.\n");

}

int main(int argc, char **argv) {

  FILE *err_file;
  FILE *log_file;
  char hostname[64];
  u_int procnum;
  

#if defined(__linux__)
  feenableexcept(FE_DIVBYZERO);
#endif

  log_file=stdout;
  err_file=stderr;

  if ((world=(struct volume *)malloc(sizeof(struct volume)))==NULL) {
    fprintf(err_file,"MCell: could not store world volume data structure\n");
    exit(EXIT_FAILURE);
  }
  world->log_file=log_file;
  world->err_file=err_file;

  world->procnum=0;
  procnum=world->procnum;
  gethostname(hostname,64);

  /*
   * Parse the command line arguments and print out errors if necessary.
   */
  if (argparse_init(argc,argv,world)) {
    if (world->log_file!=NULL) {
      log_file=world->log_file;
    }
    fprintf(log_file,"\n");
    fprintf(log_file,"MCell");
    fprintf(log_file,"  Running on %s\n\n",hostname);
    if (procnum == 0) {
      if (world->info_opt) {
        fprintf(log_file,"  -info option selected\n");
      }
      fprintf(log_file,"Usage: %s [options] mdl_file_name\n\n",argv[0]);
      fprintf(log_file,"    options:\n");
      fprintf(log_file,"       [-help]                   print this help message\n");
      fprintf(log_file,"       [-info]                   print MCell info message\n");
      fprintf(log_file,"       [-seed n]                 choose random sequence number (default: 1)\n");
      fprintf(log_file,"       [-iterations n]           override iterations in mdl_file_name\n");
      fprintf(log_file,"       [-logfile log_file_name]  send output log to file (default: stderr)\n");
      fprintf(log_file,"       [-logfreq n]              output log frequency (default: 100)\n");
      fprintf(log_file,"       [-checkpoint_infile checkpoint_file_name]  read checkpoint file \n\n");
    }

    exit(1);
  }

  log_file=world->log_file;
  fprintf(log_file,"\n");
  fprintf(log_file,"MCell");
  fprintf(log_file,"  Running on %s\n\n",hostname);
  if (world->info_opt) {
    fprintf(log_file,"  -info option selected\n");
  }

/*
  no_printf("Parsing MDL file: %s\n",world->mdl_infile_name);
  fflush(stderr);
  if (mdlparse_init(world)) {
    fprintf(log_file,"MCell: error parsing file: %s\n",world->curr_file);
    return(1);
  }
  no_printf("Done parsing MDL file: %s\n",world->mdl_infile_name);
  fflush(stderr);
*/

  if (init_sim()) {
    fprintf(log_file,"MCell: error initializing simulation\n");
    exit(EXIT_FAILURE);
  }

  printf("Running...\n");
  run_sim();
  printf("Done running.\n");
  
  exit(0);
}
