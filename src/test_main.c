
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>

#include "rng.h"
#include "strfunc.h"
#include "argparse.h"
#include "mdlparse.h"
#include "vol_util.h"
#include "diffuse.h"
#include "init.h"

struct volume *world;

void run_sim(void)
{
  struct release_event_queue *req;
  double next_release_time;
  int i;
  
  world->fully_random = 1;
/*
  for (i=0;i<10000;i++)
  {
    struct vector3 v;
    pick_displacement(&v,1.0);
    printf("vector = %7.3f %7.3f %7.3f\n",v.x,v.y,v.z);
  }
  return;
*/  
  
  while (world->it_time < world->iterations
/*         && world->it_time < world->chkpt_iterations */
        )
  {
    req = schedule_next( world->releaser );
    while (req != NULL)
    {
      release_molecules(req);
      printf("Releasing! \n");
      req = schedule_next( world->releaser );
    }

    i = schedule_anticipate( world->releaser , &next_release_time);
    if (!i) next_release_time = world->iterations + 1;
    if (next_release_time < world->it_time+1) next_release_time = world->it_time+1;
    
    for (i=0;i<world->n_subvols;i++)
    {
      run_timestep( &(world->subvol[i]) , next_release_time , (double)world->iterations );
    }
    
    world->it_time++;
    printf("Iterations: %d of %d\n",world->it_time,world->iterations);
  }

  for (i=0;i<world->n_subvols;i++)
  {
    struct molecule *mol = world->subvol[i].mol_head;
    int j;
    while (mol != NULL)
    {
      if (mol->properties == world->species_list[0]) j = 0;
      else j = 1;
      printf("location = %7.3f %7.3f %7.3f %d\n",mol->pos.x,mol->pos.y,mol->pos.z,j);
      mol = mol->next_v;
    }
  }

}

int main(int argc, char **argv) {

  FILE *log_file;
  char hostname[64];
  u_int procnum;
  

  log_file=stderr;

  if ((world=(struct volume *)malloc(sizeof(struct volume)))==NULL) {
    fprintf(log_file,"test_main: could not store world volume data structure\n");
    exit(1);
  }
  world->log_file=log_file;

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
    fprintf(log_file,"test_main");
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
  fprintf(log_file,"test_main");
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
    return(1);
  }

  printf("Running...\n");
  run_sim();
  printf("Done running.\n");
  
  exit(0);
}
