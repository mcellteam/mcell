
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>

#include "strfunc.h"
#include "argparse.h"
#include "mdlparse.h"
#include "vol_util.h"
#include "init.h"

struct volume *world;

void run_sim(void)
{
  struct release_event_queue *req;
  while (world->it_time < world->iterations
/*         && world->it_time < world->chkpt_iterations */
        )
  {
    req = schedule_next( world->releaser );
    
    if (req==NULL) world->it_time++;
    else release_molecules(req);
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

  
  run_sim();
  
  exit(0);
}
