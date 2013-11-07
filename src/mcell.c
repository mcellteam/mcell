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

#include <stdlib.h>

#include "libmcell.h"

#define ERROR_EXIT(m) {\
    mcell_print(m); exit(1);\
  }


int main(int argc, char **argv)
{
  u_int procnum = 0;

  // initialize the mcell simulation
  MCELL_STATE *world = mcell_create();
  if (!world) 
    ERROR_EXIT("Failed to initialize MCell simulation.");

  /*
   * Parse the command line arguments and print out errors if necessary.
   */
  if (mcell_argparse(argc,argv,world))
  {
    if (procnum == 0)
    {
      mcell_print_version();
      mcell_print_usage(argv[0]);
    }
    exit(1);
  }

  if (mcell_init_state(world)) 
    ERROR_EXIT("An error occured during set up of the initial simulation state");

  if (world->notify->progress_report!=NOTIFY_NONE)
    mcell_print_version();

  if (mcell_parse_mdl(world)) 
    ERROR_EXIT("An error occured during parsing of the mdl file.");

  if (mcell_init_simulation(world))
    ERROR_EXIT("An error occured during simulation creation.");

  if (mcell_read_checkpoint(world))
    ERROR_EXIT("An error occured during reading of checkpoint.");

  if (mcell_init_output(world))
    ERROR_EXIT("An error occured during setting up of output.");

  if (mcell_run_simulation(world)) 
    ERROR_EXIT("Error running mcell simulation.");

  if (world->notify->progress_report!=NOTIFY_NONE)
    mcell_print("Done running.");

  mcell_print_stats();

  exit(0);
}
