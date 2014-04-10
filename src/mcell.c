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


#define CHECKED_CALL_EXIT(function, error_message) {\
    if (function) {\
      mcell_print(error_message);\
      exit(1);\
    }\
  }


int main(int argc, char **argv)
{
  u_int procnum = 0;

  // initialize the mcell simulation
  MCELL_STATE *world = mcell_create();
  CHECKED_CALL_EXIT(!world, "Failed to initialize MCell simulation.");

  // Parse the command line arguments and print out errors if necessary.
  if (mcell_argparse(argc,argv,world))
  {
    if (procnum == 0)
    {
      mcell_print_version();
      mcell_print_usage(argv[0]);
    }
    exit(1);
  }

  CHECKED_CALL_EXIT(mcell_init_state(world), 
      "An error occured during set up of the initial simulation state");

  if (world->notify->progress_report!=NOTIFY_NONE) {
    mcell_print_version();
  }


#if 0
  /***************************************************************************
   * the below is test code for the libmcell API
   ***************************************************************************/

  /* set timestep and number of iterations */
  CHECKED_CALL_EXIT(mcell_set_time_step(world, 1e-6), "Failed to set timestep");
  CHECKED_CALL_EXIT(mcell_set_iterations(world, 100), "Failed to set iterations");

  /* create species */
  struct mcell_species_spec molA = {"A", 1e-6, 0.0, 0, 0.0, 0, 0.0, 0.0};
  mcell_symbol *molA_ptr;
  CHECKED_CALL_EXIT(mcell_create_species(world, &molA, &molA_ptr), 
    "Failed to create species A");

  struct mcell_species_spec molB = {"B", 1e-5, 0.0, 0, 0.0, 0, 0.0, 0.0};
  mcell_symbol *molB_ptr;
  CHECKED_CALL_EXIT(mcell_create_species(world, &molB, &molB_ptr),
    "Failed to create species B");

  struct mcell_species_spec molC = {"C", 2e-5, 0.0, 0, 0.0, 0, 0.0, 0.0};
  mcell_symbol *molC_ptr;
  CHECKED_CALL_EXIT(mcell_create_species(world, &molC, &molC_ptr),
    "Failed to create species C");

  /* create reactions */
  struct mcell_species *reactants = mcell_add_to_species_list(molA_ptr, false, 0, 0, NULL);
  reactants = mcell_add_to_species_list(molB_ptr, false, 0, 0, reactants);

  struct mcell_species *products = mcell_add_to_species_list(molC_ptr, false, 0, 0, NULL);

  struct mcell_species *surfs = mcell_add_to_species_list(NULL, false, 0, 0, NULL);

  struct reaction_arrow arrow = {REGULAR_ARROW, {NULL, NULL, 0, 0, 0}};

  struct reaction_rates rates = mcell_create_reaction_rates(RATE_CONSTANT, 1e8,
    RATE_UNSET, 0.0);

  if (mcell_add_reaction(world, reactants, &arrow, surfs, products, NULL, &rates, NULL) == MCELL_FAIL) {
    mcell_print("error ");
    exit(1);
  }

  mcell_delete_species_list(reactants);
  mcell_delete_species_list(products);
  mcell_delete_species_list(surfs);


  /***************************************************************************
  * end of libmcell API code
  ***************************************************************************/
#endif

  CHECKED_CALL_EXIT(mcell_parse_mdl(world),
      "An error occured during parsing of the mdl file.");

  CHECKED_CALL_EXIT(mcell_init_simulation(world),
      "An error occured during simulation creation.");

  CHECKED_CALL_EXIT(mcell_read_checkpoint(world),
      "An error occured during reading of checkpoint.");

  CHECKED_CALL_EXIT(mcell_init_output(world),
      "An error occured during setting up of output.");

  CHECKED_CALL_EXIT(mcell_run_simulation(world), "Error running mcell simulation.");

  if (world->notify->progress_report!=NOTIFY_NONE) {
    mcell_print("Done running.");
  }

  mcell_print_stats();

  exit(0);
}
