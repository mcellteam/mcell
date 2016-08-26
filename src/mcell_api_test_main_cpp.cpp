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

extern "C" {

#include "config.h"

#include <stdlib.h>

#include "mcell_init.h"
#include "mcell_misc.h"
#include "mcell_run.h"
//#include "api_test.h"


#include "mcell_objects.h"
#include "mcell_react_out.h"
#include "mcell_reactions.h"
#include "mcell_release.h"
#include "mcell_species.h"
#include "mcell_viz.h"
#include "mcell_surfclass.h"

#include <stdio.h>

}

#include <string.h>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <unordered_map>

#define CHECKED_CALL_EXIT(function, error_message)                             \
  {                                                                            \
    if (function) {                                                            \
      mcell_print(error_message);                                              \
      exit(1);                                                                 \
    }                                                                          \
  }



using namespace std;


// This is the earlier API header
namespace MCellAPI {
  class MCellSpeciesExpression {
    /**
    * This base class provides the operator overloading for C++ species expressions.
    */
    protected:
      string name;
    public:
      MCellSpeciesExpression() {
        this->name = "";
      }
      string getName();
      MCellSpeciesExpression ( string name );
      MCellSpeciesExpression operator+ (MCellSpeciesExpression);
      MCellSpeciesExpression operator> (MCellSpeciesExpression);
      MCellSpeciesExpression operator< (MCellSpeciesExpression);
      MCellSpeciesExpression operator== (MCellSpeciesExpression);
      // MCellSpeciesExpression operator<> (MCellSpeciesExpression);
  };
  class MCellSpecies : public MCellSpeciesExpression {
    /**
    * This class extends to the MCellSpeciesExpression class to include a diffusion constant.
    */
    protected:
      double diffusion_constant;
    public:
      MCellSpecies() {
        /**
        * Species constructor without a name.
        */
        this->name = "";
        this->diffusion_constant = 0;
      }
      MCellSpecies ( string species_name ) {
        /**
        * Species constructor with a name.
        */
        this->name = species_name;
        this->diffusion_constant = 0;
      }
      string to_string ();
      void set_diffusion_constant ( double d );
  };
  class MCellReaction {
    protected:
      string name;
    public:
      MCellReaction ( string name );
      void set_rate ( double d );
  };
  class MCellSim {
    protected:
      string name;
    public:
      MCellSim ( string name );
      MCellSim();
      MCellSpecies new_species ( string name );
      MCellReaction new_reaction ( MCellSpeciesExpression exp );
      void run(int n);
  };
}



// This is the earlier API implementation
namespace MCellAPI {

  MCellSpeciesExpression::MCellSpeciesExpression ( string species_name ) {
    /**
    * SpeciesExpression constructor with a name.
    */
    this->name = species_name;
    cout << "  SpeciesExpression Constructor for " << species_name << endl;
  }

  string MCellSpeciesExpression::getName() {
    /**
    * Access function to return the name.
    */
    return (this->name);
  }

  MCellSpeciesExpression MCellSpeciesExpression::operator+ (MCellSpeciesExpression rhs) {
    /**
    * Overload the "+" operator to combine species. Returns a species joined with "+".
    */
    MCellSpeciesExpression *result = new MCellSpeciesExpression ( this->name + " + " + rhs.name );
    cout << "  Species Expression: " << result->name << endl;
    return (*result);
  }

  MCellSpeciesExpression MCellSpeciesExpression::operator> (MCellSpeciesExpression mce) {
    /**
    * Overload the ">" operator to be the forward reaction symbol.
    */
    this->name = this->name + " > " + mce.name;
    cout << "  Species Expression: " << this->name << endl;
    return (*this);
  }

  MCellSpeciesExpression MCellSpeciesExpression::operator< (MCellSpeciesExpression mce) {
    /**
    * Overload the "<" operator to be the reverse reaction symbol.
    */
    this->name = this->name + " < " + mce.name;
    cout << "  Species Expression: " << this->name << endl;
    return (*this);
  }

  MCellSpeciesExpression MCellSpeciesExpression::operator== (MCellSpeciesExpression mce) {
    /**
    * Overload the "==" operator to be the bidirectional reaction symbol.
    */
    this->name = this->name + " <==> " + mce.name;
    cout << "  Species Expression: " << this->name << endl;
    return (*this);
  }

  string MCellSpecies::to_string() {
    /**
    * Produce a string representation of the species.
    */
    // string s = this->name + " has dc=" + std::to_string(this->diffusion_constant);
    // For more precision and control, use a stringstream to output the value
    std::stringstream ss;
    ss << std::setprecision(20) << this->diffusion_constant;
    string s = this->name + " has dc=" + ss.str();
    return ( s );
  }

  void MCellSpecies::set_diffusion_constant ( double d ) {
    this->diffusion_constant = d;
    cout << "Simulation is setting Diffusion Constant for " << name << " to " << d << endl;
  }

  MCellReaction::MCellReaction ( string reaction_name ) {
    this->name = reaction_name;
    cout << "  Creating a new reaction: " << reaction_name << endl;
  }

  void MCellReaction::set_rate ( double d ) {
    cout << "  Setting reaction rate for " << name << " to " << d << endl;
  }

  MCellSim::MCellSim( string simulation_name ) {
    this->name = simulation_name;
    cout << "Creating a new Simulation: " << simulation_name << endl;
  }

  MCellSim::MCellSim() {
    name = "Sim";
    cout << "Creating a new Simulation: " << name << endl;
  }

  void MCellSim::run ( int n ) {
    cout << "Running simulation \"" << this->name << "\" " << n << " steps." << endl;
  }

  MCellSpecies MCellSim::new_species ( string species_name ) {
    cout << "Simulation is creating a new species: " << species_name << endl;
    return ( MCellSpecies(species_name) );
  }

  MCellReaction MCellSim::new_reaction ( MCellSpeciesExpression exp ) {
    cout << "Simulation is creating a new Reaction Expression: " << exp.getName() << endl;
    return ( MCellReaction(exp.getName()) );
  }

}




using namespace MCellAPI;

int main(int argc, char **argv) {






  // This is an earlier C++ prototype implementation with an example dictionary (map) added

  cout << endl << endl;
  cout << "*********************************" << endl;
  cout << "***   Simple MCell API Demo   ***" << endl;
  cout << "*********************************" << endl << endl;

  MCellSim mysim = MCellSim ( "Demo" );

  MCellSpecies a = mysim.new_species ( "A" );
  MCellSpecies b = mysim.new_species ( "B" );
  MCellSpecies c = mysim.new_species ( "C" );
  MCellSpecies d = mysim.new_species ( "D" );

  a.set_diffusion_constant ( 1.2e-5 );

  cout << "Simulation is creating Reactions: " << endl;
  MCellReaction forward = mysim.new_reaction ( a + b > c + d );
  MCellReaction bidirect = mysim.new_reaction ( a + b == c + d );

  forward.set_rate ( 1e8 );
  bidirect.set_rate ( 1e7 );

  mysim.run ( 10 );

  cout << endl << "**********  Done!!  **********" << endl << endl << endl;


  // Test out an example dictionary (map or unordered_map)

  unordered_map<string, MCellSpecies*> molecules;

  cout << endl << "Assigning molecules to dictionary" << endl;

  molecules["a"] = &a;
  molecules["b"] = &b;
  molecules["c"] = &c;
  molecules["d"] = &d;

  cout << "molecules[b]=" << molecules.at("b")->to_string() << endl << endl;

  cout << "Map size: " << molecules.size() << endl;

  for( unordered_map<string,MCellSpecies*>::iterator ii=molecules.begin(); ii!=molecules.end(); ++ii)
  {
    cout << "  " << (*ii).first << ": " << (*ii).second->to_string() << endl;
  }

  // End of earlier C++ prototype







  u_int procnum = 0;

  // initialize the mcell simulation
  MCELL_STATE *world = mcell_create();
  CHECKED_CALL_EXIT(!world, "Failed to initialize MCell simulation.");

  /*
  // Parse the command line arguments and print out errors if necessary.
  if (mcell_argparse(argc, argv, world)) {
    if (procnum == 0) {
      mcell_print_version();
      mcell_print_usage(argv[0]);
    }
    exit(1);
  }
  */

  CHECKED_CALL_EXIT( mcell_init_state(world), "An error occured during set up of the initial simulation state");

  if (world->notify->progress_report != NOTIFY_NONE) {
    mcell_print_version();
  }

  cout << "\n\n================= MCell C++ API TEST EXAMPLE =================\n\n";


  /* set timestep and number of iterations */
  CHECKED_CALL_EXIT(mcell_set_time_step(world, 1e-6), "Failed to set timestep");
  CHECKED_CALL_EXIT(mcell_set_iterations(world, 1000), "Failed to set iterations");

  /* create range for partitions */
  struct num_expr_list_head list = { NULL, NULL, 0, 1 };
  mcell_generate_range(&list, -0.5, 0.5, 0.05);
  list.shared = 1;

  /* set partitions */
  CHECKED_CALL_EXIT(mcell_set_partition(world, X_PARTS, &list), "Failed to set X partition");
  CHECKED_CALL_EXIT(mcell_set_partition(world, Y_PARTS, &list), "Failed to set Y partition");
  CHECKED_CALL_EXIT(mcell_set_partition(world, Z_PARTS, &list), "Failed to set Z partition");

  /* create species */
  struct mcell_species_spec molA = { "A", 1e-6, 1, 0.0, 0, 0.0, 0.0 };
  mcell_symbol *molA_ptr;
  CHECKED_CALL_EXIT(mcell_create_species(world, &molA, &molA_ptr), "Failed to create species A");

  struct mcell_species_spec molB = { "B", 1e-5, 0, 0.0, 0, 0.0, 0.0 };
  mcell_symbol *molB_ptr;
  CHECKED_CALL_EXIT(mcell_create_species(world, &molB, &molB_ptr), "Failed to create species B");

  struct mcell_species_spec molC = { "C", 2e-5, 0, 0.0, 0, 0.0, 0.0 };
  mcell_symbol *molC_ptr;
  CHECKED_CALL_EXIT(mcell_create_species(world, &molC, &molC_ptr), "Failed to create species C");
  
  struct mcell_species_spec molD = { "D", 1e-6, 1, 0.0, 0, 0.0, 0.0 };
  mcell_symbol *molD_ptr;
  CHECKED_CALL_EXIT(mcell_create_species(world, &molD, &molD_ptr), "Failed to create species D");

  /* create reactions */
  struct mcell_species *reactants = mcell_add_to_species_list(molA_ptr, true,  1, NULL);
  reactants                       = mcell_add_to_species_list(molB_ptr, true, -1, reactants);

  struct mcell_species *products  = mcell_add_to_species_list(molC_ptr, true, -1, NULL);

  struct mcell_species *surfs     = mcell_add_to_species_list(NULL, false, 0, NULL);

  struct reaction_arrow arrow     = { REGULAR_ARROW, { NULL, NULL, 0, 0 } };

  struct reaction_rates rates = mcell_create_reaction_rates(RATE_CONSTANT, 1e7, RATE_UNSET, 0.0);

  if (mcell_add_reaction(world->notify, &world->r_step_release,
                         world->rxn_sym_table, world->radial_subdivisions,
                         world->vacancy_search_dist2, reactants, &arrow, surfs,
                         products, NULL, &rates, NULL, NULL) == MCELL_FAIL) {
    mcell_print("error ");
    exit(1);
  }

  mcell_delete_species_list(reactants);
  mcell_delete_species_list(products);
  mcell_delete_species_list(surfs);

  // create surface class
  mcell_symbol *sc_ptr;
  CHECKED_CALL_EXIT(mcell_create_surf_class(world, "SC_test", &sc_ptr), "Failed to create surface class SC_test");

  // create releases using a surface class (i.e. not a release object)
  
  // mdl equivalent: MOLECULE_DENSITY {A' = 1000}
  struct mcell_species *A = mcell_add_to_species_list(molA_ptr, true, 1, NULL);
  struct sm_dat *smd = mcell_add_mol_release_to_surf_class(world, sc_ptr, A, 1000, 0, NULL);

  // mdl equivalent: MOLECULE_NUMBER {D, = 1000}
  struct mcell_species *D = mcell_add_to_species_list(molD_ptr, true, -1, NULL);
  mcell_add_mol_release_to_surf_class(world, sc_ptr, D, 1000, 1, smd);


  // mdl equivalent: ABSORPTIVE = D
  CHECKED_CALL_EXIT( mcell_add_surf_class_properties(world, SINK, sc_ptr, molD_ptr, 0), "Failed to add surface class property");

  // mdl equivalent: REFLECTIVE = D
  /*CHECKED_CALL_EXIT(*/
  /*  mcell_add_surf_class_properties(world, RFLCT, sc_ptr, molD_ptr, 0),*/
  /*  "Failed to add surface class property");*/
  
  mcell_delete_species_list(A);
  mcell_delete_species_list(D);


  /*****************************************************************************
   * create world meta object
   *****************************************************************************/
  struct object *world_object = NULL;
  CHECKED_CALL_EXIT(mcell_create_instance_object(world, "world", &world_object), "could not create meta object");


  /****************************************************************************
   * begin code for creating a polygon mesh
   ****************************************************************************/

  struct vertex_list *verts = mcell_add_to_vertex_list(0.5, 0.5, -0.5, NULL);
  verts = mcell_add_to_vertex_list(0.5, -0.5, -0.5, verts);
  verts = mcell_add_to_vertex_list(-0.5, -0.5, -0.5, verts);
  verts = mcell_add_to_vertex_list(-0.5, 0.5, -0.5, verts);
  verts = mcell_add_to_vertex_list(0.5, 0.5, 0.5, verts);
  verts = mcell_add_to_vertex_list(0.5, -0.5, 0.5, verts);
  verts = mcell_add_to_vertex_list(-0.5, -0.5, 0.5, verts);
  verts = mcell_add_to_vertex_list(-0.5, 0.5, 0.5, verts);

  struct element_connection_list *elems =
      mcell_add_to_connection_list(1, 2, 3, NULL);
  elems = mcell_add_to_connection_list(7, 6, 5, elems);
  elems = mcell_add_to_connection_list(0, 4, 5, elems);
  elems = mcell_add_to_connection_list(1, 5, 6, elems);
  elems = mcell_add_to_connection_list(6, 7, 3, elems);
  elems = mcell_add_to_connection_list(0, 3, 7, elems);
  elems = mcell_add_to_connection_list(0, 1, 3, elems);
  elems = mcell_add_to_connection_list(4, 7, 5, elems);
  elems = mcell_add_to_connection_list(1, 0, 5, elems);
  elems = mcell_add_to_connection_list(2, 1, 6, elems);
  elems = mcell_add_to_connection_list(2, 6, 3, elems);
  elems = mcell_add_to_connection_list(4, 0, 7, elems);

  struct poly_object polygon = { "aBox", verts, 8, elems, 12 };
  struct object *new_mesh = NULL;
  CHECKED_CALL_EXIT( mcell_create_poly_object(world, world_object, &polygon, &new_mesh), "could not create polygon_object")


  /****************************************************************************
   * begin code for creating a region
   ****************************************************************************/
  struct region *test_region = mcell_create_region(world, new_mesh, "reg");
  struct element_list *region_list = mcell_add_to_region_list(NULL, 0);
  region_list = mcell_add_to_region_list(region_list, 1);
  CHECKED_CALL_EXIT(mcell_set_region_elements(test_region, region_list, 1), "could not finish creating region");


  /****************************************************************************
   * Assign surface class to "test_region" 
   ****************************************************************************/
  mcell_assign_surf_class_to_region(sc_ptr, test_region);

  /***************************************************************************
   * begin code for creating release sites
   ***************************************************************************/
  /*struct object *A_releaser = NULL;*/
  /*struct mcell_species *A =*/
  /*    mcell_add_to_species_list(molA_ptr, true, 1, 0, NULL);*/
  /*CHECKED_CALL_EXIT(mcell_create_region_release(world, world_object, new_mesh,*/
  /*                                              "A_releaser", "reg", A, 1000, 1,*/
  /*                                              NULL, &A_releaser),*/
  /*                  "could not create A_releaser");*/
  /*mcell_delete_species_list(A);*/

  struct vector3 position = { 0.0, 0.0, 0.0 };
  struct vector3 diameter = { 0.00999, 0.00999, 0.00999 };

  struct object *B_releaser = NULL;
  struct mcell_species *B =
      mcell_add_to_species_list(molB_ptr, false, 0, NULL);
  CHECKED_CALL_EXIT(mcell_create_geometrical_release_site(
                        world, world_object, "B_releaser", SHAPE_SPHERICAL,
                        &position, &diameter, B, 3000, 1, NULL, &B_releaser),
                    "could not create B_releaser");
  mcell_delete_species_list(B);


  /***************************************************************************
   * begin code for creating count statements
   ***************************************************************************/
  // struct sym_entry *where = NULL;   // we count in the world
  struct sym_entry *where = new_mesh->sym;
  // byte report_flags = REPORT_WORLD;
  // report_flags |= REPORT_CONTENTS;
  byte report_flags = REPORT_CONTENTS;

  struct output_column_list count_list;
  CHECKED_CALL_EXIT(mcell_create_count(world, molA_ptr, ORIENT_NOT_SET, where,
                                       report_flags, NULL, &count_list),
                    "Failed to create COUNT expression");

  struct output_set *os =
      mcell_create_new_output_set(NULL, 0, count_list.column_head,
                                  FILE_SUBSTITUTE, "test_files/mcell/react_data/foobar.dat");

  struct output_times_inlist outTimes;
  outTimes.type = OUTPUT_BY_STEP;
  outTimes.step = 1e-5;

  struct output_set_list output;
  output.set_head = os;
  output.set_tail = os;

  CHECKED_CALL_EXIT(
      mcell_add_reaction_output_block(world, &output, 10000, &outTimes),
      "Error setting up the reaction output block");

  struct mcell_species *mol_viz_list = mcell_add_to_species_list(molA_ptr, false, 0, NULL);
  mol_viz_list = mcell_add_to_species_list(molB_ptr, false, 0, mol_viz_list);
  mol_viz_list = mcell_add_to_species_list(molC_ptr, false, 0, mol_viz_list);
  mol_viz_list = mcell_add_to_species_list(molD_ptr, false, 0, mol_viz_list);
  CHECKED_CALL_EXIT(mcell_create_viz_output(world, "test_files/mcell/viz_data/test",
                                            mol_viz_list, 0, 1000, 2),
                    "Error setting up the viz output block");
  mcell_delete_species_list(mol_viz_list);




  CHECKED_CALL_EXIT(mcell_init_simulation(world), "An error occured during simulation creation.");

  CHECKED_CALL_EXIT(mcell_init_read_checkpoint(world), "An error occured during initialization and reading of checkpoint.");

  CHECKED_CALL_EXIT(mcell_init_output(world), "An error occured during setting up of output.");

  CHECKED_CALL_EXIT(mcell_run_simulation(world), "Error running mcell simulation.");

  if (world->notify->progress_report != NOTIFY_NONE) {
    mcell_print("Done running.");
  }

  mcell_print_stats();

  exit(0);
}
