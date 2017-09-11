#include <stdio.h>
#include <stdlib.h>

#include "mcell_structs.h"
#include "mcell_misc.h"
#include "mcell_objects.h"
#include "mcell_react_out.h"
#include "mcell_reactions.h"
#include "mcell_release.h"
#include "mcell_species.h"
#include "mcell_viz.h"
#include "mcell_surfclass.h"
#include "mcell_run.h"

#define CHECKED_CALL_EXIT(function, error_message)                             \
  {                                                                            \
    if (function) {                                                            \
      mcell_print(error_message);                                              \
      exit(1);                                                                 \
    }                                                                          \
  }

int main(void) {
  struct volume *state = mcell_create();
  CHECKED_CALL_EXIT(
      mcell_init_state(state),
      "An error occured during set up of the initial simulation state");
  /* set timestep and number of iterations */
  CHECKED_CALL_EXIT(mcell_set_time_step(state, 1e-6), "Failed to set timestep");
  CHECKED_CALL_EXIT(mcell_set_iterations(state, 1000),
                    "Failed to set iterations");

  /* create range for partitions */
  struct num_expr_list_head list = { NULL, NULL, 0, 1 };
  mcell_generate_range(&list, -0.5, 0.5, 0.05);
  list.shared = 1;
  /* set partitions */
  CHECKED_CALL_EXIT(mcell_set_partition(state, X_PARTS, &list),
                    "Failed to set X partition");
  CHECKED_CALL_EXIT(mcell_set_partition(state, Y_PARTS, &list),
                    "Failed to set Y partition");
  CHECKED_CALL_EXIT(mcell_set_partition(state, Z_PARTS, &list),
                    "Failed to set Z partition");

  /* create species */
  struct mcell_species_spec molA = { "A", 1e-6, 1, 0.0, 0, 0.0, 0.0 };
  mcell_symbol *molA_ptr;
  CHECKED_CALL_EXIT(mcell_create_species(state, &molA, &molA_ptr),
                    "Failed to create species A");

  struct mcell_species_spec molB = { "B", 1e-5, 0, 0.0, 0, 0.0, 0.0 };
  mcell_symbol *molB_ptr;
  CHECKED_CALL_EXIT(mcell_create_species(state, &molB, &molB_ptr),
                    "Failed to create species B");

  struct mcell_species_spec molC = { "C", 2e-5, 0, 0.0, 0, 0.0, 0.0 };
  mcell_symbol *molC_ptr;
  CHECKED_CALL_EXIT(mcell_create_species(state, &molC, &molC_ptr),
                    "Failed to create species C");
  
  struct mcell_species_spec molD = { "D", 1e-6, 1, 0.0, 0, 0.0, 0.0 };
  mcell_symbol *molD_ptr;
  CHECKED_CALL_EXIT(mcell_create_species(state, &molD, &molD_ptr),
                    "Failed to create species D");

  /* create reactions */
  struct mcell_species *reactants =
      mcell_add_to_species_list(molA_ptr, true, 1, NULL);
  reactants = mcell_add_to_species_list(molB_ptr, true, -1, reactants);

  struct mcell_species *products =
      mcell_add_to_species_list(molC_ptr, true, -1, NULL);

  struct mcell_species *surfs =
      mcell_add_to_species_list(NULL, false, 0, NULL);

  struct reaction_arrow arrow = { REGULAR_ARROW, { NULL, NULL, 0, 0 } };

  struct reaction_rates rates =
      mcell_create_reaction_rates(RATE_CONSTANT, 1e7, RATE_UNSET, 0.0);

  if (mcell_add_reaction(state->notify, &state->r_step_release,
                         state->rxn_sym_table, state->radial_subdivisions,
                         state->vacancy_search_dist2, reactants, &arrow, surfs,
                         products, NULL, &rates, NULL, NULL) == MCELL_FAIL) {
    mcell_print("error ");
    exit(1);
  }

  mcell_delete_species_list(reactants);
  mcell_delete_species_list(products);
  mcell_delete_species_list(surfs);

  // create surface class
  mcell_symbol *sc_ptr;
  CHECKED_CALL_EXIT(mcell_create_surf_class(state, "SC_test", &sc_ptr),
                    "Failed to create surface class SC_test");

  // create releases using a surface class (i.e. not a release object)
  
  // mdl equivalent: MOLECULE_DENSITY {A' = 1000}
  struct mcell_species *A =
      mcell_add_to_species_list(molA_ptr, true, 1, NULL);
  struct sm_dat *smd = mcell_add_mol_release_to_surf_class(
      state, sc_ptr, A, 1000, 0, NULL);

  // mdl equivalent: MOLECULE_NUMBER {D, = 1000}
  struct mcell_species *D =
      mcell_add_to_species_list(molD_ptr, true, -1, NULL);
  mcell_add_mol_release_to_surf_class(state, sc_ptr, D, 1000, 1, smd);

  // mdl equivalent: ABSORPTIVE = D
  CHECKED_CALL_EXIT(
    mcell_add_surf_class_properties(state, SINK, sc_ptr, molD_ptr, 0),
    "Failed to add surface class property");

  // mdl equivalent: REFLECTIVE = D
  /*CHECKED_CALL_EXIT(*/
  /*  mcell_add_surf_class_properties(state, RFLCT, sc_ptr, molD_ptr, 0),*/
  /*  "Failed to add surface class property");*/
  
  mcell_delete_species_list(A);
  mcell_delete_species_list(D);

  /*****************************************************************************
   * create world meta object
   *****************************************************************************/
  struct object *world_object = NULL;
  CHECKED_CALL_EXIT(mcell_create_instance_object(state, "world", &world_object),
                    "could not create meta object");

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
  CHECKED_CALL_EXIT(
      mcell_create_poly_object(state, world_object, &polygon, &new_mesh),
      "could not create polygon_object")

  /****************************************************************************
   * begin code for creating a region
   ****************************************************************************/
  struct region *test_region = mcell_create_region(state, new_mesh, "reg");
  struct element_list *region_list = mcell_add_to_region_list(NULL, 0);
  region_list = mcell_add_to_region_list(region_list, 1);
  CHECKED_CALL_EXIT(mcell_set_region_elements(test_region, region_list, 1),
                    "could not finish creating region");

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
  /*CHECKED_CALL_EXIT(mcell_create_region_release(state, world_object, new_mesh,*/
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
                        state, world_object, "B_releaser", SHAPE_SPHERICAL,
                        &position, &diameter, B, 5000, 0, 1, NULL, &B_releaser),
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
  CHECKED_CALL_EXIT(mcell_create_count(state, molA_ptr, ORIENT_NOT_SET, where,
                                       report_flags, NULL, &count_list),
                    "Failed to create COUNT expression");

  struct output_set *os =
      mcell_create_new_output_set(NULL, 0, count_list.column_head,
                                  FILE_SUBSTITUTE, "react_data/foobar.dat");

  struct output_times_inlist outTimes;
  outTimes.type = OUTPUT_BY_STEP;
  outTimes.step = 1e-5;

  struct output_set_list output;
  output.set_head = os;
  output.set_tail = os;

  CHECKED_CALL_EXIT(
      mcell_add_reaction_output_block(state, &output, 10000, &outTimes),
      "Error setting up the reaction output block");

  struct mcell_species *mol_viz_list =
      mcell_add_to_species_list(molA_ptr, false, 0, NULL);
  mol_viz_list = mcell_add_to_species_list(molB_ptr, false, 0, mol_viz_list);
  mol_viz_list = mcell_add_to_species_list(molC_ptr, false, 0, mol_viz_list);
  mol_viz_list = mcell_add_to_species_list(molD_ptr, false, 0, mol_viz_list);
  CHECKED_CALL_EXIT(mcell_create_viz_output(state, "./viz_data/test",
                                            mol_viz_list, 0, 1000, 2),
                    "Error setting up the viz output block");
  mcell_delete_species_list(mol_viz_list);
  
  CHECKED_CALL_EXIT(mcell_init_simulation(state),
                    "An error occured during simulation creation.");

  CHECKED_CALL_EXIT(
      mcell_init_read_checkpoint(state),
      "An error occured during initialization and reading of checkpoint.");

  CHECKED_CALL_EXIT(mcell_init_output(state),
                    "An error occured during setting up of output.");

  CHECKED_CALL_EXIT(mcell_run_simulation(state),
                    "Error running mcell simulation.");

  mcell_print_stats();

  return 0;
}
