#include <stdlib.h>

#include "mcell_misc.h"
#include "mcell_objects.h"
#include "mcell_react_out.h"
#include "mcell_reactions.h"
#include "mcell_release.h"
#include "mcell_species.h"
#include "mcell_viz.h"

#define CHECKED_CALL_EXIT(function, error_message)                             \
  {                                                                            \
    if (function) {                                                            \
      mcell_print(error_message);                                              \
      exit(1);                                                                 \
    }                                                                          \
  }

/***************************************************************************
 * Test code for the libmcell API
 ***************************************************************************/
void test_api(MCELL_STATE *state) {
  /* set timestep and number of iterations */
  CHECKED_CALL_EXIT(mcell_set_time_step(state, 1e-6), "Failed to set timestep");
  CHECKED_CALL_EXIT(mcell_set_iterations(state, 100), "Failed to set iterations");

  /* create range for partitions */
  struct num_expr_list_head list = {NULL, NULL, 0, 1};
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
  struct mcell_species_spec molA = {"A", 1e-6, 0, 0.0, 0, 0.0, 0.0};
  mcell_symbol *molA_ptr;
  CHECKED_CALL_EXIT(mcell_create_species(state, &molA, &molA_ptr),
    "Failed to create species A");

  struct mcell_species_spec molB = {"B", 1e-5, 0, 0.0, 0, 0.0, 0.0};
  mcell_symbol *molB_ptr;
  CHECKED_CALL_EXIT(mcell_create_species(state, &molB, &molB_ptr),
    "Failed to create species B");

  struct mcell_species_spec molC = {"C", 2e-5, 0, 0.0, 0, 0.0, 0.0};
  mcell_symbol *molC_ptr;
  CHECKED_CALL_EXIT(mcell_create_species(state, &molC, &molC_ptr),
    "Failed to create species C");

  /* create reactions */
  struct mcell_species *reactants = mcell_add_to_species_list(
    molA_ptr, false, 0, 0, NULL);
  reactants = mcell_add_to_species_list(molB_ptr, false, 0, 0, reactants);

  struct mcell_species *products = mcell_add_to_species_list(
    molC_ptr, false, 0, 0, NULL);

  struct mcell_species *surfs = mcell_add_to_species_list(
    NULL, false, 0, 0, NULL);

  struct reaction_arrow arrow = {REGULAR_ARROW, {NULL, NULL, 0, 0, 0}};

  struct reaction_rates rates = mcell_create_reaction_rates(RATE_CONSTANT, 1e6,
    RATE_UNSET, 0.0);

  if (mcell_add_reaction(state, reactants, &arrow, surfs, products, NULL, &rates, NULL) == MCELL_FAIL) {
    mcell_print("error ");
    exit(1);
  }

  mcell_delete_species_list(reactants);
  mcell_delete_species_list(products);
  mcell_delete_species_list(surfs);

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

  struct element_connection_list *elems = mcell_add_to_connection_list(1, 2, 3, NULL);
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

  struct poly_object polygon = {"aBox", verts, 8, elems, 12};
  struct object *new_mesh = NULL;
  CHECKED_CALL_EXIT(mcell_create_poly_object(state, world_object, &polygon, &new_mesh),
    "could not create polygon_object")

  /***************************************************************************
   * begin code for creating release sites
   ***************************************************************************/
  struct vector3 position = {0.0, 0.0, 0.0};
  struct vector3 diameter = {0.00999, 0.00999, 0.00999};
  struct object *A_releaser = NULL;
  struct mcell_species *A = mcell_add_to_species_list(molA_ptr, false, 0, 0, NULL);
  CHECKED_CALL_EXIT(mcell_create_geometrical_release_site(state, world_object,
    "A_releaser", SHAPE_SPHERICAL, &position, &diameter, A, 10000, 1, NULL,
    &A_releaser), "could not create A_releaser");
  mcell_delete_species_list(A);

  struct mcell_species *B = mcell_add_to_species_list(molB_ptr, false, 0, 0, NULL);
  CHECKED_CALL_EXIT(mcell_create_geometrical_release_site(state, world_object,
    "B_releaser", SHAPE_SPHERICAL, &position, &diameter, B, 5000, 1, NULL,
    &A_releaser), "could not create A_releaser");
  mcell_delete_species_list(B);

  /***************************************************************************
   * begin code for creating count statements
   ***************************************************************************/
  //struct sym_table *where = NULL;   // we count in the world
  struct sym_table *where = new_mesh->sym;
  //byte report_flags = REPORT_WORLD;
  //report_flags |= REPORT_CONTENTS;
  byte report_flags = REPORT_CONTENTS;

  struct output_column_list count_list;
  CHECKED_CALL_EXIT(mcell_create_count(state, molA_ptr, ORIENT_NOT_SET, where,
  report_flags, NULL, &count_list), "Failed to create COUNT expression");

  struct output_set *os = mcell_create_new_output_set(state, NULL, 0,
    count_list.column_head, FILE_SUBSTITUTE, "react_data/foobar.dat");

  struct output_times_inlist outTimes;
  outTimes.type = OUTPUT_BY_STEP;
  outTimes.step = 1e-5;

  struct output_set_list output;
  output.set_head = os;
  output.set_tail = os;

  CHECKED_CALL_EXIT(mcell_add_reaction_output_block(state, &output, 10000,
    &outTimes), "Error setting up the reaction output block");

  struct mcell_species *mol_viz_list = mcell_add_to_species_list(molA_ptr, false, 0, 0, NULL);
  mol_viz_list = mcell_add_to_species_list(molB_ptr, false, 0, 0, mol_viz_list);
  mol_viz_list = mcell_add_to_species_list(molC_ptr, false, 0, 0, mol_viz_list);
  CHECKED_CALL_EXIT(
    mcell_create_viz_output(state, "./viz_data/test", mol_viz_list, 0, 50, 2),
    "Error setting up the viz output block");
  mcell_delete_species_list(mol_viz_list);
}
