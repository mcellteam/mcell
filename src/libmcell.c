/***********************************************************************************
 *                                                                                 *
 * Copyright (C) 2006-2014 by *
 * The Salk Institute for Biological Studies and *
 * Pittsburgh Supercomputing Center, Carnegie Mellon University *
 *                                                                                 *
 * This program is free software; you can redistribute it and/or *
 * modify it under the terms of the GNU General Public License *
 * as published by the Free Software Foundation; either version 2 *
 * of the License, or (at your option) any later version. *
 *                                                                                 *
 * This program is distributed in the hope that it will be useful, *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the *
 * GNU General Public License for more details. *
 *                                                                                 *
 * You should have received a copy of the GNU General Public License *
 * along with this program; if not, write to the Free Software *
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301,
 *USA. *
 *                                                                                 *
 ***********************************************************************************/

#if defined(__linux__)
#define _GNU_SOURCE 1
#endif

#ifndef _WIN32
#include <sys/resource.h>
#endif
#include <stdlib.h>
#if defined(__linux__)
#include <fenv.h>
#endif

#include <assert.h>
#include <signal.h>
#include <string.h>
#include <time.h>
#include <math.h>

#include "argparse.h"
#include "chkpt.h"
#include "count_util.h"
#include "diffuse_util.h"
#include "init.h"
#include "libmcell.h"
#include "logging.h"
#include "react_output.h"
#include "react_util.h"
#include "sym_table.h"
#include "version_info.h"
#include "create_species.h"
#include "create_object.h"
#include "create_geometry.h"

#include "mcell_reactions.h"
#include "mcell_release.h"

/* simple wrapper for executing the supplied function call. In case
 * of an error returns with MCELL_FAIL and prints out error_message */
#define CHECKED_CALL(function, error_message)                                  \
  {                                                                            \
    if (function) {                                                            \
      mcell_log(error_message);                                                \
      return MCELL_FAIL;                                                       \
    }                                                                          \
  }

/* declaration of static functions */
static int install_usr_signal_handlers(void);
void swap_double(double *x, double *y);


/************************************************************************
 *
 * function for initializing the main mcell simulator. MCELL_STATE
 * keeps track of the state of the simulation.
 *
 * Returns NULL on error and a pointer to MCELL_STATE otherwise
 *
 ************************************************************************/
MCELL_STATE *
mcell_create() {
  // signal handlers
  if (install_usr_signal_handlers()) {
    return NULL;
  }

  // logging
  mcell_set_log_file(stdout);
  mcell_set_error_file(stderr);

  MCELL_STATE *state = CHECKED_MALLOC_STRUCT_NODIE(struct volume, "world");
  if (state == NULL) {
    return NULL;
  }
  memset(state, 0, sizeof(struct volume));

#if defined(__linux__)
  feenableexcept(FE_DIVBYZERO);
#endif

  state->procnum = 0;
  state->rx_hashsize = 0;
  state->iterations = INT_MIN; /* indicates iterations not set */
  state->chkpt_infile = NULL;
  state->chkpt_outfile = NULL;
  state->chkpt_init = 1;
  state->log_freq =
      ULONG_MAX; /* Indicates that this value has not been set by user */
  state->seed_seq = 1;
  state->with_checks_flag = 1;

  time_t begin_time_of_day;
  time(&begin_time_of_day);
  state->begin_timestamp = begin_time_of_day;
  state->initialization_state = "initializing";

  if (!(state->var_sym_table = init_symtab(1024))) {
    mcell_log("Failed to initialize MDL variable symbol table.");
    return NULL;
  }

  return state;
}

/************************************************************************
 *
 * function for initializing the intial simulation state (variables,
 * notifications, data structures)
 *
 * Returns 1 on error and 0 on success
 *
 ************************************************************************/
MCELL_STATUS
mcell_init_state(MCELL_STATE *state) {
  CHECKED_CALL(
      init_notifications(state),
      "Unknown error while initializing user-notification data structures.");

  CHECKED_CALL(init_variables(state),
               "Unknown error while initializing system variables.");

  CHECKED_CALL(init_data_structures(state),
               "Unknown error while initializing system data structures.");

  return MCELL_SUCCESS;
}

/************************************************************************
 *
 * function for parsing the models underlying mdl file. The function
 * updates the state accordingly.
 *
 * Returns 0 on sucess and 1 on error
 *
 * NOTE: This is currently just a very thin wrapper around parse_input()
 *
 ************************************************************************/
MCELL_STATUS
mcell_parse_mdl(MCELL_STATE *state) { return parse_input(state); }

/************************************************************************
 *
 * function for setting up all the internal data structure to get the
 * simulation into a runnable state.
 *
 * NOTE: Before this function can be called the engine user code
 *       either needs to call
 *       - mcell_parse_mdl() to parse a valid MDL file or
 *       - the individiual API functions for adding model elements
 *         (molecules, geometry, ...)
 *         XXX: These functions don't exist yet!
 *
 * Returns 0 on sucess and 1 on error
 *
 * NOTE: This is currently just a very thin wrapper around parse_input()
 *
 ************************************************************************/
MCELL_STATUS
mcell_init_simulation(MCELL_STATE *state) {
  CHECKED_CALL(init_reactions(state), "Error initializing reactions.");

  CHECKED_CALL(init_species(state), "Error initializing species.");

  if (state->notify->progress_report != NOTIFY_NONE)
    mcell_log("Creating geometry (this may take some time)");

  CHECKED_CALL(init_geom(state), "Error initializing geometry.");
  CHECKED_CALL(init_partitions(state), "Error initializing partitions.");
  CHECKED_CALL(init_vertices_walls(state),
               "Error initializing vertices and walls.");
  CHECKED_CALL(init_regions(state), "Error initializing regions.");

  if (state->place_waypoints_flag) {
    CHECKED_CALL(place_waypoints(state), "Error while placing waypoints.");
  }

  if (state->with_checks_flag) {
    CHECKED_CALL(check_for_overlapped_walls(state->n_subvols, state->subvol),
                 "Error while checking for overlapped walls.");
  }

  CHECKED_CALL(init_surf_mols(state),
               "Error while placing surface molecules on regions.");

  CHECKED_CALL(init_releases(state), "Error while initializing release sites.");

  CHECKED_CALL(init_counter_name_hash(state),
               "Error while initializing counter name hash.");

  return MCELL_SUCCESS;
}

/************************************************************************
 *
 * function for reading and initializing the checkpoint if requested
 *
 * Returns 1 on error and 0 on success
 *
 ************************************************************************/
MCELL_STATUS
mcell_init_read_checkpoint(MCELL_STATE *state) {

  if (state->chkpt_flag == 1) {
    long long exec_iterations;
    CHECKED_CALL(init_checkpoint_state(state, &exec_iterations),
                 "Error while initializing checkpoint.");

    /* XXX This is a hack to be backward compatible with the previous
     * MCell behaviour. Basically, as soon as exec_iterations <= 0
     * MCell will stop and we emulate this by returning 1 even though
     * this is not an error (as implied by returning 1). */
    if (exec_iterations <= 0) {
      mem_dump_stats(mcell_get_log_file());
      return MCELL_FAIL;
    }
  } else {
    state->chkpt_seq_num = 1;
  }

  if (state->chkpt_infile) {
    CHECKED_CALL(load_checkpoint(state),
                 "Error while loading previous checkpoint.");
  }

  // set the iteration time to the start time of the checkpoint
  state->it_time = state->start_time;

  return MCELL_SUCCESS;
}

/************************************************************************
 *
 * function for initializing the viz and reaction data output
 *
 * XXX: This function has to be called last, i.e. after the
 *      simulation has been initialized and checkpoint information
 *      been read.
 *
 * Returns 1 on error and 0 on success
 *
 ************************************************************************/
MCELL_STATUS
mcell_init_output(MCELL_STATE *state) {
  CHECKED_CALL(init_viz_data(state), "Error while initializing viz data.");
  CHECKED_CALL(init_reaction_data(state),
               "Error while initializing reaction data.");
  CHECKED_CALL(init_timers(state), "Error initializing the simulation timers.");

  // signal successful end of simulation
  state->initialization_state = NULL;

  return MCELL_SUCCESS;
}

/**************************************************************************
 * What follows are API functions for adding model elements independent of the
 * parser
 **************************************************************************/

/*************************************************************************
 mcell_set_partition:
    Set the partitioning in a particular dimension.

 In:  state: the simulation state
      dim: the dimension whose partitions we'll set
      head: the partitioning
 Out: 0 on success, 1 on failure
*************************************************************************/
MCELL_STATUS mcell_set_partition(MCELL_STATE *state, int dim,
                                 struct num_expr_list_head *head) {
  /* Allocate array for partitions */
  double *dblp = CHECKED_MALLOC_ARRAY(double, (head->value_count + 2),
                                      "volume partitions");
  if (dblp == NULL)
    return MCELL_FAIL;

  /* Copy partitions in sorted order to the array */
  unsigned int num_values = 0;
  dblp[num_values++] = -GIGANTIC;
  struct num_expr_list *nel;
  for (nel = head->value_head; nel != NULL; nel = nel->next)
    dblp[num_values++] = nel->value * state->r_length_unit;
  dblp[num_values++] = GIGANTIC;
  qsort(dblp, num_values, sizeof(double), &double_cmp);

  /* Copy the partitions into the model */
  switch (dim) {
  case X_PARTS:
    if (state->x_partitions != NULL)
      free(state->x_partitions);
    state->nx_parts = num_values;
    state->x_partitions = dblp;
    break;

  case Y_PARTS:
    if (state->y_partitions != NULL)
      free(state->y_partitions);
    state->ny_parts = num_values;
    state->y_partitions = dblp;
    break;

  case Z_PARTS:
    if (state->z_partitions != NULL)
      free(state->z_partitions);
    state->nz_parts = num_values;
    state->z_partitions = dblp;
    break;

  default:
    UNHANDLED_CASE(dim);
  }

  if (!head->shared)
    mcell_free_numeric_list(head->value_head);

  return MCELL_SUCCESS;
}

/*************************************************************************
 mcell_create_species:
    Create a new species. This uses the same helper functions as the parser,
    but is meant to be used independent of the parser.

 In: state: the simulation state
     name:  molecule name
     D:     diffusion constant
     is_2d: 1 if the species is a 2D molecule, 0 if 3D
     custom_time_step: time_step for the molecule (< 0.0 for a custom space
                       step, >0.0 for custom timestep, 0.0 for default
                       timestep)
     target_only: 1 if the molecule cannot initiate reactions
     max_step_length:
 Out: Returns 0 on sucess and 1 on error
*************************************************************************/
MCELL_STATUS
mcell_create_species(MCELL_STATE *state, struct mcell_species_spec *species,
                     mcell_symbol **species_ptr) {
  struct sym_table *sym =
      CHECKED_MALLOC_STRUCT(struct sym_table, "sym table entry");
  int error_code = new_mol_species(state, species->name, sym);
  if (error_code) {
    return error_code;
  }

  assemble_mol_species(state, sym, species);

  error_code = ensure_rdstep_tables_built(state);
  if (error_code) {
    return error_code;
  }

  if (species_ptr != NULL) {
    *species_ptr = sym;
  }

  return MCELL_SUCCESS;
}

/*************************************************************************
 mcell_set_iterations:
    Set the number of iterations for the simulation.

 In: state: the simulation state
     iterations: number of iterations to run
 Out: 0 on success; 1 on failure.
      number of iterations is set.
*************************************************************************/
MCELL_STATUS
mcell_set_iterations(MCELL_STATE *state, long long iterations) {
  if (iterations < 0) {
    return MCELL_FAIL;
  }
  state->iterations = iterations;
  return MCELL_SUCCESS;
}

/*************************************************************************
 mcell_set_time_step:
    Set the global timestep for the simulation.

 In: state: the simulation state
      step: timestep to set
 Out: 0 on success; any other integer value is a failure.
      global timestep is updated.
*************************************************************************/
MCELL_STATUS
mcell_set_time_step(MCELL_STATE *state, double step) {
  if (step <= 0) {
    return 2;
  }
  // Timestep was already set. Could introduce subtle problems if we let it
  // change after defining the species, since it is used in calculations there.
  if (state->time_unit != 0) {
    return 3;
  }
  state->time_unit = step;
  return MCELL_SUCCESS;
}

/*************************************************************************
 mcell_create_meta_object:
  Create a new meta object.

 In: state:    the simulation state
               object pointer to store created meta object
 Out: 0 on success; any other integer value is a failure.
      A mesh is created.
*************************************************************************/
MCELL_STATUS
mcell_create_instance_object(MCELL_STATE *state, char *name,
                             struct object **new_object) {
  // Create the symbol, if it doesn't exist yet.
  struct object *obj_ptr = make_new_object(state, name);
  if (obj_ptr == NULL) {
    return MCELL_FAIL;
  }
  obj_ptr->last_name = name;
  obj_ptr->object_type = META_OBJ;

  // instantiate object
  obj_ptr->parent = state->root_instance;
  add_child_objects(state->root_instance, obj_ptr, obj_ptr);

  *new_object = obj_ptr;

  return MCELL_SUCCESS;
}

/*************************************************************************
 mcell_create_poly_object:
  Create a new polygon object.

 In: state:    the simulation state
     poly_obj: all the information needed to create the polygon object (name,
               vertices, connections)
 Out: 0 on success; any other integer value is a failure.
      A mesh is created.
*************************************************************************/
MCELL_STATUS
mcell_create_poly_object(MCELL_STATE *state, struct object *parent,
                         struct poly_object *poly_obj,
                         struct object **new_object) {
  // create qualified object name
  char *qualified_name =
      CHECKED_SPRINTF("%s.%s", parent->sym->name, poly_obj->obj_name);

  // Create the symbol, if it doesn't exist yet.
  struct object *obj_ptr = make_new_object(state, qualified_name);
  if (obj_ptr == NULL) {
    return MCELL_FAIL;
  }
  obj_ptr->last_name = qualified_name;

  // Create the actual polygon object
  new_polygon_list(state, obj_ptr, poly_obj->num_vert, poly_obj->vertices,
                   poly_obj->num_conn, poly_obj->connections);

  // Do some clean-up.
  remove_gaps_from_regions(obj_ptr);
  if (check_degenerate_polygon_list(obj_ptr)) {
    return MCELL_FAIL;
  }

  // Set the parent of the object to be the root object. Not reciprocal until
  // add_child_objects is called.
  obj_ptr->parent = parent;
  add_child_objects(parent, obj_ptr, obj_ptr);

  *new_object = obj_ptr;

  return MCELL_SUCCESS;
}

/*************************************************************************
 make_new_object:
    Create a new object, adding it to the global symbol table.

 In:  state: system state
      obj_name: fully qualified object name
 Out: the newly created object
*************************************************************************/
struct object *
make_new_object(MCELL_STATE *state, char *obj_name) {
  if ((retrieve_sym(obj_name, state->obj_sym_table)) != NULL) {
    // mdlerror_fmt(parse_state,"Object '%s' is already defined", obj_name);
    return NULL;
  }

  struct sym_table *symbol;
  if ((symbol = store_sym(obj_name, OBJ, state->obj_sym_table, NULL)) == NULL) {
    return NULL;
  }

  return (struct object *)symbol->value;
}

/**************************************************************************
 *
 * The following functions are likely too low-level to be a part of the API.
 * However they are currently needed by the parser. Eventually, we should
 * try to merge these into other higher-level functions.
 *
 **************************************************************************/

/*************************************************************************
 start_object:
    Create a new object, adding it to the global symbol table. The object must
    not be defined yet. The qualified name of the object will be built by
    adding to the object_name_list, and the object is made the "current_object"
    in the mdl parser state. Because of these side effects, it is vital to call
    finish_object at the end of the scope of the object created here.

 In:  state: the simulation state
      obj_creation: information about object being created
      name: unqualified object name
 Out: the newly created object
 NOTE: This is very similar to mdl_start_object, but there is no parse state.
*************************************************************************/
struct object *
start_object(MCELL_STATE *state,
             struct object_creation *obj_creation, char *name) {
  // Create new fully qualified name.
  char *new_name;
  if ((new_name = push_object_name(obj_creation, name)) == NULL) {
    free(name);
    return NULL;
  }

  // Create the symbol, if it doesn't exist yet.
  struct object *obj_ptr = make_new_object(state, new_name);
  if (obj_ptr == NULL) {
    free(name);
    free(new_name);
    return NULL;
  }

  obj_ptr->last_name = name;
  no_printf("Creating new object: %s\n", new_name);

  // Set parent object, make this object "current".
  obj_ptr->parent = obj_creation->current_object;

  return obj_ptr;
}

/**************************************************************************
 new_polygon_list:
    Create a new polygon list object.

 In: state: the simulation state
     obj_ptr: contains information about the object (name, etc)
     n_vertices: count of vertices
     vertices: list of vertices
     n_connections: count of walls
     connections: list of walls
 Out: polygon object, or NULL if there was an error
 NOTE: This is similar to mdl_new_polygon_list
**************************************************************************/
struct polygon_object *
new_polygon_list(MCELL_STATE *state, struct object *obj_ptr, int n_vertices,
                 struct vertex_list *vertices, int n_connections,
                 struct element_connection_list *connections) {

  struct polygon_object *poly_obj_ptr =
      allocate_polygon_object("polygon list object");
  if (poly_obj_ptr == NULL) {
    goto failure;
  }

  obj_ptr->object_type = POLY_OBJ;
  obj_ptr->contents = poly_obj_ptr;

  poly_obj_ptr->n_walls = n_connections;
  poly_obj_ptr->n_verts = n_vertices;

  // Allocate and initialize removed sides bitmask
  poly_obj_ptr->side_removed = new_bit_array(poly_obj_ptr->n_walls);
  if (poly_obj_ptr->side_removed == NULL) {
    goto failure;
  }
  set_all_bits(poly_obj_ptr->side_removed, 0);

  // Keep temporarily information about vertices in the form of
  // "parsed_vertices"
  poly_obj_ptr->parsed_vertices = vertices;

  // Copy in vertices and normals
  struct vertex_list *vert_list = poly_obj_ptr->parsed_vertices;
  for (int i = 0; i < poly_obj_ptr->n_verts; i++) {
    // Rescale vertices coordinates
    vert_list->vertex->x *= state->r_length_unit;
    vert_list->vertex->y *= state->r_length_unit;
    vert_list->vertex->z *= state->r_length_unit;

    vert_list = vert_list->next;
  }

  // Allocate wall elements
  struct element_data *elem_data_ptr = NULL;
  if ((elem_data_ptr =
           CHECKED_MALLOC_ARRAY(struct element_data, poly_obj_ptr->n_walls,
                                "polygon list object walls")) == NULL) {
    goto failure;
  }
  poly_obj_ptr->element = elem_data_ptr;

  // Copy in wall elements
  for (int i = 0; i < poly_obj_ptr->n_walls; i++) {
    if (connections->n_verts != 3) {
      // mdlerror(parse_state, "All polygons must have three vertices.");
      goto failure;
    }

    struct element_connection_list *elem_conn_list_temp = connections;
    memcpy(elem_data_ptr[i].vertex_index, connections->indices,
           3 * sizeof(int));
    connections = connections->next;
    free(elem_conn_list_temp->indices);
    free(elem_conn_list_temp);
  }

  // Create object default region on polygon list object:
  struct region *reg_ptr = NULL;
  if ((reg_ptr = create_region(state, obj_ptr, "ALL")) == NULL) {
    goto failure;
  }
  if ((reg_ptr->element_list_head =
           new_element_list(0, poly_obj_ptr->n_walls - 1)) == NULL) {
    goto failure;
  }

  obj_ptr->n_walls = poly_obj_ptr->n_walls;
  obj_ptr->n_verts = poly_obj_ptr->n_verts;
  if (normalize_elements(reg_ptr, 0)) {
    // mdlerror_fmt(parse_state,
    //             "Error setting up elements in default 'ALL' region in the "
    //             "polygon object '%s'.", sym->name);
    goto failure;
  }

  return poly_obj_ptr;

failure:
  free_connection_list(connections);
  free_vertex_list(vertices);
  if (poly_obj_ptr) {
    if (poly_obj_ptr->element) {
      free(poly_obj_ptr->element);
    }
    if (poly_obj_ptr->side_removed) {
      free_bit_array(poly_obj_ptr->side_removed);
    }
    free(poly_obj_ptr);
  }
  return NULL;
}

/**************************************************************************
 finish_polygon_list:
    Finalize the polygon list, cleaning up any state updates that were made
    when we started creating the polygon.

 In: obj_ptr: contains information about the object (name, etc)
     obj_creation: information about object being created
 Out: 1 on failure, 0 on success
 NOTE: This function call might be too low-level for what we want from the API,
       but it is needed to create polygon objects for now.
**************************************************************************/
int finish_polygon_list(struct object *obj_ptr,
                        struct object_creation *obj_creation) {
  pop_object_name(obj_creation);
  remove_gaps_from_regions(obj_ptr);
  // no_printf(" n_verts = %d\n", mpvp->current_polygon->n_verts);
  // no_printf(" n_walls = %d\n", mpvp->current_polygon->n_walls);
  if (check_degenerate_polygon_list(obj_ptr)) {
    return 1;
  }
  return 0;
}

/**************************************************************************
 *
 * what follows are helper functions *not* part of the actual API.
 *
 * XXX: Many of these functions should not be called from client
 *      code and need to be removed eventually.
 *
 **************************************************************************/

/***********************************************************************
 * install_usr_signal_handlers:
 *
 *   Set signal handlers for checkpointing on SIGUSR signals.
 *
 *   In:  None
 *   Out: 0 on success, 1 on failure.
 ***********************************************************************/
static int install_usr_signal_handlers(void) {
#ifndef _WIN32 /* fixme: Windows does not support USR signals */
  struct sigaction sa, saPrev;
  sa.sa_sigaction = NULL;
  sa.sa_handler = &chkpt_signal_handler;
  sa.sa_flags = SA_RESTART;
  sigfillset(&sa.sa_mask);

  if (sigaction(SIGUSR1, &sa, &saPrev) != 0) {
    mcell_error("Failed to install USR1 signal handler.");
    return 1;
  }
  if (sigaction(SIGUSR2, &sa, &saPrev) != 0) {
    mcell_error("Failed to install USR2 signal handler.");
    return 1;
  }
#endif

  return 0;
}

/************************************************************************
 *
 * mcell_print_version prints the version string
 *
 ************************************************************************/
void mcell_print_version() { print_version(mcell_get_log_file()); }

/************************************************************************
 *
 * mcell_print_usage prints the usage information
 *
 ************************************************************************/
void mcell_print_usage(const char *executable_name) {
  print_usage(mcell_get_log_file(), executable_name);
}

/************************************************************************
 *
 * mcell_print_stats prints the simulation stats
 *
 ************************************************************************/
void mcell_print_stats() { mem_dump_stats(mcell_get_log_file()); }

/************************************************************************
 *
 * function for printing a string
 *
 * XXX: This is a temporary hack to be able to print in mcell.c
 *      since mcell disables regular printf
 *
 ************************************************************************/
void mcell_print(const char *message) { mcell_log("%s", message); }

/************************************************************************
 *
 * mcell_argparse parses the commandline and sets the
 * corresponding parts of the state (seed #, logging, ...)
 *
 ************************************************************************/
int mcell_argparse(int argc, char **argv, MCELL_STATE *state) {
  return argparse_init(argc, argv, state);
}

/*****************************************************************************
 *
 * mcell_add_to_vertex_list creates a linked list of mesh vertices belonging
 * to a polygon object
 *
 * During the first invocation of this function, NULL should be provided for
 * vertices to initialize a new vertex list. On subsecquent invocations the
 * current vertex_list should be provided as parameter vertices to which the
 * new vertex will be appended.
 *
 *****************************************************************************/
struct vertex_list *mcell_add_to_vertex_list(double x, double y, double z,
                                             struct vertex_list *vertices) {
  struct vertex_list *verts = (struct vertex_list *)CHECKED_MALLOC_STRUCT(
      struct vertex_list, "vertex list");
  if (verts == NULL) {
    return NULL;
  }

  struct vector3 *v =
      (struct vector3 *)CHECKED_MALLOC_STRUCT(struct vector3, "vector");
  if (v == NULL) {
    return NULL;
  }
  v->x = x;
  v->y = y;
  v->z = z;

  verts->vertex = v;
  verts->next = vertices;

  return verts;
}

/*****************************************************************************
 *
 * mcell_add_to_connection_list creates a linked list of element connections
 * describing a polygon object.
 *
 * During the first invocation of this function, NULL should be provided for
 * elements to initialize a new element connection list. On subsecquent
 * invocations the current element_connection_list should be provided as
 * parameter elements to which the new element connection will be appended.
 *
 *****************************************************************************/
struct element_connection_list *
mcell_add_to_connection_list(int v1, int v2, int v3,
                             struct element_connection_list *elements) {
  struct element_connection_list *elems =
      (struct element_connection_list *)CHECKED_MALLOC_STRUCT(
          struct element_connection_list, "element connection list");
  if (elems == NULL) {
    return NULL;
  }

  int *e = (int *)CHECKED_MALLOC_ARRAY(int, 3, "element connections");
  if (e == NULL) {
    return NULL;
  }
  e[0] = v1;
  e[1] = v2;
  e[2] = v3;

  elems->n_verts = 3;
  elems->indices = e;
  elems->next = elements;

  return elems;
}

/*****************************************************************************
 *
 * mcell_add_to_species_list creates a linked list of mcell_species from
 * mcell_symbols.
 *
 * The list of mcell_species is for example used to provide the list
 * of reactants, products and surface classes needed for creating
 * reactions.
 *
 * During the first invocation of this function, NULL should be provided for
 * the species_list to initialize a new mcell_species list with mcell_symbol.
 * On subsecquent invocations the current mcell_species list should
 * be provided as species_list to which the new mcell_symbol will be appended
 * with the appropriate flags for orientation and subunit status.
 *
 *****************************************************************************/
struct mcell_species *
mcell_add_to_species_list(mcell_symbol *species_ptr, bool is_oriented,
                          int orientation, bool is_subunit,
                          struct mcell_species *species_list) {
  struct mcell_species *species = (struct mcell_species *)CHECKED_MALLOC_STRUCT(
      struct mcell_species, "species list");
  if (species == NULL) {
    return NULL;
  }

  species->next = NULL;
  species->mol_type = species_ptr;
  species->orient_set = 1 ? is_oriented : 0;
  species->orient = orientation;
  species->is_subunit = 1 ? is_subunit : 0;

  if (species_list != NULL) {
    species->next = species_list;
  }

  return species;
}

/*****************************************************************************
 *
 * mcell_delete_species_list frees all memory associated with a list of
 * mcell_species
 *
 *****************************************************************************/
void mcell_delete_species_list(struct mcell_species *species) {
  struct mcell_species *tmp = species;
  while (species) {
    tmp = species->next;
    free(species);
    species = tmp;
  }
}

/*****************************************************************************
 *
 * mcell_create_reaction_rates list creates a struct reaction_rates used
 * for creating reactions from a forward and backward reaction rate.
 * The backward rate is only needed for catalytic arrow and should be
 * RATE_UNUSED otherwise.
 *
 *****************************************************************************/
struct reaction_rates mcell_create_reaction_rates(int forwardRateType,
                                                  int forwardRateConstant,
                                                  int backwardRateType,
                                                  int backwardRateConstant) {
  struct reaction_rate forwardRate;
  forwardRate.rate_type = forwardRateType;
  forwardRate.v.rate_constant = forwardRateConstant;

  struct reaction_rate backwardRate;
  backwardRate.rate_type = backwardRateType;
  backwardRate.v.rate_constant = backwardRateConstant;

  struct reaction_rates rates = { forwardRate, backwardRate };

  return rates;
}

/*************************************************************************
 mcell_copysort_numeric_list:
    Copy and sort a num_expr_list in ascending numeric order.

 In:  parse_state:  parser state
      head:  the list to sort
 Out: list is sorted
*************************************************************************/
struct num_expr_list * mcell_copysort_numeric_list(struct num_expr_list *head) {
  struct num_expr_list_head new_head;
  if (mcell_generate_range_singleton(&new_head, head->value))
    return NULL;

  head = head->next;
  while (head != NULL) {
    struct num_expr_list *insert_pt, **prev;
    for (insert_pt = new_head.value_head, prev = &new_head.value_head;
         insert_pt != NULL;
         prev = &insert_pt->next, insert_pt = insert_pt->next) {
      if (insert_pt->value >= head->value)
        break;
    }

    struct num_expr_list *new_item =
        CHECKED_MALLOC_STRUCT(struct num_expr_list, "numeric array");
    if (new_item == NULL) {
      mcell_free_numeric_list(new_head.value_head);
      return NULL;
    }

    new_item->next = insert_pt;
    new_item->value = head->value;
    *prev = new_item;
    if (insert_pt == NULL)
      new_head.value_tail = new_item;
    head = head->next;
  }

  return new_head.value_head;
}

/*************************************************************************
 mcell_sort_numeric_list:
    Sort a num_expr_list in ascending numeric order.  N.B. This uses bubble
    sort, which is O(n^2).  Don't use it if you expect your list to be very
    long.  The list is sorted in-place.

 In:  head:  the list to sort
 Out: list is sorted
*************************************************************************/
void mcell_sort_numeric_list(struct num_expr_list *head) {
  struct num_expr_list *curr, *next;
  int done = 0;
  while (!done) {
    done = 1;
    curr = head;
    while (curr != NULL) {
      next = curr->next;
      if (next != NULL) {
        if (curr->value > next->value) {
          done = 0;
          swap_double(&curr->value, &next->value);
        }
      }
      curr = next;
    }
  }
}


/*************************************************************************
 mcell_free_numeric_list:
    Free a num_expr_list.

 In:  nel:  the list to free
 Out: all elements are freed
*************************************************************************/
void mcell_free_numeric_list(struct num_expr_list *nel) {
  while (nel != NULL) {
    struct num_expr_list *n = nel;
    nel = nel->next;
    free(n);
  }
}

/*************************************************************************
 mcell_generate_range:
    Generate a num_expr_list containing the numeric values from start to end,
    incrementing by step.

 In:  state: the simulation state
      list:  destination to receive list of values
      start: start of range
      end:   end of range
      step:  increment
 Out: 0 on success, 1 on failure.  On success, list is filled in.
*************************************************************************/
MCELL_STATUS mcell_generate_range(struct num_expr_list_head *list,
                                  double start, double end, double step) {
  list->value_head = NULL;
  list->value_tail = NULL;
  list->value_count = 0;
  list->shared = 0;

  if (step > 0) {
    /* JW 2008-03-31: In the guard on the loop below, it seems to me that
     * the third condition is redundant with the second.
     */
    for (double tmp_dbl = start;
         tmp_dbl < end || !distinguishable(tmp_dbl, end, EPS_C) ||
             fabs(end - tmp_dbl) <= EPS_C;
         tmp_dbl += step) {
      if (advance_range(list, tmp_dbl))
        return MCELL_FAIL;
    }
  } else /* if (step < 0) */
  {
    /* JW 2008-03-31: In the guard on the loop below, it seems to me that
     * the third condition is redundant with the second.
     */
    for (double tmp_dbl = start;
         tmp_dbl > end || !distinguishable(tmp_dbl, end, EPS_C) ||
             fabs(end - tmp_dbl) <= EPS_C;
         tmp_dbl += step) {
      if (advance_range(list, tmp_dbl))
        return MCELL_FAIL;
    }
  }
  return MCELL_SUCCESS;
}

// Maybe move this somewhere else
int advance_range(struct num_expr_list_head *list, double tmp_dbl) {
  struct num_expr_list *nel;
  nel = CHECKED_MALLOC_STRUCT(struct num_expr_list, "numeric list");
  if (nel == NULL) {
    mcell_free_numeric_list(list->value_head);
    list->value_head = list->value_tail = NULL;
    return MCELL_FAIL;
  }
  nel->value = tmp_dbl;
  nel->next = NULL;

  ++list->value_count;
  if (list->value_tail != NULL)
    list->value_tail->next = nel;
  else
    list->value_head = nel;
  list->value_tail = nel;
  return MCELL_SUCCESS;
}

/*************************************************************************
 mcell_generate_range_singleton:
    Generate a numeric list containing a single value.

 In:  lh:   list to receive value
      value: value for list
 Out: 0 on success, 1 on failure
*************************************************************************/
int mcell_generate_range_singleton(struct num_expr_list_head *lh, double value) {

  struct num_expr_list *nel =
      CHECKED_MALLOC_STRUCT(struct num_expr_list, "numeric array");
  if (nel == NULL)
    return 1;
  lh->value_head = lh->value_tail = nel;
  lh->value_count = 1;
  lh->shared = 0;
  lh->value_head->value = value;
  lh->value_head->next = NULL;
  return 0;
}

/************************************************************************
 swap_double:
 In:  x, y: doubles to swap
 Out: Swaps references to two double values.
 ***********************************************************************/
void swap_double(double *x, double *y) {
  double temp;

  temp = *x;
  *x = *y;
  *y = temp;
}
