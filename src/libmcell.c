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

#include <signal.h>
#include <string.h>
#include <time.h>

#include "argparse.h"
#include "chkpt.h"
#include "count_util.h"
#include "init.h"
#include "libmcell.h"
#include "logging.h"
#include "mem_util.h"
#include "react_output.h"
#include "react_util.h"
#include "sym_table.h"
#include "version_info.h"
#include "create_species.h"
#include "create_object.h"
#include "create_geometry.h"


/* simple wrapper for executing the supplied function call. In case
 * of an error returns with MCELL_FAIL and prints out error_message */
#define CHECKED_CALL(function, error_message) {\
   if (function) {\
     mcell_log(error_message);\
     return MCELL_FAIL;\
   }\
 }



/* declaration of static functions */
static int install_usr_signal_handlers(void);

struct output_column* get_counter_trigger_column(MCELL_STATE* state, 
    const char *counter_name, int column_id);


/************************************************************************
 * 
 * function for initializing the main mcell simulator. MCELL_STATE 
 * keeps track of the state of the simulation.
 *
 * Returns NULL on error and a pointer to MCELL_STATE otherwise
 *
 ************************************************************************/
MCELL_STATE* 
mcell_create() 
{
  // signal handlers
  if (install_usr_signal_handlers()) 
  {
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

  state->procnum=0;
  state->iterations=INT_MIN; /* indicates iterations not set */
  state->chkpt_infile = NULL;
  state->chkpt_outfile = NULL;
  state->chkpt_init = 1;
  state->log_freq = ULONG_MAX; /* Indicates that this value has not been set by user */
  state->seed_seq = 1;
  state->with_checks_flag = 1;

  time_t begin_time_of_day;
  time(&begin_time_of_day);
  state->begin_timestamp = begin_time_of_day;
  state->initialization_state = "initializing";

  if (!(state->var_sym_table = init_symtab(1024)))
  {
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
mcell_init_state(MCELL_STATE* state)
{
  CHECKED_CALL(init_notifications(state), 
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
mcell_parse_mdl(MCELL_STATE* state)
{
  return parse_input(state);
}



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
mcell_init_simulation(MCELL_STATE* state) 
{
  CHECKED_CALL(init_species(state), "Error initializing species.");

  if (state->notify->progress_report != NOTIFY_NONE)
    mcell_log("Creating geometry (this may take some time)");

  CHECKED_CALL(init_geom(state), "Error initializing geometry.");
  CHECKED_CALL(init_partitions(state), "Error initializing partitions.");
  CHECKED_CALL(init_vertices_walls(state), "Error initializing vertices and walls.");
  CHECKED_CALL(init_regions(state), "Error initializing regions.");

  if (state->place_waypoints_flag)
  {
    CHECKED_CALL(place_waypoints(state), "Error while placing waypoints.");
  }

  if (state->with_checks_flag)
  {
    CHECKED_CALL(check_for_overlapped_walls(state->n_subvols, state->subvol),
      "Error while checking for overlapped walls.");
  }

  CHECKED_CALL(init_effectors(state), "Error while placing effectors on regions.");
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
mcell_read_checkpoint(MCELL_STATE* state)
{
  if (state->chkpt_infile)
  {
    CHECKED_CALL(load_checkpoint(state), "Error while loading previous checkpoint.");

    long long exec_iterations;
    CHECKED_CALL(init_checkpoint_state(state, &exec_iterations), 
        "Error while initializing checkpoint.");

    /* XXX This is a hack to be backward compatible with the previous
     * MCell behaviour. Basically, as soon as exec_iterations <= 0 
     * MCell will stop and we emulate this by returning 1 even though
     * this is not an error (as implied by returning 1). */
    if (exec_iterations <= 0) 
    {
      mem_dump_stats(mcell_get_log_file());
      return MCELL_FAIL;
    }
  }
  else 
  {
    state->chkpt_seq_num=1;
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
mcell_init_output(MCELL_STATE* state)
{
  CHECKED_CALL(init_viz_data(state), "Error while initializing viz data.");
  CHECKED_CALL(init_reaction_data(state), "Error while initializing reaction data.");
  CHECKED_CALL(init_timers(state), "Error initializing the simulation timers.");

  // signal successful end of simulation
  state->initialization_state = NULL;

  return MCELL_SUCCESS;
}



/************************************************************************
 * 
 * function for retrieving the current value of a given count
 * expression
 *
 * The call expects:
 *
 * - MCELL_STATE
 * - counter_name: a string containing the name of the count statement to 
 *   be retrieved. Currently, the name is identical to the full path to which 
 *   the corresponding reaction output will be written but this may change
 *   in the future
 * - column: int describing the column to be retrieved
 * - count_data: a *double which will receive the actual value
 * - count_data_type: a *count_type_t which will receive the type of the 
 *   data (for casting of count_data)
 *
 * NOTE: This function can be called anytime after the 
 *       REACTION_DATA_OUTPUT has been either parsed or
 *       set up with API calls.
 *
 * Returns 1 on error and 0 on success 
 *
 ************************************************************************/
MCELL_STATUS
mcell_get_counter_value(MCELL_STATE* state, const char *counter_name,
    int column_id, double *count_data, enum count_type_t *count_data_type)
{
  struct output_column *column = NULL;
  if ((column = get_counter_trigger_column(state, counter_name, column_id))
        == NULL)
  {
    return MCELL_FAIL;
  }
     
  // if we happen to encounter trigger data we bail
  if (column->data_type == COUNT_TRIG_STRUCT) 
  {
    return MCELL_FAIL;
  }

  // evaluate the expression and retrieve it
  eval_oexpr_tree(column->expr,1);
  *count_data = (double)column->expr->value;
  *count_data_type = column->data_type;

  return MCELL_SUCCESS;
}



/************************************************************************
 * 
 * function for changing the reaction rate constant of a given named
 * reaction.
 *
 * The call expects:
 *
 * - MCELL_STATE
 * - reaction name: const char* containing the name of reaction
 * - new rate: a double with the new reaction rate constant
 *
 * NOTE: This function can be called anytime after the 
 *       REACTION_DATA_OUTPUT has been either parsed or
 *       set up with API calls.
 *
 * Returns 1 on error and 0 on success 
 *
 ************************************************************************/
MCELL_STATUS
mcell_change_reaction_rate(MCELL_STATE* state, const char *reaction_name,
    double new_rate)
{
  // sanity check
  if (new_rate < 0.0) 
  {
    return MCELL_FAIL;
  }

  // retrive reaction corresponding to name if it exists
  struct rxn *rx = NULL;
  int path_id = 0;
  if (get_rxn_by_name(state->reaction_hash, state->rx_hashsize,
        reaction_name, &rx, &path_id)) {
    return MCELL_FAIL;
  }

  // now change the rate
  if (change_reaction_probability(state, rx, path_id, new_rate))
  {
    return MCELL_FAIL;
  }

  return MCELL_SUCCESS;
}



/**************************************************************************
 * What follows are API functions for adding model elements independent of the
 * parser
 **************************************************************************/


/*************************************************************************
 mcell_create_species:
    Create a new species. This uses the same helper functions as the parser,
    but is meant to be used independent of the parser.

 In: state: the simulation state
     name:  molecule name
     D:     diffusion constant
     D_ref: reference diffusion constant
     is_2d: 1 if the species is a 2D molecule, 0 if 3D
     custom_time_step: time_step for the molecule (< 0.0 for a custom space
                       step, >0.0 for custom timestep, 0.0 for default
                       timestep)
     target_only: 1 if the molecule cannot initiate reactions
     max_step_length:
 Out: Returns 0 on sucess and 1 on error 
*************************************************************************/
int
mcell_create_species(MCELL_STATE* state,
                     struct mcell_species *species)
{
  struct sym_table *sym = CHECKED_MALLOC_STRUCT(
    struct sym_table, "sym table entry");
  int error_code = new_mol_species(state, species->name, sym);
  if (error_code) {
    return error_code;
  }


  // Perhaps we should consider getting rid of D_ref. It doesn't seem to be
  // used for anything important. Need to rip it out of test suite first.
  assemble_mol_species(state, sym, species);

  error_code = ensure_rdstep_tables_built(state);
  if (error_code) {
    return error_code;
  }

  return 0;
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
mcell_set_iterations(MCELL_STATE* state, long long iterations)
{
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
mcell_set_time_step(MCELL_STATE* state, double step)
{
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
 mcell_create_poly_object:
  Create a new polygon object.

 In: state:    the simulation state
     poly_obj: all the information needed to create the polygon object (name,
               vertices, connections)
 Out: 0 on success; any other integer value is a failure.
      A mesh is created.
*************************************************************************/
MCELL_STATUS
mcell_create_poly_object(MCELL_STATE* state, struct poly_object *poly_obj)
{
  struct object *obj_ptr = start_object(
    state, poly_obj->obj_creation, poly_obj->obj_name);

  // Set the parent of the object to be the root object. Not reciprocal until
  // add_child_objects is called.
  obj_ptr->parent = state->root_object;

  // Create the actual polygon object
  new_polygon_list(
    state, obj_ptr, poly_obj->num_vert, poly_obj->vertices,
    poly_obj->num_conn, poly_obj->connections);

  // Do some clean-up. 
  if (finish_polygon_list(obj_ptr, poly_obj->obj_creation)) {
    return MCELL_FAIL;
  }

  // Set the polygon object to be a child object of the root object (not the
  // root instance object). The object still needs to be instantiated.
  add_child_objects(state->root_object, obj_ptr, obj_ptr);

  return MCELL_SUCCESS;
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
start_object(MCELL_STATE* state,
             struct object_creation *obj_creation,
             char *name)
{
  // Create new fully qualified name.
  char *new_name;
  if ((new_name = push_object_name(obj_creation, name)) == NULL)
  {
    free(name);
    return NULL;
  }

  // Create the symbol, if it doesn't exist yet.
  struct object *obj_ptr = make_new_object(state, new_name);
  if (obj_ptr == NULL)
  {
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
new_polygon_list(MCELL_STATE* state,
                 struct object *obj_ptr,
                 int n_vertices,
                 struct vertex_list *vertices,
                 int n_connections,
                 struct element_connection_list *connections)
{

  struct polygon_object *poly_obj_ptr = allocate_polygon_object("polygon list object");
  if (poly_obj_ptr == NULL) {
    goto failure;
  }

  obj_ptr->object_type = POLY_OBJ;
  obj_ptr->contents = poly_obj_ptr;

  poly_obj_ptr->n_walls = n_connections;
  poly_obj_ptr->n_verts = n_vertices;

  // Allocate and initialize removed sides bitmask
  poly_obj_ptr->side_removed = new_bit_array(poly_obj_ptr->n_walls);
  if (poly_obj_ptr->side_removed==NULL) {
    goto failure;
  }
  set_all_bits(poly_obj_ptr->side_removed, 0);

  // Keep temporarily information about vertices in the form of
  // "parsed_vertices"
  poly_obj_ptr->parsed_vertices = vertices;

  // Copy in vertices and normals
  struct vertex_list *vert_list = poly_obj_ptr->parsed_vertices;
  for (int i = 0; i < poly_obj_ptr->n_verts; i++)
  {
    // Rescale vertices coordinates
    vert_list->vertex->x *= state->r_length_unit;
    vert_list->vertex->y *= state->r_length_unit;
    vert_list->vertex->z *= state->r_length_unit;

    vert_list = vert_list->next;
  }

  // Allocate wall elements
  struct element_data *elem_data_ptr = NULL;
  if ((elem_data_ptr = CHECKED_MALLOC_ARRAY(
      struct element_data,
      poly_obj_ptr->n_walls,
      "polygon list object walls")) == NULL) {
    goto failure;
  }
  poly_obj_ptr->element = elem_data_ptr;

  // Copy in wall elements 
  for (int i = 0; i<poly_obj_ptr->n_walls; i++)
  {
    if (connections->n_verts != 3)
    {
      //mdlerror(parse_state, "All polygons must have three vertices.");
      goto failure;
    }

    struct element_connection_list *elem_conn_list_temp = connections;
    memcpy(elem_data_ptr[i].vertex_index, connections->indices, 3*sizeof(int));
    connections = connections->next;
    free(elem_conn_list_temp->indices);
    free(elem_conn_list_temp);
  }

  // Create object default region on polygon list object: 
  struct region *reg_ptr = NULL;
  if ((reg_ptr = create_region(state, obj_ptr, "ALL")) == NULL) {
    goto failure;
  }
  if ((reg_ptr->element_list_head = new_element_list(0, poly_obj_ptr->n_walls - 1)) == NULL) {
    goto failure;
  }

  obj_ptr->n_walls = poly_obj_ptr->n_walls;
  obj_ptr->n_verts = poly_obj_ptr->n_verts;
  if (normalize_elements(reg_ptr, 0))
  {
    //mdlerror_fmt(parse_state,
    //             "Error setting up elements in default 'ALL' region in the "
    //             "polygon object '%s'.", sym->name);
    goto failure;
  }

  return poly_obj_ptr;

failure:
  free_connection_list(connections);
  free_vertex_list(vertices);
  if (poly_obj_ptr)
  {
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
int
finish_polygon_list(struct object *obj_ptr, struct object_creation *obj_creation)
{
  pop_object_name(obj_creation);
  remove_gaps_from_regions(obj_ptr);
  //no_printf(" n_verts = %d\n", mpvp->current_polygon->n_verts);
  //no_printf(" n_walls = %d\n", mpvp->current_polygon->n_walls);
  if (check_degenerate_polygon_list(obj_ptr)) {
    return 1;
  }
  return 0;
}



/**************************************************************************
 *
 * what follows are helper functions *not* part of the actual API.
 *
 * XXX: These functions should absolutely not be called from client
 *      code and will be removed eventually.
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
static int 
install_usr_signal_handlers(void)
{
#ifndef _WIN32 /* fixme: Windows does not support USR signals */
  struct sigaction sa, saPrev;
  sa.sa_sigaction = NULL;
  sa.sa_handler = &chkpt_signal_handler;
  sa.sa_flags = SA_RESTART;
  sigfillset(&sa.sa_mask);

  if (sigaction(SIGUSR1, &sa, &saPrev) != 0)
  {
    mcell_error("Failed to install USR1 signal handler.");
    return 1;
  }
  if (sigaction(SIGUSR2, &sa, &saPrev) != 0)
  {
    mcell_error("Failed to install USR2 signal handler.");
    return 1;
  }
#endif

  return 0;
}



/************************************************************************
 * 
 * helper function printing the version string
 *
 ************************************************************************/
void 
mcell_print_version()
{
  print_version(mcell_get_log_file());
}



/************************************************************************
 * 
 * helper function printing the usage information
 *
 ************************************************************************/
void 
mcell_print_usage(const char *executable_name)
{
  print_usage(mcell_get_log_file(), executable_name);
}



/************************************************************************
 * 
 * helper function printing the simulation stats
 *
 ************************************************************************/
void 
mcell_print_stats()
{
  mem_dump_stats(mcell_get_log_file());
}



/************************************************************************
 * 
 * function for printing a string
 *
 * XXX: This is a temporary hack to be able to print in mcell.c
 *      since mcell disables regular printf
 *
 ************************************************************************/
void 
mcell_print(const char *message)
{
  mcell_log("%s", message);
}


/************************************************************************
 * 
 * helper function for parsing the commandline and setting the
 * relevant parts of the state (seed #, logging, ...)
 *
 ************************************************************************/
int
mcell_argparse(int argc, char **argv, MCELL_STATE* state)
{
  return argparse_init(argc, argv, state);
}



/************************************************************************
 * 
 * helper function for retrieving the output_column corresponding
 * to a given count or trigger statement.
 *
 ************************************************************************/
struct output_column*
get_counter_trigger_column(MCELL_STATE* state, const char *counter_name,
    int column_id)
{
  // retrieve the counter for the requested counter_name
  struct sym_table *counter_sym = retrieve_sym(counter_name, 
      state->counter_by_name);
  if (counter_sym == NULL) {
    mcell_log("Failed to retrieve symbol for counter %s.", counter_name);
    return NULL;
  }
  struct output_set *counter = (struct output_set*)(counter_sym->value);
 
  // retrieve the requested column
  struct output_column *column = counter->column_head;
  int count = 0;
  while (count < column_id && column != NULL)
  {
    count++;
    column = column->next;
  }
  if (count != column_id || column == NULL)
  {
    return NULL;
  }

  return column;
}
