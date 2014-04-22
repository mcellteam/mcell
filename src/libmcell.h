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

#ifndef LIBMCELL_H
#define LIBMCELL_H

#include <stdbool.h>

#include "config.h"

#include "mcell_engine.h"
#include "mcell_structs.h"


/**********************************************************************
 * type declarations
 **********************************************************************/
#define REGULAR_ARROW 0x00
#define ARROW_BIDIRECTIONAL 0x01
#define ARROW_CATALYTIC     0x02

/* status of libMCell API calls */
typedef int MCELL_STATUS;

#define MCELL_SUCCESS 0
#define MCELL_FAIL 1


typedef struct sym_table mcell_symbol;

/* state of mcell simulation */
typedef struct volume MCELL_STATE;

enum {
  RATE_UNSET    = -1,
  RATE_CONSTANT = 0,
  RATE_FILE     = 1,
  RATE_COMPLEX  = 2
};


/* Special pathway types. */
enum special_pathway_t
{
  RFLCT,      /* Special pathway: reflective surface */
  TRANSP,     /* Special pathway: transparent surface */
  SINK        /* Special pathway: absorptive surface */
};


struct mcell_species_spec
{
  char *name;
  double D;
  double D_ref;              // default is 0.0
  int is_2d;                 // 3D = 0; 2D = 1
  double custom_time_step;   // default is 0.0
  int target_only;           // default is 0
  double max_step_length;    // default is 0.0
  double space_step;
};


struct object_creation
{
  struct name_list *object_name_list;
  struct name_list *object_name_list_end;
  struct object *current_object;
};

struct poly_object
{
  char *obj_name;
  struct vertex_list *vertices;
  int num_vert;
  struct element_connection_list *connections;
  int num_conn;
};

struct reaction_def {
  struct sym_table *sym;
};

struct mcell_species
{
  struct mcell_species *next;
  struct sym_table *mol_type;
  short orient_set;
  short orient;
  short is_subunit;
};

struct mcell_species_list
{
  struct mcell_species *mol_type_head;
  struct mcell_species *mol_type_tail;
};

struct release_single_molecule_list
{
  struct release_single_molecule *rsm_head;
  struct release_single_molecule *rsm_tail;
  int rsm_count;
};

struct reaction_arrow
{
  int                  flags;
  struct mcell_species catalyst;
};

struct reaction_rate
{
  int rate_type;
  union
  {
    double              rate_constant;
    char                *rate_file;
    struct complex_rate *rate_complex;
  } v;
};

struct reaction_rates
{
  struct reaction_rate forward_rate;
  struct reaction_rate backward_rate;
};


/****************************************************************
 * setup routines
 ****************************************************************/
MCELL_STATE* mcell_create();

MCELL_STATUS mcell_init_state();

MCELL_STATUS mcell_parse_mdl(MCELL_STATE *state);

MCELL_STATUS mcell_init_simulation(MCELL_STATE *state);

MCELL_STATUS mcell_read_checkpoint(MCELL_STATE *state);

MCELL_STATUS mcell_init_output(MCELL_STATE *state);


/****************************************************************
 * routines for running simulations
 ****************************************************************/

/* this function runs the whole simulations */
MCELL_STATUS mcell_run_simulation(MCELL_STATE *state);

/* returns the recommended output frequence either based on
 * a user request in the MDL or via some heuristics */
long long mcell_determine_output_frequency(MCELL_STATE *state);

/* this function runs a single iteration of simulations */
MCELL_STATUS mcell_run_iteration(MCELL_STATE *state,
    long long output_frequency, int *restarted_from_checkpoint);

/* flush all output buffers to disk to disk after the simulation
 * run is complete */
MCELL_STATUS mcell_flush_data(MCELL_STATE *state);

/* print any warnings that were gererated during the simulation
 * run */
MCELL_STATUS mcell_print_final_warnings(MCELL_STATE *state);

/* print the final simulation statistics */
MCELL_STATUS mcell_print_final_statistics(MCELL_STATE *state);


/****************************************************************
 * API functions for adding model elements independent of the parser
 ****************************************************************/

MCELL_STATUS mcell_set_time_step(MCELL_STATE* state,
  double step);

MCELL_STATUS mcell_set_iterations(MCELL_STATE* state,
  long long iterations);

MCELL_STATUS mcell_create_species(MCELL_STATE* state,
  struct mcell_species_spec *species,
  mcell_symbol **species_ptr);

MCELL_STATUS mcell_add_reaction(MCELL_STATE* state,
  struct mcell_species *reactants,
  struct reaction_arrow *arrow,
  struct mcell_species *surf_class,
  struct mcell_species *products,
  struct sym_table *pathname,
  struct reaction_rates *rates,
  const char *rate_filename);

MCELL_STATUS mcell_add_surface_reaction(MCELL_STATE *state,
  int reaction_type,
  struct species *surface_class,
  struct sym_table *reactant_sym,
  short orient);

MCELL_STATUS mcell_add_concentration_clamp(MCELL_STATE *state,
  struct species *surface_class,
  struct sym_table *mol_sym,
  short orient,
  double conc);

/****************************************************************
 * API functions for manipulating model objects
 ****************************************************************/
MCELL_STATUS mcell_create_instance_object(MCELL_STATE *state,
  char *name,
  struct object **new_object);

MCELL_STATUS mcell_create_poly_object(MCELL_STATE *state,
  struct object *parent,
  struct poly_object *poly_obj,
  struct object **new_object);


/****************************************************************
 * routines for manipulating release sites
 ****************************************************************/

MCELL_STATUS mcell_create_geometrical_release_site(MCELL_STATE *state,
  struct object *parent, char *site_name, int shape, struct vector3 *position,
  struct vector3 *diameter, struct mcell_species *mol, double num_molecules,
  double release_prob, char *pattern_name, struct object **new_object);

MCELL_STATUS mcell_start_release_site(MCELL_STATE *state,
  struct sym_table *sym_ptr, struct object **obj);

MCELL_STATUS mcell_finish_release_site(struct sym_table *sym_ptr,
  struct object **obj);


/****************************************************************
 * routines for retrieving information
 ****************************************************************/

/* this function retrieves the current value of a column in
 * count expression counter_name */
MCELL_STATUS mcell_get_counter_value(MCELL_STATE* state,
    const char *counter_name, int column, double *count_data,
    enum count_type_t *count_data_type);



/****************************************************************
 * routines for changing the state of a running simulation
 ****************************************************************/

/* this function changes the reaction rate constant of the
 * given named reaction. The change happens instantaneously,
 * e.g. within the given iteration */
MCELL_STATUS mcell_change_reaction_rate(MCELL_STATE* state,
    const char *reaction_name, double new_rate);


/*****************************************************************
 * non API helper functions
 *
 * NOTE: These functions should *not* be called directly and are
 *       not part of the API. These functions are currently used by
 *       both libmcell and the parser and need to eventually be
 *       internalized by libmcell once the parser is fully API
 *       compliant.
 *****************************************************************/
void mcell_print_version();
void mcell_print_usage(const char *executable_name);
void mcell_print_stats();
int mcell_argparse(int argc, char **argv, MCELL_STATE *state);

// helper functions for dealing with polygon lists and objects
int finish_polygon_list(struct object *obj_ptr,
  struct object_creation *obj_creation);

struct polygon_object * new_polygon_list(
  MCELL_STATE* state,
  struct object *obj_ptr,
  int n_vertices,
  struct vertex_list *vertices,
  int n_connections,
  struct element_connection_list *connections);

struct object *make_new_object(MCELL_STATE *state,
  char *obj_name);

struct object *start_object(MCELL_STATE* state,
  struct object_creation *obj_creation,
  char *name);

/* helper functions for creating vertex and element_connection lists */
struct vertex_list* mcell_add_to_vertex_list(double x, double y, double z,
  struct vertex_list *vertices);

struct element_connection_list* mcell_add_to_connection_list(int v1, int v2,
  int v3, struct element_connection_list* elements);

/* helper functions for creating reactions */
struct mcell_species* mcell_add_to_species_list(mcell_symbol *species_ptr,
  bool is_oriented, int orientation, bool is_subunit,
  struct mcell_species *species_list);

struct reaction_rates mcell_create_reaction_rates(int forwardRateType,
  int forwardRate, int backwardRateType, int backwardRate);

void mcell_delete_species_list(struct mcell_species* species);

/* helper functions for release sites */
void set_release_site_location(MCELL_STATE *state,
  struct release_site_obj *rel_site_obj_ptr, struct vector3 *location);

/***********************************************************************
 * helper function for release sites
 ***********************************************************************/

int mcell_set_release_site_geometry_region(MCELL_STATE *state,
  struct release_site_obj *rel_site_obj_ptr,
  struct object *objp,
  struct release_evaluator *re);

/* Set a constant release quantity from this release site, in units of
 * molecules. */
void set_release_site_constant_number(struct release_site_obj *rel_site_obj_ptr,
  double num);

/* Set a gaussian-distributed release quantity from this release site, in units
 * of molecules. */
void set_release_site_gaussian_number(
  struct release_site_obj *rel_site_obj_ptr,
  double mean,
  double stdev);

// Create a new "release on region" expression term.
struct release_evaluator *new_release_region_expr_term(
  struct sym_table *my_sym);

// Set the geometry for a particular release site to be a region expression.
struct release_evaluator *new_release_region_expr_binary(
  struct release_evaluator *reL,
  struct release_evaluator *reR,
  int op);

int check_release_regions(struct release_evaluator *rel,
  struct object *parent,
  struct object *instance);

int is_release_site_valid(struct release_site_obj *rel_site_obj_ptr);

/* Set a release quantity from this release site based on a fixed concentration
 * within the release-site's area. */
int set_release_site_concentration(struct release_site_obj *rel_site_obj_ptr,
                                   double conc);


struct release_evaluator *
new_release_region_expr_term(struct sym_table *my_sym);

/* helper functions for IO */

// XXX this is a temporary hack to be able to print in mcell.c
// since mcell disables regular printf
void mcell_print(const char *message);

#endif
