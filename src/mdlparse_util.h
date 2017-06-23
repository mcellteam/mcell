/******************************************************************************
 *
 * Copyright (C) 2006-2017 by
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

#pragma once

#include "vector.h"
#include "mcell_structs.h"
#include "mcell_react_out.h"
#include "mcell_reactions.h"
#include "mcell_objects.h"
#include "mdlparse_aux.h"

/* ====================================
 * Caveat lector: Most of the functions in this file (and the corresponding .c
 * file) are designed to be called from the grammar, and take care of cleaning
 * up memory allocation.  If you are in doubt about whether you can still use
 * any pointers you've passed in check the function source code to be certain.
 * In most cases, the answer is "no".
 * ====================================
 */

/* Strips the quotes from a quoted string, freeing the original. */
char *mdl_strip_quotes(char *in);

/* Concatenates two strings, freeing the original strings. */
char *mdl_strcat(char *s1, char *s2);

/* Duplicates a string. */
char *mdl_strdup(char const *s1);

/* Display a warning message about something encountered during the
 * parse process. */
void mdl_warning(struct mdlparse_vars *parse_state, char const *fmt, ...)
    PRINTF_FORMAT(2);

/* Check that the speficied file mode string is valid for an fopen statement.
 */
int mdl_valid_file_mode(struct mdlparse_vars *parse_state, char *mode);

/* Error-reporting version of log */
int mdl_expr_log(struct mdlparse_vars *parse_state, double in, double *out);

/* Error-reporting version of log10 */
int mdl_expr_log10(struct mdlparse_vars *parse_state, double in, double *out);

/* Error-reporting version of MOD */
int mdl_expr_mod(struct mdlparse_vars *parse_state, double in, double divisor,
                 double *out);

/* Error-reporting version of division operator */
int mdl_expr_div(struct mdlparse_vars *parse_state, double in, double divisor,
                 double *out);

/* Error-reporting version of exponentiation operator */
int mdl_expr_pow(struct mdlparse_vars *parse_state, double in, double exponent,
                 double *out);

/* Get a uniform random number */
double mdl_expr_rng_uniform(struct mdlparse_vars *parse_state);

/* Round a value off to n digits */
double mdl_expr_roundoff(double in, int ndigits);

/* Turn a string into a double */
int mdl_expr_string_to_double(struct mdlparse_vars *parse_state, char *str,
                              double *out);

/* Add a new file handle to the symbol table. */
struct sym_entry *mdl_new_filehandle(struct mdlparse_vars *parse_state,
                                     char *name);

/* Process an fopen statement, opening a new file handle. */
int mdl_fopen(struct mdlparse_vars *parse_state, struct sym_entry *filesym,
              char *name, char *mode);

/* Process an fclose statement, closing an existing file handle. */
int mdl_fclose(struct mdlparse_vars *parse_state, struct sym_entry *filesym);

/* Create a new double argument for a printf argument list. */
struct arg *mdl_new_printf_arg_double(double d);

/* Create a new string argument for a printf argument list. */
struct arg *mdl_new_printf_arg_string(char const *str);

/* Expands C-style escape sequences in the string. */
char *mdl_expand_string_escapes(char *in);

/* printf-like formatting of MDL arguments.  Prints to the defined err_file. */
int mdl_printf(struct mdlparse_vars *parse_state, char *fmt,
               struct arg *arg_head);

/* fprintf-like formatting of MDL arguments. */
int mdl_fprintf(struct mdlparse_vars *parse_state, struct file_stream *filep,
                char *fmt, struct arg *arg_head);

/* (expression-form) sprintf-like formatting of MDL arguments. */
char *mdl_string_format(struct mdlparse_vars *parse_state, char *fmt,
                        struct arg *arg_head);

/* sprintf-like formatting of MDL arguments. */
int mdl_sprintf(struct mdlparse_vars *parse_state, struct sym_entry *assign_var,
                char *fmt, struct arg *arg_head);

/* strtime-like formatting of current time. */
int mdl_fprint_time(struct mdlparse_vars *parse_state, struct sym_entry *filep,
                    char *fmt);

/* strtime-like formatting of current time.  Prints to the defined err_file. */
void mdl_print_time(struct mdlparse_vars *parse_state, char *fmt);

/* Generate a num_expr_list containing the numeric values from start to end,
 * incrementing by step. */
int mdl_generate_range(struct mdlparse_vars *parse_state,
                       struct num_expr_list_head *list, double start,
                       double end, double step);

/* Generate a numeric list containing a single value. */
int mdl_generate_range_singleton(struct num_expr_list_head *lh, double value);

/* Add a value to a numeric list. */
int mdl_add_range_value(struct num_expr_list_head *lh, double value);

#ifdef DEBUG
/* Display a human-readable representation of the specified array. */
void mdl_debug_dump_array(struct num_expr_list *nel);
#else
#define mdl_debug_dump_array(x) /* do nothing */
#endif

/* Create a 3-D vector from a numeric array. */
struct vector3 *mdl_point(struct mdlparse_vars *parse_state,
                          struct num_expr_list_head *vals);

/* Create a 3-D vector equal to s*[1, 1, 1] for some scalar s. */
struct vector3 *mdl_point_scalar(double val);

/* Get a named variable if it exists, or create it if it doesn't. */
struct sym_entry *mdl_get_or_create_variable(struct mdlparse_vars *parse_state,
                                             char *name);

/* Assign a "double" value to a variable, freeing any previous value. */
int mdl_assign_variable_double(struct mdlparse_vars *parse_state,
                               struct sym_entry *sym, double value);

/* Assign a string value to a variable, freeing any previous value. */
int mdl_assign_variable_string(struct mdlparse_vars *parse_state,
                               struct sym_entry *sym, char *value);

/* Assign an array value to a variable, freeing any previous value. */
int mdl_assign_variable_array(struct mdlparse_vars *parse_state,
                              struct sym_entry *sym,
                              struct num_expr_list *value);

/* Assign one variable value to another variable, freeing any previous value.
 */
int mdl_assign_variable(struct mdlparse_vars *parse_state,
                        struct sym_entry *sym, struct sym_entry *value);

/* Set all notification levels to a particular value. */
void mdl_set_all_notifications(struct volume *vol, byte notify_value);

/* Set the iteration report frequency. */
int mdl_set_iteration_report_freq(struct mdlparse_vars *parse_state,
                                  long long interval);

/* Set all warning levels to a particular value. */
void mdl_set_all_warnings(struct volume *vol, byte warning_level);

/* Set the lifetime warning threshold. */
int mdl_set_lifetime_warning_threshold(struct mdlparse_vars *parse_state,
                                       long long lifetime);

/* Set the missed reaction warning threshold. */
int mdl_set_missed_reaction_warning_threshold(struct mdlparse_vars *parse_state,
                                              double rxfrac);

/* Set the global timestep for the simulation. */
int mdl_set_time_step(struct mdlparse_vars *parse_state, double step);

/* Set the maximum timestep for the simulation. */
int mdl_set_max_time_step(struct mdlparse_vars *parse_state, double step);

/* Set the global space step for the simulation. */
int mdl_set_space_step(struct mdlparse_vars *parse_state, double step);

/* Set the number of iterations for the simulation. */
int mdl_set_num_iterations(struct mdlparse_vars *parse_state,
                           long long numiters);

/* Set the number of radial directions. */
int mdl_set_num_radial_directions(struct mdlparse_vars *parse_state,
                                  int numdirs);

/* Set the number of radial subdivisions. */
int mdl_set_num_radial_subdivisions(struct mdlparse_vars *parse_state,
                                    int numdivs);

/* Set the interaction radius. */
int mdl_set_interaction_radius(struct mdlparse_vars *parse_state,
                               double interaction_radius);

/* Set the effector grid density. */
int mdl_set_grid_density(struct mdlparse_vars *parse_state, double density);

/* Schedule an asynchronous checkpoint. */
int mdl_set_realtime_checkpoint(struct mdlparse_vars *parse_state,
                                long duration, int cont_after_cp);

/* Set the input checkpoint file to use. */
int mdl_set_checkpoint_infile(struct mdlparse_vars *parse_state, char *name);

/* Set the output checkpoint file to use. */
int mdl_set_checkpoint_outfile(struct mdlparse_vars *parse_state, char *name);

/* Set if intermediate checkpoint files should be kept */
int mdl_keep_checkpoint_files(struct mdlparse_vars *parse_state, int keepFiles);

/* Set the number of iterations between checkpoints. */
int mdl_set_checkpoint_interval(struct mdlparse_vars *parse_state,
                                long long iters, int continueAfterChkpt);

/* Set the partitioning in a particular dimension. */
int mdl_set_partition(struct mdlparse_vars *parse_state, int dim,
                      struct num_expr_list_head *head);

/* Starts creating a new object.  This function has side effects, and must be
 * paired with a call to mdl_finish_object(parse_state) */
struct sym_entry *mdl_start_object(struct mdlparse_vars *parse_state,
                                   char *name);

/* "Finishes" a new object, undoing the state changes that occurred when the
 * object was "started".
 */
void mdl_finish_object(struct mdlparse_vars *parse_state);

/* Adds the first element to an empty object list. */
void mdl_object_list_singleton(struct object_list *head, struct object *objp);

/* Adds an element to an object list. */
void mdl_add_object_to_list(struct object_list *head, struct object *objp);

/* Find an existing object or print an error message if the object isn't found.
 */
struct sym_entry *mdl_existing_object(struct mdlparse_vars *parse_state,
                                      char *name);

/* Find a list of existing objects matching a particular wildcard. */
struct sym_table_list *
mdl_existing_objects_wildcard(struct mdlparse_vars *parse_state,
                              char *wildcard);

/* Find an existing region or print an error message if it isn't found. */
struct sym_entry *mdl_existing_region(struct mdlparse_vars *parse_state,
                                      struct sym_entry *obj_symp, char *name);

/* Find an existing molecule species, or print an error message if it isn't
 * found. */
struct sym_entry *mdl_existing_molecule(struct mdlparse_vars *parse_state,
                                        char *name);

/* Turn a single symbol into a singleton symbol list. */
struct sym_table_list *
mdl_singleton_symbol_list(struct mdlparse_vars *parse_state,
                          struct sym_entry *sym);

/* Find an existing molecule species, and return it in a singleton list, or
 * print an error message if it isn't found. */
struct sym_table_list *
mdl_existing_molecule_list(struct mdlparse_vars *parse_state, char *name);

/* Find a list of all molecule species matching the specified wildcard.  Print
 * an error message if it doesn't match any. */
struct sym_table_list *
mdl_existing_molecules_wildcard(struct mdlparse_vars *parse_state,
                                char *wildcard);

/* Find an existing surface molecule species, or print an error message if it
 * isn't found, or isn't a surface molecule. */
struct sym_entry *
mdl_existing_surface_molecule(struct mdlparse_vars *parse_state, char *name);

/* Find an existing surface class species, or print an error message if it
 * isn't found, or isn't a surface class. */
struct sym_entry *mdl_existing_surface_class(struct mdlparse_vars *parse_state,
                                             char *name);

/* Find a named variable if it exists, or print an error if it does not. */
struct sym_entry *mdl_existing_variable(struct mdlparse_vars *parse_state,
                                        char *name);

/* Find an existing array symbol, or print an error message if it isn't found.
 */
struct sym_entry *mdl_existing_array(struct mdlparse_vars *parse_state,
                                     char *name);

/* Find a named numeric variable if it exists, or print an error message if it
 * isn't found. */
struct sym_entry *mdl_existing_double(struct mdlparse_vars *parse_state,
                                      char *name);

/* Find a named string variable if it exists, or print an error message if it
 * isn't found. */
struct sym_entry *mdl_existing_string(struct mdlparse_vars *parse_state,
                                      char *name);

/* Get a named numeric or array variable if it exists.  Print an error message
 * if it isn't found. */
struct sym_entry *mdl_existing_num_or_array(struct mdlparse_vars *parse_state,
                                            char *name);

/* Find an existing named reaction pathway or molecule, or print an error
 * message if it isn't found. */
struct sym_entry *
mdl_existing_rxn_pathname_or_molecule(struct mdlparse_vars *parse_state,
                                      char *name);

/* Find an existing reaction pathway or release pattern, or print an error
 * message if it isn't found, or if the name could refer to either a release
 * pattern or a reaction pathway.
 */
struct sym_entry *
mdl_existing_release_pattern_or_rxn_pathname(struct mdlparse_vars *parse_state,
                                             char *name);

/* Find an existing file stream or print an error message if it isn't found. */
struct sym_entry *mdl_existing_file_stream(struct mdlparse_vars *parse_state,
                                           char *name);

/* Find all mesh objects (polygons and boxes) that match a given wildcard. */
struct sym_table_list *mdl_meshes_by_wildcard(struct mdlparse_vars *parse_state,
                                              char *wildcard);

/* Apply a rotation to the given transformation matrix. */
int mdl_transform_rotate(struct mdlparse_vars *parse_state, double (*mat)[4],
                         struct vector3 *axis, double angle);

/* Deep copy an object. */
int mdl_deep_copy_object(struct mdlparse_vars *parse_state,
                         struct object *dst_obj, struct object *src_obj);

/* Prepare a region description for use in the simulation by creating a
 * membership bitmask on the region object. */
int mdl_normalize_elements(struct mdlparse_vars *parse_state,
                           struct region *reg, int existing);

/* Finalizes the polygonal structure of the box, normalizing all regions. */
int mdl_triangulate_box_object(struct mdlparse_vars *parse_state,
                               struct sym_entry *box_sym,
                               struct polygon_object *pop,
                               double box_aspect_ratio);

/* Check that the specified diffusion constant is valid, correcting it if
 * appropriate. */
int mdl_check_diffusion_constant(struct mdlparse_vars *parse_state, double *d);

/* Finish the creation of a series of molecules, undoing any state changes we
 * made during the creation of the molecules. */
void mdl_print_species_summary(MCELL_STATE *state,
                               struct mcell_species_spec *species);
void mdl_print_species_summaries(struct volume *state,
                                 struct parse_mcell_species_list_item *mols);

int mdl_add_to_species_list(struct parse_mcell_species_list *list,
                            struct mcell_species_spec *spec);

/* Start parsing the innards of a release site. */
int mdl_start_release_site(struct mdlparse_vars *parse_state,
                           struct sym_entry *symp, int shape);

/* Finish parsing the innards of a release site. */
struct object *mdl_finish_release_site(struct mdlparse_vars *parse_state,
                                       struct sym_entry *symp);

/* Validate a release site. */
int mdl_is_release_site_valid(struct mdlparse_vars *parse_state,
                              struct release_site_obj *rsop);

/* Set the geometry for a particular release site to be a region expression. */
int mdl_set_release_site_geometry_region(struct mdlparse_vars *parse_state,
                                         struct release_site_obj *rsop,
                                         struct object *objp,
                                         struct release_evaluator *re);

/* Set the geometry for a particular release site to be an entire object. */
int
mdl_set_release_site_geometry_object(struct mdlparse_vars *parse_state,
                                     struct release_site_obj *rel_site_obj_ptr,
                                     struct object *obj_ptr);

/* Set the molecule to be released from this release site. */
int mdl_set_release_site_molecule(struct mdlparse_vars *parse_state,
                                  struct release_site_obj *rsop,
                                  struct mcell_species *mol_type);

/* Set the diameter of a release site. */
int mdl_set_release_site_diameter(struct mdlparse_vars *parse_state,
                                  struct release_site_obj *rsop, double diam);

/* Set the diameter of the release site along the X, Y, and Z axes. */
int mdl_set_release_site_diameter_array(struct mdlparse_vars *parse_state,
                                        struct release_site_obj *rsop,
                                        int n_diams,
                                        struct num_expr_list *diams,
                                        double factor);

/* Set the diameters of the release site along the X, Y, and Z axes from a
 * variable, either scalar or vector. */
int mdl_set_release_site_diameter_var(struct mdlparse_vars *parse_state,
                                      struct release_site_obj *rsop,
                                      double factor, struct sym_entry *symp);

int mdl_set_release_site_periodic_box(struct mdlparse_vars *parse_state,
                                      struct release_site_obj *rel_site_obj_ptr,
                                      struct vector3 *periodic_box);

/* Set the release probability for a release site. */
int mdl_set_release_site_probability(struct mdlparse_vars *parse_state,
                                     struct release_site_obj *rsop,
                                     double prob);

/* Set the release pattern to be used by a particular release site. */
int mdl_set_release_site_pattern(struct mdlparse_vars *parse_state,
                                 struct release_site_obj *rsop,
                                 struct sym_entry *pattern);

/* Set the molecule positions for a LIST release. */
int mdl_set_release_site_molecule_positions(
    struct mdlparse_vars *parse_state, struct release_site_obj *rsop,
    struct release_single_molecule_list *list);

/* Create a mew single molecule release position for a LIST release site. */
struct release_single_molecule *
mdl_new_release_single_molecule(struct mdlparse_vars *parse_state,
                                struct mcell_species *mol_type,
                                struct vector3 *pos);

/* Set a release quantity from this release site based on a fixed concentration
 * within the release-site's area. */
int mdl_set_release_site_concentration(struct mdlparse_vars *parse_state,
                                       struct release_site_obj *rsop,
                                       double conc);

/* Set an item to be the sole element of a vertex list. */
void mdl_vertex_list_singleton(struct vertex_list_head *head,
                               struct vertex_list *item);

/* Append a vertex to a list. */
void mdl_add_vertex_to_list(struct vertex_list_head *head,
                            struct vertex_list *item);

/* Allocate an item for a vertex list. */
struct vertex_list *mdl_new_vertex_list_item(struct vector3 *vertex);

/* Set an item to be the sole element of an element connection list. */
void
mdl_element_connection_list_singleton(struct element_connection_list_head *head,
                                      struct element_connection_list *item);

/* Append an element connection to a list. */
void
mdl_add_element_connection_to_list(struct element_connection_list_head *head,
                                   struct element_connection_list *item);

/* Create an element connection (essentially a triplet of vertex indices). */
struct element_connection_list *
mdl_new_element_connection(struct mdlparse_vars *parse_state,
                           struct num_expr_list_head *indices);

/* Create a tetrahedral element connection (essentially a quadruplet of vertex
 * indices). */
struct element_connection_list *
mdl_new_tet_element_connection(struct mdlparse_vars *parse_state,
                               struct num_expr_list_head *indices);

/* Create a new polygon list object. */
struct object *
mdl_new_polygon_list(struct mdlparse_vars *parse_state, char *obj_name,
                     int n_vertices, struct vertex_list *vertices,
                     int n_connections,
                     struct element_connection_list *connections);

/* Finalize the polygon list, cleaning up any state updates that were made when
 * we started creating the polygon. */
int mdl_finish_polygon_list(struct mdlparse_vars *parse_state,
                            struct object *obj_ptr);

/* Create a new voxel list object. */
struct voxel_object *
mdl_new_voxel_list(struct mdlparse_vars *parse_state, struct sym_entry *sym,
                   int n_vertices, struct vertex_list *vertices,
                   int n_connections,
                   struct element_connection_list *connections);

struct polygon_object *mdl_create_periodic_box(
    struct mdlparse_vars *parse_state,
    struct vector3 *llf,
    struct vector3 *urb,
    bool isPeriodicX,
    bool isPeriodicY,
    bool isPeriodicZ);

int mdl_finish_periodic_box(struct mdlparse_vars *parse_state);

/* Create a new box object, with particular corners. */
struct polygon_object *mdl_new_box_object(struct mdlparse_vars *parse_state,
                                          struct sym_entry *sym,
                                          struct vector3 *llf,
                                          struct vector3 *urb);

/* Finalize the box object, cleaning up any state updates that were made when
 * we started creating the box. */
int mdl_finish_box_object(struct mdlparse_vars *parse_state,
                          struct sym_entry *symp);

/* Create a named region on an object. */
struct region *mdl_create_region(struct mdlparse_vars *parse_state,
                                 struct object *objp, char *name);

/* Get a region on an object, creating it if it does not exist yet. */
struct region *mdl_get_region(struct mdlparse_vars *parse_state,
                              struct object *objp, char *name);

/* Begin construction of a region on an existing object. */
int mdl_start_existing_obj_region_def(struct mdlparse_vars *parse_state,
                                      struct sym_entry *obj_symp);

/* Append elements to an element list. */
void mdl_add_elements_to_list(struct element_list_head *list,
                              struct element_list *head,
                              struct element_list *tail);

/* Marks elements as being excluded, rather than included (the default). */
void mdl_set_elements_to_exclude(struct element_list *els);

/* Create a new element list for a region description. */
struct element_list *mdl_new_element_list(struct mdlparse_vars *parse_state,
                                          unsigned int begin, unsigned int end);

/* Create a new element list for a region description based on a side name. */
struct element_list *mdl_new_element_side(struct mdlparse_vars *parse_state,
                                          unsigned int side);

/* Create a new element list for a "previous region" include/exclude statement.
 */
struct element_list *mdl_new_element_previous_region(
    struct mdlparse_vars *parse_state, struct object *objp,
    struct region *rp_container, char *name_region_referent, int exclude);

/* Allocate a new region element list item for an include/exclude PATCH
 * statement. */
struct element_list *mdl_new_element_patch(struct mdlparse_vars *parse_state,
                                           struct polygon_object *poly,
                                           struct vector3 *llf,
                                           struct vector3 *urb, int exclude);

/* Set the elements for a region, normalizing the region if it's on a polygon
 * list object. */
int mdl_set_region_elements(struct mdlparse_vars *parse_state,
                            struct region *rgn, struct element_list *elements,
                            int normalize_now);

/* Create a new named reaction pathway name structure. */
struct sym_entry *mdl_new_rxn_pathname(struct mdlparse_vars *parse_state,
                                       char *name);

/* Adds an effector (or list of effectors) to a region.  These effectors will
 * be placed on the surface at initialization time. */
void mdl_add_surf_mol_to_region(struct region *rgn, struct sm_dat_list *lst);

/* Set the surface class of this region, possibly inheriting the viz_value.  */
void mdl_set_region_surface_class(struct mdlparse_vars *parse_state,
                                  struct region *rgn, struct sym_entry *scsymp);

/****************************************************************
 * Reaction output
 ***************************************************************/

/* Finalizes a reaction data output block, checking for errors, and allocating
 * the output buffer. */
int mdl_output_block_finalize(struct mdlparse_vars *parse_state,
                              struct output_block *obp);

/* Populate an output set. */
struct output_set *mdl_populate_output_set(struct mdlparse_vars *parse_state,
                                           char *comment, int exact_time,
                                           struct output_column *col_head,
                                           int file_flags, char *outfile_name);

/* Construct and add an output block to the world. */
int mdl_add_reaction_output_block_to_world(struct mdlparse_vars *parse_state,
                                           int buffer_size,
                                           struct output_times_inlist *otimes,
                                           struct output_set_list *osets);

/* Joins two subtrees into a reaction data output expression tree, with a
 * specified operation. */
struct output_expression *mdl_join_oexpr_tree(struct mdlparse_vars *parse_state,
                                              struct output_expression *left,
                                              struct output_expression *right,
                                              char oper);

/* Converts an output expression tree generated from a wildcard into a
 * summation expression tree. */
struct output_expression *mdl_sum_oexpr(struct output_expression *expr);

/* Creates a constant output expression for reaction data output. */
struct output_expression *
mdl_new_oexpr_constant(struct mdlparse_vars *parse_state, double value);

/* Generates a reaction data output expression from the first count syntax form
 * (simple molecule, unquoted, no orientation). */
struct output_expression *mdl_count_syntax_1(struct mdlparse_vars *parse_state,
                                             struct sym_entry *what,
                                             struct sym_entry *where,
                                             int hit_spec, int count_flags);

/* Generates a reaction data output expression from the second count syntax
 * form (simple molecule, unquoted, orientation in braces). */
struct output_expression *mdl_count_syntax_2(struct mdlparse_vars *parse_state,
                                             struct sym_entry *mol_type,
                                             short orient,
                                             struct sym_entry *where,
                                             int hit_spec, int count_flags);

/* Generates a reaction data output expression from the third count syntax form
 * (quoted string, possibly a wildcard, possibly an oriented molecule). */
struct output_expression *mdl_count_syntax_3(struct mdlparse_vars *parse_state,
                                             char *what,
                                             struct sym_entry *where,
                                             int hit_spec, int count_flags);


struct output_expression *mdl_count_syntax_periodic_1(struct mdlparse_vars *parse_state,
  struct sym_entry *what, struct sym_entry *where, struct vector3 *periodicBox,
  int hit_spec, int count_flags);

struct output_expression *mdl_count_syntax_periodic_2(
    struct mdlparse_vars *parse_state,
    struct sym_entry *mol_type,
    short orient,
    struct sym_entry *where,
    struct vector3 *periodicBox,
    int hit_spec,
    int count_flags); 

struct output_expression *mdl_count_syntax_periodic_3(
    struct mdlparse_vars *parse_state,
    char *what,
    struct sym_entry *where,
    struct vector3 *periodicBox,
    int hit_spec,
    int count_flags);

/* Prepare a single count expression for inclusion in an output set. */
int mdl_single_count_expr(struct mdlparse_vars *parse_state,
                          struct output_column_list *list,
                          struct output_expression *expr, char *custom_header);

/****************************************************************
 * Viz output
 ***************************************************************/

/* Build a new VIZ output block, containing parameters for an output set for
 * visualization. */
int mdl_new_viz_output_block(struct mdlparse_vars *parse_state);

/* Set the mode for a new VIZ output block. */
int mdl_set_viz_mode(struct viz_output_block *vizblk, int mode);

/* Set the molecule format for a new VIZ output block. */
int mdl_set_viz_molecule_format(struct mdlparse_vars *parse_state,
                                struct viz_output_block *vizblk, int format);

/* Set the filename prefix for a new VIZ output block. */
int mdl_set_viz_filename_prefix(struct mdlparse_vars *parse_state,
                                struct viz_output_block *vizblk,
                                char *filename);

/* Error-checking wrapper for a specified visualization state. */
int mdl_viz_state(struct mdlparse_vars *parse_state, int *target, double value);


/* Sets a flag on all of the listed species, requesting that they be visualized.
 */
int mdl_set_viz_include_molecules(struct mdlparse_vars *parse_state,
                                  struct viz_output_block *vizblk,
                                  struct sym_table_list *list, int viz_state);

/* Sets a flag on a viz block, requesting that all species be visualized. */
int mdl_set_viz_include_all_molecules(struct viz_output_block *vizblk,
                                      int viz_state);

/* Adds some new molecule output frames to a list. */
int mdl_new_viz_mol_frames(struct mdlparse_vars *parse_state,
                           struct viz_output_block *vizblk,
                           struct frame_data_list_head *frames, int time_type,
                           int mol_item_type,
                           struct num_expr_list_head *timelist);

/* Build a list of times for VIZ output, one timepoint per iteration in the
 * simulation. */
int mdl_new_viz_all_times(struct mdlparse_vars *parse_state,
                          struct num_expr_list_head *list);

/* Build a list of iterations for VIZ output, one for each iteration in the
 * simulation. */
int mdl_new_viz_all_iterations(struct mdlparse_vars *parse_state,
                               struct num_expr_list_head *list);

/* Set the viz_state value for a molecular species. */
int mdl_set_molecule_viz_state(struct viz_output_block *vizblk,
                               struct species *specp, int viz_state);

/* Set the viz_state for a particular region. */
int mdl_set_region_viz_state(struct mdlparse_vars *parse_state,
                             struct viz_output_block *vizblk, struct region *rp,
                             int viz_state);

/****************************************************************
 * Volume output
 ***************************************************************/

/* Create a new volume output request. */
struct volume_output_item *mdl_new_volume_output_item(
    struct mdlparse_vars *parse_state, char *filename_prefix,
    struct species_list *molecules, struct vector3 *location,
    struct vector3 *voxel_size, struct vector3 *voxel_count,
    struct output_times *ot);

/* Create new default output timing for volume output. */
struct output_times *
mdl_new_output_times_default(struct mdlparse_vars *parse_state);

/* Create new "step" output timing for volume output. */
struct output_times *
mdl_new_output_times_step(struct mdlparse_vars *parse_state, double step);

/* Create new "iteration list" output timing for volume output. */
struct output_times *
mdl_new_output_times_iterations(struct mdlparse_vars *parse_state,
                                struct num_expr_list_head *iters);

/* Create new "time list" output timing for volume output. */
struct output_times *
mdl_new_output_times_time(struct mdlparse_vars *parse_state,
                          struct num_expr_list_head *times);

/****************************************************************
 * Release patterns
 ***************************************************************/

/* Create a new release pattern.  There must not yet be a release pattern with
 * the given name. */
struct sym_entry *mdl_new_release_pattern(struct mdlparse_vars *parse_state,
                                          char *name);

/* Fill in the details of a release pattern. */
int mdl_set_release_pattern(struct mdlparse_vars *parse_state,
                            struct sym_entry *rpat_sym,
                            struct release_pattern *rpat_data);

/****************************************************************
 * Molecules
 ***************************************************************/

/* Create a new species.  There must not yet be a molecule or named reaction
 * pathway with the supplied name. */
struct sym_entry *mdl_new_mol_species(struct mdlparse_vars *parse_state,
                                      char *name);

/* Assemble a molecule species from its component pieces. */
struct mcell_species_spec *mdl_create_species(struct mdlparse_vars *parse_state,
                                              char *name, double D, int is_2d,
                                              double custom_time_step,
                                              int target_only,
                                              double max_step_length);

/****************************************************************
 * Reactions, surface classes
 ***************************************************************/

/* Check whether the reaction rate is valid. */
int mdl_valid_rate(struct mdlparse_vars *parse_state,
                   struct reaction_rate *rate);

/* Set the reactant/product list to contain a single item. */
int mdl_reaction_player_singleton(struct mdlparse_vars *parse_state,
                                  struct mcell_species_list *list,
                                  struct mcell_species *spec);

/* Add a single item to a reactant/product player list. */
int mdl_add_reaction_player(struct mdlparse_vars *parse_state,
                            struct mcell_species_list *list,
                            struct mcell_species *spec);

/* Set a reaction rate from a variable. */
int mdl_reaction_rate_from_var(struct mdlparse_vars *parse_state,
                               struct reaction_rate *rate,
                               struct sym_entry *symp);

/* Assemble a standard reaction from its component parts. */
struct mdlparse_vars *mdl_assemble_reaction(struct mdlparse_vars *parse_state,
                                            struct mcell_species *reactants,
                                            struct mcell_species *surface_class,
                                            struct reaction_arrow *react_arrow,
                                            struct mcell_species *products,
                                            struct reaction_rates *rate,
                                            struct sym_entry *pathname);

/* Assemble a surface reaction from its component parts. */
struct mdlparse_vars *
mdl_assemble_surface_reaction(struct mdlparse_vars *parse_state,
                              int reaction_type, struct species *surface_class,
                              struct sym_entry *reactant_sym, short orient);

/* Assemble a concentration clamp reaction from its component parts. */
struct mdlparse_vars *mdl_assemble_concentration_clamp_reaction(
    struct mdlparse_vars *parse_state, struct species *surface_class,
    struct sym_entry *mol_sym, short orient, double conc);

/* Start a surface class declaration. */
void mdl_start_surface_class(struct mdlparse_vars *parse_state,
                             struct sym_entry *symp);

/* Finish a surface class declaration.  Undoes side effects from
 * mdl_start_surface_class. */
void mdl_finish_surface_class(struct mdlparse_vars *parse_state);

/* Create a new effector data for surface molecule initialization. */
struct sm_dat *mdl_new_surf_mol_data(struct mdlparse_vars *parse_state,
                                     struct mcell_species *eff_info,
                                     double quant);

int warn_about_high_rates(struct mdlparse_vars *parse_state, FILE *warn_file,
                          int rate_warn, int print_once);

void alphabetize_pathway(struct pathway *path, struct rxn *reaction);

void check_duplicate_special_reactions(struct pathway *path);

int set_product_geometries(struct pathway *path, struct rxn *rx,
                           struct product *prod);

int scale_probabilities(struct pathway *path, struct rxn *rx,
                        struct mdlparse_vars *parse_state, double pb_factor);

void add_surface_reaction_flags(struct mdlparse_vars *parse_state);

void free_vertex_list(struct vertex_list *vlp);

/**********************************************************************
 ***  helper functions for release sites creation
 **********************************************************************/

// Adds a release molecule descriptor to a list.
void
add_release_single_molecule_to_list(struct release_single_molecule_list *list,
                                    struct release_single_molecule *mol);

// Populates a list with a single LIST release molecule descriptor.
void
release_single_molecule_singleton(struct release_single_molecule_list *list,
                                  struct release_single_molecule *mol);

/* Set a release quantity from this release site based on a fixed density
 * within the release-site's area. */
int set_release_site_density(struct release_site_obj *rel_site_obj_ptr,
                             double dens);

/* Set a release quantity from this release site based on a fixed concentration
 * in a sphere of a gaussian-distributed diameter with a particular mean and
 * std. deviation. */
void set_release_site_volume_dependent_number(
    struct release_site_obj *rel_site_obj_ptr, double mean, double stdev,
    double conc);

/****************************************************************************
 *** helper function for object creation
 ****************************************************************************/
void transform_translate(MCELL_STATE *state, double (*mat)[4],
                         struct vector3 *xlat);

void transform_scale(double (*mat)[4], struct vector3 *scale);

int transform_rotate(double (*mat)[4], struct vector3 *axis, double angle);

void check_regions(struct object *rootInstance, struct object *child_head);

int finish_polygon_list(struct object *obj_ptr,
                        struct object_creation *obj_creation);

struct object *start_object(MCELL_STATE *state,
                            struct object_creation *obj_creation,
                            char *name,
                            int *error_code);
