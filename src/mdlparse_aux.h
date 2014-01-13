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

#ifndef MDLPARSE_AUX_H
#define MDLPARSE_AUX_H

#include "mcell_structs.h"
#include "macromolecule.h"

#define ARROW_BIDIRECTIONAL 0x01
#define ARROW_CATALYTIC     0x02

/* Flags for parser to indicate which axis we are partitioning */
enum partition_axis_t
{
  X_PARTS,    /* X-axis partitions */
  Y_PARTS,    /* Y-axis partitions */
  Z_PARTS     /* Z-axis partitions */
};

/* Special pathway types. */
enum special_pathway_t
{
  RFLCT,      /* Special pathway: reflective surface */
  TRANSP,     /* Special pathway: transparent surface */
  SINK        /* Special pathway: absorptive surface */
};

#define WILDCARD_PRESENT   0x1
#define TRIGGER_PRESENT    0x2
#define COUNT_PRESENT      0x4
#define EXPRESSION_PRESENT 0x8

/***************************************************************************
 * Miscellaneous parse-time structures
 ***************************************************************************/

struct arg
{ 
  struct arg *next;
  byte arg_type; /* DBL, STR */
  void *arg_value;
};

struct arg_list
{
  struct arg *arg_head;
  struct arg *arg_tail;
};

struct species_opt_orient
{
  struct species_opt_orient *next;
  struct sym_table *mol_type;
  short orient_set;
  short orient;
  short is_subunit;
};

struct species_opt_orient_list
{
  struct species_opt_orient *mol_type_head;
  struct species_opt_orient *mol_type_tail;
};

struct eff_dat_list
{
  struct eff_dat *eff_head;
  struct eff_dat *eff_tail;
};

struct element_list_head
{
  struct element_list *elml_head;
  struct element_list *elml_tail;
};

struct release_single_molecule_list
{
  struct release_single_molecule *rsm_head;
  struct release_single_molecule *rsm_tail;
  int rsm_count;
};

struct frame_data_list_head
{
  struct frame_data_list *frame_head;
  struct frame_data_list *frame_tail;
};

struct output_column_list
{
  struct output_column *column_head;
  struct output_column *column_tail;
};

struct output_set_list
{
  struct output_set *set_head;
  struct output_set *set_tail;
};

struct num_expr_list_head
{
  struct num_expr_list *value_head;
  struct num_expr_list *value_tail;
  int value_count;
  int shared;
};

struct output_times_inlist
{
  enum output_timer_type_t    type;
  double                      step;
  struct  num_expr_list_head  values;
};

struct reaction_arrow
{
  int                           flags;
  struct species_opt_orient     catalyst;
};

struct macro_subunit_assignment_list
{
  struct macro_subunit_assignment *assign_head;
  struct macro_subunit_assignment *assign_tail;
};

enum {
  RATE_UNSET    = -1,
  RATE_CONSTANT = 0,
  RATE_FILE     = 1,
  RATE_COMPLEX  = 2
};
struct reaction_rate
{
  int rate_type;
  union
  {
    double                  rate_constant;
    char                   *rate_file;
    struct complex_rate    *rate_complex;
  } v;
};

struct reaction_rates
{
  struct reaction_rate      forward_rate;
  struct reaction_rate      backward_rate;
};

struct diffusion_constant
{
  double    D;
  int       is_2d;
};

struct object_list
{
  struct object *obj_head;
  struct object *obj_tail;
};

struct output_times
{
  enum output_timer_type_t  timer_type;
  double  step_time;
  int     num_times;
  double *times;
};

struct element_connection_list_head
{
  struct element_connection_list *connection_head;
  struct element_connection_list *connection_tail;
  int connection_count;
};

struct vertex_list_head
{
  struct vertex_list *vertex_head;
  struct vertex_list *vertex_tail;
  int vertex_count;
};



/***************************************************************************
 * Parser state structure
 ***************************************************************************/

struct mdlparse_vars
{
  /* Line number where last top-level (i.e. non-nested) multi-line (C-style)
   * comment was started in the current MDL file. */
  int           comment_started;

  /* Line numbers and filenames for all of the currently parsing files */
  u_int         line_num[MAX_INCLUDE_DEPTH];
  char const   *include_filename[MAX_INCLUDE_DEPTH];

  /* Stack pointer for filename/line number stack */
  u_int         include_stack_ptr;

  /* The world we are constructing */
  struct volume *vol;

  /* Stack of object names used to build up qualified names as objects are
   * created */
  struct name_list *object_name_list;
  struct name_list *object_name_list_end;

  /* --------------------------------------------- */
  /* Pointers to objects being created or modified */
  
  /* Object currently being created or modified. */
  struct object *current_object;

  /* Pointer to surface class currently being created or modified */
  struct species *current_surface_class;

  /* Release site currently being created or modified */
  struct release_site_obj *current_release_site;

  /* Current polygon object being created or modified */
  struct polygon_object *current_polygon;

  /* Current region object being created or modified */
  struct region *current_region;

  /* Current macromolecular complex being created or modified */
  struct complex_species *current_complex;

  /* --------------------------------------------- */
  /* Intermediate state for counting */

  /* Tracks the count statement types seen in the current output block
   * (currently only COUNT and TRIGGER). */
  int count_flags;

  /* Custom header for reaction output */
  char *header_comment;

  /* Flag indicating whether to display the exact time */
  byte exact_time_flag;

  /* --------------------------------------------- */
  /* Intermediate state for macromolecules */

  /* Name of the current complex (used for error reporting */
  char *complex_name;

  /* Is the complex a surface complex or a volume complex? */
  int complex_type;

  /* Topology of the macromolecule -- required mostly for error checking */
  struct macro_topology *complex_topo;

  /* Relationships for this macromolecule */
  struct macro_relationship *complex_relations;

  /* --------------------------------------------- */
  /* Intermediate state for regions */
  int allow_patches;

  /* --------------------------------------------- */
  /* Temporary allocators */
  struct mem_helper *species_list_mem;
  struct mem_helper *mol_data_list_mem;
  struct mem_helper *output_times_mem;
  struct mem_helper *sym_list_mem;
  struct mem_helper *path_mem;
  struct mem_helper *prod_mem;
};


/***************************************************************************
 * Macromolecule related parse-time structures
 ***************************************************************************/
struct macro_topology_element
{
  struct macro_topology_element *next;
  char *name;
  int   dimensionality;
  int  *dimensions;
};

struct macro_topology
{
  struct macro_topology_element *head;
  int total_subunits;
};

struct macro_subunit_spec
{
  struct macro_subunit_spec *next;
  int from, to;
};


struct macro_subunit_assignment
{
  struct macro_subunit_assignment *next;
  struct macro_subunit_spec  *head;
  struct species             *what;
  short                       orient;
};

struct complex_species;

struct macro_geometry
{
  struct macro_geometry *next;
  struct num_expr_list  *index;
  struct vector3         location;
};

struct macro_relationship
{
  struct macro_relationship *next;
  char *name;
  struct num_expr_list *indices;
};

struct macro_rate_clause
{
  struct macro_rate_clause   *next;
  char                       *name;
  int                         invert;
  struct species             *species;
  short                       orient;
};

struct macro_rate_rule
{
  struct macro_rate_rule     *next;
  struct macro_rate_clause   *clauses;
  double                      rate;
};

struct macro_rate_ruleset
{
  struct macro_rate_ruleset  *next;
  char                       *name;
  struct macro_rate_rule     *rules;
};

/***************************************************************************
 * Declarations for functions defined in mdlparse.y
 ***************************************************************************/

void mdlerror(struct mdlparse_vars *mpvp, char const *str);
void mdlerror_fmt(struct mdlparse_vars *mpvp, char const *fmt, ...)
  PRINTF_FORMAT(2);
int mdlparse_init(struct volume *vol);

#endif
