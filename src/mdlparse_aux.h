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

#include "mcell_structs.h"
#include "mcell_react_out.h"
#include "mcell_reactions.h"

#define WILDCARD_PRESENT 0x1
#define TRIGGER_PRESENT 0x2
#define COUNT_PRESENT 0x4
#define EXPRESSION_PRESENT 0x8

/***************************************************************************
 * Miscellaneous parse-time structures
 ***************************************************************************/

struct arg {
  struct arg *next;
  byte arg_type; /* DBL, STR */
  void *arg_value;
};

struct arg_list {
  struct arg *arg_head;
  struct arg *arg_tail;
};

struct sm_dat_list {
  struct sm_dat *sm_head;
  struct sm_dat *sm_tail;
};

struct element_list_head {
  struct element_list *elml_head;
  struct element_list *elml_tail;
};

#if 0
struct output_set_list {
  struct output_set *set_head;
  struct output_set *set_tail;
};

struct num_expr_list_head {
  struct num_expr_list *value_head;
  struct num_expr_list *value_tail;
  int value_count;
  int shared;
};

struct output_times_inlist {
  enum output_timer_type_t type;
  double step;
  struct num_expr_list_head values;
};
#endif

struct diffusion_constant {
  double D;
  int is_2d;
};

struct output_times {
  enum output_timer_type_t timer_type;
  double step_time;
  int num_times;
  double *times;
};

struct element_connection_list_head {
  struct element_connection_list *connection_head;
  struct element_connection_list *connection_tail;
  int connection_count;
};

struct vertex_list_head {
  struct vertex_list *vertex_head;
  struct vertex_list *vertex_tail;
  int vertex_count;
};

struct parse_mcell_species_list_item {
  struct parse_mcell_species_list_item *next;
  struct mcell_species_spec *spec;
};

struct parse_mcell_species_list {
  struct parse_mcell_species_list_item *species_head;
  struct parse_mcell_species_list_item *species_tail;
  int species_count;
};

struct species_list_item {
  struct species_list_item *next;
  struct species *spec;
};

struct species_list {
  struct species_list_item *species_head;
  struct species_list_item *species_tail;
  int species_count;
};

/***************************************************************************
 * Parser state structure
 ***************************************************************************/

struct mdlparse_vars {
  /* Line number where last top-level (i.e. non-nested) multi-line (C-style)
   * comment was started in the current MDL file. */
  int comment_started;

  /* Line numbers and filenames for all of the currently parsing files */
  u_int line_num[MAX_INCLUDE_DEPTH];
  char const *include_filename[MAX_INCLUDE_DEPTH];

  /* Stack pointer for filename/line number stack */
  u_int include_stack_ptr;

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
 * Declarations for functions defined in mdlparse.y
 ***************************************************************************/

void mdlerror(struct mdlparse_vars *parse_state, char const *str);
void mdlerror_fmt(struct mdlparse_vars *parse_state, char const *fmt, ...)
    PRINTF_FORMAT(2);
int mdlparse_init(struct volume *vol);
int mdlparse_file(struct mdlparse_vars *parse_state, char const *name);
