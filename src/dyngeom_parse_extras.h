#ifndef DYNGEOM_PARSE_EXTRAS_H
#define DYNGEOM_PARSE_EXTRAS_H

#include "mcell_structs.h"

#define MAX_INCLUDE_DEPTH 16

struct dyngeom_parse_vars * create_dg_parse(struct volume *state);

struct dyngeom_parse_vars {
  struct sym_table_head *reg_sym_table;
  struct sym_table_head *obj_sym_table;
  struct object *root_object;
  struct object *root_instance;
  struct object *current_object;
  struct region *current_region;
  struct name_list *object_name_list;
  struct name_list *object_name_list_end;
  char const *curr_file; /* Name of MDL file currently being parsed */

  /* Line numbers and filenames for all of the currently parsing files */
  u_int line_num[MAX_INCLUDE_DEPTH];
  char const *include_filename[MAX_INCLUDE_DEPTH];

  /* Stack pointer for filename/line number stack */
  u_int include_stack_ptr;

  /* Line number where last top-level (i.e. non-nested) multi-line (C-style)
   * comment was started in the current MDL file. */
  int comment_started;
};

#include "mcell_objects.h"

int parse_dg(struct dyngeom_parse_vars *dg_parse, char *dynamic_geometry_filename);
int parse_dg_init(struct dyngeom_parse_vars *dg_parse, char *dynamic_geometry_filename, struct volume *state);

void setup_root_obj_inst(struct dyngeom_parse_vars *dg_parse_vars, struct volume *state);
struct sym_entry *dg_start_object(
    struct dyngeom_parse_vars *dg_parse_vars,
    char *name);
struct object *dg_start_object_simple(struct dyngeom_parse_vars *dg_parse_vars,
                                      struct object_creation *obj_creation,
                                      char *name);
struct object *dg_new_polygon_list(
    struct dyngeom_parse_vars *dg_parse_vars,
    char *obj_name);
void dg_finish_object(struct dyngeom_parse_vars *dg_parse_vars);
struct region *dg_create_region(
    struct sym_table_head *reg_sym_table,
    struct object *objp,
    char *name);
struct region *dg_make_new_region(
    struct sym_table_head *reg_sym_table,
    char *obj_name,
    char *region_last_name);
int dg_copy_object_regions(
    struct dyngeom_parse_vars *dg_parse_vars,
    struct object *dst_obj,
    struct object *src_obj);
int dg_deep_copy_object(
    struct dyngeom_parse_vars *dg_parse_vars,
    struct object *dst_obj,
    struct object *src_obj);
struct sym_entry *dg_existing_object(
    struct dyngeom_parse_vars *dg_parse_vars,
    char *name);
char *find_include_file(char const *path, char const *cur_path);

#endif
