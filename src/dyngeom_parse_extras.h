#ifndef DYNGEOM_PARSE_EXTRAS_H
#define DYNGEOM_PARSE_EXTRAS_H

#include "mcell_structs.h"

int create_dg_parse();
int parse_dg(char *dynamic_geometry_filename);

struct dyngeom_parse_vars {
  struct sym_table_head *reg_sym_table;
  struct sym_table_head *obj_sym_table;
  struct object *root_object;
  struct object *root_instance;
  struct object *current_object;
  struct name_list *object_name_list;
  struct name_list *object_name_list_end;
};

int init_top_level_objs(struct dyngeom_parse_vars *dg_parse_vars);
int setup_root_obj_inst(struct dyngeom_parse_vars *dg_parse_vars);
struct sym_table *dg_start_object(
    struct dyngeom_parse_vars *dg_parse_vars,
    char *name);
void dg_finish_object(struct dyngeom_parse_vars *dg_parse_vars);
struct region *dg_create_region(
    struct sym_table_head *reg_sym_table,
    struct object *objp,
    char *name);
struct region *dg_make_new_region(
    struct sym_table_head *reg_sym_table,
    char *obj_name,
    char *region_last_name);

struct dyngeom_parse_vars *dg_parse;

#endif
