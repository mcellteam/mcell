#ifndef DYNGEOM_PARSE_H
#define DYNGEOM_PARSE_H

struct num_expr_list_head_dg {
  struct num_expr_list *value_head;
  struct num_expr_list *value_tail;
  int value_count;
  int shared;
};

struct object_list {
  struct object *obj_head;
  struct object *obj_tail;
};

struct dyngeom_parse_vars {
  struct sym_table_head *obj_sym_table;
  struct object *root_object;
  struct object *root_instance;
  struct object *current_object;
  struct name_list *object_name_list;
  struct name_list *object_name_list_end;
};

struct dyngeom_parse_vars *dg_parse;

int init_top_level_objs(struct dyngeom_parse_vars *dg_parse);
struct sym_table *dg_start_object(
    struct dyngeom_parse_vars *dg_parse, char *name);
void dg_finish_object(struct dyngeom_parse_vars *dg_parse);
int parse_dg(char *dynamic_geometry_filename);

#endif
