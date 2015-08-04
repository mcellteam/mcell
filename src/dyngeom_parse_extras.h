#ifndef DYNGEOM_PARSE_EXTRAS_H
#define DYNGEOM_PARSE_EXTRAS_H

int parse_dg();

struct num_expr_list_head_dg {
  struct num_expr_list *value_head;
  struct num_expr_list *value_tail;
  int value_count;
  int shared;
};

struct dyngeom_parse_vars {
  struct sym_table_head *obj_sym_table;
  struct object *root_object;
  struct object *root_instance;
  struct object *current_object;
  struct name_list *object_name_list;
  struct name_list *object_name_list_end;
};

struct object_list {
  struct object *obj_head;
  struct object *obj_tail;
};


int init_top_level_objs(struct dyngeom_parse_vars *dg_parse);
void object_list_singleton(struct object_list *head, struct object *objp);
void add_object_to_list(struct object_list *head, struct object *objp);
struct vector3 *point_scalar(double val);
int advance_range_dg(struct num_expr_list_head_dg *list, double tmp_dbl);
void mcell_free_numeric_list_dg(struct num_expr_list *nel);
int generate_range(
    struct num_expr_list_head_dg *list,
    double start,
    double end,
    double step);
struct sym_table *dg_start_object(
    struct dyngeom_parse_vars *dg_parse, char *name);
void dg_finish_object(struct dyngeom_parse_vars *dg_parse);

struct dyngeom_parse_vars *dg_parse;

#endif
