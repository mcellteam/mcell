#ifndef MESHPARSE
#define MESHPARSE

#include <stdio.h>

#include "mcell_structs.h"

// This is the same as num_expr_list_head in mcell_init.h, but including that
// header here causes a number of build problems that are currently difficult
// to resolve.
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

int init_top_level_objs(struct dyngeom_parse_vars *dg_parse);
void object_list_singleton(struct object_list *head, struct object *objp);
void add_object_to_list(struct object_list *head, struct object *objp);
struct vector3 *point_scalar(double val);
void mcell_free_numeric_list(struct num_expr_list *nel);
int advance_range_dg(struct num_expr_list_head_dg *list, double tmp_dbl);
int generate_range(
    struct num_expr_list_head_dg *list,
    double start,
    double end,
    double step);
struct sym_table *dg_start_object(
    struct dyngeom_parse_vars *dg_parse, char *name);
void dg_finish_object(struct dyngeom_parse_vars *dg_parse);

#endif
