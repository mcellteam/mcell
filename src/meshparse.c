#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "sym_table.h"
#include "mcell_objects.h"
#include "mcell_structs.h"
#include "mcell_misc.h"
#include "meshparse.h"

#define EPS_C 1e-12

int generate_range(
    struct num_expr_list_head_dg *list,
    double start,
    double end,
    double step) {
  list->value_head = NULL;
  list->value_tail = NULL;
  list->value_count = 0;
  list->shared = 0;

  if (step > 0) {
    for (double tmp_dbl = start;
         tmp_dbl < end || !distinguishable(tmp_dbl, end, EPS_C) ||
             fabs(end - tmp_dbl) <= EPS_C;
         tmp_dbl += step) {
      if (advance_range_dg(list, tmp_dbl))
        return 1;
    }
  } else /* if (step < 0) */
  {
    for (double tmp_dbl = start;
         tmp_dbl > end || !distinguishable(tmp_dbl, end, EPS_C) ||
             fabs(end - tmp_dbl) <= EPS_C;
         tmp_dbl += step) {
      if (advance_range_dg(list, tmp_dbl))
        return 1;
    }
  }
  return 0;
}

// This is the same as advance_range in mcell_misc.h, but including that header
// here causes a number of build problems that are currently difficult to
// resolve.
int advance_range_dg(struct num_expr_list_head_dg *list, double tmp_dbl) {
  struct num_expr_list *nel;
  nel = (struct num_expr_list *)malloc(sizeof(struct num_expr_list));
  if (nel == NULL) {
    mcell_free_numeric_list(list->value_head);
    list->value_head = list->value_tail = NULL;
    return 1;
  }
  nel->value = tmp_dbl;
  nel->next = NULL;

  ++list->value_count;
  if (list->value_tail != NULL)
    list->value_tail->next = nel;
  else
    list->value_head = nel;
  list->value_tail = nel;
  return 0;
}

void mcell_free_numeric_list(struct num_expr_list *nel) {
  while (nel != NULL) {
    struct num_expr_list *n = nel;
    nel = nel->next;
    free(n);
  }
}

struct vector3 *point_scalar(double val) {
  struct vector3 *vec;
  vec = (struct vector3 *)malloc(sizeof(struct vector3));
  if (!vec)
    return NULL;

  vec->x = val;
  vec->y = val;
  vec->z = val;
  return vec;
}

void object_list_singleton(struct object_list *head, struct object *objp) {
  objp->next = NULL;
  head->obj_tail = head->obj_head = objp;
}

void add_object_to_list(struct object_list *head, struct object *objp) {
  objp->next = NULL;
  head->obj_tail = head->obj_tail->next = objp;
}

int init_top_level_objs(struct dyngeom_parse_vars *dg_parse) {
  if ((dg_parse->obj_sym_table = init_symtab(1024)) == NULL) {
    return 1;
  }

  struct sym_table *sym;
  if ((sym = store_sym(
      "WORLD_OBJ", OBJ, dg_parse->obj_sym_table, NULL)) == NULL) {
    return 1;
  }

  dg_parse->root_object = (struct object *)sym->value;
  dg_parse->root_object->object_type = META_OBJ;
  if (!(dg_parse->root_object->last_name = CHECKED_STRDUP_NODIE("", NULL))) {
    return 1;
  }

  if ((sym = store_sym(
      "WORLD_INSTANCE", OBJ, dg_parse->obj_sym_table, NULL)) == NULL) {
    return 1;
  }

  dg_parse->root_instance = (struct object *)sym->value;
  dg_parse->root_instance->object_type = META_OBJ;
  if (!(dg_parse->root_instance->last_name = CHECKED_STRDUP("", NULL))) {
    return 1;
  }

  dg_parse->current_object = dg_parse->root_instance;

  return 0;
}

struct sym_table *dg_start_object(
    struct dyngeom_parse_vars *dg_parse,
    char *name) {
  // Create new fully qualified name.
  char *new_name;
  struct object_creation obj_creation;
  obj_creation.object_name_list = dg_parse->object_name_list;
  obj_creation.object_name_list_end = dg_parse->object_name_list_end;
  if ((new_name = push_object_name(&obj_creation, name)) == NULL) {
    free(name);
    return NULL;
  }
  dg_parse->object_name_list = obj_creation.object_name_list;
  dg_parse->object_name_list_end = obj_creation.object_name_list_end;

  // Create the symbol, if it doesn't exist yet.
  struct object *obj_ptr = make_new_object(
      dg_parse->obj_sym_table, new_name, 0);
  if (obj_ptr == NULL) {
    if (name != new_name) {
      free(name);
    }
    free(new_name);
    return NULL;
  }

  struct sym_table *sym_ptr = obj_ptr->sym;
  obj_ptr->last_name = name;

  // Set parent object, make this object "current".
  obj_ptr->parent = dg_parse->current_object;
  dg_parse->current_object = obj_ptr;

  return sym_ptr;
}

void dg_finish_object(struct dyngeom_parse_vars *dg_parse) {
  struct object_creation obj_creation;
  obj_creation.object_name_list_end = dg_parse->object_name_list_end;

  pop_object_name(&obj_creation);
  dg_parse->object_name_list_end = obj_creation.object_name_list_end;
  dg_parse->current_object = dg_parse->current_object->parent;
}
