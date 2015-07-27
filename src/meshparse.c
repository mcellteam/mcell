#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "meshparse.h"

#define EPS_C 1e-12

int distinguishable(double a, double b, double eps) {
  double c = fabs(a - b);
  a = fabs(a);
  if (a < 1) {
    a = 1;
  }
  b = fabs(b);

  if (b < a) {
    eps *= a;
  } else {
    eps *= b;
  }
  return (c > eps);
}

int generate_range(struct num_expr_list_head *list, double start, double end, double step) {
  list->value_head = NULL;
  list->value_tail = NULL;
  list->value_count = 0;
  list->shared = 0;

  if (step > 0) {
    for (double tmp_dbl = start;
         tmp_dbl < end || !distinguishable(tmp_dbl, end, EPS_C) ||
             fabs(end - tmp_dbl) <= EPS_C;
         tmp_dbl += step) {
      if (advance_range(list, tmp_dbl))
        return 1;
    }
  } else /* if (step < 0) */
  {
    for (double tmp_dbl = start;
         tmp_dbl > end || !distinguishable(tmp_dbl, end, EPS_C) ||
             fabs(end - tmp_dbl) <= EPS_C;
         tmp_dbl += step) {
      if (advance_range(list, tmp_dbl))
        return 1;
    }
  }
  return 0;
}

int advance_range(struct num_expr_list_head *list, double tmp_dbl) {
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
