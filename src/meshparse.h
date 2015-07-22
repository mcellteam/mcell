#ifndef MESHPARSE
#define MESHPARSE

#include <stdio.h>

struct num_expr_list {
  struct num_expr_list *next;
  double value; /* Value of one element of the expression */
};

struct num_expr_list_head {
  struct num_expr_list *value_head;
  struct num_expr_list *value_tail;
  int value_count;
  int shared;
};

enum partition_axis_t {
  X_PARTS, /* X-axis partitions */
  Y_PARTS, /* Y-axis partitions */
  Z_PARTS  /* Z-axis partitions */
};

struct sym_table {
  struct sym_table *next; /* Chain to next symbol in this bin of the hash */
  int sym_type;           /* Symbol Type */
  char *name;             /* Name of symbol*/
  void *value;            /* Stored value, cast by sym_type */
};

void mcell_free_numeric_list(struct num_expr_list *nel);
int advance_range(struct num_expr_list_head *list, double tmp_dbl);
int generate_range(struct num_expr_list_head *list, double start, double end, double step);

#endif
