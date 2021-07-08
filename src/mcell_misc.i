/******************************************************************************
 *
 * Copyright (C) 2006-2015 by
 * The Salk Institute for Biological Studies and
 * Pittsburgh Supercomputing Center, Carnegie Mellon University
 *
 * Use of this source code is governed by an MIT-style
 * license that can be found in the LICENSE file or at
 * https://opensource.org/licenses/MIT.
 *
******************************************************************************/

void mcell_print_version();

void mcell_print_usage(const char *executable_name);

void mcell_print_stats();

void mcell_dump_state(MCELL_STATE *state);

//int mcell_argparse(int argc, char **argv, MCELL_STATE *state); // swig-generated code cannot be built on macos

struct num_expr_list *mcell_copysort_numeric_list(struct num_expr_list *head);

void mcell_sort_numeric_list(struct num_expr_list *head);

void mcell_free_numeric_list(struct num_expr_list *nel);

MCELL_STATUS mcell_generate_range(struct num_expr_list_head *list, double start,
                                  double end, double step);

int mcell_generate_range_singleton(struct num_expr_list_head *lh, double value);

// XXX this is a temporary hack to be able to print in mcell.c
// since mcell disables regular printf
void mcell_print(const char *message);
