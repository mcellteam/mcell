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

#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "config.h"
#include "mcell_structs.h"

#include "argparse.h"
#include "version_info.h"
#include "logging.h"
#include "version_info.h"
#include "mcell_misc.h"


/* declaration of static functions */
static void swap_double(double *x, double *y);

/************************************************************************
 *
 * mcell_print_version prints the version string
 *
 ************************************************************************/
void mcell_print_version() { print_version(mcell_get_log_file()); }

/************************************************************************
 *
 * mcell_print_usage prints the usage information
 *
 ************************************************************************/
void mcell_print_usage(const char *executable_name) {
  print_usage(mcell_get_log_file(), executable_name);
}

/************************************************************************
 *
 * mcell_print_stats prints the simulation stats
 *
 ************************************************************************/
void mcell_print_stats() { mem_dump_stats(mcell_get_log_file()); }

/************************************************************************
 *
 * function for printing a string
 *
 * XXX: This is a temporary hack to be able to print in mcell.c
 *      since mcell disables regular printf
 *
 ************************************************************************/
void mcell_print(const char *message) { mcell_log("%s", message); }

/************************************************************************
 *
 * mcell_argparse parses the commandline and sets the
 * corresponding parts of the state (seed #, logging, ...)
 *
 ************************************************************************/
int mcell_argparse(int argc, char **argv, MCELL_STATE *state) {
  return argparse_init(argc, argv, state);
}

/*************************************************************************
 mcell_copysort_numeric_list:
    Copy and sort a num_expr_list in ascending numeric order.

 In:  parse_state:  parser state
      head:  the list to sort
 Out: list is sorted
*************************************************************************/
struct num_expr_list *mcell_copysort_numeric_list(struct num_expr_list *head) {
  struct num_expr_list_head new_head;
  if (mcell_generate_range_singleton(&new_head, head->value))
    return NULL;

  head = head->next;
  while (head != NULL) {
    struct num_expr_list *insert_pt, **prev;
    for (insert_pt = new_head.value_head, prev = &new_head.value_head;
         insert_pt != NULL;
         prev = &insert_pt->next, insert_pt = insert_pt->next) {
      if (insert_pt->value >= head->value)
        break;
    }

    struct num_expr_list *new_item =
        CHECKED_MALLOC_STRUCT(struct num_expr_list, "numeric array");
    if (new_item == NULL) {
      mcell_free_numeric_list(new_head.value_head);
      return NULL;
    }

    new_item->next = insert_pt;
    new_item->value = head->value;
    *prev = new_item;
    if (insert_pt == NULL)
      new_head.value_tail = new_item;
    head = head->next;
  }

  return new_head.value_head;
}

/*************************************************************************
 mcell_sort_numeric_list:
    Sort a num_expr_list in ascending numeric order.  N.B. This uses bubble
    sort, which is O(n^2).  Don't use it if you expect your list to be very
    long.  The list is sorted in-place.

 In:  head:  the list to sort
 Out: list is sorted
*************************************************************************/
void mcell_sort_numeric_list(struct num_expr_list *head) {
  struct num_expr_list *curr, *next;
  int done = 0;
  while (!done) {
    done = 1;
    curr = head;
    while (curr != NULL) {
      next = curr->next;
      if (next != NULL) {
        if (curr->value > next->value) {
          done = 0;
          swap_double(&curr->value, &next->value);
        }
      }
      curr = next;
    }
  }
}

/*************************************************************************
 mcell_free_numeric_list:
    Free a num_expr_list.

 In:  nel:  the list to free
 Out: all elements are freed
*************************************************************************/
void mcell_free_numeric_list(struct num_expr_list *nel) {
  free_numeric_list(nel);
}

/*************************************************************************
 mcell_generate_range:
    Generate a num_expr_list containing the numeric values from start to end,
    incrementing by step.

 In:  state: the simulation state
      list:  destination to receive list of values
      start: start of range
      end:   end of range
      step:  increment
 Out: 0 on success, 1 on failure.  On success, list is filled in.
*************************************************************************/
MCELL_STATUS mcell_generate_range(struct num_expr_list_head *list, double start,
                                  double end, double step) {
  if (generate_range(list, start, end, step)) {
    return MCELL_FAIL; 
  } else {
    return MCELL_SUCCESS; 
  }
}

/*************************************************************************
 mcell_generate_range_singleton:
    Generate a numeric list containing a single value.

 In:  lh:   list to receive value
      value: value for list
 Out: 0 on success, 1 on failure
*************************************************************************/
int mcell_generate_range_singleton(struct num_expr_list_head *lh,
                                   double value) {

  struct num_expr_list *nel =
      CHECKED_MALLOC_STRUCT(struct num_expr_list, "numeric array");
  if (nel == NULL)
    return 1;
  lh->value_head = lh->value_tail = nel;
  lh->value_count = 1;
  lh->shared = 0;
  lh->value_head->value = value;
  lh->value_head->next = NULL;
  return 0;
}

/************************************************************************
 swap_double:
 In:  x, y: doubles to swap
 Out: Swaps references to two double values.
 ***********************************************************************/
void swap_double(double *x, double *y) {
  double temp;

  temp = *x;
  *x = *y;
  *y = temp;
}

/************************************************************************
 mcell_find_include_file:
      Find the path for an include file.  For an absolute include path, the
      file path is unmodified, but for a relative path, the resultant path will
      be relative to the currently parsed file.

 In:  path: path from include statement
      cur_path: path of current include file
 Out: allocated buffer containing path of the include file, or NULL if the file
      path couldn't be allocated.  If we ever use a more complex mechanism for
      locating include files, we may also return NULL if no file could be
      located.
 ***********************************************************************/
char *mcell_find_include_file(char const *path, char const *cur_path) {
  char *candidate = NULL;
  if (path[0] == '/')
    candidate = strdup(path);
  else {
    char *last_slash = strrchr(cur_path, '/');
#ifdef _WIN32
    char *last_bslash = strrchr(cur_path, '\\');
    if (last_bslash > last_slash)
      last_slash = last_bslash;
#endif
    if (last_slash == NULL)
      candidate = strdup(path);
    else
      candidate = CHECKED_SPRINTF("%.*s/%s", (int)(last_slash - cur_path),
                                  cur_path, path);
  }

  return candidate;
}
