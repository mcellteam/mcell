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

#pragma once

/* Everything managed by scheduler must begin as if it were derived from
 * abstract_element */
struct abstract_element {
  struct abstract_element *next;
  double t; /* Time at which the element is scheduled */
};

/* Implements a multi-scale, discretized event scheduler */
struct schedule_helper {
  struct schedule_helper *next_scale; /* Next coarser time scale */

  double dt;   /* Timestep per slot */
  double dt_1; /* dt_1 = 1/dt */
  double now;  /* Start time of the scheduler */

  /* Items scheduled now or after now */
  int count;           /* Total number of items scheduled now or after */
  int buf_len;         /* Number of slots in the scheduler */
  int index;           /* Index of the next time block */
  int *circ_buf_count; /* How many items are scheduled in each slot */
  // Array of linked lists of scheduled items for each slot
  struct abstract_element **circ_buf_head; 
  // Array of tails of the linked lists
  struct abstract_element **circ_buf_tail; 

  /* Items scheduled before now */
  /* These events must be serviced before simulation can advance to now */
  int current_count;                     /* Number of current items */
  struct abstract_element *current;      /* List of items scheduled now */
  struct abstract_element *current_tail; /* Tail of list of items */

  int defunct_count; /* Number of defunct items (set by user)*/
  int error;         /* Error code (1 - on error, 0 - no errors) */
  int depth;         /* "Tier" of scheduler in timescale hierarchy, 0-based */
};

struct abstract_element *ae_list_sort(struct abstract_element *ae);

struct schedule_helper *create_scheduler(double dt_min, double dt_max,
                                         int maxlen, double start_iterations);

int schedule_insert(struct schedule_helper *sh, void *data,
                    int put_neg_in_current);
int schedule_deschedule(struct schedule_helper *sh, void *data);
int schedule_reschedule(struct schedule_helper *sh, void *data, double new_t);
/*void schedule_excert(struct schedule_helper *sh,void *data,void *blank,int
 * size);*/
int schedule_advance(struct schedule_helper *sh, struct abstract_element **head,
                     struct abstract_element **tail);

void *schedule_next(struct schedule_helper *sh);
void *schedule_peak(struct schedule_helper *sh);
#define schedule_add(x, y) schedule_insert((x), (y), 1)

int schedule_anticipate(struct schedule_helper *sh, double *t);
struct abstract_element *
schedule_cleanup(struct schedule_helper *sh,
                 int (*is_defunct)(struct abstract_element *e));

void delete_scheduler(struct schedule_helper *sh);
