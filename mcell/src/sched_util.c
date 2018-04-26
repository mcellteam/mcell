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

#include "config.h"

#include <float.h>
#include <string.h>
#include <stdlib.h>

#include "sched_util.h"

/*************************************************************************
ae_list_sort:
  In: head of a linked list of abstract_elements
  Out: head of the newly sorted list
  Note: uses mergesort
*************************************************************************/

struct abstract_element *ae_list_sort(struct abstract_element *ae) {
  struct abstract_element *stack[64];
  int stack_n[64];
  struct abstract_element *left = NULL, *right = NULL, *merge = NULL,
                          *tail = NULL;
  int si = 0;

  if (ae == NULL)
    return NULL;

  while (ae != NULL) {
    if (ae->next == NULL) {
      stack[si] = ae;
      stack_n[si] = 1;
      ae = NULL;
      si++;
    } else if (ae->t <= ae->next->t) {
      stack[si] = ae;
      stack_n[si] = 2;
      ae = ae->next->next;
      stack[si]->next->next = NULL;
      si++;
    } else {
      stack[si] = ae->next;
      stack_n[si] = 2;
      left = ae;
      ae = ae->next->next;
      stack[si]->next = left;
      left->next = NULL;
      si++;
    }
    while (si > 1 && stack_n[si - 1] * 2 >= stack_n[si - 2]) {
      stack_n[si - 2] += stack_n[si - 1];

      left = stack[si - 2];
      right = stack[si - 1];
      if (left->t <= right->t) {
        merge = left;
        left = left->next;
      } else {
        merge = right;
        right = right->next;
      }
      merge->next = NULL;
      tail = merge;

      while (1) {
        if (left == NULL) {
          tail->next = right;
          break;
        }
        if (right == NULL) {
          tail->next = left;
          break;
        }

        if (left->t <= right->t) {
          tail->next = left;
          tail = left;
          left = left->next;
        } else {
          tail->next = right;
          tail = right;
          right = right->next;
        }
      }

      stack[si - 2] = merge;
      si--;
    }
  }

  while (si > 1) /* Exact duplicate of code in loop--keep it this way! */
  {
    stack_n[si - 2] += stack_n[si - 1];

    left = stack[si - 2];
    right = stack[si - 1];
    if (left->t <= right->t) {
      merge = left;
      left = left->next;
    } else {
      merge = right;
      right = right->next;
    }
    merge->next = NULL;
    tail = merge;

    while (1) {
      if (left == NULL) {
        tail->next = right;
        break;
      }
      if (right == NULL) {
        tail->next = left;
        break;
      }

      if (left->t <= right->t) {
        tail->next = left;
        tail = left;
        left = left->next;
      } else {
        tail->next = right;
        tail = right;
        right = right->next;
      }
    }

    stack[si - 2] = merge;
    si--;
  }

  return stack[0];
}

/*************************************************************************
create_scheduler:
  In: timestep per slot in this scheduler
      time for all slots in this scheduler
      maximum number of slots in this scheduler
      the current time
  Out: pointer to a new instance of schedule_helper; pass this to later
       functions.  (Dispose of with delete_scheduler.)  Returns NULL
       if out of memory.
*************************************************************************/

struct schedule_helper *create_scheduler(double dt_min, double dt_max,
                                         int maxlen, double start_iterations) {
  double n_slots = dt_max / dt_min;
  int len;

  if (n_slots < (double)(maxlen - 1))
    len = (int)n_slots + 1;
  else
    len = maxlen;

  if (len < 2)
    len = 2;

  struct schedule_helper *sh = NULL;
  sh = (struct schedule_helper *)malloc(sizeof(struct schedule_helper));
  if (sh == NULL)
    return NULL;
  memset(sh, 0, sizeof(struct schedule_helper));

  sh->dt = dt_min;
  sh->dt_1 = 1 / dt_min;

  sh->now = start_iterations;
  sh->buf_len = len;

  sh->circ_buf_count = (int *)calloc(len, sizeof(int));
  if (sh->circ_buf_count == NULL)
    goto failure;

  sh->circ_buf_head = (struct abstract_element **)calloc(
      len * 2, sizeof(struct abstract_element*));
  if (sh->circ_buf_head == NULL)
    goto failure;
  sh->circ_buf_tail = sh->circ_buf_head + len;

  if (sh->dt * sh->buf_len < dt_max) {
    sh->next_scale =
        create_scheduler(dt_min * len, dt_max, maxlen, sh->now + dt_min * len);
    if (sh->next_scale == NULL)
      goto failure;
    sh->next_scale->depth = sh->depth + 1;
  }

  return sh;

failure:
  if (sh != NULL)
    delete_scheduler(sh);
  return NULL;
}

/*************************************************************************
schedule_insert:
  In: scheduler that we are using
      data to schedule (assumed to start with abstract_element struct)
      flag to indicate whether times in the "past" go into the list
         of current events (if 0, go into next event, not current).
  Out: 0 on success, 1 on memory allocation failure.  Data item is
       placed in scheduler at end of list for its time slot.
*************************************************************************/

int schedule_insert(struct schedule_helper *sh, void *data,
                    int put_neg_in_current) {
  struct abstract_element *ae = (struct abstract_element *)data;

  if (put_neg_in_current && ae->t < sh->now) {
    /* insert item into current list */

    sh->current_count++;
    if (sh->current_tail == NULL) {
      sh->current = sh->current_tail = ae;
      ae->next = NULL;
    } else {
      sh->current_tail->next = ae;
      sh->current_tail = ae;
      ae->next = NULL;
    }
    return 0;
  }

  /* insert item into future lists */
  sh->count++;
  double nsteps = (ae->t - sh->now) * sh->dt_1;

  if (nsteps < ((double)sh->buf_len)) {
    /* item fits in array for this scale */

    int i;
    if (nsteps < 0.0)
      i = sh->index;
    else
      i = (int)nsteps + sh->index;
    if (i >= sh->buf_len)
      i -= sh->buf_len;

    if (sh->circ_buf_tail[i] == NULL) {
      sh->circ_buf_count[i] = 1;
      sh->circ_buf_head[i] = sh->circ_buf_tail[i] = ae;
      ae->next = NULL;
    } else {
      sh->circ_buf_count[i]++;

      /* For schedulers other than the first tier, maintain a LIFO ordering */
      if (sh->depth) {
        ae->next = sh->circ_buf_head[i];
        sh->circ_buf_head[i] = ae;
      }

      /* For first-tier scheduler, maintain FIFO ordering */
      else {
        sh->circ_buf_tail[i]->next = ae;
        ae->next = NULL;
        sh->circ_buf_tail[i] = ae;
      }
    }
  } else {
    /* item fits in array for coarser scale */

    if (sh->next_scale == NULL) {
      sh->next_scale = create_scheduler(
          sh->dt * sh->buf_len, sh->dt * sh->buf_len * sh->buf_len, sh->buf_len,
          sh->now + sh->dt * (sh->buf_len - sh->index));
      if (sh->next_scale == NULL)
        return 1;
      sh->next_scale->depth = sh->depth + 1;
    }

    /* insert item at coarser scale and insist that item is not placed in
     * "current" list */
    return schedule_insert(sh->next_scale, data, 0);
  }

  return 0;
}

/*************************************************************************
unlink_list_item:
  Removes a specific item from the linked list.

  In: struct abstract_element **hd - pointer to head of list
      struct abstract_element **tl - pointer to tail of list
      struct abstract_element *it - pointer to item to unlink
  Out: 0 on success, 1 if the item was not found
*************************************************************************/
static int unlink_list_item(struct abstract_element **hd,
                            struct abstract_element **tl,
                            struct abstract_element *it) {
  struct abstract_element *prev = NULL;

  while (*hd != NULL) {
    if (*hd == it) {
      (*hd) = (*hd)->next;
      if (*tl == it)
        *tl = prev;
      return 0;
    }

    prev = *hd;
    hd = &prev->next;
  }

  return 1;
}

/*************************************************************************
schedule_deschedule:
  Removes an item from the schedule.

  In: struct schedule_helper *sh - the scheduler from which to remove
      void  *data - the item to remove
  Out: 0 on success, 1 if the item was not found
*************************************************************************/
int schedule_deschedule(struct schedule_helper *sh, void *data) {
  struct abstract_element *ae = (struct abstract_element *)data;

  /* If the item is in "current" */
  if (sh->current && ae->t < sh->now) {
    if (unlink_list_item(&sh->current, &sh->current_tail, ae))
      return 1;

    --sh->current_count;
    return 0;
  }

  double nsteps = (ae->t - sh->now) * sh->dt_1;
  if (nsteps < ((double)sh->buf_len)) {
    int list_idx;
    if (nsteps < 0.0)
      list_idx = sh->index;
    else
      list_idx = (int)nsteps + sh->index;
    if (list_idx >= sh->buf_len)
      list_idx -= sh->buf_len;

    if (unlink_list_item(&sh->circ_buf_head[list_idx],
                         &sh->circ_buf_tail[list_idx], ae)) {
      /* If we fail to find it in this level, it may be in the next level.
       * Note that when we are descheduling, we may need to look in more than
       * one place, depending upon how long ago the item to be descheduled was
       * originally scheduled. */
      if (sh->next_scale) {
        if (!schedule_deschedule(sh->next_scale, data)) {
          --sh->count;
          return 0;
        } else
          return 1;
      } else
        return 1;
    }

    --sh->count;
    --sh->circ_buf_count[list_idx];
    return 0;
  } else {
    if (!sh->next_scale)
      return 1;

    if (!schedule_deschedule(sh->next_scale, data)) {
      --sh->count;
      return 0;
    } else
      return 1;
  }
}

/*************************************************************************
schedule_reschedule:
  Moves an item from one time to another in the schedule.

  In: struct schedule_helper *sh - the scheduler from which to remove
      void  *data - the item to remove
      double new_t - the new time for the item
  Out: 0 on success, 1 if the item was not found
*************************************************************************/
int schedule_reschedule(struct schedule_helper *sh, void *data, double new_t) {
  if (!schedule_deschedule(sh, data)) {
    struct abstract_element *ae = (struct abstract_element *)data;
    ae->t = new_t;
    return schedule_insert(sh, data, 1);
  } else
    return 1;
}

/*************************************************************************
schedule_advance:
  In: scheduler that we are using
      a pointer to the head-pointer for the list of the next time block
      a pointer to the tail-pointer for the list of the next time block
  Out: Number of items in the next block of time.  These items start
      with *head, and end with *tail.  Returns -1 on memory error.
*************************************************************************/

int schedule_advance(struct schedule_helper *sh, struct abstract_element **head,
                     struct abstract_element **tail) {
  int n;
  struct abstract_element *p, *nextp;

  if (head != NULL)
    *head = sh->circ_buf_head[sh->index];
  if (tail != NULL)
    *tail = sh->circ_buf_tail[sh->index];

  sh->circ_buf_head[sh->index] = sh->circ_buf_tail[sh->index] = NULL;
  sh->count -= n = sh->circ_buf_count[sh->index];
  sh->circ_buf_count[sh->index] = 0;

  sh->index++;
  sh->now += sh->dt;

  if (sh->index >= sh->buf_len) {
    /* Move events from coarser time scale to this time scale */

    sh->index = 0;
    if (sh->next_scale != NULL) {
      /* Save our depth */
      int old_depth = sh->depth;
      int conservecount = sh->count;

      /* Hack: Toggle the non-zero-ness of our depth to toggle FIFO/LIFO
       * behavior
       */
      sh->depth = old_depth ? 0 : -1;

      if (schedule_advance(sh->next_scale, &p, NULL) == -1) {
        sh->depth = old_depth;
        return -1;
      }
      while (p != NULL) {
        nextp = p->next;
        if (schedule_insert(sh, (void *)p, 0)) {
          sh->depth = old_depth;
          return -1;
        }
        p = nextp;
      }

      /* moved items were already counted when originally scheduled so don't
       * count again */
      sh->count = conservecount;

      /* restore our depth */
      sh->depth = old_depth;
    }
  }

  return n;
}

/*************************************************************************
schedule_next:
  In: scheduler that we are using
  Out: Next item to deal with.  If we are out of items for the current
       timestep, NULL is returned and the time is advanced to the next
       timestep.  If there is a memory error, NULL is returned and
       sh->error is set to 1.
*************************************************************************/

void *schedule_next(struct schedule_helper *sh) {
  void *data;

  if (sh->current == NULL) {
    sh->current_count = schedule_advance(sh, &sh->current, &sh->current_tail);
    if (sh->current_count == -1)
      sh->error = 1;
    return NULL;
  } else {
    sh->current_count--;
    data = sh->current;
    sh->current = sh->current->next;
    if (sh->current == NULL)
      sh->current_tail = NULL;
    return data;
  }
}


/*************************************************************************
schedule_peak:
  In: scheduler that we are using
  Out: Next item (e.g. molecule) to deal with. NULL if scheduler is empty

  This is very similar to schedule_next, but the idea here is that we don't
  want to change the state of anything (specifically current_count).
  XXX: The caller needs to reset sh->current when it's done "peaking"
*************************************************************************/
void *schedule_peak(struct schedule_helper *sh) {
  void *data;
  // Nothing to see here; move on
  if (sh->current == NULL) {
    return NULL; 
  }
  // Ooh. Found something!
  else {
    data = sh->current;
    sh->current = sh->current->next;
    return data;
  }
}


/*************************************************************************
schedule_anticipate:
  In: scheduler that we are using
      pointer to double to store the anticipated time of next event
  Out: 1 if there is an event anticipated, 0 otherwise
*************************************************************************/

int schedule_anticipate(struct schedule_helper *sh, double *t) {
  int i, j;
  double earliest_t = DBL_MAX;

  if (sh->current != NULL) {
    *t = sh->now;
    return 1;
  } else if (sh->count == 0)
    return 0;

  while (sh->next_scale != NULL && sh->count == sh->next_scale->count)
    sh = sh->next_scale;

  for (; sh != NULL; sh = sh->next_scale) {
    if (earliest_t < sh->now)
      break;

    for (i = 0; i < sh->buf_len; i++) {
      j = i + sh->index;
      if (j >= sh->buf_len)
        j -= sh->buf_len;
      if (sh->circ_buf_count[j] > 0) {
        earliest_t = sh->now + sh->dt * i;
        break;
      }
    }
  }

  if (earliest_t < DBL_MAX) {
    *t = earliest_t;
    return 1;
  } else
    return 0;
}

/*************************************************************************
schedule_cleanup:
  In: scheduler that we are using
      pointer to a function that will return 0 if an abstract_element is
        okay, or 1 if it is defunct
  Out: all defunct items are removed from the scheduler and returned as
       a linked list (so appropriate action can be taken, such as
       deallocation)
*************************************************************************/

struct abstract_element *
schedule_cleanup(struct schedule_helper *sh,
                 int (*is_defunct)(struct abstract_element*)) {
  struct abstract_element *defunct_list;
  struct abstract_element *ae;
  struct abstract_element *temp;
  struct schedule_helper *top;
  struct schedule_helper *shp;
  int i;

  defunct_list = NULL;

  top = sh;
  for (; sh != NULL; sh = sh->next_scale) {
    sh->defunct_count = 0;

    for (i = 0; i < sh->buf_len; i++) {
      /* Remove defunct elements from beginning of list */
      while (sh->circ_buf_head[i] != NULL &&
             (*is_defunct)(sh->circ_buf_head[i])) {
        temp = sh->circ_buf_head[i]->next;
        sh->circ_buf_head[i]->next = defunct_list;
        defunct_list = sh->circ_buf_head[i];
        sh->circ_buf_head[i] = temp;
        sh->circ_buf_count[i]--;
        sh->count--;
        for (shp = top; shp != sh; shp = shp->next_scale)
          shp->count--;
      }

      if (sh->circ_buf_head[i] == NULL) {
        sh->circ_buf_tail[i] = NULL;
      } else {
        /* Now remove defunct elements from later in list */
        for (ae = sh->circ_buf_head[i]; ae != NULL; ae = ae->next) {
          while (ae->next != NULL && (*is_defunct)(ae->next)) {
            temp = ae->next->next;
            ae->next->next = defunct_list;
            defunct_list = ae->next;
            ae->next = temp;
            sh->circ_buf_count[i]--;
            sh->count--;
            for (shp = top; shp != sh; shp = shp->next_scale)
              shp->count--;
          }
          if (ae->next == NULL) {
            sh->circ_buf_tail[i] = ae;
            break;
          }
        }
      }
    }
  }

  return defunct_list;
}

/*************************************************************************
delete_scheduler:
  In: scheduler that we are using
  Out: No return value.  The scheduler is freed from dynamic memory.
*************************************************************************/

void delete_scheduler(struct schedule_helper *sh) {
  if (sh) {
    if (sh->next_scale != NULL)
      delete_scheduler(sh->next_scale);
    if (sh->circ_buf_head)
      free(sh->circ_buf_head);
    if (sh->circ_buf_count)
      free(sh->circ_buf_count);
    free(sh);
  }
}
