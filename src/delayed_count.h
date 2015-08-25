/******************************************************************************
 *
 * Copyright (C) 2006-2015 by
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

#ifndef DELAYED_COUNT_H
#define DELAYED_COUNT_H

#include "config.h"

#include "util.h"

#define UPDATE_COUNT(state, ptr, amt) do {                                         \
  if (state->sequential) (ptr) += (amt);                                    \
  else {                                                                    \
    thread_state_t *tstate_ = (thread_state_t *) pthread_getspecific(state->thread_data); \
    delayed_count_add(& tstate_->count_updates, & (ptr), (amt));            \
  }                                                                         \
} while (0)

#define UPDATE_COUNT_DBL(state, ptr, amt) do {                                     \
  if (state->sequential) (ptr) += (amt);                                    \
  else {                                                                    \
    thread_state_t *tstate_ = (thread_state_t *) pthread_getspecific(state->thread_data); \
    delayed_count_add_double(& tstate_->count_updates, & (ptr), (amt));     \
  }                                                                         \
} while (0)


typedef struct delayed_count_buffer {
    struct pointer_hash     counts;
    struct pointer_hash     counts_dbl;
    double                 *buffer;
    double                 *fill;
    double                 *limit;
} delayed_count_buffer_t;


int delayed_count_init(delayed_count_buffer_t *buf,
                       int init_size);

void delayed_count_destroy(delayed_count_buffer_t *buf);

int delayed_count_add(delayed_count_buffer_t *buf, int *target, int  amount);

int delayed_count_add_double(delayed_count_buffer_t *buf, double *target,
                             double  amount);

void delayed_count_play(delayed_count_buffer_t *buf);

#endif
