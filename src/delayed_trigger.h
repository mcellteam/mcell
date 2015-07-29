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

#ifndef DELAYED_TRIGGER_H
#define DELAYED_TRIGGER_H

#include "config.h"

#include <pthread.h>

#include "vector.h"
#include "mcell_structs.h"


#define UPDATE_TRIGGER(state, event, n, where, what) do {                          \
  if (world->sequential) fire_count_event((state), (event), (n), (where), (what));   \
  else {                                                                    \
    thread_state_t *tstate_ = (thread_state_t *) pthread_getspecific(state->thread_data); \
    delayed_trigger_fire(state, & tstate_->triggers, (event), (n), (where), (what)); \
  }                                                                         \
} while (0)

#if 0
typedef struct delayed_trigger
{
  struct counter *event;
  int             n;
  struct vector3  where;
  unsigned char   what;
} delayed_trigger_t;

typedef struct delayed_trigger_buffer
{
  delayed_trigger_t      *triggers;
  int                     fill;
  int                     length;
} delayed_trigger_buffer_t;
#endif

void delayed_trigger_init(delayed_trigger_buffer_t *buf,
                          int capacity);

void delayed_trigger_destroy(delayed_trigger_buffer_t *buf);

void delayed_trigger_fire(struct volume *world, delayed_trigger_buffer_t *buf,
                          struct counter *event,
                          int n,
                          struct vector3 const *where,
                          unsigned char what);

void delayed_trigger_flush(struct volume *world, delayed_trigger_buffer_t *buf,
                           int nolock);

#endif
