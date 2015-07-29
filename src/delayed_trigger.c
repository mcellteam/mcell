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

#include "config.h"

#include <stdlib.h>

#include "mcell_structs.h"
#include "delayed_trigger.h"
#include "mem_util.h"
#include "count_util.h"


void delayed_trigger_init(delayed_trigger_buffer_t *buf,
                          int capacity)
{
  buf->fill = 0;
  buf->length = capacity;
  buf->triggers = CHECKED_MALLOC_ARRAY(delayed_trigger_t,
                                       capacity,
                                       "delayed triggers");
}

void delayed_trigger_destroy(delayed_trigger_buffer_t *buf)
{
  if (buf->triggers != NULL)
    free(buf->triggers);
  buf->triggers = NULL;
  buf->fill = buf->length = 0;
}

void delayed_trigger_fire(struct volume *state,
  delayed_trigger_buffer_t *buf, struct counter *event, int n,
  struct vector3 const *where, unsigned char what) {

  if (buf->fill == buf->length)
    delayed_trigger_flush(state, buf, 0);

  int const idx = buf->fill ++;
  buf->triggers[idx].event = event;
  buf->triggers[idx].n     = n;
  buf->triggers[idx].where = *where;
  buf->triggers[idx].what  = what;
}

void delayed_trigger_flush(struct volume *state, delayed_trigger_buffer_t *buf,
  int nolock)
{
  /* Acquire the global triggers lock. */
  if (! nolock)
    pthread_mutex_lock(&state->trig_lock);

  /* Spill the triggers. */
  for (int i=0; i<buf->fill; ++i)
    fire_count_event(state, buf->triggers[i].event, buf->triggers[i].n,
      &buf->triggers[i].where, buf->triggers[i].what);
  buf->fill = 0;

  /* Release the global triggers lock. */
  if (! nolock)
    pthread_mutex_unlock(&state->trig_lock);
}
