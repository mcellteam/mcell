#ifndef DELAYED_TRIGGER_H
#define DELAYED_TRIGGER_H

#include "vector.h"
#include <pthread.h>

#define UPDATE_TRIGGER(event, n, where, what) do {                          \
  if (world->non_parallel) fire_count_event((event), (n), (where), (what)); \
  else {                                                                    \
    thread_state_t *tstate_ = (thread_state_t *) pthread_getspecific(world->thread_data); \
    delayed_trigger_fire(& tstate_->triggers, (event), (n), (where), (what)); \
  }                                                                         \
} while (0)

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

void delayed_trigger_init(delayed_trigger_buffer_t *buf,
                          int capacity);

void delayed_trigger_destroy(delayed_trigger_buffer_t *buf);

void delayed_trigger_fire(delayed_trigger_buffer_t *buf,
                          struct counter *event,
                          int n,
                          struct vector3 const *where,
                          unsigned char what);

void delayed_trigger_flush(delayed_trigger_buffer_t *buf,
                           int nolock);

#endif
