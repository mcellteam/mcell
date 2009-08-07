#ifndef DELAYED_COUNT_H
#define DELAYED_COUNT_H

#include "util.h"

#define UPDATE_COUNT(ptr, amt) do {                                         \
  if (world->non_parallel) (ptr) += (amt);                                  \
  else {                                                                    \
    thread_state_t *tstate_ = (thread_state_t *) pthread_getspecific(world->thread_data); \
    delayed_count_add(& tstate_->count_updates, & (ptr), (amt));            \
  }                                                                         \
} while (0)

#define UPDATE_COUNT_DBL(ptr, amt) do {                                     \
  if (world->non_parallel) (ptr) += (amt);                                  \
  else {                                                                    \
    thread_state_t *tstate_ = (thread_state_t *) pthread_getspecific(world->thread_data); \
    delayed_count_add_double(& tstate_->count_updates, & (ptr), (amt));     \
  }                                                                         \
} while (0)

typedef struct delayed_count_buffer
{
    struct pointer_hash     counts;
    struct pointer_hash     counts_dbl;
    double                 *buffer;
    double                 *fill;
    double                 *limit;
} delayed_count_buffer_t;

int delayed_count_init(delayed_count_buffer_t *buf,
                       int init_size);

void delayed_count_destroy(delayed_count_buffer_t *buf);

int delayed_count_add(delayed_count_buffer_t *buf,
                      int *target,
                      int  amount);
int delayed_count_add_double(delayed_count_buffer_t *buf,
                             double *target,
                             double  amount);

void delayed_count_play(delayed_count_buffer_t *buf);

#endif
