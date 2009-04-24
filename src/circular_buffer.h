#ifndef INCLUDED_CIRCULAR_BUFFER_H
#define INCLUDED_CIRCULAR_BUFFER_H

#include <stdlib.h>

/* Default capacity if requested capacity is <= 0. */
enum {
  CIRCULAR_BUFFER_DEFAULT_CAPACITY = 4096
};

typedef struct circular_buffer
{
  int     size;
  int     fill;
  int     read_pos;
  void   *entries[1];
} circular_buffer_t;

/* Allocate a new circular buffer with a given capacity. */
circular_buffer_t *cbuf_create(int capacity);

/* Deallocate a circular buffer, calling an optional cleanup function on each
 * leftover item in the buffer. */
void cbuf_destroy(circular_buffer_t *cb,
                  void (*cleanup)(void *entry, void *ctx),
                  void *ctx);

/* Empty a circular buffer, calling an optional cleanup function on each
 * leftover item in the buffer. */
void cbuf_clear(circular_buffer_t *cb,
                void (*cleanup)(void *entry, void *ctx),
                void *ctx);

/* Take the next item out of the circular buffer.  Returns NULL if no item was
 * found. */
void *cbuf_take_item(circular_buffer_t *cb);

/* Peek at the next item in the circular buffer.  Returns NULL if no item was
 * found. */
void *cbuf_peek_item(circular_buffer_t *cb);

/* Stash an item into the circular buffer. Returns 0 if the buffer was full
 * (i.e. the write failed), and 1 otherwise. */
int cbuf_store_item(circular_buffer_t *cb, void *item);

/* Check if this buffer is empty. */
#define cbuf_empty(cb) ((cb)->fill == 0)

/* Get the maximum number of items in the circular buffer. */
#define cbuf_size(cb) ((cb)->size)

/* Get the number of items in the circular buffer. */
#define cbuf_count(cb) ((cb)->fill)

/* Get the number of empty slots in the circular buffer (i.e. the number of
 * additional items that may be immediately stored).  CAUTION: macro evaluates
 * 'cb' more than once. */
#define cbuf_free_space(cb) ((cb)->size - (cb)->fill)

#endif
