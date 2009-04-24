#include "circular_buffer.h"
#include "mem_util.h"
#include <string.h>

#ifndef offsetof
#define offsetof(t, f) ((int) &((t *) NULL)->f)
#endif

circular_buffer_t *cbuf_create(int capacity)
{
  if (capacity <= 0)
    capacity = CIRCULAR_BUFFER_DEFAULT_CAPACITY;

  /* Allocate the buffer. */
  size_t size = offsetof(circular_buffer_t, entries[0]) + capacity * sizeof(void *);
  circular_buffer_t *cb = CHECKED_MALLOC_NODIE(size, "circular buffer");
  if (cb == NULL)
    return NULL;

  /* Initialize the buffer. */
  memset(cb, 0, size);
  cb->size = capacity;
  return cb;
}

void cbuf_destroy(circular_buffer_t *cb,
                  void (*cleanup)(void *entry, void *ctx),
                  void *ctx)
{
  if (cleanup)
    cbuf_clear(cb, cleanup, ctx);
  free(cb);
}

void cbuf_clear(circular_buffer_t *cb,
                void (*cleanup)(void *entry, void *ctx),
                void *ctx)
{
  if (cleanup)
  {
    while (cb->fill != 0)
    {
      void *item = cbuf_take_item(cb);
      (*cleanup)(item, ctx);
    }
  }
  else
    cb->fill = cb->read_pos = 0;
}

void *cbuf_take_item(circular_buffer_t *cb)
{
  if (cb->fill == 0)
    return NULL;

  void *item = cb->entries[cb->read_pos];
  if (++ cb->read_pos == cb->size)
    cb->read_pos = 0;
  -- cb->fill;
  return item;
}

void *cbuf_peek_item(circular_buffer_t *cb)
{
  if (cb->fill == 0)
    return NULL;

  return cb->entries[cb->read_pos];
}

int cbuf_store_item(circular_buffer_t *cb, void *item)
{
  if (cb->fill >= cb->size)
    return 0;

  int write_pos = cb->read_pos + cb->fill;
  if (write_pos >= cb->size)
    write_pos -= cb->size;
  cb->entries[write_pos] = item;

  ++ cb->fill;
  return 1;
}
