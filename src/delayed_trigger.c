#include "mcell_structs.h"
#include "delayed_trigger.h"
#include "mem_util.h"
#include "count_util.h"
#include <stdlib.h>

extern struct volume *world;

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

void delayed_trigger_fire(delayed_trigger_buffer_t *buf,
                          struct counter *event,
                          int n,
                          struct vector3 const *where,
                          unsigned char what)
{
  if (buf->fill == buf->length)
    delayed_trigger_flush(buf, 0);

  int const idx = buf->fill ++;
  buf->triggers[idx].event = event;
  buf->triggers[idx].n     = n;
  buf->triggers[idx].where = *where;
  buf->triggers[idx].what  = what;
}

void delayed_trigger_flush(delayed_trigger_buffer_t *buf,
                           int nolock)
{
  /* Acquire the global triggers lock. */
  if (! nolock)
    pthread_mutex_lock(& world->trig_lock);

  /* Spill the triggers. */
  for (int i=0; i<buf->fill; ++i)
    fire_count_event(buf->triggers[i].event,
                     buf->triggers[i].n,
                     &buf->triggers[i].where,
                     buf->triggers[i].what);
  buf->fill = 0;

  /* Release the global triggers lock. */
  if (! nolock)
    pthread_mutex_unlock(& world->trig_lock);
}
