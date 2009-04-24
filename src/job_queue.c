#include "job_queue.h"
#include "mem_util.h"

job_queue_t *jobq_create(int capacity)
{
  /* Allocate the job queue itself. */
  job_queue_t *jobq = CHECKED_MALLOC_STRUCT_NODIE(job_queue_t, "job queue");
  if (jobq == NULL)
    return NULL;

  /* Allocate a circular buffer for the tasks. */
  if ((jobq->buffer = cbuf_create(capacity)) == NULL)
  {
    free(jobq);
    return NULL;
  }

  /* Initialize the mutex/condition variable. */
  pthread_mutex_init(&jobq->lock, NULL);
  pthread_cond_init(&jobq->signal, NULL);

  jobq->shutdown = 0;

  return jobq;
}

void jobq_destroy(job_queue_t *jobq,
                  void (*cleanup)(void *task, void *ctx),
                  void *ctx)
{
  /* Destroy the synchronization primitives. */
  pthread_mutex_destroy(&jobq->lock);
  pthread_cond_destroy(&jobq->signal);

  /* Destroy the circular buffer. */
  cbuf_destroy(jobq->buffer, cleanup, ctx);
  free(jobq);
}

void jobq_shutdown(job_queue_t *jobq)
{
  pthread_mutex_lock(& jobq->lock);
  jobq->shutdown = 1;
  pthread_cond_broadcast(& jobq->signal);
  pthread_mutex_unlock(& jobq->lock);
}

int jobq_add_task(job_queue_t *jobq, void *task)
{
  int ntasks = 0;

  /* Lock, store, signal, unlock. */
  pthread_mutex_lock(& jobq->lock);
  ntasks = cbuf_store_item(jobq->buffer, task);
  if (ntasks)
    pthread_cond_signal(& jobq->signal);
  pthread_mutex_unlock(& jobq->lock);

  return ntasks;
}

int jobq_add_tasks(job_queue_t *jobq, void **task, int ntasks)
{
  int ntasks_added = 0;

  /* Lock, store, signal, unlock. */
  pthread_mutex_lock(& jobq->lock);

  /* Store tasks until we fail. */
  for (int i=0; i<ntasks; ++i)
  {
    if (cbuf_store_item(jobq->buffer, task[i]))
      ++ ntasks_added;
    else
      break;
  }

  /* Signal one or many, depending on number of tasks added. */
  if (ntasks_added > 1)
    pthread_cond_broadcast(& jobq->signal);
  else if (ntasks_added == 1)
    pthread_cond_signal(& jobq->signal);

  pthread_mutex_unlock(& jobq->lock);

  return ntasks_added;
}

void *jobq_wait_task(job_queue_t *jobq)
{
  void *task = NULL;

  pthread_mutex_lock(& jobq->lock);

  /* Loop until we get a task. */
  while (1)
  {
    /* Return if shutdown was requested. */
    if (jobq->shutdown)
      break;

    /* Grab the next item and return if we got one. */
    if ((task = cbuf_take_item(jobq->buffer)) != NULL)
      break;

    /* Release the lock and wait for our turn. */
    pthread_cond_wait(& jobq->signal, & jobq->lock);
  }

  pthread_mutex_unlock(& jobq->lock);

  return task;
}
