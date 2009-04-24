#ifndef JOB_QUEUE_H
#define JOB_QUEUE_H

#include "circular_buffer.h"

#include <pthread.h>

typedef struct job_queue
{
  circular_buffer_t      *buffer;
  pthread_mutex_t         lock;
  pthread_cond_t          signal;
  int                     shutdown;
} job_queue_t;

/* Create a new job queue with a given capacity. */
job_queue_t *jobq_create(int capacity);

/* Destroy a job queue.  If a cleanup function is given, it will be called on
 * all currently outstanding tasks in the queue before they are deleted. */
void jobq_destroy(job_queue_t *jobq,
                  void (*cleanup)(void *task, void *ctx),
                  void *ctx);

/* Signal a shutdown. */
void jobq_shutdown(job_queue_t *jobq);

/* Add a task to the job queue.  Returns 1 if there was sufficient space for
 * the new task, 0 otherwise. Queue is left in a valid state even if the task
 * couldn't be added. */
int jobq_add_task(job_queue_t *jobq, void *task);

/* Add several tasks to the job queue. Returns the number of tasks which were
 * successfully added.  The queue will be left in a valid state even if not all
 * tasks could be added. */
int jobq_add_tasks(job_queue_t *jobq, void **task, int ntasks);

/* Wait for any task to become available. */
void *jobq_wait_task(job_queue_t *jobq);

#endif
