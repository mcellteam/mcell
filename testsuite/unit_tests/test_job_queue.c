#include <stdio.h>
#include <unistd.h>

#include "testutils.h"
#include "job_queue.h"

/* Mutex protecting "rand", which is not threadsafe, as well as the task
 * completion counter. */
pthread_mutex_t task_mutex = PTHREAD_MUTEX_INITIALIZER;
int task_completion_count = 0;

/* Simple task which does nothing except resubmit with a probability of 0.5.
 * If N tasks are submitted, we should expect approximately 2*N will be
 * executed before all tasks are done.
 */
static void run_task(job_queue_t *jobq, void *task)
{
  pthread_mutex_lock(& task_mutex);
  int randval = (int)((2.0 * rand()) / RAND_MAX);
  ++ task_completion_count;
  pthread_mutex_unlock(& task_mutex);
  if (randval != 0)
    jobq_add_task(jobq, "dummy task");
}

/* Worker parameters structure to contain thread startup state. */
typedef struct worker_params
{
  job_queue_t *jobq;
  int          sleepytime;
} worker_params_t;

/* The worker main loop. */
static void *worker_thread(void *param)
{
  worker_params_t *params = (worker_params_t *) param;
  job_queue_t *jobq = params->jobq;
  void *task = NULL;
  while ((task = jobq_wait_task(jobq)) != NULL)
  {
    sleep(params->sleepytime);
    run_task(jobq, task);
  }
  return NULL;
}

/* Start up the requested number of workers with the provided parameters. */
static void start_workers(pthread_t *threads,
                          int nthreads,
                          int sleepytime,
                          job_queue_t *jobq)
{
  static worker_params_t params;
  params.jobq = jobq;
  params.sleepytime = sleepytime;
  for (int i=0; i<nthreads; ++i)
    pthread_create(&threads[i], NULL, worker_thread, (void *) &params);
}

/* Wait for the shutdown of the given number of threads */
static void join_workers(pthread_t *threads,
                         int nthreads)
{
  void *dummy = NULL;
  for (int i=0; i<nthreads; ++i)
    pthread_join(threads[i], &dummy);
}

int test_job_queue_1(void)
{
  task_completion_count = 0;

  job_queue_t *jobq = jobq_create(1024);

  /* Fill queue with tasks. */
  for (size_t i=0; i<1024; ++ i)
    jobq_add_task(jobq, "initial task");

  /* Start 100 workers, tasks take 1 second. */
  pthread_t workers[100];
  start_workers(workers, 100, 1, jobq);

  /* Wait for the workers to finish. */
  join_workers(workers, 100);

  /* Signal a shutdown. */
  jobq_shutdown(jobq);
  join_workers(workers, 100);

  jobq_destroy(jobq, NULL, NULL);

  /* We expect at least 1024 tasks.  It's a pretty safe bet that this will be
   * around 2048, actually, but if this test fails, we're likely to never get
   * here anyway...
   */
  FAIL_UNLESS(task_completion_count >= 1024);
  return 0;
}

int test_job_queue_2(void)
{
  task_completion_count = 0;

  job_queue_t *jobq = jobq_create(1024);

  /* Fill queue with tasks. */
  jobq_add_task(jobq, "initial task");

  /* Start 100 workers, tasks take 5 seconds. */
  pthread_t workers[100];
  start_workers(workers, 100, 5, jobq);

  /* Wait for the workers to finish. */
  join_workers(workers, 100);

  /* Signal a shutdown. */
  jobq_shutdown(jobq);
  join_workers(workers, 100);

  jobq_destroy(jobq, NULL, NULL);

  /* We expect at least 1 task completed.
   */
  FAIL_UNLESS(task_completion_count >= 1);
  return 0;
}

int suite_job_queue(void)
{
  TEST_PROLOGUE("Job Queue Tests");
  RUNTEST(test_job_queue_1);
  RUNTEST(test_job_queue_2);
  TEST_EPILOGUE;
}
