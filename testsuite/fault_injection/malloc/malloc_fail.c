#include <stdlib.h>
#include <stdio.h>
#define __USE_GNU 1
#include <dlfcn.h>
#include <string.h>
#include <time.h>
#include <sys/time.h>
#include <unistd.h>
#include <stdarg.h>

/* Pointer to next malloc function */
typedef void *(*mallocfunc)(size_t);
static mallocfunc  next_malloc = NULL;

/* Pointer to next realloc function */
typedef void *(*reallocfunc)(void *, size_t);
static reallocfunc next_realloc = NULL;

/* Pointer to next calloc function */
typedef void *(*callocfunc)(size_t, size_t);
static callocfunc  next_calloc = NULL;

/* Interval at which to fail */
static int fail_interval = 0;

/* Fail randomly or deterministically? */
static int fail_random = 0;

/* Counter tracking how many successfully allocations we've had since the last
 * injected failure */
static int ntimes = 0;

/* Break to debugger */
static void break_to_gdb()
{
#ifdef MALLOC_INJECTOR_BREAK_TO_GDB
  int tmp;
  tmp = fail_interval;
  fail_interval = -1;

  char buffer[2048];
  int pid = getpid();
  snprintf(buffer, 2048, "gdb /proc/%d/exe %d", pid, pid);
  unsetenv("LD_PRELOAD");
  system(buffer);
  fail_interval = tmp;
#endif
}

/* Grab our "randomness" setting from the environment */
static int get_fail_random()
{
  char *str_itvl = getenv("INJECT_MALLOC_FAIL_RANDOM");
  if (str_itvl == NULL)
    return 0;
  else if (! strcmp(str_itvl, "0"))
    return 0;
  else
  {
    srand(time(NULL));
    return 1;
  }
}

/* Grab our "interval" setting from the environment */
static int get_fail_interval()
{
  int value = 0;
  char *str_itvl = getenv("INJECT_MALLOC_FAIL_INTERVAL");
  if (str_itvl != NULL)
    value = atoi(str_itvl);
  if (value == 0)
    return -1;
  return value;
}

static void prepare()
{
  if (next_malloc == NULL)
  {
    next_malloc  = dlsym(RTLD_NEXT, "malloc");
    next_realloc = dlsym(RTLD_NEXT, "realloc");
    next_calloc  = dlsym(RTLD_NEXT, "calloc");
    fail_interval = get_fail_interval();
    fail_random = get_fail_random();
    printf("Initialized malloc failure injector (%d, %d)\n", fail_interval, fail_random);
  }
}

static int time_to_fail()
{
  if (fail_interval > 0)
  {
    if (fail_random)
    {
      if ((rand() % fail_interval) == 0)
        return 1;
    }
    else
    {
      if (++ ntimes == fail_interval)
      {
        ntimes = 0;
        return 1;
      }
    }
  }
  return 0;
}

/* Intercept malloc, injecting failures as appropriate */
void *malloc(size_t size)
{
  prepare();

  if (time_to_fail())
  {
    if (fail_random)
      printf("Injected random malloc failure on request for %u bytes\n", size);
    else
      printf("Injected malloc failure on request for %u bytes\n", size);
    break_to_gdb();
    return NULL;
  }

  return (*next_malloc)(size);
}

/* Intercept realloc, injecting failures as appropriate */
void *realloc(void *old, size_t size)
{
  prepare();

  if (time_to_fail())
  {
    if (fail_random)
      printf("Injected random realloc failure on request for %u bytes\n", size);
    else
      printf("Injected realloc failure on request for %u bytes\n", size);
    break_to_gdb();
    return NULL;
  }

  return (*next_realloc)(old, size);
}

/* Intercept calloc, injecting failures as appropriate */
void *calloc(size_t nmemb, size_t size)
{
  prepare();

  if (time_to_fail())
  {
    if (fail_random)
      printf("Injected random calloc failure on request for %u*%u=%u bytes\n", nmemb, size, nmemb*size);
    else
      printf("Injected calloc failure on request for %u*%u=%u bytes\n", nmemb, size, nmemb*size);
    break_to_gdb();
    return NULL;
  }

  return (*next_calloc)(nmemb, size);
}
