#include <stdlib.h>
#include <stdio.h>
#define __USE_GNU 1
#include <dlfcn.h>
#include <string.h>
#include <time.h>
#include <sys/time.h>
#include <unistd.h>

/* Pointer to next malloc function */
typedef void *(*mallocfunc)(size_t);

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

/* Intercept malloc, injecting failures as appropriate */
void *malloc(size_t size)
{
  static mallocfunc next_malloc = NULL;
  static int fail_interval = 0;
  static int fail_random = 0;
  static int ntimes = 0;

  if (next_malloc == NULL)
  {
    next_malloc = dlsym(RTLD_NEXT, "malloc");
    fail_interval = get_fail_interval();
    fail_random = get_fail_random();
    printf("Initialized malloc failure injector (%d, %d)\n", fail_interval, fail_random);
  }

  if (fail_random)
  {
    if ((rand() % fail_interval) == 0)
    {
      printf("Injected random malloc failure on request for %u bytes\n", size);
#ifdef MALLOC_INJECTOR_BREAK_TO_GDB
      char buffer[2048];
      int pid = getpid();
      snprintf(buffer, 2048, "gdb /proc/%d/exe %d", pid, pid);
      unsetenv("LD_PRELOAD");
      system(buffer);
#endif
      return NULL;
    }
  }
  else if (fail_interval > 0)
  {
    if (++ ntimes == fail_interval)
    {
      ntimes = 0;
      printf("Injected malloc failure on request for %u bytes\n", size);
#ifdef MALLOC_INJECTOR_BREAK_TO_GDB
      char buffer[2048];
      int pid = getpid();
      snprintf(buffer, 2048, "gdb /proc/%d/exe %d", pid, pid);
      unsetenv("LD_PRELOAD");
      system(buffer);
#endif
      return NULL;
    }
  }

  return (*next_malloc)(size);
}
