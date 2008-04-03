#include <stdlib.h>
#include <stdio.h>
#define __USE_GNU 1
#include <dlfcn.h>
#include <string.h>
#include <time.h>
#include <sys/time.h>

typedef void *(*mallocfunc)(size_t);

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

static int get_fail_interval()
{
  int value = 0;
  char *str_itvl = getenv("INJECT_MALLOC_FAIL_INTERVAL");
  if (str_itvl != NULL)
    value = atoi(str_itvl);
  if (value == 0)
    return -1;
}

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
      return NULL;
  }
  else if (fail_interval > 0)
  {
    if (++ ntimes == fail_interval)
    {
      ntimes = 0;
      return NULL;
    }
  }

  return (*next_malloc)(size);
}
