#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <sys/types.h>
#include <sys/wait.h>

/* Steps to add a new test suite:
 *    1. Add the test suite filename to Makefile.am
 *    2. Add a forward declaration for the top-level test suite function
 *    3. Add the function to the array below.
 *
 * See, for instance, test_circular_buffer.c for details on how to add unit
 * tests.
 */

/* Add test suites to this list. */
int (*ALL_TEST_SUITES[])(void) = {
};

/********************************************************************
 * Private implementation below.  You shouldn't need to modify below this line
 * unless there is a bug in the testing machinery.
 *******************************************************************/

/* Test suite runner: Fork a new process, redirect stderr to a file, and run
 * the test function, returning 1 if any error is encountered.
 */
static int run_suite(int (*testfunc)(void))
{
  /* Fork a child process to run the tests. */
  pid_t pid = fork();
  if (pid < 0)
  {
    printf("Failure: couldn't fork to run test suite...\n");
    return 1;
  }

  /* If we are the parent, wait for the child's response. */
  if (pid != 0)
  {
    /* Wait for the test to finish. */
    int status = 0;
    if (waitpid(pid, & status, 0) <= 0)
    {
      perror("Warning: waitpid failed on test suite");
      return 1;
    }

    /* If the test exited normally... */
    if (WIFEXITED(status))
    {
      if (WEXITSTATUS(status) == 0)
        return 0;
      else
        return 1;
    }

    /* If the test died due to a signal... */
    else if (WIFSIGNALED(status))
    {
      printf("Warning: A test suite failed with signal %d\n",
             WTERMSIG(status));
      return 1;
    }

    /* You are now entering... The Twilight Zone. */
    else
    {
      printf("Warning: A test suite finished without exiting or "
             "dying due to a signal.\nThis is highly irregular.\n");
      return 1;
    }
  }

  /* Open the logfile, and attach it to stderr. */
  FILE *f = fopen("unittest.log", "a");
  if (f == NULL)
  {
    fprintf(stderr, "Failed to open \"unittest.log\".\nPlease make sure the"
            "test suite is being run from a writable directory.\n");
    exit(1);
  }
  dup2(fileno(f), 2);

  /* Call the test. */
  if ((*testfunc)())
  {
    fclose(f);   /* Unnecessary, but safe. */
    exit(1);
  }
  else
  {
    fclose(f);   /* Unnecessary, but safe. */
    exit(0);
  }
}

#define COUNT_OF(a) (sizeof(a) / sizeof(a[0]))

int main(int c, char **v)
{
  /* Remove the old unittest.log */
  unlink("unittest.log");

  /* Loop over the tests, running them. */
  int nfailed = 0;
  for (unsigned int i=0; i<COUNT_OF(ALL_TEST_SUITES); ++i)
    nfailed += run_suite(ALL_TEST_SUITES[i]);

  /* Print a final summary. */
  if (nfailed == 0)
  {
    printf("All suites passed successfully.\n");
    return 0;
  }
  else
  {
    printf("%d test suites exhibited failure.\n", nfailed);
    return 1;
  }
}
