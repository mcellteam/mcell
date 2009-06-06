#ifndef TESTUTILS_H
#define TESTUTILS_H

#include <stdio.h>
#include <string.h>

/* Creating a test suite:
 *    1. create test functions like:
 *
 *          int test_circular_buffer_1(void)
 *          {
 *            FAIL_UNLESS_EQUAL(1, 2);
 *            return 0;
 *          }
 *
 *       NOTE: the function must return an integer, and you must include a
 *       return 0 at the bottom of the function.  The rest of the function
 *       should include at least one test for failure, which is done using the
 *       FAIL* macros:
 *           FAIL(msg)
 *           FAIL_IF(cond)
 *           FAIL_UNLESS(cond)
 *           FAIL_IF_EQUAL(a, b)
 *           FAIL_UNLESS_EQUAL(a, b)
 *           FAIL_IF_STREQUAL(a, b)
 *           FAIL_UNLESS_STREQUAL(a, b)
 *
 *       See below for more details on these macros.
 *
 *    2. Create a top-level function like:
 *          int suite_circular_buffer(void)
 *          {
 *            TEST_PROLOGUE("Circular buffer tests");
 *            RUNTEST(test_circular_buffer_1);
 *            RUNTEST(test_circular_buffer_2);
 *            RUNTEST(test_circular_buffer_3);
 *            ..
 *            RUNTEST(test_circular_buffer_N);
 *            TEST_EPILOGUE;
 *          }
 *
 *       NOTE: the function must return an integer, but you should NOT
 *       explicitly return from it.  The return value will be handled via
 *       TEST_EPILOGUE.
 */

/* Prologue: Put this at the beginning of a test suite. */
#define TEST_PROLOGUE(label)                                                \
      char const *_tu_testname = (label);                                   \
      int _tu_numtests = 0, _tu_numerrs = 0;                                \
      fprintf(stderr, "Running suite: %s\n", _tu_testname);                 \
      fflush(stderr);

/* Epilogue: Put this at the end of a test suite. */
#define TEST_EPILOGUE do {                                                  \
      if (_tu_numerrs == 0) {                                               \
        printf("\n%s: Success (%d tests)\n", _tu_testname, _tu_numtests);   \
        return 0;                                                           \
      }                                                                     \
      else {                                                                \
        printf("\n%s: Ran %d tests, %d failures.\n", _tu_testname,          \
               _tu_numtests, _tu_numerrs);                                  \
        return 1;                                                           \
      }                                                                     \
} while(0)

/* Content: Put one of these for each test in a test suite. */
#define RUNTEST(func) do {                                                  \
  ++ _tu_numtests;                                                          \
  if (func()) {                                                             \
    putchar('F');                                                           \
    _tu_numerrs += 1;                                                       \
  }                                                                         \
  else {                                                                    \
    putchar('.');                                                           \
  }                                                                         \
  fflush(stdout);                                                           \
} while(0)

#define STRINGIFY(a) #a

/* Unconditionally fail a test. */
#define FAIL(msg) do {                                                      \
  fprintf(stderr, "Test %s failed: %s", __func__, msg);                     \
  return 1;                                                                 \
} while (0)

/* Fail a test if the condition is true. */
#define FAIL_IF(a) do {                                                     \
  if ((a)) {                                                                \
    fprintf(stderr, "Test %s failed: %s", __func__, STRINGIFY(a));          \
    return 1;                                                               \
  }                                                                         \
} while (0)

/* Fail a test if the condition is false. */
#define FAIL_UNLESS(a) do {                                                 \
  if (! (a)) {                                                              \
    fprintf(stderr, "Test %s failed: ! %s", __func__, STRINGIFY(a));        \
    return 1;                                                               \
  }                                                                         \
} while (0)

/* Fail a test if the arguments are equal. */
#define FAIL_IF_EQUAL(a, b) do {                                            \
  if ((a) == (b)) {                                                         \
    fprintf(stderr, "Test %s failed: %s == %s",                             \
            __func__, STRINGIFY(a), STRINGIFY(b));                          \
    return 1;                                                               \
  }                                                                         \
} while (0)

/* Fail a test if the arguments are unequal. */
#define FAIL_UNLESS_EQUAL(a, b) do {                                        \
  if ((a) != (b)) {                                                         \
    fprintf(stderr, "Test %s failed: %s != %s",                             \
            __func__, STRINGIFY(a), STRINGIFY(b));                          \
    return 1;                                                               \
  }                                                                         \
} while (0)

/* Fail a test if the arguments are equal according to strcmp. (This is
 * NULL-safe, so you may call it with NULL values for either or both of the
 * arguments, and it will "do the right thing."
 */
#define FAIL_IF_STREQUAL(a, b) do {                                         \
  char *atmp = (a), *btmp = (b);                                            \
  if (btmp == NULL) {                                                       \
    if (atmp == NULL) {                                                     \
      fprintf(stderr, "Test %s failed: %s == %s",                           \
              __func__, STRINGIFY(a), STRINGIFY(b));                        \
      return 1;                                                             \
    }                                                                       \
  }                                                                         \
  else {                                                                    \
    if (atmp != NULL  &&  ! strcmp(atmp, btmp)) {                           \
      fprintf(stderr, "Test %s failed: %s == %s",                           \
              __func__, STRINGIFY(a), STRINGIFY(b));                        \
      return 1;                                                             \
    }                                                                       \
  }                                                                         \
} while (0)

/* Fail a test if the arguments are unequal according to strcmp. (This is
 * NULL-safe, so you may call it with NULL values for either or both of the
 * arguments, and it will "do the right thing."
 */
#define FAIL_UNLESS_STREQUAL(a, b) do {                                     \
  char *atmp = (a), *btmp = (b);                                            \
  if (btmp == NULL) {                                                       \
    if (atmp != NULL) {                                                     \
      fprintf(stderr, "Test %s failed: %s != %s",                           \
              __func__, STRINGIFY(a), STRINGIFY(b));                        \
      return 1;                                                             \
    }                                                                       \
  }                                                                         \
  else {                                                                    \
    if (atmp == NULL  ||  strcmp(atmp, btmp)) {                             \
      fprintf(stderr, "Test %s failed: %s != %s",                           \
              __func__, STRINGIFY(a), STRINGIFY(b));                        \
      return 1;                                                             \
    }                                                                       \
  }                                                                         \
} while (0)

#endif
