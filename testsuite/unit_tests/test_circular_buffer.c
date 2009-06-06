#include <stdio.h>

#include "testutils.h"
#include "circular_buffer.h"

int test_circular_buffer_empty(void)
{
  /* Create empty buffer */
  circular_buffer_t *cb = cbuf_create(10);

  /* Verify vital stats */
  FAIL_UNLESS_EQUAL(cbuf_size(cb), 10);
  FAIL_UNLESS_EQUAL(cbuf_count(cb), 0);
  FAIL_UNLESS_EQUAL(cbuf_free_space(cb), 10);
  FAIL_UNLESS(cbuf_empty(cb));

  /* Verify that taking/peeking fail. */
  FAIL_UNLESS_EQUAL(cbuf_take_item(cb), NULL);
  FAIL_UNLESS_EQUAL(cbuf_peek_item(cb), NULL);

  /* Cleanup. */
  cbuf_destroy(cb, NULL, NULL);
  return 0;
}

int test_circular_buffer_store(void)
{
  /* Create buffer with one item. */
  circular_buffer_t *cb = cbuf_create(10);
  FAIL_UNLESS(cbuf_store_item(cb, "apple"));

  /* Verify vital stats */
  FAIL_IF(cbuf_empty(cb));
  FAIL_UNLESS_EQUAL(cbuf_size(cb), 10);
  FAIL_UNLESS_EQUAL(cbuf_count(cb), 1);
  FAIL_UNLESS_EQUAL(cbuf_free_space(cb), 9);
  FAIL_IF_EQUAL(cbuf_peek_item(cb), NULL);

  /* Verify taking/peeking. */
  FAIL_UNLESS_STREQUAL(cbuf_peek_item(cb), "apple");
  FAIL_IF(cbuf_empty(cb));
  FAIL_UNLESS_STREQUAL(cbuf_take_item(cb), "apple");
  FAIL_UNLESS(cbuf_empty(cb));

  /* Cleanup. */
  cbuf_destroy(cb, NULL, NULL);
  return 0;
}

int test_circular_buffer_full(void)
{
  /* Create buffer with ten items. */
  circular_buffer_t *cb = cbuf_create(10);
  FAIL_UNLESS(cbuf_store_item(cb, "apple"));
  FAIL_UNLESS(cbuf_store_item(cb, "banana"));
  FAIL_UNLESS(cbuf_store_item(cb, "cherry"));
  FAIL_UNLESS(cbuf_store_item(cb, "durian"));
  FAIL_UNLESS(cbuf_store_item(cb, "elderberry"));
  FAIL_UNLESS(cbuf_store_item(cb, "fig"));
  FAIL_UNLESS(cbuf_store_item(cb, "guava"));
  FAIL_UNLESS(cbuf_store_item(cb, "huckleberry"));
  FAIL_UNLESS(cbuf_store_item(cb, "imbe"));
  FAIL_UNLESS(cbuf_store_item(cb, "jujube"));

  /* Verify that the next store fails. */
  FAIL_IF(    cbuf_store_item(cb, "kumquat"));

  /* Verify vital stats. */
  FAIL_UNLESS_EQUAL(cbuf_size(cb), 10);
  FAIL_UNLESS_EQUAL(cbuf_count(cb), 10);
  FAIL_UNLESS_EQUAL(cbuf_free_space(cb), 0);

  /* Verify peeking/taking */
  FAIL_UNLESS_STREQUAL(cbuf_peek_item(cb), "apple");
  FAIL_UNLESS_STREQUAL(cbuf_take_item(cb), "apple");
  FAIL_UNLESS_STREQUAL(cbuf_peek_item(cb), "banana");
  FAIL_UNLESS_STREQUAL(cbuf_take_item(cb), "banana");

  /* Try writing a few more items. */
  FAIL_UNLESS(cbuf_store_item(cb, "kumquat"));
  FAIL_UNLESS(cbuf_store_item(cb, "loquat"));

  /* Verify taking */
  FAIL_UNLESS_STREQUAL(cbuf_take_item(cb), "cherry");
  FAIL_UNLESS_STREQUAL(cbuf_take_item(cb), "durian");
  FAIL_UNLESS_STREQUAL(cbuf_take_item(cb), "elderberry");
  FAIL_UNLESS_STREQUAL(cbuf_take_item(cb), "fig");
  FAIL_UNLESS_STREQUAL(cbuf_take_item(cb), "guava");
  FAIL_UNLESS_EQUAL(cbuf_size(cb), 10);
  FAIL_UNLESS_EQUAL(cbuf_count(cb), 5);
  FAIL_UNLESS_EQUAL(cbuf_free_space(cb), 5);
  FAIL_UNLESS_STREQUAL(cbuf_take_item(cb), "huckleberry");
  FAIL_UNLESS_STREQUAL(cbuf_take_item(cb), "imbe");
  FAIL_UNLESS_STREQUAL(cbuf_take_item(cb), "jujube");
  FAIL_UNLESS_STREQUAL(cbuf_take_item(cb), "kumquat");
  FAIL_UNLESS_STREQUAL(cbuf_take_item(cb), "loquat");
  FAIL_UNLESS_STREQUAL(cbuf_take_item(cb), NULL);
  FAIL_UNLESS_EQUAL(cbuf_take_item(cb), NULL);
  FAIL_UNLESS_EQUAL(cbuf_size(cb), 10);
  FAIL_UNLESS_EQUAL(cbuf_count(cb), 0);
  FAIL_UNLESS_EQUAL(cbuf_free_space(cb), 10);

  /* Cleanup. */
  cbuf_destroy(cb, NULL, NULL);
  return 0;
}

int suite_circular_buffer(void)
{
  TEST_PROLOGUE("Circular Buffer Tests");
  RUNTEST(test_circular_buffer_empty);
  RUNTEST(test_circular_buffer_store);
  RUNTEST(test_circular_buffer_full);
  TEST_EPILOGUE;
}
