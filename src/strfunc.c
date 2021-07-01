/******************************************************************************
 *
 * Copyright (C) 2006-2017 by
 * The Salk Institute for Biological Studies and
 * Pittsburgh Supercomputing Center, Carnegie Mellon University
 *
 * Use of this source code is governed by an MIT-style
 * license that can be found in the LICENSE file or at
 * https://opensource.org/licenses/MIT.
 *
******************************************************************************/

/*  Character string handling functions */

#include "config.h"

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>

#include "strfunc.h"

/*************************************************************************
my_strcat:
  In: two strings
  Out: the two strings concatenated, in a newly malloced block of memory,
       or NULL if there isn't enough memory.
  Note: the calling function is responsible for freeing the memory.
*************************************************************************/
char *my_strcat(char const *s1, char const *s2) {
  char *temp = NULL;
  size_t len1, len2;

  len1 = (s1 == NULL) ? 0 : strlen(s1);
  len2 = (s2 == NULL) ? 0 : strlen(s2);
  if ((temp = (char *)malloc(len1 + len2 + 1)) != NULL) {
    if (len1)
      strcpy(temp, s1);
    if (len2)
      strcpy(temp + len1, s2);
    temp[len1 + len2] = '\0';
  }

  return (temp);
}

/*************************************************************************
strip_quotes:
  In: a string that must be at least two characters long
  Out: a copy of the string, newly malloced, missing the first and
       last characters (might even be quotes!); NULL is returned if
       malloc fails.
  Note: this function does NOT do any error checking!
*************************************************************************/
char *strip_quotes(char const *s) {
  char *temp = NULL;
  size_t len = strlen(s);

  if ((temp = (char *)malloc(len - 1)) != NULL) {
    strncpy(temp, s + 1, len - 2);
    temp[len - 2] = '\0';
  }

  return (temp);
}

/*
 * Format a string into an allocated buffer.
 */
char *alloc_vsprintf(char const *fmt, va_list args) {
  char stack_buffer[256];
  int len;
  char *retval = NULL;
  va_list saved_args;
  va_copy(saved_args, args);

  len = vsnprintf(stack_buffer, sizeof(stack_buffer), fmt, args);
  if (len >= (int)sizeof(stack_buffer)) {
    retval = (char *)malloc(len + 1);
    if (retval != NULL)
      vsnprintf(retval, len + 1, fmt, saved_args);
  } else
    retval = strdup(stack_buffer);

  va_end(saved_args);
  return retval;
}

/*
 * Format a string into an allocated buffer.
 */
char *alloc_sprintf(char const *fmt, ...) {
  char *retval;
  va_list args;
  va_start(args, fmt);
  retval = alloc_vsprintf(fmt, args);
  va_end(args);
  return retval;
}
