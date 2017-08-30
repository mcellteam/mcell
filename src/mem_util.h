/******************************************************************************
 *
 * Copyright (C) 2006-2017 by
 * The Salk Institute for Biological Studies and
 * Pittsburgh Supercomputing Center, Carnegie Mellon University
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301,
 * USA.
 *
******************************************************************************/

#pragma once

#ifdef MEM_UTIL_KEEP_STATS
#include <stdio.h>
#include <stdlib.h>
char *mem_util_tracking_strdup(char const *in);
void *mem_util_tracking_malloc(unsigned int size);
void *mem_util_tracking_realloc(void *data, unsigned int size);
void mem_util_tracking_free(void *data);
#undef malloc
#undef free
#undef strdup
#define malloc mem_util_tracking_malloc
#define free mem_util_tracking_free
#define strdup mem_util_tracking_strdup
#define realloc mem_util_tracking_realloc
#endif

char *checked_strdup(char const *s, char const *file, unsigned int line,
                     char const *desc, int onfailure);
void *checked_malloc(size_t size, char const *file, unsigned int line,
                     char const *desc, int onfailure);
char *checked_alloc_sprintf(char const *file, unsigned int line, int onfailure,
                            char const *fmt, ...) PRINTF_FORMAT(4);

struct mem_helper;
void *checked_mem_get(struct mem_helper *mh, char const *file,
                      unsigned int line, char const *desc, int onfailure);

#define CM_EXIT (1)
#define CHECKED_STRDUP_NODIE(s, desc)                                          \
  checked_strdup((s), __FILE__, __LINE__, desc, 0)
#define CHECKED_STRDUP(s, desc)                                                \
  checked_strdup((s), __FILE__, __LINE__, desc, CM_EXIT)
#define CHECKED_MALLOC_NODIE(sz, desc)                                         \
  checked_malloc((sz), __FILE__, __LINE__, desc, 0)
#define CHECKED_MALLOC(sz, desc)                                               \
  checked_malloc((sz), __FILE__, __LINE__, desc, CM_EXIT)
#define CHECKED_MALLOC_STRUCT_NODIE(tp, desc)                                  \
  (tp *) checked_malloc(sizeof(tp), __FILE__, __LINE__, desc, 0)
#define CHECKED_MALLOC_STRUCT(tp, desc)                                        \
  (tp *) checked_malloc(sizeof(tp), __FILE__, __LINE__, desc, CM_EXIT)
#define CHECKED_MALLOC_ARRAY_NODIE(tp, num, desc)                              \
  (tp *) checked_malloc((num) * sizeof(tp), __FILE__, __LINE__, desc, 0)
#define CHECKED_MALLOC_ARRAY(tp, num, desc)                                    \
  (tp *) checked_malloc((num) * sizeof(tp), __FILE__, __LINE__, desc, CM_EXIT)
#define CHECKED_MEM_GET_NODIE(mh, desc)                                        \
  checked_mem_get((mh), __FILE__, __LINE__, desc, 0)
#define CHECKED_MEM_GET(mh, desc)                                              \
  checked_mem_get((mh), __FILE__, __LINE__, desc, CM_EXIT)
#define CHECKED_SPRINTF_NODIE(fmt, ...)                                        \
  checked_alloc_sprintf(__FILE__, __LINE__, 0, fmt, ##__VA_ARGS__)
#define CHECKED_SPRINTF(fmt, ...)                                              \
  checked_alloc_sprintf(__FILE__, __LINE__, CM_EXIT, fmt, ##__VA_ARGS__)

/* Everything allocated by a mem_helper must start with a next pointer
   as if it were derived from an abstract_list */
struct abstract_list {
  struct abstract_list *next;
};

/* Data structure to allocate blocks of memory for a specific size of struct */
struct mem_helper {
  int buf_len;               /* Number of elements to allocate at once  */
  int buf_index;             /* Index of the next unused element in the array */
  size_t record_size;           /* Size of the element to allocate */
  unsigned char *heap_array; /* Block of memory for elements */
  struct abstract_list *defunct; /* Linked list of elements that may be reused
                                    for next memory request */
  struct mem_helper *next_helper; /* Next (fully-used) mem_helper */
#ifdef MEM_UTIL_KEEP_STATS
  struct mem_stats *stats;
#endif
};

#ifdef MEM_UTIL_KEEP_STATS
void mem_dump_stats(FILE *out);
#else
#define mem_dump_stats(out)                                                    \
  do { /* do nothing */                                                        \
  } while (0)
#endif

struct mem_helper *create_mem_named(size_t size, int length, char const *name);
struct mem_helper *create_mem(size_t size, int length);
void *mem_get(struct mem_helper *mh);
void mem_put(struct mem_helper *mh, void *defunct);
void mem_put_list(struct mem_helper *mh, void *defunct);
void delete_mem(struct mem_helper *mh);

#define stack_nonempty(sh) ((sh)->index > 0 || (sh)->next != NULL)
