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

/**************************************************************************\
 ** File: mem_util.c                                                     **
 **                                                                      **
 ** Purpose: Handles various common memory-allocation tasks (stack,etc.) **
 **                                                                      **
 ** Testing status: counter broken.  Everything else validated.          **
 **   (See validate_mem_util.c.)                                         **
\**************************************************************************/

#include "config.h"

#include <inttypes.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#include "strfunc.h"
#include "logging.h"
#include "mem_util.h"

#include "mcell_structs.h"

#ifdef MEM_UTIL_KEEP_STATS
#undef malloc
#undef free
#undef strdup
#undef realloc
#endif

#ifdef DEBUG
#define Malloc count_malloc
#else
#define Malloc malloc
#endif

#ifdef DEBUG
int howmany_count_malloc = 0;

void catch_me() {
  printf("Allocating unreasonably many memory blocks--what are you doing?!\n");
  exit(1);
}

void *count_malloc(size_t n) {
  howmany_count_malloc++;

  if (howmany_count_malloc > 200000) {
    printf("Nuts!\n");
    catch_me();
    return NULL;
  }

  return malloc(n);
}
#endif

/*************************************************************************
 * Checked malloc helpers
 *************************************************************************/

static void memalloc_failure(char const *file, unsigned int line,
                             unsigned int size, char const *desc,
                             int onfailure) {
#ifdef MEM_UTIL_LINE_NUMBERS
  if (desc)
    mcell_error_nodie("File '%s', Line %u: Failed to allocate %u bytes for %s.",
                      file, line, size, desc);
  else
    mcell_error_nodie("File '%s', Line %u: Failed to allocate %u bytes.", file,
                      line, size);
#else
  UNUSED(file);
  UNUSED(line);
  if (desc)
    mcell_error_nodie("Failed to allocate %u bytes for %s.", size, desc);
  else
    mcell_error_nodie("Failed to allocate %u bytes.", size);
#endif

  if (onfailure & CM_EXIT)
    mcell_error("Out of memory.\n"); /* extra newline */
  else
    mcell_error_nodie("Out of memory.\n"); /* extra newline */
}

char *checked_strdup(char const *s, char const *file, unsigned int line,
                     char const *desc, int onfailure) {
  if (s == NULL)
    return NULL;

  char *data = strdup(s);
  if (data == NULL)
    memalloc_failure(file, line, 1 + strlen(s), desc, onfailure);
  return data;
}

void *checked_malloc(size_t size, char const *file, unsigned int line,
                     char const *desc, int onfailure) {
  void *data = malloc(size);
  if (data == NULL)
    memalloc_failure(file, line, size, desc, onfailure);
  return data;
}

char *checked_alloc_sprintf(char const *file, unsigned int line, int onfailure,
                            char const *fmt, ...) {
  va_list args, saved_args;
  va_start(args, fmt);
  va_copy(saved_args, args);

  char *data = alloc_vsprintf(fmt, args);
  if (data == NULL) {
    int needlen = vsnprintf(NULL, 0, fmt, saved_args);
    memalloc_failure(file, line, needlen + 1, "formatted string", onfailure);
  }
  va_end(args);
  va_end(saved_args);
  return data;
}

void *checked_mem_get(struct mem_helper *mh, char const *file,
                      unsigned int line, char const *desc, int onfailure) {
  void *data = mem_get(mh);
  if (data == NULL)
    memalloc_failure(file, line, mh->record_size, desc, onfailure);
  return data;
}


/**************************************************************************\
 ** mem section: reusable block storage for small records of linked      **
 **   lists. Each element is assumed to start with a "next" pointer.     **
\**************************************************************************/

#ifdef MEM_UTIL_KEEP_STATS
struct mem_stats {
  struct mem_stats *next;
  char const *name;
  long long int size;
  long long int max_arenas;
  long long int total_arenas;
  long long int num_arenas_unfreed;
  long long int non_head_arenas;
  long long int max_non_head_arenas;
  long long int total_non_head_arenas;
  long long int unfreed_length;
  long long int max_length;
  long long int cur_free;
  long long int max_free;
  long long int cur_alloc;
  long long int max_alloc;
};

static struct mem_stats *mem_stats_root = NULL;
static long long int mem_stats_cur_mallocs = 0;
static long long int mem_stats_cur_malloc_space = 0;
static long long int mem_stats_max_mallocs = 0;
static long long int mem_stats_max_malloc_space = 0;
static long long mem_stats_total_mallocs = 0;

static long long int mem_cur_overall_allocation = 0;
static long long int mem_max_overall_allocation = 0;
static long long int mem_cur_overall_wastage = 0;
static long long int mem_max_overall_wastage = 0;

void *mem_util_tracking_malloc(size_t size) {
  unsigned char *bl = (unsigned char *)malloc(size + 16);
  *(unsigned int *)bl = size;
  if (bl != NULL) {
    if (++mem_stats_cur_mallocs > mem_stats_max_mallocs)
      mem_stats_max_mallocs = mem_stats_cur_mallocs;
    if ((mem_stats_cur_malloc_space += size) > mem_stats_max_malloc_space)
      mem_stats_max_malloc_space = mem_stats_cur_malloc_space;
    if ((mem_cur_overall_allocation += size) > mem_max_overall_allocation)
      mem_max_overall_allocation = mem_cur_overall_allocation;
    ++mem_stats_total_mallocs;
  }
  return (void *)(bl + 16);
}

void mem_util_tracking_free(void *data) {
  unsigned char *bl = (unsigned char *)data;
  bl -= 16;
  --mem_stats_cur_mallocs;
  mem_stats_cur_malloc_space -= *(unsigned int *)bl;
  free(bl);
}

void *mem_util_tracking_realloc(void *data, size_t size) {
  unsigned char *bl = (unsigned char *)data;
  bl -= 16;
  unsigned int oldsize = *(unsigned int *)bl;
  if (oldsize >= size)
    return data;
  void *block = mem_util_tracking_malloc(size);
  memcpy(block, data, oldsize);
  mem_util_tracking_free(data);
  return block;
}

char *mem_util_tracking_strdup(char const *in) {
  size_t len = strlen(in);
  void *block = mem_util_tracking_malloc(len + 1);
  memcpy(block, in, len + 1);
  return block;
}

static struct mem_stats *create_stats(char const *name, int size) {
  struct mem_stats *s = (struct mem_stats *)malloc(sizeof(struct mem_stats));
  s->next = mem_stats_root;
  s->name = name;
  s->size = size;
  s->total_arenas = 0;
  s->max_arenas = 0;
  s->num_arenas_unfreed = 0;
  s->non_head_arenas = 0;
  s->max_non_head_arenas = 0;
  s->total_non_head_arenas = 0;
  s->unfreed_length = 0;
  s->max_length = 0;
  s->cur_free = 0;
  s->max_free = 0;
  s->cur_alloc = 0;
  s->max_alloc = 0;
  mem_stats_root = s;
  return s;
}

static struct mem_stats *get_stats(char const *name, int size) {
  struct mem_stats *s;
  for (s = mem_stats_root; s != NULL; s = s->next) {
    if (s->size != size)
      continue;

    if (name == NULL) {
      if (s->name == NULL)
        return s;
    } else {
      if (s->name != NULL && !strcmp(name, s->name))
        return s;
    }
  }

  return create_stats(name, size);
}

void mem_dump_stats(FILE *out) {
  struct mem_stats *s;
  for (s = mem_stats_root; s != NULL; s = s->next) {
    char const *name = s->name;
    if (name == NULL)
      name = "(unnamed)";
    fprintf(out, "%s [%lld]\n", s->name, s->size);
    fprintf(out, "---------------------\n");
    fprintf(out, "Num Arenas:          %lld/%lld (%lld)\n",
            s->num_arenas_unfreed, s->max_arenas, s->total_arenas);
    fprintf(out, "Non-head arenas:     %lld/%lld (%lld)\n", s->non_head_arenas,
            s->max_non_head_arenas, s->total_non_head_arenas);
    fprintf(out, "Total Items:         %lld/%lld\n", s->unfreed_length,
            s->max_length);
    fprintf(out, "Free Items:          %lld/%lld\n", s->cur_free, s->max_free);
    fprintf(out, "Alloced Items:       %lld/%lld\n", s->cur_alloc,
            s->max_alloc);
    fprintf(out, "Space usage:         %lld/%lld\n", s->size * s->cur_alloc,
            s->size * s->max_alloc);
    fprintf(out, "Space wastage:       %lld/%lld\n", s->size * s->cur_free,
            s->size * s->max_free);
    fprintf(out, "Arena overhead:      %lld/%lld\n",
            (int)sizeof(struct mem_helper) * s->num_arenas_unfreed,
            (int)sizeof(struct mem_helper) * s->max_arenas);
    fprintf(out, "\n");
  }

  fprintf(out, "Malloc stats:\n");
  fprintf(out, "---------------------\n");
  fprintf(out, "Mallocs:               %lld/%lld (%lld)\n",
          mem_stats_cur_mallocs, mem_stats_max_mallocs,
          mem_stats_total_mallocs);
  fprintf(out, "Space used:            %lld/%lld\n", mem_stats_cur_malloc_space,
          mem_stats_max_malloc_space);
  fprintf(out, "\n");

  fprintf(out, "Summary stats:\n");
  fprintf(out, "---------------------\n");
  fprintf(out, "Allocation:            %lld/%lld\n", mem_cur_overall_allocation,
          mem_max_overall_allocation);
  fprintf(out, "Waste:                 %lld/%lld\n", mem_cur_overall_wastage,
          mem_max_overall_wastage);
  fprintf(out, "\n");
}

#endif

/*************************************************************************
create_mem_named:
   In: Size of a single element (including the leading "next" pointer)
       Number of elements to allocate at once
       Name of "arena" (used for statistics)
   Out: Pointer to a new mem_helper struct.
*************************************************************************/

struct mem_helper *create_mem_named(size_t size, int length, char const *name) {
  struct mem_helper *mh;
  mh = (struct mem_helper *)Malloc(sizeof(struct mem_helper));

  if (mh == NULL)
    return NULL;

  mh->buf_len = (length > 0) ? length : 128;
  mh->record_size =
      (size > (int)sizeof(void *)) ? (size_t)size : sizeof(void *);
  mh->buf_index = 0;
  mh->defunct = NULL;
  mh->next_helper = NULL;

#ifndef MEM_UTIL_NO_POOLING
#ifdef MEM_UTIL_TRACK_FREED
  mh->heap_array =
      (unsigned char *)Malloc(mh->buf_len * (mh->record_size + sizeof(int)));
  memset(mh->heap_array, 0, mh->buf_len * (mh->record_size + sizeof(int)));
#else
  mh->heap_array = (unsigned char *)Malloc(mh->buf_len * mh->record_size);
#endif

  if (mh->heap_array == NULL) {
    free(mh);
    return NULL;
  }
#else
  mh->heap_array = NULL;
#endif

#ifdef MEM_UTIL_KEEP_STATS
  struct mem_stats *s = mh->stats = get_stats(name, size);
  ++s->num_arenas_unfreed;
  ++s->total_arenas;
  if (s->num_arenas_unfreed > s->max_arenas)
    s->max_arenas = s->num_arenas_unfreed;
  s->unfreed_length += length;
  if (s->unfreed_length > s->max_length)
    s->max_length = s->unfreed_length;
  s->cur_free += length;
  if (s->cur_free > s->max_free)
    s->max_free = s->cur_free;
  if ((mem_cur_overall_wastage += size * length) > mem_max_overall_wastage)
    mem_max_overall_wastage = mem_cur_overall_wastage;
#else
  UNUSED(name);
#endif

  return mh;
}

/*************************************************************************
create_mem:
   In: Size of a single element (including the leading "next" pointer)
       Number of elements to allocate at once
   Out: Pointer to a new mem_helper struct.
*************************************************************************/
struct mem_helper *create_mem(size_t size, int length) {
  return create_mem_named(size, length, NULL);
}

/*************************************************************************
mem_get:
   In: A mem_helper
   Out: A pointer to new storage of the appropriate size.
   Note: Use this instead of "Malloc".
*************************************************************************/

void *mem_get(struct mem_helper *mh) {
#ifdef MEM_UTIL_NO_POOLING
  return malloc(mh->record_size);
#else
  if (mh->defunct != NULL) {
    struct abstract_list *retval;
    retval = mh->defunct;
    mh->defunct = retval->next;
#ifdef MEM_UTIL_KEEP_STATS
    struct mem_stats *s = mh->stats;
    --s->cur_free;
    ++s->cur_alloc;
    if (s->cur_alloc > s->max_alloc)
      s->max_alloc = s->cur_alloc;
    if ((mem_cur_overall_allocation += mh->record_size) >
        mem_max_overall_allocation)
      mem_max_overall_allocation = mem_cur_overall_allocation;
    mem_cur_overall_wastage -= mh->record_size;
#endif
#ifdef MEM_UTIL_ZERO_FREED
    {
      unsigned char *thisData = (unsigned char *)retval;
      size_t i;
      thisData += sizeof(struct abstract_element *);
      for (i = 0; i < mh->record_size - sizeof(struct abstract_element *);) {
        if (thisData[i] != '\0') {
          mcell_warn("Memory block at %08lx: non-zero at byte %d: %02x %02x "
                     "%02x %02x...",
                     retval, i, thisData[i], thisData[i + 1], thisData[i + 2],
                     thisData[i + 3]);
          i += 4;
        } else
          ++i;
      }
      memset(thisData, 0, mh->record_size - sizeof(struct abstract_element *));
    }
#endif
#ifdef MEM_UTIL_TRACK_FREED
    if (((int *)retval)[-1]) {
      mcell_warn("Duplicate allocation of ptr '%p'.", retval);
      return NULL;
    }
    ((int *)retval)[-1] = 1;
#endif
    return (void *)retval;
  } else if (mh->buf_index < mh->buf_len) {
#ifdef MEM_UTIL_TRACK_FREED
    size_t offset = mh->buf_index * (mh->record_size + sizeof(int));
#else
    size_t offset = mh->buf_index * mh->record_size;
#endif
    mh->buf_index++;
#ifdef MEM_UTIL_KEEP_STATS
    struct mem_stats *s = mh->stats;
    --s->cur_free;
    ++s->cur_alloc;
    if (s->cur_alloc > s->max_alloc)
      s->max_alloc = s->cur_alloc;
    if ((mem_cur_overall_allocation += mh->record_size) >
        mem_max_overall_allocation)
      mem_max_overall_allocation = mem_cur_overall_allocation;
    mem_cur_overall_wastage -= mh->record_size;
#endif
#ifdef MEM_UTIL_TRACK_FREED
    int *ptr = (int *)(mh->heap_array + offset);
    if (*ptr) {
      mcell_warn("Duplicate allocation of ptr '%p'.", ptr + 1);
      return NULL;
    }
    *ptr++ = 1;
    return (void *)ptr;
#else
    return (void *)(mh->heap_array + offset);
#endif
  } else {
    struct mem_helper *mhnext;
    unsigned char *temp;
#ifdef MEM_UTIL_KEEP_STATS
    struct mem_stats *s = mh->stats;
    mhnext = create_mem_named(mh->record_size, mh->buf_len, s->name);
    ++s->non_head_arenas;
    if (s->non_head_arenas > s->max_non_head_arenas)
      s->max_non_head_arenas = s->non_head_arenas;
    ++s->total_non_head_arenas;
#else
    mhnext = create_mem(mh->record_size, mh->buf_len);
#endif
    if (mhnext == NULL)
      return NULL;

    /* Swap contents of this mem_helper with new one */
    /* Keeps mh at top of list but with freshly allocated space */
    mhnext->next_helper = mh->next_helper;
    temp = mhnext->heap_array;
    mhnext->heap_array = mh->heap_array;
    mh->heap_array = temp;
    mhnext->buf_index = mh->buf_index;
    mh->next_helper = mhnext;

    mh->buf_index = 0;
    return mem_get(mh);
  }
#endif
}

/*************************************************************************
mem_put:
   In: A mem_helper
       A pointer to the memory you no longer need
   Out: No return value. The defunct memory is stored for later use.
   Note: Defunct memory need not have come from the same mem_helper.
         Use this instead of "free".
*************************************************************************/

void mem_put(struct mem_helper *mh, void *defunct) {
#ifdef MEM_UTIL_NO_POOLING
  free(defunct);
  return;
#else
  struct abstract_list *data = (struct abstract_list *)defunct;
#ifdef MEM_UTIL_TRACK_FREED
  int *ptr = (int *)data;
  if (ptr[-1] == 0) {
    mcell_warn("Duplicate free of ptr '%p'.", defunct);
    return;
  }
  ptr[-1] = 0;
#endif
#ifdef MEM_UTIL_ZERO_FREED
  memset(data, 0, mh->record_size);
#endif
  data->next = mh->defunct;
  mh->defunct = data;
#ifdef MEM_UTIL_KEEP_STATS
  struct mem_stats *s = mh->stats;
  ++s->cur_free;
  --s->cur_alloc;
  if (s->cur_free > s->max_free)
    s->max_free = s->cur_free;
  mem_cur_overall_allocation -= mh->record_size;
  if ((mem_cur_overall_wastage += mh->record_size) > mem_max_overall_wastage)
    mem_max_overall_wastage = mem_cur_overall_wastage;
#endif
#endif
}

/*************************************************************************
mem_put_list:
   In: A mem_helper
       A pointer to the memory you no longer need
   Out: No return value.  The input is taken to be a NULL-terminated
        singly-linked list.  Everything in that list is "freed" and
        held for later use.
*************************************************************************/

void mem_put_list(struct mem_helper *mh, void *defunct) {
  struct abstract_list *data = (struct abstract_list *)defunct;
  struct abstract_list *alp;

#ifdef MEM_UTIL_NO_POOLING
  struct abstract_list *alpNext;
  for (alp = data; alp != NULL; alp = alpNext) {
    alpNext = alp->next;
    free(alp);
  }
#else
#ifdef MEM_UTIL_ZERO_FREED
  for (alp = data; alp != NULL; alp = alp->next) {
    unsigned char *thisData = (unsigned char *)alp;
    thisData += sizeof(struct abstract_element *);
    memset(thisData, 0, mh->record_size - sizeof(struct abstract_element *));
  }
#endif
#ifdef MEM_UTIL_TRACK_FREED
  for (alp = data; alp != NULL; alp = alp->next) {
    int *ptr = (int *)alp;
    if (ptr[-1] == 0) {
      mcell_warn("Duplicate free of ptr '%p'.", defunct);
      return;
    }
    ptr[-1] = 0;
  }
#endif
#ifdef MEM_UTIL_KEEP_STATS
  int count = 1;
  for (alp = data; alp->next != NULL; alp = alp->next)
    ++count;
  struct mem_stats *s = mh->stats;
  s->cur_free += count;
  s->cur_alloc -= count;
  if (s->cur_free > s->max_free)
    s->max_free = s->cur_free;
  mem_cur_overall_allocation -= mh->record_size * count;
  if ((mem_cur_overall_wastage += mh->record_size * count) >
      mem_max_overall_wastage)
    mem_max_overall_wastage = mem_cur_overall_wastage;
#else
  for (alp = data; alp->next != NULL; alp = alp->next) {
  }
#endif

  alp->next = mh->defunct;
  mh->defunct = data;
#endif
}

/*************************************************************************
delete_mem:
   In: A mem_helper
   Out: All memory allocated by this mem_helper is freed, and the
        mem_helper struct itself is also.  However, defunct memory
        from other mem_helpers (or elsewhere) is not freed.
*************************************************************************/

void delete_mem(struct mem_helper *mh) {
  if (mh == NULL)
    return;
#ifndef MEM_UTIL_NO_POOLING
#ifdef MEM_UTIL_KEEP_STATS
  struct mem_stats *s = mh->stats;
  --s->num_arenas_unfreed;
  s->cur_alloc -= mh->buf_index;
  s->cur_free -= (mh->buf_len - mh->buf_index);
  s->unfreed_length -= mh->buf_len;
  mem_cur_overall_allocation -= mh->record_size * mh->buf_len;
  mem_cur_overall_wastage -= mh->record_size * (mh->buf_len - mh->buf_index);
  if (mh->next_helper) {
    delete_mem(mh->next_helper);
    --s->max_non_head_arenas;
  }
#else
  if (mh->next_helper)
    delete_mem(mh->next_helper);
#endif
  free(mh->heap_array);
#endif
  free(mh);
}
