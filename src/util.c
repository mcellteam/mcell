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

#include "config.h"

#include <assert.h>
#include <float.h>
#include <math.h>
#include <inttypes.h>
#include <memory.h>
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <sys/stat.h>
#include <unistd.h>
#include <stdbool.h>

#include "logging.h"
#include "util.h"
#include "mcell_structs.h"

/*******************************************************************
new_bit_array: mallocs an array of the desired number of bits

 In:
    bits: how many bits to place in the array

 Out:
    A pointer to a newly allocated bit_array struct, or NULL
    on memory error.
*******************************************************************/
struct bit_array *new_bit_array(int bits) {
  int n = (bits + 8 * sizeof(int) - 1) / (8 * sizeof(int));
  /* Allocate contiguous memory for struct bit_array and its associated bits */
  struct bit_array *ba = (struct bit_array *)malloc(sizeof(struct bit_array) +
    sizeof(int) * n);
  if (ba == NULL) {
    return NULL;
  }

  ba->nbits = bits;
  ba->nints = n;
  return ba;
}

/*************************************************************************
duplicate_bit_array: mallocs an array and copies an existing bit array

 In:
    old: existing bit array to duplicate (in newly malloced memory)

 Out:
    A pointer to a newly allocated bit_array struct, or NULL
    on memory error.
*************************************************************************/
struct bit_array *duplicate_bit_array(struct bit_array *old) {
  struct bit_array *ba = (struct bit_array *)malloc(sizeof(struct bit_array) +
    sizeof(int) * old->nints);
  if (ba == NULL) {
    return NULL;
  }

  memcpy(ba, old, sizeof(struct bit_array) + sizeof(int) * old->nints);
  return ba;
}

/*******************************************************************
get_bit: returns the value of a bit in a bit_array

 In:
    ba: pointer to a bit_array struct
    idx: the index of the bit to return

 Out:
    0 or 1, depending on whether the idx'th bit is set.  No
    bounds checking is performed.
*******************************************************************/
int get_bit(struct bit_array *ba, int idx) {
  int *data = &(ba->nints);
  data++; /* At start of bit array memory */

  size_t ofs = idx & (8 * sizeof(int) - 1);
  idx = idx / (8 * sizeof(int));
  ofs = 1u << ofs;

  if ((data[idx] & ofs) != 0) {
    return 1;
  } else {
    return 0;
  }
}

/*******************************************************************
set_bit: set a value in a bit array

 In:
    ba: pointer to a bit_array struct
    idx: the index of the bit to set
    value: 0 = turn bit off; nonzero = turn bit on

 Out:
    Nothing
*******************************************************************/
void set_bit(struct bit_array *ba, int idx, int value) {
  int *data = &(ba->nints);
  data++; /* At start of bit array memory */

  size_t ofs = idx & (8 * sizeof(int) - 1);
  idx = idx / (8 * sizeof(int));
  ofs = (1u << ofs);

  if (value) {
    value = ofs;
  } else {
    value = 0;
  }

  data[idx] = (data[idx] & ~ofs) | value;
}

/*******************************************************************
set_bit_range: set many bits to a value in a bit array

 In:
    ba: pointer to a bit_array struct
    idx1: the index of the first bit to set
    idx2: the index of the last bit to set
    value: 0 = turn bits off; nonzero = turn bits on

 Out:
    Nothing
*******************************************************************/
void set_bit_range(struct bit_array *ba, int idx1, int idx2, int value) {
  int *data = &(ba->nints);
  data++; /* At start of bit array memory */

  int ofs1 = idx1 & (8 * sizeof(int) - 1);
  int ofs2 = idx2 & (8 * sizeof(int) - 1);
  idx1 = idx1 / (8 * sizeof(int));
  idx2 = idx2 / (8 * sizeof(int));

  unsigned int mask, cmask;
  if (idx1 == idx2) {
    mask = 0;
    for (int i = ofs1; i <= ofs2; i++) {
      mask |= (1u << i);
    }
    cmask = ~mask;

    if (value) {
      data[idx1] = (data[idx1] & cmask) | mask;
    } else {
      data[idx1] = data[idx1] & cmask;
    }
  } else {
    if (value) {
      value = ~0;
    } else {
      value = 0;
    }
    for (int i = idx1 + 1; i < idx2; i++) {
      data[i] = value;
    }

    mask = 0;
    for (unsigned int i = ofs1; i < 8 * sizeof(int); i++) {
      mask |= (1u << i);
    }
    cmask = ~mask;
    if (value) {
      data[idx1] = (data[idx1] & cmask) | mask;
    } else {
      data[idx1] = data[idx1] & cmask;
    }

    mask = 0;
    for (int i = 0; i <= ofs2; i++) {
      mask |= (1u << i);
    }
    cmask = ~mask;
    if (value) {
      data[idx2] = (data[idx2] & cmask) | mask;
    } else {
      data[idx2] = data[idx2] & cmask;
    }
  }
}

/*******************************************************************
set_all_bits: sets all values in a bit array

 In:
    ba: pointer to a bit_array struct
    value: 0 = turn bits off; nonzero = turn bits on

 Out:
    Nothing
*******************************************************************/
void set_all_bits(struct bit_array *ba, int value) {
  if (value) {
    value = ~0;
  }

  int *data = &(ba->nints);
  data++; /* At start of bit array memory */

  for (int i = 0; i < ba->nints; i++) {
    data[i] = value;
  }
}

/*******************************************************************
bit operation: performs a logical operation on two bit arrays

 In:
    ba: pointer to a bit_array struct
    bb: pointer to another bit_array struct
    op: character that determines which operation to perform
        '!' -- ba = NOT ba
        '~' -- same
        '|' -- ba = ba OR bb
        '+' -- same
        '&' -- ba = ba AND bb
        '^' -- ba = ba XOR bb
        '-' -- ba = ba AND NOT bb

 Out:
    Nothing
*******************************************************************/
void bit_operation(struct bit_array *ba, struct bit_array *bb, char op) {
  int *da, *db;
  if (op == '!' || op == '~') {
    da = &(ba->nints);
    da++;
    for (int i = 0; i < ba->nints; i++) {
      da[i] = ~da[i];
    }
    return;
  }

  if (ba->nbits != bb->nbits) {
    return;
  }

  da = &(ba->nints);
  da++;
  db = &(bb->nints);
  db++;

  switch (op) {
  case '^':
    for (int i = 0; i < ba->nints; i++) {
      da[i] ^= db[i];
    }
    break;
  case '|':
  case '+':
    for (int i = 0; i < ba->nints; i++) {
      da[i] |= db[i];
    }
    break;
  case '-':
    for (int i = 0; i < ba->nints; i++) {
      da[i] &= ~db[i];
    }
    break;
  case '&':
    for (int i = 0; i < ba->nints; i++) {
      da[i] &= db[i];
    }
    break;
  default:
    break;
  }
}

/**********************************************************************
count_bits: count how many bits are set in a bit array

 In:
    ba: pointer to a bit_array struct

 Out:
    int containing number of nonzero bits
**********************************************************************/
int count_bits(struct bit_array *ba) {
  static const int cb_table[256] = {
    0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4, 1, 2, 2, 3, 2, 3, 3, 4,
    2, 3, 3, 4, 3, 4, 4, 5, 1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
    2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 1, 2, 2, 3, 2, 3, 3, 4,
    2, 3, 3, 4, 3, 4, 4, 5, 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
    2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 3, 4, 4, 5, 4, 5, 5, 6,
    4, 5, 5, 6, 5, 6, 6, 7, 1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
    2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 2, 3, 3, 4, 3, 4, 4, 5,
    3, 4, 4, 5, 4, 5, 5, 6, 3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
    2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 3, 4, 4, 5, 4, 5, 5, 6,
    4, 5, 5, 6, 5, 6, 6, 7, 3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
    4, 5, 5, 6, 5, 6, 6, 7, 5, 6, 6, 7, 6, 7, 7, 8
  };

  int *dd = &(ba->nints);
  dd++;
  unsigned char *d = (unsigned char *)dd;

  int n = (ba->nints - 1) * sizeof(int);
  int cnt = 0;
  for (int i = 0; i < n; i++) {
    cnt += cb_table[(*d++)];
  }

  n = ba->nbits - n * 8;
  if (n == 0)
    return cnt;

  int j = dd[ba->nints - 1];
  while (n >= 8) {
    cnt += cb_table[j & 0xFF];
    n -= 8;
    j >>= 8;
  }
  if (n > 0) {
    cnt += cb_table[j & 0xFF] - cb_table[(j & 0xFF) >> n];
  }
  return cnt;
}

/**********************************************************************
free_bit_array: frees a bit array (just a wrapper to free() for now)

 In:
    ba: pointer to a bit_array struct

 Out:
    Nothing
**********************************************************************/
void free_bit_array(struct bit_array *ba) { free(ba); }

/*************************************************************************
bisect:
  In: array of doubles, sorted low to high
      int saying how many doubles there are
      double we are using to bisect the array
  Out: index of the largest element in the array smaller than the bisector
*************************************************************************/
int bisect(double *list, int n, double val) {
  int lo = 0;
  int hi = n;
  int mid = 0;
  while (hi - lo > 1) {
    mid = (hi + lo) / 2;
    if (list[mid] > val)
    {
      hi = mid;
    } else {
      lo = mid;
    }
  }
  return lo;
}

/*************************************************************************
bisect_near:
  In: array of doubles, sorted low to high
      int saying how many doubles there are
      double we are using to bisect the array
  Out: index of the element closest to val
*************************************************************************/
int bisect_near(double *list, int n, double val) {
  int lo = 0;
  int hi = n - 1;
  int mid = 0;
  while (hi - lo > 1) {
    mid = (hi + lo) / 2;
    if (list[mid] > val)
    {
      hi = mid;
    } else {
      lo = mid;
    }
  }
  if (val > list[hi]) {
    return hi;
  } else if (val < list[lo]) {
    return lo;
  } else if (val - list[lo] < list[hi] - val) {
    return lo;
  } else {
    return hi;
  }
}

/*************************************************************************
bisect_high:
  In: array of doubles, sorted low to high
      int saying how many doubles there are
      double we are using to bisect the array
  Out: index of the smallest element in the array larger than the bisector
*************************************************************************/
int bisect_high(double *list, int n, double val) {
  int lo = 0;
  int hi = n - 1;
  int mid = 0;
  while (hi - lo > 1) {
    mid = (hi + lo) / 2;
    if (list[mid] > val) {
      hi = mid;
    } else {
      lo = mid;
    }
  }
  if (list[lo] > val)
  {
    return lo;
  } else {
    return hi;
  }
}

/**********************************************************************
distinguishable: reports whether two doubles are measurably different

 In:
    a: first double
    b: second double
    eps: fractional difference that we think is different

 Out:
    1 if the numbers are different, 0 otherwise
**********************************************************************/
inline int distinguishable(double a, double b, double eps) {
  double c = fabs(a - b);
  a = fabs(a);
  if (a < 1) {
    a = 1;
  }
  b = fabs(b);

  if (b < a) {
    eps *= a;
  } else {
    eps *= b;
  }
  return (c > eps);
}

/**********************************************************************
is_reverse_abbrev: reports whether the first string is a reverse
  abbreviation of the second, i.e. whether it matches the end of
  the second string.
**********************************************************************/
int is_reverse_abbrev(char *abbrev, char *full) {
  size_t na = strlen(abbrev);
  size_t nf = strlen(full);
  if (na > nf) {
    return 0;
  }
  return (strcmp(abbrev, full + (nf - na)) == 0);
}

/*************************************************************************
void_list_sort:
  In: linked list contining void pointers
  Out: linked list is mergesorted by memory address
*************************************************************************/
struct void_list *void_list_sort(struct void_list *vl) {
  struct void_list *stack[64];
  int stack_n[64];
  struct void_list *left, *right, *merge, *tail;
  int si = 0;

  /* HACK: If vl == NULL, we return stack[0] unmodified, so initialize it to
   *       NULL. */
  stack[0] = NULL;

  while (vl != NULL) {
    if (vl->next == NULL) {
      stack[si] = vl;
      stack_n[si] = 1;
      vl = NULL;
      si++;
    } else if ((intptr_t)vl->data <= (intptr_t)vl->next->data) {
      stack[si] = vl;
      stack_n[si] = 2;
      vl = vl->next->next;
      stack[si]->next->next = NULL;
      si++;
    } else {
      stack[si] = vl->next;
      stack_n[si] = 2;
      left = vl;
      vl = vl->next->next;
      stack[si]->next = left;
      left->next = NULL;
      si++;
    }
    while (si > 1 && stack_n[si - 1] * 2 >= stack_n[si - 2]) {
      stack_n[si - 2] += stack_n[si - 1];

      left = stack[si - 2];
      right = stack[si - 1];
      if ((intptr_t)left->data <= (intptr_t)right->data) {
        merge = left;
        left = left->next;
      } else {
        merge = right;
        right = right->next;
      }
      merge->next = NULL;
      tail = merge;

      while (1) {
        if (left == NULL) {
          tail->next = right;
          break;
        }
        if (right == NULL) {
          tail->next = left;
          break;
        }

        if ((intptr_t)left->data <= (intptr_t)right->data) {
          tail->next = left;
          tail = left;
          left = left->next;
        } else {
          tail->next = right;
          tail = right;
          right = right->next;
        }
      }

      stack[si - 2] = merge;
      si--;
    }
  }

  while (si > 1) /* Exact duplicate of code in loop--keep it this way! */
  {
    stack_n[si - 2] += stack_n[si - 1];

    left = stack[si - 2];
    right = stack[si - 1];
    if ((intptr_t)left->data <= (intptr_t)right->data) {
      merge = left;
      left = left->next;
    } else {
      merge = right;
      right = right->next;
    }
    merge->next = NULL;
    tail = merge;

    while (1) {
      if (left == NULL) {
        tail->next = right;
        break;
      }
      if (right == NULL) {
        tail->next = left;
        break;
      }

      if ((intptr_t)left->data <= (intptr_t)right->data) {
        tail->next = left;
        tail = left;
        left = left->next;
      } else {
        tail->next = right;
        tail = right;
        right = right->next;
      }
    }

    stack[si - 2] = merge;
    si--;
  }

  return stack[0];
}

/*************************************************************************
void_list_sort_by:
  In: linked list containing void pointers
      comparison function that compares two void pointers
  Out: linked list is mergesorted according to function
  Note: function should implement "less than or equal to", i.e., it
        should return a nonzero value if the first pointer is considered
        to be less than or equal to the second (based on contents or
        whatever), and it should return zero otherwise.
*************************************************************************/
struct void_list *void_list_sort_by(struct void_list *vl,
                                    int (*leq)(void *, void *)) {
  struct void_list *stack[64];
  int stack_n[64];
  struct void_list *left, *right, *merge, *tail;
  int si = 0;

  /* HACK: If vl == NULL, we return stack[0] unmodified, so initialize it to
   *       NULL. */
  stack[0] = NULL;

  while (vl != NULL) {
    if (vl->next == NULL) {
      stack[si] = vl;
      stack_n[si] = 1;
      vl = NULL;
      si++;
    } else if ((*leq)(vl->data, vl->next->data)) {
      stack[si] = vl;
      stack_n[si] = 2;
      vl = vl->next->next;
      stack[si]->next->next = NULL;
      si++;
    } else {
      stack[si] = vl->next;
      stack_n[si] = 2;
      left = vl;
      vl = vl->next->next;
      stack[si]->next = left;
      left->next = NULL;
      si++;
    }
    while (si > 1 && stack_n[si - 1] * 2 >= stack_n[si - 2]) {
      stack_n[si - 2] += stack_n[si - 1];

      left = stack[si - 2];
      right = stack[si - 1];
      if ((*leq)(left->data, right->data)) {
        merge = left;
        left = left->next;
      } else {
        merge = right;
        right = right->next;
      }
      merge->next = NULL;
      tail = merge;

      while (1) {
        if (left == NULL) {
          tail->next = right;
          break;
        }
        if (right == NULL) {
          tail->next = left;
          break;
        }

        if ((*leq)(left->data, right->data)) {
          tail->next = left;
          tail = left;
          left = left->next;
        } else {
          tail->next = right;
          tail = right;
          right = right->next;
        }
      }

      stack[si - 2] = merge;
      si--;
    }
  }

  while (si > 1) /* Exact duplicate of code in loop--keep it this way! */
  {
    stack_n[si - 2] += stack_n[si - 1];

    left = stack[si - 2];
    right = stack[si - 1];
    if ((*leq)(left->data, right->data)) {
      merge = left;
      left = left->next;
    } else {
      merge = right;
      right = right->next;
    }
    merge->next = NULL;
    tail = merge;

    while (1) {
      if (left == NULL) {
        tail->next = right;
        break;
      }
      if (right == NULL) {
        tail->next = left;
        break;
      }

      if ((*leq)(left->data, right->data)) {
        tail->next = left;
        tail = left;
        left = left->next;
      } else {
        tail->next = right;
        tail = right;
        right = right->next;
      }
    }

    stack[si - 2] = merge;
    si--;
  }

  return stack[0];
}

/*************************************************************************
void_array_search:
  In: array of void pointers sorted by memory address
      length of the array
      void pointer we're trying to find
  Out: index of the void pointer in the array, or -1 if there is no
       matching pointer in the list
*************************************************************************/
int void_array_search(void **array, int n, void *to_find) {
  int lo = 0;
  int hi = n - 1;
  int m;
  while (hi - lo > 1) {
    m = (hi + lo) / 2;
    if (to_find == array[m]) {
      return m;
    } else if ((intptr_t)to_find > (intptr_t)array[m]) {
      lo = m;
    } else {
      hi = m;
    }
  }

  if (to_find == array[lo]) {
    return lo;
  }
  if (to_find == array[hi]) {
    return hi;
  }
  return -1;
}

/*************************************************************************
void_ptr_compare:
    Utility function to allow sorting an array of pointers by address.
    Conventions are appropriate for use with qsort.

  In: void const *v1 - first pointer
      void const *v2 - second pointer
  Out: -1, 0, or 1 as *(void **)v1 <, =, or > *(void **)v2 resp.
*************************************************************************/
int void_ptr_compare(void const *v1, void const *v2) {
  void const **v1p = (void const **)v1;
  void const **v2p = (void const **)v2;
  intptr_t i1 = (intptr_t) * v1p;
  intptr_t i2 = (intptr_t) * v2p;
  if (i1 < i2) {
    return -1;
  } else if (i1 > i2) {
    return 1;
  }
  return 0;
}

/*********************************************************************
allocate_uint_array:
   In: int size - length of the array to allocate
       u_int value - value with which to initialize elements
   Out: the newly allocated array, with all elements initialized to 'value'
***********************************************************************/
u_int *allocate_uint_array(int size, u_int value) {
  u_int *arr;
  if ((arr = CHECKED_MALLOC_ARRAY_NODIE(u_int, size, NULL)) == NULL) {
    return NULL;
  }
  for (int i = 0; i < size; ++i) {
    arr[i] = value;
  }

  return arr;
}

/*********************************************************************
allocate_ptr_array:
    Allocate an array of pointers.  Use free_ptr_array to free if you want the
    pointers in the array to be freed as well.

        In: int size - length of the array to allocate
        Out: the newly allocated array, with all elements initialized to NULL.
***********************************************************************/
void **allocate_ptr_array(int size) {
  if (size == 0) {
    size = 1;
  }

  void **arr;
  if ((arr = CHECKED_MALLOC_ARRAY_NODIE(void *, size, NULL)) == NULL) {
    return NULL;
  }

  memset(arr, 0, size * sizeof(void *));
  return arr;
}

/*************************************************************************
free_ptr_array:
    Free an array of pointers, freeing any non-NULL pointers within the array.

        In:  void **pa - pointer array to free
             int count - length of pointer array
        Out: All non-NULL pointers in the array are freed, as is the array
             itself.
**************************************************************************/
void free_ptr_array(void **pa, int count) {
  for (int i = 0; i < count; ++i) {
    if (pa[i] != NULL) {
      free(pa[i]);
    }
  }
  free(pa);
}

/*************************************************************************
free_num_expr_list:
    Free a num_expr_list.

        In:  struct num_expr_list *nlist - the list to free
        Out: All elements in the list are freed.
**************************************************************************/
void free_num_expr_list(struct num_expr_list *nlist) {
  struct num_expr_list *nnext;
  while (nlist != NULL) {
    nnext = nlist->next;
    free(nlist);
    nlist = nnext;
  }
}

/*************************************************************************
dir_exists:
    Utility to check if a given directory exists.

        In:  char const *path - absolute or relative path of dir
        Out: 1 if it's a directory, 0 if not
**************************************************************************/
int dir_exists(char const *path) {
  struct stat sb;
  if (stat(path, &sb) == 0 && S_ISDIR(sb.st_mode)) {
    return 1;
  }
  return 0;
}

/*************************************************************************
is_writable_dir:
    Utility to check if a given directory exists and is readable/writable.

        In:  char const *path - absolute or relative path of dir
        Out: 1 if writable, 0 if not
**************************************************************************/
int is_writable_dir(char const *path) {
  if (dir_exists(path) && !access(path, R_OK | W_OK | X_OK)) {
    return 1;
  }
  return 0;
}

/*************************************************************************
make_parent_dir:
    Utility to make the (possibly nested) parent directory of a file.  Will
    attempt to create a directory with full rwx permission.  If the directory
    already exists and has rwx permission for the user, this function will
    return success.  Essentially, this works like mkdirs, but it strips off the
    last path element first.

        In:  char const *path - absolute or relative path of file
        Out: 0 on success, 1 on failure
**************************************************************************/
int make_parent_dir(char const *path) {
  char *pathtmp = CHECKED_STRDUP(path, "directory path");
  char *last_slash = strrchr(pathtmp, '/');
  if (last_slash) {
    *last_slash = '\0';
    if (mkdirs(pathtmp)) {
      free(pathtmp);
      return 1;
    }
  }

  free(pathtmp);
  return 0;
}

/*************************************************************************
mkdirs:
    Utility to make a (possibly nested) directory.  Will attempt to create a
    directory with full rwx permission.  If the directory already exists and
    has rwx permission for the user, this function will return success.

        In:  char const *path - absolute or relative path for dir
        Out: 0 on success, 1 on failure
**************************************************************************/
int mkdirs(char const *path) {
  char *pathtmp = CHECKED_STRDUP(path, "directory path");
  char *curpos = pathtmp;

  /* we need to skip leading '/' in case we have absolute paths */
  while (curpos != NULL && *curpos == '/') {
    ++curpos;
  }

  while (curpos != NULL) {
    /* Find next '/', turn it into '\0' */
    char *nextel = strchr(curpos, '/');
    if (nextel != NULL) {
      *nextel = '\0';
    }

    /* if this directory exists */
    if (dir_exists(pathtmp)) {
      /* Turn '\0' back to '/' */
      if (nextel) {
        *nextel = '/';
        curpos = nextel + 1;
      } else {
        curpos = NULL;
      }
      continue;
    }

    /* Make the directory */
    if (!is_writable_dir(pathtmp) && mkdir(pathtmp, 0777) != 0) {
      mcell_perror_nodie(errno, "Failed to create directory '%s'", path);
      free(pathtmp);
      return 1;
    }

    /* Turn '\0' back to '/' */
    if (nextel) {
      *nextel = '/';
      curpos = nextel + 1;
    } else {
      curpos = NULL;
    }
  }
  free(pathtmp);

  if (access(path, R_OK | W_OK | X_OK) != 0) {
    return 1;
  }
  return 0;
}

/*************************************************************************
open_file:
    Utility to open a file, printing a sensible error message if opening
    fails.
        In: char const *fname - filename for new file
            char const *mode - mode for file access
        Out: file handle for file, NULL on error
**************************************************************************/
FILE *open_file(const char *fname, const char *mode) {
  FILE *f;
  if ((f = fopen(fname, mode)) == NULL) {
    mcell_perror_nodie(errno, "Failed to open file %s.", fname);
    return NULL;
  }

  return f;
}

/*************************************************************************
erfcinv:

  Fast rational function approximation to inverse of the error function,
  based upon algorithm for inverse of Normal cumulative distribution
  function at http://home.online.no/~pjacklam/notes/invnorm/index.html
  by Peter J. Acklam. Accurate to about 4e-9 in absolute value.

  In: a value between 0 and 1 (not including endpoints)
  Out: the value y such that erfc(y) = input value
*************************************************************************/
double erfcinv(double x) {
  /* Misc constants */
  static const double tail_cutoff = 0.0485;
  static const double neg_twice_log_half = 1.386294361119891;
  // static const double sqrt_half_pi = 1.253314137315501;  /* For refinement */
  static const double scaling_const = -0.7071067811865475;

  /* Tail numerator */
  static const double tn0 = 2.938163982698783;
  static const double tn1 = 4.374664141464968;
  static const double tn2 = -2.549732539343734;
  static const double tn3 = -2.400758277161838;
  static const double tn4 = -3.223964580411365e-1;
  static const double tn5 = -7.784894002430293e-3;
  /* Tail denominator */
  static const double td1 = 3.754408661907416;
  static const double td2 = 2.445134137142996;
  static const double td3 = 3.224671290700398e-1;
  static const double td4 = 7.784695709041462e-3;

  /* Central numerator */
  static const double cn0 = 2.506628277459239;
  static const double cn1 = -3.066479806614716e1;
  static const double cn2 = 1.383577518672690e2;
  static const double cn3 = -2.759285104469687e2;
  static const double cn4 = 2.209460984245205e2;
  static const double cn5 = -3.969683028665376e1;
  /* Central denominator */
  static const double cd1 = -1.328068155288572e1;
  static const double cd2 = 6.680131188771972e1;
  static const double cd3 = -1.556989798598866e2;
  static const double cd4 = 1.615858368580409e2;
  static const double cd5 = -5.447609879822406e1;

  double p, q, r;

  if (x < tail_cutoff) {
    p = sqrt(-2 * log(x) + neg_twice_log_half);
    r = (tn0 + p * (tn1 + p * (tn2 + p * (tn3 + p * (tn4 + p * tn5))))) /
        (1.0 + p * (td1 + p * (td2 + p * (td3 + p * td4))));
  } else {
    p = 0.5 * x - 0.5;
    q = p * p;
    r = p * (cn0 + q * (cn1 + q * (cn2 + q * (cn3 + q * (cn4 + q * cn5))))) /
        (1.0 + q * (cd1 + q * (cd2 + q * (cd3 + q * (cd4 + q * cd5)))));
  }
  return scaling_const * r;
  /*
  Use the code below to refine to macine precision.  Rather slow, though.
  p = (erfc(scaling_const*r)-x)*sqrt_half_pi*exp(0.5*r*r);
  return scaling_const*(r - p/(1 + r*p/2));
  */
}

/*************************************************************************
poisson_dist:
  In: mean value
      random number distributed uniformly between 0 and 1
  Out: integer sampled from the Poisson distribution.
  Note: This does not sample the CDF.  Instead, it works its way outwards
        from the peak of the PDF.  Kinda weird.  It is not the case
        that low values of the random number will give low values.  It
        is also not super-efficient, but it works.
*************************************************************************/
int poisson_dist(double lambda, double p) {
  int i, lo, hi;
  double plo, phi, pctr;
  double lambda_i;

  i = (int)lambda;
  pctr = exp(-lambda + i * log(lambda) -
             lgamma(i + 1)); /* Highest probability bin */

  if (p < pctr)
    return i;

  lo = hi = i;
  plo = phi = pctr; /* Start at highest-probability bin and work outwards */

  p -= pctr;
  lambda_i = 1.0 / lambda;
  while (p > 0) /* Keep going until we exhaust probabilities */
  {
    if (lo > 0) /* We still have a low tail, test it */
    {
      plo *= lo * lambda_i; /* Recursive formula for p for this bin */
      lo--;
      if (p < plo)
        return lo;
      p -= plo;
    }
    /* Always test the high tail (it's infinite) */
    hi++;
    phi = phi * lambda / hi; /* Recursive formula for p for this bin */
    if (p < phi)
      return hi;
    p -= phi + DBL_EPSILON; /* Avoid infinite loop from poor roundoff */
  }

  /* should never get here */
  assert(false);
}

/*************************************************************************
byte_swap:
  In: array of bytes to be swapped
      size of this array
  Out: array of bytes swapped so that the last byte becomes the first one, etc.
       No return value
*************************************************************************/
void byte_swap(void *data, int size) {
  if (size < 2){
    return;
  }

  unsigned char temp;
  unsigned char *c = (unsigned char *)data;
  for (int i = 0, j = size - 1; i < j; i++, j--) {
    temp = c[i];
    c[i] = c[j];
    c[j] = temp;
  }
}

/************************************************************************\
                   Begin Rex's string matching code
\************************************************************************/

/*************************************************************************
  wild strings have wildcard characters * ? [...] and \ as an escape char
  feral strings have the same except no * character
  tame strings don't have any wildcard characters
*************************************************************************/

/* Measure the length of the tame string matched by a feral string of
 * length<=n*/
int feral_strlenn(char *feral, int n) {
  int real_n = 0;
  int i;
  for (i = 0; i < n; i++) {
    if (feral[i] == '\\') {
      i++;
      if (feral[i] == '\0')
        return real_n;
    } else if (feral[i] == '[') {
      while (i < n && feral[i] != ']') {
        if (feral[i] == '\0')
          return real_n;
        if (feral[i] == '\\') {
          i += 2;
          if (i > n || feral[i - 1] == '\0')
            return real_n;
        } else if (feral[i] == '-') {
          i += 2;
          if (i > n || feral[i - 1] == '\0')
            return real_n;
        } else
          i++;
      }
    } else if (feral[i] == '\0')
      return real_n;
    real_n++;
  }
  return real_n;
}

/* Check if the first n characters in the feral string is
an abbreviation for the tame string (i.e. matches the first
part of the tame string); return 0 if not found or the number
of matched characters if they are found */
int is_feral_nabbrev(char *feral, int n, char *tame) {
  char c, cc;
  int i = 0;
  int nfound = 0;
  int ok;

  if (n <= 0)
    return 0;

  while (*tame != '\0') {
    if (feral[i] == '[') /* Try to match character set */
    {
      i++;
      ok = 0;
      while (i < n && feral[i] != ']') {
        c = feral[i++];
        if (c == '\0')
          return 0; /* Malformed feral string */
        if (c == '\\') {
          if (i >= n)
            return 0; /* Malformed feral string */
          c = feral[i++];
          if (c == '\0')
            return 0; /* Malformed feral string */
        }
        if (i < n && feral[i] == '-') {
          i++;
          if (i >= n)
            return 0; /* Malformed feral string */
          cc = feral[i++];
          if (cc == '\0')
            return 0; /* Malformed feral string */
          if (cc == '\\') {
            if (i >= n)
              return 0; /* Malformed feral string */
            cc = feral[i++];
            if (cc == '\0')
              return 0; /* Malformed feral string */
          }
          if (c <= *tame && *tame <= cc) {
            ok = 1;
            break;
          }
        } else if (c == *tame) {
          ok = 1;
          break;
        }
      }
      if (i >= n)
        return 0; /* Malformed feral string */
      if (!ok)
        return 0;                      /* Set never matched */
      tame++;                          /* Matched */
      while (i < n && feral[i] != ']') /* Find trailing ] */
      {
        if (feral[i] == '\0')
          return 0; /* Malformed feral string */
        if (feral[i] == '\\') {
          i += 2;
          if (i > n || feral[i - 1] == '\0')
            return 0; /* Malformed feral string */
        } else
          i++;
      }
      if (i >= n)
        return 0; /* Malformed feral string */
      i++;
    } else /* Match single possibly escaped character */
    {
      c = feral[i++];
      if (c == '\\') {
        if (i >= n)
          return 0; /* Malformed feral string */
        c = feral[i++];
        if (c != *tame++)
          return 0; /* Mismatch */
      } else if (c != *tame++ && c != '?')
        return 0; /* Mismatch */
    }
    nfound++;
    if (i >= n)
      return nfound; /* Ran out of feral string--it's an abbreviation! */
  }

  return 0; /* Ran out of tame string with feral string left--not abbrev */
}

/* Find a substring of a tame haystack string that matches the first n
characters of the feral string needle (same syntax as strstr except using
a feral string with a length delimiter). Returns NULL if matching substring not
found. */

char *feral_strstrn(char *tame_haystack, char *feral_needle, int n) {
  char c = 0;
  char cc;
  char set[256];
  int isset = 0;
  int i, j;
  int scoot = 0;

  for (i = 0; i < n; i++)
    if (feral_needle[i] == '\0')
      break;
  n = i;

  /* Toss leading ?'s */
  i = 0;
  while (feral_needle[i] == '?' && i < n && *tame_haystack != '\0') {
    i++;
    tame_haystack++;
    scoot++;
  }

  if (i >= n)
    return tame_haystack - scoot;

  /* Beginning of needle is either a single character to match or a set of
   * characters */
  /* Efficiently search character set if it's first */
  if (feral_needle[i] == '[') {
    isset = 1;
    memset(set, 0, 256);
    set[0] = 1;
    i++;
    while (i < n && feral_needle[i] != ']') {
      c = feral_needle[i++];
      if (feral_needle[i] == '\0')
        return NULL; /* Can't match broken pattern */
      if (c == '\\') {
        if (i >= n)
          return NULL; /* Can't match broken pattern */
        c = feral_needle[i++];
      }
      if (i < n && feral_needle[i] == '-') {
        i++;
        if (i >= n)
          return NULL; /* Can't match broken pattern */
        cc = feral_needle[i++];
        if (cc == '\0')
          return NULL; /* Can't match broken pattern */
        if (cc == '\\') {
          if (i >= n)
            return NULL; /* Can't match broken pattern */
          cc = feral_needle[i++];
          if (cc == '\0')
            return NULL; /* Can't match broken pattern */
        }
        for (j = (int)c; j <= (int)cc; j++)
          set[j] = 1;
      } else
        set[(int)c] = 1;
    }
    if (i >= n)
      return NULL; /* Can't match broken pattern */
    i++;           /* Skip ] */
  } else {
    c = feral_needle[i++];
    if (c == '\\') {
      if (i >= n)
        return NULL; /* Can't match broken pattern */
      c = feral_needle[i++];
    }
    if (c == '\0')
      return NULL; /* Can't match broken pattern */
  }

  /* Match needle with haystack */
  while (*tame_haystack != '\0') {
    /* Try to match the first non-'?' character in needle with haystack */
    if (isset) /* Find next position in haystack that matches a set of
                  characters */
    {
      while (!set[(int)*tame_haystack])
        tame_haystack++;
      if (*tame_haystack == '\0')
        return NULL;
    } else /* Find next position in haystack that matches a single character */
    {
      while (*tame_haystack != c && *tame_haystack != '\0')
        tame_haystack++;
      if (*tame_haystack == '\0')
        return NULL;
    }

    if (i == n)
      return tame_haystack - scoot;
    else if (is_feral_nabbrev(feral_needle + i, n - i,
                              tame_haystack + 1)) /* Try to match the rest of
                                                     the needle */
    {
      return tame_haystack - scoot;
    }

    tame_haystack++;
  }

  return NULL;
}

/* Returns 1 if the wildcard string wild matches the tame string tame */
int is_wildcard_match(char *wild, char *tame) {
  int nstars;
  int n;

  if (*wild == '\0' && *tame == '\0')
    return 1;

  for (n = 0, nstars = 0; wild[n] != '\0'; n++) {
    if (wild[n] == '[') {
      n++;
      while (wild[n] != '\0' && wild[n] != ']') {
        if (wild[n] == '\\') {
          n++;
          if (wild[n] == '\0')
            return 0; /* Malformed wild string */
        }
        n++;
      }
      if (wild[n] == '\0')
        return 0; /* Malformed wild string */
    } else if (wild[n] == '\\') {
      n++;
      if (wild[n] == '\0')
        return 0; /* Malformed wild string */
    } else if (wild[n] == '*')
      nstars++;
  }

  if (nstars == 0)
    return (is_feral_nabbrev(wild, n, tame) == (int)strlen(tame));
  else {
    int staridx[nstars];
    int idxA[nstars + 1], idxB[nstars + 1];
    char *m;
    int nidx;
    int i, j;
    int tail_len;
    int old_length;

    for (i = n = 0; wild[n] != '\0'; n++) {
      if (wild[n] == '[') {
        do {
          n++;
          if (wild[n] == '\\')
            n++;
        } while (wild[n] != ']');
      } else if (wild[n] == '\\')
        n++;
      else if (wild[n] == '*')
        staridx[i++] = n;
    }

    for (i = 0; i < nstars && staridx[i] == i; i++) {
    } /* Skip over '*'s at the beginning of wild string */

    if (i >= nstars)
      return 1; /* All stars, of course it matches */

    if (i == 0) /* First character is not a star */
    {
      j = is_feral_nabbrev(wild, staridx[0], tame);
      if (j == 0)
        return 0; /* Didn't match start string */

      tame += j; /* Matched first j characters, toss them */

      wild += staridx[0]; /* And advance to star */
      n -= staridx[0];
      for (i = nstars - 1; i >= 0; i--)
        staridx[i] -= staridx[0];
    }

    if (staridx[nstars - 1] < n - 1) /* Last character is not a star */
    {
      j = staridx[nstars - 1] + 1;
      tail_len = feral_strlenn(wild + j, n - j);

      j = is_feral_nabbrev(wild + j, n - j, tame + (strlen(tame) - tail_len));
      if (j == 0)
        return 0; /* Didn't match tail string */
    } else
      tail_len = 0;

    /* Head and tail are matched, if any.  Now build rest of list to match */
    nidx = 0;
    for (i = 1; i < nstars; i++) {
      idxA[nidx] = staridx[i - 1] + 1;
      idxB[nidx] = staridx[i];
      nidx++;
    }

    /* And now we match all the pieces */
    old_length = 0;
    m = tame;
    for (i = 0; i < nidx; i++) {
      idxB[i] -= idxA[i]; /* Calculate length of feral string */

      if (idxB[i] == 0)
        continue; /* Just more stars */

      m = m + old_length;

      m = feral_strstrn(m, wild + idxA[i], idxB[i]);
      if (m == NULL)
        return 0; /* Couldn't find appropriate substring */
      old_length = feral_strlenn(wild + idxA[i], idxB[i]);
    }

    m = m + old_length;

    if (strlen(m) < (size_t)tail_len)
      return 0; /* Ran over tail string--no good */

    return 1;
  }
}

/************************************************************************\
                    End Rex's string matching code
\************************************************************************/

/*************************************************************************
initialize_string_buffer:
  Initialize the state of a string buffer.

  In: struct string_buffer *sb - the string buffer
      int maxstr - the maximum number of strings to add to this buffer
  Out: 0 on success; 1 if allocation fails.  All fields in the buffer are
       filled in.
**************************************************************************/
int initialize_string_buffer(struct string_buffer *sb, int maxstr) {
  if (maxstr > 0) {
    if ((sb->strings = (char **)allocate_ptr_array(maxstr)) == NULL)
      mcell_allocfailed("Failed to allocate buffer of %d strings.", maxstr);
  }

  sb->max_strings = maxstr;
  sb->n_strings = 0;
  return 0;
}

/*************************************************************************
destroy_string_buffer:
  Destroy the state of a string buffer.

  In: struct string_buffer *sb - the string buffer
  Out: The fields of the string buffer are freed, as are all strings
       added to the string buffer.  The string buffer itself is not freed.
**************************************************************************/
int destroy_string_buffer(struct string_buffer *sb) {
  if (sb->strings) {
    free_ptr_array((void **)sb->strings, sb->max_strings);
  }
  sb->strings = NULL;
  sb->max_strings = 0;
  sb->n_strings = 0;
  return 0;
}

/*************************************************************************
add_string_to_buffer:
    Add a string to the string buffer.  The string becomes "owned" by the
    buffer if this function is successful.  If this function fails, it is the
    caller's responsibility to free it.

        In: struct string_buffer *sb - the string buffer
            char *str - the string to add
        Out: 0 on success; 1 if the string is not stored because the buffer is
             full already.
**************************************************************************/
int add_string_to_buffer(struct string_buffer *sb, char *str) {
  if (sb->n_strings >= sb->max_strings) {
    mcell_internal_error("Attempt to overrun string buffer (max fill is %d).",
                         sb->max_strings);
    /*return 1;*/
  }

  sb->strings[sb->n_strings++] = str;
  return 0;
}

/*******************************************************************
 Pointer hashes implementation
*******************************************************************/

/*************************************************************************
  Initialize a pointer hash to a given initial size.  Returns 0 on success.
  Note that the desired table size may be exceeded.  Presently, the
  implementation always rounds the table size up to the nearest integer power
  of two.

  In:  struct pointer_hash *ht - the hash table to initialize
       int size - the desired table size
  Out: 0 on success; 1 if memory allocation fails.
**************************************************************************/
int pointer_hash_init(struct pointer_hash *ht, int size) {
  assert(size >= 0);

  memset(ht, 0, sizeof(struct pointer_hash));

  /* Don't allow size 0 tables.  Otherwise, the key & (tablesize-1) trick fails
   * spectacularly. */
  if (size == 0) {
    ++size;
  }

  /* Fold size to next larger power of 2 if it isn't a power of 2 already */
  if ((size) & (size - 1)) {
    size |= (size >> 1);
    size |= (size >> 2);
    size |= (size >> 4);
    size |= (size >> 8);
    size |= (size >> 16);
    ++size;
  }

  /* Initialize fields */
  ht->num_items = 0;
  ht->table_size = size;
  if ((ht->hashes = (unsigned int *)malloc(sizeof(unsigned int) *size)) ==
          NULL ||
      (ht->keys = (void const **)malloc(sizeof(void const *) *size)) == NULL ||
      (ht->values = (void **)malloc(sizeof(void *) * size)) == NULL)
    goto failure;

  /* Make sure our table starts out empty */
  memset(ht->hashes, 0, sizeof(unsigned int) * size);
  memset(ht->keys, 0, sizeof(void *) * size);
  memset(ht->values, 0, sizeof(void *) * size);

  return 0;

failure:
  pointer_hash_destroy(ht);
  return 1;
}

/*************************************************************************
  Destroy a pointer hash, freeing all memory associated with it.

  In:  struct pointer_hash *ht - the hash table to destroy
  Out: hash table is destroyed and associated memory is freed
**************************************************************************/
void pointer_hash_destroy(struct pointer_hash *ht) {
  if (ht->hashes)
    free(ht->hashes);
  if (ht->keys)
    free(ht->keys);
  if (ht->values)
    free(ht->values);
  ht->num_items = 0;
  ht->table_size = 0;
  ht->hashes = NULL;
  ht->keys = NULL;
  ht->values = NULL;
}

/*************************************************************************
  Manually resize a pointer hash to have at least 'new_size' bins.  New size
  may exceed requested size.  Works exactly like pointer_hash_init, except
  values are copied from the old table to the new one.

  In:  struct pointer_hash *ht - the hash table to resize
       int new_size - the desired table size
  Out: 0 on success; 1 if memory allocation fails.

  On failure, the hash table is left unchanged from its prior state.
**************************************************************************/
int pointer_hash_resize(struct pointer_hash *ht, int new_size) {
  if (new_size == ht->table_size)
    return 0;
  if (new_size < ht->num_items)
    return 1;

  /* Save old hash, allocate new hash */
  struct pointer_hash old = *ht;
  if (pointer_hash_init(ht, new_size)) {
    *ht = old;
    return 1;
  }

  /* Copy items over to new hash */
  for (int old_item_idx = 0; old_item_idx < old.table_size; ++old_item_idx) {
    if (old.keys[old_item_idx] == NULL)
      continue;

    if (pointer_hash_add(ht, old.keys[old_item_idx], old.hashes[old_item_idx],
                         old.values[old_item_idx]))
      goto failure;
  }

  /* Destroy the old hash */
  pointer_hash_destroy(&old);
  return 0;

failure:
  pointer_hash_destroy(ht);
  *ht = old;
  return 1;
}

/*************************************************************************
  Add a value to a pointer hash.  If a previous item was added for that key,
  the new value will replace the old value.

  In:  struct pointer_hash *ht - the hash table to receive a new item
       void const *key - the key
       unsigned int keyhash - the hash value associated with the key
       void *value - the value to store
  Out: 0 on success; 1 if memory allocation fails.
**************************************************************************/
int pointer_hash_add(struct pointer_hash *ht, void const *key,
                     unsigned int keyhash, void *value) {
  /* In case pointer hash was initialized using memset, we'll allocate space
   * on-demand.  */
  if (ht->table_size == 0) {
    if (pointer_hash_resize(ht, 2))
      return 1;
  }

  /* Make sure our table always has free space */
  if (ht->num_items >= (ht->table_size >> 1)) {
    if (pointer_hash_resize(ht, ht->table_size << 1))
      return 1;
  }

  /* Scan over entries until the end of the table */
  unsigned int start_index = keyhash & (ht->table_size - 1);
  for (unsigned int i = 0; i < (unsigned int)ht->table_size; i++) {
    unsigned int cur_index = (start_index + i) & (ht->table_size - 1);

    /* Found an old value for this key.  Replace it.  Do not increment the item
     * count. */
    if (ht->keys[cur_index] == key) {
      ht->values[cur_index] = value;
      return 0;
    }

    /* Found an empty slot.  Fill it. */
    if (ht->keys[cur_index] == NULL) {
      ht->hashes[cur_index] = keyhash;
      ht->keys[cur_index] = key;
      ht->values[cur_index] = value;
      goto done;
    }
  }

  /* This shouldn't happen unless the pointer hash is corrupted, since we've
   * already ensured that we have enough space */
  return 1;

done:
  ++ht->num_items;
  return 0;
}

/*************************************************************************
  Look up a value in a pointer hash.  Returns default_value if no item was
  found, or if the value associated with the key was default_value.

  In:  struct pointer_hash *ht - the hash table to search
       void const *key - the key to find
       unsigned int keyhash - the hash value associated with the key
       void const *default_value - value to return if item is not found
  Out: the value, or default_value if not found
**************************************************************************/
void *pointer_hash_lookup_ext(struct pointer_hash const *ht, void const *key,
                              unsigned int keyhash, void *default_value) {
  /* Empty table.  Not found. */
  if (ht->table_size == 0)
    return default_value;

  /* Search from start position to end of table */
  unsigned int start_index = keyhash & (ht->table_size - 1);
  for (unsigned int i = 0; i < (unsigned int)ht->table_size; i++) {
    unsigned int cur_index = (start_index + i) & (ht->table_size - 1);

    /* Empty slot - key not found. */
    if (ht->keys[cur_index] == NULL)
      return default_value;

    /* Found our value. */
    if (ht->keys[cur_index] == key)
      return ht->values[cur_index];
  }

  /* Worst case scenario.  We searched the entire table.  Shouldn't happen with
   * the current tuning parameters, which do not permit the table to be full.
   */
  return default_value;
}

/*************************************************************************
  Remove a value from a pointer hash.  Returns 0 if the item was
  successfully removed, or 0 if the item was not found.  Note that a
  NULL key is not allowed in a pointer hash.

  In:  struct pointer_hash *ht - the hash table to remove from
       void const *key - the key to remove
       unsigned int keyhash - the hash value associated with the key
  Out: 0 on success; 1 on failure

  Regardless of success or failure, the table will remain in a valid
  state.
**************************************************************************/
int pointer_hash_remove(struct pointer_hash *ht, void const *key,
                        unsigned int keyhash) {
  /* If the key is NULL, we're done */
  if (key == NULL)
    return 1;

  /* If the table is empty, we're done */
  if (ht->table_size == 0)
    return 1;

  int start = keyhash & (ht->table_size - 1); /* Where do we start? */
  int next_slot = -1;                         /* Next vacant slot to fill */

  /* Scan through the table to remove the item, and resituate any
   * entries which would, due to the streaming collision-avoidance
   * approach, be orphaned.
   */
  int cur = start;
  do {
    /* If we found the item to remove, remove it */
    if (key != NULL && ht->keys[cur] == key) {
      /* Clear this value */
      ht->keys[cur] = NULL;
      ht->values[cur] = NULL;
      ht->hashes[cur] = 0;
      --ht->num_items;

      if (ht->table_size > (ht->num_items << 2)) {
        /* If resizing the pointer hash succeeded, we don't need to do
         * any more work, since it will have prevented orphans. */
        if (!pointer_hash_resize(ht, ht->num_items << 1))
          return 0;
      }

      /* We may need to move something into this empty slot, so
       * remember where we were. */
      next_slot = cur;

      /* This indicates that we've already removed the value and are
       * now just fixing the table contents */
      key = NULL;
    }

    /* We found a NULL entry, so stop scanning -- the table is in a
     * valid state */
    else if (ht->keys[cur] == NULL) {
      if (key == NULL)
        return 0;
      else
        return 1;
    }

    /* This entry may need to be reseated */
    else if (next_slot != -1) {
      /* This is the slot where we'd start looking for this entry */
      int desiredSlot = ht->hashes[cur] & (ht->table_size - 1);

      /* If we would have hit our newly introduced NULL entry after we
       * started searching for this entry, but before we found what we
       * were looking for, resituate this entry.
       */
      if ((next_slot < cur &&
           (desiredSlot <= next_slot || desiredSlot > cur)) ||
          (next_slot > cur &&
           (desiredSlot <= next_slot && desiredSlot > cur))) {
        /* Move this entry to its new location */
        ht->hashes[next_slot] = ht->hashes[cur];
        ht->keys[next_slot] = ht->keys[cur];
        ht->values[next_slot] = ht->values[cur];

        /* Free this slot */
        ht->hashes[cur] = 0;
        ht->keys[cur] = NULL;
        ht->values[cur] = NULL;

        /* We may need to move something into this empty slot, so
         * remember where we were. */
        next_slot = cur;
      }
    }

    /* Advance to the next entry */
    if (++cur == ht->table_size)
      cur = 0;
  } while (cur != start);

  /* If we found our entry, success, else failure. */
  if (key == NULL)
    return 0;
  else
    return 1;
}

/*************************************************************************
remove_both_duplicates:
  In: sorted linked list containing void pointers
  Out: if linked list contains two duplicates they are both removed
       Returns number of unique items in the linked list
       after removal
  Note: linked list should be sorted in advance.
        Also this function is used in the way that we do not expect more
        than two duplicates.
        Example: input = (a,b,c,d,d, e,f,g)
                 output = (a,b,c,e,f,g)
*************************************************************************/
int remove_both_duplicates(struct void_list **head) {
  struct void_list *curr = *head, *tmp, *prev, *next_Next;
  int count = 0;

  tmp = curr;
  prev = NULL;
  while ((tmp != NULL) && (tmp->next != NULL)) {
    if (tmp->data == tmp->next->data) {
      next_Next = tmp->next->next;
      free(tmp->next);
      free(tmp);
      tmp = next_Next;
      if (prev != NULL) {
        prev->next = tmp;
      } else {
        curr = tmp;
      }
    } else {
      prev = tmp;
      tmp = tmp->next; /* only advance if there is no deletion */
    }
  }

  *head = curr;

  for (tmp = *head; tmp != NULL; tmp = tmp->next) {
    count++;
  }

  return count;
}

/*********************************************************************
delete_void_list:
  In: linked list
  Out: none.  The memory is freed.
*********************************************************************/
void delete_void_list(struct void_list *head) {
  struct void_list *nnext;
  while (head != NULL) {
    nnext = head->next;
    free(head);
    head = nnext;
  }
}


/*************************************************************************
 double_cmp:
    Comparison function for doubles, to be passed to qsort.

 In:  i1: first item for comparison
      i2: second item for comparison
 Out: -1, 0, or 1 if the first item is less than, equal to, or greater than the
      second, resp.
*************************************************************************/
int double_cmp(void const *i1, void const *i2) {
  double const *d1 = (double const *)i1;
  double const *d2 = (double const *)i2;
  if (*d1 < *d2)
    return -1;
  else if (*d1 > *d2)
    return 1;
  else
    return 0;
}

/**************************************************************************
is_string_present_in_string_array:
  In: string "str"
      array of strings "strings"
      length of the array of strings "length"
  Out: Return 1 if string "str" is present in the array of strings "strings",
       and 0 otherwise.
**************************************************************************/
int is_string_present_in_string_array(char * str, char ** strings, int length)
{
  int found = 0, i;

  for (i = 0; i < length; i++) {
    if (strcmp(str, strings[i]) == 0) {
      found = 1;
      break;
    }
  }

  return found;
}

/******************************************************************************
 convert_seconds_to_iterations:

 Converts time in seconds to iterations. This is offset by the number of
 iterations at the start of the simulation, which is usually 0 unless you're
 checkpointing. Alternatively, this function could be thought of as converting
 seconds to the internal/scaled time (i.e. the time in units of the timestep),
 but this doesn't really make sense if one changes the timestep during a
 checkpointing event (e.g. 1e-6 to 1e-7 seconds).
 *****************************************************************************/
double convert_seconds_to_iterations(
    long long start_iterations,
    double time_step_seconds,
    double simulation_start_seconds,
    double seconds) {
  double delta_iterations =
    (seconds - simulation_start_seconds) / time_step_seconds;
  return (start_iterations + delta_iterations);
}

/******************************************************************************
 convert_iterations_to_seconds:
 
 As you might imagine, this is essentially the inverse of
 convert_seconds_to_iterations.

 NOTE: Do not use iteration values that happened in the past or delta iteration
 values. Only input values that are greater the starting number of iterations
 In other words, use either the current number of iterations or some number
 corresponding to some time in the future.
 *****************************************************************************/
double convert_iterations_to_seconds(
    long long start_iterations,
    double time_step_seconds,
    double simulation_start_seconds,
    double iterations) {
  // Probably should add an assert here
  double delta_time = (iterations - start_iterations) * time_step_seconds;
  return (simulation_start_seconds + delta_time);
}

/*************************************************************************
 generate_range:
    Generate a num_expr_list containing the numeric values from start to end,
    incrementing by step.

 In:  state: the simulation state
      list:  destination to receive list of values
      start: start of range
      end:   end of range
      step:  increment
 Out: 0 on success, 1 on failure.  On success, list is filled in.
*************************************************************************/
int generate_range(struct num_expr_list_head *list, double start, double end,
                   double step) {
  list->value_head = NULL;
  list->value_tail = NULL;
  list->value_count = 0;
  list->shared = 0;

  if (step > 0) {
    /* JW 2008-03-31: In the guard on the loop below, it seems to me that
     * the third condition is redundant with the second.
     */
    for (double tmp_dbl = start;
         tmp_dbl < end || !distinguishable(tmp_dbl, end, EPS_C) ||
             fabs(end - tmp_dbl) <= EPS_C;
         tmp_dbl += step) {
      if (advance_range(list, tmp_dbl))
        return 1;
    }
  } else /* if (step < 0) */
  {
    /* JW 2008-03-31: In the guard on the loop below, it seems to me that
     * the third condition is redundant with the second.
     */
    for (double tmp_dbl = start;
         tmp_dbl > end || !distinguishable(tmp_dbl, end, EPS_C) ||
             fabs(end - tmp_dbl) <= EPS_C;
         tmp_dbl += step) {
      if (advance_range(list, tmp_dbl))
        return 1;
    }
  }
  return 0;
}

// Maybe move this somewhere else
int advance_range(struct num_expr_list_head *list, double tmp_dbl) {
  struct num_expr_list *nel;
  nel = CHECKED_MALLOC_STRUCT(struct num_expr_list, "numeric list");
  if (nel == NULL) {
    free_numeric_list(list->value_head);
    list->value_head = list->value_tail = NULL;
    return 1;
  }
  nel->value = tmp_dbl;
  nel->next = NULL;

  ++list->value_count;
  if (list->value_tail != NULL)
    list->value_tail->next = nel;
  else
    list->value_head = nel;
  list->value_tail = nel;
  return 0;
}

/*************************************************************************
 free_numeric_list:
    Free a num_expr_list.

 In:  nel:  the list to free
 Out: all elements are freed
*************************************************************************/
void free_numeric_list(struct num_expr_list *nel) {
  while (nel != NULL) {
    struct num_expr_list *n = nel;
    nel = nel->next;
    free(n);
  }
}

