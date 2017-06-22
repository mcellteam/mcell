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

#include <stdio.h>

#define COUNT_OF(arr) (sizeof((arr)) / sizeof((arr[0])))

#define BLOCK_SIZE 10000

struct num_expr_list_head {
  struct num_expr_list *value_head;
  struct num_expr_list *value_tail;
  int value_count;
  int shared;
};

struct iteration_counter {
  long long *iterations; /* array of iteration numbers, should be
                            memory allocated */
  int max_iterations; /* size of the above array */
  int n_iterations;   /* number of the filled positions in an above array */
};

struct string_buffer {
  char **strings;  /* array of strings, should be memory allocated */
  int max_strings; /* size of the above array */
  int n_strings;   /* number of the filled positions in an above array */
};

struct bit_array {
  int nbits;
  int nints;
  /* Bit array data runs off the end of this struct */
};

struct bit_array *new_bit_array(int bits);
struct bit_array *duplicate_bit_array(struct bit_array *old);
int get_bit(struct bit_array *ba, int idx);
void set_bit(struct bit_array *ba, int idx, int value);
void set_bit_range(struct bit_array *ba, int idx1, int idx2, int value);
void set_all_bits(struct bit_array *ba, int value);
void bit_operation(struct bit_array *ba, struct bit_array *bb, char op);
int count_bits(struct bit_array *ba);
void free_bit_array(struct bit_array *ba);

int bisect(double *list, int n, double val);
int bisect_near(double *list, int n, double val);
int bisect_high(double *list, int n, double val);

int distinguishable(double a, double b, double eps);
int is_reverse_abbrev(char *abbrev, char *full);

struct void_list {
  struct void_list *next;
  void *data;
};

struct void_list *void_list_sort(struct void_list *vl);
struct void_list *void_list_sort_by(struct void_list *vl,
                                    int (*leq)(void *, void *));
int remove_both_duplicates(struct void_list **head);
void delete_void_list(struct void_list *head);

int void_array_search(void **array, int n, void *to_find);
int void_ptr_compare(void const *v1, void const *v2);

unsigned int *allocate_uint_array(int size, unsigned int value);
void **allocate_ptr_array(int size);
void free_ptr_array(void **pa, int count);

struct num_expr_list;
void free_num_expr_list(struct num_expr_list *nel);

int dir_exists(char const *path);
int is_writable_dir(char const *path);
int mkdirs(char const *path);
int make_parent_dir(char const *path);

FILE *open_file(char const *fname, char const *mode);

double erfcinv(double v);

int poisson_dist(double lambda, double p);

void byte_swap(void *data, int size);

int feral_strlenn(char *feral, int n);
int is_feral_nabbrev(char *feral, int n, char *tame);
char *feral_strstrn(char *tame_haystack, char *feral_needle, int n);
int is_wildcard_match(char *wild, char *tame);

int initialize_string_buffer(struct string_buffer *sb, int maxstr);
int destroy_string_buffer(struct string_buffer *sb);
int add_string_to_buffer(struct string_buffer *sb, char *str);

double convert_seconds_to_iterations(
    long long start_iterations,
    double time_step_seconds,
    double simulation_start_seconds,
    double seconds);

double convert_iterations_to_seconds(
    long long start_iterations,
    double time_step_seconds,
    double simulation_start_seconds,
    double iterations);

/*******************************************************************
 Pointer hashes

   Essentially, pointer hashes are pointer -> pointer hash tables. There is no
   restriction on the type of pointer used for the key, but a hash value must
   be supplied along with the pointer whenever performing any operation that
   requires a key. Neither the key, nor the value pointer are ever dereferenced
   or freed by the pointer hash code.

      Usage:

        struct pointer_hash hash;
        if (pointer_hash_init(&hash)) { fail }
        pointer_hash_add(&hash, key1, hash1, value1);
        pointer_hash_add(&hash, key2, hash2, value2);
        ..
        pointer_hash_add(&hash, keyN, hashN, valueN);

        ..

        void *value1 = pointer_hash_lookup(&hash, key1, hash1);

        ..

        pointer_hash_destroy(&hash);

  Note: The pointer hash allocates some memory, which will be orphaned
  if pointer_hash_destroy is not called on the hash before it goes out
  of scope.

  The hash table collision strategy implemented by this structure is a
  streaming strategy, rather than a chaining strategy.  This means that
  when we try to insert a value into the n-th bin, if the n-th bin is
  already occupied, we will move to the n+1-th bin (modulo table size)
  until we find an insertion point.  When we do a lookup, therefore, we
  need to scan forward until we either find our key, or until we find
  an empty bin.  We can ensure that at least one of these conditions is
  met by always keeping the table size greater than the number of
  entries.

  Currently, the pointer hash is tuned to keep the table size between 2
  and 4 times the number of items it contains.
*******************************************************************/
struct pointer_hash {
  int num_items;        /* num items in table */
  int table_size;       /* size of table */
  unsigned int *hashes; /* hash values for each entry */
  void const **keys;    /* keys for each entry */
  void **values;        /* values for each entry */
};

/* Initialize a pointer hash to a given initial size.  Returns 0 on
 * success. */
int pointer_hash_init(struct pointer_hash *ht, int size);

/* Quickly clear all values from a pointer hash.  Does not free any
 * memory. */
void pointer_hash_clear(struct pointer_hash *ht);

/* Destroy a pointer hash, freeing all memory associated with it. */
void pointer_hash_destroy(struct pointer_hash *ht);

/* Manually resize a pointer hash to have at least 'new_size' bins.
 * New size may exceed requested size. */
int pointer_hash_resize(struct pointer_hash *ht, int new_size);

/* Add a value to a pointer hash.  If a previous item was added for
 * that key, the new value will replace the old value. */
int pointer_hash_add(struct pointer_hash *ht, void const *key,
                     unsigned int keyhash, void *value);

/* Look up a value in a pointer hash.  Returns NULL if no item was
 * found, or if the value associated with the key was NULL. */
#define pointer_hash_lookup(ht, key, keyhash)                                  \
  pointer_hash_lookup_ext(ht, key, keyhash, NULL)

/* Look up a value in a pointer hash.  Returns NULL if no item was
 * found, or if the value associated with the key was NULL. */
void *pointer_hash_lookup_ext(struct pointer_hash const *ht, void const *key,
                              unsigned int keyhash, void *default_value);

/* Remove a value from a pointer hash.  Returns 0 if the item was
 * successfully removed, or 0 if the item was not found.
 */
int pointer_hash_remove(struct pointer_hash *ht, void const *key,
                        unsigned int keyhash);

int double_cmp(void const *i1, void const *i2);

int is_string_present_in_string_array(char * str, char ** strings, int length);

int generate_range(struct num_expr_list_head *list, double start, double end,
                   double step);

int advance_range(struct num_expr_list_head *list, double tmp_dbl);

void free_numeric_list(struct num_expr_list *nel);

/*******************************************************************
 Inline min/max functions
*******************************************************************/

static inline double min2d(double x, double y) { return (x < y) ? x : y; }

static inline int min2i(int x, int y) { return (x < y) ? x : y; }

static inline double max2d(double x, double y) { return (x > y) ? x : y; }

static inline int max2i(int x, int y) { return (x > y) ? x : y; }

static inline double min3d(double x, double y, double z) {
  return (z < y) ? ((z < x) ? z : x) : ((y < x) ? y : x);
}

static inline int min3i(int x, int y, int z) {
  return (z < y) ? ((z < x) ? z : x) : ((y < x) ? y : x);
}

static inline double max3d(double x, double y, double z) {
  return (z > y) ? ((z > x) ? z : x) : ((y > x) ? y : x);
}

static inline int max3i(int x, int y, int z) {
  return (z > y) ? ((z > x) ? z : x) : ((y > x) ? y : x);
}

static inline long long min3ll(long long x, long long y, long long z) {
  return (z < y) ? ((z < x) ? z : x) : ((y < x) ? y : x);
}

static inline long long max3ll(long long x, long long y, long long z) {
  return (z > y) ? ((z > x) ? z : x) : ((y > x) ? y : x);
}

/* Return minimum value from the array of N doubles */
static inline double minNd(double *array, int N) {
  double smallest;
  N -= 2;
  for (smallest = array[N + 1]; N >= 0; N--) {
    if (array[N] < smallest)
      smallest = array[N];
  }
  return smallest;
}

/* Return minimum value from the array of N integers */
static inline int minNi(int *array, int N) {
  int smallest;
  N -= 2;
  for (smallest = array[N + 1]; N >= 0; N--) {
    if (array[N] < smallest)
      smallest = array[N];
  }
  return smallest;
}
