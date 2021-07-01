/******************************************************************************
 *
 * Copyright (C) 2019 by
 * The Salk Institute for Biological Studies and
 * Pittsburgh Supercomputing Center, Carnegie Mellon University
 *
 * Use of this source code is governed by an MIT-style
 * license that can be found in the LICENSE file or at
 * https://opensource.org/licenses/MIT.
 *
**/


#ifndef __MAP_C_H__
#define __MAP_C_H__

#define MAP_MISSING -3  /* No such element */
//#define MAP_FULL -2   /* Hashmap is full */
//#define MAP_OMEM -1   /* Out of Memory */
#define MAP_OK 0  /* OK */

/*
 * any_t is a pointer.  This allows you to put arbitrary structures in
 * the hashmap.
 */
typedef void *any_t;

/*
 * PFany is a pointer to a function that can take two any_t arguments
 * and return an integer. Returns status code..
 */
typedef int (*PFany)(any_t, any_t);

/*
 * map_t is a pointer to an internally maintained data structure.
 * Clients of this package do not need to know how hashmaps are
 * represented.  They see and manipulate only map_t's.
 */
typedef any_t map_t;

/*
 * Return an empty hashmap. Returns NULL if empty.
*/
map_t hashmap_new(void);

void hashmap_clear(map_t in);


//JJT: simplified functions without hash calculation (it is precalculated)
int hashmap_put_nohash(map_t in, unsigned long key, unsigned long key_hash, any_t value);
int hashmap_get_nohash(map_t in, unsigned long key, unsigned long key_hash, any_t *arg);

unsigned long crc32(const unsigned char *s, unsigned int len);

#endif //__MAP_C_H__
