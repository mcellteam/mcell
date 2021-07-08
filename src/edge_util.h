/******************************************************************************
 *
 * Copyright (C) 2006-2017,2020 by
 * The Salk Institute for Biological Studies and
 * Pittsburgh Supercomputing Center, Carnegie Mellon University
 *
 * Use of this source code is governed by an MIT-style
 * license that can be found in the LICENSE file or at
 * https://opensource.org/licenses/MIT.
 *
******************************************************************************/

#pragma once

#include "mcell_structs.h"

/* Temporary data stored about an edge of a polygon */
struct poly_edge {
  struct poly_edge *next; /* Next edge in a hash table. */

  double v1x; /* X coord of starting point */
  double v1y; /* Y coord of starting point */
  double v1z; /* Z coord of starting point */
  double v2x; /* X coord of ending point */
  double v2y; /* Y coord of ending point */
  double v2z; /* Z coord of ending point */

  int face[2]; /* wall indices on side of edge */
  int edge[2]; /* which edge of wall1/2 are we? */
  int n;     /* How many walls share this edge? */
};

/* Hash table for rapid order-invariant lookup of edges. */
struct edge_hashtable {
  struct poly_edge *data; /* Array of polygon edges */

  int nkeys;    /* Length of array */
  int stored;   /* How many things do we have in the table? */
  int distinct; /* How many of those are distinct? */
};

int edge_equals(struct poly_edge *e1, struct poly_edge *e2);
int edge_hash(struct poly_edge *pe, int nkeys);

int ehtable_init(struct edge_hashtable *eht, int nkeys);
int ehtable_add(struct edge_hashtable *eht, struct poly_edge *pe);
void ehtable_kill(struct edge_hashtable *eht);
