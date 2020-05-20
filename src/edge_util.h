/******************************************************************************
 *
 * Copyright (C) 2006-2017,2020 by
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
