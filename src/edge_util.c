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

#include <stdlib.h>
#include <iostream>

#include "edge_util.h"
#include "mcell_structs.h"
#include "sym_table.h"

#include "debug_config.h"

/**************************************************************************\
 ** Edge hash table section--finds common edges in polygons              **
\**************************************************************************/

/***************************************************************************
edge_equals:
  In: pointers to two poly_edge structs
  Out: Returns 1 if the edges are the same, 0 otherwise.
  Note: Orientation invariant, so an edge between vertex 1 and 2
        is the same as an edge between vertex 2 and 1.
***************************************************************************/

int edge_equals(struct poly_edge *e1, struct poly_edge *e2) {
  if ((!distinguishable(e1->v1x, e2->v1x, EPS_C)) &&
      (!distinguishable(e1->v1y, e2->v1y, EPS_C)) &&
      (!distinguishable(e1->v1z, e2->v1z, EPS_C)) &&
      (!distinguishable(e1->v2x, e2->v2x, EPS_C)) &&
      (!distinguishable(e1->v2y, e2->v2y, EPS_C)) &&
      (!distinguishable(e1->v2z, e2->v2z, EPS_C))) {
    return 1;
  }
  if ((!distinguishable(e1->v1x, e2->v2x, EPS_C)) &&
      (!distinguishable(e1->v1y, e2->v2y, EPS_C)) &&
      (!distinguishable(e1->v1z, e2->v2z, EPS_C)) &&
      (!distinguishable(e1->v2x, e2->v1x, EPS_C)) &&
      (!distinguishable(e1->v2y, e2->v1y, EPS_C)) &&
      (!distinguishable(e1->v2z, e2->v1z, EPS_C))) {
    return 1;
  }
  return 0;
}


/***************************************************************************
edge_hash:
  In: pe: pointer to a poly_edge struct
      nkeys: number of keys in the hash table
  Out: Returns a hash value between 0 and nkeys-1.
  Note: Orientation invariant, so a hash with the two endpoints swapped
        will be the same.
***************************************************************************/

int edge_hash(struct poly_edge *pe, int nkeys) {
  /* Get hash of X,Y,Z set of doubles for 1st and 2nd points */
  /* (Assume they're laid out consecutively in memory) */
  /* FIXME: This seems like a hack. Since poly_edge is a struct there's no
   * guarantee what the memory layout will be and the compiler may pad */
  unsigned int hashL = jenkins_hash((ub1 *)&(pe->v1x), 3 * sizeof(double));
  unsigned int hashR = jenkins_hash((ub1 *)&(pe->v2x), 3 * sizeof(double));

  return (hashL + hashR) %
         nkeys; /* ^ is symmetric so doesn't matter which is L and which is R */
}


/***************************************************************************
 *
 * create_new_poly_edge creates a new poly_edge and attaches it to the
 * past pointer to linked list of poly_edges
 *
 ***************************************************************************/
static struct poly_edge* create_new_poly_edge(struct poly_edge* list) {

  struct poly_edge *pei = CHECKED_MALLOC_STRUCT_NODIE(struct poly_edge, "polygon edge");
  if (pei == NULL) {
    return NULL;
  }

  pei->next = list->next;
  list->next = pei;
  pei->n = 0;
  pei->face[0] = -1;
  pei->face[1] = -1;
  pei->edge[0] = -1;
  pei->edge[1] = -1;
  return pei;
}


/***************************************************************************
ehtable_init:
  In: eht: pointer to an edge_hashtable struct
      nkeys: number of keys that the hash table uses
  Out: Returns 0 on success, 1 on failure.
       Hash table is initialized.
***************************************************************************/

int ehtable_init(struct edge_hashtable *eht, int nkeys) {
  eht->nkeys = nkeys;
  eht->stored = 0;
  eht->distinct = 0;
  eht->data =
      CHECKED_MALLOC_ARRAY_NODIE(struct poly_edge, nkeys, "edge hash table");
  if (eht->data == NULL)
    return 1;

  for (int i = 0; i < nkeys; i++) {
    eht->data[i].next = NULL;
    eht->data[i].n = 0;
    eht->data[i].face[0] = eht->data[i].face[1] = -1;
  }

  return 0;
}


/***************************************************************************
ehtable_add:
  In: pointer to an edge_hashtable struct
      pointer to the poly_edge to add
  Out: Returns 0 on success, 1 on failure.
       Edge is added to the hash table.
***************************************************************************/
int ehtable_add(struct edge_hashtable *eht, struct poly_edge *pe) {

#ifdef DEBUG_GEOM_OBJ_INITIALIZATION
    std::cout << "ehtable_add:\n";
    dump_poly_edge(0, pe);
#endif

  int i = edge_hash(pe, eht->nkeys);
  struct poly_edge *pep = &(eht->data[i]);

  while (pep != NULL) {
    if (pep->n == 0) /* New entry */
    {
      pep->n = 1;
      pep->face[0] = pe->face[0];
      pep->edge[0] = pe->edge[0];
      pep->v1x = pe->v1x;
      pep->v1y = pe->v1y;
      pep->v1z = pe->v1z;
      pep->v2x = pe->v2x;
      pep->v2y = pe->v2y;
      pep->v2z = pe->v2z;
      eht->stored++;
      eht->distinct++;
      return 0;
    }

    if (edge_equals(pep, pe)) /* This edge exists already ... */
    {
      if (pep->face[1] == -1) /* ...and we're the 2nd one */
      {
        pep->face[1] = pe->face[0];
        pep->edge[1] = pe->edge[0];
        pep->n++;
        eht->stored++;
        return 0;
      } else /* ...or we're 3rd and need more space */
      {
        if (pep->next != NULL) {
          if (edge_equals(pep->next, pe)) /* Space already there */
          {
            pep->n++;
            pep = pep->next;
            continue; /* Use space on next loop */
          }
        }

        struct poly_edge *pei = create_new_poly_edge(pep);
        if (pei == NULL) {
          return 1;
        }
        pep->n++;
        pep = pei;
        eht->distinct--; /* Not really distinct, just need more space */
      }
    } else if (pep->next != NULL) {
      pep = pep->next;
    } else { /* Hit end of list, so make space for use next loop. */
      struct poly_edge *pei = create_new_poly_edge(pep);
      if (pei == NULL) {
        return 1;
      }
      pep = pei;
    }
  }

  return 0;
}

/***************************************************************************
ehtable_kill:
  In: eht: pointer to an edge_hashtable struct
  Out: No return value.  Hashtable data is deallocated.
  Note: eht itself is not freed, since it isn't created with ehtable_init.
***************************************************************************/
void ehtable_kill(struct edge_hashtable *eht) {
  struct poly_edge *pe;

  for (int i = 0; i < eht->nkeys; i++) {
    while (eht->data[i].next != NULL) {
      pe = eht->data[i].next;
      eht->data[i].next = pe->next;
      free(pe);
    }
  }
  free(eht->data);
  eht->data = NULL;
  eht->nkeys = 0;
}
