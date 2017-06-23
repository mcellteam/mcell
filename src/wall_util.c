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

#include <math.h>
#include <string.h>
#include <stdlib.h>

#include "rng.h"
#include "logging.h"
#include "vector.h"
#include "util.h"
#include "init.h"
#include "sym_table.h"
#include "vol_util.h"
#include "mdlparse_util.h"
#include "grid_util.h"
#include "count_util.h"
#include "wall_util.h"
#include "react.h"
#include "strfunc.h"

/* tetrahedralVol returns the (signed) volume of the tetrahedron spanned by
 * the vertices a, b, c, and d.
 * The formula was taken from "Computational Geometry" (2nd Ed) by J. O'Rourke
 */
static double tetrahedralVol(struct vector3 *a, struct vector3 *b,
                             struct vector3 *c, struct vector3 *d) {
  return 1.0 / 6.0 * (-1 * (a->z - d->z) * (b->y - d->y) * (c->x - d->x) +
                      (a->y - d->y) * (b->z - d->z) * (c->x - d->x) +
                      (a->z - d->z) * (b->x - d->x) * (c->y - d->y) -
                      (a->x - d->x) * (b->z - d->z) * (c->y - d->y) -
                      (a->y - d->y) * (b->x - d->x) * (c->z - d->z) +
                      (a->x - d->x) * (b->y - d->y) * (c->z - d->z));
}

/* abs_max_2vec picks out the largest (absolute) value found among two vectors
 * (useful for properly handling floating-point rounding error). */
static inline double abs_max_2vec(struct vector3 *v1, struct vector3 *v2) {
  return max2d(max3d(fabs(v1->x), fabs(v1->y), fabs(v1->z)),
               max3d(fabs(v2->x), fabs(v2->y), fabs(v2->z)));
}

// create_new_poly_edge creates a new poly_edge and attaches it to the
// past pointer to linked list of poly_edges
static struct poly_edge* create_new_poly_edge(struct poly_edge* list);

// have_common_region checks if wall1 and wall2 located on the (same) object
// are part of a common region or not
static bool have_common_region(struct object *obj, int wall1, int wall2);


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
 *
 * create_new_poly_edge creates a new poly_edge and attaches it to the
 * past pointer to linked list of poly_edges
 *
 ***************************************************************************/
struct poly_edge* create_new_poly_edge(struct poly_edge* list) {

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
ehtable_add:
  In: pointer to an edge_hashtable struct
      pointer to the poly_edge to add
  Out: Returns 0 on success, 1 on failure.
       Edge is added to the hash table.
***************************************************************************/
int ehtable_add(struct edge_hashtable *eht, struct poly_edge *pe) {

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

/**************************************************************************\
 ** Edge construction section--builds permanent edges from hash table    **
\**************************************************************************/

/***************************************************************************
compatible_edges:
  In: array of pointers to walls
      index of first wall
      index of edge in first wall
      index of second wall
      index of edge in second wall
  Out: 1 if the edge joins the two walls
       0 if not (i.e. the wall doesn't contain the edge or the edge is
       traversed in the same direction in each or the two walls are
       actually the same wall)
***************************************************************************/
static int compatible_edges(struct wall **faces, int wA, int eA, int wB,
                            int eB) {
  struct vector3 *vA0, *vA1, *vA2, *vB0, *vB1, *vB2;

  if ((wA < 0) || (eA < 0) || (wB < 0) || (eB < 0))
    return 0;

  vA0 = faces[wA]->vert[eA];
  if (eA == 2)
    vA1 = faces[wA]->vert[0];
  else
    vA1 = faces[wA]->vert[eA + 1];
  if (eA == 0)
    vA2 = faces[wA]->vert[2];
  else
    vA2 = faces[wA]->vert[eA - 1];

  vB0 = faces[wB]->vert[eB];
  if (eB == 2)
    vB1 = faces[wB]->vert[0];
  else
    vB1 = faces[wB]->vert[eB + 1];
  if (eB == 0)
    vB2 = faces[wB]->vert[2];
  else
    vB2 = faces[wB]->vert[eB - 1];

  return ((vA0 == vB1 && vA1 == vB0 && vA2 != vB2) ||
          (vA0->x == vB1->x && vA0->y == vB1->y && vA0->z == vB1->z &&
           vA1->x == vB0->x && vA1->y == vB0->y && vA1->z == vB0->z &&
           !(vA2->x == vB2->x && vA2->y == vB2->y && vA2->z == vB2->z)));
}

/*****************************************************************************
 have_common_region checks if wall1 and wall2 located on the (same) object
 are part of a common region or not
******************************************************************************/
bool have_common_region(struct object *obj, int wall1, int wall2) {

  struct region_list *rl = obj->regions;
  bool common_region = false;
  while (rl != NULL) {
    struct region *r = rl->reg;
    if (strcmp(r->region_last_name, "ALL") == 0) {
      rl = rl->next;
      continue;
    }
    if (get_bit(r->membership, wall1) && get_bit(r->membership, wall2)) {
      common_region = true;
      break;
    }
    rl = rl->next;
  }
  return common_region;
}

/***************************************************************************
refine_edge_pairs:
  In: the head of a linked list of shared edges
      array of pointers to walls
  Out: No return value.  The best-matching pair of edges percolates up
       to be first in the list.  "Best-matching" means that the edge
       is traversed in different directions by each face, and that the
       normals of the two faces are as divergent as possible.
***************************************************************************/
static void refine_edge_pairs(struct poly_edge *p, struct wall **faces) {
#define TSWAP(x, y)                                                            \
  temp = (x);                                                                  \
  (x) = (y);                                                                   \
  (y) = temp

  int temp;

  double best_align = 2;
  bool share_region = false;
  struct poly_edge *best_p1 = p, *best_p2 = p;
  int best_n1 = 1;
  int best_n2 = 2;

  struct poly_edge *p1 = p;
  int n1 = 1;
  while (p1 != NULL && p1->n >= n1) {
    int wA, eA;
    if (n1 == 1) {
      wA = p1->face[0];
      eA = p1->edge[0];
    } else {
      wA = p1->face[1];
      eA = p1->edge[1];
    }

    struct poly_edge *p2;
    int n2;
    if (n1 == 1) {
      n2 = n1 + 1;
      p2 = p1;
    } else {
      n2 = 1;
      p2 = p1->next;
    }
    while (p2 != NULL && p2->n >= n2) {
      int wB, eB;
      if (n2 == 1) {
        wB = p2->face[0];
        eB = p2->edge[0];
      } else {
        wB = p2->face[1];
        eB = p2->edge[1];
      }

      // as soon as we hit an incompatible edge we can break out of the p2 loop
      // and continue scanning the next p1
      if (compatible_edges(faces, wA, eA, wB, eB)) {
        double align = faces[wA]->normal.x * faces[wB]->normal.x +
                       faces[wA]->normal.y * faces[wB]->normal.y +
                       faces[wA]->normal.z * faces[wB]->normal.z;

        // as soon as two walls have a common region we only consider walls who
        // share (any) region. We need to reset the best_align to make sure we
        // don't pick any wall that don't share a region discovered previously
        bool common_region = have_common_region(faces[wA]->parent_object, wA, wB);
        if (common_region) {
          if (!share_region) {
            best_align = 2;
          }
          share_region = true;
        }

        if (common_region || !share_region) {
          if (align < best_align) {
            best_p1 = p1;
            best_p2 = p2;
            best_n1 = n1;
            best_n2 = n2;
            best_align = align;
          }
        }
      } else {
        break;
      }

      if (n2 == 1)
        n2++;
      else {
        p2 = p2->next;
        n2 = 1;
      }
    }

    if (n1 == 1)
      n1++;
    else {
      p1 = p1->next;
      n1 = 1;
    }
  }

  /* swap best match into top spot */
  if (best_align > 1.0)
    return; /* No good pairs. */

  TSWAP(best_p1->face[best_n1-1], p->face[0]);
  TSWAP(best_p1->edge[best_n1-1], p->edge[0]);
  TSWAP(best_p2->face[best_n2-1], p->face[1]);
  TSWAP(best_p2->edge[best_n2-1], p->edge[1]);

#undef TSWAP
}

/***************************************************************************
surface_net:
  In: array of pointers to walls
      integer length of array
  Out: -1 if the surface is a manifold, 0 if it is not, 1 on malloc failure
       Walls end up connected across their edges.
  Note: Two edges must have their vertices listed in opposite order (i.e.
        connect two faces pointing the same way) to be linked.  If more than
        two faces share the same edge and can be linked, the faces with
        normals closest to each other will be linked.  We do not assume that
        the object is connected.  All pieces must be a manifold, however,
        for the entire object to be a manifold.  (That is, there must not
        be any free edges anywhere.)  It is possible to build weird, twisty
        self-intersecting things.  The behavior of these things during a
        simulation is not guaranteed to be well-defined.
***************************************************************************/
int surface_net(struct wall **facelist, int nfaces) {
  struct edge *e;
  int is_closed = 1;

  struct edge_hashtable eht;
  int nkeys = (3 * nfaces) / 2;
  if (ehtable_init(&eht, nkeys))
    return 1;

  for (int i = 0; i < nfaces; i++) {
    if (facelist[i] == NULL)
      continue;

    int k;
    int nedge = 3;
    for (int j = 0; j < nedge; j++) {

      if (j + 1 < nedge)
        k = j + 1;
      else
        k = 0;

      struct poly_edge pe;
      pe.v1x = facelist[i]->vert[j]->x;
      pe.v1y = facelist[i]->vert[j]->y;
      pe.v1z = facelist[i]->vert[j]->z;
      pe.v2x = facelist[i]->vert[k]->x;
      pe.v2y = facelist[i]->vert[k]->y;
      pe.v2z = facelist[i]->vert[k]->z;
      pe.face[0] = i;
      pe.edge[0] = j;

      if (ehtable_add(&eht, &pe))
        return 1;
    }
  }

  for (int i = 0; i < nkeys; i++) {
    struct poly_edge *pep = (eht.data + i);
    while (pep != NULL) {
      if (pep->n > 2) {
        refine_edge_pairs(pep, facelist);
      }
      if (pep->n >= 2) {
        if (pep->face[0] != -1 && pep->face[1] != -1) {
          if (compatible_edges(facelist, pep->face[0], pep->edge[0], pep->face[1],
                               pep->edge[1])) {
            facelist[pep->face[0]]->nb_walls[pep->edge[0]] = facelist[pep->face[1]];
            facelist[pep->face[1]]->nb_walls[pep->edge[1]] = facelist[pep->face[0]];
            e = (struct edge *)CHECKED_MEM_GET_NODIE(
                facelist[pep->face[0]]->birthplace->join, "edge");
            if (e == NULL)
              return 1;

            e->forward = facelist[pep->face[0]];
            e->backward = facelist[pep->face[1]];
            init_edge_transform(e, pep->edge[0]);
            facelist[pep->face[0]]->edges[pep->edge[0]] = e;
            facelist[pep->face[1]]->edges[pep->edge[1]] = e;
          }

        } else {
          is_closed = 0;
        }
      } else if (pep->n == 1) {
        is_closed = 0;
        e = (struct edge *)CHECKED_MEM_GET_NODIE(
            facelist[pep->face[0]]->birthplace->join, "edge");
        if (e == NULL)
          return 1;

        e->forward = facelist[pep->face[0]];
        e->backward = NULL;
        /* Don't call init_edge_transform unless both edges are set */
        facelist[pep->face[0]]->edges[pep->edge[0]] = e;
      }
      pep = pep->next;
    }
  }

  ehtable_kill(&eht);
  return -is_closed; /* We use 1 to indicate malloc failure so return 0/-1 */
}

/***************************************************************************
init_edge_transform
  In: e: pointer to an edge
      edgenum: integer telling which edge (0-2) of the "forward" face we are
  Out: No return value.  Coordinate transform in edge struct is set.
  Note: Don't call this on a non-shared edge.
***************************************************************************/

void init_edge_transform(struct edge *e, int edgenum) {
  struct vector2 O_f, O_b;
  struct vector2 ehat_f, ehat_b;
  struct vector2 fhat_f, fhat_b;

  struct wall *wf = e->forward;
  struct wall *wb = e->backward;
  int i = edgenum;
  int j = i + 1;
  if (j > 2)
    j = 0;

  /* Intermediate basis from the perspective of the forward frame */

  struct vector3 temp3d;
  temp3d.x = wf->vert[i]->x - wf->vert[0]->x;
  temp3d.y = wf->vert[i]->y - wf->vert[0]->y;
  temp3d.z = wf->vert[i]->z - wf->vert[0]->z;
  O_f.u = dot_prod(&temp3d, &(wf->unit_u));
  O_f.v = dot_prod(&temp3d, &(wf->unit_v)); /* Origin */

  temp3d.x = wf->vert[j]->x - wf->vert[0]->x;
  temp3d.y = wf->vert[j]->y - wf->vert[0]->y;
  temp3d.z = wf->vert[j]->z - wf->vert[0]->z;
  struct vector2 temp;
  temp.u = dot_prod(&temp3d, &(wf->unit_u)) - O_f.u;
  temp.v = dot_prod(&temp3d, &(wf->unit_v)) - O_f.v; /* Far side of e */

  double d = 1.0 / sqrt(temp.u * temp.u + temp.v * temp.v);
  ehat_f.u = temp.u * d;
  ehat_f.v = temp.v * d; /* ehat along edge */
  fhat_f.u = -ehat_f.v;
  fhat_f.v = ehat_f.u; /* fhat 90 degrees CCW */

  /* Intermediate basis from the perspective of the backward frame */

  temp3d.x = wf->vert[i]->x - wb->vert[0]->x;
  temp3d.y = wf->vert[i]->y - wb->vert[0]->y;
  temp3d.z = wf->vert[i]->z - wb->vert[0]->z;
  O_b.u = dot_prod(&temp3d, &(wb->unit_u));
  O_b.v = dot_prod(&temp3d, &(wb->unit_v)); /* Origin */

  temp3d.x = wf->vert[j]->x - wb->vert[0]->x;
  temp3d.y = wf->vert[j]->y - wb->vert[0]->y;
  temp3d.z = wf->vert[j]->z - wb->vert[0]->z;
  temp.u = dot_prod(&temp3d, &(wb->unit_u)) - O_b.u;
  temp.v = dot_prod(&temp3d, &(wb->unit_v)) - O_b.v; /* Far side of e */

  d = 1.0 / sqrt(temp.u * temp.u + temp.v * temp.v);
  ehat_b.u = temp.u * d;
  ehat_b.v = temp.v * d; /* ehat along edge */
  fhat_b.u = -ehat_b.v;
  fhat_b.v = ehat_b.u; /* fhat 90 degrees CCW */

  /* Calculate transformation matrix */

  double mtx[2][2];
  mtx[0][0] = ehat_f.u * ehat_b.u + fhat_f.u * fhat_b.u;
  mtx[0][1] = ehat_f.v * ehat_b.u + fhat_f.v * fhat_b.u;
  mtx[1][0] = ehat_f.u * ehat_b.v + fhat_f.u * fhat_b.v;
  mtx[1][1] = ehat_f.v * ehat_b.v + fhat_f.v * fhat_b.v;

  /* Calculate translation vector */

  struct vector2 q;
  q.u = O_b.u;
  q.v = O_b.v;

  q.u -= mtx[0][0] * O_f.u + mtx[0][1] * O_f.v;
  q.v -= mtx[1][0] * O_f.u + mtx[1][1] * O_f.v;

  /* Store the results */

  e->cos_theta = mtx[0][0];
  e->sin_theta = mtx[0][1];
  e->translate.u = q.u;
  e->translate.v = q.v;
}

/***************************************************************************
sharpen_object:
  In: parent: pointer to an object
  Out: 0 on success, 1 on failure.
       Adds edges to the object and all its children.
***************************************************************************/

int sharpen_object(struct object *parent) {
  if (parent->object_type == POLY_OBJ || parent->object_type == BOX_OBJ) {
    int i = surface_net(parent->wall_p, parent->n_walls);

    if (i == 1) {
      mcell_allocfailed(
          "Failed to connect walls of object %s along shared edges.",
          parent->sym->name);
    }
    else {
      parent->is_closed = -i; 
    }
  } else if (parent->object_type == META_OBJ) {
    for (struct object *o = parent->first_child; o != NULL; o = o->next) {
      if (sharpen_object(o))
        return 1;
    }
  }

  return 0;
}

/***************************************************************************
sharpen_world:
  In: nothing.  Assumes if there are polygon objects then they have been
      initialized and placed in the world in their correct memory locations.
  Out: 0 on success, 1 on failure.  Adds edges to every object.
***************************************************************************/
int sharpen_world(struct volume *world) {
  for (struct object *o = world->root_instance; o != NULL; o = o->next) {
    if (sharpen_object(o))
      return 1;
  }
  return 0;
}

/**************************************************************************\
 ** Geometry section--report on and use geometrical properties of object **
\**************************************************************************/

/***************************************************************************
closest_interior_point:
  In: a point in 3D
      a wall
      the surface coordinates of the closest interior point on the wall
      how far away the point can be before we give up
  Out: return the distance^2 between the input point and closest point.
       Sets closest interior point.
  Note: the search distance currently isn't used.  This function is just
        a wrapper for closest_pt_point_triangle.  If the closest point is
        on an edge or corner, we scoot the point towards the centroid of
        the triangle so we're contained fully within the triangle.
***************************************************************************/

double closest_interior_point(struct vector3 *pt, struct wall *w,
                              struct vector2 *ip, double r2) {
  UNUSED(r2);

  struct vector3 v;

  closest_pt_point_triangle(pt, w->vert[0], w->vert[1], w->vert[2], &v);
  xyz2uv(&v, w, ip);

  /* Check to see if we're lying on an edge; if so, scoot towards centroid. */
  /* ip lies on edge of wall if cross products are zero */

  int give_up_ctr = 0;
  int give_up = 10;
  double a1 = ip->u * w->uv_vert2.v - ip->v * w->uv_vert2.u;
  double a2 = w->uv_vert1_u * ip->v;
  struct vector2 vert_0 = {.u = 0, .v = 0};
  struct vector2 vert_1 = {.u = w->uv_vert1_u, .v = 0};
  while (give_up_ctr < give_up &&
         (!distinguishable(ip->v, 0, EPS_C) ||
          !distinguishable(a1, 0, EPS_C) ||
          !distinguishable(a1 + a2, 2.0 * w->area, EPS_C) ||
          !point_in_triangle_2D(ip, &vert_0, &vert_1, &w->uv_vert2))) {
    /* Move toward centroid. It's possible for this movement to be so small
     * that we are essentially stuck in this loop, so bail out after a set
     * number of tries. The number chosen is somewhat arbitrary. In most cases,
     * one try is sufficent. */
    ip->u = (1.0 - 5 * EPS_C) * ip->u +
            5 * EPS_C * 0.333333333333333 * (w->uv_vert1_u + w->uv_vert2.u);
    ip->v = (1.0 - 5 * EPS_C) * ip->v +
            5 * EPS_C * 0.333333333333333 * w->uv_vert2.v;

    a1 = ip->u * w->uv_vert2.v - ip->v * w->uv_vert2.u;
    a2 = w->uv_vert1_u * ip->v;

    give_up_ctr++;
  }
  return (v.x - pt->x) * (v.x - pt->x) + (v.y - pt->y) * (v.y - pt->y) +
         (v.z - pt->z) * (v.z - pt->z);
}

/***************************************************************************
find_edge_point:
  In: here: a wall
      loc: a point in the coordinate system of that wall where we are now
           (assumed to be on or inside triangle)
      disp: a 2D displacement vector to move
      edgept: a place to store the coordinate on the edge, if we hit it
  Out: index of the edge we hit (0, 1, or 2), or -1 if the new location
       is within the wall, or -2 if we can't tell.  If the result is
       0, 1, or 2, edgept is set to the new location.
***************************************************************************/
int find_edge_point(struct wall *here,
                    struct vector2 *loc,
                    struct vector2 *disp,
                    struct vector2 *edgept) {
  double f, s, t;

  double lxd = loc->u * disp->v - loc->v * disp->u;

  double lxc1 = -loc->v * here->uv_vert1_u;
  double dxc1 = -disp->v * here->uv_vert1_u;

  // Make sure that the displacement vector isn't on v0v1
  if (dxc1 < -EPS_C || dxc1 > EPS_C) {
    f = 1.0 / dxc1; /* f>0 is passing outwards */
    s = -lxd * f;
    if (0.0 < s && s < 1.0 && f > 0.0) {
      t = -lxc1 * f;
      if (EPS_C < t && t < 1.0) {
        edgept->u = loc->u + t * disp->u;
        edgept->v = loc->v + t * disp->v;
        return 0;
      } else if (t > 1.0 + EPS_C)
        return -1;
      /* else can't tell if we hit this edge, assume not */
    }
  }

  double lxc2 = loc->u * here->uv_vert2.v - loc->v * here->uv_vert2.u;
  double dxc2 = disp->u * here->uv_vert2.v - disp->v * here->uv_vert2.u;

  // Make sure that the displacement vector isn't on v1v2
  if (dxc2 < -EPS_C || dxc2 > EPS_C) {
    f = 1.0 / dxc2; /* f<0 is passing outwards */
    s = 1.0 + lxd * f;
    if (0.0 < s && s < 1.0 && f < 0.0) {
      t = -lxc2 * f;
      if (EPS_C < t && t < 1.0) {
        edgept->u = loc->u + t * disp->u;
        edgept->v = loc->v + t * disp->v;
        return 2;
      } else if (t > 1.0 + EPS_C)
        return -1;
      /* else can't tell */
    }
  }

  f = dxc2 - dxc1;

  if (f < -EPS_C || f > EPS_C) {
    f = 1.0 / f; /* f>0 is passing outwards */
    s = -(lxd + dxc1) * f;
    if (0.0 < s && s < 1.0 && f > 0.0) {
      t = (here->uv_vert1_u * here->uv_vert2.v + lxc1 - lxc2) * f;
      if (EPS_C < t && t < 1.0) {
        edgept->u = loc->u + t * disp->u;
        edgept->v = loc->v + t * disp->v;
        return 1;
      } else if (t > 1.0 + EPS_C)
        return -1;
      /* else can't tell */
    }
  }

  return -2; /* Couldn't tell whether we hit or not--calling function should
                pick another displacement */
}

/***************************************************************************
traverse_surface:
  In: here: a wall
      loc: a point in the coordinate system of that wall
      which: which edge to travel off of
      newloc: a vector to set for the new wall
  Out: NULL if the edge is not shared, or a pointer to the wall in that
       direction if it is shared. newloc is set to loc in the coordinate system
       of the new wall (after flattening the walls along their shared edge)
***************************************************************************/
struct wall *traverse_surface(struct wall *here, struct vector2 *loc, int which,
  struct vector2 *newloc) {

  struct wall *there;
  double u, v;

  struct edge *e = here->edges[which];

  if (e == NULL)
    return NULL;

  if (e->forward == here) {
    /* Apply forward transform to loc */
    there = e->backward;

    /* rotation */
    u = e->cos_theta * loc->u + e->sin_theta * loc->v;
    v = -e->sin_theta * loc->u + e->cos_theta * loc->v;

    /* translation */
    newloc->u = u + e->translate.u;
    newloc->v = v + e->translate.v;

    return there;
  } else {
    /* Apply inverse transform to loc */
    there = e->forward;

    /* inverse translation */
    u = loc->u - e->translate.u;
    v = loc->v - e->translate.v;

    /* inverse rotation */
    newloc->u = e->cos_theta * u - e->sin_theta * v;
    newloc->v = e->sin_theta * u + e->cos_theta * v;

    return there;
  }
}

/***************************************************************************
is_manifold:
  In: r: A region. This region must already be painted on walls. The edges must
         have already been added to the object (i.e. sharpened).
      count_regions_flag: This is usually set, unless we are only checking
                          volumes for dynamic geometries.
  Out: 1 if the region is a manifold, 0 otherwise.
  Note: by "manifold" we mean "orientable compact two-dimensional
        manifold without boundaries embedded in R3"
***************************************************************************/
int is_manifold(struct region *r, int count_regions_flag) {
  struct wall **wall_array = NULL, *w = NULL;
  struct region_list *rl = NULL;

  if (r->bbox == NULL)
    r->bbox = create_region_bbox(r);

  wall_array = r->parent->wall_p;

  if (wall_array == NULL) {
    mcell_internal_error("Region '%s' has NULL wall array!", r->sym->name);
    /*return 1;*/
  }

  // use the center of the region bounding box as reference point for
  // computing the volume
  struct vector3 llc = r->bbox[0];
  struct vector3 urc = r->bbox[1];
  struct vector3 d;
  d.x = (llc.x + 0.5 * urc.x);
  d.y = (llc.y + 0.5 * urc.y);
  d.z = (llc.z + 0.5 * urc.z);

  r->volume = 0.0;
  for (int n_wall = 0; n_wall < r->parent->n_walls; n_wall++) {
    if (!get_bit(r->membership, n_wall))
      continue; /* Skip wall not in region */
    w = wall_array[n_wall];
    if (count_regions_flag) {
      for (int nb = 0; nb < 3; nb++) {
        if (w->nb_walls[nb] == NULL) {
          mcell_error_nodie("BARE EDGE on wall %u edge %d.", n_wall, nb);
          return 0; /* Bare edge--not a manifold */
        }

        for (rl = w->nb_walls[nb]->counting_regions; rl != NULL; rl = rl->next) {
          if (rl->reg == r)
            break;
        }
        if (rl == NULL) {
          mcell_error_nodie("Wall %u edge %d leaves region!", n_wall, nb);
          return 0; /* Can leave region--not a manifold */
        }
      }
    }
    // compute volume of tetrahedron with w as its face
    r->volume += tetrahedralVol(w->vert[0], w->vert[1], w->vert[2], &d);
  }
  return 1;
}

/**************************************************************************\
 ** Collision section--detect whether rays intersect walls               **
\**************************************************************************/

/***************************************************************************
jump_away_line:
  In: starting coordinate
      vector we were going to move along and need to change
      fraction of way we moved before noticing we were hitting a edge
      location of the first vertex of the edge
      location of the second vertex of the edge
      normal vector to the surface containing our edge
  Out: No return value.  Movement vector is slightly changed.
***************************************************************************/
void jump_away_line(struct vector3 *p, struct vector3 *v, double k,
                    struct vector3 *A, struct vector3 *B, struct vector3 *n,
                    struct rng_state *rng) {
  struct vector3 e, f;
  double le_1, tiny;

  e.x = B->x - A->x;
  e.y = B->y - A->y;
  e.z = B->z - A->z;

  le_1 = 1.0 / sqrt(e.x * e.x + e.y * e.y + e.z * e.z);

  e.x *= le_1;
  e.y *= le_1;
  e.z *= le_1;

  f.x = n->y * e.z - n->z * e.y;
  f.y = n->z * e.x - n->x * e.z;
  f.z = n->x * e.y - n->y * e.x;

  tiny = EPS_C * (abs_max_2vec(p, v) + 1.0) /
         (k * max3d(fabs(f.x), fabs(f.y), fabs(f.z)));
  if ((rng_uint(rng) & 1) == 0) {
    tiny = -tiny;
  }
  v->x -= tiny * f.x;
  v->y -= tiny * f.y;
  v->z -= tiny * f.z;
}

/***************************************************************************
collide_wall:
  In: point: starting coordinate
      move: vector to move along
      face: wall we're checking for a collision
      t: double to store time of collision
      hitpt: vector to store the location of the collision
      update_move: flag to signal whether we should modify the movement vector
        in an ambiguous case (i.e. if we hit an edge or corner); if not, any
        ambiguous cases are treated as a miss.
  Out: Integer value indicating what happened
         COLLIDE_MISS  missed
         COLLIDE_FRONT hit the front face (face normal points out of)
         COLLIDE_BACK  hit the back face
         COLLIDE_REDO  hit an edge and modified movement vector; redo
  Note: t and/or hitpt may be modified even if there is no collision
        Not highly optimized yet.  May want to project to Cartesian
        coordinates for speed (as MCell2 did, and Rex implemented
        in pre-40308 backups in vol_utils.c).  When reflecting, use
        the value of t returned, not hitpt (reflections happen slightly
        early to avoid rounding errors with superimposed planes).
***************************************************************************/
int collide_wall(struct vector3 *point, struct vector3 *move, struct wall *face,
                 double *t, struct vector3 *hitpt, int update_move,
                 struct rng_state *rng, struct notifications *notify,
                 long long *ray_polygon_tests) {
  double dp, dv, dd;
  double nx, ny, nz;
  double a, b, c;
  double f, g, h;
  double d_eps;
  struct vector3 local;

  if (notify->final_summary == NOTIFY_FULL) {
    (*ray_polygon_tests)++;
  }

  nx = face->normal.x;
  ny = face->normal.y;
  nz = face->normal.z;

  dp = nx * point->x + ny * point->y + nz * point->z;
  dv = nx * move->x + ny * move->y + nz * move->z;
  dd = dp - face->d;

  if (dd > 0.0) {
    d_eps = EPS_C;
    if (dd < d_eps)
      d_eps = 0.5 * dd;

    /* Start & end above plane */
    if (dd + dv > d_eps)
      return COLLIDE_MISS;
  } else {
    d_eps = -EPS_C;
    if (dd > d_eps)
      d_eps = 0.5 * dd;

    /* Start & end below plane */
    if (dd < 0.0 && dd + dv < d_eps)
      return COLLIDE_MISS;
  }

  if (dd == 0.0) {
    /* Start beside plane, end above or below */
    if (dv != 0.0)
      return COLLIDE_MISS;

    if (update_move) {
      a = (abs_max_2vec(point, move) + 1.0) * EPS_C;
      if ((rng_uint(rng) & 1) == 0)
        a = -a;
      if (dd == 0.0) {
        move->x -= a * nx;
        move->y -= a * ny;
        move->z -= a * nz;
      } else {
        move->x *= (1.0 - a);
        move->y *= (1.0 - a);
        move->z *= (1.0 - a);
      }
      return COLLIDE_REDO;
    } else
      return COLLIDE_MISS;
  }

  a = 1.0 / dv;
  a *= -dd; /* Time we actually hit */
  *t = a;

  hitpt->x = point->x + a * move->x;
  hitpt->y = point->y + a * move->y;
  hitpt->z = point->z + a * move->z;

  local.x = hitpt->x - face->vert[0]->x;
  local.y = hitpt->y - face->vert[0]->y;
  local.z = hitpt->z - face->vert[0]->z;

  b = local.x * face->unit_u.x + local.y * face->unit_u.y +
      local.z * face->unit_u.z;
  c = local.x * face->unit_v.x + local.y * face->unit_v.y +
      local.z * face->unit_v.z;

  if (face->uv_vert2.v < 0.0) {
    c = -c;
    f = -face->uv_vert2.v;
  } else
    f = face->uv_vert2.v;

  if (c > 0) {
    g = b * f;
    h = c * face->uv_vert2.u;
    if (g > h) {
      if (c * face->uv_vert1_u + g < h + face->uv_vert1_u * face->uv_vert2.v) {
        if (dv > 0)
          return COLLIDE_BACK;
        else
          return COLLIDE_FRONT;
      } else if ((!distinguishable(
          c * face->uv_vert1_u + g,
          h + face->uv_vert1_u * face->uv_vert2.v,
          EPS_C))) {
        if (update_move) {
          jump_away_line(point, move, a, face->vert[1], face->vert[2],
                         &(face->normal), rng);
          return COLLIDE_REDO;
        } else
          return COLLIDE_MISS;
      } else
        return COLLIDE_MISS;
    } else if (!distinguishable(g, h, EPS_C)) {
      if (update_move) {
        jump_away_line(point, move, a, face->vert[2], face->vert[0],
                       &(face->normal), rng);
        return COLLIDE_REDO;
      } else
        return COLLIDE_MISS;
    } else
      return COLLIDE_MISS;
  } else if (!distinguishable(c, 0, EPS_C)) /* Hit first edge! */
  {
    if (update_move) {
      jump_away_line(point, move, a, face->vert[0], face->vert[1],
                     &(face->normal), rng);
      return COLLIDE_REDO;
    } else
      return COLLIDE_MISS;
  } else
    return COLLIDE_MISS;
}

/***************************************************************************
collide_mol:
  In: starting coordinate
      vector to move along
      molecule we're checking for a collision
      double to store time of collision
      vector to store the location of the collision
  Out: Integer value indicating what happened
         COLLIDE_MISS   missed
         COLLIDE_VOL_M  hit
  Note: t and/or hitpt may be modified even if there is no collision
        Not highly optimized yet.
***************************************************************************/
int collide_mol(struct vector3 *point, struct vector3 *move,
                struct abstract_molecule *a, double *t, struct vector3 *hitpt,
                double rx_radius_3d) {
  struct vector3 dir;  /* From starting point of moving molecule to target */
  struct vector3 *pos; /* Position of target molecule */

  double movelen2; /* Square of distance the moving molecule travels */
  double dirlen2;  /* Square of distance between moving and target molecules */
  double d;        /* Dot product of movement vector and vector to target */
  double sigma2;   /* Square of interaction radius */

  if ((a->properties->flags & ON_GRID) != 0)
    return COLLIDE_MISS; /* Should never call on surface molecule! */

  pos = &(((struct volume_molecule *)a)->pos);

  sigma2 = rx_radius_3d * rx_radius_3d;

  dir.x = pos->x - point->x;
  dir.y = pos->y - point->y;
  dir.z = pos->z - point->z;

  d = dir.x * move->x + dir.y * move->y + dir.z * move->z;

  /* Miss the molecule if it's behind us */
  if (d < 0)
    return COLLIDE_MISS;

  movelen2 = move->x * move->x + move->y * move->y + move->z * move->z;

  /* check whether the test molecule is further than the displacement. */
  if (d > movelen2)
    return COLLIDE_MISS;

  dirlen2 = dir.x * dir.x + dir.y * dir.y + dir.z * dir.z;

  /* check whether the moving molecule will miss interaction disk of the
     test molecule.*/
  if (movelen2 * dirlen2 - d * d > movelen2 * sigma2)
    return COLLIDE_MISS;

  *t = d / movelen2;
  //  *t = d/sqrt(movelen2*dirlen2);
  hitpt->x = point->x + (*t) * move->x;
  hitpt->y = point->y + (*t) * move->y;
  hitpt->z = point->z + (*t) * move->z;
  return COLLIDE_VOL_M;
}

/***************************************************************************
wall_in_box:
  In: array of pointers to vertices for wall (should be 3)
      normal vector for wall
      distance from wall to origin (point normal form)
      first corner of bounding box
      opposite corner of bounding box
  Out: 1 if the wall intersects the box.  0 otherwise.
***************************************************************************/
static int wall_in_box(struct vector3 **vert, struct vector3 *normal, double d,
                       struct vector3 *b0, struct vector3 *b1) {
#define n_vert 3
  int temp;
  int i, j, k;
  struct vector3 *v1, *v2;
  struct vector3 n, u, v;
  struct vector3 ba, bb, c;
  double r, a1, a2, a3, a4, cu, cv;
  double vu_[6]; /* Assume wall has 3 vertices */
  double *vv_;
  double d_box[8];
  int n_opposite;

  /* Lookup table for vertex-edge mapping for a cube */
  int which_x1[12] = { 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 1 };
  int which_y1[12] = { 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 0 };
  int which_z1[12] = { 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0 };
  int which_x2[12] = { 0, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 0 };
  int which_y2[12] = { 0, 1, 1, 1, 1, 0, 0, 0, 1, 0, 1, 1 };
  int which_z2[12] = { 1, 1, 0, 0, 1, 1, 0, 0, 0, 1, 1, 0 };

  int edge1_vt[12] = { 0, 1, 3, 2, 6, 7, 5, 4, 0, 1, 3, 4 };
  int edge2_vt[12] = { 1, 3, 2, 6, 7, 5, 4, 0, 2, 5, 7, 2 };

  /* Check if any vertex of the wall is in the box. */
  for (i = 0; i < n_vert; i++) {
    v2 = vert[i];
    if (v2->x >= b0->x && v2->x <= b1->x && v2->y >= b0->y && v2->y <= b1->y &&
        v2->z >= b0->z && v2->z <= b1->z)
      return 1;
  }

  /* Check if any wall edge intersects any face of the box */
  for (i = 0; i < n_vert; i++) {
    v2 = vert[i];
    v1 = (i == 0) ? vert[n_vert - 1] : vert[i - 1];

    /* x-faces */
    if ((v1->x <= b0->x && b0->x < v2->x) ||
        (v1->x > b0->x && b0->x >= v2->x)) {
      r = (b0->x - v1->x) / (v2->x - v1->x);
      a3 = v1->y + r * (v2->y - v1->y);
      a4 = v1->z + r * (v2->z - v1->z);
      if (b0->y <= a3 && a3 <= b1->y && b0->z <= a4 && a4 <= b1->z)
        return 2;
    }
    if ((v1->x <= b1->x && b1->x < v2->x) ||
        (v1->x > b1->x && b1->x >= v2->x)) {
      r = (b1->x - v1->x) / (v2->x - v1->x);
      a3 = v1->y + r * (v2->y - v1->y);
      a4 = v1->z + r * (v2->z - v1->z);
      if (b0->y <= a3 && a3 <= b1->y && b0->z <= a4 && a4 <= b1->z)
        return 3;
    }

    /* y-faces */
    if ((v1->y <= b0->y && b0->y < v2->y) ||
        (v1->y > b0->y && b0->y >= v2->y)) {
      r = (b0->y - v1->y) / (v2->y - v1->y);
      a3 = v1->x + r * (v2->x - v1->x);
      a4 = v1->z + r * (v2->z - v1->z);
      if (b0->x <= a3 && a3 <= b1->x && b0->z <= a4 && a4 <= b1->z)
        return 4;
    }
    if ((v1->y <= b1->y && b1->y < v2->y) ||
        (v1->y > b1->y && b1->y >= v2->y)) {
      r = (b1->y - v1->y) / (v2->y - v1->y);
      a3 = v1->x + r * (v2->x - v1->x);
      a4 = v1->z + r * (v2->z - v1->z);
      if (b0->x <= a3 && a3 <= b1->x && b0->z <= a4 && a4 <= b1->z)
        return 5;
    }

    /* z-faces */
    if ((v1->z <= b0->z && b0->z < v2->z) ||
        (v1->z > b0->z && b0->z >= v2->z)) {
      r = (b0->z - v1->z) / (v2->z - v1->z);
      a3 = v1->y + r * (v2->y - v1->y);
      a4 = v1->x + r * (v2->x - v1->x);
      if (b0->y <= a3 && a3 <= b1->y && b0->x <= a4 && a4 <= b1->x)
        return 6;
    }
    if ((v1->z <= b1->z && b1->z < v2->z) ||
        (v1->z > b1->z && b1->z >= v2->z)) {
      r = (b1->z - v1->z) / (v2->z - v1->z);
      a3 = v1->y + r * (v2->y - v1->y);
      a4 = v1->x + r * (v2->x - v1->x);
      if (b0->y <= a3 && a3 <= b1->y && b0->x <= a4 && a4 <= b1->x)
        return 7;
    }
  }

  /* Check if any box edge intersects the wall */

  n_opposite = 0;
  vv_ = &(vu_[n_vert]);

  /* Wall coordinate system n,u,v */
  n.x = normal->x;
  n.y = normal->y;
  n.z = normal->z;
  u.x = vert[1]->x - vert[0]->x;
  u.y = vert[1]->y - vert[0]->y;
  u.z = vert[1]->z - vert[0]->z;
  r = 1 / sqrt(u.x * u.x + u.y * u.y + u.z * u.z);
  u.x *= r;
  u.y *= r;
  u.z *= r;
  v.x = n.y * u.z - n.z * u.y;
  v.y = -(n.x * u.z - n.z * u.x);
  v.z = n.x * u.y - n.y * u.x;
  for (j = 0; j < n_vert; j++) {
    vu_[j] = vert[j]->x * u.x + vert[j]->y * u.y + vert[j]->z * u.z;
    vv_[j] = vert[j]->x * v.x + vert[j]->y * v.y + vert[j]->z * v.z;
  }

  /* Test every edge. */
  bb.x = b0->x;
  bb.y = b0->y;
  bb.z = b0->z;
  d_box[0] = bb.x * n.x + bb.y * n.y + bb.z * n.z;
  for (i = 0; i < 12; i++) {
    if (i < 7) /* Visiting new vertices in order */
    {
      ba.x = bb.x;
      ba.y = bb.y;
      ba.z = bb.z;
      bb.x = (which_x2[i]) ? b1->x : b0->x;
      bb.y = (which_y2[i]) ? b1->y : b0->y;
      bb.z = (which_z2[i]) ? b1->z : b0->z;
      a2 = d_box[edge2_vt[i]] = bb.x * n.x + bb.y * n.y + bb.z * n.z;
      a1 = d_box[edge1_vt[i]];

      if ((a1 - d < 0 && a2 - d < 0) || (a1 - d > 0 && a2 - d > 0))
        continue;
      else
        n_opposite++;
    } else /* Revisiting old vertices out of order */
    {
      /*      if (!n_opposite) return 0; */
      a1 = d_box[edge1_vt[i]];
      a2 = d_box[edge2_vt[i]];

      if ((a1 - d < 0 && a2 - d < 0) || (a1 - d > 0 && a2 - d > 0))
        continue;

      n_opposite++;
      ba.x = (which_x1[i]) ? b1->x : b0->x;
      ba.y = (which_y1[i]) ? b1->y : b0->y;
      ba.z = (which_z1[i]) ? b1->z : b0->z;
      bb.x = (which_x2[i]) ? b1->x : b0->x;
      bb.y = (which_y2[i]) ? b1->y : b0->y;
      bb.z = (which_z2[i]) ? b1->z : b0->z;
    }
    /* Now ba,bb = box edge endpoints ; a1,a2 = distances along wall normal */
    r = (d - a1) / (a2 - a1);
    c.x = ba.x + r * (bb.x - ba.x);
    c.y = ba.y + r * (bb.y - ba.y);
    c.z = ba.z + r * (bb.z - ba.z);
    cu = c.x * u.x + c.y * u.y + c.z * u.z;
    cv = c.x * v.x + c.y * v.y + c.z * v.z;
    /* Test for internal intersection point in wall coordinate space */
    temp = 0;
    for (j = 0; j < n_vert; j++) {
      k = (j == 0) ? n_vert - 1 : j - 1;
      if ((vu_[k] < cu && cu <= vu_[j]) || (vu_[k] >= cu && cu > vu_[j])) {
        r = (cu - vu_[k]) / (vu_[j] - vu_[k]);
        if ((vv_[k] + r * (vv_[j] - vv_[k])) > cv)
          temp++;
      }
    }
    if (temp & 1)
      return 8 + i;
  }

  return 0;
#undef n_vert
}

/***************************************************************************
init_tri_wall:
  In: object to which the wall belongs
      index of the wall within that object
      three vectors defining the vertices of the wall.
  Out: No return value.  The wall is properly initialized with normal
       vectors, local coordinate vectors, and so on.
***************************************************************************/

void init_tri_wall(struct object *objp, int side, struct vector3 *v0,
                   struct vector3 *v1, struct vector3 *v2) {
  struct wall *w; /* The wall we're working with */
  double f, fx, fy, fz;
  struct vector3 vA, vB, vX;

  w = &objp->walls[side];
  w->next = NULL;
  w->surf_class_head = NULL;
  w->num_surf_classes = 0;

  w->side = side;

  w->vert[0] = v0;
  w->vert[1] = v1;
  w->vert[2] = v2;

  w->edges[0] = NULL;
  w->edges[1] = NULL;
  w->edges[2] = NULL;
  w->nb_walls[0] = NULL;
  w->nb_walls[1] = NULL;
  w->nb_walls[2] = NULL;

  vectorize(v0, v1, &vA);
  vectorize(v0, v2, &vB);
  cross_prod(&vA, &vB, &vX);
  w->area = 0.5 * vect_length(&vX);

  if (!distinguishable(w->area, 0, EPS_C)) {
    /* this is a degenerate polygon.
    * perform initialization and quit. */
    w->unit_u.x = 0;
    w->unit_u.y = 0;
    w->unit_u.z = 0;

    w->normal.x = 0;
    w->normal.y = 0;
    w->normal.z = 0;

    w->unit_v.x = 0;
    w->unit_v.y = 0;
    w->unit_v.z = 0;
    w->d = 0;
    w->uv_vert1_u = 0;
    w->uv_vert2.u = 0;
    w->uv_vert2.v = 0;

    w->grid = NULL;

    w->parent_object = objp;
    w->flags = 0;
    w->counting_regions = NULL;

    return;
  }

  fx = (v1->x - v0->x);
  fy = (v1->y - v0->y);
  fz = (v1->z - v0->z);
  f = 1 / sqrt(fx * fx + fy * fy + fz * fz);

  w->unit_u.x = fx * f;
  w->unit_u.y = fy * f;
  w->unit_u.z = fz * f;

  fx = (v2->x - v0->x);
  fy = (v2->y - v0->y);
  fz = (v2->z - v0->z);

  w->normal.x = w->unit_u.y * fz - w->unit_u.z * fy;
  w->normal.y = w->unit_u.z * fx - w->unit_u.x * fz;
  w->normal.z = w->unit_u.x * fy - w->unit_u.y * fx;
  f = 1 / sqrt(w->normal.x * w->normal.x + w->normal.y * w->normal.y +
               w->normal.z * w->normal.z);
  w->normal.x *= f;
  w->normal.y *= f;
  w->normal.z *= f;
  w->unit_v.x = w->normal.y * w->unit_u.z - w->normal.z * w->unit_u.y;
  w->unit_v.y = w->normal.z * w->unit_u.x - w->normal.x * w->unit_u.z;
  w->unit_v.z = w->normal.x * w->unit_u.y - w->normal.y * w->unit_u.x;
  w->d = v0->x * w->normal.x + v0->y * w->normal.y + v0->z * w->normal.z;

  w->uv_vert1_u = (w->vert[1]->x - w->vert[0]->x) * w->unit_u.x +
                  (w->vert[1]->y - w->vert[0]->y) * w->unit_u.y +
                  (w->vert[1]->z - w->vert[0]->z) * w->unit_u.z;
  w->uv_vert2.u = (w->vert[2]->x - w->vert[0]->x) * w->unit_u.x +
                  (w->vert[2]->y - w->vert[0]->y) * w->unit_u.y +
                  (w->vert[2]->z - w->vert[0]->z) * w->unit_u.z;
  w->uv_vert2.v = (w->vert[2]->x - w->vert[0]->x) * w->unit_v.x +
                  (w->vert[2]->y - w->vert[0]->y) * w->unit_v.y +
                  (w->vert[2]->z - w->vert[0]->z) * w->unit_v.z;

  w->grid = NULL;

  w->parent_object = objp;
  w->flags = 0;
  w->counting_regions = NULL;
}

/***************************************************************************
wall_bounding_box:
  In: a wall
      vector to store one corner of the bounding box for that wall
      vector to store the opposite corner
  Out: No return value.  The vectors are set to define the smallest box
       that contains the wall.
***************************************************************************/
static void wall_bounding_box(struct wall *w, struct vector3 *llf,
                              struct vector3 *urb) {
  llf->x = urb->x = w->vert[0]->x;
  llf->y = urb->y = w->vert[0]->y;
  llf->z = urb->z = w->vert[0]->z;

  if (w->vert[1]->x < llf->x)
    llf->x = w->vert[1]->x;
  else if (w->vert[1]->x > urb->x)
    urb->x = w->vert[1]->x;
  if (w->vert[2]->x < llf->x)
    llf->x = w->vert[2]->x;
  else if (w->vert[2]->x > urb->x)
    urb->x = w->vert[2]->x;

  if (w->vert[1]->y < llf->y)
    llf->y = w->vert[1]->y;
  else if (w->vert[1]->y > urb->y)
    urb->y = w->vert[1]->y;
  if (w->vert[2]->y < llf->y)
    llf->y = w->vert[2]->y;
  else if (w->vert[2]->y > urb->y)
    urb->y = w->vert[2]->y;

  if (w->vert[1]->z < llf->z)
    llf->z = w->vert[1]->z;
  else if (w->vert[1]->z > urb->z)
    urb->z = w->vert[1]->z;
  if (w->vert[2]->z < llf->z)
    llf->z = w->vert[2]->z;
  else if (w->vert[2]->z > urb->z)
    urb->z = w->vert[2]->z;
}

/***************************************************************************
wall_to_vol:
  In: a wall
      the subvolume to which the wall belongs
  Out: The updated list of walls for that subvolume that now contains the
       wall requested.
***************************************************************************/
struct wall_list *wall_to_vol(struct wall *w, struct subvolume *sv) {
  struct wall_list *wl =
      CHECKED_MEM_GET_NODIE(sv->local_storage->list, "wall list");
  if (wl == NULL)
    return NULL;

  wl->this_wall = w;
  wl->next = sv->wall_head;
  sv->wall_head = wl;

  return wl;
}

/***************************************************************************
localize_wall:
  In: a wall
      the local memory storage area where this wall should be stored
  Out: A pointer to the copy of that wall in local memory, or NULL on
       memory allocation failure.
***************************************************************************/
struct wall *localize_wall(struct wall *w, struct storage *stor) {
  struct wall *ww;
  ww = CHECKED_MEM_GET_NODIE(stor->face, "wall");
  if (ww == NULL)
    return NULL;

  memcpy(ww, w, sizeof(struct wall));
  ww->next = stor->wall_head;
  stor->wall_head = ww;
  stor->wall_count++;
  ww->birthplace = stor;

  return ww;
}

/***************************************************************************
distribute_wall:
  In: a wall belonging to an object
  Out: A pointer to the wall as copied into appropriate local memory, or
       NULL on memory allocation error.  Also, the wall is added to the
       appropriate wall lists for all subvolumes it intersects; if this
       fails due to memory allocation errors, NULL is also returned.
***************************************************************************/
static struct wall *distribute_wall(struct volume *world, struct wall *w) {
  struct wall *where_am_i;       /* Version of the wall in local memory */
  struct vector3 llf, urb, cent; /* Bounding box for wall */
  int x_max, x_min, y_max, y_min, z_max,
      z_min;           /* Enlarged box to avoid rounding */
  int h, i, j, k;      /* Iteration variables for subvolumes */
  double leeway = 1.0; /* Margin of error */

  wall_bounding_box(w, &llf, &urb);

  if (llf.x < -leeway)
    leeway = -llf.x;
  if (llf.y < -leeway)
    leeway = -llf.y;
  if (llf.z < -leeway)
    leeway = -llf.z;
  if (urb.x > leeway)
    leeway = urb.x;
  if (urb.y > leeway)
    leeway = urb.y;
  if (urb.z > leeway)
    leeway = urb.z;
  leeway = EPS_C + leeway * EPS_C;
  if (world->use_expanded_list) {
    leeway += world->rx_radius_3d;
  }

  llf.x -= leeway;
  llf.y -= leeway;
  llf.z -= leeway;
  urb.x += leeway;
  urb.y += leeway;
  urb.z += leeway;

  cent.x = 0.33333333333 * (w->vert[0]->x + w->vert[1]->x + w->vert[2]->x);
  cent.y = 0.33333333333 * (w->vert[0]->y + w->vert[1]->y + w->vert[2]->y);
  cent.z = 0.33333333333 * (w->vert[0]->z + w->vert[1]->z + w->vert[2]->z);

  x_min = bisect(world->x_partitions, world->nx_parts, llf.x);
  if (urb.x < world->x_partitions[x_min + 1])
    x_max = x_min + 1;
  else
    x_max = bisect(world->x_partitions, world->nx_parts, urb.x) + 1;

  y_min = bisect(world->y_partitions, world->ny_parts, llf.y);
  if (urb.y < world->y_partitions[y_min + 1])
    y_max = y_min + 1;
  else
    y_max = bisect(world->y_partitions, world->ny_parts, urb.y) + 1;

  z_min = bisect(world->z_partitions, world->nz_parts, llf.z);
  if (urb.z < world->z_partitions[z_min + 1])
    z_max = z_min + 1;
  else
    z_max = bisect(world->z_partitions, world->nz_parts, urb.z) + 1;

  if ((z_max - z_min) * (y_max - y_min) * (x_max - x_min) == 1) {
    h = z_min + (world->nz_parts - 1) * (y_min + (world->ny_parts - 1) * x_min);
    where_am_i = localize_wall(w, world->subvol[h].local_storage);
    if (where_am_i == NULL)
      return NULL;

    if (wall_to_vol(where_am_i, &(world->subvol[h])) == NULL)
      return NULL;

    return where_am_i;
  }

  for (i = x_min; i < x_max; i++) {
    if (cent.x < world->x_partitions[i])
      break;
  }
  for (j = y_min; j < y_max; j++) {
    if (cent.y < world->y_partitions[j])
      break;
  }
  for (k = z_min; k < z_max; k++) {
    if (cent.z < world->z_partitions[k])
      break;
  }

  h = (k - 1) +
      (world->nz_parts - 1) * ((j - 1) + (world->ny_parts - 1) * (i - 1));
  where_am_i = localize_wall(w, world->subvol[h].local_storage);
  if (where_am_i == NULL)
    return NULL;

  for (k = z_min; k < z_max; k++) {
    for (j = y_min; j < y_max; j++) {
      for (i = x_min; i < x_max; i++) {
        h = k + (world->nz_parts - 1) * (j + (world->ny_parts - 1) * i);
        llf.x = world->x_fineparts[world->subvol[h].llf.x] - leeway;
        llf.y = world->y_fineparts[world->subvol[h].llf.y] - leeway;
        llf.z = world->z_fineparts[world->subvol[h].llf.z] - leeway;
        urb.x = world->x_fineparts[world->subvol[h].urb.x] + leeway;
        urb.y = world->y_fineparts[world->subvol[h].urb.y] + leeway;
        urb.z = world->z_fineparts[world->subvol[h].urb.z] + leeway;

        if (wall_in_box(w->vert, &(w->normal), w->d, &llf, &urb)) {
          if (wall_to_vol(where_am_i, &(world->subvol[h])) == NULL)
            return NULL;
        }
      }
    }
  }

  return where_am_i;
}

/***************************************************************************
distribute_object:
  In: an object
  Out: 0 on success, 1 on memory allocation failure.  The object's walls
       are copied to local memory and the wall lists in the appropriate
       subvolumes are set to refer to that wall.  The object's own copy
       of the wall is deallocated and it is set to point to the new version.
  Note: this function is recursive and is called on any children of the
        object passed to it.
***************************************************************************/
int distribute_object(struct volume *world, struct object *parent) {
  struct object *o; /* Iterator for child objects */
  int i;
  long long vert_index; /* index of the vertex in the global array
                     "world->all_vertices" */

  if (parent->object_type == BOX_OBJ || parent->object_type == POLY_OBJ) {
    for (i = 0; i < parent->n_walls; i++) {
      if (parent->wall_p[i] == NULL)
        continue; /* Wall removed. */

      parent->wall_p[i] = distribute_wall(world, parent->wall_p[i]);

      if (parent->wall_p[i] == NULL)
        mcell_allocfailed("Failed to distribute wall %d on object %s.", i,
                          parent->sym->name);

      /* create information about shared vertices */
      if (world->create_shared_walls_info_flag) {
        vert_index = (long long)(parent->wall_p[i]->vert[0] - world->all_vertices);
        push_wall_to_list(&(world->walls_using_vertex[vert_index]),
                          parent->wall_p[i]);
        vert_index = (long long)(parent->wall_p[i]->vert[1] - world->all_vertices);
        push_wall_to_list(&(world->walls_using_vertex[vert_index]),
                          parent->wall_p[i]);
        vert_index = (long long)(parent->wall_p[i]->vert[2] - world->all_vertices);
        push_wall_to_list(&(world->walls_using_vertex[vert_index]),
                          parent->wall_p[i]);
      }
    }
    if (parent->walls != NULL) {
      free(parent->walls);
      parent->walls = NULL; /* Use wall_p from now on! */
    }
  } else if (parent->object_type == META_OBJ) {
    for (o = parent->first_child; o != NULL; o = o->next) {
      if (distribute_object(world, o) != 0)
        return 1;
    }
  }

  return 0;
}

/***************************************************************************
distribute_world:
  In: No arguments.
  Out: 0 on success, 1 on memory allocation failure.  Every geometric object
       is distributed to local memory and into appropriate subvolumes.
***************************************************************************/
int distribute_world(struct volume *world) {
  struct object *o; /* Iterator for objects in the world */

  for (o = world->root_instance; o != NULL; o = o->next) {
    if (distribute_object(world, o) != 0)
      return 1;
  }

  return 0;
}

/***************************************************************************
closest_pt_point_triangle:
  In:  p - point
       a,b,c - vectors defining the vertices of the triangle.
  Out: final_result - closest point on triangle ABC to a point p.
       The code is adapted from "Real-time Collision Detection" by Christer
Ericson, ISBN 1-55860-732-3, p.141.

***************************************************************************/
void closest_pt_point_triangle(struct vector3 *p, struct vector3 *a,
                               struct vector3 *b, struct vector3 *c,
                               struct vector3 *final_result) {
  struct vector3 ab, ac, ap, bp, cp, result1;
  double d1, d2, d3, d4, vc, d5, d6, vb, va, denom, v, w;

  /* Check if P in vertex region outside A */
  vectorize(a, b, &ab);
  vectorize(a, c, &ac);
  vectorize(a, p, &ap);
  d1 = dot_prod(&ab, &ap);
  d2 = dot_prod(&ac, &ap);
  if (d1 <= 0.0f && d2 <= 0.0f) {
    memcpy(final_result, a,
           sizeof(struct vector3)); /* barycentric coordinates (1,0,0) */
    return;
  }

  /* Check if P in vertex region outside B */
  vectorize(b, p, &bp);
  d3 = dot_prod(&ab, &bp);
  d4 = dot_prod(&ac, &bp);
  if (d3 >= 0.0f && d4 <= d3) {
    memcpy(final_result, b,
           sizeof(struct vector3)); /* barycentric coordinates (0,1,0) */
    return;
  }

  /* Check if P in edge region of AB, if so return projection of P onto AB */
  vc = d1 * d4 - d3 * d2;
  if (vc <= 0.0f && d1 >= 0.0f && d3 <= 0.0f) {
    v = d1 / (d1 - d3);
    scalar_prod(&ab, v, &result1);
    vect_sum(a, &result1, final_result);
    return; /* barycentric coordinates (1-v,v,0) */
  }

  /* Check if P in vertex region outside C */
  vectorize(c, p, &cp);
  d5 = dot_prod(&ab, &cp);
  d6 = dot_prod(&ac, &cp);
  if (d6 >= 0.0f && d5 <= d6) {
    memcpy(final_result, c,
           sizeof(struct vector3)); /* barycentric coordinates (0,0,1) */
    return;
  }

  /* Check if P in edge region of AC, if so return projection of P onto AC */
  vb = d5 * d2 - d1 * d6;
  if (vb <= 0.0f && d2 >= 0.0f && d6 <= 0.0f) {
    w = d2 / (d2 - d6);
    scalar_prod(&ac, w, &result1);
    vect_sum(a, &result1, final_result);
    return; /* barycentric coordinates (0, 1-w,w) */
  }

  /* Check if P in edge region of BC, if so return projection of P onto BC */
  va = d3 * d6 - d5 * d4;
  if (va <= 0.0f && (d4 - d3) >= 0.0f && (d5 - d6) >= 0.0f) {
    w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
    vectorize(b, c, &result1);
    scalar_prod(&result1, w, &result1);
    vect_sum(b, &result1, final_result);
    return; /*barycentric coordinates (0,1-w, w) */
  }

  /* P inside face region. Compute Q through its barycentric
     coordinates (u,v,w) */
  denom = 1.0f / (va + vb + vc);
  v = vb * denom;
  w = vc * denom;
  scalar_prod(&ab, v, &ab);
  scalar_prod(&ac, w, &ac);
  vect_sum(&ab, &ac, &result1);
  vect_sum(a, &result1, final_result);
  return; /* = u*a + v*b + w*c, u = va * denom = 1.0f - v -w */
}

/***************************************************************************
test_bounding_boxes:
  In:  llf1 - lower left corner of the 1st box
       urb1 - upper right back corner of the 1st box
       llf2 - lower left corner of the 2nd box
       urb2 - upper right back corner of the 2nd box
  Out: Returns 1 if boxes intersect, 0 - otherwise
       The code is adapted from "Real-time Collision Detection"
       by Christer Ericson, ISBN 1-55860-732-3, p.79.

***************************************************************************/
int test_bounding_boxes(struct vector3 *llf1, struct vector3 *urb1,
                        struct vector3 *llf2, struct vector3 *urb2) {
  /* Two boxes overlap only if they overlap on all three axes
     while their extent along each dimension is seen as an interval
     on the corresponding axis. */

  /* exit with no intersection is separated along axis */
  if ((urb1->x < llf2->x) || (llf1->x > urb2->x))
    return 0;
  if ((urb1->y < llf2->y) || (llf1->y > urb2->y))
    return 0;
  if ((urb1->z < llf2->z) || (llf1->z > urb2->z))
    return 0;
  /* Overlapping on all axis means that boxes are intersecting. */
  return 1;
}

/* Helper struct for release_onto_regions and vacuum_from_regions */
struct reg_rel_helper_data {
  struct reg_rel_helper_data *next;
  struct surface_grid *grid;
  unsigned int index;
  double my_area;
};

/***************************************************************************
vacuum_from_regions:
  In: a release site object
      a template surface molecule we're going to remove
      the number of molecules to remove
  Out: 0 on success, 1 on failure.  Molecules of the specified type are
       removed uniformly at random from the free area in the regions
       specified by the release site object.
  Note: if the user requests to remove more molecules than actually exist,
        the function will return success and not give a warning.  The only
        reason to return failure is an out of memory condition.
***************************************************************************/
static int vacuum_from_regions(struct volume *world,
                               struct release_site_obj *rso,
                               struct surface_molecule *sm, int n) {
  struct release_region_data *rrd;
  struct mem_helper *mh;
  struct reg_rel_helper_data *rrhd_head, *p;
  int n_rrhd;
  struct wall *w;
  struct surface_molecule *smp;

  rrd = rso->region_data;

  mh = create_mem(sizeof(struct reg_rel_helper_data), 1024);
  if (mh == NULL)
    return 1;

  rrhd_head = NULL;
  n_rrhd = 0;

  for (int n_object = 0; n_object < rrd->n_objects; n_object++) {
    if (rrd->walls_per_obj[n_object] == 0)
      continue;
    for (int n_wall = 0; n_wall < rrd->in_release[n_object]->nbits; n_wall++) {
      if (!get_bit(rrd->in_release[n_object], n_wall))
        continue;

      w = rrd->owners[n_object]->wall_p[n_wall];

      if (w->grid == NULL)
        continue;

      for (unsigned int n_tile = 0; n_tile < w->grid->n_tiles; n_tile++) {
        struct surface_molecule_list *sm_list = w->grid->sm_list[n_tile];
        if (sm_list && sm_list->sm) {
          smp = w->grid->sm_list[n_tile]->sm;
          if (smp != NULL) {
            if (smp->properties == sm->properties) {
              p = CHECKED_MEM_GET_NODIE(mh, "release region helper data");
              if (p == NULL)
                return 1;

              p->next = rrhd_head;
              p->grid = w->grid;
              p->index = n_tile;
              rrhd_head = p;

              n_rrhd++;
            }
          }
        }
      }
    }
  }

  for (p = rrhd_head; n < 0 && n_rrhd > 0 && p != NULL; p = p->next, n_rrhd--) {
    if (rng_dbl(world->rng) < ((double)(-n)) / ((double)n_rrhd)) {
      smp = p->grid->sm_list[p->index]->sm;
      smp->properties->population--;
      if ((smp->properties->flags & (COUNT_CONTENTS | COUNT_ENCLOSED)) != 0)
        count_region_from_scratch(world, (struct abstract_molecule *)smp, NULL,
                                  -1, NULL, smp->grid->surface, smp->t, NULL);
      smp->properties = NULL;
      p->grid->sm_list[p->index]->sm = NULL;
      p->grid->n_occupied--;
      if (smp->flags & IN_SCHEDULE) {
        smp->grid->subvol->local_storage->timer->defunct_count++; /* Tally for
                                                                    garbage
                                                                    collection
                                                                    */
      }

      n++;
    }
  }

  delete_mem(mh);

  return 0;
}

/***************************************************************************
release_onto_regions:
  In: a release site object
      a template surface molecule we're going to release
      the number of molecules to release
  Out: 0 on success, 1 on failure.  Molecules are released uniformly at
       random onto the free area in the regions specified by the release
      site object.
  Note: if the CCNNUM method is used, the number passed in is ignored.
***************************************************************************/
int release_onto_regions(struct volume *world, struct release_site_obj *rso,
                         struct surface_molecule *sm, int n) {
  struct mem_helper *mh;
  int i;
  unsigned int grid_index;
  double A, num_to_release;
  struct wall *w;

  struct release_region_data *rrd = rso->region_data;

  int success = 0, failure = 0;
  double seek_cost = 0;

  double max_A = rrd->cum_area_list[rrd->n_walls_included - 1];
  double est_sites_avail = (int)max_A;
  const double rel_list_gen_cost = 10.0; /* Just a guess */
  double pick_cost = rel_list_gen_cost * est_sites_avail;

  if (rso->release_number_method == DENSITYNUM) {
    num_to_release = rso->concentration * est_sites_avail / world->grid_density;
    if (num_to_release > (double)INT_MAX)
      mcell_error("Release site \"%s\" tries to release more than INT_MAX "
                  "(2147483647) molecules.",
                  rso->name);
    n = (int)(num_to_release);
  }

  if (n < 0)
    return vacuum_from_regions(world, rso, sm, n);

  const int too_many_failures = 10; /* Just a guess */
  while (n > 0) {
    if (failure >= success + too_many_failures) {
      seek_cost =
          n * (((double)(success + failure + 2)) / ((double)(success + 1)));
    }
    if (seek_cost < pick_cost) {
      A = rng_dbl(world->rng) * max_A;
      i = bisect_high(rrd->cum_area_list, rrd->n_walls_included, A);
      w = rrd->owners[rrd->obj_index[i]]->wall_p[rrd->wall_index[i]];

      if (w->grid == NULL) {
        if (create_grid(world, w, NULL))
          return 1;
      }
      if (i)
        A -= rrd->cum_area_list[i - 1];
      grid_index = (unsigned int)((w->grid->n * w->grid->n) * (A / w->area));
      if (grid_index >= w->grid->n_tiles)
        grid_index = w->grid->n_tiles - 1;

        struct surface_molecule_list *sm_list = w->grid->sm_list[grid_index];
        if (sm_list && sm_list->sm) {
          failure++;
        }
        else {
          struct vector3 pos3d = {.x = 0, .y = 0, .z = 0};
          if (place_single_molecule(world, w, grid_index, sm->properties,
                                    sm->flags, rso->orientation, sm->t, sm->t2,
                                    sm->birthday, sm->periodic_box, &pos3d) == NULL) {
            struct vector3 llf, urb;
            if (world->periodic_box_obj) {
              struct polygon_object *p = (struct polygon_object*)(world->periodic_box_obj->contents);
              struct subdivided_box *sb = p->sb;
              llf = (struct vector3) {sb->x[0], sb->y[0], sb->z[0]};
              urb = (struct vector3) {sb->x[1], sb->y[1], sb->z[1]};
            }
            if (world->periodic_box_obj && !point_in_box(&llf, &urb, &pos3d)) {
              mcell_log("Cannot release '%s' outside of periodic boundaries.",
                        sm->properties->sym->name);
              failure++;
              continue;
            }
            return 1;
          }
          success++;
          n--;
        }
    } else {
      if (world->periodic_box_obj) {
        return 1;
      }
      mh = create_mem(sizeof(struct reg_rel_helper_data), 1024);
      if (mh == NULL)
        return 1;

      struct reg_rel_helper_data *rrhd_head = NULL;
      int n_rrhd = 0;
      max_A = 0;
      for (int n_object = 0; n_object < rrd->n_objects; n_object++) {
        if (rrd->walls_per_obj[n_object] == 0)
          continue;
        for (int n_wall = 0; n_wall < rrd->in_release[n_object]->nbits;
             n_wall++) {
          if (!get_bit(rrd->in_release[n_object], n_wall))
            continue;

          w = rrd->owners[n_object]->wall_p[n_wall];

          if (w->grid == NULL) {
            if (create_grid(world, w, NULL)) {
              delete_mem(mh);
              return 1;
            }
          } else if (w->grid->n_occupied == w->grid->n_tiles)
            continue;

          A = w->area / (w->grid->n_tiles);

          for (unsigned int n_tile = 0; n_tile < w->grid->n_tiles; n_tile++) {
            if (w->grid->sm_list[n_tile] == NULL || w->grid->sm_list[n_tile]->sm == NULL) {
              struct reg_rel_helper_data *new_rrd =
                  CHECKED_MEM_GET_NODIE(mh, "release region helper data");
              if (new_rrd == NULL)
                return 1;

              new_rrd->next = rrhd_head;
              new_rrd->grid = w->grid;
              new_rrd->index = n_tile;
              new_rrd->my_area = A;
              max_A += A;

              rrhd_head = new_rrd;
              n_rrhd++;
            }
          }
        }
      }

      for (struct reg_rel_helper_data *this_rrd = rrhd_head;
           this_rrd != NULL && n > 0; this_rrd = this_rrd->next) {
        if (n >= n_rrhd ||
            rng_dbl(world->rng) < (this_rrd->my_area / max_A) * ((double)n)) {
          struct vector3 pos3d = {.x = 0, .y = 0, .z = 0};
          if (place_single_molecule(world, this_rrd->grid->surface,
                                    this_rrd->index, sm->properties, sm->flags,
                                    rso->orientation, sm->t, sm->t2,
                                    sm->birthday, sm->periodic_box, &pos3d) == NULL) {
            return 1;
          }

          n--;
          n_rrhd--;
        }
        max_A -= this_rrd->my_area;
      }

      delete_mem(mh);

      if (n > 0) {
        switch (world->notify->mol_placement_failure) {
        case WARN_COPE:
          break;

        case WARN_WARN:
          mcell_warn("Could not release %d of %s (surface full).", n,
                     sm->properties->sym->name);
          break;

        case WARN_ERROR:
          mcell_error("Could not release %d of %s (surface full).", n,
                      sm->properties->sym->name);
          /*return 1;*/

        default:
          UNHANDLED_CASE(world->notify->mol_placement_failure);
        }
        break;
      }
    }
  }

  return 0;
}

/****************************************************************************
place_single_molecule:
   In: state: the simulation state
       w: the wall to receive the surface molecule
       grid_index:
       spec: the molecule
       flags:
       orientation: the orientation of the molecule
   Out: the new surface molecule is returned on success, NULL otherwise
*****************************************************************************/
struct surface_molecule *place_single_molecule(struct volume *state,
                                               struct wall *w,
                                               unsigned int grid_index,
                                               struct species *spec,
                                               short flags, short orientation,
                                               double t, double t2,
                                               double birthday,
                                               struct periodic_image *periodic_box,
                                               struct vector3 *pos3d) {

  struct vector2 s_pos;

  if (state->randomize_smol_pos)
    grid2uv_random(w->grid, grid_index, &s_pos, state->rng);
  else
    grid2uv(w->grid, grid_index, &s_pos);
  uv2xyz(&s_pos, w, pos3d);

  struct vector3 llf, urb;
  if (state->periodic_box_obj) {
    struct polygon_object *p = (struct polygon_object*)(state->periodic_box_obj->contents);
    struct subdivided_box *sb = p->sb;
    llf = (struct vector3) {sb->x[0], sb->y[0], sb->z[0]};
    urb = (struct vector3) {sb->x[1], sb->y[1], sb->z[1]};
  }

  if (state->periodic_box_obj && !point_in_box(&llf, &urb, pos3d)) {
    mcell_log("Cannot release '%s' outside of periodic boundaries.",
              spec->sym->name);
    return NULL;
  }

  struct subvolume *gsv = NULL;
  gsv = find_subvolume(state, pos3d, gsv);

  struct surface_molecule *new_sm;
  new_sm = (struct surface_molecule *)CHECKED_MEM_GET(gsv->local_storage->smol,
                                                      "surface molecule");
  if (new_sm == NULL)
    return NULL;
  new_sm->t = t;
  new_sm->t2 = t2;
  new_sm->birthday = birthday;
  new_sm->birthplace = w->birthplace->smol;
  new_sm->id = state->current_mol_id++;
  new_sm->grid_index = grid_index;
  new_sm->s_pos.u = s_pos.u;
  new_sm->s_pos.v = s_pos.v;
  new_sm->properties = spec;
  new_sm->periodic_box = CHECKED_MALLOC_STRUCT(struct periodic_image,
    "periodic image descriptor");
  new_sm->periodic_box->x = periodic_box->x;
  new_sm->periodic_box->y = periodic_box->y;
  new_sm->periodic_box->z = periodic_box->z;

  if (orientation == 0)
    new_sm->orient = (rng_uint(state->rng) & 1) ? 1 : -1;
  else
    new_sm->orient = orientation;

  new_sm->grid = w->grid;

  if (w->grid->sm_list[grid_index] == NULL) {
    struct surface_molecule_list *sm_entry = CHECKED_MALLOC_STRUCT(
      struct surface_molecule_list, "surface molecule list");
    sm_entry->next = NULL;
    w->grid->sm_list[grid_index] = sm_entry;
  }
  w->grid->sm_list[grid_index]->sm = new_sm;
  w->grid->n_occupied++;
  new_sm->properties->population++;

  new_sm->flags = flags;

  if (new_sm->properties->space_step > 0)
    new_sm->flags |= ACT_DIFFUSE;

  if ((new_sm->properties->flags & COUNT_ENCLOSED) != 0)
    new_sm->flags |= COUNT_ME;

  if (new_sm->properties->flags & (COUNT_CONTENTS | COUNT_ENCLOSED))
    count_region_from_scratch(state, (struct abstract_molecule *)new_sm, NULL,
                              1, NULL, new_sm->grid->surface, new_sm->t,
                              new_sm->periodic_box);

  if (schedule_add(gsv->local_storage->timer, new_sm)) {
    mcell_allocfailed("Failed to add volume molecule '%s' to scheduler.",
                      new_sm->properties->sym->name);
    return NULL;
  }

  return new_sm;
}

/****************************************************************************
push_wall_to_list:
   In: head of the linked list of walls
       wall to be added to the list
   Out: none. The linked list is updated.
   Note: No check for duplicates is performed
*****************************************************************************/
void push_wall_to_list(struct wall_list **wall_nbr_head, struct wall *w) {
  struct wall_list *old_head, *wlp;

  old_head = *wall_nbr_head;

  wlp = CHECKED_MALLOC_STRUCT(struct wall_list, "wall_list");
  wlp->this_wall = w;

  if (old_head == NULL) {
    wlp->next = NULL;
    old_head = wlp;
  } else {
    wlp->next = old_head;
    old_head = wlp;
  }

  *wall_nbr_head = old_head;
}

/*************************************************************************
delete_wall_list:
   In: linked list of walls
   Out: none. The memory is freed.
*************************************************************************/
void delete_wall_list(struct wall_list *wl_head) {
  struct wall_list *nnext;
  while (wl_head != NULL) {
    nnext = wl_head->next;
    free(wl_head);
    wl_head = nnext;
  }
}

/*************************************************************************
find_nbr_walls_shared_one_vertex:
   In: the origin wall
       array with information about which vertices of the origin wall
          are shared with neighbor wall (they are indices in the
          global "world->walls_using_vertex" array).
   Out: linked list of the neighbor walls that have only one common
        vertex with the origin wall (not edge-to-edge walls, but
        vertex-to-vertex walls).
   Note: the "origin" wall is not included in the list
**************************************************************************/
struct wall_list *find_nbr_walls_shared_one_vertex(struct volume *world,
                                                   struct wall *origin,
                                                   long long int *shared_vert) {
  int i;
  struct wall_list *wl;
  struct wall_list *head = NULL;

  if (!world->create_shared_walls_info_flag)
    mcell_internal_error("Function 'find_nbr_walls_shared_one_vertex()' is "
                         "called but shared walls information is not created.");

  for (i = 0; i < 3; i++) {
    if (shared_vert[i] >= 0) {
      for (wl = world->walls_using_vertex[shared_vert[i]]; wl != NULL;
           wl = wl->next) {
        if (wl->this_wall == origin)
          continue;

        if (!walls_share_full_edge(origin, wl->this_wall)) {
          push_wall_to_list(&head, wl->this_wall);
        }
      }
    }
  }

  return head;
}

/***********************************************************************
walls_share_full_edge:
  In: two walls
  Out: 1 if the walls share a full edge, 0 - otherwise.
       Here by "full" we mean that the shared edge has two endpoints
       that are the vertices of both walls w1 and w2.
************************************************************************/
int walls_share_full_edge(struct wall *w1, struct wall *w2) {
  int i, k;
  int count = 0; /* count number of shared vertices between two walls */

  for (i = 0; i < 3; i++) {
    for (k = 0; k < 3; k++) {
      if (!distinguishable_vec3(w1->vert[i], w2->vert[k], EPS_C))
        count++;
    }
  }

  if (count == 2)
    return 1;

  return 0;
}

/***********************************************************************
find_region_by_wall:
  In:  wall
  Out: an object's region list if the wall belongs to one, NULL - otherwise.
  Note: regions called "ALL" or the ones that have ALL_ELEMENTS are not
        included in the return "region list".  This is done intentionally
        since the function is used to determine region border and the
        above regions do not have region borders.
************************************************************************/
struct region_list *find_region_by_wall(struct wall *this_wall) {

  struct region_list *rlp_head = NULL;
  for (struct region_list *rlp = this_wall->parent_object->regions; rlp != NULL;
    rlp = rlp->next) {

    struct region *rp = rlp->reg;
    if ((strcmp(rp->region_last_name, "ALL") == 0) ||
        (rp->region_has_all_elements))
      continue;

    if (rp->membership == NULL)
      mcell_internal_error("Missing region membership for '%s'.",
                           rp->sym->name);

    if (get_bit(rp->membership, this_wall->side)) {
      struct region_list *rlps = CHECKED_MALLOC_STRUCT(struct region_list,
        "region_list");
      rlps->reg = rp;

      if (rlp_head == NULL) {
        rlps->next = NULL;
        rlp_head = rlps;
      } else {
        rlps->next = rlp_head;
        rlp_head = rlps;
      }
    }
  }

  return rlp_head;
}

/************************************************************************
find_regions_names_by_wall:
  In:  wall:
       ignore_regs:
  Out: linked list of wall's regions names if wall belongs to region, 
       NULL - otherwise.  Also number of regions is set.
  Note: regions called "ALL" or the ones that have ALL_ELEMENTS are not 
        included in the return regions names list. This is done intentionally
        since the function is used with dynamic geometries to place a surface
        molecule.
************************************************************************/
struct name_list *find_regions_names_by_wall(
    struct wall *w, struct string_buffer *ignore_regs)
{
  struct name_list *nl_head = NULL;
  struct region_list *reg_list_ptr_head = find_region_by_wall(w);

  struct region_list *reg_list_ptr;
  for (reg_list_ptr = reg_list_ptr_head; reg_list_ptr != NULL;) {
    struct region *reg_ptr = reg_list_ptr->reg;

    // Disregard regions which were just added during a dynamic geometry event
    if ((ignore_regs) && (is_string_present_in_string_array(
        reg_ptr->sym->name, ignore_regs->strings, ignore_regs->n_strings))) {
      reg_list_ptr_head = reg_list_ptr->next;
      free(reg_list_ptr);
      reg_list_ptr = reg_list_ptr_head;
      continue;
    }

    struct name_list *nl = CHECKED_MALLOC_STRUCT(struct name_list, "name_list");
    nl->name = alloc_sprintf("%s", reg_ptr->sym->name);   
    nl->prev = NULL;

    if (nl_head == NULL) {
      nl->next = NULL;
      nl_head = nl;
    }
    else {
      nl->next = nl_head;
      nl_head = nl;
    }
    reg_list_ptr_head = reg_list_ptr->next;
    free(reg_list_ptr);
    reg_list_ptr = reg_list_ptr_head;
  }

  return nl_head;
}

/***********************************************************************
find_restricted_regions_by_wall:
  In: wall
      surface molecule
  Out: an object's region list if the wall belongs to the region
          that is restrictive (REFL/ABSORB) to the surface molecule
       NULL - if no such regions found
  Note: regions called "ALL" or the ones that have ALL_ELEMENTS are not
        included in the return "region list".
************************************************************************/
struct region_list *
find_restricted_regions_by_wall(struct volume *world, struct wall *this_wall,
                                struct surface_molecule *sm) {

  if ((sm->properties->flags & CAN_REGION_BORDER) == 0)
    return NULL;

  struct rxn *matching_rxns[MAX_MATCHING_RXNS];
  for (int kk = 0; kk < MAX_MATCHING_RXNS; kk++) {
    matching_rxns[kk] = NULL;
  }

  int num_matching_rxns = trigger_intersect(
      world->reaction_hash, world->rx_hashsize, world->all_mols,
      world->all_volume_mols, world->all_surface_mols, sm->properties->hashval,
      (struct abstract_molecule *)sm, sm->orient, this_wall, matching_rxns, 1,
      1, 1);

  struct species *restricted_surf_class[MAX_MATCHING_RXNS];
  int num_res = 0;
  for (int kk = 0; kk < num_matching_rxns; kk++) {
    if ((matching_rxns[kk]->n_pathways == RX_REFLEC) ||
        (matching_rxns[kk]->n_pathways == RX_ABSORB_REGION_BORDER)) {
      restricted_surf_class[num_res++] = matching_rxns[kk]->players[1];
    }
  }

  struct region *rp;
  struct region_list *rlp_head = NULL;
  for (struct region_list *rlp = this_wall->parent_object->regions; rlp != NULL;
    rlp = rlp->next) {
    rp = rlp->reg;

    if (rp->membership == NULL) {
      mcell_internal_error("Missing region membership for '%s'.", rp->sym->name);
    }

    if (rp->surf_class == NULL) {
      continue;
    }

    if ((strcmp(rp->region_last_name, "ALL") == 0) ||
        (rp->region_has_all_elements)) {
      continue;
    }

    if (get_bit(rp->membership, this_wall->side)) {
      /* is this region's boundary restricted for surface molecule? */
      for (int i = 0; i < num_res; ++i) {
        if (rp->surf_class == restricted_surf_class[i]) {
          struct region_list *rlps = CHECKED_MALLOC_STRUCT(struct region_list,
            "region_list");
          rlps->reg = rp;
          if (rlp_head == NULL) {
            rlps->next = NULL;
            rlp_head = rlps;
          } else {
            rlps->next = rlp_head;
            rlp_head = rlps;
          }
          break;
        }
      }
    }
  }
  return rlp_head;
}

/***********************************************************************
find_restricted_regions_by_object:
  In: object
      surface molecule
  Out: an object's region list that are restrictive (REFL/ABSORB)
       to the surface molecule
       NULL - if no such regions found
  Note: regions called "ALL" or the ones that have ALL_ELEMENTS are not
        included in the return "region list".
************************************************************************/
struct region_list *
find_restricted_regions_by_object(struct volume *world, struct object *obj,
                                  struct surface_molecule *sm) {
  struct region *rp;
  struct region_list *rlp, *rlps, *rlp_head = NULL;
  int kk, i, wall_idx = INT_MIN;
  struct rxn *matching_rxns[MAX_MATCHING_RXNS];

  if ((sm->properties->flags & CAN_REGION_BORDER) == 0)
    return NULL;

  for (kk = 0; kk < MAX_MATCHING_RXNS; kk++) {
    matching_rxns[kk] = NULL;
  }

  for (rlp = obj->regions; rlp != NULL; rlp = rlp->next) {
    rp = rlp->reg;
    if ((strcmp(rp->region_last_name, "ALL") == 0) ||
        (rp->region_has_all_elements)) {
      continue;
    }

    /* find any wall that belongs to this region */
    for (i = 0; i < obj->n_walls; i++) {
      if (get_bit(rp->membership, i)) {
        wall_idx = i;
        break;
      }
    }

    if (wall_idx < 0)
      mcell_internal_error("Cannot find wall in the region.");

    int num_matching_rxns = 0;
    if (rp->surf_class) {
      num_matching_rxns = find_unimol_reactions_with_surf_classes(
          world->reaction_hash, world->rx_hashsize,
          (struct abstract_molecule *)sm, obj->wall_p[wall_idx],
          sm->properties->hashval, sm->orient, num_matching_rxns, 1, 1, 1,
          matching_rxns);
      num_matching_rxns = find_surface_mol_reactions_with_surf_classes(
          world->reaction_hash, world->rx_hashsize, world->all_mols,
          world->all_surface_mols, sm->orient, rp->surf_class,
          num_matching_rxns, 1, 1, 1, matching_rxns);
    }

    for (kk = 0; kk < num_matching_rxns; kk++) {
      if ((matching_rxns[kk]->n_pathways == RX_REFLEC) ||
          (matching_rxns[kk]->n_pathways == RX_ABSORB_REGION_BORDER)) {
        rlps = CHECKED_MALLOC_STRUCT(struct region_list, "region_list");
        rlps->reg = rp;

        if (rlp_head == NULL) {
          rlps->next = NULL;
          rlp_head = rlps;
        } else {
          rlps->next = rlp_head;
          rlp_head = rlps;
        }
      }
    }
  }

  return rlp_head;
}

/***********************************************************************
are_restricted_regions_for_species_on_object:
  In: object
      surface molecule
  Out: 1 if there are regions that are restrictive (REFL/ABSORB)
       to the surface molecule on this object
       0 - if no such regions found
************************************************************************/
int are_restricted_regions_for_species_on_object(struct volume *world,
                                                 struct object *obj,
                                                 struct surface_molecule *sm) {
  int wall_idx = INT_MIN;
  struct rxn *matching_rxns[MAX_MATCHING_RXNS];

  if ((sm->properties->flags & CAN_REGION_BORDER) == 0)
    return 0;

  for (int kk = 0; kk < MAX_MATCHING_RXNS; kk++) {
    matching_rxns[kk] = NULL;
  }

  for (struct region_list *rlp = obj->regions; rlp != NULL; rlp = rlp->next) {
    struct region *rp = rlp->reg;
    if ((strcmp(rp->region_last_name, "ALL") == 0) ||
        (rp->region_has_all_elements)) {
      continue;
    }

    /* find any wall that belongs to this region */
    for (int i = 0; i < obj->n_walls; i++) {
      if (get_bit(rp->membership, i)) {
        wall_idx = i;
        break;
      }
    }

    if (wall_idx < 0) {
      mcell_internal_error("Cannot find wall in the region.");
    }

    int num_matching_rxns = trigger_intersect(
        world->reaction_hash, world->rx_hashsize, world->all_mols,
        world->all_volume_mols, world->all_surface_mols,
        sm->properties->hashval, (struct abstract_molecule *)sm, sm->orient,
        obj->wall_p[wall_idx], matching_rxns, 1, 1, 1);

    if (num_matching_rxns > 0) {
      for (int kk = 0; kk < num_matching_rxns; kk++) {
        if ((matching_rxns[kk]->n_pathways == RX_REFLEC) ||
            (matching_rxns[kk]->n_pathways == RX_ABSORB_REGION_BORDER)) {
          return 1;
        }
      }
    }
  }

  return 0;
}

/***********************************************************************
is_wall_edge_region_border:
  In: wall
      wall's edge
  Out: 1 if the edge is a region's border, and 0 - otherwise.
  Note: we do not specify any particular region here, any region will
        suffice
************************************************************************/
int is_wall_edge_region_border(struct wall *this_wall, struct edge *this_edge) {
  struct region_list *rlp, *rlp_head;
  struct region *rp;
  void *key;
  unsigned int keyhash;

  int is_region_border = 0; /* flag */

  rlp_head = find_region_by_wall(this_wall);

  /* If this wall is not a part of any region (note that we do not consider
     region called ALL here) */
  if (rlp_head == NULL)
    return is_region_border;

  for (rlp = rlp_head; rlp != NULL; rlp = rlp->next) {
    rp = rlp->reg;

    if (rp->boundaries == NULL)
      mcell_internal_error("Region '%s' of the object '%s' has no boundaries.",
                           rp->region_last_name,
                           this_wall->parent_object->sym->name);

    keyhash = (unsigned int)(intptr_t)(this_edge);
    key = (void *)(this_edge);

    if (pointer_hash_lookup(rp->boundaries, key, keyhash)) {
      is_region_border = 1;
      break;
    }
  }

  if (rlp_head != NULL)
    delete_void_list((struct void_list *)rlp_head);

  return is_region_border;
}

/***********************************************************************
is_wall_edge_restricted_region_border:
  In: wall
      wall's edge
      surface molecule
  Out: 1 if the edge is a restricted region's border for above surface molecule
       0 - otherwise.
  Note: we do not specify any particular region here, any region will
        suffice for which special reactions (REFL/ABSORB) are defined.
************************************************************************/
int is_wall_edge_restricted_region_border(struct volume *world,
                                          struct wall *this_wall,
                                          struct edge *this_edge,
                                          struct surface_molecule *sm) {
  struct region_list *rlp, *rlp_head;
  struct region *rp;
  void *key;
  unsigned int keyhash;

  int is_region_border = 0; /* flag */

  rlp_head = find_restricted_regions_by_wall(world, this_wall, sm);

  /* If this wall is not a part of any region (note that we do not consider
     region called ALL here) */
  if (rlp_head == NULL)
    return is_region_border;

  for (rlp = rlp_head; rlp != NULL; rlp = rlp->next) {
    rp = rlp->reg;
    if (rp->boundaries == NULL)
      mcell_internal_error("Region '%s' of the object '%s' has no boundaries.",
                           rp->region_last_name,
                           this_wall->parent_object->sym->name);

    keyhash = (unsigned int)(intptr_t)(this_edge);
    key = (void *)(this_edge);

    if (pointer_hash_lookup(rp->boundaries, key, keyhash)) {
      is_region_border = 1;
      break;
    }
  }

  if (rlp_head != NULL)
    delete_void_list((struct void_list *)rlp_head);

  return is_region_border;
}

/*************************************************************************
find_shared_edge_index_of_neighbor_wall:
  In: original wall
      neighbor wall
  Out: index of the shared edge in the coordinate system of neighbor wall.

**************************************************************************/
int find_shared_edge_index_of_neighbor_wall(struct wall *orig_wall,
                                            struct wall *nbr_wall) {
  int nbr_edge_ind = -1;
  int shared_vert_ind_1 = -1, shared_vert_ind_2 = -1;

  find_shared_vertices_for_neighbor_walls(
      orig_wall, nbr_wall, &shared_vert_ind_1, &shared_vert_ind_2);

  if ((shared_vert_ind_1 + shared_vert_ind_2) == 1) {
    nbr_edge_ind = 0;
  } else if ((shared_vert_ind_1 + shared_vert_ind_2) == 2) {
    nbr_edge_ind = 2;
  } else if ((shared_vert_ind_1 + shared_vert_ind_2) == 3) {
    nbr_edge_ind = 1;
  } else {
    mcell_internal_error(
        "Error in the function 'find_shared_edge_index_of_neighbor_wall()");
  }

  return nbr_edge_ind;
}

/****************************************************************************
find_neighbor_wall_and_edge:
  In: orig_wall: wall
      orig_edge_ind: wall edge index (in the coordinate system of "wall")
      nbr_wall: neighbor wall (return value)
      nbr_edge_ind: index of the edge in the coordinate system of "neighbor
                    wall" that is shared with "wall" and coincides with the
                    edge with "wall edge index" (return value)

****************************************************************************/
void find_neighbor_wall_and_edge(struct wall *orig_wall, int orig_edge_ind,
                                 struct wall **nbr_wall, int *nbr_edge_ind) {
  struct wall *w;
  struct vector3 *vert_A = NULL, *vert_B = NULL;

  switch (orig_edge_ind) {
  case 0:
    vert_A = orig_wall->vert[0];
    vert_B = orig_wall->vert[1];
    break;
  case 1:
    vert_A = orig_wall->vert[1];
    vert_B = orig_wall->vert[2];
    break;
  case 2:
    vert_A = orig_wall->vert[2];
    vert_B = orig_wall->vert[0];
    break;
  default:
    mcell_internal_error("Error in function 'find_neighbor_wall_and_edge()'.");
    /*break;*/
  }

  for (int ii = 0; ii < 3; ii++) {
    w = orig_wall->nb_walls[ii];
    if (w == NULL)
      continue;

    if (wall_contains_both_vertices(w, vert_A, vert_B)) {
      *nbr_wall = w;
      *nbr_edge_ind = find_shared_edge_index_of_neighbor_wall(orig_wall, w);
      break;
    }
  }
}

/***************************************************************************
wall_contains_both_vertices:
  In: wall
      two vertices
  Out: Returns 1 if the wall contains both vertices above, and 0 otherwise.
***************************************************************************/
int wall_contains_both_vertices(struct wall *w, struct vector3 *vert_A,
                                struct vector3 *vert_B) {
  int count = 0;

  for (int ii = 0; ii < 3; ii++) {
    struct vector3 *v = w->vert[ii];

    if ((!distinguishable_vec3(v, vert_A, EPS_C)) ||
        (!(distinguishable_vec3(v, vert_B, EPS_C)))) {
      count++;
    }
  }

  if (count == 2)
    return 1;
  else
    return 0;
}

/*******************************************************************
are_walls_coincident:
  In: first wall
      second wall
      accuracy of the comparison
  Out: 0 if the walls are not coincident
       1 if the walls are coincident
*******************************************************************/
int are_walls_coincident(struct wall *w1, struct wall *w2, double eps) {
  if ((w1 == NULL) || (w2 == NULL))
    return 0;

  int count = 0;

  if (!distinguishable_vec3(w1->vert[0], w2->vert[0], eps))
    count++;
  if (!distinguishable_vec3(w1->vert[0], w2->vert[1], eps))
    count++;
  if (!distinguishable_vec3(w1->vert[0], w2->vert[2], eps))
    count++;
  if (!distinguishable_vec3(w1->vert[1], w2->vert[0], eps))
    count++;
  if (!distinguishable_vec3(w1->vert[1], w2->vert[1], eps))
    count++;
  if (!distinguishable_vec3(w1->vert[1], w2->vert[2], eps))
    count++;
  if (!distinguishable_vec3(w1->vert[2], w2->vert[0], eps))
    count++;
  if (!distinguishable_vec3(w1->vert[2], w2->vert[1], eps))
    count++;
  if (!distinguishable_vec3(w1->vert[2], w2->vert[2], eps))
    count++;

  if (count >= 3)
    return 1;

  return 0;
}

/******************************************************************
are_walls_coplanar:
  In: first wall
      second wall
      accuracy of the comparison
  Out: 1 if the walls are coplanar
       0 if walls are not coplanar
  Note: see "Real-time rendering" 2nd Ed., by Tomas Akenine-Moller and
        Eric Haines, pp. 590-591
******************************************************************/
int are_walls_coplanar(struct wall *w1, struct wall *w2, double eps) {

  /* find the plane equation of the second wall in the form (n*x + d2 = 0) */

  double d2, d1_0, d1_1, d1_2;

  d2 = -dot_prod(&(w2->normal), w2->vert[0]);

  /* check whether all vertices of the first wall satisfy
     plane equation of the second wall */
  d1_0 = dot_prod(&(w2->normal), w1->vert[0]) + d2;
  d1_1 = dot_prod(&(w2->normal), w1->vert[1]) + d2;
  d1_2 = dot_prod(&(w2->normal), w1->vert[2]) + d2;

  if ((!distinguishable(d1_0, 0, eps)) && (!distinguishable(d1_1, 0, eps)) &&
      (!distinguishable(d1_2, 0, eps))) {
    return 1;
  }

  return 0;
}

/**********************************************************************
* sorted_insert_wall_aux_list:
* In: linked list
*     new node
* Out: new node is added to the linked list in the sorted order
***********************************************************************/
void sorted_insert_wall_aux_list(struct wall_aux_list **headRef,
                                 struct wall_aux_list *newNode) {
  /* special case for the head end */
  if (*headRef == NULL || (*headRef)->d_prod >= newNode->d_prod) {
    newNode->next = *headRef;
    *headRef = newNode;
  } else {
    /* Locate the node before the point of insertion */
    struct wall_aux_list *curr = *headRef;
    while (curr->next != NULL && curr->next->d_prod < newNode->d_prod) {
      curr = curr->next;
    }
    newNode->next = curr->next;
    curr->next = newNode;
  }
}

/*******************************************************************
delete_wall_aux_list:
  In: linked list
  Out: None.  The linked list is deleted and memory is freed.
*******************************************************************/
void delete_wall_aux_list(struct wall_aux_list *head) {
  struct wall_aux_list *nnext;
  while (head != NULL) {
    nnext = head->next;
    free(head);
    head = nnext;
  }
}

/*****************************************************************
walls_belong_to_at_least_one_different_restricted_region:
  In: wall and surface molecule on it
      wall and surface molecule on it
  Out: 1 if both walls belong to at least one different restricted region
       relative to the properties of surface molecule,
       0 otherwise.
  Note: Wall can be belong to several regions simultaneously.
        Restricted region is one for the boundary of which the reactions
        REFL/ABSORB are declared.
******************************************************************/
int walls_belong_to_at_least_one_different_restricted_region(
    struct volume *world, struct wall *w1, struct surface_molecule *sm1,
    struct wall *w2, struct surface_molecule *sm2) {

  if ((w1 == NULL) || (w2 == NULL))
    return 0;

  struct region_list *rl_1 = find_restricted_regions_by_wall(world, w1, sm1);
  struct region_list *rl_2 = find_restricted_regions_by_wall(world, w2, sm2);

  if ((rl_1 == NULL) && (rl_2 == NULL))
    return 0;

  int error_code = 0;
  if (rl_1 == NULL) {
    /* Is wall 1 part of all restricted regions rl_2, then these restricted
     * regions just encompass wall 1 */
    if (wall_belongs_to_all_regions_in_region_list(w1, rl_2))
      error_code = 0;
    else
      error_code = 1;
    delete_region_list(rl_2);
    return error_code;
  }

  if (rl_2 == NULL) {
    /* Is wall 2 part of all restricted regions rl_1, then these restricted
     * regions just encompass wall 2 */
    if (wall_belongs_to_all_regions_in_region_list(w2, rl_1))
      error_code = 0;
    else
      error_code = 1;
    delete_region_list(rl_1);
    return error_code;
  }

  for (struct region_list *rl_t1 = rl_1; rl_t1 != NULL; rl_t1 = rl_t1->next) {
    struct region *rp_1 = rl_t1->reg;

    if (!region_belongs_to_region_list(rp_1, rl_2)) {
      error_code = 1;
      break;
    }
  }

  delete_region_list(rl_1);
  delete_region_list(rl_2);
  return error_code;
}

/**************************************************************************
region_belongs_to_region_list:
  In: region
      linked list of regions
  Out: 1 if region is present in the linked list of regions
       0 otherwise
***************************************************************************/
int region_belongs_to_region_list(struct region *rp, struct region_list *head) {
  int found = 0;

  for (struct region_list *rlp = head; rlp != NULL; rlp = rlp->next) {
    if (rlp->reg == rp)
      found = 1;
  }

  if (!found)
    return 0;

  return 1;
}

/*****************************************************************
wall_belongs_to_all_regions_in_region_list:
  In: wall
      region_list
  Out: 1 if wall belongs to all regions in the region list
       0 otherwise.
  Note: Wall can belong to several regions simultaneously.
******************************************************************/
int wall_belongs_to_all_regions_in_region_list(struct wall *this_wall,
                                               struct region_list *rlp_head) {
  struct region_list *rlp;
  struct region *rp;

  if (rlp_head == NULL)
    return 0;

  for (rlp = rlp_head; rlp != NULL; rlp = rlp->next) {
    rp = rlp->reg;

    if (!get_bit(rp->membership, this_wall->side))
      return 0;
  }

  return 1;
}

/*****************************************************************
wall_belongs_to_any_region_in_region_list:
  In: wall
      region_list
  Out: 1 if wall belongs to any region in the region list
       0 otherwise.
  Note: Wall can be belong to several regions simultaneously.
  Note: It is assumed that both wall and region list are defined for
        the same object.
******************************************************************/
int wall_belongs_to_any_region_in_region_list(struct wall *this_wall,
                                              struct region_list *rlp_head) {
  struct region_list *rlp;
  struct region *rp;

  if (rlp_head == NULL)
    return 0;

  for (rlp = rlp_head; rlp != NULL; rlp = rlp->next) {
    rp = rlp->reg;

    if (get_bit(rp->membership, this_wall->side))
      return 1;
  }

  return 0;
}

/*********************************************************************
* find_wall_center:
* In: wall
*     wall center (return value)
* Out: the center of the wall is found and assigned to vector "center"
**********************************************************************/
void find_wall_center(struct wall *w, struct vector3 *center)
{
  if(center == NULL) {
    mcell_internal_error("Error in function 'find_wall_center()'.");
  }

  center->x = (w->vert[0]->x + w->vert[1]->x + w->vert[2]->x) / 3;
  center->y = (w->vert[0]->y + w->vert[1]->y + w->vert[2]->y) / 3;
  center->z = (w->vert[0]->z + w->vert[1]->z + w->vert[2]->z) / 3;

}
