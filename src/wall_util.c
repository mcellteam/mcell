/**************************************************************************\
 ** File: wall_util.c                                                    **
 **                                                                      **
 ** Purpose: Build walls and surfaces, create edges from vertices and    **
 **    polygons.  All wall elements are assumed to be triangles.         **
 **                                                                      **
\**************************************************************************/


#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#include "rng.h"
#include "logging.h"
#include "vector.h"
#include "util.h"
#include "sym_table.h"
#include "mem_util.h"
#include "vol_util.h"
#include "mcell_structs.h"
#include "react_output.h"
#include "mdlparse_util.h"
#include "grid_util.h"
#include "count_util.h"
#include "wall_util.h"
#include "macromolecule.h"
#include <float.h>

extern struct volume *world;

/**************************************************************************\
 ** Internal utility function section--max/min stuff                     **
\**************************************************************************/

/**** This function is used only in this file.  It picks out the   ****/
/**** largest (absolute) value found among two vectors (useful for ****/
/**** properly handling floating-point rounding error).            ****/

static inline double abs_max_2vec(struct vector3 *v1,struct vector3 *v2)
{
  return max2d(max3d(fabs(v1->x), fabs(v1->y), fabs(v1->z)),
               max3d(fabs(v2->x), fabs(v2->y), fabs(v2->z)));
}

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

int edge_equals(struct poly_edge *e1,struct poly_edge *e2)
{
  if ( (e1->v1x == e2->v1x) && (e1->v1y == e2->v1y) && (e1->v1z == e2->v1z) &&
       (e1->v2x == e2->v2x) && (e1->v2y == e2->v2y) && (e1->v2z == e2->v2z) )
  {
    return 1;
  }
  if ( (e1->v1x == e2->v2x) && (e1->v1y == e2->v2y) && (e1->v1z == e2->v2z) &&
       (e1->v2x == e2->v1x) && (e1->v2y == e2->v1y) && (e1->v2z == e2->v1z) )
  {
    return 1;
  }
  return 0;
}


/***************************************************************************
edge_hash:
  In: pointer to a poly_edge struct
      number of keys in the hash table
  Out: Returns a hash value between 0 and nkeys-1.
  Note: Orientation invariant, so a hash with the two endpoints swapped
        will be the same.
***************************************************************************/

int edge_hash (struct poly_edge *pe,int nkeys)
{
  /* Get hash of X,Y,Z set of doubles for 1st and 2nd points */
  /* (Assume they're laid out consecutively in memory) */
  unsigned int hashL = jenkins_hash( (ub1*) &(pe->v1x) , 3*sizeof(double) );
  unsigned int hashR = jenkins_hash( (ub1*) &(pe->v2x) , 3*sizeof(double) );
  
  return (hashL+hashR)%nkeys;       /* ^ is symmetric so doesn't matter which is L and which is R */
}

/***************************************************************************
ehtable_init:
  In: pointer to an edge_hashtable struct
      number of keys that the hash table uses
  Out: Returns 0 on success, 1 on failure.
       Hash table is initialized.
***************************************************************************/

int ehtable_init(struct edge_hashtable *eht,int nkeys)
{
  int i;
  
  no_printf("Using %d keys to find edges.\n",nkeys);
  eht->nkeys = nkeys;
  eht->stored = 0;
  eht->distinct = 0;
  eht->data = CHECKED_MALLOC_ARRAY_NODIE(struct poly_edge,
                                         nkeys,
                                         "edge hash table");
  if (eht->data == NULL) return 1;
  
  for (i=0;i<nkeys;i++)
  {
    eht->data[i].next = NULL;
    eht->data[i].n = 0;
    eht->data[i].face1 = eht->data[i].face2 = -1;
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
int ehtable_add(struct edge_hashtable *eht,struct poly_edge *pe)
{
  int i;
  struct poly_edge *pep,*pei;
  
  i = edge_hash( pe , eht->nkeys );
  pep = &(eht->data[i]);
  
  while (pep != NULL)
  {
    if (pep->n==0)   /* New entry */
    {
      pep->n = 1;
      pep->face1 = pe->face1;
      pep->edge1 = pe->edge1;
      pep->v1x = pe->v1x; pep->v1y = pe->v1y; pep->v1z = pe->v1z;
      pep->v2x = pe->v2x; pep->v2y = pe->v2y; pep->v2z = pe->v2z;
      eht->stored++;
      eht->distinct++;
      return 0;
    }
    
    if (edge_equals(pep,pe))  /* This edge exists already ... */
    {
      if (pep->face2 == -1)   /* ...and we're the 2nd one */
      {
        pep->face2 = pe->face1;
        pep->edge2 = pe->edge1;
        pep->n++;
        eht->stored++;
        return 0;
      }
      else                    /* ...or we're 3rd and need more space */
      {
        if (pep->next != NULL)
        {
          if (edge_equals(pep->next,pe))  /* Space already there */
          {
            pep->n++;
            pep = pep->next;
            continue;                     /* Use space on next loop */
          }
        }
        
        pei = CHECKED_MALLOC_STRUCT_NODIE(struct poly_edge,
                                          "polygon edge");
        if (pei==NULL) return 1;

        pep->n++;
        pei->next = pep->next;
        pep->next = pei;
        pei->n = 0;
        pei->face1 = -1;
        pei->face2 = -1;
        pei->edge1 = -1;
        pei->edge2 = -1;
        pep = pei;
        eht->distinct--;  /* Not really distinct, just need more space */
      }
    }

    else if (pep->next != NULL) pep = pep->next;

    else  /* Hit end of list, so make space for use next loop. */
    {
      pei = CHECKED_MALLOC_STRUCT_NODIE(struct poly_edge,
                                        "polygon edge");
      if (pei==NULL) return 1;
      pei->next = pep->next;
      pep->next = pei;
      pei->n = 0;
      pei->face1 = -1;
      pei->face2 = -1;
      pei->edge1 = -1;
      pei->edge2 = -1;
      pep = pei;
    }
  }
  
  return 0;
}


/***************************************************************************
ehtable_kill:
  In: pointer to an edge_hashtable struct
  Out: No return value.  Hashtable data is deallocated.
  Note: eht itself is not freed, since it isn't created with ehtable_init.
***************************************************************************/
void ehtable_kill(struct edge_hashtable *eht)
{
  struct poly_edge *pe;
  int i;

  for (i=0;i<eht->nkeys;i++)
  {
    while (eht->data[i].next != NULL)
    {
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
static int compatible_edges(struct wall **faces,int wA,int eA,int wB,int eB)
{
  struct vector3 *vA0,*vA1,*vA2,*vB0,*vB1,*vB2;

  if((wA < 0) || (eA < 0) || (wB < 0) || (eB < 0)) return 0;

  vA0 = faces[wA]->vert[eA];
  if (eA==2) vA1 = faces[wA]->vert[0];
  else vA1 = faces[wA]->vert[eA+1];
  if (eA==0) vA2 = faces[wA]->vert[2];
  else vA2 = faces[wA]->vert[eA-1];

  vB0 = faces[wB]->vert[eB];
  if (eB==2) vB1 = faces[wB]->vert[0];
  else vB1 = faces[wB]->vert[eB+1];
  if (eB==0) vB2 = faces[wB]->vert[2];
  else vB2 = faces[wB]->vert[eB-1];

  return ( (vA0==vB1 && vA1==vB0 && !(vA2==vB2) ) ||
           ( vA0->x==vB1->x && vA0->y==vB1->y && vA0->z==vB1->z &&
             vA1->x==vB0->x && vA1->y==vB0->y && vA1->z==vB0->z &&
             !(vA2->x==vB2->x && vA2->y==vB2->y && vA2->z==vB2->z)
           )
         );
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
static void refine_edge_pairs(struct poly_edge *p,struct wall **faces)
{
#define TSWAP(x,y) temp=(x); (x)=(y); (y)=temp
  struct poly_edge *p1,*p2,*best_p1,*best_p2;
  int n1,n2,best_n1,best_n2;
  double align,best_align;
  int wA,wB,eA,eB;
  int temp;

  best_align = 2;
  best_p1 = best_p2 = p;
  best_n1 = 1;
  best_n2 = 2;

  p1 = p;
  n1 = 1;
  while (p1 != NULL && p1->n >= n1)
  {
    if (n1==1)
    {
      wA = p1->face1;
      eA = p1->edge1;
    }
    else
    {
      wA = p1->face2;
      eA = p1->edge2;
    }

    if (n1==1)
    {
      n2 = n1+1;
      p2 = p1;
    }
    else
    {
      n2 = 1;
      p2 = p1->next;
    }
    while (p2 != NULL && p2->n >= n2)
    {
      if (n2==1)
      {
        wB = p2->face1;
        eB = p2->edge1;
      }
      else
      {
        wB = p2->face2;
        eB = p2->edge2;
      }

      if (compatible_edges(faces,wA,eA,wB,eB))
      {
        align = faces[wA]->normal.x * faces[wB]->normal.x +
                faces[wA]->normal.y * faces[wB]->normal.y +
                faces[wA]->normal.z * faces[wB]->normal.z;

        if (align < best_align)
        {
          best_p1 = p1;
          best_p2 = p2;
          best_n1 = n1;
          best_n2 = n2;
          best_align = align;
        }
      }

      if (n2==1) n2++;
      else
      {
        p2=p2->next;
        n2=1;
      }
    }

    if (n1==1) n1++;
    else
    {
      p1=p1->next;
      n1=1;
    }
  }

  /* Now lots of boring logic to swap the values into the right spots.  Yawn. */

  if (best_align > 1.0) return;  /* No good pairs. */

  if (best_p1 == best_p2)
  {
    if (best_p1==p) return;  /* Best pair is already first */

    TSWAP(best_p1->face1,p->face1);
    TSWAP(best_p1->face2,p->face2);
    TSWAP(best_p1->edge1,p->edge1);
    TSWAP(best_p1->edge2,p->edge2);

    return;
  }

  if (best_p1==p)
  {
    if (best_n1==1)
    {
      if (best_n2==1)
      {
        TSWAP(best_p2->face1,p->face2);
        TSWAP(best_p2->edge1,p->edge2);
      }
      else
      {
        TSWAP(best_p2->face2,p->face2);
        TSWAP(best_p2->edge2,p->edge2);
      }
    }
    else
    {
      if (best_n2==1)
      {
        TSWAP(best_p2->face1,p->face1);
        TSWAP(best_p2->edge1,p->edge1);
      }
      else
      {
        TSWAP(best_p2->face2,p->face1);
        TSWAP(best_p2->edge2,p->edge1);
      }
    }
  }
  else if (best_p2==p)
  {
    if (best_n1==1)
    {
      if (best_n2==1)
      {
        TSWAP(best_p1->face1,p->face2);
        TSWAP(best_p1->edge1,p->edge2);
      }
      else
      {
        TSWAP(best_p1->face2,p->face2);
        TSWAP(best_p1->edge2,p->edge2);
      }
    }
    else
    {
      if (best_n2==1)
      {
        TSWAP(best_p1->face1,p->face1);
        TSWAP(best_p1->edge1,p->edge1);
      }
      else
      {
        TSWAP(best_p1->face2,p->face1);
        TSWAP(best_p1->edge2,p->edge1);
      }
    }
  }
  else
  {
    if (best_n1==1)
    {
      TSWAP(best_p1->face1,p->face1);
      TSWAP(best_p1->edge1,p->edge1);
    }
    else
    {
      TSWAP(best_p1->face2,p->face1);
      TSWAP(best_p1->edge2,p->edge1);
    }
    if (best_n2==1)
    {
      TSWAP(best_p2->face1,p->face2);
      TSWAP(best_p2->edge1,p->edge2);
    }
    else
    {
      TSWAP(best_p2->face2,p->face2);
      TSWAP(best_p2->edge2,p->edge2);
    }
  }
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

int surface_net( struct wall **facelist, int nfaces )
{
  struct poly_edge pe,*pep;
  struct edge *e;
  struct edge_hashtable eht;
  int i,j,k;
  int nedge;
  int nkeys;
  int is_closed = 1;
  
  nkeys = (3*nfaces)/2;
  if ( ehtable_init(&eht,nkeys) ) return 1;

  for (i=0;i<nfaces;i++)
  {
    if (facelist[i]==NULL) continue;

    nedge = 3;
    for (j=0;j<nedge;j++)
    {
      
      if (j+1 < nedge) k = j+1;
      else k = 0;
      
      pe.v1x = facelist[i]->vert[j]->x;
      pe.v1y = facelist[i]->vert[j]->y;
      pe.v1z = facelist[i]->vert[j]->z;
      pe.v2x = facelist[i]->vert[k]->x;
      pe.v2y = facelist[i]->vert[k]->y;
      pe.v2z = facelist[i]->vert[k]->z;
      pe.face1 = i;
      pe.edge1 = j;
      
      if ( ehtable_add(&eht,&pe) ) return 1;
    }
  }
  
  for (i=0;i<nkeys;i++)
  {
    pep = (eht.data + i);
    while (pep!=NULL)
    {
      if (pep->n > 2)
      {
        no_printf("Edge with more than two faces attached! Refining.\n");
        refine_edge_pairs(pep,facelist);
      }
      if (pep->n >= 2)
      {
        if (pep->face1 != -1 && pep->face2 != -1)
        { 
              if(compatible_edges(facelist,pep->face1,pep->edge1,pep->face2,pep->edge2))
              {
          	facelist[pep->face1]->nb_walls[pep->edge1] = facelist[pep->face2];
          	facelist[pep->face2]->nb_walls[pep->edge2] = facelist[pep->face1];
          	e = (struct edge*) CHECKED_MEM_GET_NODIE( facelist[pep->face1]->birthplace->join, "edge" );
          	if (e==NULL) return 1;

          	e->forward = facelist[pep->face1];
          	e->backward = facelist[pep->face2];
          	init_edge_transform(e,pep->edge1);
          	facelist[pep->face1]->edges[pep->edge1] = e;
          	facelist[pep->face2]->edges[pep->edge2] = e;
          	no_printf("  Edge: %d on %d and %d on %d\n",pep->edge1,pep->face1,pep->edge2,pep->face2);
              }

        }
        else is_closed = 0;
      }
      else if (pep->n==1)
      {
        is_closed = 0;
        e = (struct edge*) CHECKED_MEM_GET_NODIE( facelist[pep->face1]->birthplace->join, "edge" );
        if (e==NULL) return 1;

        e->forward = facelist[pep->face1];
        e->backward = NULL;
        /* Don't call init_edge_transform unless both edges are set */
        facelist[pep->face1]->edges[pep->edge1] = e;
        no_printf("  Edge: %d on %d\n",pep->edge1,pep->face1);
      }
      pep = pep->next;
    }
  }
  
  ehtable_kill(&eht);
  return -is_closed;  /* We use 1 to indicate malloc failure so return 0/-1 */

}


/***************************************************************************
init_edge_transform
  In: pointer to an edge
      integer telling which edge (0-2) of the "forward" face we are
  Out: No return value.  Coordinate transform in edge struct is set.
  Note: Don't call this on a non-shared edge.
***************************************************************************/

void init_edge_transform(struct edge *e,int edgenum)
{
  struct vector2 O_f,O_b;
  struct vector2 ehat_f,ehat_b;
  struct vector2 fhat_f,fhat_b;
  struct wall *wf,*wb;
  struct vector2 temp;
  struct vector3 temp3d;
  double d;
  int i,j;
  double mtx[2][2];
  struct vector2 q;
  
  wf = e->forward;
  wb = e->backward;
  i=edgenum;
  j=i+1;
  if (j>2) j=0;
  
  /* Intermediate basis from the perspective of the forward frame */
    
  temp3d.x = wf->vert[i]->x - wf->vert[0]->x;
  temp3d.y = wf->vert[i]->y - wf->vert[0]->y;
  temp3d.z = wf->vert[i]->z - wf->vert[0]->z;
  O_f.u = dot_prod(&temp3d,&(wf->unit_u));
  O_f.v = dot_prod(&temp3d,&(wf->unit_v));                        /* Origin */
  
  temp3d.x = wf->vert[j]->x - wf->vert[0]->x;
  temp3d.y = wf->vert[j]->y - wf->vert[0]->y;
  temp3d.z = wf->vert[j]->z - wf->vert[0]->z;
  temp.u = dot_prod(&temp3d,&(wf->unit_u)) - O_f.u;
  temp.v = dot_prod(&temp3d,&(wf->unit_v)) - O_f.v;        /* Far side of e */
  
  d = 1.0/sqrt(temp.u*temp.u + temp.v*temp.v);
  ehat_f.u = temp.u*d;
  ehat_f.v = temp.v*d;                                   /* ehat along edge */
  fhat_f.u = -ehat_f.v;
  fhat_f.v = ehat_f.u;                               /* fhat 90 degrees CCW */
  
  /* Intermediate basis from the perspective of the backward frame */
  
  temp3d.x = wf->vert[i]->x - wb->vert[0]->x;
  temp3d.y = wf->vert[i]->y - wb->vert[0]->y;
  temp3d.z = wf->vert[i]->z - wb->vert[0]->z;
  O_b.u = dot_prod(&temp3d,&(wb->unit_u));
  O_b.v = dot_prod(&temp3d,&(wb->unit_v));                        /* Origin */
  
  temp3d.x = wf->vert[j]->x - wb->vert[0]->x;
  temp3d.y = wf->vert[j]->y - wb->vert[0]->y;
  temp3d.z = wf->vert[j]->z - wb->vert[0]->z;
  temp.u = dot_prod(&temp3d,&(wb->unit_u)) - O_b.u;
  temp.v = dot_prod(&temp3d,&(wb->unit_v)) - O_b.v;        /* Far side of e */
  
  d = 1.0/sqrt(temp.u*temp.u + temp.v*temp.v);
  ehat_b.u = temp.u*d;
  ehat_b.v = temp.v*d;                                   /* ehat along edge */
  fhat_b.u = -ehat_b.v;
  fhat_b.v = ehat_b.u;                               /* fhat 90 degrees CCW */
  
  /* Calculate transformation matrix */
  
  mtx[0][0] = ehat_f.u*ehat_b.u + fhat_f.u*fhat_b.u;
  mtx[0][1] = ehat_f.v*ehat_b.u + fhat_f.v*fhat_b.u;
  mtx[1][0] = ehat_f.u*ehat_b.v + fhat_f.u*fhat_b.v;
  mtx[1][1] = ehat_f.v*ehat_b.v + fhat_f.v*fhat_b.v;
  
  /* Calculate translation vector */
  
  q.u = O_b.u;
  q.v = O_b.v;
  
  q.u -= mtx[0][0]*O_f.u + mtx[0][1]*O_f.v;
  q.v -= mtx[1][0]*O_f.u + mtx[1][1]*O_f.v;
  
  /* Store the results */
  
  e->cos_theta = mtx[0][0];
  e->sin_theta = mtx[0][1];
  e->translate.u = q.u;
  e->translate.v = q.v;
}



/***************************************************************************
sharpen_object:
  In: pointer to an object
  Out: 0 on success, 1 on failure.
       Adds edges to the object and all its children.
***************************************************************************/

int sharpen_object(struct object *parent)
{
  struct object *o;
  int i;
 
  if (parent->object_type == POLY_OBJ || parent->object_type == BOX_OBJ)
  {
    i = surface_net(parent->wall_p , parent->n_walls);
    if (i == 1)
      mcell_allocfailed("Failed to connect walls of object %s along shared edges.",
                        parent->sym->name);
  }
  else if (parent->object_type == META_OBJ)
  {
    for ( o = parent->first_child ; o != NULL ; o = o->next )
    {
      if ( sharpen_object( o ) ) return 1;
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

int sharpen_world(void)
{
  struct object *o;
  
  for (o = world->root_instance ; o != NULL ; o = o->next)
  {
    if (sharpen_object(o)) return 1;
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

double closest_interior_point(struct vector3 *pt,struct wall *w,struct vector2 *ip,double r2)
{
  UNUSED(r2);

  struct vector3 v;
  double a1,a2;
  
  closest_pt_point_triangle(pt , w->vert[0] , w->vert[1] , w->vert[2] , &v);
  xyz2uv(&v,w,ip);
  
  /* Check to see if we're lying on an edge; if so, scoot towards centroid. */
  /* ip lies on edge of wall if cross products are zero */

  a1 = ip->u*w->uv_vert2.v-ip->v*w->uv_vert2.u;
  a2 = w->uv_vert1_u*ip->v;
  while (!distinguishable(ip->v,0,EPS_C) ||
         !distinguishable(a1,0,EPS_C) ||
	 !distinguishable(a1+a2,2.0*w->area,EPS_C) )
  {
    /* Need to move centrally by a fraction larger than EPS_C or we'll have to do this many times! */
    ip->u = (1.0-5*EPS_C)*ip->u + 5*EPS_C*0.333333333333333*(w->uv_vert1_u+w->uv_vert2.u);
    ip->v = (1.0-5*EPS_C)*ip->v + 5*EPS_C*0.333333333333333*w->uv_vert2.v;
    
    a1 = ip->u*w->uv_vert2.v-ip->v*w->uv_vert2.u;
    a2 = w->uv_vert1_u*ip->v;
  }
  return (v.x-pt->x)*(v.x-pt->x) + (v.y-pt->y)*(v.y-pt->y) + (v.z-pt->z)*(v.z-pt->z);
}


/***************************************************************************
find_edge_point:
  In: a wall
      a point in the coordinate system of that wall where we are now
         (assumed to be on or inside triangle)
      a 2D displacement vector to move
      a place to store the coordinate on the edge, if we hit it
  Out: index of the edge we hit (0, 1, or 2), or -1 if the new location
       is within the wall, or -2 if we can't tell.  If the result is
       0, 1, or 2, edgept is set to the new location.
***************************************************************************/

int find_edge_point(struct wall *here,struct vector2 *loc,struct vector2 *disp,struct vector2 *edgept)
{
  double lxd;
  double lxc1,lxc2;
  double dxc1,dxc2;
  double f,s,t;
  
  lxd = loc->u*disp->v - loc->v*disp->u;
  
  lxc1 = -loc->v*here->uv_vert1_u;
  dxc1 = -disp->v*here->uv_vert1_u;
  
  if (dxc1 < -EPS_C || dxc1 > EPS_C)
  {
    f = 1.0/dxc1; /* f>0 is passing outwards */
    s = -lxd*f;
    if (0.0<s && s<1.0 && f>0.0)
    {
      t = -lxc1*f;
      if (EPS_C<t && t<1.0)
      {
	edgept->u = loc->u + t*disp->u;
	edgept->v = loc->v + t*disp->v;
	return 0;
      }
      else if (t > 1.0+EPS_C) return -1;
      /* else can't tell if we hit this edge, assume not */
    }
  }
  
  lxc2 = loc->u*here->uv_vert2.v - loc->v*here->uv_vert2.u;
  dxc2 = disp->u*here->uv_vert2.v - disp->v*here->uv_vert2.u;
  
  if (dxc2 < -EPS_C || dxc2 > EPS_C)
  {
    f = 1.0/dxc2; /* f<0 is passing outwards */
    s = 1.0 + lxd*f;
    if (0.0<s && s<1.0 && f<0.0)
    {
      t = -lxc2*f;
      if (EPS_C<t && t<1.0)
      {
	edgept->u = loc->u + t*disp->u;
	edgept->v = loc->v + t*disp->v;
	return 2;
      }
      else if (t > 1.0+EPS_C) return -1;
      /* else can't tell */
    }
  }
  
  f = dxc2-dxc1;
  
  if (f < -EPS_C || f > EPS_C)
  {
    f = 1.0/f; /* f>0 is passing outwards */
    s = -(lxd + dxc1)*f;
    if (0.0<s && s<1.0 && f>0.0)
    {
      t = (here->uv_vert1_u*here->uv_vert2.v + lxc1 - lxc2) * f;
      if (EPS_C<t && t<1.0)
      {
	edgept->u = loc->u + t*disp->u;
	edgept->v = loc->v + t*disp->v;
	return 1;
      }
      else if (t > 1.0+EPS_C) return -1;
      /* else can't tell */
    }
  }
  
  return -2;  /* Couldn't tell whether we hit or not--calling function should pick another displacement */
}


/***************************************************************************
traverse_surface:
  In: a wall
      a point in the coordinate system of that wall
      which edge to travel off of
      a vector to set for the new wall
  Out: NULL if the edge is not shared, or a pointer to the wall in that direction
       if it is shared.  newloc is set to loc in the coordinate system of the new
       wall (after flattening the walls along their shared edge)
***************************************************************************/

struct wall* traverse_surface(struct wall *here,struct vector2 *loc,int which,struct vector2 *newloc)
{
  struct edge *e;
  struct wall *there;
  double u,v;
  
  e = here->edges[which];

  if (e==NULL) return NULL;
  
  if (e->forward == here)
  {
    /* Apply forward transform to loc */
    there = e->backward;
    
    /* rotation */
    u =  e->cos_theta*loc->u + e->sin_theta*loc->v;
    v = -e->sin_theta*loc->u + e->cos_theta*loc->v;
    
    /* translation */
    newloc->u = u + e->translate.u;
    newloc->v = v + e->translate.v;

    return there;
  }
  else
  {
    /* Apply inverse transform to loc */
    there = e->forward;
    
    /* inverse translation */
    u = loc->u - e->translate.u;
    v = loc->v - e->translate.v;
    
    /* inverse rotation */
    newloc->u =  e->cos_theta*u - e->sin_theta*v;
    newloc->v =  e->sin_theta*u + e->cos_theta*v;
    
    return there;
  }
}


/***************************************************************************
is_manifold:
  In: a region.  This region must already be painted on walls.  The edges
      must have already been added to the object (i.e. sharpened).
  Out: 1 if the region is a manifold, 0 otherwise.
  Note: by "manifold" we mean "orientable compact two-dimensional
        manifold without boundaries embedded in R3"
***************************************************************************/

int is_manifold(struct region *r)
{
  struct wall **wall_array = NULL, *w = NULL;
  struct region_list *rl = NULL;

  wall_array = r->parent->wall_p;
  
  if (wall_array == NULL)
  {
    mcell_internal_error("Region '%s' has NULL wall array!", r->sym->name);
    return 0;
  }
   
  for (int n_wall=0; n_wall<r->parent->n_walls; n_wall++)
  {
    if (!get_bit(r->membership, n_wall)) continue;  /* Skip wall not in region */
    w = wall_array[n_wall];
    for (int nb=0; nb<2; nb++)
    {
      if (w->nb_walls[nb] == NULL)
      {
        mcell_error_nodie("BARE EDGE on wall %u edge %d.", n_wall, nb);
	return 0; /* Bare edge--not a manifold */
      }
      
      for (rl = w->nb_walls[nb]->counting_regions ; rl != NULL ; rl = rl->next)
      {
	if (rl->reg == r) break;
      }
      if (rl==NULL)
      {
	mcell_error_nodie("Wall %u edge %d leaves region!", n_wall, nb);
	return 0;  /* Can leave region--not a manifold */
      }
    }
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

void jump_away_line(struct vector3 *p,struct vector3 *v,double k,
                    struct vector3 *A,struct vector3 *B,struct vector3 *n)
{
  struct vector3 e,f;
  double le_1,tiny;
  
  e.x = B->x - A->x;
  e.y = B->y - A->y;
  e.z = B->z - A->z;
  
  le_1 = 1.0/sqrt(e.x*e.x + e.y*e.y + e.z*e.z);
  
  e.x *= le_1;
  e.y *= le_1;
  e.z *= le_1;

  f.x = n->y*e.z - n->z*e.y;
  f.y = n->z*e.x - n->x*e.z;
  f.z = n->x*e.y - n->y*e.x;
  
  tiny = EPS_C * (abs_max_2vec(p,v) + 1.0) / (k * max3d(fabs(f.x),fabs(f.y),fabs(f.z)));
  if ( (rng_uint(world->rng) & 1) == 0 ) {
     tiny = -tiny;
  }
  v->x -= tiny*f.x;
  v->y -= tiny*f.y;
  v->z -= tiny*f.z;
}


/***************************************************************************
touch_wall:
  In: starting coordinate
      vector to move along (forwards and backwards)
      wall we're checking for a collision
  Out: Double value between -1.0 and 1.0 if the movement ray intersected
         the wall within range.  Returns 1.0 if out of range or missed.
  Note: This code is used to estimate probabilities in constrained spaces.
        Use collide_wall to detect collisions between molecules and
        surfaces.
***************************************************************************/

double touch_wall(struct vector3 *point,struct vector3 *move,struct wall *face)
{
  double dp,dv,dd;
  double nx,ny,nz;
  double b,c,t;
  double f,g,h;
  struct vector3 local;
  
  nx = face->normal.x;
  ny = face->normal.y;
  nz = face->normal.z;
  
  dp = nx*point->x + ny*point->y + nz*point->z;
  dv = nx*move->x + ny*move->y + nz*move->z;
  dd = dp - face->d;

  if (dd==0.0 || dd*dd >= dv*dv) return 1.0;

  t = -dd/dv;
  
  local.x = point->x + t*move->x - face->vert[0]->x;
  local.y = point->y + t*move->y - face->vert[0]->y;
  local.z = point->z + t*move->z - face->vert[0]->z;
  
  b = local.x*face->unit_u.x + local.y*face->unit_u.y + local.z*face->unit_u.z;
  c = local.x*face->unit_v.x + local.y*face->unit_v.y + local.z*face->unit_v.z;
  
  if (face->uv_vert2.v < 0.0)
  {
    c = -c;
    f = -face->uv_vert2.v;
  }
  else f = face->uv_vert2.v;
    
  if (c > 0)
  {
    g = b*f;
    h = c*face->uv_vert2.u;
    if (g > h)
    {
      if ( c*face->uv_vert1_u + g < h + face->uv_vert1_u*face->uv_vert2.v ) return t;
    }
  }
  
  return 1.0;  
}


/***************************************************************************
collide_wall:
  In: starting coordinate
      vector to move along
      wall we're checking for a collision
      double to store time of collision
      vector to store the location of the collision
      flag to signal whether we should modify the movement vector in an
        ambiguous case (i.e. if we hit an edge or corner); if not, any
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

int collide_wall(struct vector3 *point,struct vector3 *move,struct wall *face,
                 double *t,struct vector3 *hitpt,int update_move)
{
  double dp,dv,dd;
  double nx,ny,nz;
  double a,b,c;
  double f,g,h;
  double d_eps;
  struct vector3 local;

  if(world->notify->final_summary == NOTIFY_FULL){
      world->ray_polygon_tests++;
  }
  
  nx = face->normal.x;
  ny = face->normal.y;
  nz = face->normal.z;
  
  dp = nx*point->x + ny*point->y + nz*point->z;
  dv = nx*move->x + ny*move->y + nz*move->z;
  dd = dp - face->d;

  if (dd > 0.0)
  {
    d_eps = EPS_C;
    if (dd < d_eps) d_eps = 0.5*dd;

    /* Start & end above plane */
    if (dd+dv > d_eps) return COLLIDE_MISS;
  }
  else
  {
    d_eps = -EPS_C;
    if (dd > d_eps) d_eps = 0.5*dd;

    /* Start & end below plane */
    if (dd<0.0 && dd+dv<d_eps) return COLLIDE_MISS;
  }
  
  if (dd==0.0)
  {
    /* Start beside plane, end above or below */
    if (dv != 0.0) return COLLIDE_MISS;

    if (update_move)
    {
      a = (abs_max_2vec( point , move ) + 1.0) * EPS_C;
      if ((rng_uint(world->rng)&1)==0) a = -a;
      if (dd==0.0)
      {
        move->x -= a*nx;
        move->y -= a*ny;
        move->z -= a*nz;
      }
      else
      {
        move->x *= (1.0-a);
        move->y *= (1.0-a);
        move->z *= (1.0-a);
      }
      return COLLIDE_REDO;
    }
    else return COLLIDE_MISS;
  }
  
  a = 1.0/dv;
  a *= -dd;         /* Time we actually hit */
  *t = a;
  
  hitpt->x = point->x + a*move->x;
  hitpt->y = point->y + a*move->y;
  hitpt->z = point->z + a*move->z;
  
  local.x = hitpt->x - face->vert[0]->x;
  local.y = hitpt->y - face->vert[0]->y;
  local.z = hitpt->z - face->vert[0]->z;
  
  b = local.x*face->unit_u.x + local.y*face->unit_u.y + local.z*face->unit_u.z;
  c = local.x*face->unit_v.x + local.y*face->unit_v.y + local.z*face->unit_v.z;
  
  if (face->uv_vert2.v < 0.0)
  {
    c = -c;
    f = -face->uv_vert2.v;
  }
  else f = face->uv_vert2.v;
    
  if (c > 0)
  {
    g = b*f;
    h = c*face->uv_vert2.u;
    if (g > h)
    {
      if ( c*face->uv_vert1_u + g < h + face->uv_vert1_u*face->uv_vert2.v )
      {
        if (dv>0) return COLLIDE_BACK;
        else return COLLIDE_FRONT;
      }
      else if (c*face->uv_vert1_u + g == h + face->uv_vert1_u*face->uv_vert2.v)
      {
        if (update_move)
        {
          jump_away_line(point,move,a,face->vert[1],face->vert[2],&(face->normal));
          return COLLIDE_REDO;
        }
        else return COLLIDE_MISS;
      }
      else return COLLIDE_MISS;
    }
    else if (g == h)
    {
      if (update_move)
      {
        jump_away_line(point,move,a,face->vert[2],face->vert[0],&(face->normal));
        return COLLIDE_REDO;
      }
      else return COLLIDE_MISS;
    }
    else return COLLIDE_MISS;
  }
  else if (c == 0) /* Hit first edge! */
  {
    if (update_move)
    {
      jump_away_line(point,move,a,face->vert[0],face->vert[1],&(face->normal));
      return COLLIDE_REDO;
    }
    else return COLLIDE_MISS;
  }
  else return COLLIDE_MISS;
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
         COLLIDE_MOL_M  hit
  Note: t and/or hitpt may be modified even if there is no collision
        Not highly optimized yet.
***************************************************************************/
int collide_mol(struct vector3 *point,struct vector3 *move,
                struct abstract_molecule *a,double *t,struct vector3 *hitpt)
{
  struct vector3 dir; /* From starting point of moving molecule to target */
  struct vector3 *pos; /* Position of target molecule */
  
  double movelen2; /* Square of distance the moving molecule travels */
  double dirlen2;  /* Square of distance between moving and target molecules */
  double d;        /* Dot product of movement vector and vector to target */
  double sigma2;   /* Square of interaction radius */
  
  if ((a->properties->flags & ON_GRID)!=0) return COLLIDE_MISS; /* Should never call on grid molecule! */
  
  pos = &( ((struct volume_molecule*)a)->pos );
  
  sigma2 = world->rx_radius_3d*world->rx_radius_3d; 

  dir.x = pos->x - point->x;
  dir.y = pos->y - point->y;
  dir.z = pos->z - point->z;

  d = dir.x*move->x + dir.y*move->y + dir.z*move->z;
  
  /* Miss the molecule if it's behind us */
  if (d<0) return COLLIDE_MISS; 
  
  movelen2 = move->x*move->x + move->y*move->y + move->z*move->z;

  /* check whether the test molecule is futher than the displacement. */
  if (d > movelen2) return COLLIDE_MISS;
  
  dirlen2 = dir.x*dir.x + dir.y*dir.y + dir.z*dir.z;
  
  /* check whether the moving molecule will miss interaction disk of the
     test molecule.*/
  if (movelen2*dirlen2 - d*d > movelen2*sigma2) return COLLIDE_MISS;

  *t = d/movelen2;
//  *t = d/sqrt(movelen2*dirlen2);
  hitpt->x = point->x + (*t)*move->x;  
  hitpt->y = point->y + (*t)*move->y;  
  hitpt->z = point->z + (*t)*move->z;  
  return COLLIDE_MOL_M;
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
static int wall_in_box(struct vector3 **vert,struct vector3 *normal,
                double d,struct vector3 *b0,struct vector3 *b1)
{
#define n_vert 3
  int temp;
  int i,j,k;
  struct vector3 *v1,*v2;
  struct vector3 n,u,v;
  struct vector3 ba,bb,c;
  double r,a1,a2,a3,a4,cu,cv;
  double vu_[6]; /* Assume wall has 3 vertices */
  double *vv_;
  int v_set;
  double d_box[8];
  int n_opposite;
  
/* Lookup table for vertex-edge mapping for a cube */
  int which_x1[12] = {0,0,0,0,1,1,1,1,0,0,0,1};
  int which_y1[12] = {0,0,1,1,1,1,0,0,0,0,1,0};
  int which_z1[12] = {0,1,1,0,0,1,1,0,0,1,1,0};
  int which_x2[12] = {0,0,0,1,1,1,1,0,0,1,1,0};
  int which_y2[12] = {0,1,1,1,1,0,0,0,1,0,1,1};
  int which_z2[12] = {1,1,0,0,1,1,0,0,0,1,1,0};
  
  int edge1_vt[12] = {0,1,3,2,6,7,5,4,0,1,3,4};
  int edge2_vt[12] = {1,3,2,6,7,5,4,0,2,5,7,2};
  
/* Check if any vertex of the wall is in the box. */
  for (i=0;i<n_vert;i++)
  {
    v2 = vert[i];
    if (v2->x >= b0->x && v2->x <= b1->x && 
        v2->y >= b0->y && v2->y <= b1->y &&
        v2->z >= b0->z && v2->z <= b1->z) return 1;
  }
  
  
/* Check if any wall edge intersects any face of the box */
  for (i=0;i<n_vert;i++)
  {
    v2 = vert[i];
    v1 = (i==0) ? vert[n_vert-1] : vert[i-1];
    
/* x-faces */
    if ((v1->x <= b0->x && b0->x < v2->x) || (v1->x > b0->x && b0->x >= v2->x))
    {
      r = (b0->x - v1->x)/(v2->x - v1->x);
      a3 = v1->y + r*(v2->y - v1->y);
      a4 = v1->z + r*(v2->z - v1->z);
      if (b0->y <= a3 && a3 <= b1->y && b0->z <= a4 && a4 <= b1->z) return 2;
    }
    if ((v1->x <= b1->x && b1->x < v2->x) || (v1->x > b1->x && b1->x >= v2->x))
    {
      r = (b1->x - v1->x)/(v2->x - v1->x);
      a3 = v1->y + r*(v2->y - v1->y);
      a4 = v1->z + r*(v2->z - v1->z);
      if (b0->y <= a3 && a3 <= b1->y && b0->z <= a4 && a4 <= b1->z) return 3;
    }

/* y-faces */
    if ((v1->y <= b0->y && b0->y < v2->y) || (v1->y > b0->y && b0->y >= v2->y))
    {
      r = (b0->y - v1->y)/(v2->y - v1->y);
      a3 = v1->x + r*(v2->x - v1->x);
      a4 = v1->z + r*(v2->z - v1->z);
      if (b0->x <= a3 && a3 <= b1->x && b0->z <= a4 && a4 <= b1->z) return 4;
    }
    if ((v1->y <= b1->y && b1->y < v2->y) || (v1->y > b1->y && b1->y >= v2->y))
    {
      r = (b1->y - v1->y)/(v2->y - v1->y);
      a3 = v1->x + r*(v2->x - v1->x);
      a4 = v1->z + r*(v2->z - v1->z);
      if (b0->x <= a3 && a3 <= b1->x && b0->z <= a4 && a4 <= b1->z) return 5;
    }

/* z-faces */
    if ((v1->z <= b0->z && b0->z < v2->z) || (v1->z > b0->z && b0->z >= v2->z))
    {
      r = (b0->z - v1->z)/(v2->z - v1->z);
      a3 = v1->y + r*(v2->y - v1->y);
      a4 = v1->x + r*(v2->x - v1->x);
      if (b0->y <= a3 && a3 <= b1->y && b0->x <= a4 && a4 <= b1->x) return 6;
    }
    if ((v1->z <= b1->z && b1->z < v2->z) || (v1->z > b1->z && b1->z >= v2->z))
    {
      r = (b1->z - v1->z)/(v2->z - v1->z);
      a3 = v1->y + r*(v2->y - v1->y);
      a4 = v1->x + r*(v2->x - v1->x);
      if (b0->y <= a3 && a3 <= b1->y && b0->x <= a4 && a4 <= b1->x) return 7;
    }

  }


/* Check if any box edge intersects the wall */

  n_opposite = 0;
  vv_ = &(vu_[n_vert]);
  v_set = 0;

/* Wall coordinate system n,u,v */  
  n.x = normal->x; n.y = normal->y; n.z = normal->z;
  u.x = vert[1]->x - vert[0]->x;
  u.y = vert[1]->y - vert[0]->y;
  u.z = vert[1]->z - vert[0]->z;
  r = 1/sqrt(u.x*u.x + u.y*u.y + u.z*u.z);
  u.x *= r; u.y *=r; u.z *= r;
  v.x = n.y*u.z - n.z*u.y;
  v.y = - (n.x*u.z - n.z*u.x);
  v.z = n.x*u.y - n.y*u.x;

  
/* Test every edge. */
  bb.x = b0->x; bb.y = b0->y; bb.z = b0->z;
  d_box[0] = bb.x*n.x + bb.y*n.y + bb.z*n.z;
  for (i=0;i<12;i++)
  {
    if (i<7) /* Visiting new vertices in order */
    {
      ba.x = bb.x; ba.y = bb.y; ba.z = bb.z;
      bb.x = (which_x2[i]) ? b1->x : b0->x;
      bb.y = (which_y2[i]) ? b1->y : b0->y;
      bb.z = (which_z2[i]) ? b1->z : b0->z;
      a2 = d_box[ edge2_vt[i] ] = bb.x*n.x + bb.y*n.y + bb.z*n.z;
      a1 = d_box[ edge1_vt[i] ];
      
      if ( (a1 - d < 0 && a2 - d < 0) ||
           (a1 - d > 0 && a2 - d > 0) ) continue;
      else n_opposite++;
    }
    else /* Revisiting old vertices out of order */
    {
/*      if (!n_opposite) return 0; */
      a1 = d_box[ edge1_vt[i] ];
      a2 = d_box[ edge2_vt[i] ];

      if ( (a1 - d < 0 && a2 - d < 0) ||
           (a1 - d > 0 && a2 - d > 0) ) continue;
      
      n_opposite++;
      ba.x = (which_x1[i]) ? b1->x : b0->x;
      ba.y = (which_y1[i]) ? b1->y : b0->y;
      ba.z = (which_z1[i]) ? b1->z : b0->z;
      bb.x = (which_x2[i]) ? b1->x : b0->x;
      bb.y = (which_y2[i]) ? b1->y : b0->y;
      bb.z = (which_z2[i]) ? b1->z : b0->z;
    }
/* Now ba,bb = box edge endpoints ; a1,a2 = distances along wall normal */
    r = (d - a1)/(a2-a1);
    c.x = ba.x + r*(bb.x-ba.x);
    c.y = ba.y + r*(bb.y-ba.y);
    c.z = ba.z + r*(bb.z-ba.z);
    cu = c.x*u.x + c.y*u.y + c.z*u.z;
    cv = c.x*v.x + c.y*v.y + c.z*v.z;
    if (!v_set)
    {
      v_set=1;
      for (j=0;j<n_vert;j++)
      {
        vu_[j] = vert[j]->x*u.x + vert[j]->y*u.y + vert[j]->z*u.z;
        vv_[j] = vert[j]->x*v.x + vert[j]->y*v.y + vert[j]->z*v.z;
      }
    }
/* Test for internal intersection point in wall coordinate space */
    temp=0;    
    for (j=0;j<n_vert;j++)
    {
      k = (j==0) ? n_vert-1 : j-1;
      if ( (vu_[k] < cu && cu <= vu_[j]) ||
           (vu_[k] >= cu && cu > vu_[j]) )
      {
        r = (cu - vu_[k])/(vu_[j] - vu_[k]);
        if ( (vv_[k] + r*(vv_[j]-vv_[k])) > cv ) temp++;
      }
    }
    if (temp & 1) return 8+i;
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

void init_tri_wall(struct object *objp, int side, struct vector3 *v0, struct vector3 *v1, struct vector3 *v2)
{
  struct wall *w;            /* The wall we're working with */
  double f,fx,fy,fz;
  struct vector3 vA,vB,vX;
 
  w=&objp->walls[side];
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
  cross_prod(&vA , &vB , &vX);
  w->area = 0.5 * vect_length(&vX);

  if(w->area == 0)
  {
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
  	w->flags=0;
  	w->counting_regions = NULL;

	return;
  }

  fx = (v1->x - v0->x);
  fy = (v1->y - v0->y);
  fz = (v1->z - v0->z);
  f = 1 / sqrt( fx*fx + fy*fy + fz*fz );
  
  w->unit_u.x = fx * f;
  w->unit_u.y = fy * f;
  w->unit_u.z = fz * f;
  
  fx = (v2->x - v0->x);
  fy = (v2->y - v0->y);
  fz = (v2->z - v0->z);

  w->normal.x = w->unit_u.y * fz - w->unit_u.z * fy;
  w->normal.y = w->unit_u.z * fx - w->unit_u.x * fz;
  w->normal.z = w->unit_u.x * fy - w->unit_u.y * fx;
  f = 1 / sqrt( w->normal.x*w->normal.x + w->normal.y*w->normal.y + w->normal.z*w->normal.z );
  w->normal.x *= f;
  w->normal.y *= f;
  w->normal.z *= f;
  w->unit_v.x = w->normal.y * w->unit_u.z - w->normal.z * w->unit_u.y;
  w->unit_v.y = w->normal.z * w->unit_u.x - w->normal.x * w->unit_u.z;
  w->unit_v.z = w->normal.x * w->unit_u.y - w->normal.y * w->unit_u.x;
  w->d = v0->x * w->normal.x + v0->y * w->normal.y + v0->z * w->normal.z;
  
  w->uv_vert1_u = (w->vert[1]->x - w->vert[0]->x)*w->unit_u.x + 
                  (w->vert[1]->y - w->vert[0]->y)*w->unit_u.y +
                  (w->vert[1]->z - w->vert[0]->z)*w->unit_u.z;
  w->uv_vert2.u = (w->vert[2]->x - w->vert[0]->x)*w->unit_u.x + 
                  (w->vert[2]->y - w->vert[0]->y)*w->unit_u.y +
                  (w->vert[2]->z - w->vert[0]->z)*w->unit_u.z;
  w->uv_vert2.v = (w->vert[2]->x - w->vert[0]->x)*w->unit_v.x + 
                  (w->vert[2]->y - w->vert[0]->y)*w->unit_v.y +
                  (w->vert[2]->z - w->vert[0]->z)*w->unit_v.z;
  
  w->grid = NULL;

  w->parent_object = objp;
  w->flags=0;
  w->counting_regions = NULL;
  no_printf("Created wall %d on object %s at:\n",w->side,w->parent_object->sym->name);
  no_printf("  vertex 0: %.9g, %.9g, %.9g\n",w->vert[0]->x,w->vert[0]->y,w->vert[0]->z);
  no_printf("  vertex 1: %.9g, %.9g, %.9g\n",w->vert[1]->x,w->vert[1]->y,w->vert[1]->z);
  no_printf("  vertex 2: %.9g, %.9g, %.9g\n",w->vert[2]->x,w->vert[2]->y,w->vert[2]->z);
}


/***************************************************************************
wall_bounding_box:
  In: a wall
      vector to store one corner of the bounding box for that wall
      vector to store the opposite corner
  Out: No return value.  The vectors are set to define the smallest box
       that contains the wall.
***************************************************************************/
static void wall_bounding_box(struct wall *w , struct vector3 *llf, struct vector3 *urb)
{
  llf->x = urb->x = w->vert[0]->x;
  llf->y = urb->y = w->vert[0]->y;
  llf->z = urb->z = w->vert[0]->z;
  
  if (w->vert[1]->x < llf->x) llf->x = w->vert[1]->x;
  else if (w->vert[1]->x > urb->x) urb->x = w->vert[1]->x;
  if (w->vert[2]->x < llf->x) llf->x = w->vert[2]->x;
  else if (w->vert[2]->x > urb->x) urb->x = w->vert[2]->x;

  if (w->vert[1]->y < llf->y) llf->y = w->vert[1]->y;
  else if (w->vert[1]->y > urb->y) urb->y = w->vert[1]->y;
  if (w->vert[2]->y < llf->y) llf->y = w->vert[2]->y;
  else if (w->vert[2]->y > urb->y) urb->y = w->vert[2]->y;

  if (w->vert[1]->z < llf->z) llf->z = w->vert[1]->z;
  else if (w->vert[1]->z > urb->z) urb->z = w->vert[1]->z;
  if (w->vert[2]->z < llf->z) llf->z = w->vert[2]->z;
  else if (w->vert[2]->z > urb->z) urb->z = w->vert[2]->z;
}


/***************************************************************************
wall_to_vol:
  In: a wall
      the subvolume to which the wall belongs
  Out: The updated list of walls for that subvolume that now contains the
       wall requested.  
***************************************************************************/
struct wall_list* wall_to_vol(struct wall *w, struct subvolume *sv)
{
  struct wall_list *wl = CHECKED_MEM_GET_NODIE(sv->local_storage->list, "wall list");
  if(wl == NULL) return NULL;
  
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
struct wall* localize_wall(struct wall *w, struct storage *stor)
{
  struct wall *ww;
  ww = CHECKED_MEM_GET_NODIE(stor->face, "wall");
  if (ww==NULL) return NULL;
  
  memcpy(ww , w , sizeof(struct wall));
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
static struct wall* distribute_wall(struct wall *w)
{
  struct wall *where_am_i;            /* Version of the wall in local memory */
  struct vector3 llf,urb,cent;                      /* Bounding box for wall */
  int x_max,x_min,y_max,y_min,z_max,z_min; /* Enlarged box to avoid rounding */
  int h,i,j,k;                         /* Iteration variables for subvolumes */
  double leeway = 1.0;                                    /* Margin of error */
  
  wall_bounding_box(w,&llf,&urb);
  
  if (llf.x<-leeway) leeway=-llf.x;
  if (llf.y<-leeway) leeway=-llf.y;
  if (llf.z<-leeway) leeway=-llf.z;
  if (urb.x>leeway) leeway=urb.x;
  if (urb.y>leeway) leeway=urb.y;
  if (urb.z>leeway) leeway=urb.z;
  leeway = EPS_C + leeway*EPS_C;
  if (world->use_expanded_list) 
  {
    leeway += world->rx_radius_3d;
  }

  llf.x -= leeway;
  llf.y -= leeway;
  llf.z -= leeway;
  urb.x += leeway;
  urb.y += leeway;
  urb.z += leeway;

  cent.x = 0.33333333333*(w->vert[0]->x + w->vert[1]->x + w->vert[2]->x);
  cent.y = 0.33333333333*(w->vert[0]->y + w->vert[1]->y + w->vert[2]->y);
  cent.z = 0.33333333333*(w->vert[0]->z + w->vert[1]->z + w->vert[2]->z);
  
  x_min = bisect( world->x_partitions , world->nx_parts , llf.x );
  if (urb.x < world->x_partitions[x_min+1]) x_max = x_min+1;
  else x_max = bisect( world->x_partitions , world->nx_parts , urb.x ) + 1;

  y_min = bisect( world->y_partitions , world->ny_parts , llf.y );
  if (urb.y < world->y_partitions[y_min+1]) y_max = y_min+1;
  else y_max = bisect( world->y_partitions , world->ny_parts , urb.y ) + 1;

  z_min = bisect( world->z_partitions , world->nz_parts , llf.z );
  if (urb.z < world->z_partitions[z_min+1]) z_max = z_min+1;
  else z_max = bisect( world->z_partitions , world->nz_parts , urb.z ) + 1;
  
  if ( (z_max-z_min)*(y_max-y_min)*(x_max-x_min) == 1 )
  {
    h = z_min + (world->nz_parts - 1)*(y_min + (world->ny_parts - 1)*x_min);
    where_am_i = localize_wall( w , world->subvol[h].local_storage );
    if(where_am_i == NULL) return NULL;
     
    if (wall_to_vol( where_am_i , &(world->subvol[h]) ) == NULL) return NULL;

    return where_am_i;
  }

  for (i=x_min;i<x_max;i++) { if (cent.x < world->x_partitions[i]) break; }
  for (j=y_min;j<y_max;j++) { if (cent.y < world->y_partitions[j]) break; }
  for (k=z_min;k<z_max;k++) { if (cent.z < world->z_partitions[k]) break; }
  
  h = (k-1) + (world->nz_parts - 1)*((j-1) + (world->ny_parts - 1)*(i-1));
  where_am_i = localize_wall( w , world->subvol[h].local_storage );
  if(where_am_i == NULL) return NULL;
  
  for (k=z_min;k<z_max;k++)
  {
    for (j=y_min;j<y_max;j++)
    {
      for (i=x_min;i<x_max;i++)
      {
        h = k + (world->nz_parts - 1)*(j + (world->ny_parts - 1)*i);
        llf.x = world->x_fineparts[ world->subvol[h].llf.x ] - leeway;
        llf.y = world->y_fineparts[ world->subvol[h].llf.y ] - leeway;
        llf.z = world->z_fineparts[ world->subvol[h].llf.z ] - leeway;
        urb.x = world->x_fineparts[ world->subvol[h].urb.x ] + leeway;
        urb.y = world->y_fineparts[ world->subvol[h].urb.y ] + leeway;
        urb.z = world->z_fineparts[ world->subvol[h].urb.z ] + leeway;

        if (wall_in_box(w->vert,&(w->normal),w->d,&llf,&urb))
	{
	  if (wall_to_vol(where_am_i,&(world->subvol[h])) == NULL) return NULL;
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

int distribute_object(struct object *parent)
{
  struct object *o;   /* Iterator for child objects */
  int i;
  int vert_index; /* index of the vertex in the global array 
                     "world->all_vertices" */
 
  if (parent->object_type == BOX_OBJ || parent->object_type == POLY_OBJ)
  {
    for (i=0;i<parent->n_walls;i++)
    {
      if (parent->wall_p[i]==NULL) continue;  /* Wall removed. */
      
      parent->wall_p[i] = distribute_wall(parent->wall_p[i]);

      if (parent->wall_p[i]==NULL)
        mcell_allocfailed("Failed to distribute wall %d on object %s.", i, parent->sym->name);

      /* create information about shared vertices */
      if(world->create_shared_walls_info_flag)
      {
         vert_index = parent->wall_p[i]->vert[0] - world->all_vertices;
         push_wall_to_list(&(world->walls_using_vertex[vert_index]), parent->wall_p[i]);
         vert_index = parent->wall_p[i]->vert[1] - world->all_vertices;
         push_wall_to_list(&(world->walls_using_vertex[vert_index]), parent->wall_p[i]);
         vert_index = parent->wall_p[i]->vert[2] - world->all_vertices;
         push_wall_to_list(&(world->walls_using_vertex[vert_index]), parent->wall_p[i]);
      
      }
    }
    if (parent->walls!=NULL)
    {
      free(parent->walls);
      parent->walls = NULL;  /* Use wall_p from now on! */
    }
  }
  else if (parent->object_type == META_OBJ)
  {
    for (o = parent->first_child; o != NULL; o = o->next)
    {
      if (distribute_object(o) != 0) return 1;
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

int distribute_world(void)
{
  struct object *o;     /* Iterator for objects in the world */

  for (o = world->root_instance ; o != NULL ; o = o->next)
  {
    if (distribute_object(o) != 0) return 1;
  }
  
  return 0;
}

/***************************************************************************
closest_pt_point_triangle:
  In:  p - point 
       a,b,c - vectors defining the vertices of the triangle. 
  Out: final_result - closest point on triangle ABC to a point p. 
       The code is adapted from "Real-time Collision Detection" by Christer Ericson, ISBN 1-55860-732-3, p.141.
       
***************************************************************************/
void closest_pt_point_triangle(struct vector3 *p, struct vector3 *a, struct vector3 *b, struct vector3 *c, struct vector3 *final_result)
{
   struct vector3 ab, ac, ap, bp, cp, result1;
   double d1, d2, d3, d4, vc, d5, d6, vb, va, denom, v, w;

   /* Check if P in vertex region outside A */
   vectorize(a, b, &ab);
   vectorize(a, c, &ac);
   vectorize(a, p, &ap);
   d1 = dot_prod(&ab, &ap);
   d2 = dot_prod(&ac, &ap);
   if(d1 <= 0.0f && d2 <= 0.0f) {
       memcpy(final_result,a,sizeof(struct vector3)); /* barycentric coordinates (1,0,0) */
       return;
   }

   /* Check if P in vertex region outside B */
   vectorize(b, p, &bp);
   d3 = dot_prod(&ab, &bp);
   d4 = dot_prod(&ac, &bp);
   if(d3 >= 0.0f && d4 <= d3) {
      memcpy(final_result,b,sizeof(struct vector3)); /* barycentric coordinates (0,1,0) */
      return;
   }

   /* Check if P in edge region of AB, if so return projection of P onto AB */
   vc = d1*d4 - d3*d2;
   if(vc <= 0.0f && d1 >= 0.0f && d3 <= 0.0f) {
        v = d1 / (d1 - d3);
        scalar_prod(&ab, v, &result1);
        vect_sum(a, &result1, final_result);
        return;  /* barycentric coordinates (1-v,v,0) */
   }

   /* Check if P in vertex region outside C */
   vectorize(c, p, &cp);
   d5 = dot_prod(&ab, &cp);
   d6 = dot_prod(&ac, &cp);
   if(d6 >=0.0f && d5 <= d6) {
        memcpy(final_result,c,sizeof(struct vector3));  /* barycentric coordinates (0,0,1) */
        return;
   }

   /* Check if P in edge region of AC, if so return projection of P onto AC */
   vb = d5*d2 - d1*d6;
   if(vb <= 0.0f && d2 >= 0.0f && d6 <= 0.0f){
      w = d2/ (d2 - d6);
      scalar_prod(&ac, w, &result1);
      vect_sum(a, &result1, final_result);
      return;      /* barycentric coordinates (0, 1-w,w) */
   }
  
   /* Check if P in edge region of BC, if so return projection of P onto BC */
   va = d3*d6 - d5*d4;
   if(va <= 0.0f && (d4 - d3) >= 0.0f && (d5 - d6) >= 0.0f) {
	w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        vectorize(b, c, &result1);
        scalar_prod(&result1, w, &result1);
        vect_sum(b, &result1, final_result);
        return;  /*barycentric coordinates (0,1-w, w) */
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
   return;   /* = u*a + v*b + w*c, u = va * denom = 1.0f - v -w */

}
/***************************************************************************
test_sphere_triangle:
  In:  s - center of the sphere
       radius - radius of the sphere
       a,b,c - vectors to the vertices of the triangle.  
  Out: Returns 1 if sphere intersects triangle ABC, 0 - otherwise.
       The point p on ABC closest to the sphere center is also returned.
       The code is adapted from "Real-time Collision Detection" by Christer Ericson, ISBN 1-55860-732-3, p.167.
       
***************************************************************************/
int test_sphere_triangle(struct vector3 *s, double radius, struct vector3 *a, struct vector3 *b, struct vector3 *c, struct vector3 *p)
{
   struct vector3 v;
   v.x = 0;
   v.y = 0;
   v.z = 0;

   /* Find point P on triangle ABC closest to the sphere center. */
   closest_pt_point_triangle(s,a,b,c,p);

   /* Sphere and triangle intersect if the (squared) distance from the sphere
      center to point p is less than the (squared) sphere radius. */
      
    vectorize(s, p, &v);
    return (dot_prod(&v,&v) <= radius*radius);

}

/***************************************************************************
compute_plane:
  In:  a,b,c - vectors to the three noncollinear points.
       p - pointer to the struct plane.  
  Out: Computes plane equation. 
       Returnes plane.
       The code is adapted from "Real-time Collision Detection" by Christer Ericson, ISBN 1-55860-732-3, p.55.
       
***************************************************************************/

void compute_plane(struct vector3 *a, struct vector3 *b, struct vector3 *c, struct plane *p)
{
	struct vector3 ba, ca;

        vectorize(a, b, &ba);
        vectorize(a, c, &ca);

        cross_prod(&ba, &ca, &(p->n));
        /* normalize the plane normal */
        normalize(&(p->n));

        p->d = dot_prod(&(p->n), a);

        return;
}

/***************************************************************************
test_sphere_plane:
  In:  s - center of the sphere
       radius - radius of the sphere
       p - struct plane.  
  Out: Returns 1 if sphere intersects plane p, 0 - otherwise.
       The code is adapted from "Real-time Collision Detection" by Christer Ericson, ISBN 1-55860-732-3, p.160.
       
***************************************************************************/
int test_sphere_plane(struct vector3 *s, double radius, struct plane *p)
{
	/* For a normalized plane (|p.n = 1|), evaluating the plane equation
           for a point gives the signed distance of the point to the plane */

        double dist;
    
        dist = dot_prod(s, &(p->n)) - p->d;
        /* If sphere center within +/- radius from the plane, plane 
           intersects sphere */
        return fabs(dist) <= radius;

}
/***************************************************************************
test_sphere_ray:
  In:  p - start point of the ray
       d - unit vector of the ray
       s - center of the sphere
       radius - radius of the sphere
       t - parameter in the ray equation (r = p + t*d)
       q - point of the intersection 
  Out: Returns 1 if sphere intersects ray, 0 - otherwise.
       The code is adapted from "Real-time Collision Detection" by Christer Ericson, ISBN 1-55860-732-3, p.178.
       
***************************************************************************/
int test_sphere_ray(struct vector3 *p, struct vector3 *d, struct vector3 *s,
       double radius, double *t, struct vector3 *q)
{
	struct vector3 m, result;
        double b,c, discr;

        vectorize(s, p, &m);
        b = dot_prod(&m, d);
        c = dot_prod(&m,&m) - radius*radius;
        /* exit if ray's origin outside sphere (c > 0) and ray pointing 
           away from sphere ( b > 0) */   
        if (c > 0.0 && b > 0.0) return 0;
        
        discr = b*b - c;
        /* A negative discriminant corresponds to the ray missing sphere */
        if (discr < 0.0) return 0;
        /* ray now found to intersect sphere, compute smallest t value of 
           intersection */         
        *t = - b - sqrt(discr);
        /* If t is negative, ray started inside sphere so clamp t to zero */
        if (*t < 0.0) *t = 0.0;
        scalar_prod(d, *t, &result);
        vect_sum(p, &result, q);
        return 1;

}

/***************************************************************************
test_segment_plane:
  In:  a - start point of the segment
       b - end point of the segment
       p - plane
       t - parameter in the ray equation (r = p + t*d)
       q - point of the intersection 
  Out: Returns 1 if segment intersects plane, 0 - otherwise.
       The code is adapted from "Real-time Collision Detection" by Christer Ericson, ISBN 1-55860-732-3, p.176.
       
***************************************************************************/
int test_segment_plane(struct vector3 *a, struct vector3 *b, struct plane *p, double *t, struct vector3 *q)
{
   struct vector3 ab, t_ab; 
   double n_a, n_ab;

    /* Compute the t value for the directed line ab intersecting the plane */
    vectorize(a, b, &ab);
    n_a = dot_prod(&(p->n), a);
    n_ab = dot_prod(&(p->n), &ab);

    if(n_ab == 0) return 0; /* segment is parallel to the plane */

    *t = (p->d - n_a) / n_ab;

    /* If t in [0..1] compute and return intersection point */

   if((*t >= 0.0) && (*t <= 1.0)) {
	scalar_prod(&ab, *t, &t_ab);
        vect_sum(a, &t_ab, q);
        return 1;
   }
   /* Else no intersection */
   return 0;
  
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
int test_bounding_boxes(struct vector3 *llf1, struct vector3 *urb1, struct vector3 *llf2, struct vector3 *urb2)
{
  /* Two boxes overlap only if they overlap on all three axes
     while their extent along each dimension is seen as an interval
     on the corresponding axis. */

  /* exit with no intersection is separated along axis */
  if((urb1->x <  llf2->x) || (llf1->x > urb2->x)) return 0;
  if((urb1->y <  llf2->y) || (llf1->y > urb2->y)) return 0;
  if((urb1->z <  llf2->z) || (llf1->z > urb2->z)) return 0;
  /* Overlapping on all axis means that boxes are intersecting. */
  return 1;
}



/***************************************************************************
surface_point_in_region:
  In: object upon which the point is
      index of the wall on that object
      vector of the actual point
      expression for where we should release
  Out: Returns 1 if the surface point is in the release specification,
       0 otherwise.
***************************************************************************/
int surface_point_in_region(struct object *ob,int wall_n,struct vector3 *v,struct release_evaluator *expr)
{
  struct subvolume *sv = find_subvolume(v,NULL);
  struct waypoint *wp = &(world->waypoints[sv - world->subvol]);
  struct wall_list *wl;
  struct wall_list pre_wall,my_wall;
  struct vector3 delta,hit;
  struct region_list *irl,*trl,*ttrl,*rl,*arl,**qrl,**prl;
  double t;
  int i;
  
  rl = arl = NULL;
  delta.x = v->x - wp->loc.x;
  delta.y = v->y - wp->loc.y;
  delta.z = v->z - wp->loc.z;
  pre_wall.next = &my_wall;
  my_wall.next = NULL;
  my_wall.this_wall = ob->wall_p[wall_n];
  
  for (wl=sv->wall_head ; wl!=NULL ; wl=wl->next)
  {
    if (wl!=&my_wall)
    {
      if (wl->this_wall==my_wall.this_wall) continue;  /* Don't try to collide with the wall the point is on */
      i = collide_wall(&(wp->loc) , &delta , wl->this_wall , &t , &hit , 0);
      if (i==COLLIDE_MISS || i==COLLIDE_REDO || !(t >= 0 && t < 1.0)) continue;
    }
    else i = COLLIDE_FRONT;
    
    for (irl = wl->this_wall->parent_object->regions ; irl!=NULL ; irl=irl->next)
    {
      if (!get_bit(irl->reg->membership,wl->this_wall->side)) continue;
      if (i==COLLIDE_FRONT) { prl=&rl; qrl=&arl; }
      else { qrl=&rl; prl=&arl; }

      if (*prl==irl) { trl=(*prl)->next; mem_put(sv->local_storage->regl,*prl); (*prl)=trl; }
      else if (*prl!=NULL)
      {
	for (trl=*prl; trl->next!=NULL && trl->next!=irl ; trl=trl->next) {}
	if (trl->next!=NULL) { ttrl = trl->next; trl->next=ttrl->next; mem_put(sv->local_storage->regl,ttrl); }
	else { trl = (struct region_list*) CHECKED_MEM_GET(sv->local_storage->regl, "region list"); trl->reg=irl->reg; trl->next=*qrl; *qrl=trl; }
      }
    }
    if (wl->next==NULL && wl!=&my_wall) wl=&pre_wall; /* Cheat to go through loop one extra time with the wall on which the point is */
  }
  
  i = eval_rel_region_3d(expr,wp,rl,arl);
  
  if (rl!=NULL) mem_put_list(sv->local_storage->regl,rl);
  if (arl!=NULL) mem_put_list(sv->local_storage->regl,arl);
  
  return i;
}


/* Helper struct for release_onto_regions and vacuum_from_regions */
struct reg_rel_helper_data
{
  struct reg_rel_helper_data *next;
  struct surface_grid *grid;
  int index;
  double my_area;
};


/***************************************************************************
vacuum_from_regions:
  In: a release site object
      a template grid molecule we're going to remove
      the number of molecules to remove
  Out: 0 on success, 1 on failure.  Molecules of the specified type are
       removed uniformly at random from the free area in the regions
       specified by the release site object.
  Note: if the user requests to remove more molecules than actually exist,
        the function will return success and not give a warning.  The only
	reason to return failure is an out of memory condition.
***************************************************************************/
static int vacuum_from_regions(struct release_site_obj *rso,struct grid_molecule *g,int n)
{
  struct release_region_data *rrd;
  struct mem_helper *mh;
  struct reg_rel_helper_data *rrhd_head,*p;
  int n_rrhd;
  struct wall *w;
  struct grid_molecule *gp;  
  
  rrd = rso->region_data;

  mh = create_mem( sizeof(struct reg_rel_helper_data) , 1024 );
  if (mh==NULL) return 1;
  
  rrhd_head = NULL;
  n_rrhd=0;
  
  for (int n_object=0; n_object<rrd->n_objects; n_object++)
  {
    if (rrd->walls_per_obj[n_object]==0) continue;
    for (int n_wall=0; n_wall<rrd->in_release[n_object]->nbits; n_wall++)
    {
      if (!get_bit(rrd->in_release[n_object], n_wall)) continue;
      
      w = rrd->owners[n_object]->wall_p[n_wall];
      
      if (w->grid==NULL) continue;
      
      for (unsigned int n_tile=0; n_tile<w->grid->n_tiles; n_tile++)
      {
        gp = w->grid->mol[n_tile];
        if (gp!=NULL)
        {
          if (gp->properties == g->properties)
          {
            if (rrd->refinement && !grid_release_check(rrd, n_object, n_wall, n_tile, NULL)) continue;
            p = CHECKED_MEM_GET_NODIE(mh, "release region helper data");
            if (p==NULL) return 1;
            
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

  for (p=rrhd_head ; n<0 && n_rrhd>0 && p!=NULL ; p=p->next , n_rrhd--)
  {
    if (rng_dbl(world->rng) < ((double)(-n))/((double)n_rrhd))
    {
      gp = p->grid->mol[ p->index ];
      gp->properties->population--;
      if ((gp->properties->flags & (COUNT_CONTENTS|COUNT_ENCLOSED)) != 0)
        count_region_from_scratch((struct abstract_molecule*)gp,NULL,-1,NULL,gp->grid->surface,gp->t);
      gp->properties = NULL;
      p->grid->mol[ p->index ] = NULL;
      p->grid->n_occupied--;
      if (gp->flags & IN_SCHEDULE)
      {
        gp->grid->subvol->local_storage->timer->defunct_count++; /* Tally for garbage collection */
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
      a template grid molecule we're going to release
      the number of molecules to release
  Out: 0 on success, 1 on failure.  Molecules are released uniformly at
       random onto the free area in the regions specified by the release
      site object.
  Note: if the CCNNUM method is used, the number passed in is ignored.
***************************************************************************/
int release_onto_regions(struct release_site_obj *rso,struct grid_molecule *g,int n)
{
  int success,failure;
  double est_sites_avail;
  double seek_cost,pick_cost;
  const double rel_list_gen_cost = 10.0;  /* Just a guess */
  const int too_many_failures = 10;       /* Also a guess */
  struct release_region_data *rrd;
  struct mem_helper *mh;
  int i;
  unsigned int grid_index;
  double A,max_A, num_to_release;
  struct wall *w;
  struct grid_molecule *new_g;
  struct subvolume *gsv = NULL;
  struct vector3 pos3d;

  long long skipped_placements = 0;
  int is_complex = 0;
  if (g->properties->flags & IS_COMPLEX)
    is_complex = 1;
  
  rrd = rso->region_data;
  
  success = failure = 0;
  seek_cost = 0;
  
  max_A = rrd->cum_area_list[rrd->n_walls_included-1];
  est_sites_avail = (int)max_A;
  pick_cost = rel_list_gen_cost * est_sites_avail;
  
  if (rso->release_number_method == DENSITYNUM)
  {
    num_to_release = rso->concentration * est_sites_avail / world->grid_density;
    if (num_to_release > (double) INT_MAX)
      mcell_error("Release site \"%s\" tries to release more than INT_MAX (2147483647) molecules.", rso->name);
    n = (int)(num_to_release);

  }
  
  if (n<0) return vacuum_from_regions(rso,g,n);
  if(world->notify->release_events==NOTIFY_FULL) 
  {
    if(n > 0)
      mcell_log_raw("Releasing %d molecules %s ...", n, g->properties->sym->name);
  }
  
  while (n>0)
  {
    if (! is_complex  &&  failure >= success+too_many_failures)
    {
      seek_cost = n*( ((double)(success+failure+2))/((double)(success+1)) );
    }
    if (seek_cost < pick_cost)
    {
      A = rng_dbl( world->rng )*max_A;
      i = bisect_high( rrd->cum_area_list , rrd->n_walls_included , A );
      w = rrd->owners[rrd->obj_index[i]]->wall_p[ rrd->wall_index[i] ];
      
      if (w->grid==NULL)
      {
        if (create_grid(w, NULL))
          return 1;
      }
      if (i) A -= rrd->cum_area_list[i-1];
      grid_index = (int)((w->grid->n*w->grid->n)*(A/w->area));
      if (grid_index>=w->grid->n_tiles) grid_index=w->grid->n_tiles-1;

      if (is_complex)
      {
        short orient = 0;
        if (rso->orientation > 0) orient = 1;
        else if (rso->orientation < 0) orient = -1;
        else
        {
          orient = (rng_uint(world->rng)&1)?1:-1;
        }
        struct grid_molecule *gp = macro_insert_molecule_grid_2(g->properties, orient, w, grid_index, g->t, NULL, rrd);
        if (gp == NULL)
        {
          ++ failure;
          if (failure == world->complex_placement_attempts)
          {
            -- n;
            if (++ skipped_placements >= world->notify->complex_placement_failure_threshold)
            {
              switch (world->notify->complex_placement_failure)
              {
                case WARN_COPE:
                  break;

                case WARN_WARN:
                  mcell_warn("Could not release %lld of %s (surface full).",
                             skipped_placements + n,
                             g->properties->sym->name);
                  break;

                case WARN_ERROR:
                  mcell_error("Could not release %lld of %s (surface full).",
                              skipped_placements + n,
                              g->properties->sym->name);
                  return 1;

                default: UNHANDLED_CASE(world->notify->complex_placement_failure);
              }
              break;
            }
          }
        }
        else
        {
          failure = 0;
          ++ success;
          -- n;
        }
      }
      else
      {
        if (w->grid->mol[grid_index] != NULL || (rrd->refinement && !grid_release_check(rrd,rrd->obj_index[i],rrd->wall_index[i],grid_index,NULL))) failure++;
        else
        {
          struct vector2 s_pos;
          if (world->randomize_gmol_pos) grid2uv_random(w->grid,grid_index,&s_pos);
          else grid2uv(w->grid,grid_index,&s_pos);
          uv2xyz(&s_pos, w, &pos3d);
          gsv = find_subvolume(&pos3d, gsv);

          new_g = (struct grid_molecule*)CHECKED_MEM_GET( gsv->local_storage->gmol, "grid molecule" );
          if (new_g==NULL) return 1;
          memcpy(new_g,g,sizeof(struct grid_molecule));
          new_g->birthplace = w->grid->subvol->local_storage->gmol;
          new_g->grid_index = grid_index;
          new_g->s_pos.u = s_pos.u;
          new_g->s_pos.v = s_pos.v;

          if (rso->orientation > 0) new_g->orient = 1;
          else if (rso->orientation < 0) new_g->orient = -1;
          else {
            new_g->orient = (rng_uint(world->rng)&1)?1:-1;
          }

          new_g->grid = w->grid;

          w->grid->mol[grid_index] = new_g;

          w->grid->n_occupied++;
          new_g->properties->population++;
          if ((new_g->properties->flags&COUNT_ENCLOSED) != 0) new_g->flags |= COUNT_ME;
          if (new_g->properties->flags & (COUNT_CONTENTS|COUNT_ENCLOSED))
            count_region_from_scratch((struct abstract_molecule*)new_g,NULL,1,NULL,new_g->grid->surface,new_g->t);

          if (schedule_add( gsv->local_storage->timer, new_g))
            return 1;

          success++;
          n--;
        }
      }
    }
    else
    {
      mh = create_mem( sizeof(struct reg_rel_helper_data) , 1024 );
      if (mh==NULL) return 1;

      struct reg_rel_helper_data *rrhd_head = NULL;
      int n_rrhd=0;
      max_A=0;
      for (int n_object=0; n_object<rrd->n_objects; n_object++)
      {
        if (rrd->walls_per_obj[n_object]==0) continue;
        for (int n_wall=0; n_wall<rrd->in_release[n_object]->nbits; n_wall++)
        {
          if (!get_bit(rrd->in_release[n_object], n_wall)) continue;
          
          w = rrd->owners[n_object]->wall_p[n_wall];
          
          if (w->grid==NULL)
          {
            if (create_grid(w,NULL))
              return 1;
          }
          else if (w->grid->n_occupied == w->grid->n_tiles) continue;
          
          A = w->area / (w->grid->n_tiles);
          
          for (unsigned int n_tile=0; n_tile<w->grid->n_tiles; n_tile++)
          {
            if (w->grid->mol[n_tile]==NULL && !(rrd->refinement && !grid_release_check(rrd, n_object, n_wall, n_tile, NULL)))
            {
              struct reg_rel_helper_data *new_rrd = CHECKED_MEM_GET_NODIE(mh, "release region helper data");
              if (new_rrd == NULL) return 1;

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
      
      
      for (struct reg_rel_helper_data *this_rrd = rrhd_head; this_rrd != NULL && n>0; this_rrd = this_rrd->next)
      {
        if (n>=n_rrhd || rng_dbl(world->rng)<(this_rrd->my_area/max_A)*((double)n))
        {
          struct vector2 s_pos;
          if (world->randomize_gmol_pos) grid2uv_random(this_rrd->grid,this_rrd->index,&s_pos);
          else grid2uv(this_rrd->grid,this_rrd->index,&s_pos);
          uv2xyz(&s_pos, this_rrd->grid->surface, &pos3d);
          gsv = find_subvolume(&pos3d, gsv);

          new_g = (struct grid_molecule*)CHECKED_MEM_GET( gsv->local_storage->gmol, "grid molecule" );
          if (new_g==NULL) return 1;
          memcpy(new_g,g,sizeof(struct grid_molecule));
          new_g->birthplace = this_rrd->grid->subvol->local_storage->gmol;
          new_g->grid_index = this_rrd->index;
          new_g->s_pos.u = s_pos.u;
          new_g->s_pos.v = s_pos.v;
	  
	  if (rso->orientation>0) new_g->orient=1;
	  else if (rso->orientation<0) new_g->orient=-1;
	  else{ 
            new_g->orient = (rng_uint(world->rng)&1)?1:-1;
	  }
          new_g->grid = this_rrd->grid;
          
          this_rrd->grid->mol[ this_rrd->index ] = new_g;

          this_rrd->grid->n_occupied++;
          new_g->properties->population++;
          if ((new_g->properties->flags&COUNT_ENCLOSED) != 0) new_g->flags |= COUNT_ME;
          if (new_g->properties->flags & (COUNT_CONTENTS|COUNT_ENCLOSED))
            count_region_from_scratch((struct abstract_molecule*)new_g,NULL,1,NULL,NULL,new_g->t);

          if (schedule_add(gsv->local_storage->timer, new_g))
            return 1;
          
          n--;
          n_rrhd--;
        }
        max_A -= this_rrd->my_area;
      }
      
      delete_mem(mh);
      
      if (n>0)
      {
        switch (world->notify->mol_placement_failure)
        {
          case WARN_COPE:
            break;

          case WARN_WARN:
            mcell_warn("Could not release %d of %s (surface full).", n, g->properties->sym->name);
            break;

          case WARN_ERROR:
            mcell_error("Could not release %d of %s (surface full).", n, g->properties->sym->name);
            return 1;

          default: UNHANDLED_CASE(world->notify->mol_placement_failure);
        }
        break;
      }
    }
  }
  
  return 0;
}


/****************************************************************************
push_wall_to_list:
   In: head of the linked list of walls
       wall to be added to the list
   Out: none. The linked list is updated.
   Note: No check for duplicates is performed
*****************************************************************************/
void push_wall_to_list(struct wall_list **wall_nbr_head, struct wall *w)
{
   struct wall_list *old_head, *wlp;

   old_head = *wall_nbr_head;

   wlp = CHECKED_MALLOC_STRUCT(struct wall_list, "wall_list");
   wlp->this_wall = w;

   if(old_head == NULL)
   {
      wlp->next = NULL;
      old_head = wlp;
   }else{
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
void delete_wall_list(struct wall_list *wl_head)
{
   struct wall_list *nnext;
   while(wl_head != NULL)
   {
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
struct wall_list* find_nbr_walls_shared_one_vertex(struct wall *origin, int *shared_vert)
{
  int i;
  struct wall_list *wl;
  struct wall_list *head = NULL;

  if(!world->create_shared_walls_info_flag) mcell_internal_error("Function 'find_nbr_walls_shared_one_vertex()' is called but shared walls information is not created.");

  for(i = 0; i < 3; i++)
  {
     if(shared_vert[i] >= 0)
     {
        for(wl = world->walls_using_vertex[shared_vert[i]]; wl != NULL; wl = wl->next)
        {
           if(wl->this_wall == origin) continue;

           if(!walls_share_full_edge(origin, wl->this_wall))
           { 
              push_wall_to_list(&head, wl->this_wall);
           } 
        }

     }
  }

  return head;

}


/************************************************************************
wall_share_vertex:
   In:  wall
        vertex of another wall
   Out: 1 if wall share the argument vertex, 0 - otherwise
************************************************************************/
int wall_share_vertex(struct wall *w, struct vector3 *vert)
{
   int i;

   for(i = 0; i < 3; i++)
   {
     if(!distinguishable_vec3(w->vert[i], vert, EPS_C)) return 1;
   }
   return 0;

}

/***********************************************************************
walls_share_full_edge:
  In: two walls
  Out: 1 if the walls share a full edge, 0 - otherwise.
       Here by "full" we mean that the shared edge has two endpoints
       that are the vertices of both walls w1 and w2.
************************************************************************/
int walls_share_full_edge(struct wall *w1, struct wall *w2)
{
   int i, k;
   int  count = 0; /* count number of shared vertices between two walls */
    
   for(i = 0; i < 3; i++)
   {
     for(k = 0; k < 3; k++)
     {
        if(!distinguishable_vec3(w1->vert[i], w2->vert[k], EPS_C)) count++;
     }
   }

   if(count == 2) return 1;

   return 0;
}

/***********************************************************************
find_region_by_wall:
  In: object
      wall
  Out: an object's region list if the wall belongs to one, NULL - otherwise.
  Note: regions called "ALL" or the ones that have ALL_ELEMENTS are not 
        included in the return "region list".  This is done intentionally
        since the function is used to determine region border and the 
        above regions do not have region borders.
************************************************************************/
struct region_list * find_region_by_wall(struct object *parent, struct wall *this_wall)
{
  struct region *rp;
  struct region_list *rlp, *rlps, *rlp_head = NULL;
  int this_wall_idx = -1;

  for(int i = 0; i < parent->n_walls; i++)
  {
    if(parent->wall_p[i] == this_wall)
    {
       this_wall_idx = i;
       break;
    }
  }

  for(rlp = parent->regions; rlp != NULL; rlp = rlp->next)
  {
    rp = rlp->reg;
    if((strcmp(rp->region_last_name,"ALL") == 0) || (rp->region_has_all_elements))  continue;

    if(rp->membership == NULL)
       mcell_internal_error("Missing region membership for '%s'.", rp->sym->name);

    if(get_bit(rp->membership, this_wall_idx))
    {
      rlps = CHECKED_MALLOC_STRUCT(struct region_list, "region_list");
      rlps->reg = rp;

      if(rlp_head == NULL)
      {
        rlps->next = NULL;
        rlp_head = rlps;
      }else{
        rlps->next = rlp_head;
        rlp_head = rlps;

      }
    }
  }

  return rlp_head;

}

/***********************************************************************
is_wall_edge_region_border:
  In: wall
      wall's edge
  Out: 1 if the edge is a region's border, and 0 - otherwise.
  Note: we do not specify any particular region here, any region will 
        suffice
************************************************************************/
int is_wall_edge_region_border(struct wall *this_wall, struct edge *this_edge)
{
  struct region_list *rlp, *rlp_head;
  struct region *rp;
  void *key;
  unsigned int keyhash;
 
  int is_region_border = 0;  /* flag */

  rlp_head = find_region_by_wall(this_wall->parent_object, this_wall);

  /* If this wall is not a part of any region (note that we do not consider 
     region called ALL here) */
  if(rlp_head == NULL) return is_region_border;
  
  for(rlp = rlp_head; rlp != NULL; rlp = rlp->next)
  {
    rp = rlp->reg;

    if(rp->boundaries == NULL) mcell_internal_error("Region '%s' of the object '%s' has no boundaries.", rp->region_last_name, this_wall->parent_object->sym->name);

    keyhash = (unsigned int)(intptr_t)(this_edge);
    key = (void *)(this_edge);
  
    if(pointer_hash_lookup(rp->boundaries, key, keyhash))
    {
      is_region_border = 1;
      break;
    }
  }

  if(rlp_head != NULL) delete_void_list((struct void_list *)rlp_head);

  return is_region_border;

}

/*************************************************************************
find_shared_edge_index_of_neighbor_wall:
  In: original wall
      neighbor wall
  Out: index of the shared edge in the coordinate system of neighbor wall.

**************************************************************************/
int find_shared_edge_index_of_neighbor_wall(struct wall *orig_wall, struct wall *nbr_wall)
{
   int nbr_edge_ind = -1;
   int shared_vert_ind_1 = -1, shared_vert_ind_2 = -1;

   find_shared_vertices_for_neighbor_walls(orig_wall, nbr_wall, &shared_vert_ind_1, &shared_vert_ind_2);
   
   if((shared_vert_ind_1 + shared_vert_ind_2) == 1)
   {
      nbr_edge_ind = 0;
   }else if((shared_vert_ind_1 + shared_vert_ind_2) == 2){
      nbr_edge_ind = 2;
   }else if((shared_vert_ind_1 + shared_vert_ind_2) == 3){
      nbr_edge_ind = 1;
   }else{
      mcell_internal_error("Error in the function 'find_shared_edge_index_of_neighbor_wall()");
   }

   return nbr_edge_ind;

}

/****************************************************************************
find_neighbor_wall_and_edge:
  In: wall
      wall edge index ( in the coordinate system of "wall")
      neighbor wall (return value)
      index of the edge in the coordinate system of 
      "neighbor wall" that is shared with "wall" and 
      coincides with the edge with "wall edge index" (return value)

****************************************************************************/
void find_neighbor_wall_and_edge(struct wall *orig_wall, int orig_edge_ind, struct wall **nbr_wall, int *nbr_edge_ind)
{
  int ii;
  struct wall *w;
  struct vector3 *vert_A, *vert_B;

  switch(orig_edge_ind)
  {
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
       break;
  }

  for(ii = 0; ii < 3; ii++)
  {
    w = orig_wall->nb_walls[ii];
    if(w == NULL) continue;  
  
   
    if(wall_contains_both_vertices(w, vert_A, vert_B))
    {
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
int wall_contains_both_vertices(struct wall *w, struct vector3 *vert_A, struct vector3 *vert_B)
{
  int count = 0, ii;
  struct vector3 *v;
  
  for(ii = 0; ii < 3; ii++)
  {
    v = w->vert[ii];
  
    if((!distinguishable_vec3(v, vert_A, EPS_C)) || (!(distinguishable_vec3(v, vert_B, EPS_C))))
    {
      count++;
    }

  }
  
  if(count == 2) return 1;
  else return 0;

}

/*******************************************************************
are_walls_coincident:
  In: first wall
      second wall
      accuracy of the comparison
  Out: 0 if the walls are not coincident
       1 if the walls are coincident
*******************************************************************/
int are_walls_coincident(struct wall *w1, struct wall *w2, double eps)
{
  if((w1 == NULL) || (w2 == NULL)) return 0;

  int count = 0;

  if(!distinguishable_vec3(w1->vert[0], w2->vert[0], eps)) count ++; 
  if(!distinguishable_vec3(w1->vert[0], w2->vert[1], eps)) count ++; 
  if(!distinguishable_vec3(w1->vert[0], w2->vert[2], eps)) count ++; 
  if(!distinguishable_vec3(w1->vert[1], w2->vert[0], eps)) count ++; 
  if(!distinguishable_vec3(w1->vert[1], w2->vert[1], eps)) count ++; 
  if(!distinguishable_vec3(w1->vert[1], w2->vert[2], eps)) count ++; 
  if(!distinguishable_vec3(w1->vert[2], w2->vert[0], eps)) count ++; 
  if(!distinguishable_vec3(w1->vert[2], w2->vert[1], eps)) count ++; 
  if(!distinguishable_vec3(w1->vert[2], w2->vert[2], eps)) count ++; 

  if(count >= 3) return 1;

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
int are_walls_coplanar(struct wall *w1, struct wall *w2, double eps)
{

  /* find the plane equation of the second wall in the form (n*x + d2 = 0) */
  
  double d2, d1_0, d1_1, d1_2;
 
  d2 = -dot_prod(&(w2->normal), w2->vert[0]);

  /* check whether all vertices of the first wall satisfy
     plane equation of the second wall */  
  d1_0 = dot_prod(&(w2->normal), w1->vert[0]) + d2;
  d1_1 = dot_prod(&(w2->normal), w1->vert[1]) + d2;
  d1_2 = dot_prod(&(w2->normal), w1->vert[2]) + d2;

  if((!distinguishable(d1_0, 0, eps)) && (!distinguishable(d1_1, 0, eps)) 
     && (!distinguishable(d1_2, 0, eps))) 
  {
     return 1;
  }   

  return 0;

}


/**********************************************************************
* overlap_coplanar_walls:
* In: first wall
*     second wall
* Out: 1 if the walls overlap
*      0 if the walls do not overlap
* Note: the walls are assumed to be coplanar and no special check is
*       performed here for coplanarity.
***********************************************************************/
int overlap_coplanar_walls(struct wall *w1, struct wall *w2)
{
  /* check whether each of the vertices of w1 lie inside w2
     and vice versa */
  if(point_inside_triangle(w1->vert[0], w2->vert[0], w2->vert[1], w2->vert[2]))
  {
     return 1;
  }
  if(point_inside_triangle(w1->vert[1], w2->vert[0], w2->vert[1], w2->vert[2])) 
  {
      return 1;
  }
  if(point_inside_triangle(w1->vert[2], w2->vert[0], w2->vert[1], w2->vert[2])) 
  {
     return 1;
  }


  if(point_inside_triangle(w2->vert[0], w1->vert[0], w1->vert[1], w1->vert[2]))
  {
     return 1;
  }
  if(point_inside_triangle(w2->vert[1], w1->vert[0], w1->vert[1], w1->vert[2]))
  {
     return 1;
  }
  if(point_inside_triangle(w2->vert[2], w1->vert[0], w1->vert[1], w1->vert[2]))
  {
     return 1;
  }

  return 0;
}

/***********************************************************************
* overlap_tri_tri_3d:
*  In: arrays of doubles representing coordinates of vertices
*      and normals of the triangles
*  Out: 1 if triangles overlap
*       0 if triangles do not overlap
  Note: see "Real-time rendering" 2nd Ed., by Tomas Akenine-Moller and 
        Eric Haines, pp. 582, 592.  Also based on "Fast and Robust 
        Triangle-Triangle Overlap Test Using Orientation Predicates"
        by P. Guigue and O. Devillers, Journal of Graphic Tools,
        8(1), 2003.
  Note: walls are assumed to be coplanar. No separate check is done
        for coplanarity inside this function.
***********************************************************************/
int overlap_tri_tri_3d(double p1[3], double q1[3], double r1[3],
		       double p2[3], double q2[3], double r2[3],
		       double normal_1[3], double normal_2[3])
{
  /* Since triangles are coplanar they are projected onto 
     the axis-aligned plane where the areas of the triangles
     are maximized. Then a simple two-dimensional triangle-triangle
     overlap test is performed. Let the normal to the wall n
     has coordinates (n_x, n_y, n_z). The coordinate 
     component that corresponds to max(n_x, n_y, n_z) can be skipped
     and the others are kept as two-dimensional coordinates. */
  
  double P1[2],Q1[2],R1[2];
  double P2[2],Q2[2],R2[2];

  double n_x, n_y, n_z;

  n_x = ((normal_1[0]<0)?-normal_1[0]:normal_1[0]);
  n_y = ((normal_1[1]<0)?-normal_1[1]:normal_1[1]);
  n_z = ((normal_1[2]<0)?-normal_1[2]:normal_1[2]);


  /* Projection of the triangle in 3D onto 2D such that the area
     of the projection is maximized */

  if (( n_x > n_z ) && ( n_x >= n_y )) {
    // Project onto plane YZ

      P1[0] = q1[2]; P1[1] = q1[1];
      Q1[0] = p1[2]; Q1[1] = p1[1];
      R1[0] = r1[2]; R1[1] = r1[1]; 
    
      P2[0] = q2[2]; P2[1] = q2[1];
      Q2[0] = p2[2]; Q2[1] = p2[1];
      R2[0] = r2[2]; R2[1] = r2[1]; 
  } else if (( n_y > n_z ) && ( n_y >= n_x )) {
    // Project onto plane XZ

    P1[0] = q1[0]; P1[1] = q1[2];
    Q1[0] = p1[0]; Q1[1] = p1[2];
    R1[0] = r1[0]; R1[1] = r1[2]; 
 
    P2[0] = q2[0]; P2[1] = q2[2];
    Q2[0] = p2[0]; Q2[1] = p2[2];
    R2[0] = r2[0]; R2[1] = r2[2]; 

  } else {
    // Project onto plane XY

    P1[0] = p1[0]; P1[1] = p1[1]; 
    Q1[0] = q1[0]; Q1[1] = q1[1]; 
    R1[0] = r1[0]; R1[1] = r1[1]; 
    
    P2[0] = p2[0]; P2[1] = p2[1]; 
    Q2[0] = q2[0]; Q2[1] = q2[1]; 
    R2[0] = r2[0]; R2[1] = r2[1]; 
  }

  return tri_tri_overlap_test_2d(P1,Q1,R1,P2,Q2,R2);
    
}


 /* some 2D macros */
#define ORIENT_2D(a, b, c)  ((a[0]-c[0])*(b[1]-c[1])-(a[1]-c[1])*(b[0]-c[0]))

#define INTERSECTION_TEST_VERTEX(P1, Q1, R1, P2, Q2, R2) {\
  if (ORIENT_2D(R2,P2,Q1) >= 0.0f)\
    if (ORIENT_2D(R2,Q2,Q1) <= 0.0f)\
      if (ORIENT_2D(P1,P2,Q1) > 0.0f) {\
	if (ORIENT_2D(P1,Q2,Q1) <= 0.0f) return 1; \
	else return 0;} else {\
	if (ORIENT_2D(P1,P2,R1) >= 0.0f)\
	  if (ORIENT_2D(Q1,R1,P2) >= 0.0f) return 1; \
	  else return 0;\
	else return 0;}\
    else \
      if (ORIENT_2D(P1,Q2,Q1) <= 0.0f)\
	if (ORIENT_2D(R2,Q2,R1) <= 0.0f)\
	  if (ORIENT_2D(Q1,R1,Q2) >= 0.0f) return 1; \
	  else return 0;\
	else return 0;\
      else return 0;\
  else\
    if (ORIENT_2D(R2,P2,R1) >= 0.0f) \
      if (ORIENT_2D(Q1,R1,R2) >= 0.0f)\
	if (ORIENT_2D(P1,P2,R1) >= 0.0f) return 1;\
	else return 0;\
      else \
	if (ORIENT_2D(Q1,R1,Q2) >= 0.0f) {\
	  if (ORIENT_2D(R2,R1,Q2) >= 0.0f) return 1; \
	  else return 0; }\
	else return 0; \
    else  return 0; \
 };

#define INTERSECTION_TEST_EDGE(P1, Q1, R1, P2, Q2, R2) { \
  if (ORIENT_2D(R2,P2,Q1) >= 0.0f) {\
    if (ORIENT_2D(P1,P2,Q1) >= 0.0f) { \
        if (ORIENT_2D(P1,Q1,R2) >= 0.0f) return 1; \
        else return 0;} else { \
      if (ORIENT_2D(Q1,R1,P2) >= 0.0f){ \
	if (ORIENT_2D(R1,P1,P2) >= 0.0f) return 1; else return 0;} \
      else return 0; } \
  } else {\
    if (ORIENT_2D(R2,P2,R1) >= 0.0f) {\
      if (ORIENT_2D(P1,P2,R1) >= 0.0f) {\
	if (ORIENT_2D(P1,R1,R2) >= 0.0f) return 1;  \
	else {\
	  if (ORIENT_2D(Q1,R1,R2) >= 0.0f) return 1; else return 0;}}\
      else  return 0; }\
    else return 0; }}

int ccw_tri_tri_intersection_2d(double p1[2], double q1[2], double r1[2], 
				double p2[2], double q2[2], double r2[2]) {
  if ( ORIENT_2D(p2,q2,p1) >= 0.0f ) {
    if ( ORIENT_2D(q2,r2,p1) >= 0.0f ) {
      if ( ORIENT_2D(r2,p2,p1) >= 0.0f ) return 1;
      else INTERSECTION_TEST_EDGE(p1,q1,r1,p2,q2,r2)
    } else {  
      if ( ORIENT_2D(r2,p2,p1) >= 0.0f ) 
	INTERSECTION_TEST_EDGE(p1,q1,r1,r2,p2,q2)
      else INTERSECTION_TEST_VERTEX(p1,q1,r1,p2,q2,r2)}}
  else {
    if ( ORIENT_2D(q2,r2,p1) >= 0.0f ) {
      if ( ORIENT_2D(r2,p2,p1) >= 0.0f ) 
	INTERSECTION_TEST_EDGE(p1,q1,r1,q2,r2,p2)
      else  INTERSECTION_TEST_VERTEX(p1,q1,r1,q2,r2,p2)}
    else INTERSECTION_TEST_VERTEX(p1,q1,r1,r2,p2,q2)}
};

/**********************************************************************
* tri_tri_overlap_test_2d:
*  In: coordinates of the vertices of the two triangles
*  Out: 1 if triangles overlap
*       0 if triangles do not overlap
*  Note: triangles are assumed to be coplanar
*  Note:  Code based on "Fast and Robust Triangle-Triangle Overlap Test 
*         Using Orientation Predicates" by P. Guigue and O. Devillers, 
*         Journal of Graphic Tools, 8(1), 2003.
* http://jgt.akpeters.com/papers/GuigueDevillers03/triangle_triangle_intersectio* n.html
**********************************************************************/
int tri_tri_overlap_test_2d(double p1[2], double q1[2], double r1[2], 
			    double p2[2], double q2[2], double r2[2]) 
{
  if ( ORIENT_2D(p1,q1,r1) < 0.0f )
    if ( ORIENT_2D(p2,q2,r2) < 0.0f )
      return ccw_tri_tri_intersection_2d(p1,r1,q1,p2,r2,q2);
    else
      return ccw_tri_tri_intersection_2d(p1,r1,q1,p2,q2,r2);
  else
    if ( ORIENT_2D(p2,q2,r2) < 0.0f )
      return ccw_tri_tri_intersection_2d(p1,q1,r1,p2,r2,q2);
    else
      return ccw_tri_tri_intersection_2d(p1,q1,r1,p2,q2,r2);

}

/**********************************************************************
* sorted_insert_wall_aux_list:
* In: linked list
*     new node
* Out: new node is added to the linked list in the sorted order
***********************************************************************/
void sorted_insert_wall_aux_list(struct wall_aux_list **headRef, struct wall_aux_list *newNode)
{
  /* special case for the head end */
  if(*headRef == NULL || (*headRef)->d_prod >= newNode->d_prod)
  {
    newNode->next = *headRef;
    *headRef = newNode;
  }
  else {
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
void delete_wall_aux_list(struct wall_aux_list *head)
{
  struct wall_aux_list *nnext;
  while(head != NULL)
  {
    nnext = head->next;
    free(head);
    head = nnext;
  }
}

/*****************************************************************
walls_belong_to_same_region:
  In: two walls
  Out: 1 if both walls belong to the same region,
       0 otherwise.
  Note: Wall can be belong to several regions simultaneously.
        It is important that both walls belong to at least one
        same region.
******************************************************************/
int walls_belong_to_same_region(struct wall *w1, struct wall *w2)
{
  struct region_list *rl_1, *rl_2, *rl_t1, *rl_t2;
  struct region *rp_1, *rp_2;

  if(w1 == w2) return 1;

  rl_1 = find_region_by_wall(w1->parent_object, w1);
  if(rl_1 == NULL) return 0;

  rl_2 = find_region_by_wall(w2->parent_object, w2);
  if(rl_2 == NULL) return 0;

  for(rl_t1 = rl_1; rl_t1 != NULL; rl_t1 = rl_t1->next)
  {
    rp_1 = rl_t1->reg;
    for(rl_t2 = rl_2; rl_t2 != NULL; rl_t2 = rl_t2->next)
    {
       rp_2 = rl_t2->reg;
       if(rp_1 == rp_2) return 1;
    }
  }

  return 0;

}
