/**************************************************************************\
 ** File: wall_util.c                                                    **
 **                                                                      **
 ** Purpose: Build walls and surfaces, create edges from vertices and    **
 **    polygons.  All wall elements are assumed to be triangles.         **
 **                                                                      **
 ** Testing status: previously tested, compiles after changes.           **
\**************************************************************************/


#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "rng.h"
#include "wall_util.h"
#include "mcell_structs.h"

#ifdef DEBUG
#define no_printf printf
#endif


extern struct volume *world;



/**************************************************************************\
 ** Internal utility function section--max/min stuff                     **
\**************************************************************************/

#define max(x,y) ((x)>(y)) ? (x): (y)
#define min(x,y) ((x)<(y)) ? (x): (y)

inline double abs_max_2vec(struct vector3 *v1,struct vector3 *v2)
{
    return max( max( abs(v1->x) , abs(v1->y) ) ,
                max( max( abs(v1->z) , abs(v2->z) ) ,
                     max( abs(v2->y) , abs(v2->x) ) ) );
}


inline double max3(double f1, double f2, double f3)
{
  return (max(f1,max(f2,f3)));
}


inline double min3(double f1, double f2, double f3)
{
  return (min(f1,min(f2,f3)));
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

int edge_hash(struct poly_edge *pe,int nkeys)
{
  unsigned short *a = (unsigned short*) &(pe->v1x);
  unsigned int hashL=1 , hashR=1;
  int i,j;
  for (i=j=0;i<12;i++)
  {
    j += 3;
    if (j>=14) j-=14;
    hashL += ((int)a[i])<<j;
  }
  for (j=0;i<24;i++)
  {
    j += 3;
    if (j>=14) j-=14;
    hashR += ((int)a[i])<<j;
  }
  return ( (hashL ^ hashR) % nkeys );
}


/***************************************************************************
ehtable_init:
  In: pointer to an edge_hashtable struct
      number of keys that the hash table uses
  Out: No return value.  Hash table is initialized.
***************************************************************************/

void ehtable_init(struct edge_hashtable *eht,int nkeys)
{
  int i;
  
  printf("Using %d keys.\n",nkeys);
  eht->nkeys = nkeys;
  eht->stored = 0;
  eht->distinct = 0;
  eht->data = (struct poly_edge*) malloc( nkeys * sizeof(struct poly_edge) );
  if (eht->data == NULL)
  { 
    printf("Malloc failed in ehtable_init\n");
    exit(1);
  }
  
  for (i=0;i<nkeys;i++)
  {
    eht->data[i].next = NULL;
    eht->data[i].n = 0;
    eht->data[i].face1 = eht->data[i].face2 = -1;
  }
}


/***************************************************************************
ehtable_add:
  In: pointer to an edge_hashtable struct
      pointer to the poly_edge to add
  Out: No return value.  Edge is added to the hash table.
***************************************************************************/

void ehtable_add(struct edge_hashtable *eht,struct poly_edge *pe)
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
      printf(" (new)\n");
      return;
    }
    
    if (edge_equals(pep,pe))  /* This edge exists already ... */
    {
      if (pep->face2 == -1)   /* ...and we're the 2nd one */
      {
        pep->face2 = pe->face1;
        pep->edge2 = pe->edge1;
        pep->n++;
        eht->stored++;
        printf(" (join)\n");
        return;
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
        
        pei = (struct poly_edge*) malloc( sizeof(struct poly_edge) );

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
      pei = (struct poly_edge*) malloc( sizeof(struct poly_edge) );

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
surface_net:
  In: array of pointers to walls
      integer length of array
  Out: 1 if the surface is closed, 0 otherwise
  Note: Assumes no triply-connected edges.
***************************************************************************/

int surface_net( struct wall **facelist, int nfaces )
{
  struct poly_edge pe,*pep;
  struct edge *e;
  struct edge_hashtable eht;
  int i,j,k;
  int nedge;
  int nkeys;
  int same;
  int is_closed = 1;
  
  nkeys = (3*nfaces)/2;
  ehtable_init(&eht,nkeys);
  
  for (i=0;i<nfaces;i++)
  {
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
      
      ehtable_add(&eht,&pe);
    }
  }
  
  for (i=0;i<nkeys;i++)
  {
    pep = (eht.data + i);
    while (pep!=NULL)
    {
      if (pep->n > 2)
      {
        printf("Edge with more than two faces attached! Ignoring extra.\n");
        same=1;
      }
      if (pep->n==2 || pep->n==3)
      {
        if (pep->face1 != -1 && pep->face2 != -1)
        {
          facelist[pep->face1]->nb_walls[pep->edge1] = facelist[pep->face2];
          facelist[pep->face2]->nb_walls[pep->edge2] = facelist[pep->face1];
          e = (struct edge*) malloc( sizeof(struct edge) );
          e->forward = facelist[pep->face1];
          e->backward = facelist[pep->face2];
          init_edge_transform(e,pep->edge1);
          facelist[pep->face1]->edges[pep->edge1] = e;
          facelist[pep->face2]->edges[pep->edge2] = e;
        }
      }
      else if (pep->n==1)
      {
        if (!same) is_closed = 0;
        same=0;
        e = (struct edge*) malloc( sizeof(struct edge) );
        e->forward = facelist[pep->face1];
        e->backward = NULL;
        init_edge_transform(e,pep->edge1);
        facelist[pep->face1]->edges[pep->edge1] = e;
      }
      pep = pep->next;
    }
  }
  
  ehtable_kill(&eht);
  return is_closed;
}


/***************************************************************************
init_edge_transform
  In: pointer to an edge
      integer telling the which edge (0-2) of the "forward" face we are
  Out: No return value.  Coordinate transform in edge struct is set.
***************************************************************************/

void init_edge_transform(struct edge *e,int edgenum)
{
  struct vector3 v,*vp;
  struct vector2 ehatf,ehatb;
  int i,j;
  
  if (edgenum+1 < 3)
  { 
    i = edgenum;
    j = i+1; 
  }
  else
  {
    i = edgenum;
    j = 0;
  }
  v.x = e->forward->vert[j]->x - e->forward->vert[i]->x;
  v.y = e->forward->vert[j]->y - e->forward->vert[i]->y;
  v.z = e->forward->vert[j]->z - e->forward->vert[i]->z;
  
  e->length = sqrt(v.x*v.x + v.y*v.y + v.z*v.z);
  e->length_1 = 1 / e->length;
  
  v.x *= e->length_1;
  v.y *= e->length_1;
  v.z *= e->length_1;
  
  if (e->backward == NULL) return;
  
  vp = &(e->forward->unit_u);
  ehatf.u = v.x * vp->x + v.y * vp->y + v.z * vp->z;
  vp = &(e->forward->unit_v);
  ehatf.v = v.x * vp->x + v.y * vp->y + v.z * vp->z;
  vp = &(e->backward->unit_u);
  ehatb.u = v.x * vp->x + v.y * vp->y + v.z * vp->z;
  vp = &(e->backward->unit_v);
  ehatb.v = v.x * vp->x + v.y * vp->y + v.z * vp->z;
  
  e->cos_theta = ehatf.u*ehatb.u + ehatf.v*ehatb.v;
  e->sin_theta = ehatf.v*ehatb.u - ehatf.u*ehatb.v;

/*  
  if (abs(e->cos_theta*e->cos_theta + e->sin_theta*e->sin_theta - 1.0) > EPS_C)
  {
    printf("Linear transformation error.\n");
    exit(1);
  }
*/
  if (e->forward->vert[0] == e->backward->vert[0])
  {
    e->translate.u = e->translate.v = 0.0;
  }
  else
  {
    vp = e->forward->vert[0];
    v.x = e->forward->vert[edgenum]->x - vp->x;
    v.y = e->forward->vert[edgenum]->y - vp->y;
    v.z = e->forward->vert[edgenum]->z - vp->z;
    ehatf.u = v.x*e->forward->unit_u.x + v.y*e->forward->unit_u.y + v.z*e->forward->unit_u.z;
    ehatf.v = v.x*e->forward->unit_v.x + v.y*e->forward->unit_v.y + v.z*e->forward->unit_v.z;
    vp = e->backward->vert[0];
    v.x = e->forward->vert[edgenum]->x - vp->x;
    v.y = e->forward->vert[edgenum]->y - vp->y;
    v.z = e->forward->vert[edgenum]->z - vp->z;
    ehatb.u = v.x*e->forward->unit_u.x + v.y*e->forward->unit_u.y + v.z*e->forward->unit_u.z;
    ehatb.v = v.x*e->forward->unit_v.x + v.y*e->forward->unit_v.y + v.z*e->forward->unit_v.z;
    e->translate.u = -ehatf.u + e->cos_theta*ehatb.u - e->sin_theta*ehatb.v;
    e->translate.v = -ehatf.v + e->sin_theta*ehatb.u + e->cos_theta*ehatb.v;
  }  
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
  WARNING: broken--not picking tiny values properly!
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
  
  e.x = k*v->x;
  e.y = k*v->y;
  e.z = k*v->z;
  
  tiny = (abs_max_2vec( p , &e ) + 1.0) * EPS_C;
  if ( (rng_uint(world->rng_idx) & 1) == 0 ) tiny = -tiny;
  
  v->x -= tiny*f.x;
  v->y -= tiny*f.y;
  v->z -= tiny*f.z;
}


/***************************************************************************
collide_wall:
  In: starting coordinate
      vector to move along
      wall we're checking for a collision
      double to store time of collision
      vector to store the location of the collision
  Out: Integer value indicating what happened
         COLLIDE_MISS  missed
         COLLIDE_FRONT hit the front face (face normal points out of)
         COLLIDE_BACK  hit the back face
         COLLIDE_REDO  hit an edge and modified movement vector; redo
  Note: t and/or hitpt may be modified even if there is no collision
        Not highly optimized yet.  May want to project to Cartesian
        coordinates for speed (as MCell2 did, and Rex implemented
        in pre-40308 backups in vol_utils.c).
***************************************************************************/

int collide_wall(struct vector3 *point,struct vector3 *move,struct wall *face,
                 double *t,struct vector3 *hitpt)
{
  double dp,dv,dd;
  double nx,ny,nz;
  double a,b,c;
  double f,g,h;
  struct vector3 local;
  
  nx = face->normal.x;
  ny = face->normal.y;
  nz = face->normal.z;
  
  dp = nx*point->x + ny*point->y + nz*point->z;
  dv = nx*move->x + ny*move->y + nz*move->z;
  dd = dp - face->d;
  
  if ( (dd>0.0 && dd+dv>0.0) ||  /* Start & end above plane */
       (dd<0.0 && dd+dv<0.0) ||  /* Start & end below plane */
       (dd==0.0 && dv!=0.0) )    /* Start beside plane, end above or below */
  {
    return COLLIDE_MISS;
  }
  
  if (dd+dv == 0.0)
  {
    a = (abs_max_2vec( point , move ) + 1.0) * EPS_C;
    if ((rng_uint(world->rng_idx++)&1)==0) a = -a;
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
  
  *t = a = -dd/dv;
  
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
        jump_away_line(point,move,a,face->vert[1],face->vert[2],&(face->normal));
        return COLLIDE_REDO;
      }
      else return COLLIDE_MISS;
    }
    else if (g == h)
    {
      jump_away_line(point,move,a,face->vert[2],face->vert[0],&(face->normal));
      return COLLIDE_REDO;
    }
    else return COLLIDE_MISS;
  }
  else if (c == 0) /* Hit first edge! */
  {
    jump_away_line(point,move,a,face->vert[0],face->vert[1],&(face->normal));
    return COLLIDE_REDO;
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
  struct vector3 dir;
  struct vector3 *pos;
  double movelen2,d,dirlen2;
  double sigma2;
  
  if ((a->properties->flags & ON_GRID)!=0) return COLLIDE_MISS; /* Should never call on grid molecule! */
  
  if ((a->properties->flags & ON_SURFACE)==0) pos = &( ((struct molecule*)a)->pos );
  else pos = &( ((struct surface_molecule*)a)->pos );
  
  sigma2 = a->properties->radius * a->properties->radius;
  
  dir.x = pos->x - point->x;
  dir.y = pos->y - point->y;
  dir.z = pos->z - point->z;
  
  d = dir.x*move->x + dir.y*move->y + dir.z*move->z;
  
  if (d<0) return COLLIDE_MISS;
  
  movelen2 = move->x*move->x + move->y*move->y + move->z*move->z;
  
  if (d > movelen2) return COLLIDE_MISS;
  
  dirlen2 = dir.x*dir.x + dir.y*dir.y + dir.z*dir.z;
  
  if (movelen2*dirlen2 - d*d > movelen2*sigma2) return COLLIDE_MISS;
  
  return COLLIDE_MOL_M;
}


/***************************************************************************
intersect_box:
  In: llf corner of box
      urb corner of box
      wall we're checking
  Out: Integer value indicating what happened
         0    missed
         1    tri corner inside box
         2    intersect with -x face
  Note: t and/or hitpt may be modified even if there is no collision
        Not highly optimized yet.
***************************************************************************/

#if 0
/* New version, not completed, probably worse than old version.*/
int intersect_box(struct vector3 *llf,struct vector3 *urb,struct wall *w)
{
  int i,j;
  double t,loc;
  struct vector3 wall_llf,wall_urb,vec,hit;
  
  wall_llf.x = GIGANTIC;
  wall_llf.y = GIGANTIC;
  wall_llf.z = GIGANTIC;
  wall_urb.x = -GIGANTIC;
  wall_urb.y = -GIGANTIC;
  wall_urb.z = -GIGANTIC;
  
  for (i=0;i<3;i++)
  {
    if (w->vert[i].x < wall_llf.x) wall_llf.x = w->vert[i].x;
    if (w->vert[i].x > wall_urb.x) wall_urb.x = w->vert[i].x;
    if (w->vert[i].y < wall_llf.y) wall_llf.y = w->vert[i].y;
    if (w->vert[i].y > wall_urb.y) wall_urb.y = w->vert[i].y;
    if (w->vert[i].z < wall_llf.z) wall_llf.z = w->vert[i].z;
    if (w->vert[i].z > wall_urb.z) wall_urb.z = w->vert[i].z;
  }
  
  if (wall_llf.x > urb->x || wall_urb.x < llf->x ||
      wall_llf.y > urb->y || wall_urb.y < llf->y ||
      wall_llf.z > urb->z || wall_urb.z < llf->z)
  {
    return 0;  /* Bounding boxes don't intersect! */
  }
  
  for (i=0;i<3;i++)
  {
    if ( w->vert[i].x < urb->x && w->vert[i].x > llf->x &&
         w->vert[i].y < urb->y && w->vert[i].y > llf->y && 
         w->vert[i].z < urb->z && w->vert[i].z > llf->z)
    {
      return 1;  /* Corner of wall is in box */
    }
  }
  
  for (i=0;i<3;i++)
  {
    j = i+1;
    if (j>=3) j=0;
    
    if ( (w->vert[i].x < llf->x && w->vert[j].x > llf->x) ||
         (w->vert[i].x > llf->x && w->vert[j].x < llf->x) )
    {
      t = (llf->x - w->vert[i].x) / (w->vert[j].x - w->vert[i].x);
      loc = w->vert[i].y + t*(w->vert[j]->y - w->vert[i]->y);
      if (llf->y <= loc && urb->y >= loc)
      {
        loc = w->vert[i].z + t*(w->vert[j]->z - w->vert[j]->z);
        if (llf->z <= loc && urb->y >= loc)
        {
          return 1;  /* Edge of triangle intersects box */
        }
      }
    }
    if ( (w->vert[i].x < urb->x && w->vert[j].x > urb->x) ||
         (w->vert[i].x > urb->x && w->vert[j].x < urb->x) )
    {
      t = (urb->x - w->vert[i].x) / (w->vert[j].x - w->vert[i].x);
      loc = w->vert[i].y + t*(w->vert[j]->y - w->vert[i]->y);
      if (llf->y <= loc && urb->y >= loc)
      {
        loc = w->vert[i].z + t*(w->vert[j]->z - w->vert[j]->z);
        if (llf->z <= loc && urb->y >= loc)
        {
          return 1;  /* Edge of triangle intersects box */
        }
      }
    }
      
    if ( (w->vert[i].y < llf->y && w->vert[j].y > llf->y) ||
         (w->vert[i].y > llf->y && w->vert[j].y < llf->y) )
    {
      t = (llf->y - w->vert[i].y) / (w->vert[j].y - w->vert[i].y);
      loc = w->vert[i].x + t*(w->vert[j]->x - w->vert[i]->x);
      if (llf->x <= loc && urb->x >= loc)
      {
        loc = w->vert[i].z + t*(w->vert[j]->z - w->vert[j]->z);
        if (llf->z <= loc && urb->y >= loc)
        {
          return 1;  /* Edge of triangle intersects box */
        }
      }
    }
    if ( (w->vert[i].y < urb->y && w->vert[j].y > urb->y) ||
         (w->vert[i].y > urb->y && w->vert[j].y < urb->y) )
    {
      t = (urb->y - w->vert[i].y) / (w->vert[j].y - w->vert[i].y);
      loc = w->vert[i].x + t*(w->vert[j]->x - w->vert[i]->x);
      if (llf->x <= loc && urb->x >= loc)
      {
        loc = w->vert[i].z + t*(w->vert[j]->z - w->vert[j]->z);
        if (llf->z <= loc && urb->y >= loc)
        {
          return 1;  /* Edge of triangle intersects box */
        }
      }
    }
      
    if ( (w->vert[i].z < llf->z && w->vert[j].z > llf->z) ||
         (w->vert[i].z > llf->z && w->vert[j].z < llf->z) )
    {
      t = (llf->z - w->vert[i].z) / (w->vert[j].z - w->vert[i].z);
      loc = w->vert[i].y + t*(w->vert[j]->y - w->vert[i]->y);
      if (llf->y <= loc && urb->y >= loc)
      {
        loc = w->vert[i].x + t*(w->vert[j]->x - w->vert[j]->x);
        if (llf->x <= loc && urb->x >= loc)
        {
          return 1;  /* Edge of triangle intersects box */
        }
      }
    }
    if ( (w->vert[i].z < urb->z && w->vert[j].z > urb->z) ||
         (w->vert[i].z > urb->z && w->vert[j].z < urb->z) )
    {
      t = (urb->z - w->vert[i].z) / (w->vert[j].z - w->vert[i].z);
      loc = w->vert[i].y + t*(w->vert[j]->y - w->vert[i]->y);
      if (llf->y <= loc && urb->y >= loc)
      {
        loc = w->vert[i].x + t*(w->vert[j]->x - w->vert[j]->x);
        if (llf->x <= loc && urb->x >= loc)
        {
          return 1;  /* Edge of triangle intersects box */
        }
      }
    }
  }
  
  vec.x = urb->x - llf->x; vec.y = 0.0; vec.z = 0.0;
  if (vec.x != 0.0)
  {
    i = collide_wall(llf,&vec,w,&t,&hit);
    if (i) return 1;  /* Edge of wall intersects triangle */
    vec.x = -vec.x;
    i = collide_wall(urb,&vec,w,&t,&hit);
    if (i) return 1;  /* Edge of wall intersects triangle */
  }

  vec.y = urb->y - llf->y; vec.y = 0.0; vec.z = 0.0;
  if (vec.y != 0.0)
  {
    i = collide_wall(llf,&vec,w,&t,&hit);
    if (i) return 1;  /* Edge of wall intersects triangle */
    vec.y = -vec.y;
    i = collide_wall(urb,&vec,w,&t,&hit);
    if (i) return 1;  /* Edge of wall intersects triangle */
  }

  vec.z = urb->z - llf->z; vec.y = 0.0; vec.z = 0.0;
  if (vec.z != 0.0)
  {
    i = collide_wall(llf,&vec,w,&t,&hit);
    if (i) return 1;  /* Edge of wall intersects triangle */
    vec.z = -vec.z;
    i = collide_wall(urb,&vec,w,&t,&hit);
    if (i) return 1;  /* Edge of wall intersects triangle */
  }
  
  return 0;  /* Near miss! */
}
#endif


void init_tri_wall(struct object *objp, int side, struct vector3 *v0, struct vector3 *v1, struct vector3 *v2)
{
  struct wall *w;
  double f,fx,fy,fz;
  
  w=&objp->walls[side];
  w->next = NULL;
  w->surf_class = world->g_surf;
  w->side = side;
  
/*
  fx = max3( abs(v0->x - v1->x) , abs(v0->x - v2->x) , abs(v1->x - v2->x) );
  fy = max3( abs(v0->y - v1->y) , abs(v0->y - v2->y) , abs(v1->y - v2->y) );
  fz = max3( abs(v0->z - v1->z) , abs(v0->z - v2->z) , abs(v1->z - v2->z) );
  
  if (fx == min3(fx,fy,fz)) w->projection = 0;
  else if (fy == min3(fx,fy,fz)) w->projection = 1;
  else w->projection = 2;
*/
  
  w->vert[0] = v0;
  w->vert[1] = v1;
  w->vert[2] = v2;
  w->edges[0] = NULL;
  w->edges[1] = NULL;
  w->edges[2] = NULL;
  w->nb_walls[0] = NULL;
  w->nb_walls[1] = NULL;
  w->nb_walls[2] = NULL;
  
  w->area = sqrt(0.5 * ( ((v1->x - v0->x)*(v1->x - v0->x)+
                          (v1->y - v0->y)*(v1->y - v0->y)+
                          (v1->z - v0->z)*(v1->z - v0->z)) *
                         ((v2->x - v0->x)*(v2->x - v0->x)+
                          (v2->y - v0->y)*(v2->y - v0->y)+
                          (v2->z - v0->z)*(v2->z - v0->z)) +
                         ((v1->x - v0->x)*(v2->x - v0->x)+
                          (v1->y - v0->y)*(v2->y - v0->y)+
                          (v1->z - v0->z)*(v2->z - v0->z)) *
                         ((v1->x - v0->x)*(v2->x - v0->x)+
                          (v1->y - v0->y)*(v2->y - v0->y)+
                          (v1->z - v0->z)*(v2->z - v0->z)) ) );

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
  
  w->mol = NULL;
  w->mol_count = 0;
  w->effectors = NULL;
  w->viz_state = EXCLUDE_OBJ; 
  w->parent_object = objp;
  w->regions = NULL;
  no_printf("Created wall %d on object %s at:\n",w->side,w->parent_object->sym->name);
  no_printf("  vertex 0: %.9g, %.9g, %.9g\n",w->vert[0]->x,w->vert[0]->y,w->vert[0]->z);
  no_printf("  vertex 1: %.9g, %.9g, %.9g\n",w->vert[1]->x,w->vert[1]->y,w->vert[1]->z);
  no_printf("  vertex 2: %.9g, %.9g, %.9g\n",w->vert[2]->x,w->vert[2]->y,w->vert[2]->z);
}

