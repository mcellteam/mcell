/**************************************************************************\
** File: grid_util.c                                                      **
**                                                                        **
** Purpose: Translates between 3D world coordinates and surface grid index**
**                                                                        **
** Testing status: partially tested (validate_grid_util.c).               **
\**************************************************************************/


#include <math.h>
#include <string.h>
#include <stdlib.h>

#include "rng.h"
#include "grid_util.h"
#include "vol_util.h"
#include "wall_util.h"
#include "mcell_structs.h"


extern struct volume *world;


/*************************************************************************
xyz2uv and uv2xyz:
  In: 2D and 3D vectors and a wall
  Out: first vector is converted to 2nd vector
       WARNING: no error checking--point assumed to be valid!
*************************************************************************/
void xyz2uv(struct vector3 *a,struct wall *w,struct vector2 *b)
{
  if (w->effectors)
  {
    b->u = a->x * w->unit_u.x + a->y * w->unit_u.y + a->z * w->unit_u.z - w->effectors->vert0.u;
    b->v = a->x * w->unit_v.x + a->y * w->unit_v.y + a->z * w->unit_v.z - w->effectors->vert0.v;
  }
  else
  {
    struct vector3 p;
    p.x = a->x - w->vert[0]->x;
    p.y = a->y - w->vert[0]->y;
    p.z = a->z - w->vert[0]->z;
    b->u = p.x * w->unit_u.x + p.y * w->unit_u.y + p.z * w->unit_u.z;
    b->v = p.x * w->unit_v.x + p.y * w->unit_v.y + p.z * w->unit_v.z;
  }
}

void uv2xyz(struct vector2 *a,struct wall *w,struct vector3 *b)
{
  b->x = a->u * w->unit_u.x + a->v * w->unit_v.x + w->vert[0]->x;
  b->y = a->u * w->unit_u.y + a->v * w->unit_v.y + w->vert[0]->y;
  b->z = a->u * w->unit_u.z + a->v * w->unit_v.z + w->vert[0]->z;
}


/*************************************************************************
xyz2grid and uv2grid:
  In: a vector and a surface grid
  Out: int containing the index on the grid of that vector
       WARNING: no error checking--point assumed to be valid!
  Note: xyz2grid just does a dot-product to uv coordinates first.
*************************************************************************/

int xyz2grid(struct vector3 *v,struct surface_grid *g)
{
  struct vector3 *unit_u = &(g->surface->unit_u);
  struct vector3 *unit_v = &(g->surface->unit_v);
  double i,j;
  double u0,u1_u0;
  double striploc,striprem,stripeloc,striperem;
  int strip,stripe,flip;
  
  i = v->x * unit_u->x + v->y * unit_u->y + v->z * unit_u->z - g->vert0.u;
  j = v->x * unit_v->x + v->y * unit_v->y + v->z * unit_v->z - g->vert0.v;
  
  striploc = j*g->inv_strip_wid;
  strip = (int)striploc;
  striprem = striploc - strip;
  
  strip = g->n - strip - 1;

  u0 = j*g->vert2_slope;
  u1_u0 = g->surface->uv_vert1_u - j*g->fullslope;
  
  stripeloc = ((i-u0)/u1_u0)*(((double)strip)+(1.0-striprem));
  stripe = (int)( stripeloc );
  striperem = stripeloc - stripe;
  
  flip = (striperem < 1.0-striprem) ? 0 : 1;

  return strip*strip + 2*stripe + flip;
}

int uv2grid(struct vector2 *v,struct surface_grid *g)
{
  double i,j;
  double u0,u1_u0;
  double striploc,striprem,stripeloc,striperem;
  int strip,stripe,flip;
  
  i = v->u;
  j = v->v;
  
  striploc = j*g->inv_strip_wid;
  strip = (int)striploc;
  striprem = striploc - strip;
  
  strip = g->n - strip - 1;

  u0 = j*g->vert2_slope;
  u1_u0 = g->surface->uv_vert1_u - j*g->fullslope;
  
  stripeloc = ((i-u0)/u1_u0)*(((double)strip)+(1.0-striprem));
  stripe = (int)( stripeloc );
  striperem = stripeloc - stripe;
  
  flip = (striperem < 1.0-striprem) ? 0 : 1;

  return strip*strip + 2*stripe + flip;
}


/*************************************************************************
grid2xyz and grid2uv and grid2uv_random
  In: a surface grid
      index of a tile on that grid
      vector to store the results
  Out: vector contains the coordinates of the center of that tile, or
       a random coordinate within that tile.
       WARNING: no error checking--index assumed to be valid!
  Note: grid2xyz just multiplies by uv unit vectors at the end.
*************************************************************************/

void grid2xyz(struct surface_grid *g,int index,struct vector3 *v)
{
  struct vector3 *unit_u = &(g->surface->unit_u);
  struct vector3 *unit_v = &(g->surface->unit_v);
  int root;
  int rootrem;
  int k,j,i;
  double ucoef,vcoef,over3n;
  
  root = (int)( sqrt((double)index) );
  rootrem = index - root*root;
  k = g->n - root - 1;
  j = rootrem/2;
  i = rootrem - 2*j;
  
  over3n = 1.0 / (double) (3*g->n);
  
  ucoef = ((double)(3*j+i+1))*over3n*g->surface->uv_vert1_u + 
          ((double)(3*k+i+1))*over3n*g->surface->uv_vert2.u;
  vcoef = ((double)(3*k+i+1))*over3n*g->surface->uv_vert2.v;
  
  v->x = ucoef*unit_u->x + vcoef*unit_v->x + g->surface->vert[0]->x;
  v->y = ucoef*unit_u->y + vcoef*unit_v->y + g->surface->vert[0]->y;
  v->z = ucoef*unit_u->z + vcoef*unit_v->z + g->surface->vert[0]->z;
}

void grid2uv(struct surface_grid *g,int index,struct vector2 *v)
{
  int root;
  int rootrem;
  int k,j,i;
  double over3n;
  
  root = (int)( sqrt((double)index) );
  rootrem = index - root*root;
  k = g->n - root - 1;
  j = rootrem/2;
  i = rootrem - 2*j;
  
  over3n = 1.0 / (double) (3*g->n);
  
  v->u = ((double)(3*j+i+1))*over3n*g->surface->uv_vert1_u + 
          ((double)(3*k+i+1))*over3n*g->surface->uv_vert2.u;
  v->v = ((double)(3*k+i+1))*over3n*g->surface->uv_vert2.v;
}

void grid2uv_random(struct surface_grid *g,int index,struct vector2 *v)
{
  int root;
  int rootrem;
  int k,j,i;
  double over_n;
  double u_ran,v_ran;
  
  root = (int)( sqrt((double)index) );
  rootrem = index - root*root;
  k = g->n - root - 1;
  j = rootrem/2;
  i = rootrem - 2*j;
  
  over_n = 1.0 / (double) (g->n);
  
  u_ran = rng_dbl(world->rng);
  v_ran = 1.0 - sqrt(rng_dbl(world->rng));
  
  v->u = ((double)(j+i) + (1-2*i)*(1.0-v_ran)*u_ran)*over_n*g->surface->uv_vert1_u + 
          ((double)(k+i) + (1-2*i)*v_ran)*over_n*g->surface->uv_vert2.u;
  v->v = ((double)(k+i) + (1-2*i)*v_ran)*over_n*g->surface->uv_vert2.v;
}


/*************************************************************************
init_grid_geometry: 
  In: a surface grid with correct # of divisions/edge and wall pointer
  Out: all the precomputed geometry speedup values are properly set
*************************************************************************/

void init_grid_geometry(struct surface_grid *g)
{
  g->inv_strip_wid = 1.0 / (g->surface->uv_vert2.v / ((double)g->n));
  g->vert2_slope = g->surface->uv_vert2.u / g->surface->uv_vert2.v;
  g->fullslope = g->surface->uv_vert1_u / g->surface->uv_vert2.v;

  g->vert0.u = g->surface->vert[0]->x*g->surface->unit_u.x +
               g->surface->vert[0]->y*g->surface->unit_u.y +
               g->surface->vert[0]->z*g->surface->unit_u.z;
  g->vert0.v = g->surface->vert[0]->x*g->surface->unit_v.x +
               g->surface->vert[0]->y*g->surface->unit_v.y +
               g->surface->vert[0]->z*g->surface->unit_v.z;
  
  g->n_tiles = g->n * g->n;
}


/*************************************************************************
create_grid: 
  In: a wall pointer that needs to have its grid created
      a guess for the subvolume the center of the grid is in
  Out: integer, 0 if grid exists or was created, 1 on memory error.
       The grid is created and the wall is set to point at it.
*************************************************************************/

int create_grid(struct wall *w,struct subvolume *guess)
{
  struct surface_grid *sg;
  struct vector3 center;
  int i;

  if (w->effectors != NULL) return 0;
  
  sg = (struct surface_grid *) mem_get(w->birthplace->effs);
  if (sg == NULL) return 1;
  
  center.x = 0.33333333333*(w->vert[0]->x + w->vert[1]->x + w->vert[2]->x);
  center.y = 0.33333333333*(w->vert[0]->y + w->vert[1]->y + w->vert[2]->y);
  center.z = 0.33333333333*(w->vert[0]->z + w->vert[1]->z + w->vert[2]->z);
  
  sg->surface = w;
  sg->subvol = find_subvolume(&center , guess);
  
  /* sg->n = sqrt( w->area ); */
  /* sg->n = 0.5 + sqrt( w->area ); */
  sg->n = (int) ceil(sqrt( w->area ));
  if (sg->n<1) sg->n=1;

  sg->n_tiles = sg->n * sg->n;
  sg->n_occupied = 0;

  sg->binding_factor = ((double)sg->n_tiles) / w->area;
  init_grid_geometry(sg);
  
  sg->mol = (struct grid_molecule**)malloc(sg->n_tiles*sizeof(struct grid_molecule*));
  if (sg->mol == NULL) return 1;

  for (i=0;i<sg->n_tiles;i++) sg->mol[i] = NULL;
  
  w->effectors = sg;

  return 0;
}



/*************************************************************************
nearest_free: 
  In: a surface grid
      a vector in u,v coordinates on that surface
      the maximum distance we can search for free spots
  Out: integer containing the index of the closest unoccupied grid point
       to the vector, or -1 if no unoccupied points are found in range
  Note: we assume you've already checked the grid element that contains
        the point, so we don't bother looking there first.
*************************************************************************/

int nearest_free(struct surface_grid *g,struct vector2 *v,double max_d2,double *found_dist2)
{
  int h,i,j,k;
  int span;
  int can_flip;
  int idx;
  double d2;
  double f,ff,fff;
  double over3n = 0.333333333333333/(double)(g->n);
  
  idx = -1;
  d2 = 2*max_d2 + 1.0;
  
  for (k=0;k<g->n;k++)
  {
    f = v->v - ((double)(3*k+1))*over3n*g->surface->uv_vert2.v;
    ff = f - over3n*g->surface->uv_vert2.v;
    ff *= ff;
    f *= f;
    if (f > max_d2 && ff > max_d2) continue;  /* Entire strip is too far away */
    
    span = (g->n - k);
    for (j=0 ; j < span ; j++)
    {
      can_flip = (j!=span-1);
      for (i=0 ; i <= can_flip ; i++)
      {
	fff = v->u - over3n*( (double)(3*j+i+1)*g->surface->uv_vert1_u + (double)(3*k+i+1)*g->surface->uv_vert2.u );
	fff *= fff;
	if (i) fff += ff;
	else fff += f;
	
	if (fff<max_d2 && (idx==-1 || fff<d2))
	{
	  h = (g->n-k)-1;
	  h = h*h + 2*j + i;
	  
	  if (g->mol[h]==NULL)
	  {
	    idx = h;
	    d2 = fff;
	  }
	  else if (idx==-1)
	  {
	    if (fff < d2) d2=fff;
	  }
	}
      }
    }
  }
 
  if (found_dist2!=NULL) *found_dist2 = d2;
  
  return idx;
}


/*************************************************************************
search_nbhd_for_free: 
  In: the wall that we ought to be in
      a vector in u,v coordinates on that surface where we should go
      the maximum distance we can search for free spots
      a place to store the index of our free slot
      a function that we'll call to make sure a wall is okay
      context for that function passed in by whatever called us
  Out: pointer to the wall that has the free slot, or NULL if no wall
       exist in range.
  Note: This is not recursive.  It should be made recursive.
*************************************************************************/

struct wall *search_nbhd_for_free(struct wall *origin,struct vector2 *point,double max_d2,int *found_idx,
                             int (*ok)(void*,struct wall*),void *context)
{
  struct wall *there;
  int i,j;
  double d2;
  struct vector2 pt,ed;
  struct vector2 vurt0,vurt1;
  int best_i;
  double best_d2;
  struct wall *best_w;
  
  best_i = -1;
  best_d2 = 2.0*max_d2+1.0;
  best_w = NULL;
  
  if (origin->effectors==NULL)
  {
    if (create_grid(origin,NULL)) return NULL;  /* FIXME: handle out of memory properly */
  }
  
  i = nearest_free(origin->effectors,point,max_d2,&d2);

  if (i != -1)
  {
    best_i = i;  
    best_d2 = d2;
    best_w = origin;
  }
  
  for (j=0 ; j<3 ; j++)
  {
    if (origin->edges[j]==NULL || origin->edges[j]->backward==NULL) continue;
    
    if (origin->edges[j]->forward==origin) there = origin->edges[j]->backward;
    else there = origin->edges[j]->forward;
    
    if (ok!=NULL && !(*ok)(context,there) ) continue;  /* Calling function doesn't like this wall */
    
    switch (j)
    {
      case 0:
        vurt0.u = vurt0.v = 0.0;
	vurt1.u = origin->uv_vert1_u; vurt1.v = 0;
	break;
      case 1:
      	vurt0.u = origin->uv_vert1_u; vurt0.v = 0;
        memcpy(&vurt1,&(origin->uv_vert2),sizeof(struct vector2));
	break;
      case 2:
        memcpy(&vurt0,&(origin->uv_vert2),sizeof(struct vector2));
	vurt1.u = vurt1.v = 0.0;
        break;
      /* No default case since 0<=j<=2 */
    }
    ed.u = vurt1.u - vurt0.u;
    ed.v = vurt1.v - vurt0.v;
    pt.u = point->u - vurt0.u;
    pt.v = point->v - vurt0.v;
    
    d2 = pt.u*ed.u + pt.v*ed.v;
    d2 = (pt.u*pt.u + pt.v*pt.v) - d2*d2/(ed.u*ed.u+ed.v*ed.v); /* Distance squared to line */
    
    if (d2<best_d2)
    {
      if (there->effectors==NULL)
      {
	if (create_grid(there,NULL)) return NULL;  /* FIXME: handle out of memory properly */
      }
      traverse_surface(origin,point,j,&pt);
      i = nearest_free(there->effectors,&pt,max_d2,&d2);
      
      if (i!=-1 && d2 < best_d2)
      {
	best_i = i;
	best_d2 = d2;
	best_w = there;
      }
    }
  }
  
  *found_idx = best_i;
  return best_w;
} 
