/**************************************************************************\
** File: grid_util.c                                                      **
**                                                                        **
** Purpose: Translates between 3D world coordinates and surface grid index**
**                                                                        **
** Testing status: partially tested (validate_grid_util.c).               **
\**************************************************************************/


#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include "logging.h"
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
  if (w->grid)
  {
    b->u = a->x * w->unit_u.x + a->y * w->unit_u.y + a->z * w->unit_u.z - w->grid->vert0.u;
    b->v = a->x * w->unit_v.x + a->y * w->unit_v.y + a->z * w->unit_v.z - w->grid->vert0.v;
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

void grid2xyz(struct surface_grid *g,int idx,struct vector3 *v)
{
  struct vector3 *unit_u = &(g->surface->unit_u);
  struct vector3 *unit_v = &(g->surface->unit_v);
  int root;
  int rootrem;
  int k,j,i;
  double ucoef,vcoef,over3n;
  
  root = (int)( sqrt((double) idx) );
  rootrem = idx - root*root;
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

void grid2uv(struct surface_grid *g,int idx,struct vector2 *v)
{
  int root;
  int rootrem;
  int k,j,i;
  double over3n;
  
  root = (int)( sqrt((double) idx) );
  rootrem = idx - root*root;
  k = g->n - root - 1;
  j = rootrem/2;
  i = rootrem - 2*j;
  
  over3n = 1.0 / (double) (3*g->n);
  
  v->u = ((double)(3*j+i+1))*over3n*g->surface->uv_vert1_u + 
          ((double)(3*k+i+1))*over3n*g->surface->uv_vert2.u;
  v->v = ((double)(3*k+i+1))*over3n*g->surface->uv_vert2.v;
}

void grid2uv_random(struct surface_grid *g,int idx,struct vector2 *v)
{
  int root;
  int rootrem;
  int k,j,i;
  double over_n;
  double u_ran,v_ran;
  
  root = (int)( sqrt((double) idx) );
  rootrem = idx - root*root;
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
  struct surface_grid *sg = NULL;
  struct vector3 center;

  if (w->grid != NULL) return 0;
  
  sg = (struct surface_grid *) CHECKED_MEM_GET(w->birthplace->grids, "surface grid");
  if (sg == NULL) return 1;
  
  center.x = 0.33333333333*(w->vert[0]->x + w->vert[1]->x + w->vert[2]->x);
  center.y = 0.33333333333*(w->vert[0]->y + w->vert[1]->y + w->vert[2]->y);
  center.z = 0.33333333333*(w->vert[0]->z + w->vert[1]->z + w->vert[2]->z);
  
  sg->surface = w;
  sg->subvol = find_subvolume(&center , guess);
  
  sg->n = (int) ceil(sqrt( w->area ));
  if (sg->n<1) sg->n=1;

  sg->n_tiles = sg->n * sg->n;
  sg->n_occupied = 0;

  sg->binding_factor = ((double)sg->n_tiles) / w->area;
  init_grid_geometry(sg);
  
  sg->mol = CHECKED_MALLOC_ARRAY(struct grid_molecule *,
                                 sg->n_tiles,
                                 "surface grid");

  for (unsigned int i=0;i<sg->n_tiles;i++) sg->mol[i] = NULL;
  
  w->grid = sg;

  return 0;
}



/*************************************************************************
grid_neighbors: 
  In: a surface grid
      an index on that grid
      flag that tells whether we have to create a grid on the 
      neighbor wall if there is no grid there
      an array[3] of pointers to be filled in with neighboring grid(s)
      an array[3] of pointers to be filled in with neighboring indices
  Out: no return value.  The three nearest neighbors are returned,
       which may be on a neighboring grid if the supplied index is
       at an edge.  If there is no neighbor in one of the three
       directions, the neighboring grid pointer is set to NULL.
  Note: the three neighbors are returned in the same order as the
        edges, i.e. the 0th will be the nearest neighbor in the
        direction of the 0th edge, and so on.
  Note: If this code is used to find neighboring molecules,
        the "create_grid_flag" should be set to zero. 
        In such case if a nearby wall exists but has no grid placed
	on it, this function returns NULL for that grid, even though
	there is space there (just no molecules).
        If this code is used to find free spots, the "create_grid_flag" 
        should be set to 1 (or any positive value) and the function
        returns newly created grid for this wall.
*************************************************************************/

void grid_neighbors(struct surface_grid *grid,
                    int idx,
                    int create_grid_flag,
                    struct surface_grid **nb_grid,
                    int *nb_idx)
{
  int i,j,k,root,rootrem;
  struct vector3 loc_3d;
  struct vector2 near_2d;
  double d;
  
  /* Calculate strip (k), stripe (j), and flip (i) indices from idx */
  root = (int)( sqrt((double) idx) );
  rootrem = idx - root*root;
  k = root;
  j = rootrem/2;
  i = rootrem - 2*j;
  
  /* First look "left" (towards edge 2) */
  if (j>0 || i>0) /* all tiles except upright tiles in stripe 0 */
  {
    nb_grid[2] = grid;
    nb_idx[2] = idx-1;
  }
  else /* upright tiles in stripe 0 */
  {
    if (grid->surface->nb_walls[2]==NULL) nb_grid[2] = NULL;
    else if ((grid->surface->nb_walls[2]->grid==NULL) && (!create_grid_flag)) nb_grid[2] = NULL;  
    else
    {
      if ((grid->surface->nb_walls[2]->grid==NULL) && create_grid_flag) create_grid(grid->surface->nb_walls[2], NULL);

      if (grid->mol[idx]!=NULL) uv2xyz(&grid->mol[idx]->s_pos,grid->surface,&loc_3d);
      else grid2xyz(grid,idx,&loc_3d);
      d = closest_interior_point(&loc_3d,grid->surface->nb_walls[2],&near_2d,GIGANTIC);
      if (d==GIGANTIC) nb_grid[2]=NULL;
      else
      {
	nb_grid[2] = grid->surface->nb_walls[2]->grid;
	nb_idx[2] = uv2grid(&near_2d,nb_grid[2]);
      }
    }
  }
  
  /* Then "right" (towards edge 1) */
  if (j < k) /* all tiles except upright tiles in last stripe */
  {
    nb_grid[1] = grid;
    nb_idx[1] = idx+1;
  }
  else  /* upright tiles in last stripe */
  {
    if (grid->surface->nb_walls[1]==NULL) nb_grid[1] = NULL;
    else if ((grid->surface->nb_walls[1]->grid==NULL) && (!create_grid_flag)) nb_grid[1] = NULL;
    else
    {
      if ((grid->surface->nb_walls[1]->grid==NULL) && create_grid_flag) create_grid(grid->surface->nb_walls[1], NULL);

      if (grid->mol[idx]!=NULL) uv2xyz(&grid->mol[idx]->s_pos,grid->surface,&loc_3d);
      else grid2xyz(grid,idx,&loc_3d);
      d = closest_interior_point(&loc_3d,grid->surface->nb_walls[1],&near_2d,GIGANTIC);
      if (d==GIGANTIC) nb_grid[1]=NULL;
      else
      {
	nb_grid[1] = grid->surface->nb_walls[1]->grid;
	nb_idx[1] = uv2grid(&near_2d,nb_grid[1]);
      }
    }
  }

    
  /* Finally "up/down" (towards edge 0 if not flipped) */
  if (i || k+1 < grid->n)  /* all tiles except upright tiles in last strip */
  {
    nb_grid[0] = grid;
    if (i) nb_idx[0] = 2*j+(k-1)*(k-1);    /* unflip and goto previous strip */
    else nb_idx[0] = 1 + 2*j+(k+1)*(k+1);  /* flip and goto next strip */
  }
  else  /* upright tiles in last strip */
  {
    if (grid->surface->nb_walls[0]==NULL) nb_grid[0] = NULL;
    else if ((grid->surface->nb_walls[0]->grid==NULL) && (!create_grid_flag)) nb_grid[0] = NULL;
    else
    {
      if ((grid->surface->nb_walls[0]->grid==NULL) && create_grid_flag) create_grid(grid->surface->nb_walls[0], NULL);

      if (grid->mol[idx]!=NULL) uv2xyz(&grid->mol[idx]->s_pos,grid->surface,&loc_3d);
      else grid2xyz(grid,idx,&loc_3d);
      d = closest_interior_point(&loc_3d,grid->surface->nb_walls[0],&near_2d,GIGANTIC);
      if (d==GIGANTIC) nb_grid[0]=NULL;
      else
      {
	nb_grid[0] = grid->surface->nb_walls[0]->grid;
	nb_idx[0] = uv2grid(&near_2d,nb_grid[0]);
      }
    }
  }
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
  Note: if no unoccupied tile is found, found_dist2 contains distance to
        closest occupied tile.
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

  /* check whether the grid is fully occupied */
  if(g->n_occupied >= g->n_tiles){
    *found_dist2 = 0;
     return -1;
  }
  
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
  Note: usually the calling function will create a grid if needed and
        check the grid element at u,v; if that is not done this function
        will return the correct result but not efficiently.
  Note: This is not recursive.  It should be made recursive.
*************************************************************************/

struct wall *search_nbhd_for_free(struct wall *origin,struct vector2 *point,double max_d2,int *found_idx,
                             int (*ok)(void*,struct wall*),void *context)
{
  struct wall *there = NULL;
  int i, j;
  double d2;
  struct vector2 pt,ed;
  struct vector2 vurt0,vurt1;
  int best_i;
  double best_d2;
  struct wall *best_w = NULL;
  
  best_i = -1;
  best_d2 = 2.0*max_d2+1.0;
  
  if (origin->grid == NULL  &&  create_grid(origin, NULL))
    mcell_allocfailed("Failed to create grid for wall.");
  
  i = -1; /* default return value */

  /* Find index and distance of nearest free grid element on origin wall */
  if(origin->grid->n_occupied < origin->grid->n_tiles){
     i = nearest_free(origin->grid, point, max_d2,&d2);
  }

  if (i != -1)
  {
    best_i = i;  
    best_d2 = d2;
    best_w = origin;
  }
 
  /* if there are no free slots on the origin wall - look around */

  if(best_w == NULL)
  {
     /* Check for closer free grid elements on neighboring walls */
     for (j=0 ; j<3 ; j++)
     {
       if (origin->edges[j]==NULL || origin->edges[j]->backward==NULL) continue;
    
       if (origin->edges[j]->forward==origin) there = origin->edges[j]->backward;
       else there = origin->edges[j]->forward;
    
       if (ok!=NULL && !(*ok)(context,there) ) continue;  /* Calling function doesn't like this wall */
  
       /* check whether there are any available spots on the neighbor wall */   
         if (there->grid != NULL)
         {
             if(there->grid->n_occupied >= there->grid->n_tiles){
                   continue;
             }
         }
 
       /* Calculate distance between point and edge j of origin wall */
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
         default:
           /* default case should not occur since 0<=j<=2 */
           UNHANDLED_CASE(j);
       }
       ed.u = vurt1.u - vurt0.u;
       ed.v = vurt1.v - vurt0.v;
       pt.u = point->u - vurt0.u;
       pt.v = point->v - vurt0.v;
    
       d2 = pt.u*ed.u + pt.v*ed.v;
       d2 = (pt.u*pt.u + pt.v*pt.v) - d2*d2/(ed.u*ed.u+ed.v*ed.v); /* Distance squared to line */
    
       /* Check for free grid element on neighbor if point to edge distance is closer than best_d2  */
       if (d2<best_d2)
       {
               
         if (there->grid == NULL  &&  create_grid(there, NULL))
           mcell_allocfailed("Failed to create grid for wall.");

         traverse_surface(origin,point,j,&pt);
         i = nearest_free(there->grid,&pt,max_d2,&d2);
         
         if (i!=-1 && d2 < best_d2)
         {
	   best_i = i;
	   best_d2 = d2;
	   best_w = there;
         }
       }
     }
  }  

  *found_idx = best_i;
  return best_w;
} 


/*************************************************************************
grid_release_check: 
  In: a region release data structure
      index of object in that structure we're testing
      index of wall in that object
      index on surface grid of that wall
      expression to evaluate (or NULL if we are to evaluate full
        expression).
  Out: 1 if the identified grid element meets the criteria of the region
       release, 0 if not.
  Note: Only call this if rrd->refinement is set, meaning that you have
        to pick subsets of 2d regions based on clipping by 3d volumes.
	If you are just working on 2d, rrd already contains bit arrays
	that say whether or not a given wall is OK.
  Note: This function is recursive.
*************************************************************************/

int grid_release_check(struct release_region_data *rrd,
                       int obj_n,
                       int wall_n,
                       int grid_n,
                       struct release_evaluator *expr)
{
  struct region *r;
  int okL,okR;
  if (expr==NULL) expr = rrd->expression;
  if (expr->left==NULL) return 0;
  if (expr->op&REXP_LEFT_REGION)
  {
    r = (struct region*)expr->left;
    if (r->parent != rrd->owners[obj_n]) okL = 0;
    else okL = get_bit(r->membership,wall_n);
  }
  else okL = grid_release_check(rrd,obj_n,wall_n,grid_n,expr->left);
  if (expr->right==NULL) return okL;
  
  if (expr->op&(REXP_SUBTRACTION|REXP_INTERSECTION|REXP_INCLUSION) && !okL) return 0;  /* Don't need to check right */
  else if (expr->op&REXP_UNION && okL) return 1;  /* Don't need to check right */
  
  if (expr->op&REXP_INCLUSION)
  {
    struct wall *w = rrd->owners[obj_n]->wall_p[wall_n];
    struct vector3 pt;
    grid2xyz(w->grid,grid_n,&pt);
    okR = surface_point_in_region(rrd->owners[obj_n],wall_n,&pt,expr->right);
  }
  else if (expr->op&REXP_RIGHT_REGION)
  {
    r = (struct region*)expr->right;
    if (r->parent != rrd->owners[obj_n]) okR = 0;
    else okR = get_bit(r->membership,wall_n);
  }
  else okR = grid_release_check(rrd,obj_n,wall_n,grid_n,expr->right);
  if (expr->op&REXP_UNION) return okL || okR;
  else if (expr->op&REXP_SUBTRACTION) return okL && !okR;
  else if (expr->op&(REXP_INTERSECTION|REXP_INCLUSION)) return okL && okR;
  
  return 0;
}

/***************************************************************************
delete_tile_neighbor_list:
   In: linked list of tile_neighbors
   Out: none.  The memory is freed
****************************************************************************/
void delete_tile_neighbor_list(struct tile_neighbor *head)
{
   struct tile_neighbor *nnext;
   while(head != NULL)
   {
      nnext = head->next;
      free(head);
      head = nnext;
   }
}



/***************************************************************************
push_tile_neighbor_to_list:
   In: head of the linked list
       surface_grid of the wall the tile is on
       index of the tile
   Out: none. The linked list is expanded by one node (grid/idx).
****************************************************************************/
void push_tile_neighbor_to_list(struct tile_neighbor **head, struct surface_grid *grid, int idx)
{
    struct tile_neighbor *tile_nbr, *old_head;
    
    old_head = *head;

    tile_nbr = CHECKED_MALLOC_STRUCT(struct tile_neighbor, "tile_neighbor");
    tile_nbr->grid = grid;
    tile_nbr->flag = 0;
    tile_nbr->idx = idx;

    if(old_head == NULL)
    {
        tile_nbr->next = NULL;
        old_head = tile_nbr;
    }else{
        tile_nbr->next = old_head;
        old_head = tile_nbr;
    }

    *head = old_head;

}

/*********************************************************************
get_tile_neighbor_from_list_of_vacant_neighbors:
   In:  head of the linked list of tile_neighbors
        index of the node in the list (indexing starts from zero)
        surface_grid (return value)
        tile index (return value)
   Out: Only if the node was not previously selected, the surface_grid 
        and tile index on the surface grid are set.
        Returns number of really vacant tiles at the start of the function 
        (some tiles may be already selected by the previous calls 
           to the function).
*********************************************************************/
int get_tile_neighbor_from_list_of_vacant_neighbors(struct tile_neighbor *head, int list_index, struct surface_grid **grid, int *tile_idx)
{

   struct tile_neighbor *curr = head;
   
   int iter = 0; /* iterator through the linked list like through the array */
   int count = 0; /* number of really vacant tiles */

   while(curr != NULL)
   {
      if((curr->flag & TILE_CHECKED) == 0) count++;
      if((iter == list_index) && ((curr->flag & TILE_CHECKED) == 0)) {
         curr->flag |= TILE_CHECKED;   
         *grid = curr->grid;
          *tile_idx = curr->idx;
      }
      iter++;
      curr = curr->next;
   }
   return count;  
}

/*************************************************************************
get_tile_vertices:
   In: surface grid
       tile index
       tile flip information (return value)
       first tile vertex R   (return value)
       second tile vertex S  (return value)
       third tile vertex T   (return value)
   Out: the tile vertices (R,S,T) coordinates are defined
   Note: the vertices (R,S,T) are listed clockwise.  For the upright tile 
         (orientation similar to the wall) two vertices R and S are
         on the same line PQ parallel and closest to the u-axis.  
         For the inverted tile (orientation inverted compared to the wall) 
         two vertices S and T are on the same line XY parallel and furthest 
         to the u-axis.

*************************************************************************/
void get_tile_vertices(struct surface_grid *sg, int idx, int *flp, struct vector2 *R, struct vector2 *S, struct vector2 *T)
{
   /* indices of the barycentric subdivision */
   int strip, stripe, flip; 
   int root, rootrem;
   struct vector2 P, Q, X, Y;
   /* length of the segments PQ and XY (see below) */
   double pq, xy; 
   /* cotangent of the angle formed between the u-axis
      and the wall edge opposite to the origin of the 
      uv-coordinate system */
   double cot_angle;

  /* Calculate strip, stripe, and flip indices from idx */
   root = (int)( sqrt((double) idx) );
   rootrem = idx - root*root;
   strip = sg->n - root - 1;
   stripe = rootrem/2;
   flip = rootrem - 2*stripe;

  
  /* Let PQ to be the segment on the grid containing the vertex R
     and the one closest to the u-axis.  Let XY to be the segment
     containing the vertex T and the one furthest to the u-axis.
     Let point P to be on the left side of point Q.
     Let point X to be on the left side of point Y.
  */


   /* find v-coordinates of P, Q, X, Y */
   P.v = Q.v = strip/sg->inv_strip_wid;
   X.v = Y.v = (strip + 1)/sg->inv_strip_wid;

   /* find u-coordinates of P, Q, X, Y */
   P.u = (P.v)*(sg->vert2_slope);
   X.u = (X.v)*(sg->vert2_slope);

   cot_angle = (sg->surface->uv_vert1_u - sg->surface->uv_vert2.u)/(sg->surface->uv_vert2.v);
   Q.u = sg->surface->uv_vert1_u - (Q.v)* cot_angle;
   Y.u = sg->surface->uv_vert1_u - (Y.v)* cot_angle;
 
   pq = Q.u - P.u;
   if(idx == 0) {
      xy = 0;
   }else{
      xy = Y.u - X.u;
   }

   /* find coordinates of the tile vertices */
   /* For the upright tile the vertices R and S lie on
      the line PQ and vertex T - on the line XY.
      For the inverted tile only the vertex R lies on the line PQ
      while the vertices S and T lie on the line XY */

    if(flip == 1)
    {
       /* inverted tile */
       R->v = P.v;
       /* note: pq/(sg->n - strip) tells us 
                the number of slots on the line PQ */
       R->u = P.u + pq*(stripe + 1)/(sg->n - strip);

       S->v = T->v = X.v;
       T->u = X.u + xy * stripe / (sg->n - strip -1);
       S->u = X.u + xy * (stripe + 1) / (sg->n - strip - 1);

    }else{
      /* upright tile */
       R->v = S->v = P.v;
       T->v = X.v;

       /* note: pq/(sg->n - strip) tells us 
                the number of slots on the line PQ */
       R->u = P.u + pq * stripe / (sg->n - strip);
       S->u = P.u + pq * (stripe + 1) / (sg->n - strip);
       if(idx == 0){
          T->u = X.u;
       }else{
          T->u = X.u + xy * stripe / (sg->n - strip - 1);
       }

    }

    /* set the tile flip value */
    *flp = flip;
} 

/*************************************************************************
tile_orientation:
  In: a vector to the point on the grid and a surface grid
  Out: 0 if the triangle containg the point  is upright, 
       and 1 if it is inverted.
       WARNING: no error checking--point assumed to be valid.
*************************************************************************/
int tile_orientation(struct vector2 *v,struct surface_grid *g)
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

  return flip;
}



/*************************************************************************
get_tile_area:
  In: vertices of the tile
  Out: area of the tile
**************************************************************************/
double get_tile_area(struct vector2 *A, struct vector2 *B, struct vector2 *C)
{

  return 0.5*(-(B->u)*(A->v) + (C->u)*(A->v) + (A->u)*(B->v) - (C->u)*(B->v) - (A->u)*(C->v) + (B->u)*(C->v));

}


/*****************************************************************************
grid_all_neighbors_for_inner_tile:
  In: a surface grid
      an index on that grid
      a point on that tile
      a linked list of  neighbor tiles (return value)
      a length of the linked list above (return value)
  Out: The list and list length of nearest neighbors 
       are returned.  Neighbors should share either common edge or 
       common vertice.
  Note: The code below is valid only for the inner tile - the one that has 12 
        neighbor tiles all belonging to the same grid as the start tile.
******************************************************************************/
void grid_all_neighbors_for_inner_tile(struct surface_grid *grid, int idx, struct vector2 *pos,  struct tile_neighbor **tile_neighbor_head, int *list_length)
{
   struct tile_neighbor *tile_nbr_head = NULL;  
   int count = 0, tile_orient, temp_ind = -1;               

   int vert_nbr_ind = -1;

   /* The tile that has shape similar to the wall is called upright tile.
      The tile that has shape inverted relative to the wall is 
      called inverted tile. */ 

   int kk; 
   struct surface_grid *sg[3];/* Neighboring surface grids (edge-to-edge) */
   int si[3];  /* Indices on those grids (edge-to-edge) of neighbor molecules */

   if((u_int)idx >= grid->n_tiles){
      mcell_internal_error("Grid molecule tile index is greater than or equal of the number of tiles on the grid\n");
   }

   for(kk = 0; kk < 3; kk++)
   {
     sg[kk] = NULL;
     si[kk] = -1;
   }

   /* find neighbors to react with */
   grid_neighbors(grid, idx, 0, sg, si);

   if((grid != sg[0]) || (grid != sg[1]) || (grid != sg[2]))
   {
      mcell_internal_error("The function 'grid_all_neighbors_for_inner_tile()' is called for the tile %d that is not an inner tile.", idx);
   }

   for(kk = 0; kk < 3; kk++)
   {
      if((si[kk] != idx -1) && (si[kk] != idx + 1))
      {
         vert_nbr_ind = si[kk];
         break;
      }
   }

   /* The tile has 2 neighbors to the left and 2 neighbors to the right */
    push_tile_neighbor_to_list(&tile_nbr_head, grid, idx -1);
    count++;
    push_tile_neighbor_to_list(&tile_nbr_head, grid, idx - 2);
    count++;
    push_tile_neighbor_to_list(&tile_nbr_head, grid, idx + 1);
    count++;
    push_tile_neighbor_to_list(&tile_nbr_head, grid, idx + 2);
    count++;


    /* find the orientation of the tile */
    tile_orient = tile_orientation(pos, grid);

    if(tile_orient == 0)
    {
      /* upright tile has 5 neighbors in the row above it */
      push_tile_neighbor_to_list(&tile_nbr_head, grid, vert_nbr_ind);
      count++;
      push_tile_neighbor_to_list(&tile_nbr_head, grid, vert_nbr_ind - 1);
      count++;
      push_tile_neighbor_to_list(&tile_nbr_head, grid, vert_nbr_ind - 2);
      count++;
      push_tile_neighbor_to_list(&tile_nbr_head, grid, vert_nbr_ind + 1);
      count++;
      push_tile_neighbor_to_list(&tile_nbr_head, grid, vert_nbr_ind + 2);
      count++;

      /* upright tile has 3 neighbors in the row below it */
      temp_ind = move_strip_down(grid, idx);
      if(temp_ind == -1)
      {
         mcell_internal_error("The function 'grid_all_neighbors_for_inner_tile() is called for the tile %d that is not an inner tile.", idx);
      }
      push_tile_neighbor_to_list(&tile_nbr_head, grid, temp_ind);
      count++;
      push_tile_neighbor_to_list(&tile_nbr_head, grid, temp_ind - 1);
      count++;
      push_tile_neighbor_to_list(&tile_nbr_head, grid, temp_ind + 1);
      count++;
    }else{
        /* inverted tile has 3 neighbors in the row above it  */
      temp_ind = move_strip_up(grid, idx);
      if(temp_ind == -1)
      {
         mcell_internal_error("The function 'grid_all_neighbors_for_inner_tile() is called for the tile %d that is not an inner tile.", idx);
      }
      push_tile_neighbor_to_list(&tile_nbr_head, grid, temp_ind);
      count++;
      push_tile_neighbor_to_list(&tile_nbr_head, grid, temp_ind - 1);
      count++;
      push_tile_neighbor_to_list(&tile_nbr_head, grid, temp_ind + 1);
      count++;

       /*   inverted tile has 5 neighbors in the row below it  */
      push_tile_neighbor_to_list(&tile_nbr_head, grid, vert_nbr_ind);
      count++;
      push_tile_neighbor_to_list(&tile_nbr_head, grid, vert_nbr_ind - 1);
      count++;
      push_tile_neighbor_to_list(&tile_nbr_head, grid, vert_nbr_ind - 2);
      count++;
      push_tile_neighbor_to_list(&tile_nbr_head, grid, vert_nbr_ind + 1);
      count++;
      push_tile_neighbor_to_list(&tile_nbr_head, grid, vert_nbr_ind + 2);
      count++;
    }

    if(count != 12)
    {
       mcell_internal_error("The function 'grid_all_neighbors_for_inner_tile() is called for the tile %d that is not an inner tile.", idx);
    }else{
       *list_length = count;
    }
    
    *tile_neighbor_head = tile_nbr_head;

}

/**************************************************************************
grid_all_neighbors_across_walls_through_vertices: 
  In: a surface grid
      tile index on that grid
      linked list of the neighbor walls that share one vertex only
      flag that tells whether we need to create a grid on a neighbor wall
      a linked list of  neighbor tiles (return value)
      a length of the linked list above (return value)
  Out: The list of nearest neighbors are returned,
       Neighbors should share either common edge or common vertice.
  Note: This version allows looking for the neighbors at the neighbor walls
       that are connected to the start wall through vertices only.
****************************************************************************/
void grid_all_neighbors_across_walls_through_vertices(struct surface_grid *grid, int idx, struct wall_list *wall_nbr_head, int create_grid_flag, struct tile_neighbor **tile_neighbor_head, int *list_length)
{
   struct tile_neighbor *tile_nbr_head = NULL; 
   struct wall * origin; /* wall where the start grid is located */
   struct  vector3 *vert_0, *vert_1, *vert_2; /* vertices of the origin wall */
   struct wall_list *wl;
   struct wall *w;

   /* only one corner tile from each neighbor wall
      can be a neighbor to our start tile */
   int tile_idx_0, tile_idx_mid, tile_idx_last;

   int tiles_count = 0; /* number of tiles added */

   origin = grid->surface;
   vert_0 = origin->vert[0];
   vert_1 = origin->vert[1];
   vert_2 = origin->vert[2];

   /* since the neighbor walls are connected to the origin wall by just
      one vertex, and this code is valid only for the corner tile
      on the origin wall, from each neighbor wall we will pick up
      only one corner tile that shared a vertex with the origin wall */
   for(wl = wall_nbr_head; wl != NULL; wl = wl->next)
   {
       w = wl->this_wall;
       if(w->grid == NULL)
       {
          if(create_grid_flag)
          {
             if(create_grid(w, NULL)) 
                mcell_allocfailed("Failed to allocate grid for wall.");
          }else{
               continue;
          }
       }
       if(w->grid->n_tiles == 1)
       {
          push_tile_neighbor_to_list(&tile_nbr_head, w->grid, 0);
          tiles_count++;
       }else{
          /* find tile indices of the corner tiles */
          tile_idx_0 = 0;
          /* see function "move_strip_up()" */
          tile_idx_mid = w->grid->n_tiles - 2*(w->grid->n) + 1;
          tile_idx_last = w->grid->n_tiles - 1;

          /* find index of the neighbor tile on that wall */
          if(wall_share_vertex(w, vert_0)) 
          {
               push_tile_neighbor_to_list(&tile_nbr_head, w->grid, tile_idx_mid);
               tiles_count++;
          }else if(wall_share_vertex(w, vert_1)){
               push_tile_neighbor_to_list(&tile_nbr_head, w->grid, tile_idx_last);
               tiles_count++;
          }else if(wall_share_vertex(w, vert_2)){
               push_tile_neighbor_to_list(&tile_nbr_head, w->grid, tile_idx_0);
               tiles_count++;
          }
       }

   }
    
   *list_length = tiles_count;
   *tile_neighbor_head = tile_nbr_head;

}

/**************************************************************************
grid_all_neighbors_across_walls_through_edges: 
  In: a surface grid
      an index on that grid
      flag that tells whether we need to create a grid on a neighbor wall
      a linked list of  neighbor tiles (return value)
      a length of the linked list above (return value)
  Out: The list of nearest neighbors are returned,
       Neighbors should share either common edge or common vertice.
  Note: This version allows looking for the neighbors at the neighbor walls
        that are connected to the start wall through edges only.
****************************************************************************/
void grid_all_neighbors_across_walls_through_edges(struct surface_grid *grid, int idx, int create_grid_flag, struct tile_neighbor **tile_neighbor_head, int *list_length)
{
   struct tile_neighbor *tile_nbr_head = NULL; 
   int tiles_count = 0;
   int tiles_added = 0;  /* return value from the function 
                          "add_more_tile_neighbors_to_list()" */

   int kk; 
   struct surface_grid *sg[3];/* Neighboring surface grids (edge-to-edge) */
   int si[3];  /* Indices on those grids (edge-to-edge) of neighbor tiles */

   struct surface_grid *main_sg[3];/* Neighboring surface grids (edge-to-edge)                                        for the start tile  */
   int main_si[3];  /* Indices on those grids (edge-to-edge) of neighbor tiles                      for the start tile */

   int root, rootrem, strip, stripe, flip;
   int temp_idx;

   if((u_int)idx >= grid->n_tiles){
      mcell_internal_error("Grid molecule tile index is greater than or equal of the number of tiles on the grid\n");
   }
   
   /* find (strip, stripe, flip) coordinates of the tile */
   root  = (int)(sqrt((double) idx));
   rootrem = idx - root*root;
   strip = grid->n - root -1;
   stripe = rootrem/2;
   flip = rootrem - 2*stripe;

   for(kk = 0; kk < 3; kk++)
   {
     main_sg[kk] = NULL;
     main_si[kk] = -1;
   }

   /* find neighbors to react with */
   grid_neighbors(grid, idx, create_grid_flag, main_sg, main_si);

   if(stripe == 0)
   {
      if(flip > 0) /* inverted tile */
      {
         /* put in the list tiles that are on the same strip */
         push_tile_neighbor_to_list(&tile_nbr_head, grid, idx - 1);
         tiles_count++;
         push_tile_neighbor_to_list(&tile_nbr_head, grid, idx + 1);
         tiles_count++;
         if(strip < grid->n - 2)
         {
            push_tile_neighbor_to_list(&tile_nbr_head, grid, idx + 2);
            tiles_count++;
         }

         /* put in the list tiles that are on the row below the start tile
            but on the same grid */
         for(kk = 0; kk < 3; kk++)
         {
            if(main_sg[kk] == grid)
            {         
              if((main_si[kk] != idx -1) && (main_si[kk] != idx + 1))
              {
                push_tile_neighbor_to_list(&tile_nbr_head, grid, main_si[kk]);
                tiles_count++;
                if(strip < grid->n - 2)
                {
                   push_tile_neighbor_to_list(&tile_nbr_head, grid, main_si[kk] + 1);
                   tiles_count++;
                   push_tile_neighbor_to_list(&tile_nbr_head, grid, main_si[kk] + 2);
                   tiles_count++;
                }
                break;  
              }
            }
         }


         if(strip > 0)
         {
         /* put in the list tiles that are on the row above the start tile
            but on the same grid */
            temp_idx = move_strip_up(grid, idx);
            if(temp_idx == -1)
            {
               mcell_internal_error("Error in navigating on the grid");
            }
            push_tile_neighbor_to_list(&tile_nbr_head, grid, temp_idx);
            tiles_count++;
            push_tile_neighbor_to_list(&tile_nbr_head, grid, temp_idx - 1);
            tiles_count++;
            push_tile_neighbor_to_list(&tile_nbr_head, grid, temp_idx + 1);
            tiles_count++;
         }

         /* get the neighbors from the neighbor walls */
         for(kk = 0; kk < 3; kk++)
         {
            sg[kk] = NULL;
            si[kk] = -1;
         }
         grid_neighbors(grid, idx - 1, create_grid_flag, sg, si);
         for(kk = 0; kk < 3; kk++)
         {
            if((sg[kk] != grid) && (sg[kk] != NULL))
            {
               if(is_neighbor_tile(grid, idx, sg[kk], si[kk]))
               {
                  push_tile_neighbor_to_list(&tile_nbr_head, sg[kk], si[kk]);
                  tiles_count++;
               }
                    
               tiles_added = add_more_tile_neighbors_to_list(&tile_nbr_head, grid, idx, sg[kk], si[kk]); 
               tiles_count += tiles_added; 
                      
            }
         }
       

      }else{   /* upright tile (flip == 0)  */
        if(idx == 0) /* it is a special case */
        {
          for(kk = 0; kk < 3; kk++)
          {
            if((main_sg[kk] != NULL) && (main_sg[kk] == grid))
            {
               /* put in the list tiles that are on the row above the start tile */
               push_tile_neighbor_to_list(&tile_nbr_head, grid, main_si[kk]);
               tiles_count++;
               push_tile_neighbor_to_list(&tile_nbr_head, grid, main_si[kk] - 1);
               tiles_count++;
               push_tile_neighbor_to_list(&tile_nbr_head, grid, main_si[kk] + 1);
               tiles_count++;
            }else if((main_sg[kk] != NULL) && (main_sg[kk] != grid)){
               if(is_neighbor_tile(grid, idx, main_sg[kk], main_si[kk]))
               {
                  push_tile_neighbor_to_list(&tile_nbr_head, main_sg[kk], main_si[kk]);
                  tiles_count++;
               }
                    
               tiles_added = add_more_tile_neighbors_to_list(&tile_nbr_head, grid, idx, main_sg[kk], main_si[kk]); 
               tiles_count += tiles_added; 
                       
            }

          }
        }else{  /* if (idx != 0) */

            /* put in the list tiles that are on the same strip */
            push_tile_neighbor_to_list(&tile_nbr_head, grid, idx + 1);
            tiles_count++;
            push_tile_neighbor_to_list(&tile_nbr_head, grid, idx + 2);
            tiles_count++;
 
            /* put in the list tiles that are on the row below the start tile */
            temp_idx = move_strip_down(grid, idx + 1);
            if(temp_idx == -1)
            {
               mcell_internal_error("Error in navigating on the grid");
            }
            push_tile_neighbor_to_list(&tile_nbr_head, grid, temp_idx);
            tiles_count++;

            if(strip > 0)
            {

               /* put in the list tiles that are on the row above the start tile
                  but on the same grid */
               temp_idx = move_strip_up(grid, idx);
               if(temp_idx == -1)
               {
                  mcell_internal_error("Error in navigating on the grid");
               }
               push_tile_neighbor_to_list(&tile_nbr_head, grid, temp_idx);
               tiles_count++;
               push_tile_neighbor_to_list(&tile_nbr_head, grid, temp_idx - 1);
               tiles_count++;
               push_tile_neighbor_to_list(&tile_nbr_head, grid, temp_idx + 1);
               tiles_count++;
               push_tile_neighbor_to_list(&tile_nbr_head, grid, temp_idx + 2);
               tiles_count++;


            }else{         /* strip == 0 */
               /* put in the list tiles that are on the row above the start tile
               but on the different grid */
            
               for(kk = 0; kk < 3; kk++)
               {
                  if((main_sg[kk] != NULL) && (main_sg[kk] != grid))
                  {
                      push_tile_neighbor_to_list(&tile_nbr_head, main_sg[kk], main_si[kk]);
                      tiles_count++;
                            
                    tiles_added = add_more_tile_neighbors_to_list(&tile_nbr_head, grid, idx, main_sg[kk], main_si[kk]); 
                    tiles_count += tiles_added;  
                           
                  }
               }

            }

        } /* end if-else (idx == 0) */

      } /* end upright or inverted tile */

   } /* end if(stripe == 0) */



   if((strip == 0) && (stripe > 0))
   {
      /* put in the list tiles that are on the same row */ 
      push_tile_neighbor_to_list(&tile_nbr_head, grid, idx - 1);
      tiles_count++;
      push_tile_neighbor_to_list(&tile_nbr_head, grid, idx - 2);
      tiles_count++;

      if((stripe < grid->n - 2) || ((stripe == grid->n -2) && (flip == 0)))
      {
         push_tile_neighbor_to_list(&tile_nbr_head, grid, idx + 1);
         tiles_count++;
         push_tile_neighbor_to_list(&tile_nbr_head, grid, idx + 2);
         tiles_count++;
      }else if((stripe == grid->n - 2) && (flip == 1)){
         push_tile_neighbor_to_list(&tile_nbr_head, grid, idx + 1);
         tiles_count++;

         for(kk = 0; kk < 3; kk++)
         {
            sg[kk] = NULL;
            si[kk] = -1;
         }
         grid_neighbors(grid, idx + 1, create_grid_flag, sg, si);
         for(kk = 0; kk < 3; kk++)
         {
            if((sg[kk] != NULL) && (sg[kk] != grid))
            {
               if(is_neighbor_tile(grid, idx, sg[kk], si[kk]))
               {
                  push_tile_neighbor_to_list(&tile_nbr_head, sg[kk], si[kk]);
                  tiles_count++;
               }
               tiles_added = add_more_tile_neighbors_to_list(&tile_nbr_head, grid, idx, sg[kk], si[kk]); 
               tiles_count += tiles_added; 
            }
          }
      }else if(stripe == grid->n - 1){
         /* this is a corner tile */
         for(kk = 0; kk < 3; kk++)
         {
            sg[kk] = NULL;
            si[kk] = -1;
         }
         grid_neighbors(grid, idx, create_grid_flag, sg, si);
         for(kk = 0; kk < 3; kk++)
         {
            if((sg[kk] != NULL) && (sg[kk] != grid))
            {
               if(is_neighbor_tile(grid, idx, sg[kk], si[kk]))
               {
                  push_tile_neighbor_to_list(&tile_nbr_head, sg[kk], si[kk]);
                  tiles_count++;
               }
               tiles_added = add_more_tile_neighbors_to_list(&tile_nbr_head, grid, idx, sg[kk], si[kk]); 
               tiles_count += tiles_added; 
            }
          }

      }



      /* put in the list tiles that are on the row above */
      if(flip == 0)
      {
         for(kk = 0; kk < 3; kk++)
         {
            if((main_sg[kk] != NULL) && (main_sg[kk] != grid))
            {
               if(stripe < grid->n - 1)
               {
                  push_tile_neighbor_to_list(&tile_nbr_head, main_sg[kk], main_si[kk]);
                  tiles_count++;
                     
                  tiles_added = add_more_tile_neighbors_to_list(&tile_nbr_head, grid, idx, main_sg[kk], main_si[kk]); 
                  tiles_count += tiles_added; 
                  break;
               }      
            }         
         }
      }else{ /* flip > 0 */
         for(kk = 0; kk < 3; kk++)
         {
            sg[kk] = NULL;
            si[kk] = -1;
         }
         grid_neighbors(grid, idx - 1, create_grid_flag, sg, si);
         for(kk = 0; kk < 3; kk++)
         {
            if((sg[kk] != NULL) && (sg[kk] != grid))
            {
               if(is_neighbor_tile(grid, idx, sg[kk], si[kk]))
               {
                  push_tile_neighbor_to_list(&tile_nbr_head, sg[kk], si[kk]);
                  tiles_count++;
               }
                      
                  tiles_added = add_more_tile_neighbors_to_list(&tile_nbr_head, grid, idx, sg[kk], si[kk]); 
                  tiles_count += tiles_added; 
                  break; 
            }
         }
      }

      /* put in the list tiles that are on the row below */
      if(flip > 0)
      {
         temp_idx = move_strip_down(grid, idx);
         if(temp_idx == -1)
         {
             mcell_internal_error("Error in navigating on the grid");
         }
         push_tile_neighbor_to_list(&tile_nbr_head, grid, temp_idx);
         tiles_count++;
         push_tile_neighbor_to_list(&tile_nbr_head, grid, temp_idx - 1);
         tiles_count++;
         push_tile_neighbor_to_list(&tile_nbr_head, grid, temp_idx - 2);
         tiles_count++;
         if(stripe < grid->n - 2)
         {
            push_tile_neighbor_to_list(&tile_nbr_head, grid, temp_idx + 1);
            tiles_count++;
            push_tile_neighbor_to_list(&tile_nbr_head, grid, temp_idx + 2);
            tiles_count++;
         }
      }else{  /* (flip == 0) */
         if((unsigned int)idx < grid->n_tiles -1)
         {
            temp_idx = move_strip_down(grid, idx);
            if(temp_idx == -1)
            {
               mcell_internal_error("Error in navigating on the grid");
            }
            push_tile_neighbor_to_list(&tile_nbr_head, grid, temp_idx);
            tiles_count++;
            push_tile_neighbor_to_list(&tile_nbr_head, grid, temp_idx - 1);
            tiles_count++;
            push_tile_neighbor_to_list(&tile_nbr_head, grid, temp_idx + 1);
            tiles_count++;
         }else{
             /* this is a corner tile */
            temp_idx = move_strip_down(grid, idx -1);
            if(temp_idx == -1)
            {
               mcell_internal_error("Error in navigating on the grid");
            }
            push_tile_neighbor_to_list(&tile_nbr_head, grid, temp_idx);
            tiles_count++;
         }
      }


   } /* end if ((strip == 0) && (stripe > 0)) */



   if((strip > 0) && (stripe > 0))
   {
      /* We are guaranteed to be in the right layer of the wall.
         Note that no extra checking is done here for that
         assumption besides calling the function "is_inner_tile()"
         before calling this function. */
      if(flip > 0)
      {
         /* put in the list tiles that are on the same row */
            push_tile_neighbor_to_list(&tile_nbr_head, grid, idx - 1);
            tiles_count++;
            push_tile_neighbor_to_list(&tile_nbr_head, grid, idx - 2);
            tiles_count++;
            push_tile_neighbor_to_list(&tile_nbr_head, grid, idx + 1);
            tiles_count++;
            for(kk = 0; kk < 3; kk++)
            {
               sg[kk] = NULL;
               si[kk] = -1;
            }
            grid_neighbors(grid, idx + 1, create_grid_flag, sg, si);
            for(kk = 0; kk < 3; kk++)
            {
               if((sg[kk] != NULL) && (sg[kk] != grid))
               { 
                  if(is_neighbor_tile(grid, idx, sg[kk], si[kk]))
                  {
                     push_tile_neighbor_to_list(&tile_nbr_head, sg[kk], si[kk]);
                     tiles_count++;
                  }
                         
                  tiles_added = add_more_tile_neighbors_to_list(&tile_nbr_head, grid, idx, sg[kk], si[kk]); 
                  tiles_count += tiles_added; 
                  break;      
               }
            }

         /* put in the list tiles that are above the current row */
            temp_idx = move_strip_up(grid, idx);
            if(temp_idx == -1)
            {
               mcell_internal_error("Error in navigating on the grid");
            }
            push_tile_neighbor_to_list(&tile_nbr_head, grid, temp_idx);
            tiles_count++;
            push_tile_neighbor_to_list(&tile_nbr_head, grid, temp_idx - 1);
            tiles_count++;
            push_tile_neighbor_to_list(&tile_nbr_head, grid, temp_idx + 1);
            tiles_count++;
         /* put in the list tiles that are below the current row */
            temp_idx = move_strip_down(grid, idx);
            if(temp_idx == -1)
            {
               mcell_internal_error("Error in navigating on the grid");
            }
            push_tile_neighbor_to_list(&tile_nbr_head, grid, temp_idx);
            tiles_count++;
            push_tile_neighbor_to_list(&tile_nbr_head, grid, temp_idx - 1);
            tiles_count++;
            push_tile_neighbor_to_list(&tile_nbr_head, grid, temp_idx - 2);
            tiles_count++;

      }else{   /* (flip == 0) */

         /* put in the list tiles that are on the same row */
            push_tile_neighbor_to_list(&tile_nbr_head, grid, idx - 1);
            tiles_count++;
            push_tile_neighbor_to_list(&tile_nbr_head, grid, idx - 2);
            tiles_count++;
            for(kk = 0; kk < 3; kk++)
            {
               if((main_sg[kk] != NULL) && (main_sg[kk] != grid))
               {
                   push_tile_neighbor_to_list(&tile_nbr_head, main_sg[kk], main_si[kk]);
                   tiles_count++;
                        
                  tiles_added = add_more_tile_neighbors_to_list(&tile_nbr_head, grid, idx, main_sg[kk], main_si[kk]); 
                  tiles_count += tiles_added;  
                  break;    
               }
            }

         /* put in the list tiles that are above the current row */
            temp_idx = move_strip_up(grid, idx);
            if(temp_idx == -1)
            {
               mcell_internal_error("Error in navigating on the grid");
            }
            push_tile_neighbor_to_list(&tile_nbr_head, grid, temp_idx);
            tiles_count++;
            push_tile_neighbor_to_list(&tile_nbr_head, grid, temp_idx - 1);
            tiles_count++;
            push_tile_neighbor_to_list(&tile_nbr_head, grid, temp_idx - 2);
            tiles_count++;
            push_tile_neighbor_to_list(&tile_nbr_head, grid, temp_idx + 1);
            tiles_count++;
         /* put in the list tiles that are below the current row */
            temp_idx = move_strip_down(grid, idx - 1);
            if(temp_idx == -1)
            {
               mcell_internal_error("Error in navigating on the grid");
            }
            push_tile_neighbor_to_list(&tile_nbr_head, grid, temp_idx);
            tiles_count++;
      }

   } /* end if((strip > 0) && (stripe > 0)) */


    *list_length = tiles_count;
    *tile_neighbor_head = tile_nbr_head;

}


/***************************************************************************
add_more_tile_neighbors_to_list: 
  In: a linked list of  neighbor tiles
      a surface grid of the start tile
      an index on that grid of the start tile
      a surface grid of the neighbor tile
      an index on that neighbor grid of the neighbor tile
  Out: Returns number of new neighbor tiles added to the original
       linked list of neighbor tiles.
       Neighbors should share either common edge or common vertice.
  Note: This function finds neighbors of the tile (orig_grid/orig_idx)
        using first found neighbor (new_grid/new_idx).  The tile
        (new_grid/new_idx) is already present in "tile_nbr_head".
        The function looks only for the tiles that are the neighbors of
        "orig_grid/orig_idx" and reside on the neighbor "new_grid".
****************************************************************************/
int add_more_tile_neighbors_to_list(struct tile_neighbor **tile_nbr_head, struct surface_grid *orig_grid, int orig_idx, struct surface_grid *new_grid, int new_idx)
{
    if(orig_grid == new_grid) {
       mcell_internal_error("Function 'add_more_tile_neighbors_to_list()' should be called for different grids only");
    }
    
    if((unsigned int)new_idx > new_grid->n_tiles - 1){
       mcell_internal_error("Function 'add_more_tile_neighbors_to_list()' was called with tile index greater than number of tiles");
    }  
    
    int root, rootrem, strip, stripe, flip;
    int tiles_added = 0;
    int check_1, check_2;    /* flags */ 
    int temp_idx = -1;

   /* find (strip, stripe, flip) coordinates of the tile with "new_idx"  */
   root  = (int)(sqrt((double) new_idx));
   rootrem = new_idx - root*root;
   strip = new_grid->n - root -1;
   stripe = rootrem/2;
   flip = rootrem - 2*stripe;

   if(stripe == 0)
   {
     if(flip > 0)
     {
        if(is_neighbor_tile(orig_grid, orig_idx, new_grid, new_idx - 1))
        {
           push_tile_neighbor_to_list(tile_nbr_head, new_grid, new_idx - 1);
           tiles_added++;
        }

        /* go up */
        check_1 = 1;
        check_2 = 1;
        temp_idx = new_idx -1;
        while(check_1 || check_2)
        {
           check_1 = 0;
           check_2 = 0;
           temp_idx = move_strip_up(new_grid, temp_idx);
           if(temp_idx < 0) break; /* we have moved from the new_grid */

           check_1 = is_neighbor_tile(orig_grid, orig_idx, new_grid, temp_idx);
           if(temp_idx >= 1)
           {
             check_2 = is_neighbor_tile(orig_grid, orig_idx, new_grid, temp_idx -1);
           }
           if(check_1)
           {
              push_tile_neighbor_to_list(tile_nbr_head, new_grid, temp_idx);
              tiles_added++;
           }
           if(check_2)
           {
              push_tile_neighbor_to_list(tile_nbr_head, new_grid, temp_idx -1);
              tiles_added++;
           }
           temp_idx = temp_idx -1;
        }
        /* go down */
        temp_idx = new_idx;
        check_1 = 1;
        check_2 = 1;
        while(check_1 || check_2)
        {
           check_1 = 0;
           check_2 = 0;
           temp_idx = move_strip_down(new_grid, temp_idx);
           if(temp_idx < 0) break; /* we have moved from the new_grid */
           if((unsigned int)temp_idx > new_grid->n_tiles - 1) break;

           check_1 = is_neighbor_tile(orig_grid, orig_idx, new_grid, temp_idx);
           if((unsigned int)temp_idx < new_grid->n_tiles - 1)
           {
              check_2 = is_neighbor_tile(orig_grid, orig_idx, new_grid, temp_idx + 1);
           }
           if(check_1)
           {
              push_tile_neighbor_to_list(tile_nbr_head, new_grid, temp_idx);
              tiles_added++;
           }
           if(check_2)
           {
              push_tile_neighbor_to_list(tile_nbr_head, new_grid, temp_idx + 1);
              tiles_added++;
           }
           temp_idx = temp_idx + 1;
        }
     }else{
       /* flip == 0 */
        if(is_neighbor_tile(orig_grid, orig_idx, new_grid, new_idx + 1))
        {
           push_tile_neighbor_to_list(tile_nbr_head, new_grid, new_idx + 1);
           tiles_added++;
        }
        /* go up */
        temp_idx = new_idx;
        check_1 = 1;
        check_2 = 1;
        while(check_1 || check_2)
        {
           check_1 = 0;
           check_2 = 0;
           temp_idx = move_strip_up(new_grid, temp_idx);
           if(temp_idx < 0) break; /* we have moved from the new_grid */

           check_1 = is_neighbor_tile(orig_grid, orig_idx, new_grid, temp_idx);
           if(temp_idx >= 1)
           {
              check_2 = is_neighbor_tile(orig_grid, orig_idx, new_grid, temp_idx -1);
           }
           if(check_1)
           {
              push_tile_neighbor_to_list(tile_nbr_head, new_grid, temp_idx);
              tiles_added++;
           }
           if(check_2)
           {
              push_tile_neighbor_to_list(tile_nbr_head, new_grid, temp_idx -1);
              tiles_added++;
           }
           temp_idx = temp_idx -1;
        }
        /* go down */
        temp_idx = new_idx + 1;
        check_1 = 1;
        check_2 = 1;
        while(check_1 || check_2)
        {
           check_1 = 0;
           check_2 = 0;
           temp_idx = move_strip_down(new_grid, temp_idx);
           if(temp_idx < 0) break; /* we have moved from the new_grid */
           if((unsigned int)temp_idx > new_grid->n_tiles - 1) break;

           check_1 = is_neighbor_tile(orig_grid, orig_idx, new_grid, temp_idx);
           if((unsigned int)temp_idx < new_grid->n_tiles - 1)
           {
              check_2 = is_neighbor_tile(orig_grid, orig_idx, new_grid, temp_idx + 1);
           }
           if(check_1)
           {
              push_tile_neighbor_to_list(tile_nbr_head, new_grid, temp_idx);
              tiles_added++;
           }
           if(check_2)
           {
              push_tile_neighbor_to_list(tile_nbr_head, new_grid, temp_idx + 1);
              tiles_added++;
           }
           temp_idx = temp_idx + 1;
        }

     }

   } /* end if(stripe == 0) */

   if((strip == 0) && (stripe > 0))
   {
        /* go up */
        temp_idx = new_idx + 1;
        while(1)
        {
            if((unsigned int)temp_idx > new_grid->n_tiles -1) break;
           if(is_neighbor_tile(orig_grid, orig_idx, new_grid, temp_idx))
           {
              push_tile_neighbor_to_list(tile_nbr_head, new_grid, temp_idx);
              tiles_added++;
              temp_idx = temp_idx + 1;
           }else break;
        }
        /* go down */
        temp_idx = new_idx - 1;
        while(1)
        {
           if(temp_idx < 0) break;
           if(is_neighbor_tile(orig_grid, orig_idx, new_grid, temp_idx))
           {
              push_tile_neighbor_to_list(tile_nbr_head, new_grid, temp_idx);
              tiles_added++;
              temp_idx = temp_idx - 1;
           }else break;
        }
   }        

   if((strip > 0) && (stripe > 0))
   {
      if(flip > 0)
      {
        if(is_neighbor_tile(orig_grid, orig_idx, new_grid, new_idx + 1))
        {
           push_tile_neighbor_to_list(tile_nbr_head, new_grid, new_idx + 1);
           tiles_added++;
        }
        /* go up */
        check_1 = 1;
        check_2 = 1;
        temp_idx = new_idx + 1;
        while(check_1 || check_2)
        {
           check_1 = 0;
           check_2 = 0;
           temp_idx = move_strip_up(new_grid, temp_idx);
           if(temp_idx < 0) break; /* we have moved from the new_grid */
           if((unsigned int)temp_idx > new_grid->n_tiles -1) break;

           check_1 = is_neighbor_tile(orig_grid, orig_idx, new_grid, temp_idx);
           if((unsigned int)temp_idx < new_grid->n_tiles - 1)
           {
              check_2 = is_neighbor_tile(orig_grid, orig_idx, new_grid, temp_idx + 1);
           }
           if(check_1)
           {
              push_tile_neighbor_to_list(tile_nbr_head, new_grid, temp_idx);
              tiles_added++;
           }
           if(check_2)
           {
              push_tile_neighbor_to_list(tile_nbr_head, new_grid, temp_idx + 1);
              tiles_added++;
           }
           temp_idx = temp_idx + 1;
        }
        /* go down */
        temp_idx = new_idx;
        check_1 = 1;
        check_2 = 1;
        while(check_1 || check_2)
        {
           check_1 = 0;
           check_2 = 0;
           temp_idx = move_strip_down(new_grid, temp_idx);
           if(temp_idx < 0) break; /* we have moved from the new_grid */

           check_1 = is_neighbor_tile(orig_grid, orig_idx, new_grid, temp_idx);
           if(temp_idx >= 1)
           {
              check_2 = is_neighbor_tile(orig_grid, orig_idx, new_grid, temp_idx - 1);
           }
           if(check_1)
           {
              push_tile_neighbor_to_list(tile_nbr_head, new_grid, temp_idx);
              tiles_added++;
           }
           if(check_2)
           {
              push_tile_neighbor_to_list(tile_nbr_head, new_grid, temp_idx - 1);
              tiles_added++;
           }
           temp_idx = temp_idx - 1;
        }

      }else{
         /* flip == 0 */
        if(is_neighbor_tile(orig_grid, orig_idx, new_grid, new_idx - 1))
        {
           push_tile_neighbor_to_list(tile_nbr_head, new_grid, new_idx - 1);
           tiles_added++;
        }
        /* go up */
        temp_idx = new_idx;
        check_1 = 1;
        check_2 = 1;
        while(check_1 || check_2)
        {
           check_1 = 0;
           check_2 = 0;
           temp_idx = move_strip_up(new_grid, temp_idx);
           if(temp_idx < 0) break; /* we have moved from the new_grid */
           if((unsigned int)temp_idx > new_grid->n_tiles - 1) break;

           check_1 = is_neighbor_tile(orig_grid, orig_idx, new_grid, temp_idx);
           if((unsigned int)temp_idx < new_grid->n_tiles - 1)
           {
              check_2 = is_neighbor_tile(orig_grid, orig_idx, new_grid, temp_idx + 1);
           }
           if(check_1)
           {
              push_tile_neighbor_to_list(tile_nbr_head, new_grid, temp_idx);
              tiles_added++;
           }
           if(check_2)
           {
              push_tile_neighbor_to_list(tile_nbr_head, new_grid, temp_idx + 1);
              tiles_added++;
           }
           temp_idx = temp_idx + 1;
        }
        /* go down */
        temp_idx = new_idx - 1;
        check_1 = 1;
        check_2 = 1;
        while(check_1 || check_2)
        {
           check_1 = 0;
           check_2 = 0;
           temp_idx = move_strip_down(new_grid, temp_idx);
           if(temp_idx < 0) break; /* we have moved from the new_grid */

           check_1 = is_neighbor_tile(orig_grid, orig_idx, new_grid, temp_idx);
           if(temp_idx >= 1)
           { 
             check_2 = is_neighbor_tile(orig_grid, orig_idx, new_grid, temp_idx - 1);
           }
           if(check_1)
           {
              push_tile_neighbor_to_list(tile_nbr_head, new_grid, temp_idx);
              tiles_added++;
           }
           if(check_2)
           {
              push_tile_neighbor_to_list(tile_nbr_head, new_grid, temp_idx - 1);
              tiles_added++;
           }
           temp_idx = temp_idx - 1;
        }

      }

   } /* end if((strip > 0) && (stripe > 0)) */

   return tiles_added;
   
}

/***************************************************************************
is_neighbor_tile: 
  In: a surface grid of the start tile
      an index on that grid of the start tile
      a surface grid of the neighbor tile
      an index on that neighbor grid of the neighbor tile
  Out: Returns 1 if the tiles are neighbors, and 0 otherwise.
       Neighbors should share either common edge or common vertice.
****************************************************************************/
int is_neighbor_tile(struct surface_grid *orig_grid, int orig_idx, struct surface_grid *new_grid, int new_idx)
{
   /* corners of the original tile in surface and world coordinates */
   struct vector2 A, B, C;
   struct vector3 A_3d, B_3d, C_3d;
   /* corners of the neighbor tile in surface and world coordinates */
   struct vector2 R, S, T;
   struct vector3 R_3d, S_3d, T_3d;
   int flip_orig, flip_new;  /* flip coordinate of the tile */

   if(new_grid == NULL)
   {
      mcell_internal_error("Function 'is_neighbor_tile()' called for NULL surface_grid");
   }
   
   if(new_idx < 0) return 0;
   if((u_int)new_idx > new_grid->n_tiles - 1) return 0;
   
   get_tile_vertices(orig_grid, orig_idx, &flip_orig, &A, &B, &C);
   uv2xyz(&A, orig_grid->surface, &A_3d);
   uv2xyz(&B, orig_grid->surface, &B_3d);
   uv2xyz(&C, orig_grid->surface, &C_3d);

   get_tile_vertices(new_grid, new_idx, &flip_new, &R, &S, &T); 
   uv2xyz(&R, new_grid->surface, &R_3d);
   uv2xyz(&S, new_grid->surface, &S_3d);
   uv2xyz(&T, new_grid->surface, &T_3d);

   if(point_in_triangle(&R_3d, &A_3d, &B_3d, &C_3d)) {
         return 1; 
   }
   if(point_in_triangle(&S_3d, &A_3d, &B_3d, &C_3d)){
         return 1; 
   }
   if(point_in_triangle(&T_3d, &A_3d, &B_3d, &C_3d)){
        return 1;  
   }


   return 0;
}

/************************************************************************
find_closest_position:
  In: surface grid of the first tile
      first tile index
      surface grid of the second tile
      second (neighbor) tile index
  Out: position of the product on the first tile that is closest
       to the second tile. If the neighbor tiles have common edge
       this position happens to be very close to the center of the 
       common edge but inward the first tile.  If the neighbor tiles
       have common vertex, this position happens to be very close to
       to the vertex but inward the first tile.                               
*************************************************************************/
void find_closest_position(struct surface_grid *grid1, int idx1, struct surface_grid *grid2, int idx2, struct vector2 *p)
{

  /* vertices of the first tile */
   struct vector2 R, S, T;
   struct vector3 R_3d, S_3d, T_3d;
  
   /* vertices of the second tile */
   struct vector2 A, B, C;
   struct vector3 A_3d, B_3d, C_3d;
   /* vertices A,B,C in the coordinate system RST */
   struct vector2 A_new, B_new, C_new;
   
   /* the ratios in which we divide the segment */
   double k1 = 1e-10;  /* this is our good faith assumption */
   double k2 = 1;

   int flip1; /* flip information about first tile */
   int flip2; /* flip information about second tile */
   
   int num_shared_vertices = 0, num_exact_shared_vertices = 0;
   /* flags */
   int R_shared = 0, S_shared = 0, T_shared = 0;
   int A_shared = 0, B_shared = 0, C_shared = 0;

   /* find out the vertices of the first tile where we will put the product */
   get_tile_vertices(grid1, idx1, &flip1, &R, &S, &T);
   uv2xyz(&R, grid1->surface, &R_3d);
   uv2xyz(&S, grid1->surface, &S_3d);
   uv2xyz(&T, grid1->surface, &T_3d);

   /* find out the vertices of the second tile  */
   get_tile_vertices(grid2, idx2, &flip2, &A, &B, &C);
   uv2xyz(&A, grid2->surface, &A_3d);
   uv2xyz(&B, grid2->surface, &B_3d);
   uv2xyz(&C, grid2->surface, &C_3d);

    /* find shared vertices */
    if(grid1 == grid2)
    {
       if(!distinguishable_vec2(&R, &A, EPS_C) || (!distinguishable_vec2(&R, &B, EPS_C)) || (!distinguishable_vec2(&R, &C, EPS_C)))
       {
       
           num_exact_shared_vertices++;
           R_shared = 1;
       }
       if(!distinguishable_vec2(&S, &A, EPS_C) || (!distinguishable_vec2(&S, &B, EPS_C)) || (!distinguishable_vec2(&S, &C, EPS_C)))
       {
           num_exact_shared_vertices++;
           S_shared = 1;
       }
       if(!distinguishable_vec2(&T, &A, EPS_C) || (!distinguishable_vec2(&T, &B, EPS_C)) || (!distinguishable_vec2(&T, &C, EPS_C)))
       {  
           num_exact_shared_vertices++;
           T_shared = 1;
       }

    }else{
        /* below there are cases when the grid structures on the neighbor 
           walls are not shifted relative to one another */
     if(!distinguishable_vec3(&R_3d, &A_3d, EPS_C) || (!distinguishable_vec3(&R_3d, &B_3d, EPS_C)) || (!distinguishable_vec3(&R_3d, &C_3d, EPS_C)))
     {
              num_exact_shared_vertices++;
              R_shared = 1;
      }

     if(!distinguishable_vec3(&S_3d, &A_3d, EPS_C) || (!distinguishable_vec3(&S_3d, &B_3d, EPS_C)) || (!distinguishable_vec3(&S_3d, &C_3d, EPS_C)))
      {
             num_exact_shared_vertices++;
             S_shared = 1;
      }

     if(!distinguishable_vec3(&T_3d, &A_3d, EPS_C) || (!distinguishable_vec3(&T_3d, &B_3d, EPS_C)) || (!distinguishable_vec3(&T_3d, &C_3d, EPS_C)))
     {
             num_exact_shared_vertices++;
             T_shared = 1;
      }
    }
    

   if(num_exact_shared_vertices == 1)
   {
      if(R_shared)
      {
         place_product_shared_vertex(&R, &S, &T, p, k1, k2);
         return;
      }else if(S_shared){
         place_product_shared_vertex(&S, &R, &T, p, k1, k2);
         return;
      }else{ /*T is shared */
         place_product_shared_vertex(&T, &R, &S, p, k1, k2);
         return;
      }
   }

   if(num_exact_shared_vertices == 2)
   {
      if(R_shared && S_shared)
      {
         place_product_shared_segment(&R, &S, &T, p, k1, k2); 
         return;
      }else if(R_shared && T_shared){
          place_product_shared_segment(&R, &T, &S, p, k1, k2);
          return;
      }else{  /*S_shared and T_shared */
          place_product_shared_segment(&S, &T, &R, p, k1, k2);
          return;
      }
   }   

    if(num_exact_shared_vertices == 0)
    {
        /* below are the cases when the grids on the neighbor walls
           are shifted relative to one another */
       /* find out whether the vertices of one tile cross the sides of
          another tile */

       
       if((intersect_point_segment(&S_3d, &A_3d, &B_3d))
         ||  (intersect_point_segment(&S_3d, &B_3d, &C_3d))
         ||  (intersect_point_segment(&S_3d, &A_3d, &C_3d)))
       {
          S_shared = 1;
          num_shared_vertices++;
       }

       if((intersect_point_segment(&R_3d, &A_3d, &B_3d))
          || (intersect_point_segment(&R_3d, &B_3d, &C_3d))
          || (intersect_point_segment(&R_3d, &A_3d, &C_3d)))
       {
          R_shared = 1;
          num_shared_vertices++;
       }

       if((intersect_point_segment(&T_3d, &A_3d, &B_3d))
         || (intersect_point_segment(&T_3d, &B_3d, &C_3d))
         || (intersect_point_segment(&T_3d, &A_3d, &C_3d)))
       {
           T_shared = 1;
           num_shared_vertices++;
       }
       
       if((intersect_point_segment(&A_3d, &R_3d, &S_3d))
        ||  (intersect_point_segment(&A_3d, &S_3d, &T_3d))
        ||  (intersect_point_segment(&A_3d, &R_3d, &T_3d)))
       {
           A_shared = 1;
           num_shared_vertices++;
       }

       if((intersect_point_segment(&B_3d, &R_3d, &S_3d))
        ||  (intersect_point_segment(&B_3d, &S_3d, &T_3d))
        ||  (intersect_point_segment(&B_3d, &R_3d, &T_3d)))
       {
           B_shared = 1;
           num_shared_vertices++;
       }
       
       if((intersect_point_segment(&C_3d, &R_3d, &S_3d))
        ||  (intersect_point_segment(&C_3d, &S_3d, &T_3d))
        ||  (intersect_point_segment(&C_3d, &R_3d, &T_3d)))
       {
           C_shared = 1;
           num_shared_vertices++;
       }
 
        /* two vertices shared from the same tile */
       if(R_shared && S_shared)
       {
         place_product_shared_segment(&R, &S, &T, p, k1, k2); 
         return;
       }else if(R_shared && T_shared){
         place_product_shared_segment(&R, &T, &S, p, k1, k2); 
         return;
       }else if(S_shared && T_shared){
         place_product_shared_segment(&S, &T, &R, p, k1, k2);
         return; 
       }

      
        /* two vertices shared from the same tile */
       if(A_shared && B_shared)
       {
         if(parallel_segments(&A_3d, &B_3d, &R_3d, &S_3d))
         {
            xyz2uv(&A_3d, grid1->surface, &A_new);
            xyz2uv(&B_3d, grid1->surface, &B_new);
            place_product_shared_segment(&A_new, &B_new, &T, p, k1, k2);
            return; 
         }else if(parallel_segments(&A_3d, &B_3d, &R_3d, &T_3d)){
            xyz2uv(&A_3d, grid1->surface, &A_new);
            xyz2uv(&B_3d, grid1->surface, &B_new);
            place_product_shared_segment(&A_new, &B_new, &S, p, k1, k2);
            return; 
         }else if(parallel_segments(&A_3d, &B_3d, &S_3d, &T_3d)){
            xyz2uv(&A_3d, grid1->surface, &A_new);
            xyz2uv(&B_3d, grid1->surface, &B_new);
            place_product_shared_segment(&A_new, &B_new, &R, p, k1, k2);
            return; 
         }

       }else if(A_shared && C_shared){
         if(parallel_segments(&A_3d, &C_3d, &R_3d, &S_3d))
         {
            xyz2uv(&A_3d, grid1->surface, &A_new);
            xyz2uv(&C_3d, grid1->surface, &C_new);
            place_product_shared_segment(&A_new, &C_new, &T, p, k1, k2);
            return; 
         }else if(parallel_segments(&A_3d, &C_3d, &R_3d, &T_3d)){
            xyz2uv(&A_3d, grid1->surface, &A_new);
            xyz2uv(&C_3d, grid1->surface, &C_new);
            place_product_shared_segment(&A_new, &C_new, &S, p, k1, k2);
            return; 
         }else if(parallel_segments(&A_3d, &C_3d, &S_3d, &T_3d)){
            xyz2uv(&A_3d, grid1->surface, &A_new);
            xyz2uv(&C_3d, grid1->surface, &C_new);
            place_product_shared_segment(&A_new, &C_new, &R, p, k1, k2);
            return; 
         }

       }else if(B_shared && C_shared){
         if(parallel_segments(&B_3d, &C_3d, &R_3d, &S_3d))
         {
            xyz2uv(&B_3d, grid1->surface, &B_new);
            xyz2uv(&C_3d, grid1->surface, &C_new);
            place_product_shared_segment(&B_new, &C_new, &T, p, k1, k2);
            return; 
         }else if(parallel_segments(&B_3d, &C_3d, &R_3d, &T_3d)){
            xyz2uv(&B_3d, grid1->surface, &B_new);
            xyz2uv(&C_3d, grid1->surface, &C_new);
            place_product_shared_segment(&B_new, &C_new, &S, p, k1, k2);
            return; 
         }else if(parallel_segments(&B_3d, &C_3d, &S_3d, &T_3d)){
            xyz2uv(&B_3d, grid1->surface, &B_new);
            xyz2uv(&C_3d, grid1->surface, &C_new);
            place_product_shared_segment(&B_new, &C_new, &R, p, k1, k2);
            return; 
         }

       }

       /* one vertex shared from each tile */
       if (R_shared)
       {
           if(A_shared)
           {
              if(parallel_segments(&R_3d, &A_3d, &R_3d, &S_3d))
              {
                 xyz2uv(&A_3d, grid1->surface, &A_new);
                 place_product_shared_segment(&A_new, &R, &T, p, k1, k2);
                 return; 
              }else if(parallel_segments(&R_3d, &A_3d, &R_3d, &T_3d)){
                 xyz2uv(&A_3d, grid1->surface, &A_new);
                 place_product_shared_segment(&A_new, &R, &S, p, k1, k2);
                 return; 
              }

           }else if(B_shared){
              if(parallel_segments(&R_3d, &B_3d, &R_3d, &S_3d))
              {
                 xyz2uv(&B_3d, grid1->surface, &B_new);
                 place_product_shared_segment(&B_new, &R, &T, p, k1, k2);
                 return; 
              }else if(parallel_segments(&R_3d, &B_3d, &R_3d, &T_3d)){
                 xyz2uv(&B_3d, grid1->surface, &B_new);
                 place_product_shared_segment(&B_new, &R, &S, p, k1, k2);
                 return; 
              }

           }else if(C_shared){
              if(parallel_segments(&R_3d, &C_3d, &R_3d, &S_3d))
              {
                 xyz2uv(&C_3d, grid1->surface, &C_new);
                 place_product_shared_segment(&C_new, &R, &T, p, k1, k2);
                 return; 
              }else if(parallel_segments(&R_3d, &C_3d, &R_3d, &T_3d)){
                 xyz2uv(&C_3d, grid1->surface, &C_new);
                 place_product_shared_segment(&C_new, &R, &S, p, k1, k2);
                 return; 
              }

           }else{
              place_product_shared_vertex(&R, &S, &T, p, k1, k2);
              return;  
           }

       }else if(S_shared){
           if(A_shared)
           {
              if(parallel_segments(&S_3d, &A_3d, &S_3d, &T_3d))
              {
                 xyz2uv(&A_3d, grid1->surface, &A_new);
                 place_product_shared_segment(&A_new, &S, &R, p, k1, k2);
                 return; 
              }else if(parallel_segments(&S_3d, &A_3d, &S_3d, &R_3d)){
                 xyz2uv(&A_3d, grid1->surface, &A_new);
                 place_product_shared_segment(&A_new, &S, &T, p, k1, k2);
                 return; 
              }

           }else if(B_shared){
              if(parallel_segments(&S_3d, &B_3d, &S_3d, &T_3d))
              {
                 xyz2uv(&B_3d, grid1->surface, &B_new);
                 place_product_shared_segment(&B_new, &S, &R, p, k1, k2);
                 return; 
              }else if(parallel_segments(&S_3d, &B_3d, &S_3d, &R_3d)){
                 xyz2uv(&B_3d, grid1->surface, &B_new);
                 place_product_shared_segment(&B_new, &S, &T, p, k1, k2);
                 return; 
              }

           }else if(C_shared){
              if(parallel_segments(&S_3d, &C_3d, &S_3d, &T_3d))
              {
                 xyz2uv(&C_3d, grid1->surface, &C_new);
                 place_product_shared_segment(&C_new, &S, &R, p, k1, k2);
                 return; 
              }else if(parallel_segments(&S_3d, &C_3d, &S_3d, &R_3d)){
                 xyz2uv(&C_3d, grid1->surface, &C_new);
                 place_product_shared_segment(&C_new, &S, &T, p, k1, k2);
                 return; 
              }

           }else{
              place_product_shared_vertex(&S, &R, &T, p, k1, k2);
              return;  
           }

       }else if(T_shared){
           if(A_shared)
           {
              if(parallel_segments(&T_3d, &A_3d, &T_3d, &S_3d))
              {
                 xyz2uv(&A_3d, grid1->surface, &A_new);
                 place_product_shared_segment(&A_new, &T, &R, p, k1, k2);
                 return; 
              }else if(parallel_segments(&T_3d, &A_3d, &T_3d, &R_3d)){
                 xyz2uv(&A_3d, grid1->surface, &A_new);
                 place_product_shared_segment(&A_new, &T, &S, p, k1, k2);
                 return;
              } 
           }else if(B_shared){
              if(parallel_segments(&T_3d, &B_3d, &T_3d, &S_3d))
              {
                 xyz2uv(&B_3d, grid1->surface, &B_new);
                 place_product_shared_segment(&B_new, &T, &R, p, k1, k2);
                 return; 
              }else if(parallel_segments(&T_3d, &B_3d, &T_3d, &R_3d)){
                 xyz2uv(&B_3d, grid1->surface, &B_new);
                 place_product_shared_segment(&B_new, &T, &S, p, k1, k2);
                 return;
              } 

           }else if(C_shared){
              if(parallel_segments(&T_3d, &C_3d, &T_3d, &S_3d))
              {
                 xyz2uv(&C_3d, grid1->surface, &C_new);
                 place_product_shared_segment(&C_new, &T, &R, p, k1, k2);
                 return; 
              }else if(parallel_segments(&T_3d, &C_3d, &T_3d, &R_3d)){
                 xyz2uv(&C_3d, grid1->surface, &C_new);
                 place_product_shared_segment(&C_new, &T, &S, p, k1, k2);
                 return;
              } 

           }else{
              place_product_shared_vertex(&T, &R, &S, p, k1, k2);
              return;  
           }

       }
         
       /* only one vertex is shared */
       if(A_shared)
       {
          xyz2uv(&A_3d, grid1->surface, &A_new); 
          if(intersect_point_segment(&A_3d, &R_3d, &S_3d))
          {
             place_product_close_to_segment_endpoint(&T, &A_new, p, k1, k2); 
             return;
          }else if(intersect_point_segment(&A_3d, &R_3d, &T_3d)){
             place_product_close_to_segment_endpoint(&S, &A_new, p, k1, k2); 
             return;
          }else if(intersect_point_segment(&A_3d, &S_3d, &T_3d)){
             place_product_close_to_segment_endpoint(&R, &A_new, p, k1, k2); 
             return;
          }

       }else if(B_shared){
          xyz2uv(&B_3d, grid1->surface, &B_new); 
          if(intersect_point_segment(&B_3d, &R_3d, &S_3d))
          {
             place_product_close_to_segment_endpoint(&T, &B_new, p, k1, k2); 
             return;
          }else if(intersect_point_segment(&B_3d, &R_3d, &T_3d)){
             place_product_close_to_segment_endpoint(&S, &B_new, p, k1, k2); 
             return;
          }else if(intersect_point_segment(&B_3d, &S_3d, &T_3d)){
             place_product_close_to_segment_endpoint(&R, &B_new, p, k1, k2); 
             return;
          }
       }else if(C_shared){
          xyz2uv(&C_3d, grid1->surface, &C_new); 
          if(intersect_point_segment(&C_3d, &R_3d, &S_3d))
          {
             place_product_close_to_segment_endpoint(&T, &C_new, p, k1, k2); 
             return;
          }else if(intersect_point_segment(&C_3d, &R_3d, &T_3d)){
             place_product_close_to_segment_endpoint(&S, &C_new, p, k1, k2); 
             return;
          }else if(intersect_point_segment(&C_3d, &S_3d, &T_3d)){
             place_product_close_to_segment_endpoint(&R, &C_new, p, k1, k2); 
             return;
          }
       }

    } /* end if (num_shared_vertices == 0) */

       /* No vertices are shared but the tiles are still neighbors.
          It may mean that there is a roundoff error in calculating
          tile vertices */

    if((num_exact_shared_vertices == 0) && (num_shared_vertices == 0)) 
    {
      /* find points on the triangle RST that are closest to A, B, C */
      struct vector3 A_close_3d, B_close_3d, C_close_3d;
      double dist_A_A_close_3d, dist_B_B_close_3d, dist_C_C_close_3d, min_dist;
      struct vector3  prod_pos_3d;
      struct vector2 prod_pos;

      closest_pt_point_triangle(&A_3d, &R_3d, &S_3d, &T_3d, &A_close_3d);
      closest_pt_point_triangle(&B_3d, &R_3d, &S_3d, &T_3d, &B_close_3d);
      closest_pt_point_triangle(&C_3d, &R_3d, &S_3d, &T_3d, &C_close_3d);

      dist_A_A_close_3d = distance_vec3(&A_3d, &A_close_3d);
      dist_B_B_close_3d = distance_vec3(&B_3d, &B_close_3d);
      dist_C_C_close_3d = distance_vec3(&C_3d, &C_close_3d);

      min_dist = min3d(dist_A_A_close_3d, dist_B_B_close_3d, dist_C_C_close_3d);

      if(!distinguishable(min_dist, dist_A_A_close_3d, EPS_C))
      {
          prod_pos_3d.x = A_close_3d.x;
          prod_pos_3d.y = A_close_3d.y;
          prod_pos_3d.z = A_close_3d.z;
      }else if(!distinguishable(min_dist, dist_B_B_close_3d, EPS_C)){
          prod_pos_3d.x = B_close_3d.x;
          prod_pos_3d.y = B_close_3d.y;
          prod_pos_3d.z = B_close_3d.z;
      }else{
          prod_pos_3d.x = C_close_3d.x;
          prod_pos_3d.y = C_close_3d.y;
          prod_pos_3d.z = C_close_3d.z;
      }  
      
      xyz2uv(&prod_pos_3d, grid1->surface, &prod_pos); 

      if(intersect_point_segment(&prod_pos_3d, &R_3d, &S_3d))
      {
         place_product_close_to_segment_endpoint(&T, &prod_pos, p, k1, k2);
         return; 
      }else if(intersect_point_segment(&prod_pos_3d, &R_3d, &T_3d))
      {
         place_product_close_to_segment_endpoint(&S, &prod_pos, p, k1, k2);
         return; 
      }else if(intersect_point_segment(&prod_pos_3d, &S_3d, &T_3d))
      {
         place_product_close_to_segment_endpoint(&R, &prod_pos, p, k1, k2);
         return; 
      }else{
         p->u = prod_pos.u;
         p->v = prod_pos.v; 
         return;       
      }
    }

    if((num_exact_shared_vertices == 0) && (num_shared_vertices == 0)) 
    {
        mcell_internal_error("Error in the function 'find_closest_position()': number of shared_vertices = 0");
    }
}


/*****************************************************************************
is_inner_tile:
   In: Surface grid
       Index of the tile on that grid
   Out: Returns 1 if the tile is an inner tile 
        (not on the border with the neighbor walls).
        Returns 0 if the tile is on the border with the neighbor wall.
*****************************************************************************/
int is_inner_tile(struct surface_grid *g, int idx)
{
   int root, rootrem, strip, stripe, flip;

   root  = (int)(sqrt((double) idx));
   rootrem = idx - root*root;
   strip = g->n - root -1;
   stripe = rootrem/2;
   flip = rootrem - 2*stripe;

   if((strip == 0) || (stripe == 0)) return 0; 

   if((strip + stripe) == g->n -1) {
      return 0;  
   }

   if(((strip + stripe) == g->n - 2) && (flip == 1)) {
       return 0;  
   }

   return 1;
}

/*****************************************************************************
is_corner_tile:
   In: Surface grid
       Index of the tile on that grid
   Out: Returns 1 if the tile is a corner tile 
        (there are only three corners on each wall).
        Returns 0 if the tile is not a corner tile.
*****************************************************************************/
int is_corner_tile(struct surface_grid *g, int idx)
{
   int root, rootrem, strip, stripe, flip;

   root  = (int)(sqrt((double) idx));
   rootrem = idx - root*root;
   strip = g->n - root -1;
   stripe = rootrem/2;
   flip = rootrem - 2*stripe;

   if((idx == 0) || ((u_int)idx == (g->n_tiles -1))) return 1;

   if((stripe == 0) &&  (strip == 0) && (flip == 0)) return 1;

   return 0;
}

/*****************************************************************************
find_shared_vertices:
   In: Surface grid
       Index of the tile on that grid
       3-member array of indices of parent wall vertices that are shared 
       with other walls (return value)
   Out: Returns 3-member array of the indices of parent wall vertices 
        that are shared with other neighbor walls. 
        If the wall vertex is not shared the corresponding value
        in the return array is not set. 
   Note: Used only for corner tiles.  Here some of the tile vertices
         coincide with the wall vertices which in turn may be shared
         with the neighbor walls.
*****************************************************************************/
void find_shared_vertices(struct surface_grid *sg, int idx, int *shared_vert)
{
   struct vector2 R,S,T; /* tile vertices  in 2D coordinates */  
   struct vector3 R_3d, S_3d, T_3d;  /* tile vertices in 3D coordinates */ 
   int flip;
   
   get_tile_vertices(sg, idx, &flip, &R, &S, &T); 
   uv2xyz(&R, sg->surface, &R_3d); 
   uv2xyz(&S, sg->surface, &S_3d); 
   uv2xyz(&T, sg->surface, &T_3d); 

   /* find tile vertex that coincides with wall vertex */
   if(!distinguishable_vec3(sg->surface->vert[0], &R_3d, EPS_C) ||
      (!distinguishable_vec3(sg->surface->vert[0], &S_3d, EPS_C)) ||
      (!distinguishable_vec3(sg->surface->vert[0], &T_3d, EPS_C)))
   {
      if(sg->surface->vert_0_head != NULL) {
            shared_vert[0] = 0;
      }   
   }

   if(!distinguishable_vec3(sg->surface->vert[1], &R_3d, EPS_C) ||
      (!distinguishable_vec3(sg->surface->vert[1], &S_3d, EPS_C)) ||
      (!distinguishable_vec3(sg->surface->vert[1], &T_3d, EPS_C)))
   {
      if(sg->surface->vert_1_head != NULL) {
            shared_vert[1] = 1;
      }   
   }

   if(!distinguishable_vec3(sg->surface->vert[2], &R_3d, EPS_C) ||
      (!distinguishable_vec3(sg->surface->vert[2], &S_3d, EPS_C)) ||
      (!distinguishable_vec3(sg->surface->vert[2], &T_3d, EPS_C)))
   {
      if(sg->surface->vert_2_head != NULL) {
            shared_vert[2] = 2;
      }   
   }

}


/*****************************************************************************
move_strip_up:
   In: Wall surface grid and tile index on that grid
   Out: Returns index of the tile that is exactly above the tile 
        with index "idx".  Both "idx" and return index should be
        on the same surface grid.
        Returns (-1) if there is no such tile that satisfies the above
        condition (e.g. when the "idx" tile is in the (strip == 0).)
   Note: The direction "up" means moving to the strip where tile indices
         are greater than index "idx".
*****************************************************************************/
int move_strip_up(struct surface_grid *grid, int idx)
{
   int root, rootrem, strip, stripe, flip;
   int num_tiles_per_strip;  
   int tile_up_idx; /* return value */

 
   /* find internal coordinates (strip, stripe, flip) */
   root  = (int)(sqrt((double) idx));
   rootrem = idx - root*root;
   strip = grid->n - root -1;
   stripe = rootrem/2;
   flip = rootrem - 2*stripe;

   num_tiles_per_strip = 2*(grid->n) - 2*strip - 1;

   if(strip == 0) 
   {
      /* tile above is on another wall */
      tile_up_idx = -1;
   }else{
      tile_up_idx = idx + num_tiles_per_strip + 1;
   }
   
   return tile_up_idx;

}

/*****************************************************************************
move_strip_down:
   In: Wall surface grid and tile index on that grid
   Out: Returns index of the tile that is exactly below the tile 
        with index "idx".  Both "idx" and return index should be
        on the same surface grid. Both "idx" and "tile down"
        should share the side.
        Returns (-1) if there is no such tile that satisfies the above
        conditions.
   Note: The direction "down" means moving to the strip where tile indices
         are less than index "idx".
*****************************************************************************/
int move_strip_down(struct surface_grid *grid, int idx)
{
   int root, rootrem, strip, stripe, flip;
   int num_tiles_per_strip;  
   int tile_down_idx; /* return value */

   /* find internal coordinates (strip, stripe, flip) */
   root  = (int)(sqrt((double) idx));
   rootrem = idx - root*root;
   strip = grid->n - root -1;
   stripe = rootrem/2;
   flip = rootrem - 2*stripe;

   num_tiles_per_strip = 2*(grid->n) - 2*strip - 1;

   if(is_inner_tile(grid, idx))
   {
      tile_down_idx = idx - num_tiles_per_strip + 1;
   }else{
     if((strip == 0) && (stripe > 0))
     {
        if((unsigned int)idx == grid->n_tiles - 1)
        {
            tile_down_idx = -1;
        }else{
            tile_down_idx = idx - num_tiles_per_strip + 1;
        }
     }else{
        /* here we are at the left or right border layers */
        if(flip == 0)
        {
            tile_down_idx = -1;
        }else{
            tile_down_idx = idx - num_tiles_per_strip + 1;
        }
     }

   }

   return tile_down_idx;

}



/*****************************************************************************
intersect_point_segment:
   In: point P and segment AB
   Out: 1 if the point P lies on the segment AB, and 0 - otherwise
******************************************************************************/ 
int intersect_point_segment(struct vector3 *P, struct vector3 *A, struct vector3 *B)
{
    struct vector3 ba, pa;
    double t; /* parameter in the line parametrical equation */
    double ba_length, pa_length; /* length of the vectors */
    double cosine_angle; /* cosine of the angle between ba and pa */

    vectorize(A, B, &ba);
    vectorize(A, P, &pa); 
    
    ba_length = vect_length(&ba);
    pa_length = vect_length(&pa);
    
   /* if point intersects segment, vectors pa and ba should be collinear */
    cosine_angle = dot_prod(&ba, &pa)/(ba_length * pa_length);
    if(distinguishable(cosine_angle, 1.0, EPS_C)){
        return 0;
    }

 
   /* Project P on AB, computing parameterized position d(t) = A + t(B - A ) */
   t = dot_prod(&pa, &ba) / dot_prod(&ba, &ba);

   if ((t < 0) || (t > 1)) return 0;

   return 1;

}


/*****************************************************************************
place_product_shared_segment:
   In: segment defined by points R_shared and S_shared
       point T
       position of the product (return value)
       ratios k1 and k2
   Out: coordinates of the point "prod" dividing the distance between T
       and midpont on (R_shared, S_shared) in the ratio k1:k2, so that
       midpoint_prod/midpoint_T = k1/k2
   Note: all points are on the plane 
******************************************************************************/ 

void place_product_shared_segment(struct vector2 *R_shared, struct vector2 *S_shared, struct vector2 *T, struct vector2  *prod, double k1, double k2)
{
    struct vector2 M;   /* midpoint */ 
     
    /* find midpoint on the segment (R_shared - S_shared) */
    M.u = 0.5*(R_shared->u + S_shared->u);
    M.v = 0.5*(R_shared->v + S_shared->v);
   
    /* find coordinates of the  point such that internally
        divides the segment TM in the ratio (k1:k2)  
        We want to place the product close to the common edge.
    */
    prod->u = (k1*T->u + k2*M.u)/(k1 + k2);
    prod->v = (k1*T->v + k2*M.v)/(k1 + k2);

}



/*****************************************************************************
place_product_shared_vertex:
   In: point R_shared
       segment defined by points S and T
       position of the product (return value)
       ratios k1 and k2
   Out: coordinates of the point "prod" dividing the distance between R_shared
       and midpont on (S, T) in the ratio k1:k2, so that 
       R_shared_prod/prod_midpoint = k1/k2
   Note: all points are on the plane 
******************************************************************************/ 
void place_product_shared_vertex(struct vector2 *R_shared, struct vector2 *S, struct vector2 *T, struct vector2* prod, double k1, double k2)
{
   struct vector2 M; /* midpoint */

   /* find midpoint on the segment ST */
   M.u = 0.5*(S->u + T->u);
   M.v = 0.5*(S->v + T->v);
 
   /* find coordinates of the  point such that internally
      divides the segment RM in the ratio (k1:k2)  
   */
   prod->u = (k1*M.u + k2*R_shared->u)/(k1 + k2);
   prod->v = (k1*M.v + k2*R_shared->v)/(k1 + k2);

}


/*****************************************************************************
place_product_close_to_segment_endpoint:
   In: segment defined by points S (start) and E (end)
       position of the product (return value)
       ratios k1 and k2
   Out: coordinates of the point "prod" dividing the distance between S and E
        in the ratio k1:k2, so that E_prod/S_prod = k1/k2
   Note: all points are on the plane
******************************************************************************/ 
void place_product_close_to_segment_endpoint(struct vector2 *S, struct vector2 *E, struct vector2 *prod, double k1, double k2)
{
   prod->u = (k1*S->u + k2*E->u)/(k1 + k2);
   prod->v = (k1*S->v + k2*E->v)/(k1 + k2);
}

/***************************************************************************
append_tile_neigbhor_list:
   In: linked list of neighbor tiles head1
       linked list of neighbor tiles head2
   Out: list "head2" is appended to the list "head1"
        The original "head2" head is set to NULL
***************************************************************************/
void append_tile_neighbor_list(struct tile_neighbor **head1, struct tile_neighbor **head2)
{
   struct tile_neighbor *curr;

   if(*head1 == NULL)    /* special case if "head1" is empty */
   {
      *head1 = *head2;
   }
   else{
      curr = *head1;
      while (curr->next != NULL)  /* find the last node */
      {
         curr = curr->next;     
      }
      curr->next = *head2;
   }
   *head2 = NULL;
}
