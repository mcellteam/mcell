/**************************************************************************\
** File: vol_util.c                                                       **
**                                                                        **
** Purpose: Adds, subtracts, and moves particles around (bookkeeping).    **
**                                                                        **
** Testing status: compiles.  Worked earlier, but has been changed.       **
\**************************************************************************/

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "rng.h"
#include "mem_util.h"
#include "count_util.h"
#include "mcell_structs.h"
#include "vol_util.h"
#include "react.h"
#include "react_output.h"
#include "util.h"
#include "wall_util.h"
#include "grid_util.h"

extern struct volume *world;


/*************************************************************************
inside_subvolume:
  In: pointer to vector3
      pointer to subvolume
  Out: nonzero if the vector is inside the subvolume.
*************************************************************************/

int inside_subvolume(struct vector3 *point,struct subvolume *subvol)
{
  return ( (point->x >= world->x_fineparts[ subvol->llf.x ] ) &&
           (point->x <= world->x_fineparts[ subvol->urb.x ] ) &&
           (point->y >= world->y_fineparts[ subvol->llf.y ] ) &&
           (point->y <= world->y_fineparts[ subvol->urb.y ] ) &&
           (point->z >= world->z_fineparts[ subvol->llf.z ] ) &&
           (point->z <= world->z_fineparts[ subvol->urb.z ] ) );
}


/*************************************************************************
find_coarse_subvolume:
  In: pointer to vector3
  Out: pointer to the coarse subvolume that the vector is within
*************************************************************************/

struct subvolume* find_coarse_subvol(struct vector3 *loc)
{
  int i,j,k;
  i = bisect(world->x_partitions,world->nx_parts,loc->x);
  j = bisect(world->y_partitions,world->ny_parts,loc->y);
  k = bisect(world->z_partitions,world->nz_parts,loc->z);
  return 
    &( world->subvol
      [
        k + (world->nz_parts-1)*(j + (world->ny_parts-1)*i)
      ]
    );
}


/*************************************************************************
traverse_subvol:
  In: pointer to our current subvolume
      pointer to a vector3 of where we want to be
      which direction we're traveling to get there
  Out: subvolume that's closest to where we want to be in our direction
  Note: BSP trees traverse is not yet implemented
*************************************************************************/
struct subvolume* traverse_subvol(struct subvolume *here,struct vector3 *point,int which)
{
    int new_index;

    switch(which)
    {
      case X_NEG:
          if (here->world_edge&X_NEG_BIT) return NULL;
          new_index = here->index - (world->nz_parts - 1)*(world->ny_parts - 1);
          break;
      case X_POS:
          if (here->world_edge&X_POS_BIT) return NULL;
          new_index = here->index + (world->nz_parts - 1)*(world->ny_parts - 1);
          break;
      case Y_NEG:
          if (here->world_edge&Y_NEG_BIT) return NULL;
          new_index = here->index - (world->nz_parts - 1);
          break;
      case Y_POS:
          if (here->world_edge&Y_POS_BIT) return NULL;
          new_index = here->index + (world->nz_parts - 1);
          break;
      case Z_NEG:
          if (here->world_edge&Z_NEG_BIT) return NULL;
          new_index = here->index - 1;
          break;
      case Z_POS:
          if (here->world_edge&Z_POS_BIT) return NULL;
          new_index = here->index + 1;
          break;
      default: 
          fprintf(world->err_file, "Wrong direction caculated in %s, %ld\n", __FILE__, (long)__LINE__);
	  return NULL;
          break;

    } /* end switch */

    return &(world->subvol[new_index]);
    
    /*
  int flag = 1<<which;
  int left_path;
  struct bsp_tree *branch;
  
  if ((here->is_bsp & flag) == 0) return (struct subvolume*)here->neighbor[which]; 
  else
  {
    branch = (struct bsp_tree*) here->neighbor[which];
    while (branch != NULL)
    {
      if ( (branch->flags & X_AXIS) != 0 )
      {
        if ( point->x <= world->x_fineparts[ branch->partition ] ) left_path = 1;
        else left_path = 0;
      }
      else
      {
        if ( (branch->flags & Y_AXIS) != 0 )
        {
          if ( point->y <= world->y_fineparts[ branch->partition ] ) left_path = 1;
          else left_path = 0;
        }
        else // Must be Z_AXIS 
        {
          if ( point->z <= world->z_fineparts[ branch->partition ] ) left_path = 1;
          else left_path = 0;
        }
      }
      if (left_path)
      {
        if ((branch->flags & BRANCH_L) == 0) return (struct subvolume*) branch->left;
        else branch = (struct bsp_tree*) branch->left;
      }
      else
      {
        if ((branch->flags & BRANCH_R) == 0) return (struct subvolume*) branch->right;
        else branch = (struct bsp_tree*) branch->right;
      }
    }
  }

  return NULL;
  */
}

/*************************************************************************
collide_sv_time:
  In: pointer to a vector3 of where we are (*here)
      pointer to a vector3 of where we want to be
      our current subvolume
  Out: time to hit the closest wall of the subvolume
*************************************************************************/

double collide_sv_time(struct vector3 *here,struct vector3 *move,struct subvolume *sv)
{
  double dx,dy,dz,tx,ty,tz,t;
  int whichx,whichy,whichz,which;
  
  whichx = whichy = whichz = 1;
  if (move->x==0 && move->y==0 && move->z==0) return GIGANTIC;
  
  if (move->x > 0) dx = world->x_fineparts[ sv->urb.x ] - here->x;
  else { dx = world->x_fineparts[ sv->llf.x ] - here->x; whichx = 0; }
  
  if (move->y > 0) dy = world->y_fineparts[ sv->urb.y ] - here->y;
  else { dy = world->y_fineparts[ sv->llf.y ] - here->y; whichy = 0; }
  
  if (move->z > 0) dz = world->z_fineparts[ sv->urb.z ] - here->z;
  else { dz = world->z_fineparts[ sv->llf.z ] - here->z; whichz = 0; }
  
  tx = dx * move->y * move->z; if (tx<0) tx = -tx;
  ty = move->x * dy * move->z; if (ty<0) ty = -ty;
  tz = move->x * move->y * dz; if (tz<0) tz = -tz;
  
  if (tx<ty || move->y==0.0)
  {
    if (tx<tz || move->z==0.0) { t = dx / move->x; which = X_NEG + whichx; }
    else { t = dz / move->z; which = Z_NEG + whichz; }
  }
  else /* ty<tx */
  {
    if (ty<tz || move->z==0.0) { t = dy / move->y; which = Y_NEG + whichy; }
    else { t = dz / move->z; which = Z_NEG + whichz; }
  }
  
  return t;
}



/*************************************************************************
next_subvol:
  In: pointer to a vector3 of where we are (*here)
      pointer to a vector3 of where we want to be
      our current subvolume
  Out: next subvolume along that vector or NULL if the endpoint is 
         in the current subvolume.  *here is updated to just inside
         the next subvolume.
*************************************************************************/

struct subvolume* next_subvol(struct vector3 *here,struct vector3 *move,struct subvolume *sv)
{
  double dx,dy,dz,tx,ty,tz,t;
  int whichx,whichy,whichz,which;
  
  whichx = whichy = whichz = 1;
  if (move->x==0 && move->y==0 && move->z==0) return NULL;
  
  if (move->x > 0) dx = world->x_fineparts[ sv->urb.x ] - here->x;
  else { dx = world->x_fineparts[ sv->llf.x ] - here->x; whichx = 0; }
  
  if (move->y > 0) dy = world->y_fineparts[ sv->urb.y ] - here->y;
  else { dy = world->y_fineparts[ sv->llf.y ] - here->y; whichy = 0; }
  
  if (move->z > 0) dz = world->z_fineparts[ sv->urb.z ] - here->z;
  else { dz = world->z_fineparts[ sv->llf.z ] - here->z; whichz = 0; }
  
  tx = dx * move->y * move->z; if (tx<0) tx = -tx;
  ty = move->x * dy * move->z; if (ty<0) ty = -ty;
  tz = move->x * move->y * dz; if (tz<0) tz = -tz;
  
  if (tx<ty || move->y==0.0)
  {
    if (tx<tz || move->z==0.0) { t = dx / move->x; which = X_NEG + whichx; }
    else { t = dz / move->z; which = Z_NEG + whichz; }
  }
  else /* ty<tx */
  {
    if (ty<tz || move->z==0.0) { t = dy / move->y; which = Y_NEG + whichy; }
    else { t = dz / move->z; which = Z_NEG + whichz; }
  }
      
  if (t>=1.0)
  {
    here->x += move->x;
    here->y += move->y;
    here->z += move->z;
    
    return NULL;
  }
  else
  {
    here->x += t*move->x;
    here->y += t*move->y;
    here->z += t*move->z;
    
    t = 1.0-t;
    
    move->x *= t;
    move->y *= t;
    move->z *= t;
    
    return traverse_subvol(sv,here,which); 
  }
}
  


/*************************************************************************
find_subvolume:
  In: pointer to a vector3 of where we are
      pointer to a subvolume we might be in or near
  Out: subvolume that we are in
*************************************************************************/

struct subvolume* find_subvolume(struct vector3 *loc,struct subvolume *guess)
{
#if 1
  /* This code is faster if coarse subvolumes are always used */
  
  if (guess==NULL) return find_coarse_subvol(loc);
  else
  {
    if (world->x_fineparts[guess->llf.x] <= loc->x && loc->x <= world->x_fineparts[guess->urb.x] &&
        world->y_fineparts[guess->llf.y] <= loc->y && loc->y <= world->y_fineparts[guess->urb.y] &&
        world->z_fineparts[guess->llf.z] <= loc->z && loc->z <= world->z_fineparts[guess->urb.z])
    {
      return guess;
    }
    else return find_coarse_subvol(loc);
  }
#else
  /* This code should be used if we ever subdivide subvolumes */
  struct subvolume *sv;
  struct vector3 center;
  
  if (guess == NULL) sv = find_coarse_subvol(loc);
  else sv = guess;

  center.x = 0.5*(world->x_fineparts[ sv->llf.x ] + world->x_fineparts[ sv->urb.x ]);
  center.y = 0.5*(world->y_fineparts[ sv->llf.y ] + world->y_fineparts[ sv->urb.y ]);
  center.z = 0.5*(world->z_fineparts[ sv->llf.z ] + world->z_fineparts[ sv->urb.z ]);
  
  while (loc->x < world->x_fineparts[ sv->llf.x ] )
  {
    sv = traverse_subvol(sv , &center , X_NEG);
    center.x = 0.5*(world->x_fineparts[ sv->llf.x ] + world->x_fineparts[ sv->urb.x ]);
  }
  while (loc->x > world->x_fineparts[ sv->urb.x ] )
  {
    sv = traverse_subvol(sv , &center , X_POS);
    center.x = 0.5*(world->x_fineparts[ sv->llf.x ] + world->x_fineparts[ sv->urb.x ]);
  }
  center.x = loc->x;
  
  while (loc->y < world->y_fineparts[ sv->llf.y ] )
  {
    sv = traverse_subvol(sv , &center , Y_NEG);
    center.y = 0.5*(world->y_fineparts[ sv->llf.y ] + world->y_fineparts[ sv->urb.y ]);
  }
  while (loc->y > world->y_fineparts[ sv->urb.y ] )
  {
    sv = traverse_subvol(sv , &center , Y_POS);
    center.y = 0.5*(world->y_fineparts[ sv->llf.y ] + world->y_fineparts[ sv->urb.y ]);
  }
  center.y = loc->y;

  while (loc->z < world->z_fineparts[ sv->llf.z ] )
  {
    sv = traverse_subvol(sv , &center , Z_NEG);
    center.z = 0.5*(world->z_fineparts[ sv->llf.z ] + world->z_fineparts[ sv->urb.z ]);
  }
  while (loc->z > world->z_fineparts[ sv->urb.z ] )
  {
    sv = traverse_subvol(sv , &center , Z_POS);
    center.z = 0.5*(world->z_fineparts[ sv->llf.z ] + world->z_fineparts[ sv->urb.z ]);
  }
  center.z = loc->z;
  
  return sv;
#endif
}  



/*************************************************************************
is_defunct_molecule
  In: abstract_element that is assumed to be an abstract_molecule
  Out: 0 if the properties field is set, 1 if it is NULL
  Note: This function is passed to sched_util so it can tell which
        molecules are active and which are defunct and can be cleaned up.
*************************************************************************/
int is_defunct_molecule(struct abstract_element *e)
{
  return ((struct abstract_molecule*)e)->properties == NULL;
}


/*************************************************************************
insert_grid_molecule
  In: species for the new molecule
      3D location of the new molecule
      orientation of the new molecule
      diameter to search for a free surface spot (vector3 now, should be double!)
      schedule time for the new molecule
  Out: pointer to the new molecule, or NULL if no free spot was found.
  Note: This function halts the program if it runs out of memory.
*************************************************************************/

struct grid_molecule* insert_grid_molecule(struct species *s,struct vector3 *loc,short orient,double search_diam,double t)
{
  double search_d2,d2;
  struct vector2 s_loc;

  double best_d2;
  struct wall *best_w;
  struct vector2 best_uv;
  struct vector3 best_xyz;
  
  int i0,i1,j0,j1,k0,k1,i,j,k,h;

  struct subvolume *sv;
  struct wall_list *wl;
//  struct wall *w;
  struct grid_molecule *g;
  
  
  if (search_diam<=EPS_C) search_d2 = EPS_C*EPS_C;
  else search_d2 = search_diam * search_diam;
  
  sv = find_subvolume(loc,NULL);
  
  best_d2 = search_d2*2 + 1;
  best_w = NULL;
  for (wl=sv->wall_head ; wl!=NULL ; wl=wl->next)
  {
    d2 = closest_interior_point(loc,wl->this_wall,&s_loc,search_d2);
    if (d2 < search_d2 && d2 < best_d2)
    {
      best_d2 = d2;
      best_w = wl->this_wall;
      best_uv.u = s_loc.u;
      best_uv.v = s_loc.v;
    }
  }

  if (search_d2 > EPS_C*EPS_C)  /* Might need to look in adjacent subvolumes */
  {
    i = sv->index / ((world->ny_parts-1)*(world->nz_parts-1));
    j = sv->index/(world->nz_parts-1) - i*(world->ny_parts-1);
    k = sv->index - (world->nz_parts-1)*(j + i*(world->ny_parts-1));
    
    for (i0=i ; i0>0 ; i0--)
    {
      d2 = loc->x - world->x_partitions[ i0 ]; d2 *= d2;
      if (d2 >= best_d2 || d2 >= search_d2) break;
    }
    for (i1=i ; i1<world->nx_parts-1 ; i1++)
    {
      d2 = loc->x - world->x_partitions[ i1+1 ]; d2 *= d2;
      if (d2 >= best_d2 || d2 >= search_d2) break;
    }
    for (j0=j ; j0>0 ; j0--)
    {
      d2 = loc->y - world->y_partitions[ j0 ]; d2 *= d2;
      if (d2 >= best_d2 || d2 >= search_d2) break;
    }
    for (j1=j ; j1<world->ny_parts-1 ; j1++)
    {
      d2 = loc->y - world->y_partitions[ j1+1 ]; d2 *= d2;
      if (d2 >= best_d2 || d2 >= search_d2) break;
    }
    for (k0=k ; k0>0 ; k0--)
    {
      d2 = loc->z - world->z_partitions[ k0 ]; d2 *= d2;
      if (d2 >= best_d2 || d2 >= search_d2) break;
    }
    for (k1=k ; k1<world->nz_parts-1 ; k1++)
    {
      d2 = loc->z - world->z_partitions[ k1+1 ]; d2 *= d2;
      if (d2 >= best_d2 || d2 >= search_d2) break;
    }
    
    if (i0<i || i1>i || j0<j || j1>j || k0<k || k1>k)
    {
      for (i=i0;i<=i1;i++)
      {
	for (j=j0;j<=j1;j++)
	{
	  for (k=k0;k<=k1;k++)
	  {
	    h = k + (world->nz_parts-1)*(j + (world->ny_parts-1)*i);
	    if (h == sv->index) continue;
	    
	    for (wl=world->subvol[h].wall_head ; wl!=NULL ; wl=wl->next)
	    {
	      d2 = closest_interior_point(loc,wl->this_wall,&s_loc,search_d2);
	      if (d2 <= search_d2 && d2 < best_d2)
	      {
		best_d2 = d2;
		best_w = wl->this_wall;
		best_uv.u = s_loc.u;
		best_uv.v = s_loc.v;
	      }	    
	    }
	  }
	}
      }
      if (best_w!=NULL)
      {
	uv2xyz(&best_uv,best_w,&best_xyz);
	sv = find_subvolume(&best_xyz,sv);  /* May have switched subvolumes */
      }
    }
  }
  
  if (best_w==NULL) return NULL;
  
  d2 = search_d2 - best_d2;  /* We can look this far around the surface we hit for an empty spot */
  
  if (best_w->grid==NULL)
  {
    if (create_grid(best_w,sv))
    {
      fprintf(world->err_file,"File '%s', Line %ld: Out of memory while trying to insert molecules\n", __FILE__, (long)__LINE__);
      i = emergency_output();
      fprintf(world->err_file,"%d errors while attempting emergency output\n",i);
      exit(EXIT_FAILURE);
    }
    i = uv2grid(&best_uv,best_w->grid);
  }
  else
  {
    i = uv2grid(&best_uv,best_w->grid);
    if (best_w->grid->mol[i]!=NULL)
    {
      if (d2 <= EPS_C*EPS_C) return NULL;
      else
      {
	best_w = search_nbhd_for_free(best_w,&best_uv,d2,&i,NULL,NULL);
	if (best_w==NULL) return NULL;
	
	if (world->randomize_gmol_pos) grid2uv_random(best_w->grid,i,&best_uv);
	else grid2uv(best_w->grid,i,&best_uv);
      }
    }
  }
  
  g = mem_get(sv->local_storage->gmol);
  if (g==NULL)
  {
    fprintf(world->err_file,"File '%s', Line %ld: Out of memory while trying to insert molecules\n", __FILE__, (long)__LINE__);
    i = emergency_output();
    fprintf(world->err_file,"%d errors while attempting emergency output\n",i);
    exit(EXIT_FAILURE);
  }

  g->birthplace = sv->local_storage->gmol;
  g->birthday = t;
  g->properties = s;
  s->population++;
  g->flags = TYPE_GRID | ACT_NEWBIE | IN_SCHEDULE;
  if (s->space_step > 0) g->flags |= ACT_DIFFUSE;
  if (trigger_unimolecular(s->hashval,(struct abstract_molecule*)g)!= NULL || (s->flags&CAN_GRIDWALL)!=0 ) g->flags |= ACT_REACT;
  
  g->t = t;
  g->t2 = 0.0;
  g->grid = best_w->grid;
  g->grid_index = i;
  g->s_pos.u = best_uv.u;
  g->s_pos.v = best_uv.v;
  g->orient = orient;
  
  g->grid->n_occupied++;
  
  g->grid->mol[ g->grid_index ] = g;
  
  if (g->properties->flags & (COUNT_CONTENTS|COUNT_ENCLOSED))
    count_region_from_scratch( (struct abstract_molecule*)g , NULL , 1 , NULL , g->grid->surface , g->t );
  
  if ( schedule_add(sv->local_storage->timer,g) )
  {
    fprintf(world->err_file,"File '%s', Line %ld: Out of memory while trying to insert molecules\n", __FILE__, (long)__LINE__);
    i = emergency_output();
    fprintf(world->err_file,"%d errors while attempting emergency output\n",i);
    exit(EXIT_FAILURE);    
  }
  
  return g;
}


/*************************************************************************
insert_volume_molecule
  In: pointer to a volume_molecule that we're going to place in local storage
      pointer to a volume_molecule that may be nearby
  Out: pointer to the new volume_molecule (copies data from volume molecule            passed in), or NULL if out of memory.                                           Molecule is placed in scheduler also.
*************************************************************************/

struct volume_molecule* insert_volume_molecule(struct volume_molecule *m,struct volume_molecule *guess)
{
  struct volume_molecule *new_m;
  struct subvolume *sv;
  
  if (guess == NULL) sv = find_subvolume(&(m->pos),NULL);
  else if ( inside_subvolume(&(m->pos),guess->subvol) ) sv = guess->subvol;
  else sv = find_subvolume(&(m->pos),guess->subvol);
  
  new_m = mem_get(sv->local_storage->mol);
  if(new_m == NULL) {
	fprintf(world->err_file, "File '%s', Line %ld: Out of memory, trying to save intermediate results.\n", __FILE__, (long)__LINE__);        int i = emergency_output();
        fprintf(world->err_file, "Fatal error: out of memory during inserting %s molecule.\nAttempt to write intermediate results had %d errors.\n", m->properties->sym->name, i);
        exit(EXIT_FAILURE);
  }

  memcpy(new_m,m,sizeof(struct volume_molecule));

  new_m->birthplace = sv->local_storage->mol;
  new_m->next = NULL;
  new_m->subvol = sv;
  new_m->next_v = sv->mol_head;
  sv->mol_head = new_m;
  sv->mol_count++;
  new_m->properties->population++;

  if (new_m->properties->flags & (COUNT_CONTENTS|COUNT_ENCLOSED))
  {
    count_region_from_scratch( (struct abstract_molecule*)new_m , NULL , 1 , &(new_m->pos) , NULL , new_m->t );
  }
  
  if ( schedule_add(sv->local_storage->timer,new_m) ) {
	fprintf(world->err_file, "File '%s', Line %ld: Out of memory, trying to save intermediate results.\n", __FILE__, (long)__LINE__);
        int i = emergency_output();
        fprintf(world->err_file, "Fatal error: out of memory during inserting %s molecule.\nAttempt to write intermediate results had %d errors.\n", m->properties->sym->name, i);
        exit(EXIT_FAILURE);

  } 
  return new_m;
}


/*************************************************************************
excert_volume_molecule:
  In: pointer to a volume_molecule that we're going to remove from local storage
  Out: no return value; molecule is marked for removal.
*************************************************************************/

void excert_volume_molecule(struct volume_molecule *m)
{
  if (m->properties->flags & (COUNT_CONTENTS|COUNT_ENCLOSED))
  {
    count_region_from_scratch( (struct abstract_molecule*)m , NULL,  -1 , NULL , NULL , m->t );
  }
  m->subvol->mol_count--;
  m->properties->n_deceased++;
  m->properties->cum_lifetime += m->t - m->birthday;
  m->properties->population--;
  m->properties = NULL;
}


/*************************************************************************
insert_volume_molecule_list:
  In: pointer to a linked list of volume_molecules to copy into subvolumes.
  Out: 0 on success, 1 on memory allocation error; molecules are placed
       in their subvolumes.
*************************************************************************/

int insert_volume_molecule_list(struct volume_molecule *m)
{
  struct volume_molecule *new_m,*guess;
  
  guess=NULL;
  while (m != NULL)
  {
    new_m = insert_volume_molecule(m,guess);
    if(new_m == NULL) { 
	fprintf(world->err_file, "File '%s', Line %ld: Out of memory, trying to save intermediate results.\n", __FILE__, (long)__LINE__);
        int i = emergency_output();
        fprintf(world->err_file, "Fatal error: out of memory during inserting %s molecule.\nAttempt to write intermediate results had %d errors.\n", m->properties->sym->name, i);
        exit(EXIT_FAILURE);
    }
    guess = new_m;
    m = (struct volume_molecule*)m->next;
  }
  
  return 0;
}


/*************************************************************************
migrate_volume_molecule:
  In: pointer to a volume_molecule already in a subvolume
      pointer to the new subvolume to move it to
  Out: pointer to moved molecule.  The molecule's position is updated
       but it is not rescheduled.  Returns NULL if out of memory.
*************************************************************************/

struct volume_molecule* migrate_volume_molecule(struct volume_molecule *m,struct subvolume *new_sv)
{
  struct volume_molecule *new_m;

  new_m = mem_get(new_sv->local_storage->mol);
  if (new_m==NULL){ 
	fprintf(world->err_file, "File '%s', Line %ld: Out of memory, trying to save intermediate results.\n", __FILE__, (long)__LINE__);        int i = emergency_output();
        fprintf(world->err_file, "Fatal error: out of memory during migrating  %s molecule.\nAttempt to write intermediate results had %d errors.\n", m->properties->sym->name, i);
        exit(EXIT_FAILURE);
  }
  if (new_m==m) printf("File '%s', Line %ld: Unexpected behavior!\n", __FILE__, (long)__LINE__);
  
  memcpy(new_m,m,sizeof(struct volume_molecule));
  new_m->birthplace = new_sv->local_storage->mol;
  
  new_m->next = NULL;
  new_m->subvol = new_sv;
  new_m->next_v = new_sv->mol_head;
  new_sv->mol_head = new_m;
  new_sv->mol_count++;
 
  m->subvol->mol_count--;
  m->properties = NULL;

  return new_m;    
}


/*************************************************************************
eval_rel_region_3d:
  In: an expression tree containing regions to release on
      the waypoint for the current subvolume
      a list of regions entered from the waypoint to the release loc.
      a list of regions exited from the waypoint to the release loc.
  Out: 1 if the location chosen satisfies the expression, 0 if not.
*************************************************************************/

int eval_rel_region_3d(struct release_evaluator *expr,struct waypoint *wp,struct region_list *in_regions,struct region_list *out_regions)
{
  struct region *r;
  struct region_list *rl;
  int satisfies_l,satisfies_r;
  
  satisfies_l=0;
  if (expr->op & REXP_LEFT_REGION)
  {
    r = (struct region*)expr->left;
    for (rl=wp->regions ; rl!=NULL ; rl=rl->next)
    {
      if (rl->reg == r)
      {
        satisfies_l=1;
        break;
      }
    }
    if (satisfies_l)
    {
      for (rl=out_regions ; rl!=NULL ; rl=rl->next)
      {
        if (rl->reg==r)
        {
          satisfies_l=0;
          break;
        }
      }
    }
    else
    {
      for (rl=in_regions ; rl!=NULL ; rl=rl->next)
      {
        if (rl->reg==r)
        {
          satisfies_l=1;
          break;
        }
      }
    }
  }
  else satisfies_l = eval_rel_region_3d(expr->left,wp,in_regions,out_regions);
  
  if (expr->op & REXP_NO_OP) return satisfies_l;
  
  satisfies_r=0;
  if (expr->op & REXP_RIGHT_REGION)
  {
    r = (struct region*)expr->right;
    for (rl=wp->regions ; rl!=NULL ; rl=rl->next)
    {
      if (rl->reg == r)
      {
        satisfies_r=1;
        break;
      }
    }
    if (satisfies_r)
    {
      for (rl=out_regions ; rl!=NULL ; rl=rl->next)
      {
        if (rl->reg==r)
        {
          satisfies_r=0;
          break;
        }
      }
    }
    else
    {
      for (rl=in_regions ; rl!=NULL ; rl=rl->next)
      {
        if (rl->reg==r)
        {
          satisfies_r=1;
          break;
        }
      }
    }
  }
  else satisfies_r = eval_rel_region_3d(expr->right,wp,in_regions,out_regions);
  
  if (expr->op & REXP_UNION) return (satisfies_l || satisfies_r);
  else if (expr->op & REXP_INTERSECTION) return (satisfies_l && satisfies_r);
  else if (expr->op & REXP_SUBTRACTION) return (satisfies_l && !satisfies_r);

  return 0;
}


/*************************************************************************
vacuum_inside_regions:
  In: pointer to a release site object
      template molecule to remove
      integer number of molecules to remove (negative)
  Out: 0 on success, 1 on failure; molecule(s) are removed from the
       world as specified.
  Note: if more molecules are to be removed than actually exist, all
        existing molecules of the specified type are removed.
*************************************************************************/

int vacuum_inside_regions(struct release_site_obj *rso,struct volume_molecule *m,int n)
{
  struct volume_molecule *mp;
  struct release_region_data *rrd;
  struct region_list *extra_in,*extra_out;
  struct region_list *rl,*rl2;
  struct waypoint *wp;
  struct subvolume *sv = NULL;
  int i1,i2,j1,j2,k1,k2;
  int h,i,j,k,l;
  struct mem_helper *mh;
  struct void_list *vl;
  struct void_list *vl_head = NULL;
  int vl_num = 0;
  double t;
  struct vector3 hit,delta;
  struct vector3 *origin;
  struct wall_list *wl;
  
  rrd = rso->region_data;
  
  mh = create_mem(sizeof(struct void_list),1024);
  if (mh==NULL) return 1;
  
  i1 = bisect(world->x_partitions,world->nx_parts,rrd->llf.x);
  i2 = bisect_high(world->x_partitions,world->nx_parts,rrd->urb.x);
  j1 = bisect(world->y_partitions,world->ny_parts,rrd->llf.y);
  j2 = bisect_high(world->y_partitions,world->ny_parts,rrd->urb.y);
  k1 = bisect(world->z_partitions,world->nz_parts,rrd->llf.z);
  k2 = bisect_high(world->z_partitions,world->nz_parts,rrd->urb.z);
  
  for (i=i1;i<i2;i++)
  {
    for (j=j1;j<j2;j++)
    {
      for (k=k1;k<k2;k++)
      {
        h = k + (world->nz_parts - 1)*(j + (world->ny_parts-1)*i);
        sv = &(world->subvol[h]);
        
        for (mp=sv->mol_head ; mp!=NULL ; mp=mp->next_v)
        {
          if (mp->properties==m->properties)
          {
            extra_in=extra_out=NULL;
            wp = &(world->waypoints[sv->index]);
            origin = &(wp->loc);
            delta.x = mp->pos.x - origin->x;
            delta.y = mp->pos.y - origin->y;
            delta.z = mp->pos.z - origin->z;
            
            for (wl=sv->wall_head ; wl!=NULL ; wl=wl->next)
            {
              l = collide_wall(origin,&delta,wl->this_wall,&t,&hit);
              
              if (l!=COLLIDE_MISS)
              {
                if(world->notify->final_summary == NOTIFY_FULL) {
                   world->ray_polygon_colls++;
                }
                
                for (rl=wl->this_wall->counting_regions ; rl!=NULL ; rl=rl->next)
                {
                  if (l==COLLIDE_FRONT || l==COLLIDE_BACK)
                  {
                    rl2 = (struct region_list*)mem_get(sv->local_storage->regl);
                    rl2->reg = rl->reg;
                    
                    if (l==COLLIDE_FRONT)
                    {
                      rl2->next = extra_in;
                      extra_in = rl2;
                    }
                    else  /*l==COLLIDE_BACK*/
                    {
                      rl2->next = extra_out;
                      extra_out = rl2;
                    }
                  }
                }
              }
            }
            
            for (rl=extra_in ; rl!=NULL ; rl=rl->next)
            {
              if (rl->reg==NULL) continue;
              for (rl2=extra_out ; rl2!=NULL ; rl2=rl2->next)
              {
                if (rl2->reg==NULL) continue;
                if (rl->reg==rl2->reg)
                {
                  rl->reg = NULL;
                  rl2->reg = NULL;
                  break;
                }
              }
            }

            l = eval_rel_region_3d(rrd->expression,wp,extra_in,extra_out);
            
            if (l)
            {
              vl = (struct void_list*)mem_get(mh);
              if (vl==NULL) return 1;
              
              vl->data = mp;
              vl->next = vl_head;
              vl_head = vl;
              vl_num++;
            }
	    
	    if (extra_in!=NULL) mem_put_list(sv->local_storage->regl,extra_in);
	    if (extra_out!=NULL) mem_put_list(sv->local_storage->regl,extra_out);
          }
        }
      }
    }
  }
  
  for ( vl=vl_head ; n<0 && vl_num>0 && vl!=NULL ; vl=vl->next , vl_num--)
  {
    if ( rng_dbl(world->rng) < ((double)(-n))/((double)vl_num) )
    {
      if(world->notify->final_summary == NOTIFY_FULL){
         world->random_number_use++;
      }
      mp = (struct volume_molecule*)vl->data;
      mp->properties->population--;
      mp->subvol->mol_count--;
      if ((mp->properties->flags & (COUNT_CONTENTS|COUNT_ENCLOSED)) != 0)
        count_region_from_scratch((struct abstract_molecule*)mp,NULL,-1,&(mp->pos),NULL,mp->t);
      mp->properties = NULL;
      if (mp->flags & IN_SCHEDULE)
      {
        mp->subvol->local_storage->timer->defunct_count++; /* Tally for garbage collection */
      }
      
      n++;
    }
  }
  
  delete_mem(mh);
  return 0;
}


/*************************************************************************
release_inside_regions:
  In: pointer to a release site object
      template molecule to release
      integer number of molecules to release
  Out: 0 on success, 1 on failure; molecule(s) are released into the
       world as specified.
  Note: if the CCNNUM release method is used, the number of molecules
        passed in is ignored.
*************************************************************************/

int release_inside_regions(struct release_site_obj *rso,struct volume_molecule *m,int n)
{
  struct volume_molecule *new_m;
  struct release_region_data *rrd;
  struct region_list *extra_in,*extra_out;
  struct region_list *rl,*rl2;
  struct subvolume *sv = NULL;
  struct wall_list *wl;
  struct vector3 delta;
  struct vector3 *origin;
  struct waypoint *wp;
  double t;
  struct vector3 hit;
  int bad_location;
  int i;
  
  rrd = rso->region_data;
  new_m = NULL;
  m->previous_wall = NULL;
  m->index = -1;
  
  if (rso->release_number_method==CCNNUM)
  {
    double vol = (rrd->urb.x-rrd->llf.x)*(rrd->urb.y-rrd->llf.y)*(rrd->urb.z-rrd->llf.z);
    n = (int)(0.5 + N_AV*1e-15*rso->concentration*vol*world->length_unit*world->length_unit*world->length_unit);
  }
  
  if (n<0) return vacuum_inside_regions(rso,m,n);
  
  while (n>0)
  {
    m->pos.x = rrd->llf.x + (rrd->urb.x-rrd->llf.x)*rng_dbl(world->rng);
    if(world->notify->final_summary == NOTIFY_FULL){
       world->random_number_use++;
    }
    m->pos.y = rrd->llf.y + (rrd->urb.y-rrd->llf.y)*rng_dbl(world->rng);
    if(world->notify->final_summary == NOTIFY_FULL){
       world->random_number_use++;
    }
    m->pos.z = rrd->llf.z + (rrd->urb.z-rrd->llf.z)*rng_dbl(world->rng);
    if(world->notify->final_summary == NOTIFY_FULL){
       world->random_number_use++;
    }
    
    if (sv == NULL) sv = find_subvolume(&(m->pos),NULL);
    else if ( !inside_subvolume(&(m->pos),sv) )
    {
      sv = find_subvolume(&(m->pos),sv);
    }
    
    extra_in=extra_out=NULL;
    wp = &(world->waypoints[sv->index]);
    origin = &(wp->loc);
    delta.x = m->pos.x - origin->x;
    delta.y = m->pos.y - origin->y;
    delta.z = m->pos.z - origin->z;
    
    bad_location = 0;
    
    for (wl=sv->wall_head ; wl!=NULL ; wl=wl->next)
    {
      i = collide_wall(origin,&delta,wl->this_wall,&t,&hit);
      
      if (i!=COLLIDE_MISS)
      {
        if(world->notify->final_summary == NOTIFY_FULL) {
            world->ray_polygon_colls++;
        }

        if ( (t>-EPS_C && t<EPS_C) || (t>1.0-EPS_C && t<1.0+EPS_C) )
        {
          bad_location = 1;
          break;
        }
        for (rl=wl->this_wall->counting_regions ; rl!=NULL ; rl=rl->next)
        {
          rl2 = (struct region_list*)mem_get(sv->local_storage->regl);
          rl2->reg = rl->reg;
          
          if (i==COLLIDE_FRONT)
          {
            rl2->next = extra_in;
            extra_in = rl2;
          }
          else if (i==COLLIDE_BACK)
          {
            rl2->next = extra_out;
            extra_out = rl2;
          }
          else
          {
            bad_location = 1;
            break;
          }
        }
      }
    }
    if (bad_location)
    {
      if (extra_in!=NULL) mem_put_list(sv->local_storage->regl,extra_in);
      if (extra_out!=NULL) mem_put_list(sv->local_storage->regl,extra_out);
      continue;
    }
    
    for (rl=extra_in ; rl!=NULL ; rl=rl->next)
    {
      if (rl->reg==NULL) continue;
      for (rl2=extra_out ; rl2!=NULL ; rl2=rl2->next)
      {
        if (rl2->reg==NULL) continue;
        if (rl->reg==rl2->reg)
        {
          rl->reg = NULL;
          rl2->reg = NULL;
          break;
        }
      }
    }

    i = eval_rel_region_3d(rrd->expression,wp,extra_in,extra_out);

    if (extra_in!=NULL) mem_put_list(sv->local_storage->regl,extra_in);
    if (extra_out!=NULL) mem_put_list(sv->local_storage->regl,extra_out);
    
    if (!i)
    {
      if (rso->release_number_method==CCNNUM) n--;
      continue;
    }
    
    m->subvol = sv;
    new_m =  insert_volume_molecule(m,new_m);
    
    if (new_m==NULL) return 1;
    
    n--;
  }
  
  return 0;
}


/*************************************************************************
release_molecules:
  In: pointer to a release event
  Out: 0 on success, 1 on failure; next event is scheduled and molecule(s)
       are released into the world as specified.
  Note: if a release is triggered by a reaction, there isn't anything
        to schedule.  Also, in that case, rpat isn't really a release
	pattern (it's a rxn_pathname in disguise) so be sure to not
	dereference it!
*************************************************************************/
int release_molecules(struct release_event_queue *req)
{
  struct release_site_obj *rso;
  struct release_pattern *rpat;
  struct volume_molecule m;
  struct grid_molecule g;
  struct abstract_molecule *ap;
  struct volume_molecule *mp;
  struct grid_molecule *gp;
  struct volume_molecule *guess;
  int i,number;
  short orient;
  struct vector3 *diam_xyz;
  struct vector3 pos;
  double diam,vol;
  double t,k;
  struct release_single_molecule *rsm;
  double location[1][4];
  
  if(req == NULL) return 0;
  rso = req->release_site;
  rpat = rso->pattern;

  /* Set up canonical molecule to be released */
  /* If we have a list, assume a 3D molecule and fix later */
  if (rso->mol_list!=NULL || (rso->mol_type->flags & NOT_FREE)==0)
  {
    mp = &m;  /* Avoid strict-aliasing message */
    ap = (struct abstract_molecule*)mp;
    ap->flags = TYPE_3D | IN_VOLUME;
  }
  else
  {
    gp = &g;  /* Avoid strict-aliasing message */
    ap = (struct abstract_molecule*)gp;
    ap->flags = TYPE_GRID | IN_SURFACE;
  }
  ap->flags |= IN_SCHEDULE + ACT_NEWBIE;

  if(req->train_counter == 0)
  {
	req->train_counter++;
  }
  
  guess = NULL;

  /* Skip events that happened in the past (delay<0 or after checkpoint) */
  if( req->event_time < world->it_time && rso->release_prob!=MAGIC_PATTERN_PROBABILITY)
  {
    do
    {
      /* Schedule next release event and leave the function. 
	 This part of the code is relevant to checkpointing. */
      if (rso->release_prob < 1.0)
      {
	k = -log( 1.0 - rso->release_prob );
	t = -log( rng_dbl(world->rng) ) / k;  /* Poisson dist. */
        if(world->notify->final_summary == NOTIFY_FULL){
           world->random_number_use++;
        }
	req->event_time += rpat->release_interval * (ceil(t)-1.0); /* Rounded to integers */
      }
      else
      {
	req->event_time += rpat->release_interval;
      }
      /* we may need to move to the next train. */
      if (!distinguishable(req->event_time,req->train_high_time + rpat->train_duration,EPS_C) ||
	   req->event_time > req->train_high_time + rpat->train_duration)
      {
	req->train_high_time += rpat->train_interval;
	req->event_time = req->train_high_time;
	req->train_counter++;
      }
    } while(req->event_time <= world->start_time); 
  
    if (req->train_counter <= rpat->number_of_trains && req->event_time < FOREVER)
    {
      if ( schedule_add(world->releaser,req) )
      {
	fprintf(world->err_file, "File '%s', Line %ld: Out of memory, trying to save intermediate results.\n", __FILE__, (long)__LINE__);
	int i = emergency_output();
	fprintf(world->err_file, "Fatal error: out of memory during release molecule event.\nAttempt to write intermediate results had %d errors.\n", i);
	exit(EXIT_FAILURE);
      } 
    }
    return 0;
  }

  
  /* Set molecule characteristics. */
  ap->t = req->event_time;
  ap->properties = rso->mol_type;
  ap->t2 = 0.0;
  ap->birthday = ap->t;
  
  if (rso->mol_list==NULL)  /* All molecules are the same, so we can set flags */
  {
    if (trigger_unimolecular(rso->mol_type->hashval , ap) != NULL || (rso->mol_type->flags&CAN_GRIDWALL)!=0) ap->flags |= ACT_REACT;
    if (rso->mol_type->space_step > 0.0) ap->flags |= ACT_DIFFUSE;
  }
  
  switch(rso->release_number_method)
  {
    case CONSTNUM:
      number = rso->release_number;
      break;
    case GAUSSNUM:
      if (rso->standard_deviation > 0)
      {
	number = (int) rng_gauss(world->rng)*rso->standard_deviation + rso->release_number;
      }
      else
      {
        rso->release_number_method = CONSTNUM;
        number = rso->release_number;
      }
      break;
    case VOLNUM:
      diam = rso->mean_diameter;
      if (rso->standard_deviation > 0)
      {
	diam += rng_gauss(world->rng)*rso->standard_deviation;
      }
      vol = (MY_PI/6.0) * diam*diam*diam;
      number = (int)(N_AV * 1e-15 * rso->concentration * vol + 0.5);
      break;
    case CCNNUM:
      if (rso->diameter==NULL) number = 0;
      else
      {
        switch(rso->release_shape)
        {
          case SHAPE_SPHERICAL:
          case SHAPE_ELLIPTIC:
            vol = (1.0/6.0)*MY_PI*rso->diameter->x*rso->diameter->y*rso->diameter->z*(world->length_unit*world->length_unit*world->length_unit);
            break;
          case SHAPE_RECTANGULAR:
          case SHAPE_CUBIC:
            vol = rso->diameter->x*rso->diameter->y*rso->diameter->z*(world->length_unit*world->length_unit*world->length_unit);
            break;
          default:
            fprintf(world->err_file,"File '%s', Line %ld: Can't release a concentration on a spherical shell\n", __FILE__, (long)__LINE__);
            vol = 0;
            break;
        }
        number = (int)(N_AV * 1e-15 * rso->concentration * vol + 0.5);
      }
      break;
    
    default:
      number = 0;
      break;
  }
  
  if (rso->release_shape == SHAPE_REGION)
  {
    int pop_before = ap->properties->population;
    if (ap->flags & TYPE_3D)
    {
      i = release_inside_regions(rso,(struct volume_molecule*)ap,number);
      if (i) return 1;
      
      if (world->notify->release_events==NOTIFY_FULL)
      {
        if (number>0 || (rso->release_number_method==CCNNUM && rso->concentration>0))
        {
          fprintf(world->log_file, "Releasing %d %s at iteration %lld\n",
            ap->properties->population-pop_before,req->release_site->mol_type->sym->name,world->it_time);
        }
        else
        {
          fprintf(world->log_file, "Unreleasing %d %s at iteration %lld\n",
            ap->properties->population-pop_before,req->release_site->mol_type->sym->name,world->it_time);
        }
      }
    }
    else
    {
      i = release_onto_regions(rso,(struct grid_molecule*)ap,number);
      if (i) return 1;

      if (world->notify->release_events==NOTIFY_FULL)
      {
        if (number>0)
        {
          fprintf(world->log_file, "Releasing %d %s at iteration %lld\n",
            ap->properties->population-pop_before,req->release_site->mol_type->sym->name,world->it_time);
        }
        else
        {
          fprintf(world->log_file, "Unreleasing %d %s at iteration %lld\n",
            ap->properties->population-pop_before,req->release_site->mol_type->sym->name,world->it_time);
        }
      }
    }
  }
  else  /* Guaranteed to be 3D molecule or at least specified by 3D location if in list */
  {
    m.previous_wall = NULL;
    m.index = -1;
    
    diam_xyz = rso->diameter;
    rsm = rso->mol_list;
    if (rsm != NULL)
    {
      for (i=0 ; rsm!=NULL ; rsm=rsm->next,i++)
      {
	location[0][0] = rsm->loc.x + rso->location->x;
	location[0][1] = rsm->loc.y + rso->location->y;
	location[0][2] = rsm->loc.z + rso->location->z;
	location[0][3] = 1;
	
	mult_matrix(location,req->t_matrix,location,1,4,4);
	
	m.pos.x = location[0][0];
	m.pos.y = location[0][1];
	m.pos.z = location[0][2];
	  
	if ((rsm->mol_type->flags & NOT_FREE)==0)
	{
	  m.properties = rsm->mol_type;
	  guess = insert_volume_molecule(&m,guess);
	  if (guess==NULL) return 1;
	}
	else
	{
          if (diam_xyz==NULL) diam=0.0;
          else diam=diam_xyz->x;
	  
	  if (rsm->orient>0) orient=1;
	  else if (rsm->orient<0) orient=-1;
	  else {
             orient = (rng_uint(world->rng)&1)?1:-1;
             if(world->notify->final_summary == NOTIFY_FULL){
                world->random_number_use++;
             }
          }

	  gp = insert_grid_molecule(rsm->mol_type,&(m.pos),orient,diam,req->event_time);
	  if (gp==NULL)
	  {
	    fprintf(world->log_file,"Warning: unable to find surface upon which to place molecule %s\n",rsm->mol_type->sym->name);
	    fprintf(world->log_file,"  Perhaps you want to SITE_DIAMETER larger to increase search distance?\n");
	  }
	}
      }
      if (world->notify->release_events==NOTIFY_FULL)
      {
        fprintf(world->log_file, "Releasing %d molecules from list at iteration %lld\n",
            i,world->it_time);
      }
    }
    else if (diam_xyz != NULL)
    {
      
      for (i=0;i<number;i++)
      {
	do /* Pick values in unit square, toss if not in unit circle */
	{
	  pos.x = (rng_dbl(world->rng)-0.5);
          if(world->notify->final_summary == NOTIFY_FULL){
             world->random_number_use++;
          }
	  pos.y = (rng_dbl(world->rng)-0.5);
          if(world->notify->final_summary == NOTIFY_FULL){
             world->random_number_use++;
          }
	  pos.z = (rng_dbl(world->rng)-0.5);
          if(world->notify->final_summary == NOTIFY_FULL){
             world->random_number_use++;
          }
	} while ( (rso->release_shape == SHAPE_SPHERICAL || rso->release_shape == SHAPE_ELLIPTIC || rso->release_shape == SHAPE_SPHERICAL_SHELL)
		  && pos.x*pos.x + pos.y*pos.y + pos.z*pos.z >= 0.25 );
	
	if (rso->release_shape == SHAPE_SPHERICAL_SHELL)
	{
	  double r;
	  r = sqrt( pos.x*pos.x + pos.y*pos.y + pos.z*pos.z)*2.0;
	  if (r==0.0) { pos.x = 0.0; pos.y = 0.0; pos.z = 0.5; }
	  else { pos.x /= r; pos.y /= r; pos.z /= r; }
	}
	
	location[0][0] = pos.x*diam_xyz->x + rso->location->x;
	location[0][1] = pos.y*diam_xyz->y + rso->location->y;
	location[0][2] = pos.z*diam_xyz->z + rso->location->z;
	location[0][3] = 1;
        
        mult_matrix(location,req->t_matrix,location,1,4,4);
        
        m.pos.x = location[0][0];
        m.pos.y = location[0][1];
        m.pos.z = location[0][2];
        
        guess = insert_volume_molecule(&m,guess);  /* Insert copy of m into world */
        if (guess == NULL) return 1;
      }
      if (world->notify->release_events==NOTIFY_FULL)
      {
        fprintf(world->log_file, "Releasing %d %s at iteration %lld\n",
            number,req->release_site->mol_type->sym->name,world->it_time);
      }
    }
    else
    {
      location[0][0] = rso->location->x;
      location[0][1] = rso->location->y;
      location[0][2] = rso->location->z;
      location[0][3] = 1;
      
      mult_matrix(location,req->t_matrix,location,1,4,4);
      
      m.pos.x = location[0][0];
      m.pos.y = location[0][1];
      m.pos.z = location[0][2];
      
      for (i=0;i<number;i++)
      {
         guess = insert_volume_molecule(&m,guess);
         if (guess == NULL) return 1;
      }
      if (world->notify->release_events==NOTIFY_FULL)
      {
        fprintf(world->log_file, "Releasing %d %s at iteration %lld\n",
            number,req->release_site->mol_type->sym->name,world->it_time);
      }
    }
  }
 
  
  /* Schedule next release event. */
  if (rso->release_prob==MAGIC_PATTERN_PROBABILITY) return 0;  /* Triggered by reaction, don't schedule */
  if (rso->release_prob < 1.0)
  {
    k = -log( 1.0 - rso->release_prob );
    t = -log( rng_dbl(world->rng) ) / k;  /* Poisson dist. */
    if(world->notify->final_summary == NOTIFY_FULL){
        world->random_number_use++;
    }
    req->event_time += rpat->release_interval * (ceil(t)-1.0); /* Rounded to integers */
  }
  else if (rso->release_prob==MAGIC_PATTERN_PROBABILITY)
  {
    req->event_time += FOREVER;
  }
  else
  {
    req->event_time += rpat->release_interval;
  }
    /* we may need to move to the next train. */
  if (!distinguishable(req->event_time,req->train_high_time + rpat->train_duration,EPS_C) ||
       req->event_time > req->train_high_time + rpat->train_duration)
  {
    req->train_high_time += rpat->train_interval;
    req->event_time = req->train_high_time;
    req->train_counter++;
  }
    
  if (req->train_counter <= rpat->number_of_trains && req->event_time < FOREVER)
  {
    if ( schedule_add(world->releaser,req) ){
      fprintf(world->err_file, "File '%s', Line %ld: Out of memory, trying to save intermediate results.\n", __FILE__, (long)__LINE__);
      int i = emergency_output();
      fprintf(world->err_file, "Fatal error: out of memory during release molecule event.\nAttempt to write intermediate results had %d errors.\n", i);
      exit(EXIT_FAILURE);
    } 
  }
 
  return 0;
}



/*************************************************************************
find_exponential_params:
  In: value of f(0)
      value of f(N)
      difference between f(1) and f(0)
      number of data points
      pointer to where we store the scaling factor A
      pointer to the constant offset B
      pointer to the rate of decay k
  Out: no return value.  This is a utility function that uses bisection
       to solve for A,B,k to find an exponentially increasing function
         f(n) = A*exp(n*k)+B
       subject to the contstraints
         f(0) = c
         f(1) = c+d
         f(N) = C
*************************************************************************/

void find_exponential_params(double c,double C,double d,double N,double *A,double *B, double *k)
{
  double k_min,k_max,k_mid,f;
  int i;
  
  k_min = 0;
  k_max = log(GIGANTIC)/N;
  for (i=0;i<720;i++)
  {
    k_mid = 0.5*(k_min + k_max);
    f = c + (exp(N*k_mid)-1.0)*d/(exp(k_mid)-1.0);
    if (C > f) k_min = k_mid;
    else k_max = k_mid;
    if ((k_max-k_min)/(k_max+k_min) < EPS_C) break;
  }
  
  *k = k_mid;
  *A = d / ( exp(*k) - 1.0 );
  *B = c - *A;
}


/*************************************************************************
set_partitions:
  In: nothing.  Uses struct volume *world, assumes bounding box is set.
  Out: 0 on success, 1 on error; coarse and fine partitions are set.
*************************************************************************/

/*FIXME: I am impossible to understand.  Comprehensibilize this function!*/
int set_partitions()
{
  double f_min,f_max,f,df,dfx,dfy,dfz;
  int i,j;
  double steps_min,steps_max;
  double x_aspect,y_aspect,z_aspect;
  int x_in,y_in,z_in;
  int x_start,y_start,z_start;
  double A,B,k;
  struct vector3 part_min,part_max;
  double smallest_spacing;

  /* Set sensible bounds for spacing between fine partitions (minimum size of subdivision) */  
  smallest_spacing = 0.1/world->length_unit;  /* 100nm */

  if (2*world->rx_radius_3d > smallest_spacing) smallest_spacing=2*world->rx_radius_3d;

  /* We have 2^15 possible fine partitions; we'll use 24k of them */
  if (world->n_fineparts != 4096 + 16384 + 4096)
  {
    world->n_fineparts = 4096 + 16384 + 4096;
    world->x_fineparts = (double*)malloc(sizeof(double)*world->n_fineparts);
    world->y_fineparts = (double*)malloc(sizeof(double)*world->n_fineparts);
    world->z_fineparts = (double*)malloc(sizeof(double)*world->n_fineparts);
  }
  if((world->x_fineparts == NULL) || (world->y_fineparts == NULL) ||
        (world->z_fineparts == NULL))
  {
    fprintf(world->err_file, "File '%s', Line %ld: Out of memory while trying to create partitions.\n", __FILE__, (long)__LINE__);
    return 1;
  }

  /* Something like the maximum expected error--not sure exactly what this is */
  dfx = 1e-3 + (world->bb_urb.x - world->bb_llf.x)/8191.0;
  dfy = 1e-3 + (world->bb_urb.y - world->bb_llf.y)/8191.0;
  dfz = 1e-3 + (world->bb_urb.z - world->bb_llf.z)/8191.0;

  /* Not sure how this is supposed to work--looks like two ideas mixed, probably broken */
  /* Was supposed to make sure that the fine partitions would still obey the 2*reaction radius rule */
  f_min = world->bb_llf.x - dfx;
  f_max = world->bb_urb.x + dfx;
  if (f_max - f_min < smallest_spacing)
  {
    printf("Rescaling: was %.3f to %.3f, now ",f_min*world->length_unit,f_max*world->length_unit);
    f = smallest_spacing - (f_max-f_min);
    f_max += 0.5*f;
    f_min -= 0.5*f;
    printf("%.3f to %.3f\n",f_min*world->length_unit,f_max*world->length_unit);
  }
  /* Set bounds over which to do linear subdivision (world bounding box) */
  part_min.x = f_min;
  part_max.x = f_max;
  df = (f_max - f_min)/16383.0;
  /* Subdivide world bounding box */
  for (i=0;i<16384;i++)
  {
    world->x_fineparts[ 4096 + i ] = f_min + df*((double)i);
  }
  /* Create an exponentially increasing fine partition size as we go to -infinity */
  find_exponential_params(-f_min,1e12,df,4096,&A,&B,&k);
  for (i=1;i<=4096;i++) world->x_fineparts[4096-i] = -(A*exp(i*k)+B);
  /* And again as we go to +infinity */
  find_exponential_params(f_max,1e12,df,4096,&A,&B,&k);
  for (i=1;i<=4096;i++) world->x_fineparts[4096+16383+i] = A*exp(i*k)+B;
  dfx = df;

  /* Same thing for y as we just did for x */
  f_min = world->bb_llf.y - dfy;
  f_max = world->bb_urb.y + dfy;
  if (f_max - f_min < smallest_spacing)
  {
    printf("Rescaling: was %.3f to %.3f, now ",f_min*world->length_unit,f_max*world->length_unit);
    f = smallest_spacing - (f_max-f_min);
    f_max += 0.5*f;
    f_min -= 0.5*f;
    printf("%.3f to %.3f\n",f_min*world->length_unit,f_max*world->length_unit);
  }
  part_min.y = f_min;
  part_max.y = f_max; 
  df = (f_max - f_min)/16383.0;
  for (i=0;i<16384;i++)
  {
    world->y_fineparts[ 4096 + i ] = f_min + df*((double)i);
  }
  find_exponential_params(-f_min,1e12,df,4096,&A,&B,&k);
  for (i=1;i<=4096;i++) world->y_fineparts[4096-i] = -(A*exp(i*k)+B);
  find_exponential_params(f_max,1e12,df,4096,&A,&B,&k);
  for (i=1;i<=4096;i++) world->y_fineparts[4096+16383+i] = A*exp(i*k)+B;
  dfy = df;

  /* And same again for z */
  f_min = world->bb_llf.z - dfz;
  f_max = world->bb_urb.z + dfz;
  if (f_max - f_min < smallest_spacing)
  {
    printf("Rescaling: was %.3f to %.3f, now ",f_min*world->length_unit,f_max*world->length_unit);
    f = smallest_spacing - (f_max-f_min);
    f_max += 0.5*f;
    f_min -= 0.5*f;
    printf("%.3f to %.3f\n",f_min*world->length_unit,f_max*world->length_unit);
  }
  part_min.z = f_min;
  part_max.z = f_max;
  df = (f_max - f_min)/16383.0;
  for (i=0;i<16384;i++)
  {
    world->z_fineparts[ 4096 + i ] = f_min + df*((double)i);
  }
  find_exponential_params(-f_min,1e12,df,4096,&A,&B,&k);
  for (i=1;i<=4096;i++) world->z_fineparts[4096-i] = -(A*exp(i*k)+B);
  find_exponential_params(f_max,1e12,df,4096,&A,&B,&k);
  for (i=1;i<=4096;i++) world->z_fineparts[4096+16383+i] = A*exp(i*k)+B;
  dfz = df;
  
  /* Try to figure out how many timesteps our fastest particle can make in the whole world (along longest and shortest axes) */
  f = part_max.x - part_min.x;
  f_min = f_max = f;
  f = part_max.y - part_min.y;
  if (f < f_min) f_min = f;
  else if (f > f_max) f_max = f;
  f = part_max.z - part_min.z;
  if (f < f_min) f_min = f;
  else if (f > f_max) f_max = f;
  
  if (world->speed_limit == 0)
  {
    steps_min = f_min;
    steps_max = f_max;
  }
  else
  {
    steps_min = f_min / world->speed_limit;
    steps_max = f_max / world->speed_limit;
  }
 
  /* Verify that partitions are not closer than interaction diameter. */
  if (world->x_partitions!=NULL)
  {
    for (i=1;i<world->nx_parts;i++)
    {
      if (world->x_partitions[i] - world->x_partitions[i-1] < 2*world->rx_radius_3d)
      {
        fprintf(world->err_file,"Error: X partitions closer than interaction diameter\n");
        fprintf(world->err_file,"  X partition #%d at %g\n",i,world->length_unit*world->x_partitions[i-1]);
        fprintf(world->err_file,"  X partition #%d at %g\n",i+1,world->length_unit*world->x_partitions[i]);
        fprintf(world->err_file,"  Interaction diameter %g\n",2*world->length_unit*world->rx_radius_3d);
        return 1;
      }
    }
  }
  if (world->y_partitions!=NULL)
  {
    for (i=1;i<world->ny_parts;i++)
    {
      if (world->y_partitions[i] - world->y_partitions[i-1] < 2*world->rx_radius_3d)
      {
        fprintf(world->err_file,"Error: Y partitions closer than interaction diameter\n");
        fprintf(world->err_file,"  Y partition #%d at %g\n",i,world->length_unit*world->y_partitions[i-1]);
        fprintf(world->err_file,"  Y partition #%d at %g\n",i+1,world->length_unit*world->y_partitions[i]);
        fprintf(world->err_file,"  Interaction diameter %g\n",2*world->length_unit*world->rx_radius_3d);
        return 1;
      }
    }
  }
  if (world->z_partitions!=NULL)
  {
    for (i=1;i<world->nz_parts;i++)
    {
      if (world->z_partitions[i] - world->z_partitions[i-1] < 2*world->rx_radius_3d)
      {
        fprintf(world->err_file,"Error: Z partitions closer than interaction diameter\n");
        fprintf(world->err_file,"  Z partition #%d at %g\n",i,world->length_unit*world->z_partitions[i-1]);
        fprintf(world->err_file,"  Z partition #%d at %g\n",i+1,world->length_unit*world->z_partitions[i]);
        fprintf(world->err_file,"  Interaction diameter %g\n",2*world->length_unit*world->rx_radius_3d);
        return 1;
      }
    }
  }
  

  /* Use automatic partitioning if some of the partitions are not set */
  /* FIXME: should probably be an error, not just a warning, if only some are set */
  if (world->x_partitions == NULL ||
      world->y_partitions == NULL ||
      world->z_partitions == NULL)
  {
    if (world->x_partitions!=NULL || world->y_partitions!=NULL || world->z_partitions!=NULL)
    {
      fprintf(world->err_file,"Warning: some but not all axes are manually partitioned.\n  Using automatic partitions instead--no manual partitions used.\n");
    }
    /* Guess how big to make partitions--nothing really clever about what's done here */
    if (steps_max / MAX_TARGET_TIMESTEP > MAX_COARSE_PER_AXIS)
    {
      world->nx_parts = world->ny_parts = world->nz_parts = MAX_COARSE_PER_AXIS;
    }
    else if (steps_min / MIN_TARGET_TIMESTEP < MIN_COARSE_PER_AXIS)
    {
      world->nx_parts = world->ny_parts = world->nz_parts = MIN_COARSE_PER_AXIS;
    }
    else
    {
      world->nx_parts = steps_min / MIN_TARGET_TIMESTEP;
      if (world->nx_parts > MAX_COARSE_PER_AXIS)
        world->nx_parts = MAX_COARSE_PER_AXIS;
      if ((world->nx_parts & 1) != 0) world->nx_parts += 1;
      
      world->ny_parts = world->nz_parts = world->nx_parts;
    }
    
    /* Allocate memory for our automatically created partitions */
    world->x_partitions = (double*) malloc( sizeof(double) * world->nx_parts );
    world->y_partitions = (double*) malloc( sizeof(double) * world->ny_parts );
    world->z_partitions = (double*) malloc( sizeof(double) * world->nz_parts );
  
    if((world->x_partitions == NULL) || (world->y_partitions == NULL) ||
        (world->z_partitions == NULL))
    {
      fprintf(world->err_file, "File '%s', Line %ld: Out of memory while trying to create partitions.\n", __FILE__, (long)__LINE__);
      return 1;
    }

    /* Calculate aspect ratios so that subvolumes are approximately cubic */
    x_aspect = (part_max.x - part_min.x) / f_max;
    y_aspect = (part_max.y - part_min.y) / f_max;
    z_aspect = (part_max.z - part_min.z) / f_max;
    
    x_in = floor( (world->nx_parts - 2) * x_aspect + 0.5 );
    y_in = floor( (world->ny_parts - 2) * y_aspect + 0.5 );
    z_in = floor( (world->nz_parts - 2) * z_aspect + 0.5 );
    if (x_in < 2) x_in = 2;
    if (y_in < 2) y_in = 2;
    if (z_in < 2) z_in = 2;
    
    /* If we've violated our 2*reaction radius criterion, fix it */
    smallest_spacing = 2*world->rx_radius_3d;
    if ( (part_max.x-part_min.x)/(x_in-1) < smallest_spacing )
    {
      x_in = 1 + floor((part_max.x-part_min.x)/smallest_spacing);
    }
    if ( (part_max.y-part_min.y)/(y_in-1) < smallest_spacing )
    {
      y_in = 1 + floor((part_max.y-part_min.y)/smallest_spacing);
    }
    if ( (part_max.z-part_min.z)/(z_in-1) < smallest_spacing )
    {
      z_in = 1 + floor((part_max.z-part_min.z)/smallest_spacing);
    }
    
    /* Set up to walk symmetrically out from the center of the world, dropping partitions on the way */
    if (x_in < 2) x_in = 2;
    if (y_in < 2) y_in = 2;
    if (z_in < 2) z_in = 2;
    x_start = (world->nx_parts - x_in)/2;
    y_start = (world->ny_parts - y_in)/2;
    z_start = (world->nz_parts - z_in)/2;
    if (x_start < 1) x_start = 1;
    if (y_start < 1) y_start = 1;
    if (z_start < 1) z_start = 1;
    
    /* Now go through and drop partitions in each direction (picked from sensibly close fine partitions) */
    f = (part_max.x - part_min.x) / (x_in - 1);
    world->x_partitions[0] = world->x_fineparts[1];
    /* Dunno how this actually works! */
    for (i=x_start;i<x_start+x_in;i++)
    {
      world->x_partitions[i] = world->x_fineparts[4096 + (i-x_start)*16384/(x_in-1)];
    }
    for (i=x_start-1;i>0;i--)
    {
      for (j=0 ; world->x_partitions[i+1]-world->x_fineparts[4095-j] < f ; j++) {}
      world->x_partitions[i] = world->x_fineparts[4095-j];
    }
    for (i=x_start+x_in;i<world->nx_parts-1;i++)
    {
      for (j=0 ; world->x_fineparts[4096+16384+j]-world->x_partitions[i-1] < f ; j++) {}
      world->x_partitions[i] = world->x_fineparts[4096+16384+j];
    }
    world->x_partitions[world->nx_parts-1] = world->x_fineparts[4096+16384+4096-2];
    
    /* Same thing for y axis */
    f = (part_max.y - part_min.y) / (y_in - 1);
    world->y_partitions[0] = world->y_fineparts[1];
    for (i=y_start;i<y_start+y_in;i++)
    {
      world->y_partitions[i] = world->y_fineparts[4096 + (i-y_start)*16384/(y_in-1)];
    }
    for (i=y_start-1;i>0;i--)
    {
      for (j=0 ; world->y_partitions[i+1]-world->y_fineparts[4095-j] < f ; j++) {}
	world->y_partitions[i] = world->y_fineparts[4095-j];
    }
    for (i=y_start+y_in;i<world->ny_parts-1;i++)
    {
      for (j=0 ; world->y_fineparts[4096+16384+j]-world->y_partitions[i-1] < f ; j++) {}
      world->y_partitions[i] = world->y_fineparts[4096+16384+j];
    }
    world->y_partitions[world->ny_parts-1] = world->y_fineparts[4096+16384+4096-2];
    
    /* Again for z axis */
    f = (part_max.z - part_min.z) / (z_in - 1);
    world->z_partitions[0] = world->z_fineparts[1];
    for (i=z_start;i<z_start+z_in;i++)
    {
      world->z_partitions[i] = world->z_fineparts[4096 + (i-z_start)*16384/(z_in-1)];
    }
    for (i=z_start-1;i>0;i--)
    {
      for (j=0 ; world->z_partitions[i+1]-world->z_fineparts[4095-j] < f ; j++) {}
      world->z_partitions[i] = world->z_fineparts[4095-j];
    }
    for (i=z_start+z_in;i<world->nz_parts-1;i++)
    {
      for (j=0 ; world->z_fineparts[4096+16384+j]-world->z_partitions[i-1] < f ; j++) {}
      world->z_partitions[i] = world->z_fineparts[4096+16384+j];
    }
    world->z_partitions[world->nz_parts-1] = world->z_fineparts[4096+16384+4096-2];

  }
  else /* User-supplied partitions */
  {

    double *dbl_array;

    /* We need to keep the outermost partition away from the world bounding box */
    /* We do this by adding a larger outermost partition, calculated somehow or other */
    dfx += 1e-3;
    dfy += 1e-3;
    dfz += 1e-3;

    /* All this code just adds extra outermost partitions if they might be too close to the outermost objects in the world */
    /* Don't ask me how it actually does it (or if it does it successfully....) */    
    if (world->x_partitions[1] + dfx > world->bb_llf.x)
    {
      if (world->x_partitions[1] - dfx < world->bb_llf.x) world->x_partitions[1] = world->bb_llf.x-dfx;
      else
      {
	dbl_array = (double*) malloc( sizeof(double)*(world->nx_parts+1) );
	if (dbl_array == NULL)
	{ 
	  fprintf(world->err_file, "File '%s', Line %ld: Out of memory while trying to create partitions.\n", __FILE__, (long)__LINE__);
	  return 1;
	}
  
	dbl_array[0] = world->x_partitions[0];
	dbl_array[1] = world->bb_llf.x - dfx;
	memcpy(&(dbl_array[2]),&(world->x_partitions[1]),sizeof(double)*(world->nx_parts-1));
	free( world->x_partitions );
	world->x_partitions = dbl_array;
	world->nx_parts++;
      }
    }
    if (world->x_partitions[world->nx_parts-2] - dfx < world->bb_urb.x)
    {
      if (world->x_partitions[world->nx_parts-2] + dfx > world->bb_urb.x)
	world->x_partitions[world->nx_parts-2] = world->bb_urb.x + dfx;
      else
      {
	dbl_array = (double*) malloc( sizeof(double)*(world->nx_parts+1) );
	if (dbl_array == NULL)
	{ 
	  fprintf(world->err_file, "File '%s', Line %ld: Out of memory while trying to create partitions.\n", __FILE__, (long)__LINE__);
	  return 1;
	}
  
	dbl_array[world->nx_parts] = world->x_partitions[world->nx_parts-1];
	dbl_array[world->nx_parts-1] = world->bb_urb.x + dfx;
	memcpy(dbl_array,world->x_partitions,sizeof(double)*(world->nx_parts-1));
	free( world->x_partitions );
	world->x_partitions = dbl_array;
	world->nx_parts++;
	}
    }
     if (world->y_partitions[1] + dfy > world->bb_llf.y)
    {
      if (world->y_partitions[1] - dfy < world->bb_llf.y)
	world->y_partitions[1] = world->bb_llf.y-dfy;
      else
      {
	dbl_array = (double*) malloc( sizeof(double)*(world->ny_parts+1) );
	if (dbl_array==NULL)
	{ 
	  fprintf(world->err_file, "File '%s', Line %ld: Out of memory while trying to create partitions.\n", __FILE__, (long)__LINE__);
	  return 1;
	}
  
	dbl_array[0] = world->y_partitions[0];
	dbl_array[1] = world->bb_llf.y - dfy;
	memcpy(&(dbl_array[2]),&(world->y_partitions[1]),sizeof(double)*(world->ny_parts-1));
	free( world->y_partitions );
	world->y_partitions = dbl_array;
	world->ny_parts++;
      }
    }
    if (world->y_partitions[world->ny_parts-2] - dfy < world->bb_urb.y)
    {
      if (world->y_partitions[world->ny_parts-2] + dfy > world->bb_urb.y)
	world->y_partitions[world->ny_parts-2] = world->bb_urb.y + dfy;
      else
      {
	dbl_array = (double*) malloc( sizeof(double)*(world->ny_parts+1) );
	if (dbl_array==NULL)
	{
	  fprintf(world->err_file, "File '%s', Line %ld: Out of memory while trying to create partitions.\n", __FILE__, (long)__LINE__);
	  return 1;
	}
  
	dbl_array[world->ny_parts] = world->y_partitions[world->ny_parts-1];
	dbl_array[world->ny_parts-1] = world->bb_urb.y + dfy;
	memcpy(dbl_array,world->y_partitions,sizeof(double)*(world->ny_parts-1));
	free( world->y_partitions );
	world->y_partitions = dbl_array;
	world->ny_parts++;
      }
    }
    if (world->z_partitions[1] + dfz > world->bb_llf.z)
    {
      if (world->z_partitions[1] - dfz < world->bb_llf.z)
	world->z_partitions[1] = world->bb_llf.z-dfz;
      else
      {
	dbl_array = (double*) malloc( sizeof(double)*(world->nz_parts+1) );
	if (dbl_array==NULL)
	{
	  fprintf(world->err_file, "File '%s', Line %ld: Out of memory while trying to create partitions.\n", __FILE__, (long)__LINE__);
	  return 1;
	} 
  
	dbl_array[0] = world->z_partitions[0];
	dbl_array[1] = world->bb_llf.z - dfz;
	memcpy(&(dbl_array[2]),&(world->z_partitions[1]),sizeof(double)*(world->nz_parts-1));
	free( world->z_partitions );
	world->z_partitions = dbl_array;
	world->nz_parts++;
      }
    }
    if (world->z_partitions[world->nz_parts-2] - dfz < world->bb_urb.z)
    {
      if (world->z_partitions[world->nz_parts-2] + dfz > world->bb_urb.z)
	world->z_partitions[world->nz_parts-2] = world->bb_urb.z + dfz;
      else
      {
	dbl_array = (double*) malloc( sizeof(double)*(world->nz_parts+1) );
	if (dbl_array==NULL){
	  fprintf(world->err_file, "File '%s', Line %ld: Out of memory while trying to create partitions.\n", __FILE__, (long)__LINE__);
	  return 1;
	} 
  
	dbl_array[world->nz_parts] = world->z_partitions[world->nz_parts-1];
	dbl_array[world->nz_parts-1] = world->bb_urb.z + dfz;
	memcpy(dbl_array,world->z_partitions,sizeof(double)*(world->nz_parts-1));
	free( world->z_partitions );
	world->z_partitions = dbl_array;
	world->nz_parts++;
      }
    }
   
    /* Now that we've added outermost partitions, we find the closest fine partition along each axis */
    world->x_partitions[0] = world->x_fineparts[1];
    for (i=1;i<world->nx_parts-1;i++)
    {
      world->x_partitions[i] = 
        world->x_fineparts[ 
          bisect_near( 
            world->x_fineparts , world->n_fineparts ,
            world->x_partitions[i]
          )
        ];
    }
    world->x_partitions[world->nx_parts-1] = world->x_fineparts[4096+16384+4096-2];

    world->y_partitions[0] = world->y_fineparts[1];
    for (i=1;i<world->ny_parts-1;i++)
    {
      world->y_partitions[i] = 
        world->y_fineparts[ 
          bisect_near( 
            world->y_fineparts , world->n_fineparts ,
            world->y_partitions[i]
          )
        ];
    }
    world->y_partitions[world->ny_parts-1] = world->y_fineparts[4096+16384+4096-2];

    world->z_partitions[0] = world->z_fineparts[1];
    for (i=1;i<world->nz_parts-1;i++)
    {
      world->z_partitions[i] = 
        world->z_fineparts[ 
          bisect_near( 
            world->z_fineparts , world->n_fineparts ,
            world->z_partitions[i]
          )
        ];
    }
    world->z_partitions[world->nz_parts-1] = world->z_fineparts[4096+16384+4096-2];
  }
  
  /* And finally we tell the user what happened */
  if (world->notify->partition_location==NOTIFY_FULL)
  {
    printf("X partitions: ");
    printf("-inf ");
    for (i=1;i<world->nx_parts - 1;i++) printf("%.5f ",world->length_unit * world->x_partitions[i]);
    printf("inf");
    printf("\n");
    printf("Y partitions: ");
    printf("-inf ");
    for (i=1;i<world->ny_parts - 1;i++) printf("%.5f ",world->length_unit * world->y_partitions[i]);
    printf("inf");
    printf("\n");
    printf("Z partitions: ");
    printf("-inf ");
    for (i=1;i<world->nz_parts - 1;i++) printf("%.5f ",world->length_unit * world->z_partitions[i]);
    printf("inf");
    printf("\n");
  }

  return 0;
}
/**********************************************************************
distance_point_line -- returns distance between point and line in 3D space
              The line is defined by 2 points
              The formulas are taken from "Computer Graphics" 
              by Michael Mortenson, ISBN 0-8311-1182-8, p.222
Parameters
	q -- location of the fixed point
	v0 -- first point on the line
	v1 -- second point on the line

Returns
	distance between the point and the line
**********************************************************************/
double distance_point_line(struct vector3 *q, struct vector3 *v0, struct vector3 *v1)
{
   double x,y,z; /* coordinates of the point q */
   double x0,y0,z0; /* coordinates of the point v0 */
   double x1,y1,z1; /* coordinates of the point v1 */
   double u; /* parameter in the equation of the line
                p(u) = p0 + u(p1 - p0) */
   double p_x,p_y,p_z; /* X,Y, and Z components of the vector p (line) */
   double nominator, denominator, d_min;

   x = q->x;
   y = q->y;
   z = q->z;

   x0 = v0->x;
   y0 = v0->y;
   z0 = v0->z;

   x1 = v1->x;
   y1 = v1->y;
   z1 = v1->z;
   
   nominator = (x1 - x0)*(x - x0) + (y1 - y0)*(y - y0) + (z1 - z0)*(z - z0);
   denominator = sqrt(pow((x1-x0),2) + pow((y1 - y0),2) + pow((z1 - z0),2));
   u = nominator/denominator;

   p_x = x0 +u*(x1 - x0);
   p_y = y0 +u*(y1 - y0);
   p_z = z0 +u*(z1 - z0);

   d_min = sqrt(pow((p_x - x),2) + pow((p_y - y),2) + pow((p_z - z),2)); 
   return d_min;
} 

/**********************************************************************
navigate_world -- returns index of the destination subvolume

Parameters
	cur_index -- index of the starting subvolume 
	direction -- direction of the movement
                     Valid values are X_NEG, X_POS, Y_NEG, Y_POS,
                                      Z_NEG, Z_POS.

Returns
	index of the destination neighbor subvolume ("face-to-face")
**********************************************************************/
int navigate_world(int curr_index, int direction)
{
        switch(direction)
        {
		case(X_NEG):
		   return curr_index - (world->nz_parts-1)*(world->ny_parts-1);
                   break;
                case(X_POS):
		   return curr_index + (world->nz_parts-1)*(world->ny_parts-1);
                   break;
                case(Y_NEG):  
		   return curr_index - (world->nz_parts-1);
                   break;
                case(Y_POS):
		   return curr_index + (world->nz_parts-1);
                   break;
                case(Z_NEG):
                   return curr_index - 1;
                   break;
                case(Z_POS):
                   return curr_index + 1;
                   break;
                default:
		   return INT_MAX;
                   break;

        }
 
}

/**********************************************************************
navigate_world_by_edge -- returns index of the destination subvolume

Parameters
	cur_index -- index of the starting subvolume 
	direction1 -- direction of the first movement
	direction2 -- direction of the second movement
                     Valid values are X_NEG, X_POS, Y_NEG, Y_POS,
                                      Z_NEG, Z_POS.

Returns
	index of the destination neighbor subvolume ("edge-to-edge")
**********************************************************************/
int navigate_world_by_edge(int curr_index, int direction1, int direction2)
{
  	int first_index, final_index;
        first_index = navigate_world(curr_index, direction1);
        final_index = navigate_world(first_index, direction2);
    
        return final_index; 
}
/**********************************************************************
navigate_world_by_corner -- returns index of the destination subvolume

Parameters
	cur_index -- index of the starting subvolume 
	direction1 -- direction of the first movement
	direction2 -- direction of the second movement
	direction3 -- direction of the third movement
                     Valid values are X_NEG, X_POS, Y_NEG, Y_POS,
                                      Z_NEG, Z_POS.

Returns
	index of the destination neighbor subvolume ("corner-to-corner")
**********************************************************************/
int navigate_world_by_corner(int curr_index, int direction1, int direction2, int direction3)
{
  	int first_index, second_index, final_index;
        first_index = navigate_world(curr_index, direction1);
        second_index = navigate_world(first_index, direction2);
        final_index = navigate_world(second_index, direction3);
    
        return final_index; 
}

/************************************************************************
   In: starting position of the molecule
       displacement (random walk) vector
       vector to store one corner of the bounding box
       vector to store the opposite corner of the bounding box 
   Out: No return value. The vectors are set to define the bounding box
        of the random walk movement that extends for R_INT in all
        directions.
************************************************************************/
void path_bounding_box(struct vector3 *loc, struct vector3 * displacement,
 struct vector3 *llf, struct vector3 *urb)
{
   struct vector3 final;  /* final position of the molecule after random walk */
   double R;     /* molecule interaction radius */

  
   R = world->rx_radius_3d; 
   vect_sum(loc, displacement, &final);

   llf->x = urb->x = loc->x;
   llf->y = urb->y = loc->y;
   llf->z = urb->z = loc->z;

   if(final.x < llf->x) {
         llf->x = final.x;
   }
   if(final.x > urb->x){
       urb->x = final.x;
   }
   if(final.y < llf->y) {
         llf->y = final.y;
   }
   if(final.y > urb->y){
       urb->y = final.y;
   }
   if(final.z < llf->z) {
         llf->z = final.z;
   }
   if(final.z > urb->z){
       urb->z = final.z;
   }
   /* Extend the bounding box at the distance R. */
   llf->x -= R;
   llf->y -= R;
   llf->z -= R;

   urb->x += R;
   urb->y += R;
   urb->z += R;
}

