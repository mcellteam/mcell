/***********************************************************************************
 *                                                                                 *
 * Copyright (C) 2006-2013 by                                                      *
 * The Salk Institute for Biological Studies and                                   *
 * Pittsburgh Supercomputing Center, Carnegie Mellon University                    *
 *                                                                                 *
 * This program is free software; you can redistribute it and/or                   *
 * modify it under the terms of the GNU General Public License                     *
 * as published by the Free Software Foundation; either version 2                  *
 * of the License, or (at your option) any later version.                          *
 *                                                                                 *
 * This program is distributed in the hope that it will be useful,                 *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of                  *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                   *
 * GNU General Public License for more details.                                    *
 *                                                                                 *
 * You should have received a copy of the GNU General Public License               *
 * along with this program; if not, write to the Free Software                     *
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA. *
 *                                                                                 *
 ***********************************************************************************/

/**************************************************************************\
** File: vol_util.c                                                       **
**                                                                        **
** Purpose: Adds, subtracts, and moves particles around (bookkeeping).    **
**                                                                        **
** Testing status: compiles.  Worked earlier, but has been changed.       **
\**************************************************************************/

#include "config.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "logging.h"
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
#include "macromolecule.h"

static int release_inside_regions(struct volume *world, 
    struct release_site_obj *rso, struct volume_molecule *m,
    int n);

/*************************************************************************
inside_subvolume:
  In: pointer to vector3
      pointer to subvolume
  Out: nonzero if the vector is inside the subvolume.
*************************************************************************/
int 
inside_subvolume(struct vector3 *point, struct subvolume *subvol,
    double *x_fineparts, double *y_fineparts, double *z_fineparts)
{
  return ( (point->x >= x_fineparts[ subvol->llf.x ] ) &&
           (point->x <= x_fineparts[ subvol->urb.x ] ) &&
           (point->y >= y_fineparts[ subvol->llf.y ] ) &&
           (point->y <= y_fineparts[ subvol->urb.y ] ) &&
           (point->z >= z_fineparts[ subvol->llf.z ] ) &&
           (point->z <= z_fineparts[ subvol->urb.z ] ) );
}


/*************************************************************************
find_coarse_subvolume:
  In: pointer to vector3
  Out: pointer to the coarse subvolume that the vector is within
*************************************************************************/
struct subvolume* 
find_coarse_subvol(struct volume *world, struct vector3 *loc)
{
  int i,j,k;
  i = bisect(world->x_partitions,world->nx_parts,loc->x);
  j = bisect(world->y_partitions,world->ny_parts,loc->y);
  k = bisect(world->z_partitions,world->nz_parts,loc->z);
  return
    &(world->subvol
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
struct subvolume* 
traverse_subvol(struct subvolume *here, struct vector3 *point,
    int which, int nx_part, int ny_parts, int nz_parts)
{
  UNUSED(point);
    switch(which)
    {
      case X_NEG:
          if (here->world_edge&X_NEG_BIT) return NULL;
          return here - (nz_parts - 1)*(ny_parts - 1);
      case X_POS:
          if (here->world_edge&X_POS_BIT) return NULL;
          return here + (nz_parts - 1)*(ny_parts - 1);
      case Y_NEG:
          if (here->world_edge&Y_NEG_BIT) return NULL;
          return here - (nz_parts - 1);
      case Y_POS:
          if (here->world_edge&Y_POS_BIT) return NULL;
          return here + (nz_parts - 1);
      case Z_NEG:
          if (here->world_edge&Z_NEG_BIT) return NULL;
          return here - 1;
      case Z_POS:
          if (here->world_edge&Z_POS_BIT) return NULL;
          return here + 1;
      default:
          mcell_internal_error("Invalid direction specified in traverse_subvol (dir=%d).", which);
          return NULL;
    } /* end switch */

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
      if ((branch->flags & X_AXIS) != 0)
      {
        if (point->x <= world->x_fineparts[ branch->partition ]) left_path = 1;
        else left_path = 0;
      }
      else
      {
        if ((branch->flags & Y_AXIS) != 0)
        {
          if (point->y <= world->y_fineparts[ branch->partition ]) left_path = 1;
          else left_path = 0;
        }
        else // Must be Z_AXIS
        {
          if (point->z <= world->z_fineparts[ branch->partition ]) left_path = 1;
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
double 
collide_sv_time(struct vector3 *here, struct vector3 *move, 
    struct subvolume *sv, double *x_fineparts, double *y_fineparts,
    double *z_fineparts)
{
  double dx,dy,dz,tx,ty,tz,t;

  if (move->x==0 && move->y==0 && move->z==0) return GIGANTIC;

  if (move->x > 0) dx = x_fineparts[ sv->urb.x ] - here->x;
  else { dx = x_fineparts[ sv->llf.x ] - here->x; }

  if (move->y > 0) dy = y_fineparts[ sv->urb.y ] - here->y;
  else { dy = y_fineparts[ sv->llf.y ] - here->y; }

  if (move->z > 0) dz = z_fineparts[ sv->urb.z ] - here->z;
  else { dz = z_fineparts[ sv->llf.z ] - here->z; }

  tx = dx * move->y * move->z; if (tx<0) tx = -tx;
  ty = move->x * dy * move->z; if (ty<0) ty = -ty;
  tz = move->x * move->y * dz; if (tz<0) tz = -tz;

  if (tx<ty || move->y==0.0)
  {
    if (tx<tz || move->z==0.0)
    { t = dx / move->x; } /* Collision with X */
    else                       { t = dz / move->z; } /* Collision with Z */
  }
  else /* ty<tx */
  {
    if (ty<tz || move->z==0.0)
    { t = dy / move->y; } /* Collision with Y */
    else                       { t = dz / move->z; } /* Collision with Z */
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
struct subvolume* 
next_subvol(struct vector3 *here, struct vector3 *move, struct subvolume *sv,
    double *x_fineparts, double *y_fineparts, double *z_fineparts, 
    int nx_parts, int ny_parts, int nz_parts)
{
  double dx,dy,dz,tx,ty,tz,t;
  int whichx,whichy,whichz,which;

  whichx = whichy = whichz = 1;
  if (move->x==0 && move->y==0 && move->z==0) return NULL;

  if (move->x > 0) dx = x_fineparts[ sv->urb.x ] - here->x;
  else { dx = x_fineparts[ sv->llf.x ] - here->x; whichx = 0; }

  if (move->y > 0) dy = y_fineparts[ sv->urb.y ] - here->y;
  else { dy = y_fineparts[ sv->llf.y ] - here->y; whichy = 0; }

  if (move->z > 0) dz = z_fineparts[ sv->urb.z ] - here->z;
  else { dz = z_fineparts[ sv->llf.z ] - here->z; whichz = 0; }

  if (move->x == 0.0)
  {
    ty = dy * move->z; if (ty<0) ty = -ty;
    tz = move->y * dz; if (tz<0) tz = -tz;
    if (ty < tz)
    { t = dy / move->y; which = Y_NEG + whichy; }
    else { t = dz / move->z; which = Z_NEG + whichz; }
  }
  else if (move->y == 0.0)
  {
    tx = dx * move->z; if (tx<0) tx = -tx;
    tz = move->x * dz; if (tz<0) tz = -tz;
    if (tx < tz)
    { t = dx / move->x; which = X_NEG + whichx; }
    else { t = dz / move->z; which = Z_NEG + whichz; }
  }
  else if (move->z == 0.0)
  {
    tx = dx * move->y; if (tx<0) tx = -tx;
    ty = move->x * dy; if (ty<0) ty = -ty;
    if (tx < ty)
    { t = dx / move->x; which = X_NEG + whichx; }
    else { t = dy / move->y; which = Y_NEG + whichy; }
  }
  else
  {
    tx = dx * move->y * move->z; if (tx<0) tx = -tx;
    ty = move->x * dy * move->z; if (ty<0) ty = -ty;
    tz = move->x * move->y * dz; if (tz<0) tz = -tz;

    if (tx<ty)
    {
      if (tx<tz)
      { t = dx / move->x; which = X_NEG + whichx; }
      else { t = dz / move->z; which = Z_NEG + whichz; }
    }
    else /* ty<tx */
    {
      if (ty<tz)
      { t = dy / move->y; which = Y_NEG + whichy; }
      else { t = dz / move->z; which = Z_NEG + whichz; }
    }
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

    return traverse_subvol(sv,here,which, nx_parts, ny_parts, nz_parts);
  }
}



/*************************************************************************
find_subvolume:
  In: pointer to a vector3 of where we are
      pointer to a subvolume we might be in or near
  Out: subvolume that we are in
*************************************************************************/
struct subvolume* 
find_subvolume(struct volume *world, struct vector3 *loc,
    struct subvolume *guess)
{
#if 1
  /* This code is faster if coarse subvolumes are always used */

  if (guess==NULL) return find_coarse_subvol(world, loc);
  else
  {
    if (world->x_fineparts[guess->llf.x] <= loc->x && loc->x <= world->x_fineparts[guess->urb.x] &&
        world->y_fineparts[guess->llf.y] <= loc->y && loc->y <= world->y_fineparts[guess->urb.y] &&
        world->z_fineparts[guess->llf.z] <= loc->z && loc->z <= world->z_fineparts[guess->urb.z])
    {
      return guess;
    }
    else return find_coarse_subvol(world, loc);
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

  while (loc->x < world->x_fineparts[ sv->llf.x ])
  {
    sv = traverse_subvol(sv , &center , X_NEG);
    center.x = 0.5*(world->x_fineparts[ sv->llf.x ] + world->x_fineparts[ sv->urb.x ]);
  }
  while (loc->x > world->x_fineparts[ sv->urb.x ])
  {
    sv = traverse_subvol(sv , &center , X_POS);
    center.x = 0.5*(world->x_fineparts[ sv->llf.x ] + world->x_fineparts[ sv->urb.x ]);
  }
  center.x = loc->x;

  while (loc->y < world->y_fineparts[ sv->llf.y ])
  {
    sv = traverse_subvol(sv , &center , Y_NEG);
    center.y = 0.5*(world->y_fineparts[ sv->llf.y ] + world->y_fineparts[ sv->urb.y ]);
  }
  while (loc->y > world->y_fineparts[ sv->urb.y ])
  {
    sv = traverse_subvol(sv , &center , Y_POS);
    center.y = 0.5*(world->y_fineparts[ sv->llf.y ] + world->y_fineparts[ sv->urb.y ]);
  }
  center.y = loc->y;

  while (loc->z < world->z_fineparts[ sv->llf.z ])
  {
    sv = traverse_subvol(sv , &center , Z_NEG);
    center.z = 0.5*(world->z_fineparts[ sv->llf.z ] + world->z_fineparts[ sv->urb.z ]);
  }
  while (loc->z > world->z_fineparts[ sv->urb.z ])
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
int 
is_defunct_molecule(struct abstract_element *e)
{
  return ((struct abstract_molecule*)e)->properties == NULL;
}



/*************************************************************************
place_grid_molecule
  In: species for the new molecule
      3D location of the new molecule
      orientation of the new molecule
      diameter to search for a free surface spot
      schedule time for the new molecule
  Out: pointer to the new molecule, or NULL if no free spot was found.
  Note: This function halts the program if it runs out of memory.
        This function is similar to insert_grid_molecule, but it does
        not schedule the molecule or add it to the count.  This is done
        to simplify the logic when placing a surface macromolecule.
        (i.e. place all molecules, and once we're sure we've succeeded,
        schedule them all and count them all.)
 *************************************************************************/
struct grid_molecule* 
place_grid_molecule(struct volume *world, struct species *s, 
    struct vector3 *loc, short orient, double search_diam, double t, 
    struct subvolume **psv, struct grid_molecule **cmplx)
{
  double search_d2,d2;
  struct vector2 s_loc;

  double best_d2;
  struct wall *best_w;
  struct vector2 best_uv;
  struct vector3 best_xyz;

  struct subvolume *sv;
  struct wall_list *wl;
  struct grid_molecule *g;

  if (search_diam<=EPS_C) search_d2 = EPS_C*EPS_C;
  else search_d2 = search_diam * search_diam;

  sv = find_subvolume(world, loc, NULL);

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
    const int sv_index = sv - world->subvol;
    int sv_remain = sv_index;

    /* Turn linear sv_index into part_x, part_y, part_z triple. */
    const int part_x = sv_remain / ((world->ny_parts-1)*(world->nz_parts-1));
    sv_remain -= part_x * ((world->ny_parts-1)*(world->nz_parts-1));
    const int part_y = sv_remain / (world->nz_parts-1);
    sv_remain -= part_y * (world->nz_parts-1);
    const int part_z = sv_remain;

    /* Find min x partition. */
    int x_min;
    for (x_min=part_x; x_min>0; x_min--)
    {
      d2 = loc->x - world->x_partitions[x_min]; d2 *= d2;
      if (d2 >= best_d2 || d2 >= search_d2) break;
    }

    /* Find max x partition. */
    int x_max;
    for (x_max=part_x; x_max<world->nx_parts-1 ; x_max++)
    {
      d2 = loc->x - world->x_partitions[x_max + 1]; d2 *= d2;
      if (d2 >= best_d2 || d2 >= search_d2) break;
    }

    /* Find min y partition. */
    int y_min;
    for (y_min=part_y; y_min>0; y_min--)
    {
      d2 = loc->y - world->y_partitions[y_min]; d2 *= d2;
      if (d2 >= best_d2 || d2 >= search_d2) break;
    }

    /* Find max y partition. */
    int y_max;
    for (y_max=part_y; y_max<world->ny_parts-1; y_max++)
    {
      d2 = loc->y - world->y_partitions[y_max+1]; d2 *= d2;
      if (d2 >= best_d2 || d2 >= search_d2) break;
    }

    /* Find min z partition. */
    int z_min;
    for (z_min=part_z; z_min>0; z_min--)
    {
      d2 = loc->z - world->z_partitions[z_min]; d2 *= d2;
      if (d2 >= best_d2 || d2 >= search_d2) break;
    }

    /* Find max z partition. */
    int z_max;
    for (z_max=part_z; z_max<world->nz_parts-1; z_max++)
    {
      d2 = loc->z - world->z_partitions[z_max+1]; d2 *= d2;
      if (d2 >= best_d2 || d2 >= search_d2) break;
    }

    if (x_min<part_x || x_max>part_x || y_min<part_y || y_max>part_y || z_min<part_z || z_max>part_z)
    {
      for (int px=x_min; px<=x_max; px++)
      {
        for (int py=y_min; py<=y_max; py++)
        {
          for (int pz=z_min; pz<=z_max; pz++)
          {
            const int this_sv = pz + (world->nz_parts-1)*(py + (world->ny_parts-1)*px);
            if (this_sv == sv_index) continue;

            for (wl=world->subvol[this_sv].wall_head; wl!=NULL; wl=wl->next)
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
        sv = find_subvolume(world, &best_xyz,sv);  /* May have switched subvolumes */
      }
    }
  }

  if (best_w==NULL)
  {
    return NULL;
  }

  d2 = search_d2 - best_d2;  /* We can look this far around the surface we hit for an empty spot */

  int grid_index;
  if (best_w->grid==NULL)
  {
    if (create_grid(world, best_w,sv))
      mcell_allocfailed("Failed to create grid for wall.");
    grid_index = uv2grid(&best_uv,best_w->grid);
  }
  else if ((s->flags & IS_COMPLEX) != 0)
  {
    grid_index = uv2grid(&best_uv,best_w->grid);
  }
  else
  {
    grid_index = uv2grid(&best_uv,best_w->grid);
    if (best_w->grid->mol[grid_index]!=NULL)
    {
      if (d2 <= EPS_C*EPS_C)
      {
        return NULL;
      }
      else
      {
        best_w = search_nbhd_for_free(world, best_w,&best_uv,d2,
            &grid_index,NULL,NULL);
        if (best_w==NULL)
        {
          return NULL;
        }

        if (world->randomize_gmol_pos) grid2uv_random(best_w->grid,
            grid_index, &best_uv, world->rng);
        else grid2uv(best_w->grid,grid_index,&best_uv);
      }
    }
  }

  uv2xyz(&best_uv, best_w, &best_xyz);
  sv = find_subvolume(world, &best_xyz, sv);

  g = CHECKED_MEM_GET(sv->local_storage->gmol, "grid molecule");
  g->birthplace = sv->local_storage->gmol;
  g->birthday = t;
  g->id = world->current_mol_id++;
  g->properties = s;
  s->population++;
  g->cmplx = cmplx;
  g->flags = TYPE_GRID | ACT_NEWBIE | IN_SCHEDULE;
  if ((s->flags & IS_COMPLEX) != 0)
    g->flags |= COMPLEX_MASTER;
  else if (g->cmplx)
    g->flags |= COMPLEX_MEMBER;
  if (s->space_step > 0) g->flags |= ACT_DIFFUSE;
  if (trigger_unimolecular(world->reaction_hash, world->rx_hashsize,
        s->hashval, (struct abstract_molecule*)g)!= NULL || (s->flags&CAN_GRIDWALL)!=0 ) g->flags |= ACT_REACT;

  g->t = t;
  g->t2 = 0.0;
  g->grid = best_w->grid;
  g->grid_index = grid_index;
  g->s_pos.u = best_uv.u;
  g->s_pos.v = best_uv.v;
  g->orient = orient;

  /* Put it on the grid if it doesn't represent a macromolecular complex */
  if ((s->flags & IS_COMPLEX) == 0)
  {
    g->grid->mol[ g->grid_index ] = g;
    g->grid->n_occupied++;
    g->flags |= IN_SURFACE;
  }

  if ((s->flags&COUNT_ENCLOSED) != 0) g->flags |= COUNT_ME;

  *psv = sv;
  return g;
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
struct grid_molecule* 
insert_grid_molecule(struct volume *world, struct species *s,
    struct vector3 *loc, short orient, double search_diam, double t,
    struct grid_molecule **cmplx)
{
  struct subvolume *sv = NULL;
  struct grid_molecule *g = place_grid_molecule(world, s, loc, orient, 
      search_diam, t, &sv, cmplx);
  if (g == NULL)
    return NULL;

  if (g->properties->flags & (COUNT_CONTENTS|COUNT_ENCLOSED))
    count_region_from_scratch(world, (struct abstract_molecule*)g, NULL, 
        1, NULL, g->grid->surface, g->t);

  if (schedule_add(sv->local_storage->timer,g))
    mcell_allocfailed("Failed to add grid molecule to scheduler.");

  return g;
}



/*************************************************************************
insert_volume_molecule
  In: pointer to a volume_molecule that we're going to place in local storage
      pointer to a volume_molecule that may be nearby
  Out: pointer to the new volume_molecule (copies data from volume molecule
       passed in), or NULL if out of memory.  Molecule is placed in scheduler
       also.
*************************************************************************/
struct volume_molecule* 
insert_volume_molecule(struct volume *world, struct volume_molecule *m,
    struct volume_molecule *guess)
{
  struct volume_molecule *new_m;
  struct subvolume *sv;

  if (guess == NULL) sv = find_subvolume(world, &(m->pos),NULL);
  else if (inside_subvolume(&(m->pos), guess->subvol, world->x_fineparts,
                            world->y_fineparts, world->z_fineparts)) 
    sv = guess->subvol;
  else sv = find_subvolume(world, &(m->pos),guess->subvol);

  new_m = CHECKED_MEM_GET(sv->local_storage->mol, "volume molecule");
  memcpy(new_m,m,sizeof(struct volume_molecule));
  new_m->birthplace = sv->local_storage->mol;
  new_m->id = world->current_mol_id++;
  new_m->prev_v = NULL;
  new_m->next_v = NULL;
  new_m->next = NULL;
  new_m->subvol = sv;
  ht_add_molecule_to_list(&sv->mol_by_species, new_m);
  sv->mol_count++;
  new_m->properties->population++;

  if ((new_m->properties->flags&COUNT_SOME_MASK) != 0) new_m->flags |= COUNT_ME;
  if (new_m->properties->flags & (COUNT_CONTENTS|COUNT_ENCLOSED))
  {
    count_region_from_scratch(world, (struct abstract_molecule*)new_m, NULL, 
        1, &(new_m->pos), NULL, new_m->t);
  }

  if (schedule_add(sv->local_storage->timer,new_m))
    mcell_allocfailed("Failed to add volume molecule to scheduler.");
  return new_m;
}



/*************************************************************************
exsert_volume_molecule:
  In: pointer to a volume_molecule that we're going to remove from local storage
  Out: no return value; molecule is marked for removal.
*************************************************************************/
void 
exsert_volume_molecule(struct volume *world, struct volume_molecule *m)
{
  if (m->properties->flags & (COUNT_CONTENTS|COUNT_ENCLOSED))
  {
    count_region_from_scratch(world, (struct abstract_molecule*)m, NULL,  
        -1, NULL, NULL, m->t);
  }
  m->subvol->mol_count--;
  m->properties->n_deceased++;
  m->properties->cum_lifetime += m->t - m->birthday;
  m->properties->population--;
  collect_molecule(m);
}



/*************************************************************************
insert_volume_molecule_list:
  In: pointer to a linked list of volume_molecules to copy into subvolumes.
  Out: 0 on success, 1 on memory allocation error; molecules are placed
       in their subvolumes.
*************************************************************************/
int 
insert_volume_molecule_list(struct volume *world, struct volume_molecule *m)
{
  struct volume_molecule *new_m,*guess;

  guess=NULL;
  while (m != NULL)
  {
    new_m = insert_volume_molecule(world, m, guess);
    if (new_m == NULL)
      mcell_allocfailed("Failed to add volume molecule to world.");
    guess = new_m;
    m = (struct volume_molecule*)m->next;
  }

  return 0;
}


static int remove_from_list(struct volume_molecule *it)
{
  if (it->prev_v)
  {
#ifdef DEBUG_LIST_CHECKS
    if (*it->prev_v != it)
    {
      mcell_error_nodie("Stale previous pointer!");
    }
#endif
  *(it->prev_v) = it->next_v;
  }
  else
  {
#ifdef DEBUG_LIST_CHECKS
    mcell_error_nodie("No previous pointer.");
#endif
  }
  if (it->next_v)
  {
#ifdef DEBUG_LIST_CHECKS
    if (it->next_v->prev_v != &it->next_v)
    {
      mcell_error_nodie("Stale next pointer!");
    }
#endif
    it->next_v->prev_v = it->prev_v;
  }
  it->prev_v = NULL;
  it->next_v = NULL;
  return 1;
}

/*************************************************************************
migrate_volume_molecule:
  In: pointer to a volume_molecule already in a subvolume
      pointer to the new subvolume to move it to
  Out: pointer to moved molecule.  The molecule's position is updated
       but it is not rescheduled.  Returns NULL if out of memory.
*************************************************************************/
struct volume_molecule* 
migrate_volume_molecule(struct volume_molecule *m, struct subvolume *new_sv)
{
  struct volume_molecule *new_m;

  new_sv->mol_count++;
  m->subvol->mol_count--;

  if (m->subvol->local_storage == new_sv->local_storage)
  {
    if (remove_from_list(m))
    {
      m->subvol = new_sv;
      ht_add_molecule_to_list(&new_sv->mol_by_species, m);
      return m;
    }
  }

  new_m = CHECKED_MEM_GET(new_sv->local_storage->mol, "volume molecule");
  memcpy(new_m,m,sizeof(struct volume_molecule));
  new_m->birthplace = new_sv->local_storage->mol;
  new_m->prev_v = NULL;
  new_m->next_v = NULL;
  new_m->next = NULL;
  new_m->subvol = new_sv;

  ht_add_molecule_to_list(&new_sv->mol_by_species, new_m);

  collect_molecule(m);

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
int 
eval_rel_region_3d(struct release_evaluator *expr, struct waypoint *wp,
    struct region_list *in_regions, struct region_list *out_regions)
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
  else satisfies_r = eval_rel_region_3d(expr->right, wp,
      in_regions, out_regions);

  if (expr->op & REXP_UNION) return (satisfies_l || satisfies_r);
  else if (expr->op & (REXP_INTERSECTION|REXP_INCLUSION)) 
    return (satisfies_l && satisfies_r);
  else if (expr->op & REXP_SUBTRACTION) 
    return (satisfies_l && !satisfies_r);

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
static int 
vacuum_inside_regions(struct volume *world, struct release_site_obj *rso,
    struct volume_molecule *m, int n)
{
  struct volume_molecule *mp;
  struct release_region_data *rrd;
  struct region_list *extra_in,*extra_out;
  struct region_list *rl,*rl2;
  struct waypoint *wp;
  struct subvolume *sv = NULL;
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

  const int x_min = bisect(world->x_partitions,world->nx_parts,rrd->llf.x);
  const int x_max = bisect_high(world->x_partitions,world->nx_parts,rrd->urb.x);
  const int y_min = bisect(world->y_partitions,world->ny_parts,rrd->llf.y);
  const int y_max = bisect_high(world->y_partitions,world->ny_parts,rrd->urb.y);
  const int z_min = bisect(world->z_partitions,world->nz_parts,rrd->llf.z);
  const int z_max = bisect_high(world->z_partitions,world->nz_parts,rrd->urb.z);

  for (int px=x_min; px<x_max; px++)
  {
    for (int py=y_min; py<y_max; py++)
    {
      for (int pz=z_min; pz<z_max; pz++)
      {
        const int this_sv = pz + (world->nz_parts - 1)*(py + (world->ny_parts-1)*px);
        sv = &(world->subvol[this_sv]);

        struct per_species_list *psl = (struct per_species_list *) pointer_hash_lookup(&sv->mol_by_species, m->properties, m->properties->hashval);
        if (psl != NULL)
        {
          for (mp = psl->head; mp != NULL; mp = mp->next_v)
          {
            extra_in=extra_out=NULL;
            wp = &(world->waypoints[this_sv]);
            origin = &(wp->loc);
            delta.x = mp->pos.x - origin->x;
            delta.y = mp->pos.y - origin->y;
            delta.z = mp->pos.z - origin->z;

            for (wl=sv->wall_head ; wl!=NULL ; wl=wl->next)
            {
              int hitcode = collide_wall(origin,&delta,wl->this_wall,&t,
                  &hit, 0, world->rng, world->notify, 
                  &(world->ray_polygon_tests));
              if (hitcode != COLLIDE_MISS)
              {
                world->ray_polygon_colls++;

                for (rl=wl->this_wall->counting_regions ; rl!=NULL ; rl=rl->next)
                {
                  if (hitcode == COLLIDE_FRONT || hitcode == COLLIDE_BACK)
                  {
                    rl2 = (struct region_list*) CHECKED_MEM_GET(sv->local_storage->regl, "region list");
                    rl2->reg = rl->reg;

                    if (hitcode == COLLIDE_FRONT)
                    {
                      rl2->next = extra_in;
                      extra_in = rl2;
                    }
                    else  /*hitcode == COLLIDE_BACK*/
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

            if (eval_rel_region_3d(rrd->expression,wp,extra_in,extra_out))
            {
              vl = (struct void_list *) CHECKED_MEM_GET(mh, "temporary list");
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

  for (vl=vl_head ; n<0 && vl_num>0 && vl!=NULL ; vl=vl->next , vl_num--)
  {
    if (rng_dbl(world->rng) < ((double)(-n))/((double)vl_num))
    {
      mp = (struct volume_molecule*)vl->data;
      mp->properties->population--;
      mp->subvol->mol_count--;
      if ((mp->properties->flags & (COUNT_CONTENTS|COUNT_ENCLOSED)) != 0)
        count_region_from_scratch(world, (struct abstract_molecule*)mp, 
            NULL, -1, &(mp->pos), NULL, mp->t);
      if (mp->flags & IN_SCHEDULE)
      {
        mp->subvol->local_storage->timer->defunct_count++; /* Tally for garbage collection */
      }
      collect_molecule(mp);

      n++;
    }
  }

  delete_mem(mh);
  return 0;
}


/*************************************************************************
 is_point_inside_release_region:
    Check if a given point is inside the specified region.

*************************************************************************/
static int 
is_point_inside_region(struct volume *world, struct vector3 const *pos,
    struct release_evaluator *expression, struct subvolume *sv)
{
  struct region_list *extra_in=NULL, *extra_out=NULL, *cur_region;
  struct waypoint *wp;
  struct vector3 delta;
  struct vector3 *origin;
  struct wall_list *wl;
  int bad_location = 0;
  int result;

  /* If no subvolume hint was given, or the hint is incorrect, find the right
   * subvolume
   */
  if (sv == NULL || 
      ! inside_subvolume((struct vector3 *) pos, sv, world->x_fineparts,
        world->y_fineparts, world->z_fineparts))
    sv = find_subvolume(world, (struct vector3 *) pos, sv);

  /* Find waypoint, compute trajectory from waypoint */
  wp = &(world->waypoints[sv - world->subvol]);
  origin = &(wp->loc);
  delta.x = pos->x - origin->x;
  delta.y = pos->y - origin->y;
  delta.z = pos->z - origin->z;

  for (wl=sv->wall_head ; wl!=NULL ; wl=wl->next)
  {
    struct vector3 hit_pos;
    double hit_time;
    int hit_check = collide_wall(origin, &delta, wl->this_wall, &hit_time, 
        &hit_pos, 0, world->rng, world->notify, &(world->ray_polygon_tests));

    if (hit_check!=COLLIDE_MISS)
    {
      world->ray_polygon_colls++;

      if ((hit_time>-EPS_C && hit_time<EPS_C) || (hit_time>1.0-EPS_C && hit_time<1.0+EPS_C))
      {
        bad_location = 1;
        break;
      }

      for (cur_region = wl->this_wall->counting_regions;
           cur_region != NULL;
           cur_region = cur_region->next)
      {
        struct region_list *crossed_region = (struct region_list*) CHECKED_MEM_GET(sv->local_storage->regl, "region list");
        crossed_region->reg = cur_region->reg;

        if (hit_check==COLLIDE_FRONT)
        {
          crossed_region->next = extra_in;
          extra_in = crossed_region;
        }
        else if (hit_check==COLLIDE_BACK)
        {
          crossed_region->next = extra_out;
          extra_out = crossed_region;
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
    return 0;
  }

  for (cur_region = extra_in;
       cur_region != NULL;
       cur_region = cur_region->next)
  {
    struct region_list *out_region = NULL;
    if (cur_region->reg == NULL) continue;
    for (out_region = extra_out;
         out_region != NULL;
         out_region = out_region->next)
    {
      if (out_region->reg==NULL) continue;
      if (cur_region->reg == out_region->reg)
      {
        cur_region->reg = NULL;
        out_region->reg = NULL;
        break;
      }
    }
  }

  result = eval_rel_region_3d(expression, wp, extra_in, extra_out);

  if (extra_in!=NULL) mem_put_list(sv->local_storage->regl,extra_in);
  if (extra_out!=NULL) mem_put_list(sv->local_storage->regl,extra_out);
  return result;
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
static int 
release_inside_regions(struct volume *world, struct release_site_obj *rso, 
    struct volume_molecule *m, int n)
{
  struct volume_molecule *new_m;
  struct release_region_data *rrd;
  struct subvolume *sv = NULL;
  double num_to_release;

  rrd = rso->region_data;
  new_m = NULL;
  m->previous_wall = NULL;
  m->index = -1;

  if (rso->release_number_method==CCNNUM)
  {
    double vol = (rrd->urb.x-rrd->llf.x)*(rrd->urb.y-rrd->llf.y)*(rrd->urb.z-rrd->llf.z);
    num_to_release = (N_AV*1e-15*rso->concentration*vol*world->length_unit*world->length_unit*world->length_unit) + 0.5;
    if (num_to_release > INT_MAX)
      mcell_error("Release site \"%s\" tries to release more than INT_MAX (2147483647) molecules.", rso->name);
    n = (int)(num_to_release);
  }

  if (n<0) 
    return vacuum_inside_regions(world, rso,m,n);
  if(world->notify->release_events == NOTIFY_FULL)
  {
     if (n > 0)
       mcell_log_raw("Releasing %d molecules %s ...", n, m->properties->sym->name);
  }

  long long skipped_placements = 0;
  int can_place = 1;
  int nfailures = 0;
  while (n>0)
  {
    m->pos.x = rrd->llf.x + (rrd->urb.x-rrd->llf.x)*rng_dbl(world->rng);
    m->pos.y = rrd->llf.y + (rrd->urb.y-rrd->llf.y)*rng_dbl(world->rng);
    m->pos.z = rrd->llf.z + (rrd->urb.z-rrd->llf.z)*rng_dbl(world->rng);

    if (!is_point_inside_region(world, &m->pos, rrd->expression, NULL))
    {
      if (rso->release_number_method==CCNNUM) n--;
      continue;
    }

    can_place = 1;
    if (m->properties->flags & IS_COMPLEX)
    {
      int subunit_idx;
      struct complex_species *cspec = (struct complex_species *) m->properties;
      sv = find_subvolume(world, &m->pos, NULL);
      for (subunit_idx = 0; subunit_idx < cspec->num_subunits; ++ subunit_idx)
      {
        struct vector3 subunit_pos;
        subunit_pos.x = m->pos.x + cspec->rel_locations[ subunit_idx ].x;
        subunit_pos.y = m->pos.y + cspec->rel_locations[ subunit_idx ].y;
        subunit_pos.z = m->pos.z + cspec->rel_locations[ subunit_idx ].z;
        if (!is_point_inside_region(world, &subunit_pos, rrd->expression, sv))
        {
          can_place = 0;
          break;
        }
      }
    }

    if (! can_place)
    {
      if (++ nfailures >= world->complex_placement_attempts)
      {
        -- n;
        if (++ skipped_placements >= world->notify->complex_placement_failure_threshold)
        {
          switch (world->notify->complex_placement_failure)
          {
            case WARN_COPE:
              break;

            case WARN_WARN:
              if (world->notify->release_events == NOTIFY_FULL)
                mcell_log_raw("\n");
              mcell_warn("Failed to place volume macromolecule '%s' in region %d times in a row.\n"
                         "         Leaving %lld molecules unplaced.",
                         m->properties->sym->name,
                         nfailures,
                         n + skipped_placements);
              break;

            case WARN_ERROR:
              if (world->notify->release_events == NOTIFY_FULL)
                mcell_log_raw("\n");
              mcell_error("Failed to place volume macromolecule '%s' in region %d times in a row.",
                          m->properties->sym->name,
                          nfailures);
              return 1;

            default: UNHANDLED_CASE(world->notify->complex_placement_failure);
          }
          break;
        }
        nfailures = 0;
      }
      continue;
    }

    /* Actually place the molecule */
    nfailures = 0;
    m->subvol = sv;
    if (m->properties->flags & IS_COMPLEX)
      new_m = macro_insert_molecule_volume(world, m, new_m);
    else
      new_m = insert_volume_molecule(world, m, new_m);
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
int 
release_molecules(struct volume *world, struct release_event_queue *req)
{
  struct release_site_obj *rso;
  struct release_pattern *rpat;
  struct volume_molecule m;
  struct abstract_molecule *ap;
  struct grid_molecule *gp;
  struct volume_molecule *guess;
  int i,i_failed,number;
  double num_to_release;
  short orient;
  struct vector3 *diam_xyz;
  struct vector3 pos;
  double diam,vol;
  double k;
  struct release_single_molecule *rsm;
  double location[1][4];

  if (req == NULL) return 0;
  rso = req->release_site;
  rpat = rso->pattern;

  memset(&m, 0, sizeof(struct volume_molecule));
  ap = (struct abstract_molecule*)(&m);

  /* Set up canonical molecule to be released */
  /* If we have a list, assume a 3D molecule and fix later */
  if (rso->mol_list!=NULL || (rso->mol_type->flags & NOT_FREE)==0)
  {
    m.flags = TYPE_3D | IN_VOLUME;
  }
  else
  {
    m.flags = TYPE_GRID | IN_SURFACE;
  }
  m.flags |= IN_SCHEDULE | ACT_NEWBIE;

  if (req->train_counter == 0)
  {
    req->train_counter++;
  }

  guess = NULL;

  /* Skip events that happened in the past (delay<0 or after checkpoint) */
  if (req->event_time < world->it_time && rso->release_prob!=MAGIC_PATTERN_PROBABILITY)
  {
    do
    {
      /* Schedule next release event and leave the function.
         This part of the code is relevant to checkpointing. */
      if (rso->release_prob < 1.0)
      {
        if (rso->release_prob == 0) return 0;
        req->event_time += rpat->release_interval;
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
      if (schedule_add(world->releaser,req))
        mcell_allocfailed("Failed to add release request to scheduler.");
    }
    return 0;
  }

  /* check whether the release will happen */
  if (rso->release_prob < 1.0)
  {
     k  = rng_dbl(world->rng);
     if (rso->release_prob < k)
     {
        /* make sure we will try the release pattern again in the future */
        req->event_time += rpat->release_interval;

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
            if (schedule_add(world->releaser,req))
               mcell_allocfailed("Failed to add release request to scheduler.");
        }
        return 0;

     }

  }

  /* Set molecule characteristics. */
  m.t = req->event_time;
  m.properties = rso->mol_type;
  m.t2 = 0.0;
  m.birthday = m.t;
  m.cmplx = NULL;

  if (rso->mol_list==NULL)  /* All molecules are the same, so we can set flags */
  {
    if (trigger_unimolecular(world->reaction_hash, world->rx_hashsize,
          rso->mol_type->hashval , ap) != NULL || (rso->mol_type->flags&CAN_GRIDWALL)!=0) ap->flags |= ACT_REACT;
    if (rso->mol_type->space_step > 0.0) ap->flags |= ACT_DIFFUSE;
  }

  switch (rso->release_number_method)
  {
    case CONSTNUM:
      num_to_release = rso->release_number;
      if (num_to_release > INT_MAX)
        mcell_error("Release site \"%s\" tries to release more than INT_MAX (2147483647) molecules.", rso->name);
      number = (int)(num_to_release);
      break;

    case GAUSSNUM:
      if (rso->standard_deviation > 0)
      {
        num_to_release = (rng_gauss(world->rng)*rso->standard_deviation + rso->release_number);
        if (num_to_release > INT_MAX)
          mcell_error("Release site \"%s\" tries to release more than INT_MAX (2147483647) molecules.", rso->name);
        number = (int)(num_to_release);
      }
      else
      {
        rso->release_number_method = CONSTNUM;
        num_to_release = rso->release_number;
        if (num_to_release > INT_MAX)
          mcell_error("Release site \"%s\" tries to release more than INT_MAX (2147483647) molecules.", rso->name);
        number = (int)(num_to_release);
      }
      break;

    case VOLNUM:
      diam = rso->mean_diameter;
      if (rso->standard_deviation > 0)
      {
        diam += rng_gauss(world->rng)*rso->standard_deviation;
      }
      vol = (MY_PI/6.0) * diam*diam*diam;
      num_to_release = N_AV * 1e-15 * rso->concentration * vol + 0.5;
      if (num_to_release > INT_MAX)
        mcell_error("Release site \"%s\" tries to release more than INT_MAX (2147483647) molecules.", rso->name);
      number = (int)(num_to_release);
      break;

    case CCNNUM:
    case DENSITYNUM:
      if (rso->diameter==NULL) number = 0;
      else
      {
        switch (rso->release_shape)
        {
          case SHAPE_SPHERICAL:
          case SHAPE_ELLIPTIC:
            vol = (1.0/6.0)*MY_PI*rso->diameter->x*rso->diameter->y*rso->diameter->z;
            break;
          case SHAPE_RECTANGULAR:
          case SHAPE_CUBIC:
            vol = rso->diameter->x*rso->diameter->y*rso->diameter->z;
            break;

          case SHAPE_SPHERICAL_SHELL:
            mcell_error("Release site \"%s\" tries to release a concentration on a spherical shell.", rso->name);
            vol = 0;
            break;

          default:
            mcell_internal_error("Release by concentration on invalid release site shape (%d) for release site \"%s\".", rso->release_shape, rso->name);
            break;
        }
        num_to_release = N_AV * 1e-15 * rso->concentration * vol * world->length_unit*world->length_unit*world->length_unit + 0.5;
        if (num_to_release > INT_MAX)
          mcell_error("Release site \"%s\" tries to release more than INT_MAX (2147483647) molecules.", rso->name);
         number = (int)(num_to_release);
      }
      break;

    default:
      mcell_internal_error("Release site \"%s\" has invalid release number method (%d).", rso->name, rso->release_number_method);
      number = 0;
      break;
  }

  if (rso->release_shape == SHAPE_REGION)
  {
    u_int pop_before = ap->properties->population;
    if (ap->flags & TYPE_3D)
    {
      if (release_inside_regions(world, rso,
            (struct volume_molecule*)ap,number))
        return 1;

      if (world->notify->release_events==NOTIFY_FULL)
      {
        if (number >= 0)
        {
          mcell_log("  Released %d %s from \"%s\" at iteration %lld.",
            ap->properties->population-pop_before, rso->mol_type->sym->name, rso->name, world->it_time);
        }
        else
        {
          mcell_log("  Removed %d %s from \"%s\" at iteration %lld.",
            pop_before-ap->properties->population, rso->mol_type->sym->name, rso->name, world->it_time);
        }
      }
    }
    else
    {
      i = release_onto_regions(world, rso,(struct grid_molecule*)ap,number);
      if (i) return 1;

      if (world->notify->release_events==NOTIFY_FULL)
      {
        if (number >= 0)
        {
          mcell_log("  Released %d %s from \"%s\" at iteration %lld.",
            ap->properties->population-pop_before, rso->mol_type->sym->name, rso->name, world->it_time);
        }
        else
        {
          mcell_log("  Removed %d %s from \"%s\" at iteration %lld.",
            pop_before-ap->properties->population, rso->mol_type->sym->name, rso->name, world->it_time);
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
      i = 0; /* serves as counter for released molecules */
      i_failed = 0; /* serves as counted for the failed to release molecules */
      for (; rsm!=NULL ; rsm=rsm->next)
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
          if ((rsm->mol_type->flags & IS_COMPLEX))
          {
            guess = macro_insert_molecule_volume(world, &m, guess);
            i++;
          }
          else
          {
            m.properties = rsm->mol_type;

            /* Have to set flags, since insert_volume_molecule doesn't */
            if (trigger_unimolecular(world->reaction_hash, world->rx_hashsize,
                  ap->properties->hashval , ap) != NULL ||
                (ap->properties->flags&CAN_GRIDWALL)!=0)
            {
              ap->flags |= ACT_REACT;
            }
            if (m.properties->space_step > 0.0) ap->flags |= ACT_DIFFUSE;
            guess = insert_volume_molecule(world, &m, guess);
            i++;
          }
          if (guess==NULL) return 1;
        }
        else
        {
          if (diam_xyz==NULL) diam=0.0;
          else diam=diam_xyz->x;

          if (rsm->orient>0) orient=1;
          else if (rsm->orient<0) orient=-1;
          else 
          {
             orient = (rng_uint(world->rng)&1)?1:-1;
          }

          /* Don't have to set flags, insert_grid_molecule takes care of it */
          if ((rsm->mol_type->flags & IS_COMPLEX))
          {
            /* XXX: Retry? */
            gp = macro_insert_molecule_grid(world, rsm->mol_type, &m.pos, 
                orient, diam, req->event_time);
          }
          else
          {
            gp = insert_grid_molecule(world, rsm->mol_type, &m.pos, orient, 
                diam, req->event_time, NULL);
          }
          if (gp==NULL)
          {
            mcell_warn("Molecule release is unable to find surface upon which to place molecule %s.\n"
                       "  This could be caused by too small of a SITE_DIAMETER on the release site '%s'.",
                       rsm->mol_type->sym->name,
                       rso->name);
            i_failed++;
          }
          else
            i++;
        }
      }
      if (world->notify->release_events==NOTIFY_FULL)
      {
          mcell_log("Releasing %d molecules from list \"%s\" at iteration %lld.", i, rso->name, world->it_time);
      }
      if (i_failed > 0)
          mcell_warn("Failed to release %d molecules from list \"%s\" at iteration %lld.", i_failed, rso->name, world->it_time);
    }
    else if (diam_xyz != NULL)
    {

      if (world->notify->release_events == NOTIFY_FULL)
      {
        if (number > 0)
          mcell_log_raw("Releasing %d molecules %s ...", number, rso->mol_type->sym->name);
      }
      const int is_spheroidal =
            (rso->release_shape == SHAPE_SPHERICAL ||
             rso->release_shape == SHAPE_ELLIPTIC  ||
             rso->release_shape == SHAPE_SPHERICAL_SHELL);
      for (i=0;i<number;i++)
      {
        do /* Pick values in unit square, toss if not in unit circle */
        {
          pos.x = (rng_dbl(world->rng)-0.5);
          pos.y = (rng_dbl(world->rng)-0.5);
          pos.z = (rng_dbl(world->rng)-0.5);
        } while (is_spheroidal
                  && pos.x*pos.x + pos.y*pos.y + pos.z*pos.z >= 0.25);

        if (rso->release_shape == SHAPE_SPHERICAL_SHELL)
        {
          double r;
          r = sqrt(pos.x*pos.x + pos.y*pos.y + pos.z*pos.z)*2.0;
          if (r==0.0)
          { pos.x = 0.0; pos.y = 0.0; pos.z = 0.5; }
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
        if ((m.properties->flags & IS_COMPLEX))
          guess = macro_insert_molecule_volume(world, &m, guess);
        else
          guess = insert_volume_molecule(world, &m,guess);  /* Insert copy of m into world */
        if (guess == NULL) return 1;
      }
      if (world->notify->release_events==NOTIFY_FULL)
      {
        mcell_log("  Released %d %s from \"%s\" at iteration %lld.", number,rso->mol_type->sym->name, rso->name, world->it_time);
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

      if (world->notify->release_events == NOTIFY_FULL)
      {
        if (number > 0)
          mcell_log_raw("Releasing %d molecules %s ...", number, rso->mol_type->sym->name);
      }

      for (i=0;i<number;i++)
      {
        if ((rso->mol_type->flags & IS_COMPLEX))
          guess = macro_insert_molecule_volume(world, &m, guess);
        else
          guess = insert_volume_molecule(world, &m, guess);
        if (guess == NULL) return 1;
      }
      if (world->notify->release_events==NOTIFY_FULL)
      {
        mcell_log("  Released %d %s from \"%s\" at iteration %lld.", number,rso->mol_type->sym->name, rso->name, world->it_time);
      }
    }
  }


  /* Schedule next release event. */
  if (rso->release_prob==MAGIC_PATTERN_PROBABILITY) return 0;  /* Triggered by reaction, don't schedule */
  req->event_time += rpat->release_interval;

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
    if (schedule_add(world->releaser,req))
      mcell_allocfailed("Failed to add release request to scheduler.");
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
static void 
find_exponential_params(double c, double C, double d, double N, double *A,
    double *B, double *k)
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
  *A = d / (exp(*k) - 1.0);
  *B = c - *A;
}



/*************************************************************************
 check_partitions_against_interaction_diameter:
  In: nothing.  Uses struct volume *world, assumes partitions are set.
  Out: 0 on success, 1 on error
*************************************************************************/
static int 
check_partitions_against_interaction_diameter(struct volume *world)
{
  int i;

  if (world->x_partitions!=NULL)
  {
    for (i=1;i<world->nx_parts;i++)
    {
      if (world->x_partitions[i] - world->x_partitions[i-1] < 2*world->rx_radius_3d)
      {
        mcell_error("X partitions closer than interaction diameter\n"
                    "  X partition #%d at %g\n"
                    "  X partition #%d at %g\n"
                    "  Interaction diameter %g",
                    i,   world->length_unit*world->x_partitions[i-1],
                    i+1, world->length_unit*world->x_partitions[i],
                    2*world->length_unit*world->rx_radius_3d);
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
        mcell_error("Y partitions closer than interaction diameter\n"
                    "  Y partition #%d at %g\n"
                    "  Y partition #%d at %g\n"
                    "  Interaction diameter %g",
                    i,   world->length_unit*world->y_partitions[i-1],
                    i+1, world->length_unit*world->y_partitions[i],
                    2*world->length_unit*world->rx_radius_3d);
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
        mcell_error("Z partitions closer than interaction diameter\n"
                    "  Z partition #%d at %g\n"
                    "  Z partition #%d at %g\n"
                    "  Interaction diameter %g\n",
                    i,   world->length_unit*world->z_partitions[i-1],
                    i+1, world->length_unit*world->z_partitions[i],
                    2*world->length_unit*world->rx_radius_3d);
        return 1;
      }
    }
  }
  return 0;
}



/*************************************************************************
set_partitions:
  In: nothing.  Uses struct volume *world, assumes bounding box is set.
  Out: 0 on success, 1 on error; coarse and fine partitions are set.
*************************************************************************/
int 
set_partitions(struct volume *world)
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
  smallest_spacing = 0.1 * world->r_length_unit;  /* 100nm */

  if (2*world->rx_radius_3d > smallest_spacing) smallest_spacing=2*world->rx_radius_3d;

  /* We have 2^15 possible fine partitions; we'll use 24k of them */
  if (world->n_fineparts != 4096 + 16384 + 4096)
  {

    world->n_fineparts = 4096 + 16384 + 4096;
    world->x_fineparts = CHECKED_MALLOC_ARRAY(double, world->n_fineparts, "x fine partitions");
    world->y_fineparts = CHECKED_MALLOC_ARRAY(double, world->n_fineparts, "y fine partitions");
    world->z_fineparts = CHECKED_MALLOC_ARRAY(double, world->n_fineparts, "z fine partitions");
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
    if (world->notify->progress_report != NOTIFY_NONE)
      mcell_log_raw("Rescaling: was %.3f to %.3f, now ",f_min*world->length_unit,f_max*world->length_unit);
    f = smallest_spacing - (f_max-f_min);
    f_max += 0.5*f;
    f_min -= 0.5*f;
    if (world->notify->progress_report != NOTIFY_NONE)
      mcell_log_raw("%.3f to %.3f\n",f_min*world->length_unit,f_max*world->length_unit);
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
    if (world->notify->progress_report != NOTIFY_NONE)
      mcell_log_raw("Rescaling: was %.3f to %.3f, now ",f_min*world->length_unit,f_max*world->length_unit);
    f = smallest_spacing - (f_max-f_min);
    f_max += 0.5*f;
    f_min -= 0.5*f;
    if (world->notify->progress_report != NOTIFY_NONE)
      mcell_log_raw("%.3f to %.3f\n",f_min*world->length_unit,f_max*world->length_unit);
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
    if (world->notify->progress_report != NOTIFY_NONE)
      mcell_log_raw("Rescaling: was %.3f to %.3f, now ",f_min*world->length_unit,f_max*world->length_unit);
    f = smallest_spacing - (f_max-f_min);
    f_max += 0.5*f;
    f_min -= 0.5*f;
    if (world->notify->progress_report != NOTIFY_NONE)
      mcell_log_raw("%.3f to %.3f\n",f_min*world->length_unit,f_max*world->length_unit);
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
  if (check_partitions_against_interaction_diameter(world))
    return 1;


  /* Use automatic partitioning only when there are no user-specified partitions */
  if (world->x_partitions!=NULL || world->y_partitions!=NULL || world->z_partitions!=NULL)
  {
    if (world->x_partitions == NULL)
      mcell_error("Some axes are partitioned, but the X-axis is not.");
    if (world->y_partitions == NULL)
      mcell_error("Some axes are partitioned, but the Y-axis is not.");
    if (world->z_partitions == NULL)
      mcell_error("Some axes are partitioned, but the Z-axis is not.");
  }

  if (world->x_partitions == NULL &&
      world->y_partitions == NULL &&
      world->z_partitions == NULL)
  {
     /* perform automatic partitioning */

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
    world->x_partitions = CHECKED_MALLOC_ARRAY(double, world->nx_parts, "x partitions");
    world->y_partitions = CHECKED_MALLOC_ARRAY(double, world->ny_parts, "y partitions");
    world->z_partitions = CHECKED_MALLOC_ARRAY(double, world->nz_parts, "z partitions");

    /* Calculate aspect ratios so that subvolumes are approximately cubic */
    x_aspect = (part_max.x - part_min.x) / f_max;
    y_aspect = (part_max.y - part_min.y) / f_max;
    z_aspect = (part_max.z - part_min.z) / f_max;

    x_in = floor((world->nx_parts - 2) * x_aspect + 0.5);
    y_in = floor((world->ny_parts - 2) * y_aspect + 0.5);
    z_in = floor((world->nz_parts - 2) * z_aspect + 0.5);
    if (x_in < 2) x_in = 2;
    if (y_in < 2) y_in = 2;
    if (z_in < 2) z_in = 2;

    /* If we've violated our 2*reaction radius criterion, fix it */
    smallest_spacing = 2*world->rx_radius_3d;
    if ((part_max.x-part_min.x)/(x_in-1) < smallest_spacing)
    {
      x_in = 1 + floor((part_max.x-part_min.x)/smallest_spacing);
    }
    if ((part_max.y-part_min.y)/(y_in-1) < smallest_spacing)
    {
      y_in = 1 + floor((part_max.y-part_min.y)/smallest_spacing);
    }
    if ((part_max.z-part_min.z)/(z_in-1) < smallest_spacing)
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
        dbl_array = CHECKED_MALLOC_ARRAY(double, (world->nx_parts+1), "x partitions (expanded in -X dir)");
        dbl_array[0] = world->x_partitions[0];
        dbl_array[1] = world->bb_llf.x - dfx;
        memcpy(&(dbl_array[2]),&(world->x_partitions[1]),sizeof(double)*(world->nx_parts-1));
        free(world->x_partitions);
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
        dbl_array = CHECKED_MALLOC_ARRAY(double, (world->nx_parts+1), "x partitions (expanded in +X dir)");
        dbl_array[world->nx_parts] = world->x_partitions[world->nx_parts-1];
        dbl_array[world->nx_parts-1] = world->bb_urb.x + dfx;
        memcpy(dbl_array,world->x_partitions,sizeof(double)*(world->nx_parts-1));
        free(world->x_partitions);
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
        dbl_array = CHECKED_MALLOC_ARRAY(double, (world->ny_parts+1), " y partitions (expanded in -Y dir)");
        dbl_array[0] = world->y_partitions[0];
        dbl_array[1] = world->bb_llf.y - dfy;
        memcpy(&(dbl_array[2]),&(world->y_partitions[1]),sizeof(double)*(world->ny_parts-1));
        free(world->y_partitions);
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
        dbl_array = CHECKED_MALLOC_ARRAY(double, (world->ny_parts+1), "y partitions (expanded in +Y dir)");
        dbl_array[world->ny_parts] = world->y_partitions[world->ny_parts-1];
        dbl_array[world->ny_parts-1] = world->bb_urb.y + dfy;
        memcpy(dbl_array,world->y_partitions,sizeof(double)*(world->ny_parts-1));
        free(world->y_partitions);
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
        dbl_array = CHECKED_MALLOC_ARRAY(double, (world->nz_parts+1), "z partitions (expanded in -Z dir)");
        dbl_array[0] = world->z_partitions[0];
        dbl_array[1] = world->bb_llf.z - dfz;
        memcpy(&(dbl_array[2]),&(world->z_partitions[1]),sizeof(double)*(world->nz_parts-1));
        free(world->z_partitions);
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
        dbl_array = CHECKED_MALLOC_ARRAY(double, (world->nz_parts+1), "z partitions (expanded in +Z dir)");
        dbl_array[world->nz_parts] = world->z_partitions[world->nz_parts-1];
        dbl_array[world->nz_parts-1] = world->bb_urb.z + dfz;
        memcpy(dbl_array,world->z_partitions,sizeof(double)*(world->nz_parts-1));
        free(world->z_partitions);
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
    mcell_log_raw("X partitions: ");
    mcell_log_raw("-inf ");
    for (i=1;i<world->nx_parts - 1;i++) mcell_log_raw("%.5f ",world->length_unit * world->x_partitions[i]);
    mcell_log_raw("inf\n");
    mcell_log_raw("Y partitions: ");
    mcell_log_raw("-inf ");
    for (i=1;i<world->ny_parts - 1;i++) mcell_log_raw("%.5f ",world->length_unit * world->y_partitions[i]);
    mcell_log_raw("inf\n");
    mcell_log_raw("Z partitions: ");
    mcell_log_raw("-inf ");
    for (i=1;i<world->nz_parts - 1;i++) mcell_log_raw("%.5f ",world->length_unit * world->z_partitions[i]);
    mcell_log_raw("inf\n");
  }

  return 0;
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
void 
path_bounding_box(struct vector3 *loc, struct vector3 * displacement,
    struct vector3 *llf, struct vector3 *urb, double rx_radius_3d)
{
   struct vector3 final;  /* final position of the molecule after random walk */
   double R;     /* molecule interaction radius */


   R = rx_radius_3d;
   vect_sum(loc, displacement, &final);

   llf->x = urb->x = loc->x;
   llf->y = urb->y = loc->y;
   llf->z = urb->z = loc->z;

   if (final.x < llf->x)
   {
         llf->x = final.x;
   }
   if (final.x > urb->x)
   {
       urb->x = final.x;
   }
   if (final.y < llf->y)
   {
         llf->y = final.y;
   }
   if (final.y > urb->y)
   {
       urb->y = final.y;
   }
   if (final.z < llf->z)
   {
         llf->z = final.z;
   }
   if (final.z > urb->z)
   {
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



/***************************************************************************
 This function puts volume molecule
 in the random positions in the world.
 It is a research function that should not be called
 during regular MCell3 execution.
 In: volume_molecule
     vector pointing to the low-left-corner of the world bounding box
     sizes of the world bounding box in 3 dimensions
 Out: if molecule moves out of the subvolume a new copy of that molecule is
      created and rescheduled, otherwise the existing molecule gets random
      position in the original subvolume it belonged to.
***************************************************************************/
void 
randomize_vol_mol_position(struct volume *world, struct volume_molecule *mp, 
    struct vector3 *low_end, double size_x, double size_y, double size_z)
{
   double num; /* random number */
   struct subvolume *new_sv, *old_sv;
   struct vector3 loc;
   struct volume_molecule *new_mp;

    /* find future molecule position */
    num = rng_dbl(world->rng);
    loc.x = low_end->x + num*size_x;
    num = rng_dbl(world->rng);
    loc.y = low_end->y + num*size_y;
    num = rng_dbl(world->rng);
    loc.z = low_end->z + num*size_z;
    /* find old subvolume */
    old_sv = find_subvolume(world, &(mp->pos), NULL);


    /* now remove molecule from old subvolume
    and place it into the new location into new one */
    mp->pos.x = loc.x;
    mp->pos.y = loc.y;
    mp->pos.z = loc.z;
    if(!inside_subvolume(&(mp->pos), old_sv, world->x_fineparts,
          world->y_fineparts, world->z_fineparts))
    {
       /* find new subvolume after reshuffling */
       new_sv = find_subvolume(world, &loc, NULL);
       new_mp = migrate_volume_molecule(mp, new_sv);
       if (schedule_add(new_sv->local_storage->timer, (struct abstract_molecule *)new_mp))
         mcell_allocfailed("Failed to add volume molecule to scheduler.");
    }
}


/***************************************************************************
 collect_molecule:
    Perform garbage collection on a discarded molecule.  If the molecule is no
    longer in any lists, it will be freed.

 In: m: the molecule
 Out: Nothing.  Molecule is unlinked from its list in the subvolume, and
      possibly returned to its birthplace.
***************************************************************************/
void 
collect_molecule(struct volume_molecule *m)
{
  /* Unlink from the previous item */
  if (m->prev_v != NULL)
  {
#ifdef DEBUG_LIST_CHECKS
    if (*m->prev_v != m)
    {
      mcell_error_nodie("Stale previous pointer!  ACK!  THRBBPPPPT!");
    }
#endif
  *(m->prev_v) = m->next_v;
  }

  /* Unlink from the following item */
  if (m->next_v != NULL)
  {
#ifdef DEBUG_LIST_CHECKS
    if (m->next_v->prev_v != &m->next_v)
    {
      mcell_error_nodie("Stale next pointer!  ACK!  THRBBPPPPT!");
    }
#endif
    m->next_v->prev_v = m->prev_v;
  }

  /* Clear our next/prev pointers */
  m->prev_v = NULL;
  m->next_v = NULL;

  /* Dispose of the molecule */
  m->properties = NULL;
  m->flags &= ~IN_VOLUME;
  if ((m->flags & IN_MASK) == 0)
    mem_put(m->birthplace, m);
}



/***************************************************************************
 ht_add_molecule_to_list:
    Add a molecule to the appropriate molecule list in a subvolume's pointer
    hash.  It is assumed that the molecule's subvolume pointer is valid and
    points to the right subvolume.

    If the molecule takes part in any mol-mol interactions (including
    trimolecular reactions involving two or more volume molecules), it is added
    to a list containing only molecules of the same species.  If it does NOT
    take part in any such interactions, it is added to a single molecule list
    which keeps track of all molecules which do not interact with other volume
    molecules.

 In: h: the pointer hash to which to add the molecule
     m: the molecule
 Out: Nothing.  Molecule is added to the subvolume's molecule lists.
***************************************************************************/
void 
ht_add_molecule_to_list(struct pointer_hash *h, struct volume_molecule *m)
{
  struct per_species_list *list = NULL;

  /* If the molecule does not interact with other volume molecules... */
  if (! (m->properties->flags & (CAN_MOLMOL|CAN_MOLMOLMOL|CAN_MOLMOLGRID)))
  {
    /* If we have a list for these molecules, it's always at the head of the
     * species lists, and the list always has the species set to NULL. */
    if (m->subvol->species_head != NULL  &&
        m->subvol->species_head->properties == NULL)
    {
      list = m->subvol->species_head;
    }
    else
    {
      list = (struct per_species_list *) CHECKED_MEM_GET(m->subvol->local_storage->pslv, "per-species molecule list");
      list->properties = NULL;
      list->head = NULL;
      list->next = m->subvol->species_head;
      m->subvol->species_head = list;
    }
  }

  /* The molecule DOES interact with other volume molecules */
  else
  {
    /* See if we have a list */
    list = (struct per_species_list *) pointer_hash_lookup(h, m->properties, m->properties->hashval);

    /* If not, create one and add it in */
    if (list == NULL)
    {
      list = (struct per_species_list *) CHECKED_MEM_GET(m->subvol->local_storage->pslv, "per-species molecule list");
      list->properties = m->properties;
      list->head = NULL;
      if (pointer_hash_add(h, m->properties, m->properties->hashval, list))
        mcell_allocfailed("Failed to add species to subvolume species table.");

      /* If the first per-species list is for non-volume-interacting molecules,
       * add our new list after that. */
      if (m->subvol->species_head != NULL  &&
          m->subvol->species_head->properties == NULL)
      {
        list->next = m->subvol->species_head->next;
        m->subvol->species_head->next = list;
      }

      /* Otherwise, add it to the beginning of the list */
      else
      {
        list->next = m->subvol->species_head;
        m->subvol->species_head = list;
      }
    }
  }

  /* Link the molecule into the list */
  m->next_v = list->head;
  if (list->head)
    list->head->prev_v = &m->next_v;
  m->prev_v = &list->head;
  list->head = m;
}


/***************************************************************************
 ht_remove:
    Remove a species list from a pointer hash.  This is a fairly simple wrapper
    around pointer_hash_remove, which simply looks up the key (the species),
    and avoids trying to remove the non-interacting molecules list from the
    pointer hash, since it should never be added in the first place.

    Note that the per-species list must still be valid at this point -- it must
    not have been freed, nor its 'properties' pointer nilled.

 In: h: the pointer hash to which to add the molecule
     psl: the species list to remove
 Out: Nothing.  Molecule is added to the subvolume's molecule lists.
***************************************************************************/
void 
ht_remove(struct pointer_hash *h, struct per_species_list *psl)
{
  struct species *s = psl->properties;
  if (s == NULL)
    return;

  (void) pointer_hash_remove(h, s, s->hashval);
}
