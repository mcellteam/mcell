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

#include <assert.h>
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "diffuse.h"
#include "rng.h"
#include "util.h"
#include "logging.h"
#include "mcell_structs.h"
#include "count_util.h"
#include "grid_util.h"
#include "vol_util.h"
#include "wall_util.h"
#include "react.h"
#include "react_output.h"

/**********************************************************************
ray_trace_trimol:
  In: molecule that is moving
      linked list of potential collisions with molecules (we could react)
      subvolume that we start in
      displacement vector from current to new location
      wall we have reflected off of and should not hit again
      start time of the  molecule random walk (local
         to the molecule timestep)
  Out: collision list of walls and molecules we intersected along our ray
       (current subvolume only), plus the subvolume wall.  Will always
       return at least the subvolume wall--NULL indicates an out of
       memory eriror.
  Note: This is a version of the "ray_trace()" function adapted for
        the case when moving molecule can engage in trimolecular collisions

**********************************************************************/
struct sp_collision *ray_trace_trimol(struct volume *world,
                                      struct volume_molecule *m,
                                      struct sp_collision *c,
                                      struct subvolume *sv, struct vector3 *v,
                                      struct wall *reflectee,
                                      double walk_start_time) {
  struct sp_collision *smash, *shead;
  struct abstract_molecule *a;
  struct wall_list *wlp;
  struct wall_list fake_wlp;
  double dx, dy, dz;
  /* time, in units of of the molecule's time step, at which molecule
     will cross the x,y,z partitions, respectively. */
  double tx, ty, tz;
  int i, j, k;

  world->ray_voxel_tests++;

  shead = NULL;
  smash = (struct sp_collision *)CHECKED_MEM_GET(sv->local_storage->sp_coll,
                                                 "collision structure");

  fake_wlp.next = sv->wall_head;

  for (wlp = sv->wall_head; wlp != NULL; wlp = wlp->next) {
    if (wlp->this_wall == reflectee)
      continue;

    i = collide_wall(&(m->pos), v, wlp->this_wall, &(smash->t), &(smash->loc),
                     1, world->rng, world->notify, &(world->ray_polygon_tests));
    if (i == COLLIDE_REDO) {
      if (shead != NULL)
        mem_put_list(sv->local_storage->sp_coll, shead);
      shead = NULL;
      wlp = &fake_wlp;
      continue;
    } else if (i != COLLIDE_MISS) {
      world->ray_polygon_colls++;

      smash->what = COLLIDE_WALL + i;
      smash->moving = m->properties;
      smash->target = (void *)wlp->this_wall;
      smash->t_start = walk_start_time;
      smash->pos_start.x = m->pos.x;
      smash->pos_start.y = m->pos.y;
      smash->pos_start.z = m->pos.z;
      smash->sv_start = sv;

      smash->disp.x = v->x;
      smash->disp.y = v->y;
      smash->disp.z = v->z;

      smash->next = shead;
      shead = smash;
      smash = (struct sp_collision *)CHECKED_MEM_GET(sv->local_storage->sp_coll,
                                                     "collision structure");
    }
  }

  dx = dy = dz = 0.0;
  i = -10;
  if (v->x < 0.0) {
    dx = world->x_fineparts[sv->llf.x] - m->pos.x;
    i = 0;
  } else if (v->x > 0.0) {
    dx = world->x_fineparts[sv->urb.x] - m->pos.x;
    i = 1;
  }

  j = -10;
  if (v->y < 0.0) {
    dy = world->y_fineparts[sv->llf.y] - m->pos.y;
    j = 0;
  } else if (v->y > 0.0) {
    dy = world->y_fineparts[sv->urb.y] - m->pos.y;
    j = 1;
  }

  k = -10;
  if (v->z < 0.0) {
    dz = world->z_fineparts[sv->llf.z] - m->pos.z;
    k = 0;
  } else if (v->z > 0.0) {
    dz = world->z_fineparts[sv->urb.z] - m->pos.z;
    k = 1;
  }

  if (i + j + k < 0) /* At least one vector is zero */
  {
    if (i + j + k < -15) /* Two or three vectors are zero */
    {
      if (i >= 0) /* X is the nonzero one */
      {
        smash->t = dx / v->x;
        smash->what = COLLIDE_SUBVOL + COLLIDE_SV_NX + i;
      } else if (j >= 0) /* Y is nonzero */
      {
        smash->t = dy / v->y;
        smash->what = COLLIDE_SUBVOL + COLLIDE_SV_NY + j;
      } else if (k >= 0) /* Z is nonzero */
      {
        smash->t = dz / v->z;
        smash->what = COLLIDE_SUBVOL + COLLIDE_SV_NZ + k;
      } else {
        smash->t = FOREVER;
        smash->what =
            COLLIDE_SUBVOL; /*Wrong, but we'll never hit it, so it's ok*/
      }
    } else /* One vector is zero; throw out other two */
    {
      if (i < 0) {
        ty = fabs(dy * v->z);
        tz = fabs(v->y * dz);
        if (ty < tz) {
          smash->t = dy / v->y;
          smash->what = COLLIDE_SUBVOL + COLLIDE_SV_NY + j;
        } else {
          smash->t = dz / v->z;
          smash->what = COLLIDE_SUBVOL + COLLIDE_SV_NZ + k;
        }
      } else if (j < 0) {
        tx = fabs(dx * v->z);
        tz = fabs(v->x * dz);
        if (tx < tz) {
          smash->t = dx / v->x;
          smash->what = COLLIDE_SUBVOL + COLLIDE_SV_NX + i;
        } else {
          smash->t = dz / v->z;
          smash->what = COLLIDE_SUBVOL + COLLIDE_SV_NZ + k;
        }
      } else /* k<0 */
      {
        tx = fabs(dx * v->y);
        ty = fabs(v->x * dy);
        if (tx < ty) {
          smash->t = dx / v->x;
          smash->what = COLLIDE_SUBVOL + COLLIDE_SV_NX + i;
        } else {
          smash->t = dy / v->y;
          smash->what = COLLIDE_SUBVOL + COLLIDE_SV_NY + j;
        }
      }
    }
  } else /* No vectors are zero--use alternate method */
  {
    tx = fabs(dx * v->y * v->z);
    ty = fabs(v->x * dy * v->z);
    tz = fabs(v->x * v->y * dz);

    if (tx < ty) {
      if (tx < tz) {
        smash->t = dx / v->x;
        smash->what = COLLIDE_SUBVOL + COLLIDE_SV_NX + i;
      } else {
        smash->t = dz / v->z;
        smash->what = COLLIDE_SUBVOL + COLLIDE_SV_NZ + k;
      }
    } else {
      if (ty < tz) {
        smash->t = dy / v->y;
        smash->what = COLLIDE_SUBVOL + COLLIDE_SV_NY + j;
      } else {
        smash->t = dz / v->z;
        smash->what = COLLIDE_SUBVOL + COLLIDE_SV_NZ + k;
      }
    }
  }

  smash->loc.x = m->pos.x + smash->t * v->x;
  smash->loc.y = m->pos.y + smash->t * v->y;
  smash->loc.z = m->pos.z + smash->t * v->z;

  smash->moving = m->properties;
  smash->target = (void *)sv;
  smash->t_start = walk_start_time;
  smash->pos_start.x = m->pos.x;
  smash->pos_start.y = m->pos.y;
  smash->pos_start.z = m->pos.z;
  smash->sv_start = sv;

  smash->disp.x = v->x;
  smash->disp.y = v->y;
  smash->disp.z = v->z;

  smash->next = shead;
  shead = smash;

  for (; c != NULL; c = c->next) {
    a = (struct abstract_molecule *)c->target;
    if (a->properties == NULL)
      continue;

    i = collide_mol(&(m->pos), v, a, &(c->t), &(c->loc), world->rx_radius_3d);
    if (i != COLLIDE_MISS) {
      smash = (struct sp_collision *)CHECKED_MEM_GET(sv->local_storage->sp_coll,
                                                     "collision structure");
      memcpy(smash, c, sizeof(struct sp_collision));

      smash->t_start = walk_start_time;
      smash->pos_start.x = m->pos.x;
      smash->pos_start.y = m->pos.y;
      smash->pos_start.z = m->pos.z;

      smash->disp.x = v->x;
      smash->disp.y = v->y;
      smash->disp.z = v->z;

      smash->next = shead;
      shead = smash;
    }
  }

  return shead;
}

/****************************************************************************
expand_collision_partner_list:
  In: molecule that is moving
      displacement to the new location
      subvolume that we start in
  Out: Returns linked list of molecules from neighbor subvolumes
       that are located within "interaction_radius" from the the subvolume
       border.
       The molecules are added only when the molecule displacement
       bounding box intersects with the subvolume bounding box.
  Note:  This is a version of the function "expand_collision_list()"
        adapted for the case when molecule can engage in trimolecular
        collisions.
****************************************************************************/
static struct sp_collision *expand_collision_partner_list(
    struct volume_molecule *m, struct vector3 *mv, struct subvolume *sv,
    double rx_radius_3d, double *x_fineparts, double *y_fineparts,
    double *z_fineparts, int nx_parts, int ny_parts, int nz_parts,
    int rx_hashsize, struct rxn **reaction_hash) {
  struct sp_collision *shead1 = NULL;
  /* lower left and upper_right corners of the molecule path
     bounding box expanded by R. */
  struct vector3 path_llf, path_urb;
  double R; /* molecule interaction radius */
  R = (rx_radius_3d);

  /* find the molecule path bounding box. */
  path_bounding_box(&m->pos, mv, &path_llf, &path_urb, rx_radius_3d);

  /* Decide which directions we need to go */
  int x_neg = 0, x_pos = 0, y_neg = 0, y_pos = 0, z_neg = 0, z_pos = 0;
  if (!(sv->world_edge & X_POS_BIT) && path_urb.x + R > x_fineparts[sv->urb.x])
    x_pos = 1;
  if (!(sv->world_edge & X_NEG_BIT) && path_llf.x - R < x_fineparts[sv->llf.x])
    x_neg = 1;
  if (!(sv->world_edge & Y_POS_BIT) && path_urb.y + R > y_fineparts[sv->urb.y])
    y_pos = 1;
  if (!(sv->world_edge & Y_NEG_BIT) && path_llf.y - R < y_fineparts[sv->llf.y])
    y_neg = 1;
  if (!(sv->world_edge & Z_POS_BIT) && path_urb.z + R > z_fineparts[sv->urb.z])
    z_pos = 1;
  if (!(sv->world_edge & Z_NEG_BIT) && path_llf.z - R < z_fineparts[sv->llf.z])
    z_neg = 1;

  /* go +X */
  if (x_pos) {
    struct subvolume *newsv_x = sv + (nz_parts - 1) * (ny_parts - 1);
    shead1 = expand_collision_partner_list_for_neighbor(
        sv, m, mv, newsv_x, &path_llf, &path_urb, shead1, R, 0.0, 0.0,
        x_fineparts, y_fineparts, z_fineparts, rx_hashsize, reaction_hash);

    /* go +X, +Y */
    if (y_pos) {
      struct subvolume *newsv_y = newsv_x + (nz_parts - 1);
      shead1 = expand_collision_partner_list_for_neighbor(
          sv, m, mv, newsv_y, &path_llf, &path_urb, shead1, R, R, 0.0,
          x_fineparts, y_fineparts, z_fineparts, rx_hashsize, reaction_hash);

      /* go +X, +Y, +Z */
      if (z_pos)
        shead1 = expand_collision_partner_list_for_neighbor(
            sv, m, mv, newsv_y + 1, &path_llf, &path_urb, shead1, R, R, R,
            x_fineparts, y_fineparts, z_fineparts, rx_hashsize, reaction_hash);

      /* go +X, +Y, -Z */
      if (z_neg)
        shead1 = expand_collision_partner_list_for_neighbor(
            sv, m, mv, newsv_y - 1, &path_llf, &path_urb, shead1, R, R, -R,
            x_fineparts, y_fineparts, z_fineparts, rx_hashsize, reaction_hash);
    }

    /* go +X, -Y */
    if (y_neg) {
      struct subvolume *newsv_y = newsv_x - (nz_parts - 1);
      shead1 = expand_collision_partner_list_for_neighbor(
          sv, m, mv, newsv_y, &path_llf, &path_urb, shead1, R, -R, 0.0,
          x_fineparts, y_fineparts, z_fineparts, rx_hashsize, reaction_hash);

      /* go +X, -Y, +Z */
      if (z_pos)
        shead1 = expand_collision_partner_list_for_neighbor(
            sv, m, mv, newsv_y + 1, &path_llf, &path_urb, shead1, R, -R, R,
            x_fineparts, y_fineparts, z_fineparts, rx_hashsize, reaction_hash);

      /* go +X, -Y, -Z */
      if (z_neg)
        shead1 = expand_collision_partner_list_for_neighbor(
            sv, m, mv, newsv_y - 1, &path_llf, &path_urb, shead1, R, -R, -R,
            x_fineparts, y_fineparts, z_fineparts, rx_hashsize, reaction_hash);
    }

    /* go +X, +Z */
    if (z_pos)
      shead1 = expand_collision_partner_list_for_neighbor(
          sv, m, mv, newsv_x + 1, &path_llf, &path_urb, shead1, R, 0.0, R,
          x_fineparts, y_fineparts, z_fineparts, rx_hashsize, reaction_hash);

    /* go +X, -Z */
    if (z_neg)
      shead1 = expand_collision_partner_list_for_neighbor(
          sv, m, mv, newsv_x - 1, &path_llf, &path_urb, shead1, R, 0.0, -R,
          x_fineparts, y_fineparts, z_fineparts, rx_hashsize, reaction_hash);
  }

  /* go -X */
  if (x_neg) {
    struct subvolume *newsv_x = sv - (nz_parts - 1) * (ny_parts - 1);
    shead1 = expand_collision_partner_list_for_neighbor(
        sv, m, mv, newsv_x, &path_llf, &path_urb, shead1, -R, 0.0, 0.0,
        x_fineparts, y_fineparts, z_fineparts, rx_hashsize, reaction_hash);

    /* go -X, +Y */
    if (y_pos) {
      struct subvolume *newsv_y = newsv_x + (nz_parts - 1);
      shead1 = expand_collision_partner_list_for_neighbor(
          sv, m, mv, newsv_y, &path_llf, &path_urb, shead1, -R, R, 0.0,
          x_fineparts, y_fineparts, z_fineparts, rx_hashsize, reaction_hash);

      /* go -X, +Y, +Z */
      if (z_pos)
        shead1 = expand_collision_partner_list_for_neighbor(
            sv, m, mv, newsv_y + 1, &path_llf, &path_urb, shead1, -R, R, R,
            x_fineparts, y_fineparts, z_fineparts, rx_hashsize, reaction_hash);

      /* go -X, +Y, -Z */
      if (z_neg)
        shead1 = expand_collision_partner_list_for_neighbor(
            sv, m, mv, newsv_y - 1, &path_llf, &path_urb, shead1, -R, R, -R,
            x_fineparts, y_fineparts, z_fineparts, rx_hashsize, reaction_hash);
    }

    /* go -X, -Y */
    if (y_neg) {
      struct subvolume *newsv_y = newsv_x - (nz_parts - 1);
      shead1 = expand_collision_partner_list_for_neighbor(
          sv, m, mv, newsv_y, &path_llf, &path_urb, shead1, -R, -R, 0.0,
          x_fineparts, y_fineparts, z_fineparts, rx_hashsize, reaction_hash);

      /* go -X, -Y, +Z */
      if (z_pos)
        shead1 = expand_collision_partner_list_for_neighbor(
            sv, m, mv, newsv_y + 1, &path_llf, &path_urb, shead1, -R, -R, R,
            x_fineparts, y_fineparts, z_fineparts, rx_hashsize, reaction_hash);

      /* go -X, -Y, -Z */
      if (z_neg)
        shead1 = expand_collision_partner_list_for_neighbor(
            sv, m, mv, newsv_y - 1, &path_llf, &path_urb, shead1, -R, -R, -R,
            x_fineparts, y_fineparts, z_fineparts, rx_hashsize, reaction_hash);
    }

    /* go -X, +Z */
    if (z_pos)
      shead1 = expand_collision_partner_list_for_neighbor(
          sv, m, mv, newsv_x + 1, &path_llf, &path_urb, shead1, -R, 0.0, R,
          x_fineparts, y_fineparts, z_fineparts, rx_hashsize, reaction_hash);

    /* go -X, -Z */
    if (z_neg)
      shead1 = expand_collision_partner_list_for_neighbor(
          sv, m, mv, newsv_x - 1, &path_llf, &path_urb, shead1, -R, 0.0, -R,
          x_fineparts, y_fineparts, z_fineparts, rx_hashsize, reaction_hash);
  }

  /* go +Y */
  if (y_pos) {
    struct subvolume *newsv_y = sv + (nz_parts - 1);
    shead1 = expand_collision_partner_list_for_neighbor(
        sv, m, mv, newsv_y, &path_llf, &path_urb, shead1, 0.0, R, 0.0,
        x_fineparts, y_fineparts, z_fineparts, rx_hashsize, reaction_hash);

    /* go +Y, +Z */
    if (z_pos)
      shead1 = expand_collision_partner_list_for_neighbor(
          sv, m, mv, newsv_y + 1, &path_llf, &path_urb, shead1, 0.0, R, R,
          x_fineparts, y_fineparts, z_fineparts, rx_hashsize, reaction_hash);

    /* go +Y, -Z */
    if (z_neg)
      shead1 = expand_collision_partner_list_for_neighbor(
          sv, m, mv, newsv_y - 1, &path_llf, &path_urb, shead1, 0.0, R, -R,
          x_fineparts, y_fineparts, z_fineparts, rx_hashsize, reaction_hash);
  }

  /* go -Y */
  if (y_pos) {
    struct subvolume *newsv_y = sv - (nz_parts - 1);
    shead1 = expand_collision_partner_list_for_neighbor(
        sv, m, mv, newsv_y, &path_llf, &path_urb, shead1, 0.0, -R, 0.0,
        x_fineparts, y_fineparts, z_fineparts, rx_hashsize, reaction_hash);

    /* go -Y, +Z */
    if (z_pos)
      shead1 = expand_collision_partner_list_for_neighbor(
          sv, m, mv, newsv_y + 1, &path_llf, &path_urb, shead1, 0.0, -R, R,
          x_fineparts, y_fineparts, z_fineparts, rx_hashsize, reaction_hash);

    /* go -Y, -Z */
    if (z_neg)
      shead1 = expand_collision_partner_list_for_neighbor(
          sv, m, mv, newsv_y - 1, &path_llf, &path_urb, shead1, 0.0, -R, -R,
          x_fineparts, y_fineparts, z_fineparts, rx_hashsize, reaction_hash);
  }

  /* go +Z */
  if (z_pos)
    shead1 = expand_collision_partner_list_for_neighbor(
        sv, m, mv, sv + 1, &path_llf, &path_urb, shead1, 0.0, 0.0, R,
        x_fineparts, y_fineparts, z_fineparts, rx_hashsize, reaction_hash);

  /* go -Z */
  if (z_neg)
    shead1 = expand_collision_partner_list_for_neighbor(
        sv, m, mv, sv - 1, &path_llf, &path_urb, shead1, 0.0, 0.0, -R,
        x_fineparts, y_fineparts, z_fineparts, rx_hashsize, reaction_hash);

  return shead1;
}

/***************************************************************************
diffuse_3D_big_list:
  In:   molecule that is moving
        maximum time we can spend diffusing
  Out:  Pointer to the molecule if it still exists (may have been
        reallocated), NULL otherwise.
        Position and time are updated, but molecule is not rescheduled.
  Note: This version takes into account both 2-way and 3-way reactions
***************************************************************************/
struct volume_molecule *diffuse_3D_big_list(struct volume *world,
                                            struct volume_molecule *m,
                                            double max_time) {
  struct vector3 displacement;  /* Molecule moves along this vector */
  struct vector3 displacement2; /* Used for 3D mol-mol unbinding */
  double disp_length;           /* length of the displacement */
  struct sp_collision *smash,
      *new_smash;             /* Thing we've hit that's under consideration */
  struct sp_collision *shead; /* Things we might hit (can interact with) */
  struct sp_collision *stail; /* tail of the collision list shead */
  struct sp_collision *shead_exp; /* Things we might hit (can interact with)
                                     from neighbor subvolumes */
  struct sp_collision *shead2; /* Things that we will hit, given our motion */

  struct sp_collision *main_shead2 =
      NULL; /* Things that we will hit, given our motion */
  struct sp_collision *new_coll = NULL;
  struct tri_collision *tentative; /* Things we already hit but haven't yet
                                      counted */
  struct tri_collision *tri_smash = NULL; /* potential trimolecular collision */
  struct tri_collision *main_tri_shead =
      NULL; /* head of the main collision list */

  struct subvolume *sv;
  struct wall *w;
  struct wall *reflectee; /* Bounced off this one, don't hit it again */
  struct rxn *rx;
  struct volume_molecule *mp, *new_mp;
  struct surface_molecule *sm = NULL;
  struct abstract_molecule *am1, *am2 = NULL;
  double steps = 1.0;
  double t_steps = 1.0;
  double rate_factor = 1.0, r_rate_factor = 1.0, factor, factor1, factor2;
  double t_start = 0; /* allows to account for the collision time after
                      reflection from a wall or moving to another subvolume */
  double f;
  double t_confident; /* We're sure we can count things up til this time */
  struct vector3 *loc_certain; /* We've counted up to this location */
  /* this flag is set to 1 only after reflection from a wall and only with
   * expanded lists. */
  int redo_expand_collision_list_flag = 0;

  int i, j, k;

  int calculate_displacement = 1;

  /* array of pointers to the possible reactions */
  struct rxn *matching_rxns[MAX_MATCHING_RXNS];
  int num_matching_rxns = 0;
  int is_reflec_flag;

  /* Flags that tell whether moving and target molecules
     can participate in the MOL_MOL_MOL, MOL_MOL, MOL_MOL_GRID,
     or MOL_GRID_GRID interactions.  Here MOL means volume molecule
     and GRID means surface molecule  */
  int moving_tri_molecular_flag = 0, moving_bi_molecular_flag = 0,
      moving_mol_mol_grid_flag = 0, moving_mol_grid_grid_flag = 0;
  /* collision flags */
  int col_tri_molecular_flag = 0, col_bi_molecular_flag = 0,
      col_mol_mol_grid_flag;

  struct species *spec = m->properties;
  struct periodic_image *periodic_box = m->periodic_box;
  if (spec == NULL)
    mcell_internal_error(
        "Attempted to take a diffusion step for a defunct molecule.");
  if (spec->space_step <= 0.0) {
    m->t += max_time;
    return m;
  }

  /* volume_reversibility and surface_reversibility routines are not valid
     in case of tri-molecular reactions */
  if (world->volume_reversibility || world->surface_reversibility) {
    if (world->volume_reversibility &&
        m->index <= DISSOCIATION_MAX) /* Only set if volume_reversibility is */
    {
      m->index = -1;
    } else if (!world->surface_reversibility) {
      if (m->flags & ACT_CLAMPED) /* Pretend we were already moving */
      {
        m->birthday -= 5 * spec->time_step; /* Pretend to be old */
      }
    }
  } else {
    if (m->flags & ACT_CLAMPED) /* Pretend we were already moving */
    {
      m->birthday -= 5 * spec->time_step; /* Pretend to be old */
    } else if ((m->flags & MATURE_MOLECULE) == 0) {
      /* Newly created particles that have long time steps gradually increase */
      /* their timestep to the full value */
      if (spec->time_step > 1.0) {
        double sched_time = convert_iterations_to_seconds(
            world->start_iterations, world->time_unit,
            world->simulation_start_seconds, m->t);
        f = 1 + 0.2 * ((sched_time - m->birthday)/world->time_unit);
        if (f < 1)
          mcell_internal_error("A %s molecule is scheduled to move before it "
                               "was born [birthday=%.15g, t=%.15g]",
                               spec->sym->name, m->birthday,
                               sched_time);
        if (max_time > f)
          max_time = f;
        if (f > m->subvol->local_storage->max_timestep)
          m->flags |= MATURE_MOLECULE;
      }
    }
  }

/* Done housekeeping, now let's do something fun! */

pretend_to_call_diffuse_3D_big_list: /* Label to allow fake recursion */

  sv = m->subvol;

  shead = NULL;
  stail = NULL;
  shead_exp = NULL;
  shead2 = NULL;

  if (calculate_displacement) {
    if (m->flags &
        ACT_CLAMPED) /* Surface clamping and microscopic reversibility */
    {
      if (m->index <= DISSOCIATION_MAX) /* Volume microscopic reversibility */
      {
        pick_release_displacement(&displacement, &displacement2,
                                  spec->space_step, world->r_step_release,
                                  world->d_step, world->radial_subdivisions,
                                  world->directions_mask, world->num_directions,
                                  world->rx_radius_3d, world->rng);

        t_steps = 0;
      } else /* Clamping or surface microscopic reversibility */
      {
        pick_clamped_displacement(&displacement, m, world->r_step_surface,
                                  world->rng, world->radial_subdivisions);
        t_steps = spec->time_step;
        m->previous_wall = NULL;
        m->index = -1;
      }
      m->flags -= ACT_CLAMPED;
      r_rate_factor = rate_factor = 1.0;
      steps = 1.0;
    } else {
      /* XXX: I don't think this is safe.  We probably need to pass in a list
       * of nearby molecules... */
      if (max_time > MULTISTEP_WORTHWHILE)
        steps = safe_diffusion_step(m, NULL, world->radial_subdivisions,
                                    world->r_step, world->x_fineparts,
                                    world->y_fineparts, world->z_fineparts);
      else
        steps = 1.0;

      t_steps = steps * spec->time_step;
      if (t_steps > max_time) {
        t_steps = max_time;
        steps = max_time / spec->time_step;
      }
      if (steps < EPS_C) {
        steps = EPS_C;
        t_steps = EPS_C * spec->time_step;
      }

      if (steps == 1.0) {
        pick_displacement(&displacement, spec->space_step, world->rng);
        r_rate_factor = rate_factor = 1.0;
      } else {
        rate_factor = sqrt(steps);
        r_rate_factor = 1.0 / rate_factor;
        pick_displacement(&displacement, rate_factor * spec->space_step,
                          world->rng);
      }
    }

    if (spec->flags & SET_MAX_STEP_LENGTH) {
      disp_length = vect_length(&displacement);
      if (disp_length > spec->max_step_length) {
        /* rescale displacement to the level of MAXIMUM_STEP_LENGTH */
        displacement.x *= (spec->max_step_length / disp_length);
        displacement.y *= (spec->max_step_length / disp_length);
        displacement.z *= (spec->max_step_length / disp_length);
      }
    }

    world->diffusion_number++;
    world->diffusion_cumtime += steps;
  }

  moving_bi_molecular_flag =
      ((spec->flags & (CAN_VOLVOL | CANT_INITIATE)) == CAN_VOLVOL);
  moving_tri_molecular_flag =
      ((spec->flags & (CAN_VOLVOLVOL | CANT_INITIATE)) == CAN_VOLVOLVOL);
  moving_mol_mol_grid_flag =
      ((spec->flags & (CAN_VOLVOLSURF | CANT_INITIATE)) == CAN_VOLVOLSURF);
  moving_mol_grid_grid_flag =
      ((spec->flags & CAN_VOLSURFSURF) == CAN_VOLSURFSURF);

  if (moving_tri_molecular_flag || moving_bi_molecular_flag ||
      moving_mol_mol_grid_flag) {
    /* scan molecules from this SV */
    struct per_species_list *psl_next, *psl,
        **psl_head = &m->subvol->species_head;
    for (psl = m->subvol->species_head; psl != NULL; psl = psl_next) {
      psl_next = psl->next;
      if (psl->properties == NULL) {
        psl_head = &psl->next;
        continue;
      }

      /* Garbage collection of empty per-species lists */
      if (psl->head == NULL) {
        *psl_head = psl->next;
        ht_remove(&sv->mol_by_species, psl);
        mem_put(sv->local_storage->pslv, psl);
        continue;
      } else
        psl_head = &psl->next;

      col_bi_molecular_flag =
          moving_bi_molecular_flag &&
          ((psl->properties->flags & CAN_VOLVOL) == CAN_VOLVOL);
      col_tri_molecular_flag =
          moving_tri_molecular_flag &&
          ((psl->properties->flags & CAN_VOLVOLVOL) == CAN_VOLVOLVOL);
      col_mol_mol_grid_flag =
          moving_mol_mol_grid_flag &&
          ((psl->properties->flags & CAN_VOLVOLSURF) == CAN_VOLVOLSURF);

      if (col_bi_molecular_flag &&
          !trigger_bimolecular_preliminary(
               world->reaction_hash, world->rx_hashsize, spec->hashval,
               psl->properties->hashval, spec, psl->properties))
        col_bi_molecular_flag = 0;

      /* What types of collisions are we concerned with for this molecule type?
       */
      int what = 0;
      if (col_bi_molecular_flag)
        what |= COLLIDE_VOL;
      if (col_tri_molecular_flag)
        what |= COLLIDE_VOL_VOL;
      if (col_mol_mol_grid_flag)
        what |= COLLIDE_VOL_SURF;

      /* If we are interested in collisions with this molecule type, add all
       * local molecules to our collision list */
      if (what != 0) {
        for (mp = psl->head; mp != NULL; mp = mp->next_v) {
          if (mp == m)
            continue;

          smash = (struct sp_collision *)CHECKED_MEM_GET(
              sv->local_storage->sp_coll, "collision data");
          smash->t = 0.0;
          smash->t_start = 0.0;
          smash->pos_start.x = m->pos.x;
          smash->pos_start.y = m->pos.y;
          smash->pos_start.z = m->pos.z;
          smash->sv_start = sv;
          smash->disp.x = displacement.x;
          smash->disp.y = displacement.y;
          smash->disp.z = displacement.z;
          smash->moving = m->properties;
          smash->target = (void *)mp;
          smash->loc.x = 0.0;
          smash->loc.y = 0.0;
          smash->loc.z = 0.0;
          smash->what = what;

          smash->next = shead;
          shead = smash;
        }
      }
    }

    if (world->use_expanded_list && shead != NULL) {
      for (stail = shead; stail->next != NULL; stail = stail->next) {
      }
    }

    if (world->use_expanded_list &&
        (moving_tri_molecular_flag || moving_bi_molecular_flag ||
         moving_mol_mol_grid_flag)) {
      shead_exp = expand_collision_partner_list(
          m, &displacement, sv, world->rx_radius_3d, world->x_fineparts,
          world->y_fineparts, world->z_fineparts, world->nx_parts,
          world->ny_parts, world->nz_parts, world->rx_hashsize,
          world->reaction_hash);

      if (stail != NULL)
        stail->next = shead_exp;
      else {
        if (shead != NULL) {
          mcell_internal_error("Collision lists corrupted. While expanding the "
                               "collision lists, expected shead to be NULL, "
                               "but it wasn't.");
        }
        shead = shead_exp;
      }
    }
  }

  reflectee = NULL;

#define TRI_CLEAN_AND_RETURN(x)                                                \
  do {                                                                         \
    if (main_tri_shead != NULL)                                                \
      mem_put_list(sv->local_storage->tri_coll, main_tri_shead);               \
    if (main_shead2 != NULL)                                                   \
      mem_put_list(sv->local_storage->sp_coll, main_shead2);                   \
    return (x);                                                                \
  } while (0)

  do {
    if (world->use_expanded_list && redo_expand_collision_list_flag) {
      /* split the combined collision list into two original lists
         and remove old "shead_exp" */
      if (shead_exp != NULL) {
        if (shead == shead_exp) {
          mem_put_list(sv->local_storage->sp_coll, shead_exp);
          shead = NULL;
        } else if (shead != NULL) {
          stail->next = NULL;
          mem_put_list(sv->local_storage->sp_coll, shead_exp);
        }
        shead_exp = NULL;
      }

      if (moving_tri_molecular_flag || moving_bi_molecular_flag ||
          moving_mol_mol_grid_flag) {
        shead_exp = expand_collision_partner_list(
            m, &displacement, sv, world->rx_radius_3d, world->x_fineparts,
            world->y_fineparts, world->z_fineparts, world->nx_parts,
            world->ny_parts, world->nz_parts, world->rx_hashsize,
            world->reaction_hash);

        /* combine two collision lists */
        if (shead_exp != NULL) {
          if (shead != NULL)
            stail->next = shead_exp;
          else
            shead = shead_exp;
        }
      }

      /* reset the flag */
      redo_expand_collision_list_flag = 0;
    }

    shead2 = ray_trace_trimol(world, m, shead, sv, &displacement, reflectee,
                              t_start);

    if (shead2 == NULL)
      mcell_internal_error("ray_trace_trimol returned NULL.");

    if (shead2->next != NULL) {
      shead2 = (struct sp_collision *)ae_list_sort(
          (struct abstract_element *)shead2);
    }

    for (smash = shead2; smash != NULL; smash = smash->next) {
      if (smash->t >= 1.0 || smash->t < 0.0) {
        if ((smash->what &
             (COLLIDE_VOL | COLLIDE_VOL_VOL | COLLIDE_VOL_SURF)) != 0)
          mcell_internal_error("Detected a mol-mol[-*] collision outside of "
                               "the 0.0...1.0 time window.  Iteration %lld, "
                               "time of collision %.8e",
                               world->current_iterations, smash->t);

        smash = NULL;
        break;
      }

      /* copy the collision objects of the type COLLIDE_VOL, COLLIDE_VOL_VOL,
         COLLIDE_VOL_SURF and COLLIDE_WALL to the main collision list */

      if (((smash->what & (COLLIDE_VOL | COLLIDE_VOL_VOL | COLLIDE_VOL_SURF)) !=
           0)) {

        new_coll = (struct sp_collision *)CHECKED_MEM_GET(
            sv->local_storage->sp_coll, "collision data");
        memcpy(new_coll, smash, sizeof(struct sp_collision));

        new_coll->t += new_coll->t_start;
        new_coll->next = main_shead2;
        main_shead2 = new_coll;

      } else if ((smash->what & COLLIDE_WALL) != 0) {
        new_coll = (struct sp_collision *)CHECKED_MEM_GET(
            sv->local_storage->sp_coll, "collision data");
        memcpy(new_coll, smash, sizeof(struct sp_collision));

        new_coll->t += new_coll->t_start;
        new_coll->next = main_shead2;
        main_shead2 = new_coll;

        /* if the wall is reflective - start over new search for collision
         * partners*/

        w = (struct wall *)smash->target;

        if ((smash->what & COLLIDE_MASK) == COLLIDE_FRONT)
          k = 1;
        else
          k = -1;

        if ((spec->flags & CAN_VOLWALL) != 0) {
          /* m->index = -1; */
          is_reflec_flag = 0;

          num_matching_rxns = trigger_intersect(
              world->reaction_hash, world->rx_hashsize, world->all_mols,
              world->all_volume_mols, world->all_surface_mols, spec->hashval,
              (struct abstract_molecule *)m, k, w, matching_rxns, 1, 1, 0);

          if (num_matching_rxns > 0) {
            for (int ii = 0; ii < num_matching_rxns; ii++) {
              rx = matching_rxns[ii];
              if (rx->n_pathways == RX_REFLEC) {
                is_reflec_flag = 1;
                break;
              }
            }
          }

          /* the wall is reflective */
          /* if there is no defined reaction between molecule
             and this particular wall, it means that this wall is
             reflective for this molecule */
          if ((num_matching_rxns == 0) || is_reflec_flag) {
            m->pos.x = smash->loc.x;
            m->pos.y = smash->loc.y;
            m->pos.z = smash->loc.z;
            m->t += t_steps * smash->t;
            reflectee = w;

            t_start += t_steps * smash->t;

            t_steps *= (1.0 - smash->t);

            factor = -2.0 * (displacement.x * w->normal.x +
                             displacement.y * w->normal.y +
                             displacement.z * w->normal.z);
            displacement.x =
                (displacement.x + factor * w->normal.x) * (1.0 - smash->t);
            displacement.y =
                (displacement.y + factor * w->normal.y) * (1.0 - smash->t);
            displacement.z =
                (displacement.z + factor * w->normal.z) * (1.0 - smash->t);

            redo_expand_collision_list_flag = 1; /* Only useful if we're using
                                                    expanded lists, but easier
                                                    to always set it */

            break;
          }
        } else {
          /* the case when (spec->flags&CAN_VOLWALL) == 0) */
          /* the default property of the wall is to be REFLECTIVE.
             It works if we do not specifically describe
             the properties of the wall */

          m->pos.x = smash->loc.x;
          m->pos.y = smash->loc.y;
          m->pos.z = smash->loc.z;
          m->t += t_steps * smash->t;
          reflectee = w;

          t_start += t_steps * smash->t;

          t_steps *= (1.0 - smash->t);

          factor = -2.0 * (displacement.x * w->normal.x +
                           displacement.y * w->normal.y +
                           displacement.z * w->normal.z);
          displacement.x =
              (displacement.x + factor * w->normal.x) * (1.0 - smash->t);
          displacement.y =
              (displacement.y + factor * w->normal.y) * (1.0 - smash->t);
          displacement.z =
              (displacement.z + factor * w->normal.z) * (1.0 - smash->t);

          redo_expand_collision_list_flag =
              1; /* Only useful if we're using expanded lists, but easier to
                    always set it */

          break;
        }
      } else if ((smash->what & COLLIDE_SUBVOL) != 0) {
        struct subvolume *nsv;

        m->pos.x = smash->loc.x;
        m->pos.y = smash->loc.y;
        m->pos.z = smash->loc.z;

        displacement.x *= (1.0 - smash->t);
        displacement.y *= (1.0 - smash->t);
        displacement.z *= (1.0 - smash->t);

        m->t += t_steps * smash->t;
        t_start += t_steps * smash->t;
        t_steps *= (1.0 - smash->t);
        if (t_steps < EPS_C)
          t_steps = EPS_C;

        nsv = traverse_subvol(sv, smash->what - COLLIDE_SV_NX - COLLIDE_SUBVOL,
                              world->ny_parts, world->nz_parts);
        if (nsv == NULL) {
          mcell_internal_error(
              "A %s molecule escaped the world at [%.2f, %.2f, %.2f]",
              spec->sym->name, m->pos.x * world->length_unit,
              m->pos.y * world->length_unit, m->pos.z * world->length_unit);
          // Never get here
          /*if (world->place_waypoints_flag && (m->flags & COUNT_ME))*/
          /*  count_region_from_scratch(world, (struct abstract_molecule *)m,*/
          /*                            NULL, -1, &(m->pos), NULL, m->t);*/
          /*spec->population--;*/
          /*collect_molecule(m);*/
          /*CLEAN_AND_RETURN(NULL);*/
        } else {
          m = migrate_volume_molecule(m, nsv);
        }

        if (shead2 != NULL) {
          mem_put_list(sv->local_storage->sp_coll, shead2);
          shead2 = NULL;
        }
        if (shead != NULL) {
          mem_put_list(sv->local_storage->sp_coll, shead);
          shead = NULL;
        }
        calculate_displacement = 0;

        if (m->properties == NULL)
          mcell_internal_error("A defunct molecule is diffusing.");

        goto pretend_to_call_diffuse_3D_big_list; /* Jump to beginning of
                                                     function */
      }
    } /* end for (smash ...) */

    if (shead2 != NULL) {
      mem_put_list(sv->local_storage->sp_coll, shead2);
      shead2 = NULL;
    }
  } while (smash != NULL);

  if (shead2 != NULL) {
    mem_put_list(sv->local_storage->sp_coll, shead2);
    shead2 = NULL;
  }
  if (shead != NULL) {
    mem_put_list(sv->local_storage->sp_coll, shead);
    shead = NULL;
  }

  for (smash = main_shead2; smash != NULL; smash = smash->next) {
    smash->t += smash->t_start;
  }

  if (main_shead2 != NULL) {
    if (main_shead2->next != NULL) {
      main_shead2 = (struct sp_collision *)ae_list_sort(
          (struct abstract_element *)main_shead2);
    }
  }

  /* build main_tri_shead list */
  for (smash = main_shead2; smash != NULL; smash = smash->next) {
    if ((smash->what & (COLLIDE_VOL | COLLIDE_VOL_VOL | COLLIDE_VOL_SURF)) !=
        0) {

      mp = (struct volume_molecule *)smash->target;

      if (moving_bi_molecular_flag && ((smash->what & COLLIDE_VOL) != 0)) {
        num_matching_rxns = trigger_bimolecular(
            world->reaction_hash, world->rx_hashsize, spec->hashval,
            mp->properties->hashval, (struct abstract_molecule *)m,
            (struct abstract_molecule *)mp, 0, 0, matching_rxns);

        if (num_matching_rxns > 0) {
          for (i = 0; i < num_matching_rxns; i++) {
            tri_smash = (struct tri_collision *)CHECKED_MEM_GET(
                sv->local_storage->tri_coll, "tri_collision data");
            tri_smash->t = smash->t;
            tri_smash->target1 = (void *)mp;
            tri_smash->target2 = NULL;
            tri_smash->orient = 0; /* default value */
            tri_smash->what = 0;
            tri_smash->what |= COLLIDE_VOL;
            tri_smash->loc = smash->loc;
            tri_smash->loc1 = smash->loc;
            tri_smash->loc2 = smash->loc;
            tri_smash->last_walk_from = smash->pos_start;
            tri_smash->intermediate = matching_rxns[i];

            tri_smash->factor = exact_disk(
                world, &(smash->loc), &(smash->disp), world->rx_radius_3d,
                smash->sv_start, m, (struct volume_molecule *)smash->target,
                world->use_expanded_list, world->x_fineparts,
                world->y_fineparts, world->z_fineparts);

            tri_smash->wall = NULL;
            tri_smash->factor *= r_rate_factor; /* scaling the reaction rate */
            tri_smash->local_prob_factor = 0;
            tri_smash->next = main_tri_shead;
            main_tri_shead = tri_smash;
          }
        }
      }
      if (moving_tri_molecular_flag && ((smash->what & COLLIDE_VOL_VOL) != 0)) {
        for (new_smash = smash->next; new_smash != NULL;
             new_smash = new_smash->next) {
          if ((new_smash->what & COLLIDE_VOL_VOL) == 0)
            continue;

          new_mp = (struct volume_molecule *)new_smash->target;

          num_matching_rxns = trigger_trimolecular(
              world->reaction_hash, world->rx_hashsize, smash->moving->hashval,
              mp->properties->hashval, new_mp->properties->hashval,
              smash->moving, mp->properties, new_mp->properties, 0, 0, 0,
              matching_rxns);

          if (num_matching_rxns > 0) {
            for (i = 0; i < num_matching_rxns; i++) {
              tri_smash = (struct tri_collision *)CHECKED_MEM_GET(
                  sv->local_storage->tri_coll, "collision data");
              tri_smash->loc = new_smash->loc;
              tri_smash->t = new_smash->t;
              tri_smash->target2 = (void *)new_mp;
              tri_smash->loc2 = new_smash->loc;
              tri_smash->target1 = (void *)mp;
              tri_smash->loc1 = smash->loc;
              tri_smash->last_walk_from = new_smash->pos_start;
              tri_smash->orient = 0; /* default value */

              factor1 = exact_disk(world, &(smash->loc), &(smash->disp),
                                   world->rx_radius_3d, smash->sv_start, m,
                                   (struct volume_molecule *)smash->target,
                                   world->use_expanded_list, world->x_fineparts,
                                   world->y_fineparts, world->z_fineparts);
              factor2 = exact_disk(world, &(new_smash->loc), &(new_smash->disp),
                                   world->rx_radius_3d, new_smash->sv_start, m,
                                   (struct volume_molecule *)new_smash->target,
                                   world->use_expanded_list, world->x_fineparts,
                                   world->y_fineparts, world->z_fineparts);
              tri_smash->factor = factor1 * factor2;
              tri_smash->factor *=
                  r_rate_factor; /* scaling the reaction rate */
              tri_smash->local_prob_factor = 0;
              tri_smash->what = 0;
              tri_smash->what |= COLLIDE_VOL_VOL;
              tri_smash->intermediate = matching_rxns[i];
              tri_smash->wall = NULL;
              tri_smash->next = main_tri_shead;
              main_tri_shead = tri_smash;
            }

          } /* end if (...) */

        } /* end for (new_smash...) */
      }   /* end if (...) */
      if (moving_mol_mol_grid_flag && ((smash->what & COLLIDE_VOL_SURF) != 0)) {
        for (new_smash = smash->next; new_smash != NULL;
             new_smash = new_smash->next) {
          if ((new_smash->what & COLLIDE_WALL) == 0)
            continue;

          w = (struct wall *)new_smash->target;

          if ((new_smash->what & COLLIDE_MASK) == COLLIDE_FRONT)
            k = 1;
          else
            k = -1;

          if (w->grid != NULL) {
            j = xyz2grid(&(new_smash->loc), w->grid);
            if (w->grid->sm_list[j] && w->grid->sm_list[j]->sm) {
              if (m->index != j || m->previous_wall != w) {
                sm = w->grid->sm_list[j]->sm;
                num_matching_rxns = trigger_trimolecular(
                    world->reaction_hash, world->rx_hashsize,
                    smash->moving->hashval, mp->properties->hashval,
                    sm->properties->hashval, smash->moving, mp->properties,
                    sm->properties, k, k, sm->orient, matching_rxns);

                if (num_matching_rxns > 0) {
                  for (i = 0; i < num_matching_rxns; i++) {
                    tri_smash = (struct tri_collision *)CHECKED_MEM_GET(
                        sv->local_storage->tri_coll, "tri_collision data");
                    tri_smash->t = new_smash->t;
                    tri_smash->target1 = (void *)mp;
                    tri_smash->target2 = (void *)sm;
                    tri_smash->orient = k;
                    tri_smash->what = 0;
                    tri_smash->what |= COLLIDE_VOL_SURF;
                    tri_smash->loc = new_smash->loc;
                    tri_smash->loc1 = smash->loc;
                    tri_smash->loc2 = new_smash->loc;
                    tri_smash->last_walk_from = new_smash->pos_start;
                    tri_smash->intermediate = matching_rxns[i];

                    factor1 =
                        exact_disk(world, &(smash->loc), &(smash->disp),
                                   world->rx_radius_3d, smash->sv_start, m,
                                   (struct volume_molecule *)smash->target,
                                   world->use_expanded_list, world->x_fineparts,
                                   world->y_fineparts, world->z_fineparts);

                    factor2 = r_rate_factor / w->grid->binding_factor;
                    tri_smash->factor = factor1 * factor2;

                    tri_smash->local_prob_factor = 0;
                    tri_smash->wall = w;

                    tri_smash->next = main_tri_shead;
                    main_tri_shead = tri_smash;
                  }
                } /* end if (num_matching_rxns > 0) */
              }
              /* Matched previous wall and index--don't rebind */
              else {
                m->index = -1; // Avoided rebinding, but next time it's OK
              }

            } /* end if (w->grid->mol[j] ...) */
          }   /* end if (w->grid != NULL ...) */

        } /* end for (new_smash...) */
      }   /* end if (...) */

    } /* end if (...) */
    else if ((smash->what & COLLIDE_WALL) != 0) {
      w = (struct wall *)smash->target;
      int wall_was_accounted_for = 0; /* flag */

      if ((smash->what & COLLIDE_MASK) == COLLIDE_FRONT)
        k = 1;
      else
        k = -1;

      /* first look for the bimolecular reactions between moving and
         surface molecules */
      if (w->grid != NULL && (spec->flags & CAN_VOLSURF) != 0) {
        j = xyz2grid(&(smash->loc), w->grid);
        if (w->grid->sm_list[j] && w->grid->sm_list[j]->sm) {
          if (m->index != j || m->previous_wall != w) {
            sm = w->grid->sm_list[j]->sm;
            // look for bimolecular reactions between volume and surface mols
            num_matching_rxns = trigger_bimolecular(
                world->reaction_hash, world->rx_hashsize, spec->hashval,
                sm->properties->hashval, (struct abstract_molecule *)m,
                (struct abstract_molecule *)sm, k, sm->orient, matching_rxns);
            if (num_matching_rxns > 0) {
              for (i = 0; i < num_matching_rxns; i++) {
                tri_smash = (struct tri_collision *)CHECKED_MEM_GET(
                    sv->local_storage->tri_coll, "collision data");
                tri_smash->t = smash->t;
                tri_smash->target1 = (void *)sm;
                tri_smash->target2 = NULL;
                tri_smash->orient = k;
                tri_smash->what = 0;
                tri_smash->what |= COLLIDE_SURF;
                tri_smash->loc = smash->loc;
                tri_smash->loc1 = smash->loc;
                tri_smash->loc2 = smash->loc;
                tri_smash->last_walk_from = smash->pos_start;
                tri_smash->intermediate = matching_rxns[i];
                tri_smash->factor = r_rate_factor / w->grid->binding_factor;
                tri_smash->local_prob_factor = 0;
                tri_smash->wall = w;
                tri_smash->next = main_tri_shead;
                main_tri_shead = tri_smash;
                wall_was_accounted_for = 1;
              }
            } /* end if (num_matching_rxns > 0) */
          } else {
            m->index = -1; /* Avoided rebinding, but next time it's OK */
          }
        } /* end if (w->grid->mol[j] ...) */
      }   /* end if (w->grid != NULL ...) */

      /* now look for the trimolecular reactions */
      if (moving_mol_grid_grid_flag) {
        if (w->grid != NULL) {
          j = xyz2grid(&(smash->loc), w->grid);
          if (w->grid->sm_list[j] && w->grid->sm_list[j]->sm) {
            sm = w->grid->sm_list[j]->sm;
            if (m->index != j || m->previous_wall != w) {
              /* search for neighbors that can participate
                in 3-way reaction */
              struct surface_molecule *smp; /* Neighboring molecules */
              struct tile_neighbor *tile_nbr_head = NULL, *curr;
              int list_length = 0;

              /* find neighbor molecules to react with */
              find_neighbor_tiles(world, sm, sm->grid, sm->grid_index, 0, 1,
                                  &tile_nbr_head, &list_length);
              if (tile_nbr_head != NULL) {
                double local_prob_factor; /*local probability factor for the
                                             reaction */
                local_prob_factor = 3.0 / list_length;

                /* step through the neighbors */
                for (curr = tile_nbr_head; curr != NULL; curr = curr->next) {
                  struct surface_molecule_list *sm_list = curr->grid->sm_list[curr->idx]; 
                  if (sm_list == NULL || sm_list->sm == NULL)
                    continue;
                  smp = curr->grid->sm_list[curr->idx]->sm;

                  /* check whether any of potential partners
                are behind restrictive (REFLECTIVE/ABSORPTIVE) boundary */
                  if ((sm->properties->flags & CAN_REGION_BORDER) ||
                      (smp->properties->flags & CAN_REGION_BORDER)) {
                    if (sm->grid->surface != smp->grid->surface) {
                      /* INSIDE-OUT check */
                      if (walls_belong_to_at_least_one_different_restricted_region(
                              world, sm->grid->surface, sm,
                              smp->grid->surface, smp))
                        continue;

                      /* OUTSIDE-IN check */
                      if (walls_belong_to_at_least_one_different_restricted_region(
                              world, sm->grid->surface, smp,
                              smp->grid->surface, sm))
                        continue;
                    }
                  }

                  num_matching_rxns = trigger_trimolecular(
                      world->reaction_hash, world->rx_hashsize,
                      smash->moving->hashval, sm->properties->hashval,
                      smp->properties->hashval, smash->moving, sm->properties,
                      smp->properties, k, sm->orient, smp->orient,
                      matching_rxns);
                  if (num_matching_rxns > 0) {
                    for (i = 0; i < num_matching_rxns; i++) {
                      tri_smash = (struct tri_collision *)CHECKED_MEM_GET(
                          sv->local_storage->tri_coll, "collision data");
                      tri_smash->t = smash->t;
                      tri_smash->target1 = (void *)sm;
                      tri_smash->target2 = (void *)smp;
                      tri_smash->orient = k;
                      tri_smash->what = 0;
                      tri_smash->what |= COLLIDE_SURF_SURF;
                      grid2xyz(curr->grid, curr->idx, &(tri_smash->loc));
                      tri_smash->loc1 = smash->loc;
                      tri_smash->loc2 = tri_smash->loc;
                      tri_smash->last_walk_from = smash->pos_start;
                      tri_smash->intermediate = matching_rxns[i];
                      tri_smash->factor = r_rate_factor /
                                          (w->grid->binding_factor) *
                                          (curr->grid->binding_factor);
                      tri_smash->local_prob_factor = local_prob_factor;
                      tri_smash->wall = w;
                      tri_smash->next = main_tri_shead;
                      main_tri_shead = tri_smash;

                      wall_was_accounted_for = 1;
                    }
                  }
                }
                if (tile_nbr_head != NULL)
                  delete_tile_neighbor_list(tile_nbr_head);
              }
            }
          }
        }

      } /* end if (moving_mol_grid_grid_flag) */

      /* now look for the mol-wall interactions */
      if ((spec->flags & CAN_VOLWALL) != 0) {

        /*  m->index = -1;  */
        num_matching_rxns = trigger_intersect(
            world->reaction_hash, world->rx_hashsize, world->all_mols,
            world->all_volume_mols, world->all_surface_mols, spec->hashval,
            (struct abstract_molecule *)m, k, w, matching_rxns, 1, 1, 0);

        for (i = 0; i < num_matching_rxns; i++) {
          rx = matching_rxns[i];
          tri_smash = (struct tri_collision *)CHECKED_MEM_GET(
              sv->local_storage->tri_coll, "tri_collision data");
          tri_smash->t = smash->t;
          tri_smash->target1 = (void *)w;
          tri_smash->target2 = NULL;
          tri_smash->orient = k;
          tri_smash->what = 0;
          tri_smash->what |= COLLIDE_WALL;
          tri_smash->loc = smash->loc;
          tri_smash->loc1 = smash->loc;
          tri_smash->loc2 = smash->loc;
          tri_smash->last_walk_from = smash->pos_start;
          tri_smash->intermediate = rx;
          tri_smash->factor = r_rate_factor;
          tri_smash->local_prob_factor = 0;
          tri_smash->wall = w;

          tri_smash->next = main_tri_shead;
          main_tri_shead = tri_smash;

          wall_was_accounted_for = 1;
        } /* end for (i = 0; i < num_matching_rxns; ...) */
      }   /* end if (spec->flags & CAN_WALLMOL ...) */

      if (!wall_was_accounted_for) {

        /* This is a simple reflective wall
           (default wall behavior).
           We want to keep it in the "tri_smash"
           list just in order to account for the hits with it */
        tri_smash = (struct tri_collision *)CHECKED_MEM_GET(
            sv->local_storage->tri_coll, "tri_collision data");
        tri_smash->t = smash->t;
        tri_smash->target1 = (void *)w;
        tri_smash->target2 = NULL;
        tri_smash->orient = k;
        tri_smash->what = 0;
        tri_smash->what |= COLLIDE_WALL;
        tri_smash->loc = smash->loc;
        tri_smash->loc1 = smash->loc;
        tri_smash->loc2 = smash->loc;
        tri_smash->last_walk_from = smash->pos_start;
        tri_smash->intermediate = NULL;
        tri_smash->factor = r_rate_factor;
        tri_smash->local_prob_factor = 0;
        tri_smash->wall = w;

        tri_smash->next = main_tri_shead;
        main_tri_shead = tri_smash;
      }

    } /* end if (smash->what & COLLIDE_WALL)... */
  }   /* end for (smash ...) */

  if (main_tri_shead != NULL) {
    if (main_tri_shead->next != NULL) {
      main_tri_shead = (struct tri_collision *)ae_list_sort(
          (struct abstract_element *)main_tri_shead);
    }
  }

  tentative = main_tri_shead;
  loc_certain = NULL;

  /* now check for the reactions going through the 'main_tri_shead' list */
  for (tri_smash = main_tri_shead; tri_smash != NULL;
       tri_smash = tri_smash->next) {

    if (world->notify->molecule_collision_report == NOTIFY_FULL) {
      if (((tri_smash->what & COLLIDE_VOL) != 0) &&
          (world->rxn_flags.vol_vol_reaction_flag)) {
        world->vol_vol_colls++;
      } else if (((tri_smash->what & COLLIDE_SURF) != 0) &&
                 (world->rxn_flags.vol_surf_reaction_flag)) {
        world->vol_surf_colls++;
      } else if (((tri_smash->what & COLLIDE_VOL_VOL) != 0) &&
                 (world->rxn_flags.vol_vol_vol_reaction_flag)) {
        world->vol_vol_vol_colls++;
      } else if (((tri_smash->what & COLLIDE_VOL_SURF) != 0) &&
                 (world->rxn_flags.vol_vol_surf_reaction_flag)) {
        world->vol_vol_surf_colls++;
      } else if (((tri_smash->what & COLLIDE_SURF_SURF) != 0) &&
                 (world->rxn_flags.vol_surf_surf_reaction_flag)) {
        world->vol_surf_surf_colls++;
      }
    }

    j = INT_MIN;

    if (((tri_smash->what & (COLLIDE_VOL | COLLIDE_SURF | COLLIDE_VOL_VOL |
                             COLLIDE_VOL_SURF | COLLIDE_SURF_SURF)) != 0)) {

      if (tri_smash->t < EPS_C)
        continue;

      if ((tri_smash->factor < 0)) /* one of the targets is blocked by a wall */
        continue;                  /* Reaction blocked by a wall */

      am1 = (struct abstract_molecule *)tri_smash->target1;
      if (tri_smash->target2 != NULL) {
        am2 = (struct abstract_molecule *)tri_smash->target2;
      } else
        am2 = NULL;

      /* if one of the targets was already destroyed
         - move on  */
      if (am1 != NULL && am1->properties == NULL)
        continue;
      if (am2 != NULL && am2->properties == NULL)
        continue;

      rx = tri_smash->intermediate;

      k = tri_smash->orient;

      if ((rx != NULL) && (rx->prob_t != NULL))
        update_probs(world, rx, m->t);

      /* XXX: Change required here to support macromol+trimol */
      i = test_bimolecular(rx, tri_smash->factor, tri_smash->local_prob_factor,
                           NULL, NULL,
                           world->rng);

      if (i < RX_LEAST_VALID_PATHWAY)
        continue;

      if ((tri_smash->what & COLLIDE_VOL) != 0) {
        j = outcome_bimolecular(world, rx, i, (struct abstract_molecule *)m,
                                am1, 0, 0, m->t + tri_smash->t,
                                &(tri_smash->loc), loc_certain);
      } else if ((tri_smash->what & COLLIDE_SURF) != 0) {
        j = outcome_bimolecular(
            world, rx, i, (struct abstract_molecule *)m, am1, k,
            ((struct surface_molecule *)am1)->orient, m->t + tri_smash->t,
            &(tri_smash->loc), &(tri_smash->last_walk_from));
      } else if ((tri_smash->what & COLLIDE_VOL_VOL) != 0) {
        j = outcome_trimolecular(world, rx, i, (struct abstract_molecule *)m,
                                 am1, am2, 0, 0, 0, m->t + tri_smash->t,
                                 &(tri_smash->loc),
                                 &(tri_smash->last_walk_from));
      } else if ((tri_smash->what & COLLIDE_VOL_SURF) != 0) {
        short orient_target = 0;
        if ((am1->properties->flags & ON_GRID) != 0) {
          orient_target = ((struct surface_molecule *)am1)->orient;

        } else {
          orient_target = ((struct surface_molecule *)am2)->orient;
        }

        j = outcome_trimolecular(world, rx, i, (struct abstract_molecule *)m,
                                 am1, am2, k, k, orient_target,
                                 m->t + tri_smash->t, &(tri_smash->loc),
                                 &tri_smash->last_walk_from);

      } else if ((tri_smash->what & COLLIDE_SURF_SURF) != 0) {
        short orient1, orient2;
        orient1 = ((struct surface_molecule *)am1)->orient;
        orient2 = ((struct surface_molecule *)am2)->orient;

        j = outcome_trimolecular(world, rx, i, (struct abstract_molecule *)m,
                                 am1, am2, k, orient1, orient2,
                                 m->t + tri_smash->t, &(tri_smash->loc),
                                 &tri_smash->last_walk_from);
      }

      if (j != RX_DESTROY)
        continue;
      else {

        /* Count the hits up until we were destroyed */
        for (; tentative != NULL && tentative->t <= tri_smash->t;
             tentative = tentative->next) {
          if (tentative->wall == NULL)
            continue;
          if (!(spec->flags & (tentative->wall->flags) & COUNT_SOME_MASK))
            continue;
          count_region_update(world, spec, periodic_box,
            tentative->wall->counting_regions, tentative->orient, 0,
            &(tentative->loc), tentative->t);
          if (tentative == tri_smash)
            break;
        }

        TRI_CLEAN_AND_RETURN(NULL);
      }

    }

    if ((tri_smash->what & COLLIDE_WALL) != 0) {
      k = tri_smash->orient;

      w = (struct wall *)tri_smash->target1;

      if (tri_smash->next == NULL)
        t_confident = tri_smash->t;
      else if (tri_smash->next->t * (1.0 - EPS_C) > tri_smash->t)
        t_confident = tri_smash->t;
      else
        t_confident = tri_smash->t * (1.0 - EPS_C);

      if ((spec->flags & CAN_VOLWALL) != 0) {
        rx = tri_smash->intermediate;

        if (rx != NULL) {
          if ((rx->n_pathways > RX_SPECIAL) &&
              (world->notify->molecule_collision_report == NOTIFY_FULL)) {
            if (world->rxn_flags.vol_wall_reaction_flag)
              world->vol_wall_colls++;
          }

          if (rx->n_pathways == RX_TRANSP) {
            rx->n_occurred++;
            if ((m->flags & COUNT_ME) != 0 &&
                (spec->flags & COUNT_SOME_MASK) != 0) {
              /* Count as far up as we can unambiguously */
              for (; tentative != NULL && tentative->t <= t_confident;
                   tentative = tentative->next) {
                if (tentative->wall == NULL)
                  continue;
                if (!(spec->flags & (tentative->wall->flags) & COUNT_SOME_MASK))
                  continue;
                count_region_update(world, spec, periodic_box,
                  tentative->wall->counting_regions, tentative->orient, 1,
                  &(tentative->loc), tentative->t);
                if (tentative == tri_smash)
                  break;
              }
            }

            continue; /* Ignore this wall and keep going */
          } else if (rx->n_pathways != RX_REFLEC) {
            if (rx->prob_t != NULL)
              update_probs(world, rx, m->t);
            i = test_intersect(rx, r_rate_factor, world->rng);
            if (i > RX_NO_RX) {
              /* Save m flags in case it gets collected in outcome_intersect */
              int mflags = m->flags;
              j = outcome_intersect(
                  world, rx, i, w, (struct abstract_molecule *)m, k,
                  m->t + t_steps * tri_smash->t, &(tri_smash->loc), NULL);

              if (j == RX_FLIP) {
                if ((m->flags & COUNT_ME) != 0 &&
                    (spec->flags & COUNT_SOME_MASK) != 0) {
                  /* Count as far up as we can unambiguously */
                  for (; tentative != NULL && tentative->t <= t_confident;
                       tentative = tentative->next) {
                    if (tentative->wall == NULL)
                      continue;
                    if (!(spec->flags & (tentative->wall->flags) &
                          COUNT_SOME_MASK))
                      continue;
                    count_region_update(world, spec, periodic_box,
                      tentative->wall->counting_regions, tentative->orient, 1,
                      &(tentative->loc), tentative->t);
                    if (tentative == tri_smash)
                      break;
                  }
                }

                continue; /* pass through */
              } else if (j == RX_DESTROY) {
                if ((mflags & COUNT_ME) != 0 &&
                    (spec->flags & COUNT_HITS) != 0) {
                  /* Count the hits up until we were destroyed */
                  for (; tentative != NULL && tentative->t <= tri_smash->t;
                       tentative = tentative->next) {
                    if (tentative->wall == NULL)
                      continue;
                    if (!(spec->flags & (tentative->wall->flags) &
                          COUNT_SOME_MASK))
                      continue;
                    count_region_update(world, spec, periodic_box,
                      tentative->wall->counting_regions, tentative->orient, 0,
                      &(tentative->loc), tentative->t);
                    if (tentative == tri_smash)
                      break;
                  }
                }

                TRI_CLEAN_AND_RETURN(NULL);
              }
            }
          } else if (rx->n_pathways == RX_REFLEC) {
            /* We reflected, so we hit but did not cross things we tentatively
             * hit earlier */
            if ((m->flags & COUNT_ME) != 0 && (spec->flags & COUNT_HITS) != 0) {
              for (; tentative != NULL && tentative->t <= tri_smash->t;
                   tentative = tentative->next) {
                if (tentative->wall == NULL)
                  continue;
                if (!(spec->flags & (tentative->wall->flags) & COUNT_SOME_MASK))
                  continue;
                count_region_update(world, spec, periodic_box,
                  tentative->wall->counting_regions, tentative->orient, 0,
                  &(tentative->loc), tentative->t);
                if (tentative == tri_smash)
                  break;
              }
            }
            continue;
          }
        } else { /* (rx == NULL) - simple reflective wall */
          if ((m->flags & COUNT_ME) != 0 && (spec->flags & COUNT_HITS) != 0) {
            for (; tentative != NULL && tentative->t <= tri_smash->t;
                 tentative = tentative->next) {
              if (tentative->wall == NULL)
                continue;
              if (!(spec->flags & (tentative->wall->flags) & COUNT_SOME_MASK))
                continue;
              count_region_update(world, spec, periodic_box,
                tentative->wall->counting_regions, tentative->orient, 0,
                &(tentative->loc), tentative->t);
              if (tentative == tri_smash)
                break;
            }
          }
          continue;
        }
      } /* end if (spec->flags & CAN_WALLMOL ...) */

    } /* end if ((tri_smash->what & COLLIDE_WALL) ... */

  } /* end for (tri_smash ...) */

#undef TRI_CLEAN_AND_RETURN

  m->pos.x += displacement.x;
  m->pos.y += displacement.y;
  m->pos.z += displacement.z;
  m->t += t_steps;

  m->index = -1;
  m->previous_wall = NULL;

  if (main_tri_shead != NULL)
    mem_put_list(sv->local_storage->tri_coll, main_tri_shead);
  if (main_shead2 != NULL)
    mem_put_list(sv->local_storage->sp_coll, main_shead2);

  return m;
}

/*************************************************************************
react_2D_trimol_all_neighbors:
  In: molecule that may react
      maximum duration we have to react
  Out: Pointer to the molecule if it still exists (may have been
       destroyed), NULL otherwise.
  Note: Time is not updated--assume that's already taken care of
        elsewhere.  Only nearest neighbors can react.
  PostNote: This function is valid only for the trimolecular reaction
            involving all three neighbor surface molecules whose tiles are
            connected by edge or vertex
*************************************************************************/
struct surface_molecule *react_2D_trimol_all_neighbors(
    struct volume *world, struct surface_molecule *sm, double t,
    enum notify_level_t molecule_collision_report,
    enum notify_level_t final_summary, int grid_grid_grid_reaction_flag,
    long long *surf_surf_surf_colls) {

  struct surface_molecule *gm_f, *gm_s; /* Neighboring molecule */

  int i;     /* points to the pathway of the reaction */
  int j;     /* points to the the reaction */
  int n = 0; /* total number of possible reactions for a given molecules
                with all its neighbors */
  int k; /* return value from "outcome_trimolecular()" */
  int l = 0, jj, kk;
  int num_matching_rxns = 0;
  struct rxn *matching_rxns[MAX_MATCHING_RXNS];

  /* linked lists of the tile neighbors (first and second level) */
  struct tile_neighbor *tile_nbr_head_f = NULL, *tile_nbr_head_s = NULL,
                       *curr_f, *curr_s;
  int list_length_f, list_length_s; /* length of the linked lists above */

  int max_size = 12 * 12 * MAX_MATCHING_RXNS; /* reasonable assumption */
  struct rxn *rxn_array[max_size]; /* array of reaction objects with neighbor
                                     molecules */
  /* local probability factors for the reactions */
  double local_prob_factor_f, local_prob_factor_s;
  double local_prob_factor[max_size];
  double cf[max_size]; /* Correction factors for area for those molecules */
  /* points to the first partner in the trimol reaction */
  struct surface_molecule *first_partner[max_size];
  /* points to the second partner in the trimol reaction */
  struct surface_molecule *second_partner[max_size];

  for (kk = 0; kk < max_size; kk++) {
    rxn_array[kk] = NULL;
    first_partner[kk] = NULL;
    second_partner[kk] = NULL;
    cf[kk] = 0;
    local_prob_factor[kk] = 0;
  }

  /* find first level neighbor molecules to react with */
  find_neighbor_tiles(world, sm, sm->grid, sm->grid_index, 0, 1,
                      &tile_nbr_head_f, &list_length_f);

  if (tile_nbr_head_f == NULL)
    return sm;

  /* Calculate local_prob_factor for the reaction probability.
     Here we convert from 3 neighbor tiles (upper probability
     limit) to the real number of neighbor tiles. */
  local_prob_factor_f = 1.0 / list_length_f;

  /* step through the neighbors */
  for (curr_f = tile_nbr_head_f; curr_f != NULL; curr_f = curr_f->next) {
    struct surface_molecule_list *sm_list = curr_f->grid->sm_list[curr_f->idx]; 
    if (sm_list == NULL || sm_list->sm == NULL)
      continue;
    gm_f = sm_list->sm;

    /* check whether the neighbor molecule is behind
       the restrictive region boundary   */
    if ((sm->properties->flags & CAN_REGION_BORDER) ||
        (gm_f->properties->flags & CAN_REGION_BORDER)) {
      if (sm->grid->surface != gm_f->grid->surface) {
        /* INSIDE-OUT check */
        if (walls_belong_to_at_least_one_different_restricted_region(
                world, sm->grid->surface, sm, gm_f->grid->surface, gm_f))
          continue;
        /* OUTSIDE-IN check */
        if (walls_belong_to_at_least_one_different_restricted_region(
                world, sm->grid->surface, gm_f, gm_f->grid->surface, sm))
          continue;
      }
    }

    /* find nearest neighbor molecules to react with (2nd level) */
    find_neighbor_tiles(world, gm_f, gm_f->grid, gm_f->grid_index, 0, 1,
                        &tile_nbr_head_s, &list_length_s);

    if (tile_nbr_head_s == NULL)
      continue;

    local_prob_factor_s = 1.0 / (list_length_s - 1);

    for (curr_s = tile_nbr_head_s; curr_s != NULL; curr_s = curr_s->next) {
      sm_list = curr_s->grid->sm_list[curr_s->idx]; 
      if (sm_list == NULL || sm_list->sm == NULL)
        continue;
      gm_s = curr_s->grid->sm_list[curr_s->idx]->sm;
      if (gm_s == NULL)
        continue;
      if (gm_s == gm_f)
        continue; /* no self reaction for trimolecular reaction */
      if (gm_s == sm)
        continue;

      /* Check whether there are restrictive region boundaries between "sm" and
       * "gm_s". By now we know that there are no restrictive region boundaries
       * between "sm" and "gm_f". */

      if ((sm->properties->flags & CAN_REGION_BORDER) ||
          (gm_f->properties->flags & CAN_REGION_BORDER) ||
          (gm_s->properties->flags & CAN_REGION_BORDER)) {
        if (gm_f->grid->surface != gm_s->grid->surface) {
          /* INSIDE-OUT check */
          if (walls_belong_to_at_least_one_different_restricted_region(
                  world, gm_f->grid->surface, gm_f, gm_s->grid->surface, gm_s))
            continue;
          /* OUTSIDE-IN check */
          if (walls_belong_to_at_least_one_different_restricted_region(
                  world, gm_f->grid->surface, gm_s, gm_s->grid->surface, gm_f))
            continue;
        }
        if (sm->grid->surface != gm_s->grid->surface) {
          /* INSIDE-OUT check */
          if (walls_belong_to_at_least_one_different_restricted_region(
                  world, sm->grid->surface, sm, gm_s->grid->surface, gm_s))
            continue;
          /* OUTSIDE-IN check */
          if (walls_belong_to_at_least_one_different_restricted_region(
                  world, sm->grid->surface, gm_s, gm_s->grid->surface, sm))
            continue;
        }
      }

      num_matching_rxns = trigger_trimolecular(
          world->reaction_hash, world->rx_hashsize, sm->properties->hashval,
          gm_f->properties->hashval, gm_s->properties->hashval, sm->properties,
          gm_f->properties, gm_s->properties, sm->orient, gm_f->orient,
          gm_s->orient, matching_rxns);
      if (num_matching_rxns > 0) {
        if ((final_summary == NOTIFY_FULL) &&
            (molecule_collision_report == NOTIFY_FULL)) {
          if (grid_grid_grid_reaction_flag)
            surf_surf_surf_colls++;
        }
        for (jj = 0; jj < num_matching_rxns; jj++) {
          if (matching_rxns[jj] != NULL) {
            rxn_array[l] = matching_rxns[jj];
            cf[l] = (t / (gm_f->grid->binding_factor)) *
                    (t / (gm_s->grid->binding_factor));
            local_prob_factor[l] = local_prob_factor_f * local_prob_factor_s;

            first_partner[l] = gm_f;
            second_partner[l] = gm_s;
            l++;
          }
        }

        n += num_matching_rxns;
      }
    }
    if (tile_nbr_head_s != NULL)
      delete_tile_neighbor_list(tile_nbr_head_s);
  }

  if (tile_nbr_head_f != NULL)
    delete_tile_neighbor_list(tile_nbr_head_f);

  if (n > max_size)
    mcell_internal_error("The size of the reactions array in the function "
                         "'react_2D_trimol_all_neighbors()' is not "
                         "sufficient.");

  if (n == 0) {
    return sm; /* Nobody to react with */
  } else if (n == 1) {
    /* XXX: Change required here to support macromol+trimol */
    i = test_bimolecular(rxn_array[0], cf[0], local_prob_factor[0], NULL, NULL, world->rng);
    j = 0;
  } else {
    /* XXX: Change required here to support macromol+trimol */

    j = test_many_reactions_all_neighbors(rxn_array, cf, local_prob_factor, n,
                                          &(i), world->rng);
  }

  if ((j == RX_NO_RX) || (i < RX_LEAST_VALID_PATHWAY)) {
    return sm; /* No reaction */
  }

  /* run the reaction */
  k = outcome_trimolecular(
      world, rxn_array[j], i, (struct abstract_molecule *)sm,
      (struct abstract_molecule *)first_partner[j],
      (struct abstract_molecule *)second_partner[j], sm->orient,
      first_partner[j]->orient, second_partner[j]->orient, sm->t, NULL, NULL);

  if (k == RX_DESTROY) {
    mem_put(sm->birthplace, sm);
    return NULL;
  }

  return sm;
}
