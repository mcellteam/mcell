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

/**************************************************************************\
** File: vol_util.c                                                       **
**                                                                        **
** Purpose: Adds, subtracts, and moves particles around (bookkeeping).    **
**                                                                        **
** Testing status: compiles.  Worked earlier, but has been changed.       **
\**************************************************************************/

#include "config.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "diffuse.h"
#include "vector.h"
#include "logging.h"
#include "rng.h"
#include "mem_util.h"
#include "count_util.h"
#include "vol_util.h"
#include "react.h"
#include "wall_util.h"
#include "grid_util.h"
#include "diffuse.h"

static int test_max_release(double num_to_release, char *name);

static int check_release_probability(double release_prob, struct volume *state,
                                     struct release_event_queue *req,
                                     struct release_pattern *rpat);

static int skip_past_events(double release_prob, struct volume *state,
                            struct release_event_queue *req,
                            struct release_pattern *rpat);

static int calculate_number_to_release(struct release_site_obj *rso,
                                       struct volume *state);

static int release_inside_regions(struct volume *state,
                                  struct release_site_obj *rso,
                                  struct volume_molecule *vm, int n);

static int num_vol_mols_from_conc(struct release_site_obj *rso,
                                  double length_unit, bool *exactNumber);

/*************************************************************************
inside_subvolume:
  In: pointer to vector3
      pointer to subvolume
  Out: nonzero if the vector is inside the subvolume.
*************************************************************************/
int inside_subvolume(struct vector3 *point, struct subvolume *subvol,
                     double *x_fineparts, double *y_fineparts,
                     double *z_fineparts) {
  return ((point->x >= x_fineparts[subvol->llf.x]) &&
          (point->x <= x_fineparts[subvol->urb.x]) &&
          (point->y >= y_fineparts[subvol->llf.y]) &&
          (point->y <= y_fineparts[subvol->urb.y]) &&
          (point->z >= z_fineparts[subvol->llf.z]) &&
          (point->z <= z_fineparts[subvol->urb.z]));
}

/*************************************************************************
find_coarse_subvolume:
  In: pointer to vector3
  Out: pointer to the coarse subvolume that the vector is within
*************************************************************************/
struct subvolume *find_coarse_subvol(struct volume *state,
                                     struct vector3 *loc) {
  int i = bisect(state->x_partitions, state->nx_parts, loc->x);
  int j = bisect(state->y_partitions, state->ny_parts, loc->y);
  int k = bisect(state->z_partitions, state->nz_parts, loc->z);
  return &(state->subvol
               [k + (state->nz_parts - 1) * (j + (state->ny_parts - 1) * i)]);
}

/*************************************************************************
traverse_subvol:
  In: pointer to our current subvolume
      pointer to a vector3 of where we want to be
      which direction we're traveling to get there
  Out: subvolume that's closest to where we want to be in our direction
  Note: BSP trees traverse is not yet implemented
*************************************************************************/
struct subvolume *traverse_subvol(struct subvolume *here,
                                  int which, int ny_parts, int nz_parts) {
  switch (which) {
  case X_NEG:
    if (here->world_edge & X_NEG_BIT)
      return NULL;
    return here - (nz_parts - 1) * (ny_parts - 1);
  case X_POS:
    if (here->world_edge & X_POS_BIT)
      return NULL;
    return here + (nz_parts - 1) * (ny_parts - 1);
  case Y_NEG:
    if (here->world_edge & Y_NEG_BIT)
      return NULL;
    return here - (nz_parts - 1);
  case Y_POS:
    if (here->world_edge & Y_POS_BIT)
      return NULL;
    return here + (nz_parts - 1);
  case Z_NEG:
    if (here->world_edge & Z_NEG_BIT)
      return NULL;
    return here - 1;
  case Z_POS:
    if (here->world_edge & Z_POS_BIT)
      return NULL;
    return here + 1;
  default:
    mcell_internal_error(
        "Invalid direction specified in traverse_subvol (dir=%d).", which);
    return NULL;
  } /* end switch */
}

/*************************************************************************
collide_sv_time:
  In: pointer to a vector3 of where we are (*here)
      pointer to a vector3 of where we want to be
      our current subvolume
  Out: time to hit the closest wall of the subvolume
*************************************************************************/
double collide_sv_time(struct vector3 *here, struct vector3 *move,
                       struct subvolume *sv, double *x_fineparts,
                       double *y_fineparts, double *z_fineparts) {
  double dx, dy, dz, tx, ty, tz, t;

  if ((!distinguishable(move->x, 0, EPS_C)) &&
      (!distinguishable(move->y, 0, EPS_C)) &&
      (!distinguishable(move->z, 0, EPS_C))) {
    return GIGANTIC;
  }

  if (move->x > 0)
    dx = x_fineparts[sv->urb.x] - here->x;
  else {
    dx = x_fineparts[sv->llf.x] - here->x;
  }

  if (move->y > 0)
    dy = y_fineparts[sv->urb.y] - here->y;
  else {
    dy = y_fineparts[sv->llf.y] - here->y;
  }

  if (move->z > 0)
    dz = z_fineparts[sv->urb.z] - here->z;
  else {
    dz = z_fineparts[sv->llf.z] - here->z;
  }

  tx = GIGANTIC;
  if (distinguishable(move->x, 0, EPS_C)) {
    tx = dx / move->x;
  }

  ty = GIGANTIC;
  if (distinguishable(move->y, 0, EPS_C)) {
    ty = dy / move->y;
  }

  tz = GIGANTIC;
  if (distinguishable(move->z, 0, EPS_C)) {
    tz = dz / move->z;
  }

  if (tx < ty) {
    if (tx < tz) {
      t = tx;
    } else {
      t = tz;
    }
  } else {
    if (ty < tz) {
      t = ty;
    } else {
      t = tz;
    }
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
struct subvolume *next_subvol(struct vector3 *here, struct vector3 *move,
                              struct subvolume *sv, double *x_fineparts,
                              double *y_fineparts, double *z_fineparts,
                              int ny_parts, int nz_parts) {
  double dx, dy, dz, tx, ty, tz, t;
  int which;

  int whichx = 1, whichy = 1, whichz = 1;
  if ((!distinguishable(move->x, 0, EPS_C)) &&
      (!distinguishable(move->y, 0, EPS_C)) &&
      (!distinguishable(move->z, 0, EPS_C))) {
    return NULL;
  }

  if (move->x > 0)
    dx = x_fineparts[sv->urb.x] - here->x;
  else {
    dx = x_fineparts[sv->llf.x] - here->x;
    whichx = 0;
  }

  if (move->y > 0)
    dy = y_fineparts[sv->urb.y] - here->y;
  else {
    dy = y_fineparts[sv->llf.y] - here->y;
    whichy = 0;
  }

  if (move->z > 0)
    dz = z_fineparts[sv->urb.z] - here->z;
  else {
    dz = z_fineparts[sv->llf.z] - here->z;
    whichz = 0;
  }

  if (move->x == 0.0) {
    ty = dy * move->z;
    if (ty < 0)
      ty = -ty;
    tz = move->y * dz;
    if (tz < 0)
      tz = -tz;
    if (ty < tz) {
      t = dy / move->y;
      which = Y_NEG + whichy;
    } else {
      t = dz / move->z;
      which = Z_NEG + whichz;
    }
  } else if (move->y == 0.0) {
    tx = dx * move->z;
    if (tx < 0)
      tx = -tx;
    tz = move->x * dz;
    if (tz < 0)
      tz = -tz;
    if (tx < tz) {
      t = dx / move->x;
      which = X_NEG + whichx;
    } else {
      t = dz / move->z;
      which = Z_NEG + whichz;
    }
  } else if (move->z == 0.0) {
    tx = dx * move->y;
    if (tx < 0)
      tx = -tx;
    ty = move->x * dy;
    if (ty < 0)
      ty = -ty;
    if (tx < ty) {
      t = dx / move->x;
      which = X_NEG + whichx;
    } else {
      t = dy / move->y;
      which = Y_NEG + whichy;
    }
  } else {
    tx = dx * move->y * move->z;
    if (tx < 0)
      tx = -tx;
    ty = move->x * dy * move->z;
    if (ty < 0)
      ty = -ty;
    tz = move->x * move->y * dz;
    if (tz < 0)
      tz = -tz;

    if (tx < ty) {
      if (tx < tz) {
        t = dx / move->x;
        which = X_NEG + whichx;
      } else {
        t = dz / move->z;
        which = Z_NEG + whichz;
      }
    } else /* ty<tx */
    {
      if (ty < tz) {
        t = dy / move->y;
        which = Y_NEG + whichy;
      } else {
        t = dz / move->z;
        which = Z_NEG + whichz;
      }
    }
  }

  if (t >= 1.0) {
    here->x += move->x;
    here->y += move->y;
    here->z += move->z;

    return NULL;
  } else {
    here->x += t * move->x;
    here->y += t * move->y;
    here->z += t * move->z;

    t = 1.0 - t;

    move->x *= t;
    move->y *= t;
    move->z *= t;

    return traverse_subvol(sv, which, ny_parts, nz_parts);
  }
}

/*************************************************************************
find_subvolume:
  In: pointer to a vector3 of where we are
      pointer to a subvolume we might be in or near
  Out: subvolume that we are in
*************************************************************************/
struct subvolume *find_subvolume(struct volume *state, struct vector3 *loc,
                                 struct subvolume *guess) {
  /* This code is faster if coarse subvolumes are always used */

  if (guess == NULL)
    return find_coarse_subvol(state, loc);
  else {
    if (state->x_fineparts[guess->llf.x] <= loc->x &&
        loc->x <= state->x_fineparts[guess->urb.x] &&
        state->y_fineparts[guess->llf.y] <= loc->y &&
        loc->y <= state->y_fineparts[guess->urb.y] &&
        state->z_fineparts[guess->llf.z] <= loc->z &&
        loc->z <= state->z_fineparts[guess->urb.z]) {
      return guess;
    } else
      return find_coarse_subvol(state, loc);
  }
}

/*************************************************************************
is_defunct_molecule
  In: abstract_element that is assumed to be an abstract_molecule
  Out: 0 if the properties field is set, 1 if it is NULL
  Note: This function is passed to sched_util so it can tell which
        molecules are active and which are defunct and can be cleaned up.
*************************************************************************/
int is_defunct_molecule(struct abstract_element *e) {
  return ((struct abstract_molecule *)e)->properties == NULL;
}

/*struct surface_molecule **/
/*place_surface_molecule(struct volume *state, struct species *s,*/
/*                       struct vector3 *loc, short orient, double search_diam,*/
/*                       double t, struct subvolume **psv, char *mesh_name,*/
/*                       struct string_buffer *reg_names,*/
/*                       struct string_buffer *regions_to_ignore) {*/
struct wall* find_closest_wall(
    struct volume *state, struct vector3 *loc, double search_diam,
    struct vector2 *best_uv, int *grid_index, struct species *s, char *mesh_name,
    struct string_buffer *reg_names, struct string_buffer *regions_to_ignore) {

  double d2;
  struct vector2 s_loc;

  /*struct vector2 best_uv;*/
  struct vector3 best_xyz;

  double search_d2;
  if (search_diam <= EPS_C)
    search_d2 = EPS_C * EPS_C;
  else
    search_d2 = search_diam * search_diam;

  struct subvolume *sv = find_subvolume(state, loc, NULL);

  char *species_name = s->sym->name;
  unsigned int keyhash = (unsigned int)(intptr_t)(species_name);
  void *key = (void *)(species_name);
  struct mesh_transparency *mesh_transp = (
      struct mesh_transparency *)pointer_hash_lookup(state->species_mesh_transp,
                                                     key, keyhash);

  double best_d2 = search_d2 * 2 + 1;
  struct wall *best_w = NULL;
  struct wall_list *wl;
  for (wl = sv->wall_head; wl != NULL; wl = wl->next) {
    if (verify_wall_regions_match(
        mesh_name, reg_names, wl->this_wall, regions_to_ignore, mesh_transp,
        species_name)) {
      continue; 
    }

    d2 = closest_interior_point(loc, wl->this_wall, &s_loc, search_d2);
    if (d2 <= search_d2 && d2 < best_d2) {
      best_d2 = d2;
      best_w = wl->this_wall;
      best_uv->u = s_loc.u;
      best_uv->v = s_loc.v;
    }
  }

  if (search_d2 > EPS_C * EPS_C) /* Might need to look in adjacent subvolumes */
  {
    const int sv_index = sv - state->subvol;
    int sv_remain = sv_index;

    /* Turn linear sv_index into part_x, part_y, part_z triple. */
    const int part_x =
        sv_remain / ((state->ny_parts - 1) * (state->nz_parts - 1));
    sv_remain -= part_x * ((state->ny_parts - 1) * (state->nz_parts - 1));
    const int part_y = sv_remain / (state->nz_parts - 1);
    sv_remain -= part_y * (state->nz_parts - 1);
    const int part_z = sv_remain;

    /* Find min x partition. */
    int x_min;
    for (x_min = part_x; x_min > 0; x_min--) {
      d2 = loc->x - state->x_partitions[x_min];
      d2 *= d2;
      if (d2 >= best_d2 || d2 >= search_d2)
        break;
    }

    /* Find max x partition. */
    int x_max;
    for (x_max = part_x; x_max < state->nx_parts - 1; x_max++) {
      d2 = loc->x - state->x_partitions[x_max + 1];
      d2 *= d2;
      if (d2 >= best_d2 || d2 >= search_d2)
        break;
    }

    /* Find min y partition. */
    int y_min;
    for (y_min = part_y; y_min > 0; y_min--) {
      d2 = loc->y - state->y_partitions[y_min];
      d2 *= d2;
      if (d2 >= best_d2 || d2 >= search_d2)
        break;
    }

    /* Find max y partition. */
    int y_max;
    for (y_max = part_y; y_max < state->ny_parts - 1; y_max++) {
      d2 = loc->y - state->y_partitions[y_max + 1];
      d2 *= d2;
      if (d2 >= best_d2 || d2 >= search_d2)
        break;
    }

    /* Find min z partition. */
    int z_min;
    for (z_min = part_z; z_min > 0; z_min--) {
      d2 = loc->z - state->z_partitions[z_min];
      d2 *= d2;
      if (d2 >= best_d2 || d2 >= search_d2)
        break;
    }

    /* Find max z partition. */
    int z_max;
    for (z_max = part_z; z_max < state->nz_parts - 1; z_max++) {
      d2 = loc->z - state->z_partitions[z_max + 1];
      d2 *= d2;
      if (d2 >= best_d2 || d2 >= search_d2)
        break;
    }

    if (x_min < part_x || x_max > part_x || y_min < part_y || y_max > part_y ||
        z_min < part_z || z_max > part_z) {
      for (int px = x_min; px <= x_max; px++) {
        for (int py = y_min; py <= y_max; py++) {
          for (int pz = z_min; pz <= z_max; pz++) {
            const int this_sv =
                pz + (state->nz_parts - 1) * (py + (state->ny_parts - 1) * px);
            if (this_sv == sv_index)
              continue;

            for (wl = state->subvol[this_sv].wall_head; wl != NULL;
                 wl = wl->next) {
              if (verify_wall_regions_match(
                  mesh_name, reg_names, wl->this_wall, regions_to_ignore,
                  mesh_transp, species_name)) {
                continue; 
              }

              d2 =
                  closest_interior_point(loc, wl->this_wall, &s_loc, search_d2);
              if (d2 <= search_d2 && d2 < best_d2) {
                best_d2 = d2;
                best_w = wl->this_wall;
                best_uv->u = s_loc.u;
                best_uv->v = s_loc.v;
              }
            }
          }
        }
      }
      if (best_w != NULL) {
        uv2xyz(best_uv, best_w, &best_xyz);
        sv = find_subvolume(state, &best_xyz,
                            sv); /* May have switched subvolumes */
      }
    }
  }

  if (best_w == NULL) {
    return NULL;
  }

  /* We can look this far around the surface we hit for an empty spot */
  d2 = search_d2 - best_d2; 

  if (best_w->grid == NULL) {
    if (create_grid(state, best_w, sv))
      mcell_allocfailed("Failed to create grid for wall.");
    *grid_index = uv2grid(best_uv, best_w->grid);
  } else {
    *grid_index = uv2grid(best_uv, best_w->grid);
    struct surface_molecule_list *sm_list = best_w->grid->sm_list[*grid_index];
    if (sm_list && sm_list->sm) {
      // XXX: this isn't good enough. we should only return this if the PB of
      // sm isn't represented in the PB list.
      if (state->periodic_box_obj && !state->periodic_traditional) {
        return best_w;
      }
      if (d2 <= EPS_C * EPS_C) {
        return NULL;
      } else {
        best_w = search_nbhd_for_free(
            state, best_w, best_uv, d2, grid_index, NULL, NULL, mesh_name,
            reg_names);
        if (best_w == NULL) {
          return NULL;
        }

        if (state->randomize_smol_pos)
          grid2uv_random(best_w->grid, *grid_index, best_uv, state->rng);
        else
          grid2uv(best_w->grid, *grid_index, best_uv);
      }
    }
  }

  return best_w;
}

/*************************************************************************
place_surface_molecule
  In: species for the new molecule
      3D location of the new molecule
      orientation of the new molecule
      diameter to search for a free surface spot
      schedule time for the new molecule
  Out: pointer to the new molecule, or NULL if no free spot was found.
  Note: This function halts the program if it runs out of memory.
        This function is similar to insert_surface_molecule, but it does
        not schedule the molecule or add it to the count.  This is done
        to simplify the logic when placing a surface macromolecule.
        (i.e. place all molecules, and once we're sure we've succeeded,
        schedule them all and count them all.)
 *************************************************************************/
struct surface_molecule *
place_surface_molecule(struct volume *state, struct species *s,
                       struct vector3 *loc, short orient, double search_diam,
                       double t, struct subvolume **psv, char *mesh_name,
                       struct string_buffer *reg_names,
                       struct string_buffer *regions_to_ignore,
                       struct periodic_image *periodic_box) {

  struct vector2 best_uv;
  struct vector3 best_xyz;
  int grid_index = 0;
  int *grid_index_p = &grid_index;
  struct wall *best_w = find_closest_wall(
    state, loc, search_diam, &best_uv, grid_index_p, s, mesh_name, reg_names,
    regions_to_ignore);

  if (best_w == NULL) {
    return NULL; 
  }
  if (state->periodic_box_obj) {
    struct polygon_object *p = (struct polygon_object*)(state->periodic_box_obj->contents);
    struct subdivided_box *sb = p->sb;
    struct vector3 llf = {sb->x[0], sb->y[0], sb->z[0]};
    struct vector3 urb = {sb->x[1], sb->y[1], sb->z[1]};
    struct vector3 pos3d;
    uv2xyz(&best_uv, best_w, &pos3d);
    if (!point_in_box(&llf, &urb, &pos3d)) {
      return NULL;
    }
  }
  struct surface_molecule_list *sm_list = best_w->grid->sm_list[grid_index];
  if (state->periodic_box_obj && periodicbox_in_surfmol_list(periodic_box, sm_list)) {
    return NULL;
  }

  uv2xyz(&best_uv, best_w, &best_xyz);
  struct subvolume *sv = NULL;
  sv = find_subvolume(state, &best_xyz, sv);

  struct surface_molecule *sm;
  sm = CHECKED_MEM_GET(sv->local_storage->smol, "surface molecule");
  sm->mesh_name = NULL;
  sm->birthplace = sv->local_storage->smol;
  sm->birthday = convert_iterations_to_seconds(
      state->start_iterations, state->time_unit,
      state->simulation_start_seconds, t);
  sm->id = state->current_mol_id++;
  sm->properties = s;
  s->population++;
  sm->periodic_box = CHECKED_MALLOC_STRUCT(struct periodic_image,
    "periodic image descriptor");
  sm->periodic_box->x = periodic_box->x;
  sm->periodic_box->y = periodic_box->y;
  sm->periodic_box->z = periodic_box->z;

  sm->flags = TYPE_SURF | ACT_NEWBIE | IN_SCHEDULE;
  if (s->space_step > 0)
    sm->flags |= ACT_DIFFUSE;
  if (trigger_unimolecular(state->reaction_hash, state->rx_hashsize, s->hashval,
                           (struct abstract_molecule *)sm) != NULL ||
      (s->flags & CAN_SURFWALL) != 0)
    sm->flags |= ACT_REACT;

  sm->t = t;
  sm->t2 = 0.0;
  sm->grid = best_w->grid;
  sm->grid_index = grid_index;
  sm->s_pos.u = best_uv.u;
  sm->s_pos.v = best_uv.v;
  sm->orient = orient;

  sm_list = add_surfmol_with_unique_pb_to_list(sm_list, sm);
  if (sm_list == NULL) {
    return NULL; 
  }
  sm->grid->sm_list[sm->grid_index] = sm_list;
  
  sm->grid->n_occupied++;
  sm->flags |= IN_SURFACE;

  if ((s->flags & COUNT_ENCLOSED) != 0)
    sm->flags |= COUNT_ME;

  *psv = sv;
  return sm;
}

/*************************************************************************
insert_surface_molecule
  In: species for the new molecule
      3D location of the new molecule
      orientation of the new molecule
      diameter to search for a free surface spot (vector3 now, should be
double!)
      schedule time for the new molecule
  Out: pointer to the new molecule, or NULL if no free spot was found.
  Note: This function halts the program if it runs out of memory.
*************************************************************************/
struct surface_molecule *
insert_surface_molecule(struct volume *state, struct species *s,
                        struct vector3 *loc, short orient, double search_diam,
                        double t, char *mesh_name,
                        struct string_buffer *reg_names,
                        struct string_buffer *regions_to_ignore,
                        struct periodic_image *periodic_box) {
  struct subvolume *sv = NULL;
  struct surface_molecule *sm =
      place_surface_molecule(
          state, s, loc, orient, search_diam, t, &sv, mesh_name, reg_names,
          regions_to_ignore, periodic_box);
  if (sm == NULL)
    return NULL;

  if (periodic_box != NULL) {
    sm->periodic_box->x = periodic_box->x;
    sm->periodic_box->y = periodic_box->y;
    sm->periodic_box->z = periodic_box->z;
  }

  if (sm->properties->flags & (COUNT_CONTENTS | COUNT_ENCLOSED))
    count_region_from_scratch(state, (struct abstract_molecule *)sm, NULL, 1,
                              NULL, sm->grid->surface, sm->t, NULL);

  if (schedule_add(sv->local_storage->timer, sm))
    mcell_allocfailed("Failed to add surface molecule to scheduler.");

  return sm;
}

/*************************************************************************
insert_volume_molecule
  In: pointer to a volume_molecule that we're going to place in local storage
      pointer to a volume_molecule that may be nearby
  Out: pointer to the new volume_molecule (copies data from volume molecule
       passed in), or NULL if out of memory.  Molecule is placed in scheduler
       also.
*************************************************************************/
struct volume_molecule *insert_volume_molecule(
    struct volume *state, struct volume_molecule *vm,
    struct volume_molecule *vm_guess) {

  struct subvolume *sv;

  if (vm_guess == NULL)
    sv = find_subvolume(state, &(vm->pos), NULL);
  else if (inside_subvolume(&(vm->pos), vm_guess->subvol, state->x_fineparts,
                            state->y_fineparts, state->z_fineparts))
    sv = vm_guess->subvol;
  else
    sv = find_subvolume(state, &(vm->pos), vm_guess->subvol);

  // Make sure this molecule isn't outside of the periodic boundaries
  struct vector3 llf, urb;
  if (state->periodic_box_obj) {
    struct polygon_object *p = (struct polygon_object*)(state->periodic_box_obj->contents);
    struct subdivided_box *sb = p->sb;
    llf = (struct vector3) {sb->x[0], sb->y[0], sb->z[0]};
    urb = (struct vector3) {sb->x[1], sb->y[1], sb->z[1]};
  }
  if (state->periodic_box_obj && !point_in_box(&llf, &urb, &vm->pos)) {
    mcell_error("cannot release '%s' outside of periodic boundaries.",
              vm->properties->sym->name);
    return NULL;
  }

  struct volume_molecule *new_vm;
  new_vm = CHECKED_MEM_GET(sv->local_storage->mol, "volume molecule");
  memcpy(new_vm, vm, sizeof(struct volume_molecule));
  new_vm->mesh_name = NULL;
  new_vm->birthplace = sv->local_storage->mol;
  new_vm->id = state->current_mol_id++;
  new_vm->prev_v = NULL;
  new_vm->next_v = NULL;
  new_vm->next = NULL;
  new_vm->subvol = sv;
  ht_add_molecule_to_list(&sv->mol_by_species, new_vm);
  sv->mol_count++;
  new_vm->properties->population++;
  new_vm->periodic_box = CHECKED_MALLOC_STRUCT(struct periodic_image,
    "periodic image descriptor");
  new_vm->periodic_box->x = vm->periodic_box->x;
  new_vm->periodic_box->y = vm->periodic_box->y;
  new_vm->periodic_box->z = vm->periodic_box->z;

  if ((new_vm->properties->flags & COUNT_SOME_MASK) != 0)
    new_vm->flags |= COUNT_ME;
  if (new_vm->properties->flags & (COUNT_CONTENTS | COUNT_ENCLOSED)) {
    count_region_from_scratch(state, (struct abstract_molecule *)new_vm, NULL,
                              1, &(new_vm->pos), NULL, new_vm->t,
                              new_vm->periodic_box);
  }

  if (schedule_add(sv->local_storage->timer, new_vm))
    mcell_allocfailed("Failed to add volume molecule to scheduler.");
  return new_vm;
}

static int remove_from_list(struct volume_molecule *it) {
  if (it->prev_v) {
#ifdef DEBUG_LIST_CHECKS
    if (*it->prev_v != it) {
      mcell_error_nodie("Stale previous pointer!");
    }
#endif
    *(it->prev_v) = it->next_v;
  } else {
#ifdef DEBUG_LIST_CHECKS
    mcell_error_nodie("No previous pointer.");
#endif
  }
  if (it->next_v) {
#ifdef DEBUG_LIST_CHECKS
    if (it->next_v->prev_v != &it->next_v) {
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
struct volume_molecule *migrate_volume_molecule(struct volume_molecule *vm,
                                                struct subvolume *new_sv) {
  struct volume_molecule *new_vm;

  new_sv->mol_count++;
  vm->subvol->mol_count--;

  if (vm->subvol->local_storage == new_sv->local_storage) {
    if (remove_from_list(vm)) {
      vm->subvol = new_sv;
      ht_add_molecule_to_list(&new_sv->mol_by_species, vm);
      return vm;
    }
  }

  new_vm = CHECKED_MEM_GET(new_sv->local_storage->mol, "volume molecule");
  memcpy(new_vm, vm, sizeof(struct volume_molecule));
  new_vm->birthplace = new_sv->local_storage->mol;
  new_vm->mesh_name = NULL;
  new_vm->prev_v = NULL;
  new_vm->next_v = NULL;
  new_vm->next = NULL;
  new_vm->subvol = new_sv;

  ht_add_molecule_to_list(&new_sv->mol_by_species, new_vm);

  collect_molecule(vm);

  return new_vm;
}

/*************************************************************************
eval_rel_region_3d:
  In: an expression tree containing regions to release on
      the waypoint for the current subvolume
      a list of regions entered from the waypoint to the release loc.
      a list of regions exited from the waypoint to the release loc.
  Out: 1 if the location chosen satisfies the expression, 0 if not.
*************************************************************************/
int eval_rel_region_3d(struct release_evaluator *expr, struct waypoint *wp,
                       struct region_list *in_regions,
                       struct region_list *out_regions) {
  struct region *r;
  struct region_list *rl;
  int satisfies_l, satisfies_r;

  satisfies_l = 0;
  if (expr->op & REXP_LEFT_REGION) {
    r = (struct region *)expr->left;
    for (rl = wp->regions; rl != NULL; rl = rl->next) {
      if (rl->reg == r) {
        satisfies_l = 1;
        break;
      }
    }
    if (satisfies_l) {
      for (rl = out_regions; rl != NULL; rl = rl->next) {
        if (rl->reg == r) {
          satisfies_l = 0;
          break;
        }
      }
    } else {
      for (rl = in_regions; rl != NULL; rl = rl->next) {
        if (rl->reg == r) {
          satisfies_l = 1;
          break;
        }
      }
    }
  } else
    satisfies_l = eval_rel_region_3d(expr->left, wp, in_regions, out_regions);

  if (expr->op & REXP_NO_OP)
    return satisfies_l;

  satisfies_r = 0;
  if (expr->op & REXP_RIGHT_REGION) {
    r = (struct region *)expr->right;
    for (rl = wp->regions; rl != NULL; rl = rl->next) {
      if (rl->reg == r) {
        satisfies_r = 1;
        break;
      }
    }
    if (satisfies_r) {
      for (rl = out_regions; rl != NULL; rl = rl->next) {
        if (rl->reg == r) {
          satisfies_r = 0;
          break;
        }
      }
    } else {
      for (rl = in_regions; rl != NULL; rl = rl->next) {
        if (rl->reg == r) {
          satisfies_r = 1;
          break;
        }
      }
    }
  } else
    satisfies_r = eval_rel_region_3d(expr->right, wp, in_regions, out_regions);

  if (expr->op & REXP_UNION)
    return (satisfies_l || satisfies_r);
  else if (expr->op & (REXP_INTERSECTION))
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
       state as specified.
  Note: if more molecules are to be removed than actually exist, all
        existing molecules of the specified type are removed.
*************************************************************************/
static int vacuum_inside_regions(struct volume *state,
                                 struct release_site_obj *rso,
                                 struct volume_molecule *vm, int n) {
  struct volume_molecule *mp;
  struct release_region_data *rrd;
  struct region_list *extra_in, *extra_out;
  struct region_list *rl, *rl2;
  struct waypoint *wp;
  struct subvolume *sv = NULL;
  struct mem_helper *mh;
  struct void_list *vl;
  struct void_list *vl_head = NULL;
  int vl_num = 0;
  double t;
  struct vector3 hit, delta;
  struct vector3 *origin;
  struct wall_list *wl;

  rrd = rso->region_data;
  mh = create_mem(sizeof(struct void_list), 1024);
  if (mh == NULL)
    return 1;

  const int x_min = bisect(state->x_partitions, state->nx_parts, rrd->llf.x);
  const int x_max =
      bisect_high(state->x_partitions, state->nx_parts, rrd->urb.x);
  const int y_min = bisect(state->y_partitions, state->ny_parts, rrd->llf.y);
  const int y_max =
      bisect_high(state->y_partitions, state->ny_parts, rrd->urb.y);
  const int z_min = bisect(state->z_partitions, state->nz_parts, rrd->llf.z);
  const int z_max =
      bisect_high(state->z_partitions, state->nz_parts, rrd->urb.z);

  for (int px = x_min; px < x_max; px++) {
    for (int py = y_min; py < y_max; py++) {
      for (int pz = z_min; pz < z_max; pz++) {
        const int this_sv =
            pz + (state->nz_parts - 1) * (py + (state->ny_parts - 1) * px);
        sv = &(state->subvol[this_sv]);

        struct per_species_list *psl =
            (struct per_species_list *)pointer_hash_lookup(
                &sv->mol_by_species, vm->properties, vm->properties->hashval);

        if (psl != NULL) {
          for (mp = psl->head; mp != NULL; mp = mp->next_v) {
            extra_in = extra_out = NULL;
            wp = &(state->waypoints[this_sv]);
            origin = &(wp->loc);
            delta.x = mp->pos.x - origin->x;
            delta.y = mp->pos.y - origin->y;
            delta.z = mp->pos.z - origin->z;

            for (wl = sv->wall_head; wl != NULL; wl = wl->next) {
              int hitcode = collide_wall(origin, &delta, wl->this_wall, &t,
                                         &hit, 0, state->rng, state->notify,
                                         &(state->ray_polygon_tests));
              if (hitcode != COLLIDE_MISS) {
                state->ray_polygon_colls++;

                for (rl = wl->this_wall->counting_regions; rl != NULL;
                     rl = rl->next) {
                  if (hitcode == COLLIDE_FRONT || hitcode == COLLIDE_BACK) {
                    rl2 = (struct region_list *)CHECKED_MEM_GET(
                        sv->local_storage->regl, "region list");
                    rl2->reg = rl->reg;

                    if (hitcode == COLLIDE_FRONT) {
                      rl2->next = extra_in;
                      extra_in = rl2;
                    } else /*hitcode == COLLIDE_BACK*/
                    {
                      rl2->next = extra_out;
                      extra_out = rl2;
                    }
                  }
                }
              }
            }

            for (rl = extra_in; rl != NULL; rl = rl->next) {
              if (rl->reg == NULL)
                continue;
              for (rl2 = extra_out; rl2 != NULL; rl2 = rl2->next) {
                if (rl2->reg == NULL)
                  continue;
                if (rl->reg == rl2->reg) {
                  rl->reg = NULL;
                  rl2->reg = NULL;
                  break;
                }
              }
            }

            if (eval_rel_region_3d(rrd->expression, wp, extra_in, extra_out)) {
              vl = (struct void_list *)CHECKED_MEM_GET(mh, "temporary list");
              vl->data = mp;
              vl->next = vl_head;
              vl_head = vl;
              vl_num++;
            }

            if (extra_in != NULL)
              mem_put_list(sv->local_storage->regl, extra_in);
            if (extra_out != NULL)
              mem_put_list(sv->local_storage->regl, extra_out);
          }
        }
      }
    }
  }

  for (vl = vl_head; n < 0 && vl_num > 0 && vl != NULL;
       vl = vl->next, vl_num--) {
    if (rng_dbl(state->rng) < ((double)(-n)) / ((double)vl_num)) {
      mp = (struct volume_molecule *)vl->data;
      mp->properties->population--;
      mp->subvol->mol_count--;
      if ((mp->properties->flags & (COUNT_CONTENTS | COUNT_ENCLOSED)) != 0)
        count_region_from_scratch(state, (struct abstract_molecule *)mp, NULL,
                                  -1, &(mp->pos), NULL, mp->t, NULL);
      if (mp->flags & IN_SCHEDULE) {
        mp->subvol->local_storage->timer
            ->defunct_count++; /* Tally for garbage collection */
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
static int is_point_inside_region(struct volume *state,
                                  struct vector3 const *pos,
                                  struct release_evaluator *expression,
                                  struct subvolume *sv) {
  struct region_list *extra_in = NULL, *extra_out = NULL, *cur_region;
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
      !inside_subvolume((struct vector3 *)pos, sv, state->x_fineparts,
                        state->y_fineparts, state->z_fineparts))
    sv = find_subvolume(state, (struct vector3 *)pos, sv);

  /* Find waypoint, compute trajectory from waypoint */
  wp = &(state->waypoints[sv - state->subvol]);
  origin = &(wp->loc);
  delta.x = pos->x - origin->x;
  delta.y = pos->y - origin->y;
  delta.z = pos->z - origin->z;

  for (wl = sv->wall_head; wl != NULL; wl = wl->next) {
    struct vector3 hit_pos;
    double hit_time;
    int hit_check =
        collide_wall(origin, &delta, wl->this_wall, &hit_time, &hit_pos, 0,
                     state->rng, state->notify, &(state->ray_polygon_tests));

    if (hit_check != COLLIDE_MISS) {
      state->ray_polygon_colls++;

      if ((hit_time > -EPS_C && hit_time < EPS_C) ||
          (hit_time > 1.0 - EPS_C && hit_time < 1.0 + EPS_C)) {
        bad_location = 1;
        break;
      }

      for (cur_region = wl->this_wall->counting_regions; cur_region != NULL;
           cur_region = cur_region->next) {
        struct region_list *crossed_region =
            (struct region_list *)CHECKED_MEM_GET(sv->local_storage->regl,
                                                  "region list");
        crossed_region->reg = cur_region->reg;

        if (hit_check == COLLIDE_FRONT) {
          crossed_region->next = extra_in;
          extra_in = crossed_region;
        } else if (hit_check == COLLIDE_BACK) {
          crossed_region->next = extra_out;
          extra_out = crossed_region;
        } else {
          bad_location = 1;
          break;
        }
      }
    }
  }

  if (bad_location) {
    if (extra_in != NULL)
      mem_put_list(sv->local_storage->regl, extra_in);
    if (extra_out != NULL)
      mem_put_list(sv->local_storage->regl, extra_out);
    return 0;
  }

  for (cur_region = extra_in; cur_region != NULL;
       cur_region = cur_region->next) {
    struct region_list *out_region = NULL;
    if (cur_region->reg == NULL)
      continue;
    for (out_region = extra_out; out_region != NULL;
         out_region = out_region->next) {
      if (out_region->reg == NULL)
        continue;
      if (cur_region->reg == out_region->reg) {
        cur_region->reg = NULL;
        out_region->reg = NULL;
        break;
      }
    }
  }

  result = eval_rel_region_3d(expression, wp, extra_in, extra_out);

  if (extra_in != NULL)
    mem_put_list(sv->local_storage->regl, extra_in);
  if (extra_out != NULL)
    mem_put_list(sv->local_storage->regl, extra_out);
  return result;
}

/*************************************************************************
release_inside_regions:
  In: pointer to a release site object
      template molecule to release
      integer number of molecules to release
  Out: 0 on success, 1 on failure; molecule(s) are released into the
       state as specified.
  Note: if the CCNNUM release method is used, the number of molecules
        passed in is ignored.
*************************************************************************/
static int release_inside_regions(struct volume *state,
                                  struct release_site_obj *rso,
                                  struct volume_molecule *vm, int n) {

  struct release_region_data *rrd = rso->region_data;
  vm->previous_wall = NULL;
  vm->index = -1;

  // test if the release region is a single object (versus a CSG expression)
  // since then we can compute the number of molecules exactly.
  bool exactNumber = false;
  if (rso->release_number_method == CCNNUM) {
    n = num_vol_mols_from_conc(rso, state->length_unit, &exactNumber);
  }

  if (n < 0)
    return vacuum_inside_regions(state, rso, vm, n);

  struct volume_molecule *new_vm = NULL;
  struct subvolume *sv = NULL;
  while (n > 0) {
    vm->pos.x = rrd->llf.x + (rrd->urb.x - rrd->llf.x) * rng_dbl(state->rng);
    vm->pos.y = rrd->llf.y + (rrd->urb.y - rrd->llf.y) * rng_dbl(state->rng);
    vm->pos.z = rrd->llf.z + (rrd->urb.z - rrd->llf.z) * rng_dbl(state->rng);

    if (!is_point_inside_region(state, &vm->pos, rrd->expression, NULL)) {
      if (rso->release_number_method == CCNNUM && !exactNumber)
        n--;
      continue;
    }

    /* Actually place the molecule */
    vm->subvol = sv;
    vm->periodic_box->x = rso->periodic_box->x;
    vm->periodic_box->y = rso->periodic_box->y;
    vm->periodic_box->z = rso->periodic_box->z;
    new_vm = insert_volume_molecule(state, vm, new_vm);
    if (new_vm == NULL)
      return 1;

    n--;
  }

  return 0;
}

/*************************************************************************
release_molecules:
  In: pointer to a release event
  Out: 0 on success, 1 on failure; next event is scheduled and molecule(s)
       are released into the state as specified.
  Note: if a release is triggered by a reaction, there isn't anything
        to schedule.  Also, in that case, rpat isn't really a release
        pattern (it's a rxn_pathname in disguise) so be sure to not
        dereference it!
*************************************************************************/
int release_molecules(struct volume *state, struct release_event_queue *req) {

  if (req == NULL)
    return 0;

  struct release_site_obj *rso = req->release_site;

  struct volume_molecule vm;
  memset(&vm, 0, sizeof(struct volume_molecule));

  /* Set up canonical molecule to be released */
  /* If we have a list, assume a 3D molecule and fix later */
  if (rso->mol_list != NULL || (rso->mol_type->flags & NOT_FREE) == 0) {
    vm.flags = TYPE_VOL | IN_VOLUME;
  } else {
    vm.flags = TYPE_SURF | IN_SURFACE;
  }
  vm.flags |= IN_SCHEDULE | ACT_NEWBIE;

  if (req->train_counter == 0) {
    req->train_counter++;
  }

  struct release_pattern *rpat = rso->pattern;
  if (skip_past_events(rso->release_prob, state, req, rpat)) {
    return 0;
  }

  if (check_release_probability(rso->release_prob, state, req, rpat)) {
    return 0;
  }

  // Set molecule characteristics.
  vm.mesh_name = NULL;
  vm.t = req->event_time;
  vm.properties = rso->mol_type;
  vm.t2 = 0.0;
  vm.birthday = convert_iterations_to_seconds(
      state->start_iterations, state->time_unit,
      state->simulation_start_seconds, vm.t);
  struct periodic_image periodic_box = { .x = rso->periodic_box->x,
                                         .y = rso->periodic_box->y,
                                         .z = rso->periodic_box->z
                                       };
  vm.periodic_box = &periodic_box;

  struct abstract_molecule *ap = (struct abstract_molecule *)(&vm);

  // All molecules are the same, so we can set flags
  if (rso->mol_list == NULL) {
    if (trigger_unimolecular(state->reaction_hash, state->rx_hashsize,
                             rso->mol_type->hashval, ap) != NULL ||
        (rso->mol_type->flags & CAN_SURFWALL) != 0)
      ap->flags |= ACT_REACT;
    if (rso->mol_type->space_step > 0.0)
      ap->flags |= ACT_DIFFUSE;
  }

  int number = calculate_number_to_release(rso, state);

  if (rso->release_shape == SHAPE_REGION) {
    u_int pop_before = ap->properties->population;
    if (ap->flags & TYPE_VOL) {
      if (release_inside_regions(state, rso, (struct volume_molecule *)ap,
                                 number))
        return 1;
    } else {
      if (release_onto_regions(state, rso, (struct surface_molecule *)ap,
                               number))
        return 1;
    }
    if (state->notify->release_events == NOTIFY_FULL) {
      if (number >= 0) {
        mcell_log("Released %d %s from \"%s\" at iteration %lld.",
                  ap->properties->population - pop_before,
                  rso->mol_type->sym->name, rso->name, state->current_iterations);
      } else {
        mcell_log("Removed %d %s from \"%s\" at iteration %lld.",
                  pop_before - ap->properties->population,
                  rso->mol_type->sym->name, rso->name, state->current_iterations);
      }
    }
  }
  // Guaranteed to be 3D molec or at least specified by 3D location if in list
  else {
    vm.previous_wall = NULL;
    vm.index = -1;

    if (rso->mol_list != NULL) {
      if (release_by_list(state, req, &vm)) {
        return 1;
      }
    } else if (rso->diameter != NULL) {

      if (release_ellipsoid_or_rectcuboid(state, req, &vm, number)) {
        return 1;
      }
    } else {
      double location[1][4];
      location[0][0] = rso->location->x;
      location[0][1] = rso->location->y;
      location[0][2] = rso->location->z;
      location[0][3] = 1;

      mult_matrix(location, req->t_matrix, location, 1, 4, 4);

      vm.pos.x = location[0][0];
      vm.pos.y = location[0][1];
      vm.pos.z = location[0][2];

      struct volume_molecule *vm_guess = NULL;
      for (int i = 0; i < number; i++) {
        vm_guess = insert_volume_molecule(state, &vm, vm_guess);
        if (vm_guess == NULL)
          return 1;
        vm.periodic_box->x = rso->periodic_box->x;
        vm.periodic_box->y = rso->periodic_box->y;
        vm.periodic_box->z = rso->periodic_box->z;
      }
      if (state->notify->release_events == NOTIFY_FULL) {
        mcell_log("Released %d %s from \"%s\" at iteration %lld.", number,
                  rso->mol_type->sym->name, rso->name, state->current_iterations);
      }
    }
  }

  /* Schedule next release event. */
  if (!distinguishable(rso->release_prob, MAGIC_PATTERN_PROBABILITY, EPS_C))
    return 0; /* Triggered by reaction, don't schedule */
  req->event_time += rpat->release_interval;

  /* we may need to move to the next train. */
  if (!distinguishable(req->event_time,
                       req->train_high_time + rpat->train_duration, EPS_C) ||
      req->event_time > req->train_high_time + rpat->train_duration) {
    req->train_high_time += rpat->train_interval;
    req->event_time = req->train_high_time;
    req->train_counter++;
  }

  if (req->train_counter <= rpat->number_of_trains &&
      req->event_time < FOREVER) {
    if (schedule_add(state->releaser, req))
      mcell_allocfailed("Failed to add release request to scheduler.");
  }
  else {
    free(req);
  }

  return 0;
}

/*************************************************************************
 release_ellipsoid_or_rectcuboid:
    This function is used for CUBIC (aka RECTANGULAR), SPHERICAL (aka
    ELLIPTICAL), and SPHERICAL_SHELL release sites.

  In: state: MCell simulation state
      req:
      vm: volume molecule being released
      number: number to release
  Out: 0 on success, 1 on failure
*************************************************************************/
int release_ellipsoid_or_rectcuboid(struct volume *state,
                                    struct release_event_queue *req,
                                    struct volume_molecule *vm, int number) {

  struct release_site_obj *rso = req->release_site;
  struct vector3 *diam_xyz = rso->diameter;

  struct vector3 pos;
  const int is_spheroidal = (rso->release_shape == SHAPE_SPHERICAL ||
                             rso->release_shape == SHAPE_ELLIPTIC ||
                             rso->release_shape == SHAPE_SPHERICAL_SHELL);

  for (int i = 0; i < number; i++) {
    do /* Pick values in unit square, toss if not in unit circle */
    {
      pos.x = (rng_dbl(state->rng) - 0.5);
      pos.y = (rng_dbl(state->rng) - 0.5);
      pos.z = (rng_dbl(state->rng) - 0.5);
    } while (is_spheroidal &&
             pos.x * pos.x + pos.y * pos.y + pos.z * pos.z >= 0.25);

    if (rso->release_shape == SHAPE_SPHERICAL_SHELL) {
      double r = sqrt(pos.x * pos.x + pos.y * pos.y + pos.z * pos.z) * 2.0;
      if (r == 0.0) {
        pos.x = 0.0;
        pos.y = 0.0;
        pos.z = 0.5;
      } else {
        pos.x /= r;
        pos.y /= r;
        pos.z /= r;
      }
    }

    double location[1][4];
    location[0][0] = pos.x * diam_xyz->x + rso->location->x;
    location[0][1] = pos.y * diam_xyz->y + rso->location->y;
    location[0][2] = pos.z * diam_xyz->z + rso->location->z;
    location[0][3] = 1;

    mult_matrix(location, req->t_matrix, location, 1, 4, 4);

    vm->pos.x = location[0][0];
    vm->pos.y = location[0][1];
    vm->pos.z = location[0][2];
    struct volume_molecule *guess = NULL;
    /* Insert copy of vm into state */
    vm->periodic_box->x = rso->periodic_box->x;
    vm->periodic_box->y = rso->periodic_box->y;
    vm->periodic_box->z = rso->periodic_box->z;
    guess = insert_volume_molecule(state, vm, guess); 
    if (guess == NULL)
      return 1;
  }
  if (state->notify->release_events == NOTIFY_FULL) {
    mcell_log("Released %d %s from \"%s\" at iteration %lld.", number,
              rso->mol_type->sym->name, rso->name, state->current_iterations);
  }
  return 0;
}

/*************************************************************************
release_by_list:
    This function is used for LIST based release sites.

  In: state: MCell simulation state
      req:
      vm: volume molecule being released
  Out: 0 on success, 1 on failure
*************************************************************************/
int release_by_list(struct volume *state, struct release_event_queue *req,
                    struct volume_molecule *vm) {

  int i = 0;        /* serves as counter for released molecules */
  int i_failed = 0; /* serves as counted for the failed to release molecules */

  struct release_site_obj *rso = req->release_site;
  struct release_single_molecule *rsm = rso->mol_list;

  for (; rsm != NULL; rsm = rsm->next) {
    double location[1][4];
    location[0][0] = rsm->loc.x + rso->location->x;
    location[0][1] = rsm->loc.y + rso->location->y;
    location[0][2] = rsm->loc.z + rso->location->z;
    location[0][3] = 1;

    mult_matrix(location, req->t_matrix, location, 1, 4, 4);

    vm->pos.x = location[0][0];
    vm->pos.y = location[0][1];
    vm->pos.z = location[0][2];

    struct volume_molecule *vm_guess = NULL;
    if ((rsm->mol_type->flags & NOT_FREE) == 0) {
      struct abstract_molecule *ap = (struct abstract_molecule *)(vm);
      vm->properties = rsm->mol_type;
      // Have to set flags, since insert_volume_molecule doesn't
      if (trigger_unimolecular(state->reaction_hash, state->rx_hashsize,
                               ap->properties->hashval, ap) != NULL ||
          (ap->properties->flags & CAN_SURFWALL) != 0) {
        ap->flags |= ACT_REACT;
      }
      if (vm->properties->space_step > 0.0)
        ap->flags |= ACT_DIFFUSE;
      vm_guess = insert_volume_molecule(state, vm, vm_guess);
      i++;
      if (vm_guess == NULL)
        return 1;
      vm_guess->periodic_box->x = rso->periodic_box->x;
      vm_guess->periodic_box->y = rso->periodic_box->y;
      vm_guess->periodic_box->z = rso->periodic_box->z;
      i++;
    } else {
      double diam;
      if (rso->diameter == NULL)
        diam = 0.0;
      else
        diam = rso->diameter->x;

      short orient;
      if (rsm->orient > 0)
        orient = 1;
      else if (rsm->orient < 0)
        orient = -1;
      else {
        orient = (rng_uint(state->rng) & 1) ? 1 : -1;
      }

      // Don't have to set flags, insert_surface_molecule takes care of it
      struct surface_molecule *sm;
      sm = insert_surface_molecule(state, rsm->mol_type, &vm->pos, orient,
                                   diam, req->event_time, NULL, NULL, NULL,
                                   rso->periodic_box);
      if (sm == NULL) {
        mcell_warn("Molecule release is unable to find surface upon which "
                   "to place molecule %s.\n"
                   "  This could be caused by too small of a SITE_DIAMETER "
                   "on the release site '%s'.",
                   rsm->mol_type->sym->name, rso->name);
        i_failed++;
      } else {
        i++;
      }
    }
  }
  if (state->notify->release_events == NOTIFY_FULL) {
    mcell_log("Released %d molecules from list \"%s\" at iteration %lld.", i,
              rso->name, state->current_iterations);
  }
  if (i_failed > 0)
    mcell_warn("Failed to release %d molecules from list \"%s\" at "
               "iteration %lld.",
               i_failed, rso->name, state->current_iterations);

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
static void find_exponential_params(double c, double C, double d, double N,
                                    double *A, double *B, double *k) {

  double k_min = 0;
  double k_mid = 0;
  double k_max = log(GIGANTIC) / N;
  for (int i = 0; i < 720; i++) {
    k_mid = 0.5 * (k_min + k_max);
    double f = c + (exp(N * k_mid) - 1.0) * d / (exp(k_mid) - 1.0);
    if (C > f)
      k_min = k_mid;
    else
      k_max = k_mid;
    if ((k_max - k_min) / (k_max + k_min) < EPS_C)
      break;
  }

  *k = k_mid;
  *A = d / (exp(*k) - 1.0);
  *B = c - *A;
}

/*************************************************************************
 check_partitions_against_interaction_diameter:
  In: nothing.  Uses struct volume *state, assumes partitions are set.
  Out: 0 on success, 1 on error
*************************************************************************/
static int check_partitions_against_interaction_diameter(struct volume *state) {
  int i;

  if (state->x_partitions != NULL) {
    for (i = 1; i < state->nx_parts; i++) {
      if (state->x_partitions[i] - state->x_partitions[i - 1] <
          2 * state->rx_radius_3d) {
        mcell_error("X partitions closer than interaction diameter\n"
                    "  X partition #%d at %g\n"
                    "  X partition #%d at %g\n"
                    "  Interaction diameter %g",
                    i, state->length_unit * state->x_partitions[i - 1], i + 1,
                    state->length_unit * state->x_partitions[i],
                    2 * state->length_unit * state->rx_radius_3d);
        /*return 1;*/
      }
    }
  }
  if (state->y_partitions != NULL) {
    for (i = 1; i < state->ny_parts; i++) {
      if (state->y_partitions[i] - state->y_partitions[i - 1] <
          2 * state->rx_radius_3d) {
        mcell_error("Y partitions closer than interaction diameter\n"
                    "  Y partition #%d at %g\n"
                    "  Y partition #%d at %g\n"
                    "  Interaction diameter %g",
                    i, state->length_unit * state->y_partitions[i - 1], i + 1,
                    state->length_unit * state->y_partitions[i],
                    2 * state->length_unit * state->rx_radius_3d);
        /*return 1;*/
      }
    }
  }
  if (state->z_partitions != NULL) {
    for (i = 1; i < state->nz_parts; i++) {
      if (state->z_partitions[i] - state->z_partitions[i - 1] <
          2 * state->rx_radius_3d) {
        mcell_error("Z partitions closer than interaction diameter\n"
                    "  Z partition #%d at %g\n"
                    "  Z partition #%d at %g\n"
                    "  Interaction diameter %g\n",
                    i, state->length_unit * state->z_partitions[i - 1], i + 1,
                    state->length_unit * state->z_partitions[i],
                    2 * state->length_unit * state->rx_radius_3d);
        /*return 1;*/
      }
    }
  }
  return 0;
}

/*************************************************************************
set_partitions:
  In: nothing.  Uses struct volume *state, assumes bounding box is set.
  Out: 0 on success, 1 on error; coarse and fine partitions are set.
*************************************************************************/
int set_partitions(struct volume *state) {
  /* Set sensible bounds for spacing between fine partitions (minimum size of
   * subdivision) */
  double smallest_spacing = 0.1 * state->r_length_unit; /* 100nm */

  if (2 * state->rx_radius_3d > smallest_spacing)
    smallest_spacing = 2 * state->rx_radius_3d;

  /* We have 2^15 possible fine partitions; we'll use 24k of them */
  if (state->n_fineparts != 4096 + 16384 + 4096) {

    state->n_fineparts = 4096 + 16384 + 4096;
    state->x_fineparts =
        CHECKED_MALLOC_ARRAY(double, state->n_fineparts, "x fine partitions");
    state->y_fineparts =
        CHECKED_MALLOC_ARRAY(double, state->n_fineparts, "y fine partitions");
    state->z_fineparts =
        CHECKED_MALLOC_ARRAY(double, state->n_fineparts, "z fine partitions");
  }

  // Something like the maximum expected error--not sure exactly what this is
  double dfx = 1e-3 + (state->bb_urb.x - state->bb_llf.x) / 8191.0;
  double dfy = 1e-3 + (state->bb_urb.y - state->bb_llf.y) / 8191.0;
  double dfz = 1e-3 + (state->bb_urb.z - state->bb_llf.z) / 8191.0;

  // Make sure fine partition follow the "2 x reaction radius" rule, I guess.
  double f_min = state->bb_llf.x - dfx;
  double f_max = state->bb_urb.x + dfx;
  dfx = increase_fine_partition_size(state, state->x_fineparts, &f_min, &f_max,
                                     smallest_spacing);
  struct vector3 part_min, part_max;
  part_min.x = f_min;
  part_max.x = f_max;

  // Same thing for y as we just did for x
  f_min = state->bb_llf.y - dfy;
  f_max = state->bb_urb.y + dfy;
  dfy = increase_fine_partition_size(state, state->y_fineparts, &f_min, &f_max,
                                     smallest_spacing);
  part_min.y = f_min;
  part_max.y = f_max;

  // And same again for z
  f_min = state->bb_llf.z - dfz;
  f_max = state->bb_urb.z + dfz;
  dfz = increase_fine_partition_size(state, state->z_fineparts, &f_min, &f_max,
                                     smallest_spacing);
  part_min.z = f_min;
  part_max.z = f_max;

  /* Try to figure out how many timesteps our fastest particle can make in the
   * whole state (along longest and shortest axes) */
  double f = part_max.x - part_min.x;
  f_min = f_max = f;
  f = part_max.y - part_min.y;
  if (f < f_min)
    f_min = f;
  else if (f > f_max)
    f_max = f;
  f = part_max.z - part_min.z;
  if (f < f_min)
    f_min = f;
  else if (f > f_max)
    f_max = f;

  double steps_min, steps_max;
  if (!distinguishable(state->speed_limit, 0, EPS_C)) {
    steps_min = f_min;
    steps_max = f_max;
  } else {
    steps_min = f_min / state->speed_limit;
    steps_max = f_max / state->speed_limit;
  }

  /* Verify that partitions are not closer than interaction diameter. */
  if (check_partitions_against_interaction_diameter(state))
    return 1;

  if (state->x_partitions != NULL || state->y_partitions != NULL ||
      state->z_partitions != NULL) {
    if (state->x_partitions == NULL)
      mcell_error("Some axes are partitioned, but the X-axis is not.");
    if (state->y_partitions == NULL)
      mcell_error("Some axes are partitioned, but the Y-axis is not.");
    if (state->z_partitions == NULL)
      mcell_error("Some axes are partitioned, but the Z-axis is not.");
  }

  /* Use automatic partitioning only when there are no user-specified
   * partitions */
  if (state->x_partitions == NULL && state->y_partitions == NULL &&
      state->z_partitions == NULL) {
    set_auto_partitions(state, steps_min, steps_max, &part_min, &part_max,
                        f_max, smallest_spacing);
  } else {
    set_user_partitions(state, dfx, dfy, dfz);
  }

  /* And finally we tell the user what happened */
  if (state->notify->partition_location == NOTIFY_FULL) {
    mcell_log_raw("X partitions: ");
    mcell_log_raw("-inf ");
    for (int i = 1; i < state->nx_parts - 1; i++)
      mcell_log_raw("%.5f ", state->length_unit * state->x_partitions[i]);
    mcell_log_raw("inf\n");
    mcell_log_raw("Y partitions: ");
    mcell_log_raw("-inf ");
    for (int i = 1; i < state->ny_parts - 1; i++)
      mcell_log_raw("%.5f ", state->length_unit * state->y_partitions[i]);
    mcell_log_raw("inf\n");
    mcell_log_raw("Z partitions: ");
    mcell_log_raw("-inf ");
    for (int i = 1; i < state->nz_parts - 1; i++)
      mcell_log_raw("%.5f ", state->length_unit * state->z_partitions[i]);
    mcell_log_raw("inf\n");
  }

  return 0;
}

double increase_fine_partition_size(struct volume *state, double *fineparts,
                                    double *f_min, double *f_max,
                                    double smallest_spacing) {
  /* Not sure how this is supposed to work--looks like two ideas mixed,
   * probably broken */
  /* Was supposed to make sure that the fine partitions would still obey the
   * 2*reaction radius rule */

  if (*f_max - *f_min < smallest_spacing) {
    if (state->notify->progress_report != NOTIFY_NONE)
      mcell_log_raw("Rescaling: was %.3f to %.3f, now ",
                    (*f_min) * state->length_unit,
                    (*f_max) * state->length_unit);
    double f = smallest_spacing - (*f_max - *f_min);
    *f_max += 0.5 * f;
    *f_min -= 0.5 * f;
    if (state->notify->progress_report != NOTIFY_NONE)
      mcell_log_raw("%.3f to %.3f\n", (*f_min) * state->length_unit,
                    (*f_max) * state->length_unit);
  }
  // Set bounds over which to do linear subdivision (state bounding box)
  double df = (*f_max - *f_min) / 16383.0;
  // Subdivide state bounding box
  for (int i = 0; i < 16384; i++) {
    fineparts[4096 + i] = *f_min + df * ((double)i);
  }

  /* Create an exponentially increasing fine partition size as we go to
   * -infinity */
  double A, B, k;
  find_exponential_params(-*f_min, 1e12, df, 4096, &A, &B, &k);
  for (int i = 1; i <= 4096; i++)
    fineparts[4096 - i] = -(A * exp(i * k) + B);
  /* And again as we go to +infinity */
  find_exponential_params(*f_max, 1e12, df, 4096, &A, &B, &k);
  for (int i = 1; i <= 4096; i++)
    fineparts[4096 + 16383 + i] = A * exp(i * k) + B;
  return df;
}

void set_auto_partitions(struct volume *state, double steps_min,
                         double steps_max, struct vector3 *part_min,
                         struct vector3 *part_max, double f_max,
                         double smallest_spacing) {
  /* perform automatic partitioning */

  /* Guess how big to make partitions--nothing really clever about what's done
   * here */
  if (steps_max / MAX_TARGET_TIMESTEP > MAX_COARSE_PER_AXIS) {
    state->nx_parts = state->ny_parts = state->nz_parts = MAX_COARSE_PER_AXIS;
  } else if (steps_min / MIN_TARGET_TIMESTEP < MIN_COARSE_PER_AXIS) {
    state->nx_parts = state->ny_parts = state->nz_parts = MIN_COARSE_PER_AXIS;
  } else {
    state->nx_parts = steps_min / MIN_TARGET_TIMESTEP;
    if (state->nx_parts > MAX_COARSE_PER_AXIS)
      state->nx_parts = MAX_COARSE_PER_AXIS;
    if ((state->nx_parts & 1) != 0)
      state->nx_parts += 1;

    state->ny_parts = state->nz_parts = state->nx_parts;
  }

  /* Allocate memory for our automatically created partitions */
  state->x_partitions =
      CHECKED_MALLOC_ARRAY(double, state->nx_parts, "x partitions");
  state->y_partitions =
      CHECKED_MALLOC_ARRAY(double, state->ny_parts, "y partitions");
  state->z_partitions =
      CHECKED_MALLOC_ARRAY(double, state->nz_parts, "z partitions");

  /* Calculate aspect ratios so that subvolumes are approximately cubic */
  double x_aspect = (part_max->x - part_min->x) / f_max;
  double y_aspect = (part_max->y - part_min->y) / f_max;
  double z_aspect = (part_max->z - part_min->z) / f_max;

  int x_in = (int)floor((state->nx_parts - 2) * x_aspect + 0.5);
  int y_in = (int)floor((state->ny_parts - 2) * y_aspect + 0.5);
  int z_in = (int)floor((state->nz_parts - 2) * z_aspect + 0.5);
  if (x_in < 2)
    x_in = 2;
  if (y_in < 2)
    y_in = 2;
  if (z_in < 2)
    z_in = 2;

  /* If we've violated our 2*reaction radius criterion, fix it */
  smallest_spacing = 2 * state->rx_radius_3d;
  if ((part_max->x - part_min->x) / (x_in - 1) < smallest_spacing) {
    x_in = 1 + (int)floor((part_max->x - part_min->x) / smallest_spacing);
  }
  if ((part_max->y - part_min->y) / (y_in - 1) < smallest_spacing) {
    y_in = 1 + (int)floor((part_max->y - part_min->y) / smallest_spacing);
  }
  if ((part_max->z - part_min->z) / (z_in - 1) < smallest_spacing) {
    z_in = 1 + (int)floor((part_max->z - part_min->z) / smallest_spacing);
  }

  /* Set up to walk symmetrically out from the center of the state, dropping
   * partitions on the way */
  if (x_in < 2)
    x_in = 2;
  if (y_in < 2)
    y_in = 2;
  if (z_in < 2)
    z_in = 2;
  int x_start = (state->nx_parts - x_in) / 2;
  int y_start = (state->ny_parts - y_in) / 2;
  int z_start = (state->nz_parts - z_in) / 2;
  if (x_start < 1)
    x_start = 1;
  if (y_start < 1)
    y_start = 1;
  if (z_start < 1)
    z_start = 1;

  set_fineparts(part_min->x, part_max->x, state->x_partitions,
                state->x_fineparts, state->nx_parts, x_in, x_start);
  set_fineparts(part_min->y, part_max->y, state->y_partitions,
                state->y_fineparts, state->ny_parts, y_in, y_start);
  set_fineparts(part_min->z, part_max->z, state->z_partitions,
                state->z_fineparts, state->nz_parts, z_in, z_start);
}

void set_fineparts(double min, double max, double *partitions,
                   double *fineparts, int n_parts, int in, int start) {
  /* Now go through and drop partitions in each direction (picked from
   * sensibly close fine partitions) */
  double f = (max - min) / (in - 1);
  int j = 0;
  partitions[0] = fineparts[1];
  /* Dunno how this actually works! */
  for (int i = start; i < start + in; i++) {
    partitions[i] = fineparts[4096 + (i - start) * 16384 / (in - 1)];
  }
  for (int i = start - 1; i > 0; i--) {
    for (j = 0; partitions[i + 1] - fineparts[4095 - j] < f; j++) {
    }
    partitions[i] = fineparts[4095 - j];
  }
  for (int i = start + in; i < n_parts - 1; i++) {
    for (j = 0; fineparts[4096 + 16384 + j] - partitions[i - 1] < f; j++) {
    }
    partitions[i] = fineparts[4096 + 16384 + j];
  }
  partitions[n_parts - 1] = fineparts[4096 + 16384 + 4096 - 2];
}

void set_user_partitions(struct volume *state, double dfx, double dfy,
                         double dfz) {
  // User-supplied partitions
  // We need to keep the outermost partition away from the state bounding box.
  // We do this by adding a larger outermost partition, calculated somehow or
  // other.

  dfx += 1e-3;
  dfy += 1e-3;
  dfz += 1e-3;

  state->x_partitions =
      add_extra_outer_partitions(state->x_partitions, state->bb_llf.x,
                                 state->bb_urb.x, dfx, &state->nx_parts);
  state->y_partitions =
      add_extra_outer_partitions(state->y_partitions, state->bb_llf.y,
                                 state->bb_urb.y, dfy, &state->ny_parts);
  state->z_partitions =
      add_extra_outer_partitions(state->z_partitions, state->bb_llf.z,
                                 state->bb_urb.z, dfz, &state->nz_parts);

  find_closest_fine_part(state->x_partitions, state->x_fineparts,
                         state->n_fineparts, state->nx_parts);
  find_closest_fine_part(state->y_partitions, state->y_fineparts,
                         state->n_fineparts, state->ny_parts);
  find_closest_fine_part(state->z_partitions, state->z_fineparts,
                         state->n_fineparts, state->nz_parts);
}

void find_closest_fine_part(double *partitions, double *fineparts,
                            int n_fineparts, int n_parts) {
  /* Now that we've added outermost partitions, we find the closest fine
   * partition along each axis */
  partitions[0] = fineparts[1];
  for (int i = 1; i < n_parts - 1; i++) {
    partitions[i] =
        fineparts[bisect_near(fineparts, n_fineparts, partitions[i])];
  }
  partitions[n_parts - 1] = fineparts[4096 + 16384 + 4096 - 2];
}

double *add_extra_outer_partitions(double *partitions, double bb_llf_val,
                                   double bb_urb_val, double df_val,
                                   int *n_parts) {
  /* All this code just adds extra outermost partitions if they might be too
   * close to the outermost objects in the state */
  /* Don't ask me how it actually does it (or if it does it successfully...) */
  double *dbl_array;
  if (partitions[1] + df_val > bb_llf_val) {
    if (partitions[1] - df_val < bb_llf_val)
      partitions[1] = bb_llf_val - df_val;
    else {
      dbl_array = CHECKED_MALLOC_ARRAY(double, (*n_parts + 1),
                                       "x partitions (expanded in -X dir)");
      dbl_array[0] = partitions[0];
      dbl_array[1] = bb_llf_val - df_val;
      memcpy(&(dbl_array[2]), &(partitions[1]),
             sizeof(double) * (*n_parts - 1));
      free(partitions);
      partitions = dbl_array;
      *n_parts = *n_parts + 1;
    }
  }
  if (partitions[*n_parts - 2] - df_val < bb_urb_val) {
    if (partitions[*n_parts - 2] + df_val > bb_urb_val)
      partitions[*n_parts - 2] = bb_urb_val + df_val;
    else {
      dbl_array = CHECKED_MALLOC_ARRAY(double, (*n_parts + 1),
                                       "x partitions (expanded in +X dir)");
      dbl_array[*n_parts] = partitions[*n_parts - 1];
      dbl_array[*n_parts - 1] = bb_urb_val + df_val;
      memcpy(dbl_array, partitions, sizeof(double) * (*n_parts - 1));
      free(partitions);
      partitions = dbl_array;
      *n_parts = *n_parts + 1;
    }
  }
  return partitions;
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
void path_bounding_box(struct vector3 *loc, struct vector3 *displacement,
                       struct vector3 *llf, struct vector3 *urb,
                       double rx_radius_3d) {
  struct vector3 final; /* final position of the molecule after random walk */
  double R;             /* molecule interaction radius */

  R = rx_radius_3d;
  vect_sum(loc, displacement, &final);

  llf->x = urb->x = loc->x;
  llf->y = urb->y = loc->y;
  llf->z = urb->z = loc->z;

  if (final.x < llf->x) {
    llf->x = final.x;
  }
  if (final.x > urb->x) {
    urb->x = final.x;
  }
  if (final.y < llf->y) {
    llf->y = final.y;
  }
  if (final.y > urb->y) {
    urb->y = final.y;
  }
  if (final.z < llf->z) {
    llf->z = final.z;
  }
  if (final.z > urb->z) {
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
 collect_molecule:
    Perform garbage collection on a discarded molecule.  If the molecule is no
    longer in any lists, it will be freed.

 In: vm: the molecule
 Out: Nothing.  Molecule is unlinked from its list in the subvolume, and
      possibly returned to its birthplace.
***************************************************************************/
void collect_molecule(struct volume_molecule *vm) {
  /* Unlink from the previous item */
  if (vm->prev_v != NULL) {
#ifdef DEBUG_LIST_CHECKS
    if (*vm->prev_v != vm) {
      mcell_error_nodie("Stale previous pointer!  ACK!  THRBBPPPPT!");
    }
#endif
    *(vm->prev_v) = vm->next_v;
  }

  /* Unlink from the following item */
  if (vm->next_v != NULL) {
#ifdef DEBUG_LIST_CHECKS
    if (vm->next_v->prev_v != &vm->next_v) {
      mcell_error_nodie("Stale next pointer!  ACK!  THRBBPPPPT!");
    }
#endif
    vm->next_v->prev_v = vm->prev_v;
  }

  /* Clear our next/prev pointers */
  vm->prev_v = NULL;
  vm->next_v = NULL;

  /* Dispose of the molecule */
  vm->properties = NULL;
  vm->flags &= ~IN_VOLUME;
  if ((vm->flags & IN_MASK) == 0)
    mem_put(vm->birthplace, vm);
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
     vm: the molecule
 Out: Nothing.  Molecule is added to the subvolume's molecule lists.
***************************************************************************/
void ht_add_molecule_to_list(struct pointer_hash *h,
                             struct volume_molecule *vm) {
  struct per_species_list *list = NULL;

  /* See if we have a list */
  list = (struct per_species_list *)pointer_hash_lookup(h, vm->properties,
                                                        vm->properties->hashval);

  /* If not, create one and add it in */
  if (list == NULL) {
    list = (struct per_species_list *)CHECKED_MEM_GET(
        vm->subvol->local_storage->pslv, "per-species molecule list");
    list->properties = vm->properties;
    list->head = NULL;
    if (pointer_hash_add(h, vm->properties, vm->properties->hashval, list))
      mcell_allocfailed("Failed to add species to subvolume species table.");

    list->next = vm->subvol->species_head;
    vm->subvol->species_head = list;
  }

  /* Link the molecule into the list */
  vm->next_v = list->head;
  if (list->head)
    list->head->prev_v = &vm->next_v;
  vm->prev_v = &list->head;
  list->head = vm;
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
void ht_remove(struct pointer_hash *h, struct per_species_list *psl) {
  struct species *s = psl->properties;
  if (s == NULL)
    return;

  (void)pointer_hash_remove(h, s, s->hashval);
}

/***************************************************************************
 test_max_release:

 In: num_to_release: The number to release
     name: The of the release site
 Out: The number to release
***************************************************************************/
static int test_max_release(double num_to_release, char *name) {
  int num = (int)num_to_release;
  if (num > INT_MAX)
    mcell_error("Release site \"%s\" tries to release more than INT_MAX "
                "(2147483647) molecules.",
                name);
  return num;
}

/***************************************************************************
 check_release_probability:

 In: release_prob: the probability of release
     state: MCell state
     req: release event
     rpat: release pattern
 Out: Return 1 if release probability is < k (random number). Otherwise 0.
***************************************************************************/
static int check_release_probability(double release_prob, struct volume *state,
                                     struct release_event_queue *req,
                                     struct release_pattern *rpat) {
  /* check whether the release will happen */
  if (release_prob < 1.0) {
    double k = rng_dbl(state->rng);
    if (release_prob < k) {
      /* make sure we will try the release pattern again in the future */
      req->event_time += rpat->release_interval;

      /* we may need to move to the next train. */
      if (!distinguishable(req->event_time,
                           req->train_high_time + rpat->train_duration,
                           EPS_C) ||
          req->event_time > req->train_high_time + rpat->train_duration) {
        req->train_high_time += rpat->train_interval;
        req->event_time = req->train_high_time;
        req->train_counter++;
      }

      if (req->train_counter <= rpat->number_of_trains &&
          req->event_time < FOREVER) {
        if (schedule_add(state->releaser, req))
          mcell_allocfailed("Failed to add release request to scheduler.");
      }
      return 1;
    }
  }
  return 0;
}

/***************************************************************************
 check_past_events:
    Skip events that happened in the past (delay<0 or after checkpoint)

 In: release_prob: the probability of release
     state: MCell state
     req: release event
     rpat: release pattern
 Out: Return 1 if release pattern is not a reaction and release event hasn't
      happened yet
***************************************************************************/
static int skip_past_events(double release_prob, struct volume *state,
                            struct release_event_queue *req,
                            struct release_pattern *rpat) {

  if (req->event_time < state->current_iterations &&
      (distinguishable(release_prob, MAGIC_PATTERN_PROBABILITY, EPS_C))) {
    do {
      /* Schedule next release event and leave the function.
         This part of the code is relevant to checkpointing. */
      if (release_prob < 1.0) {
        if (!distinguishable(release_prob, 0, EPS_C))
          return 0;
        req->event_time += rpat->release_interval;
      } else {
        req->event_time += rpat->release_interval;
      }
      /* we may need to move to the next train. */
      if (!distinguishable(req->event_time,
                           req->train_high_time + rpat->train_duration,
                           EPS_C) ||
          req->event_time > req->train_high_time + rpat->train_duration) {
        req->train_high_time += rpat->train_interval;
        req->event_time = req->train_high_time;
        req->train_counter++;
      }
    } while (req->event_time <= state->start_iterations);

    if (req->train_counter <= rpat->number_of_trains &&
        req->event_time < FOREVER) {
      if (schedule_add(state->releaser, req))
        mcell_allocfailed("Failed to add release request to scheduler.");
    }
    return 1;
  }
  return 0;
}

/***************************************************************************
 calculate_number_to_release:

 In:  rso: Release site object
      state: MCell state
 Out: number: The number of molecules to be released
***************************************************************************/
static int calculate_number_to_release(struct release_site_obj *rso,
                                       struct volume *state) {
  int number;
  double vol, num_to_release;
  switch (rso->release_number_method) {

  case CONSTNUM:
    num_to_release = rso->release_number;
    number = test_max_release(num_to_release, rso->name);
    break;

  case GAUSSNUM:
    if (rso->standard_deviation > 0) {
      num_to_release = (rng_gauss(state->rng) * rso->standard_deviation +
                        rso->release_number);
      number = test_max_release(num_to_release, rso->name);
    } else {
      rso->release_number_method = CONSTNUM;
      num_to_release = rso->release_number;
      number = test_max_release(num_to_release, rso->name);
    }
    break;

  case VOLNUM: {
    double diam = rso->mean_diameter;
    if (rso->standard_deviation > 0) {
      diam += rng_gauss(state->rng) * rso->standard_deviation;
    }
    vol = (MY_PI / 6.0) * diam * diam * diam;
    num_to_release = N_AV * 1e-15 * rso->concentration * vol + 0.5;
    number = test_max_release(num_to_release, rso->name);
    break;
  }

  case CCNNUM:
  case DENSITYNUM:
    if (rso->diameter == NULL)
      number = 0;
    else {
      switch (rso->release_shape) {
      case SHAPE_SPHERICAL:
      case SHAPE_ELLIPTIC:
        vol = (1.0 / 6.0) * MY_PI * rso->diameter->x * rso->diameter->y *
              rso->diameter->z;
        break;
      case SHAPE_RECTANGULAR:
      case SHAPE_CUBIC:
        vol = rso->diameter->x * rso->diameter->y * rso->diameter->z;
        break;

      case SHAPE_SPHERICAL_SHELL:
        mcell_error("Release site \"%s\" tries to release a concentration on a "
                    "spherical shell.",
                    rso->name);
        /*vol = 0;*/
        /*break;*/

      default:
        mcell_internal_error("Release by concentration on invalid release site "
                             "shape (%d) for release site \"%s\".",
                             rso->release_shape, rso->name);
        /*break;*/
      }
      num_to_release = N_AV * 1e-15 * rso->concentration * vol *
                           state->length_unit * state->length_unit *
                           state->length_unit +
                       0.5;
      number = test_max_release(num_to_release, rso->name);
    }
    break;

  default:
    mcell_internal_error(
        "Release site \"%s\" has invalid release number method (%d).",
        rso->name, rso->release_number_method);
    /*number = 0;*/
    /*break;*/
  }
  return number;
}

/*
 * num_vol_mols_from_conc computes the number of volume molecules to be
 * released within a closed object. There are two cases:
 * - for a single closed object we know the exact volume and can thus compute
 *   the exact number of molecules required and release them by setting
 *   exactNumber to true.
 * - for a release object consisting of a boolean expression of closed objects
 *   we are currently not able to compute the volume exactly. Instead we compute
 *   the number of molecules in the bounding box and then release an approximate
 *   number by setting exactNumber to false.
 */
int num_vol_mols_from_conc(struct release_site_obj *rso, double length_unit,
                           bool *exactNumber) {

  struct release_region_data *rrd = rso->region_data;
  struct release_evaluator *eval = rrd->expression;
  double vol = 0.0;
  if (eval->left != NULL && (eval->op & REXP_LEFT_REGION) &&
      eval->right == NULL && (eval->op & REXP_NO_OP)) {
    struct region *r = (struct region *)eval->left;
    assert(r->manifold_flag !=
           MANIFOLD_UNCHECKED); // otherwise we have no volume
    vol = r->volume;
    *exactNumber = true;
  } else {
    vol = (rrd->urb.x - rrd->llf.x) * (rrd->urb.y - rrd->llf.y) *
          (rrd->urb.z - rrd->llf.z);
  }
  double num_to_release = (N_AV * 1e-15 * rso->concentration * vol *
                           length_unit * length_unit * length_unit) +
                          0.5;
  return test_max_release(num_to_release, rso->name);
}

/*************************************************************************
  periodic_boxes_are_identical() tests if two periodic boxes are identical

 In:  b1: pointer to first periodic_box struct
      b2: pointer to second periodic_box struct
 Out: true if the two boxes are identical and false otherwise.
      NOTE: If one or both of b1 or b2 are NULL we also return true.
      This behavior makes sure that e.g. a COUNT without any periodic
      box defined always matches any other periodic box.
*************************************************************************/
bool periodic_boxes_are_identical(const struct periodic_image *b1,
  const struct periodic_image *b2) {
  if (b1 == NULL || b2 == NULL) {
    return true;
  }
  return (b1->x == b2->x) && (b1->y == b2->y) && (b1->z == b2->z);
}

/*************************************************************************
  convert_relative_to_abs_PBC_coords is used to convert the PBC coordinate
  system. This is probably not the best name for this function, since the
  meaning of relative and absolute are somewhat ambiguous in this context.
  Note: This is only needed for the non-traditional form of PBCs.

 In:  periodic_box_obj: The actual periodic box object
      periodic_box: The current periodic box that the molecule is in
      periodic_traditional: A flag to indicate whether we are using traditional
                            PBCs or mirrored geometry PBCs
      pos: The position of the molecule prior to conversion
      pos_output: The position of the molecule after conversion
 Out: If 0, then coordinates were successfully converted. If 1, then
      coordinates were not or did not need to be converted.
*************************************************************************/
int convert_relative_to_abs_PBC_coords(
    struct object *periodic_box_obj,
    struct periodic_image *periodic_box,
    bool periodic_traditional,
    struct vector3 *pos,
    struct vector3 *pos_output) {
  double llx = 0.0;
  double urx = 0.0;
  double lly = 0.0;
  double ury = 0.0;
  double llz = 0.0;
  double urz = 0.0;
  double x_box_length = 0.0;
  double y_box_length = 0.0;
  double z_box_length = 0.0;
  if (periodic_box_obj && !(periodic_traditional)) {
    assert(periodic_box_obj->object_type == BOX_OBJ);
    struct polygon_object* p = (struct polygon_object*)(periodic_box_obj->contents);
    struct subdivided_box* sb = p->sb;
    x_box_length = sb->x[1] - sb->x[0];
    y_box_length = sb->y[1] - sb->y[0];
    z_box_length = sb->z[1] - sb->z[0];
    llx = sb->x[0];
    urx = sb->x[1];
    lly = sb->y[0];
    ury = sb->y[1];
    llz = sb->z[0];
    urz = sb->z[1];
  }

  if (periodic_box_obj && !(periodic_traditional)) {

    int pos_or_neg = (periodic_box->x > 0) ? 1 : -1;
    double difference = (periodic_box->x > 0) ? urx - pos->x : pos->x - llx;

    // translate X
    if (periodic_box->x == 0) {
      pos_output->x = pos->x; 
    }
    else if (periodic_box->x % 2 == 0) {
      pos_output->x = pos->x + pos_or_neg * (abs(periodic_box->x) * x_box_length);
    }
    else {
      pos_output->x = pos->x + pos_or_neg * ((abs(periodic_box->x) - 1) * x_box_length + 2 * difference);
    }

    // translate Y
    pos_or_neg = (periodic_box->y > 0) ? 1 : -1;
    difference = (periodic_box->y > 0) ? ury - pos->y : pos->y - lly;
    if (periodic_box->y == 0) {
      pos_output->y = pos->y; 
    }
    else if (periodic_box->y % 2 == 0) {
      pos_output->y = pos->y + pos_or_neg * (abs(periodic_box->y) * y_box_length);
    }
    else {
      pos_output->y = pos->y + pos_or_neg * ((abs(periodic_box->y) - 1) * y_box_length + 2 * difference);
    }

    // translate Z
    pos_or_neg = (periodic_box->z > 0) ? 1 : -1;
    difference = (periodic_box->z > 0) ? urz - pos->z : pos->z - llz;
    if (periodic_box->z == 0) {
      pos_output->z = pos->z; 
    }
    else if (periodic_box->z % 2 == 0) {
      pos_output->z = pos->z + pos_or_neg * (abs(periodic_box->z) * z_box_length);
    }
    else {
      pos_output->z = pos->z + pos_or_neg * ((abs(periodic_box->z) - 1) * z_box_length + 2 * difference);
    }

    return 0;
  }
  else {
    return 1; 
  }
}

/*************************************************************************
  add_surfmol_with_unique_pb_to_list

 In:  sm_list: a list of surface molecules
      sm: the surface molecule we want to add to the list
 Out: Return the head of the list. Also, sm should be added to sm_list if the
      periodic box it inhabits isn't already in the list.
*************************************************************************/
struct surface_molecule_list* add_surfmol_with_unique_pb_to_list(
    struct surface_molecule_list *sm_list,
    struct surface_molecule *sm) {
  struct surface_molecule_list *sm_list_head = sm_list;
  struct surface_molecule_list *sm_entry = CHECKED_MALLOC_STRUCT(
    struct surface_molecule_list, "surface molecule list");
  sm_entry->sm = sm;
  sm_entry->next = NULL;
  if (sm_list == NULL) {
    sm_list_head = sm_entry;
  }
  else if (sm_list->sm == NULL) {
    sm_list_head->sm = sm;
    free(sm_entry);
  }
  else {
    for (; sm_list != NULL; sm_list = sm_list->next) {
      if (sm && periodic_boxes_are_identical(
          sm_list->sm->periodic_box, sm->periodic_box)) {
        free(sm_entry);
        return NULL;
      }
      if (sm_list->next == NULL) {
        sm_list->next = sm_entry;
        break;
      }
    }
  }
  return sm_list_head;
}

/*************************************************************************
  remove_surfmol_from_list

 In:  sm_head: pointer to the head of a list of surface molecules
      sm: the surface molecule we want to remove from the list
 Out: Remove sm from the surface molecule list. Return 1 on failure, 0
      otherwise.
*************************************************************************/
void remove_surfmol_from_list(
    struct surface_molecule_list **sm_head,
    struct surface_molecule *sm) {

  struct surface_molecule_list *sm_list = *sm_head;
  struct surface_molecule_list *prev = *sm_head;

  if (sm_list == NULL) {
  }
  else if (sm_list->sm == sm) {
    if (sm_list->next != NULL) {
      *sm_head = sm_list->next; 
    }
    else {
      *sm_head = NULL;
    }
    free(sm_list);
  }
  else {
    for (; sm_list != NULL; sm_list = sm_list->next) {
      if (sm_list->sm == sm) {
        prev->next = sm_list->next; 
        free(sm_list);
        sm_list = NULL;
        break;
      }
      prev = sm_list;
    }
  }
  return;
}
