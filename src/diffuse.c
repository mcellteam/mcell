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
#include <stdlib.h>

#include "diffuse.h"
#include "logging.h"
#include "mcell_structs.h"
#include "count_util.h"
#include "grid_util.h"
#include "vol_util.h"
#include "wall_util.h"
#include "react.h"


#define FREE_COLLISION_LISTS()                                                 \
  do {                                                                         \
    if (shead2 != NULL)                                                        \
      mem_put_list(sv->local_storage->coll, shead2);                           \
    if (shead != NULL)                                                         \
      mem_put_list(sv->local_storage->coll, shead);                            \
  } while (0)


static const int inert_to_mol = 1;
static const int inert_to_all = 2;

/* declaration of static functions */
int move_sm_on_same_triangle(
    struct volume *state,
    struct surface_molecule *sm,
    struct vector2 *new_loc,
    struct periodic_image *previous_box,
    struct wall *new_wall,
    struct hit_data *hd_info);

static void redo_collision_list(struct volume* world, struct collision** shead,
  struct collision** stail, struct collision** shead_exp,
  struct volume_molecule* m, struct vector3* displacement, struct subvolume* sv);

static int collide_and_react_with_vol_mol(
  struct volume* world, struct collision* smash, struct volume_molecule* vm,
  struct collision** tentative, struct vector3* displacement,
  struct vector3* loc_certain, double t_steps, double r_rate_factor);

static int collide_and_react_with_surf_mol(
  struct volume* world, struct collision* smash, struct volume_molecule* vm,
  struct collision** tentative, struct vector3** loc_certain, double t_steps,
  int mol_grid_flag, int mol_mol_grid_flag, double r_rate_factor);

static int collide_and_react_with_walls(
  struct volume* world, struct collision* smash, struct volume_molecule* vm,
  struct collision** tentative, struct vector3** loc_certain, double t_steps,
  int inertness, double r_rate_factor);

static int reflect_or_periodic_bc(
  struct volume* world, struct collision* smash, struct vector3* displacement,
  struct volume_molecule** mol, struct wall** reflectee,
  struct collision** tentative, double* t_steps);

static void reflect_absorb_inside_out(
    struct volume *world, struct surface_molecule *sm, struct hit_data *hd_head,
    struct rxn **rx, struct rxn *matching_rxns[], struct vector2 boundary_pos,
    struct wall *this_wall, int index_edge_was_hit, int *reflect_now,
    int *absorb_now, int *this_wall_edge_region_border);

int reflect_absorb_outside_in( 
    struct volume *world,
    struct surface_molecule *sm,
    struct hit_data **hd_head,
    struct rxn **rx,
    struct rxn *matching_rxns[],
    struct vector2 boundary_pos,
    struct wall *target_wall,
    struct wall *this_wall,
    int *reflect_now,
    int *absorb_now,
    int this_wall_edge_region_border);

void collide_and_react_with_subvol(
  struct volume* world, struct collision *smash, struct vector3* displacement,
  struct volume_molecule** mol, struct collision** tentative, double* t_steps);

void compute_displacement(
  struct volume* world, struct collision* shead, struct volume_molecule* vm,
  struct vector3* displacement, struct vector3* displacement2,
  double* rate_factor, double* r_rate_factor, double* steps, double* t_steps,
  double max_time);

void determine_mol_mol_reactions(
  struct volume* world, struct volume_molecule* vm, struct collision** shead,
  struct collision** stail, int interness);

void set_inertness_and_maxtime(
  struct volume* world, struct volume_molecule* vm, double* maxtime,
  int* inertness);

void register_hits(
  struct volume* world, struct volume_molecule* m,
  struct collision** tentative, struct wall** reflect_w, double* reflect_t,
  struct vector3* displacement, struct collision* smash, double* t_steps);

void count_tentative_collisions(
  struct volume *world, struct collision **tc, struct collision *smash,
  struct species *spec, double t_confident, int destroy_flag,
  struct periodic_image *box);

void change_boxes_2D(
    bool periodic_traditional,
    struct surface_molecule *sm,
    struct object *periodic_box_obj,
    struct vector3 *hit_xyz,
    struct vector3 *teleport_xyz);

/*************************************************************************
pick_2D_displacement:
  In: v: vector2 to store the new displacement
      scale: scale factor to apply to the displacement
      rng:
  Out: No return value.  vector is set to a random orientation and a
         distance chosen from the probability distribution of a diffusing
         2D molecule, scaled by the scaling factor.
*************************************************************************/
void pick_2D_displacement(struct vector2 *v, double scale,
                          struct rng_state *rng) {
  static const double one_over_2_to_16th = 1.52587890625e-5;
  struct vector2 a;

  /*
   * NOTE: The below algorithm is the polar method due to Marsaglia
   * combined with a rejection method for picking uniform random
   * variates in C2.
   * Both methods are nicely described in Chapters V.4.3 and V.4.4
   * of "Non-Uniform Random Variate Generation" by Luc Devroye
   * (http://luc.devroye.org/rnbookindex.html).
   */
  double f;
  do {
    unsigned int n = rng_uint(rng);

    a.u = 2.0 * one_over_2_to_16th * (n & 0xFFFF) - 1.0;
    a.v = 2.0 * one_over_2_to_16th * (n >> 16) - 1.0;
    f = a.u * a.u + a.v * a.v;
  } while ((f < EPS_C) || (f > 1.0));

  /*
   * NOTE: The scaling factor to go from a uniform to
   * a normal distribution is sqrt(-2log(f)/f).
   * However, since we use two normally distributed
   * variates to generate a normally distributed
   * 2d vector (with variance 1) we have to normalize
   * and divide by an additional factor of sqrt(2)
   * resulting in normalFactor.
   */
  double normalFactor = sqrt(-log(f) / f);
  v->u = a.u * normalFactor * scale;
  v->v = a.v * normalFactor * scale;
}

/*************************************************************************
pick_clamped_displacement:
  In: v: vector3 to store the new displacement
      vm: molecule that just came through the surface
      rng:
      radial_subdivisions:
  Out: No return value.  vector is set to a random orientation and a
         distance chosen from the probability distribution of a diffusing
         3D molecule that has come through a surface from a uniform
         concentration on the other side.
  Note: vm->previous_wall points to the wall we're coming from, and
        vm->index is the orientation we came off with
*************************************************************************/
void pick_clamped_displacement(struct vector3 *v, struct volume_molecule *vm,
                               double *r_step_surface, struct rng_state *rng,
                               u_int radial_subdivisions) {
  static const double one_over_2_to_20th = 9.5367431640625e-7;
  struct wall *w = vm->previous_wall;

  unsigned int n = rng_uint(rng);

  /* Correct distribution along normal from surface (from lookup table) */
  double r_n = r_step_surface[n & (radial_subdivisions - 1)];

  double p = one_over_2_to_20th * ((n >> 12) + 0.5);
  double t = r_n / erfcinv(p * erfc(r_n));
  struct vector2 r_uv;
  pick_2D_displacement(&r_uv, sqrt(t) * vm->properties->space_step, rng);

  r_n *= vm->index * vm->properties->space_step;
  v->x = r_n * w->normal.x + r_uv.u * w->unit_u.x + r_uv.v * w->unit_v.x;
  v->y = r_n * w->normal.y + r_uv.u * w->unit_u.y + r_uv.v * w->unit_v.y;
  v->z = r_n * w->normal.z + r_uv.u * w->unit_u.z + r_uv.v * w->unit_v.z;
}

/*************************************************************************
pick_release_displacement:
  In: in_disk: vector3 to store the position on interaction disk to go to
      away: vector3 along which to travel away from the disk
      scale:
      r_step_release:
      d_step:
      radial_subdivisions:
      directions_mask:
      num_directions:
      rx_radius_3d:
      rng:
  Out: No return value.  Vectors are set to random orientation with
         distances chosen from the probability distribution matching
         the binding of a 3D molecule (distance and direction to
         interaction disk and distance along disk).
*************************************************************************/
void pick_release_displacement(struct vector3 *in_disk, struct vector3 *away,
                               double scale, double *r_step_release,
                               double *d_step, u_int radial_subdivisions,
                               int directions_mask, u_int num_directions,
                               double rx_radius_3d, struct rng_state *rng) {
  static const double one_over_2_to_16th = 1.52587890625e-5;
  struct vector2 disk;
  struct vector3 orth, axo;

  u_int bits = rng_uint(rng);

  u_int x_bit = (bits & 0x80000000);
  u_int y_bit = (bits & 0x40000000);
  u_int z_bit = (bits & 0x20000000);
  u_int thetaphi_bits = (bits & 0x1FFFF000) >> 12;
  u_int r_bits = (bits & 0x00000FFF);

  double r = scale * r_step_release[r_bits & (radial_subdivisions - 1)];

  u_int idx = thetaphi_bits & directions_mask;
  while (idx >= num_directions) {
    idx = rng_uint(rng) & directions_mask;
  }

  if (x_bit)
    away->x = d_step[idx];
  else
    away->x = -d_step[idx];
  if (y_bit)
    away->y = d_step[idx + 1];
  else
    away->y = -d_step[idx + 1];
  if (z_bit)
    away->z = d_step[idx + 2];
  else
    away->z = -d_step[idx + 2];

  if (d_step[idx] < d_step[idx + 1]) {
    if (d_step[idx] < d_step[idx + 2]) {
      orth.x = 0;
      orth.y = away->z;
      orth.z = -away->y;
    } else {
      orth.x = away->y;
      orth.y = -away->x;
      orth.z = 0;
    }
  } else if (d_step[idx + 1] < d_step[idx + 2]) {
    orth.x = away->z;
    orth.y = 0;
    orth.z = -away->x;
  } else {
    orth.x = away->y;
    orth.y = -away->x;
    orth.z = 0;
  }

  normalize(&orth);
  cross_prod(away, &orth, &axo);

  double f;
  do {
    bits = rng_uint(rng);

    disk.u = 2.0 * one_over_2_to_16th * (bits & 0xFFFF) - 1.0;
    disk.v = 2.0 * one_over_2_to_16th * (bits >> 16) - 1.0;
    f = disk.u * disk.u + disk.v * disk.v;
  } while (f < 0.01 || f > 1.0);

  in_disk->x = (disk.u * orth.x + disk.v * axo.x) * rx_radius_3d;
  in_disk->y = (disk.u * orth.y + disk.v * axo.y) * rx_radius_3d;
  in_disk->z = (disk.u * orth.z + disk.v * axo.z) * rx_radius_3d;

  away->x *= r;
  away->y *= r;
  away->z *= r;
}

/*************************************************************************
pick_displacement:
  In: vector3 to store the new displacement
      scale factor to apply to the displacement
  Out: No return value.  vector is set to a random orientation and a
         distance chosen from the probability distribution of a diffusing
         3D molecule, scaled by the scaling factor.
*************************************************************************/
void pick_displacement(struct vector3 *v, double scale, struct rng_state *rng) {
  v->x = scale * rng_gauss(rng) * .70710678118654752440;
  v->y = scale * rng_gauss(rng) * .70710678118654752440;
  v->z = scale * rng_gauss(rng) * .70710678118654752440;
}

/*************************************************************************
reflect_periodic_2D:
  In: state: simulation state
      index_edge_was_hit: the index of the edge that we might have hit
      origin_uv: uv coordinates of where the SM started
      curr_wall: the wall that the SM is currently on
      disp_uv: uv coordinates of the displacement
      boundary_uv: uv coordinates of the edge boundary.
      origin_xyz: uv coordinates of where the SM startedi
  Note: for origin_xyz, we used to compute it in this function from origin_uv,
        but it can be dangerous if we've just hit a PB, since the conversions
        back and forth from uv to xyz led to precision errors. 
  Out: the following things are updated if we hit the periodic box:
       disp_uv will be reversed if we hit. should do a proper reflection
       return the XYZ loc if we hit the periodic box, NULL otherwise.
*************************************************************************/
struct vector3* reflect_periodic_2D(
    struct volume *state,
    int index_edge_was_hit,
    struct vector2 *origin_uv,
    struct wall *curr_wall,
    struct vector2 *disp_uv,
    struct vector2 *boundary_uv,
    struct vector3 *origin_xyz) {

  struct vector3 target_xyz;
  struct vector2 target_uv = {
    .u = origin_uv->u + disp_uv->u,
    .v = origin_uv->v + disp_uv->v
  };
  // Still within current wall, but we might have hit the periodic box
  // set the target_xyz
  if (index_edge_was_hit == -1) {
    uv2xyz(&target_uv, curr_wall, &target_xyz);
  }
  // Hit the edge of current wall, and we might have hit the periodic box
  else if (index_edge_was_hit == 0 || 
           index_edge_was_hit == 1 ||
           index_edge_was_hit == 2) {
    uv2xyz(boundary_uv, curr_wall, &target_xyz);
  }
  // It's unclear what we hit
  else {
    return NULL;
  }
  struct vector3 delta_xyz = {target_xyz.x - origin_xyz->x,
                              target_xyz.y - origin_xyz->y,
                              target_xyz.z - origin_xyz->z};
  struct vector3 updated_xyz = *origin_xyz;
  // Go through all the subvolumes between where we are to where we want to be
  for (struct subvolume *sv = find_subvolume(state, origin_xyz, NULL);
       sv != NULL; sv = next_subvol(
          &updated_xyz, &delta_xyz, sv, state->x_fineparts, state->y_fineparts,
          state->z_fineparts, state->ny_parts,
          state->nz_parts)) {

    // Check all the walls in this subvolume
    for (struct wall_list *wl = sv->wall_head; wl != NULL; wl = wl->next) {
      // Skip it if it's not part of the periodic box
      if (wl->this_wall->parent_object != state->periodic_box_obj) {
        continue;
      }

      struct vector3 *hit_xyz = malloc(sizeof(*hit_xyz));
      double t = 0.0;
      int i = collide_wall(
          &updated_xyz, &delta_xyz, wl->this_wall, &t, hit_xyz, 0, state->rng,
          state->notify, &(state->ray_polygon_tests));
      if (i != COLLIDE_MISS &&
          (hit_xyz->x - target_xyz.x) * delta_xyz.x +
          (hit_xyz->y - target_xyz.y) * delta_xyz.y +
          (hit_xyz->z - target_xyz.z) * delta_xyz.z < 0) {
        if (!state->periodic_traditional) {
          struct vector2 hit_uv;
          xyz2uv(hit_xyz, curr_wall, &hit_uv);
          // Take the remaining displacement and reverse it. Raytrace backwards
          // from the point we hit. Not a proper reflection but okay for now.
          disp_uv->u = -(target_uv.u-hit_uv.u);
          disp_uv->v = -(target_uv.v-hit_uv.v);
          origin_uv->u = hit_uv.u;
          origin_uv->v = hit_uv.v;
        }
        return hit_xyz;
      }
      free(hit_xyz);
    }
  }
  return NULL;
}
        
/*************************************************************************
change_boxes_2D:
  In: sm: the surface molecule
      periodic_box_obj: the actual periodic box object (*not* periodic_image)
      hit_xyz: where we hit on the periodic box
  Out: If we are using traditional PBCs, then we set the teleport_xyz to the
       opposite side of the box, which is where the molecule will move to.
       Othwerise, the periodic box on the sm is updated (i.e. sm->periodic_box)
  Todo: this function should probably be renamed, as it doesn't just change the
        periodic box.
*************************************************************************/
void change_boxes_2D(
    bool periodic_traditional,
    struct surface_molecule *sm,
    struct object *periodic_box_obj,
    struct vector3 *hit_xyz,
    struct vector3 *teleport_xyz) {

  assert(periodic_box_obj->object_type == BOX_OBJ);

  struct polygon_object* p = (struct polygon_object*)(periodic_box_obj->contents);
  struct subdivided_box* sb = p->sb;

  // Lower left and upper right corners of the periodic box
  double llx = sb->x[0];
  double urx = sb->x[1];
  double lly = sb->y[0];
  double ury = sb->y[1];
  double llz = sb->z[0];
  double urz = sb->z[1];

  int x_inc = (sm->periodic_box->x % 2 == 0) ? 1 : -1;
  int y_inc = (sm->periodic_box->y % 2 == 0) ? 1 : -1;
  int z_inc = (sm->periodic_box->z % 2 == 0) ? 1 : -1;
  int box_inc_x = 0;
  int box_inc_y = 0;
  int box_inc_z = 0;
  double x_pos = 0;
  double y_pos = 0;
  double z_pos = 0;

  if (!distinguishable(hit_xyz->x, llx, EPS_C)) {
    x_pos = urx - EPS_C;
    box_inc_x = -x_inc;
  } else if (!distinguishable(hit_xyz->x, urx, EPS_C)) {
    x_pos = llx + EPS_C;
    box_inc_x = x_inc;
  }
  // Wrap molecule around to other side of box
  if (periodic_traditional && x_pos) {
    teleport_xyz->x = x_pos;
  }

  if (!distinguishable(hit_xyz->y, lly, EPS_C)) {
    y_pos = ury - EPS_C;
    box_inc_y = -y_inc;
  } else if (!distinguishable(hit_xyz->y, ury, EPS_C)) {
    y_pos = lly + EPS_C;
    box_inc_y = y_inc;
  }
  // Wrap molecule around to other side of box
  if (periodic_traditional && y_pos) {
    teleport_xyz->y = y_pos;
  }

  if (!distinguishable(hit_xyz->z, llz, EPS_C)) {
    z_pos = urz - EPS_C;
    box_inc_z = -z_inc;
  } else if (!distinguishable(hit_xyz->z, urz, EPS_C)) {
    z_pos = llz + EPS_C;
    box_inc_z = z_inc;
  }
  // Wrap molecule around to other side of box
  if (periodic_traditional && z_pos) {
    teleport_xyz->z = z_pos;
  }

  if (!(periodic_traditional) && (box_inc_x || box_inc_y || box_inc_z)) {
    sm->periodic_box->x += box_inc_x;
    sm->periodic_box->y += box_inc_y;
    sm->periodic_box->z += box_inc_z;
  }
}

/*************************************************************************
reflect_absorb_inside_out:
  In: world: simulation state
      sm: molecule that is moving
      hd_head: region border hit data information
      rx: the type of reaction if any - absorptive/reflective
      matching_rxns: an array of possible reactions
      boundary_pos: the uv coordinates where we hit
      this_wall: the wall that we are on
      index_edge_was_hit: the index of the edge we just hit (0,1,2)
      reflect_now: should the sm reflect
      absorb_now: should the sm be absorbed
      this_wall_edge_region_border:
  Out: 1 if we are about to reflect or absorb (reflect_now, absorb_now). 0
       otherwise. hd_head and this_wall_edge_region_border are updated.
*************************************************************************/
void reflect_absorb_inside_out(
    struct volume *world,
    struct surface_molecule *sm,
    struct hit_data *hd_head,
    struct rxn **rx,
    struct rxn *matching_rxns[],
    struct vector2 boundary_pos,
    struct wall *this_wall,
    int index_edge_was_hit,
    int *reflect_now,
    int *absorb_now,
    int *this_wall_edge_region_border) {

  struct edge *this_edge = this_wall->edges[index_edge_was_hit];

  if (is_wall_edge_region_border(this_wall, this_edge)) {
    *this_wall_edge_region_border = 1;
  }

  /* find neighbor wall that shares this_edge and it's index in the coordinate
   * system of neighbor wall */
  struct wall *nbr_wall = NULL;
  /* index of the shared edge with neighbor wall in the coordinate system of
   * neighbor wall */
  int nbr_edge_ind = -1;
  find_neighbor_wall_and_edge(this_wall, index_edge_was_hit, &nbr_wall, &nbr_edge_ind);

  int nbr_wall_edge_region_border = 0;
  if (nbr_wall != NULL) {
    if (is_wall_edge_region_border(nbr_wall, nbr_wall->edges[nbr_edge_ind])) {
      nbr_wall_edge_region_border = 1;
    }
  }

  if (is_wall_edge_restricted_region_border(world, this_wall, this_edge, sm)) {
    int num_matching_rxns = trigger_intersect(
        world->reaction_hash, world->rx_hashsize, world->all_mols,
        world->all_volume_mols, world->all_surface_mols,
        sm->properties->hashval, (struct abstract_molecule *)sm, sm->orient,
        this_wall, matching_rxns, 1, 1, 1);

    /* check if this wall has any reflective or absorptive region borders for
     * this molecule (aka special reactions) */
    for (int i = 0; i < num_matching_rxns; i++) {
      *rx = matching_rxns[i];

      if ((*rx)->n_pathways == RX_REFLEC) {
        /* check for REFLECTIVE border */
        *reflect_now = 1;
        break;
      } else if ((*rx)->n_pathways == RX_ABSORB_REGION_BORDER) {
        /* check for ABSORPTIVE border */
        *absorb_now = 1;
        break;
      }
    }

    /* count hits if we absorb or reflect */
    if (reflect_now || absorb_now) {
      if (this_wall->flags & sm->properties->flags & COUNT_HITS) {
        // XXX: treat hd_head like we do in reflect_absorb_outside_in?
        update_hit_data(&hd_head, this_wall, this_wall, sm, boundary_pos, 1, 0);
      }

      if (nbr_wall != NULL && nbr_wall_edge_region_border) {
        if (nbr_wall->flags & sm->properties->flags & COUNT_HITS) {
          // XXX: treat hd_head like we do in reflect_absorb_outside_in?
          update_hit_data(&hd_head, this_wall, nbr_wall, sm, boundary_pos, 0, 0);
        }
      }
    }
  }
}

/*************************************************************************
reflect_absorb_outside_in:
  In: world: simulation state
      sm: molecule that is moving
      hd_head: region border hit data information
      rx: the type of reaction if any - absorptive/reflective
      matching_rxns: an array of possible reactions
      boundary_pos: the uv coordinates where we hit
      target_wall: the wall we hit
      this_wall: the wall that we are on
      reflect_now: should the sm reflect
      absorb_now: should the sm be absorbed
      this_wall_edge_region_border:
  Out: 1 if we are about to reflect or absorb (reflect_now, absorb_now). 0
       otherwise. hd_head is updated.
*************************************************************************/
int reflect_absorb_outside_in( 
    struct volume *world,
    struct surface_molecule *sm,
    struct hit_data **hd_head,
    struct rxn **rx,
    struct rxn *matching_rxns[],
    struct vector2 boundary_pos,
    struct wall *target_wall,
    struct wall *this_wall,
    int *reflect_now,
    int *absorb_now,
    int this_wall_edge_region_border) {

  /* index of the shared edge in the coordinate system of target wall */
  int target_edge_ind = find_shared_edge_index_of_neighbor_wall(this_wall, target_wall);

  int target_wall_edge_region_border = 0;
  if (is_wall_edge_region_border(target_wall, target_wall->edges[target_edge_ind])) {
    target_wall_edge_region_border = 1;
  }

  if (is_wall_edge_restricted_region_border(world, target_wall, target_wall->edges[target_edge_ind], sm)) {
    *reflect_now = 0;
    *absorb_now = 0;
    int num_matching_rxns = trigger_intersect(
        world->reaction_hash, world->rx_hashsize, world->all_mols,
        world->all_volume_mols, world->all_surface_mols,
        sm->properties->hashval, (struct abstract_molecule *)sm,
        sm->orient, target_wall, matching_rxns, 1, 1, 1);

    for (int i = 0; i < num_matching_rxns; i++) {
      *rx = matching_rxns[i];
      if ((*rx)->n_pathways == RX_REFLEC) {
        /* check for REFLECTIVE border */
        *reflect_now = 1;
        break;
      } else if ((*rx)->n_pathways == RX_ABSORB_REGION_BORDER) {
        /* check for ABSORPTIVE border */
        *absorb_now = 1;
        break;
      }
    }

    /* count hits if we reflect or absorb */
    if (*reflect_now || *absorb_now) {
      if (target_wall->flags & sm->properties->flags & COUNT_HITS) {
        /* this is OUTSIDE IN hit */
        update_hit_data(hd_head, this_wall, target_wall, sm, boundary_pos, 0, 0);

        /* this is INSIDE OUT hit for the same region border */
        update_hit_data(hd_head, this_wall, this_wall, sm, boundary_pos, 1, 0);
      }
    }

    if (*reflect_now || *absorb_now) {
      return 1; 
    }
  }

  if (this_wall_edge_region_border) {
    /* if we get to this point in the code the molecule crossed
       the region border inside out - update hits count */
    if (this_wall->flags & sm->properties->flags & COUNT_HITS) {
      update_hit_data(hd_head, this_wall, this_wall, sm, boundary_pos, 1, 1);
    }
  }
  if (target_wall_edge_region_border) {
    /* if we get to this point in the code the molecule crossed
       the region border outside in - update hits count */
    if (target_wall->flags & sm->properties->flags & COUNT_HITS) {
      update_hit_data(hd_head, this_wall, target_wall, sm, boundary_pos, 0, 1);
    }
  }
  return 0;
}

/*************************************************************************
ray_trace_2D:
  In: world: simulation state
      sm: molecule that is moving
      disp: displacement vector from current to new location
      pos: place to store new coordinate (in coord system of new wall)
      kill_me: flag that tells that molecule hits ABSORPTIVE region border
           (value = 1)
      rxp: reaction object (valid only in case of hitting ABSORPTIVE region
         border)
      hit_data_info: region border hit data information
  Out: Return wall at endpoint of movement vector. Otherwise NULL if we hit
       ambiguous location or if SM was absorbed.
       pos: location of that endpoint in the coordinate system of the new wall.
       kill_me, rxp, and hit_data_info will all be updated if we hit absorptive
       boundary.
*************************************************************************/
struct wall *ray_trace_2D(
    struct volume *world,
    struct surface_molecule *sm,
    struct vector2 *disp,
    struct vector2 *pos,
    int *kill_me,
    struct rxn **rxp,
    struct hit_data **hit_data_info) {

  struct hit_data *hit_data_head = NULL;

  struct wall *this_wall = sm->grid->surface;

  struct vector2 orig_pos = { .u = sm->s_pos.u,
                               .v = sm->s_pos.v
                            };
  struct vector2 this_pos = { .u = sm->s_pos.u,
                              .v = sm->s_pos.v
                            };
  struct vector2 this_disp = { .u = disp->u,
                               .v = disp->v
                             };
  struct periodic_image orig_box = { .x = sm->periodic_box->x,
                                     .y = sm->periodic_box->y,
                                     .z = sm->periodic_box->z
                                   };
  struct vector3 origin_xyz;
  uv2xyz(&this_pos, this_wall, &origin_xyz);

  struct rxn *rx = NULL;
  /* Will break out with return or break when we're done traversing walls */
  while (1) {

    int this_wall_edge_region_border = 0;
    int absorb_now = 0;
    int reflect_now = 0;

    /* Index of the wall edge that the SM hits */
    struct vector2 boundary_pos;
    int index_edge_was_hit =
        find_edge_point(this_wall, &this_pos, &this_disp, &boundary_pos);

    if (world->periodic_box_obj) {
      struct vector3 *hit_xyz = reflect_periodic_2D(
          world,
          index_edge_was_hit,
          &this_pos,
          this_wall,
          &this_disp,
          &boundary_pos,
          &origin_xyz);

      // We hit the periodic box! Update PBC, "reflect", and keep moving.
      if (hit_xyz) {
        struct vector3 teleport_xyz = { .x = hit_xyz->x,
                                        .y = hit_xyz->y,
                                        .z = hit_xyz->z
                                      };
        change_boxes_2D(
          world->periodic_traditional, sm, world->periodic_box_obj, hit_xyz,
          &teleport_xyz);
        if (world->periodic_traditional) {
          // Some of these uv coords might not be valid (i.e. they fall off
          // edge of triangle), but (i think) that's okay. we're ultimately
          // just trying to get remaining uv displacement.
          struct vector2 target_uv = { .u = this_pos.u + this_disp.u,
                                       .v = this_pos.v + this_disp.v
                                     };
          struct vector3 target_xyz;
          uv2xyz(&target_uv, this_wall, &target_xyz);
          struct vector3 remaining_disp_xyz = { .x = target_xyz.x - hit_xyz->x,
                                                .y = target_xyz.y - hit_xyz->y,
                                                .z = target_xyz.z - hit_xyz->z
                                              };
          struct vector3 new_target_xyz = { .x = teleport_xyz.x + remaining_disp_xyz.x,
                                            .y = teleport_xyz.y + remaining_disp_xyz.y,
                                            .z = teleport_xyz.z + remaining_disp_xyz.z,
                                          };
          int grid_index = 0;
          int *grid_index_p = &grid_index;
          struct wall *prev_wall = this_wall;
          // this_pos is also being updated here.
          this_wall = find_closest_wall(
            world, &teleport_xyz, 0.0, &this_pos, grid_index_p, sm->properties, NULL, NULL, NULL);
          // Try again if we can't find a place
          if ((this_wall == NULL) ||
              (this_wall->parent_object != prev_wall->parent_object) ) {
            *hit_data_info = hit_data_head;
            free(hit_xyz);
            return NULL;
          }
          struct vector2 new_target_uv;
          xyz2uv(&new_target_xyz, this_wall, &new_target_uv);
          this_disp.u = new_target_uv.u - this_pos.u;
          this_disp.v = new_target_uv.v - this_pos.v;
        }
        else {
          origin_xyz.x = hit_xyz->x;
          origin_xyz.y = hit_xyz->y;
          origin_xyz.z = hit_xyz->z;
        }
        free(hit_xyz);
        continue;
      }
    }

    // Ambiguous edge collision. Give up and try again from diffuse_2D.
    if (index_edge_was_hit == -2) {
      sm->s_pos.u = orig_pos.u;
      sm->s_pos.v = orig_pos.v;
      sm->periodic_box->x = orig_box.x;
      sm->periodic_box->y = orig_box.y;
      sm->periodic_box->z = orig_box.z;
      *hit_data_info = hit_data_head;
      return NULL;
    }
    // We didn't hit the edge. Stay inside this wall. We're done!
    else if (index_edge_was_hit == -1) {
      pos->u = this_pos.u + this_disp.u;
      pos->v = this_pos.v + this_disp.v;

      sm->s_pos.u = orig_pos.u;
      sm->s_pos.v = orig_pos.v;
      *hit_data_info = hit_data_head;
      return this_wall;
    }
    // Not ambiguous (-2) or inside wall (-1), must have hit edge (0, 1, 2)

    struct vector2 old_pos = {.u = this_pos.u,
                              .v = this_pos.v
                             };
    /* We hit the edge - check for the reflection/absorption from the
       edges of the wall if they are region borders
       Note - here we test for potential collisions with the region
       border while moving INSIDE OUT */
    struct rxn *matching_rxns[MAX_MATCHING_RXNS];
    if (sm->properties->flags & CAN_REGION_BORDER) {
      reflect_absorb_inside_out(
          world, sm, hit_data_head, &rx, matching_rxns, boundary_pos, this_wall,
          index_edge_was_hit, &reflect_now, &absorb_now,
          &this_wall_edge_region_border);
      if (absorb_now) {
        *kill_me = 1;
        *rxp = rx;
        *hit_data_info = hit_data_head;
        return NULL;
      }
    }

    /* no reflection - keep going */
    struct vector2 new_disp;
    if (!reflect_now) {
      struct wall *target_wall =
          traverse_surface(this_wall, &old_pos, index_edge_was_hit, &this_pos);

      if (target_wall != NULL) {
        if (sm->properties->flags & CAN_REGION_BORDER) {

          /* We hit the edge - check for the reflection/absorption from the
             edges of the wall if they are region borders
             Note - here we test for potential collisions with the region
             border while moving OUTSIDE IN */

          if (reflect_absorb_outside_in(
              world, sm, &hit_data_head, &rx, matching_rxns, boundary_pos,
              target_wall, this_wall, &reflect_now, &absorb_now,
              this_wall_edge_region_border)) {
            if (absorb_now) {
              *kill_me = 1;
              *rxp = rx;
              *hit_data_info = hit_data_head;
              return NULL;
            }
          }
        }

        if (!reflect_now) {
          this_disp.u = old_pos.u + this_disp.u;
          this_disp.v = old_pos.v + this_disp.v;
          traverse_surface(this_wall, &this_disp, index_edge_was_hit, &new_disp);
          this_disp.u = new_disp.u - this_pos.u;
          this_disp.v = new_disp.v - this_pos.v;
          this_wall = target_wall;

          continue;
        }
      }
    }

    if (!reflect_now) {
      *hit_data_info = hit_data_head;
    }

  /* If we reach this point, assume we reflect off the edge since there is no
   * neighboring wall
   *
   * NOTE: this_pos has been corrupted by traverse_surface; use old_pos to find
   * out whether the present wall edge is a region border
   */
    new_disp.u = this_disp.u - (boundary_pos.u - old_pos.u);
    new_disp.v = this_disp.v - (boundary_pos.v - old_pos.v);

    double f;
    struct vector2 reflector;

    switch (index_edge_was_hit) {
    case 0:
      new_disp.v *= -1.0;
      break;
    case 1:
      reflector.u = -this_wall->uv_vert2.v;
      reflector.v = this_wall->uv_vert2.u - this_wall->uv_vert1_u;
      f = 1.0 / sqrt(reflector.u * reflector.u + reflector.v * reflector.v);
      reflector.u *= f;
      reflector.v *= f;
      f = 2.0 * (new_disp.u * reflector.u + new_disp.v * reflector.v);
      new_disp.u -= f * reflector.u;
      new_disp.v -= f * reflector.v;
      break;
    case 2:
      reflector.u = this_wall->uv_vert2.v;
      reflector.v = -this_wall->uv_vert2.u;
      f = 1.0 / sqrt(reflector.u * reflector.u + reflector.v * reflector.v);
      reflector.u *= f;
      reflector.v *= f;
      f = 2.0 * (new_disp.u * reflector.u + new_disp.v * reflector.v);
      new_disp.u -= f * reflector.u;
      new_disp.v -= f * reflector.v;
      break;

    default:
      UNHANDLED_CASE(index_edge_was_hit);
    }

    this_pos.u = boundary_pos.u;
    this_pos.v = boundary_pos.v;
    this_disp.u = new_disp.u;
    this_disp.v = new_disp.v;

  } /* end while(1) */

  sm->s_pos.u = orig_pos.u;
  sm->s_pos.v = orig_pos.v;

  *hit_data_info = hit_data_head;

  return NULL;
}

/*************************************************************************
ray_trace:
  In: world: simulation state
      init_pos: position of molecule that is moving
      c: linked list of potential collisions with molecules (we could react)
      sv: subvolume that we start in
      v: displacement vector from current to new location
      reflectee: wall we have reflected off of and should not hit again
  Out: collision list of walls and molecules we intersected along our ray
       (current subvolume only), plus the subvolume wall.  Will always
       return at least the subvolume wall--NULL indicates an out of
       memory error.
*************************************************************************/
struct collision *ray_trace(struct volume *world, struct vector3 *init_pos,
                            struct collision *c, struct subvolume *sv,
                            struct vector3 *v, struct wall *reflectee) {
  /* time, in units of of the molecule's time step, at which molecule
     will cross the x,y,z partitions, respectively. */
  double tx, ty, tz;

  world->ray_voxel_tests++;

  struct collision *shead = NULL;
  struct collision *smash = (struct collision *)CHECKED_MEM_GET(
      sv->local_storage->coll, "collision structure");

  struct wall_list fake_wlp;
  fake_wlp.next = sv->wall_head;

  // Check wall collisions
  for (struct wall_list *wlp = sv->wall_head; wlp != NULL; wlp = wlp->next) {
    if (wlp->this_wall == reflectee)
      continue;

    int i = collide_wall(init_pos, v, wlp->this_wall, &(smash->t), &(smash->loc),
                     1, world->rng, world->notify, &(world->ray_polygon_tests));
    if (i == COLLIDE_REDO) {
      if (shead != NULL)
        mem_put_list(sv->local_storage->coll, shead);
      shead = NULL;
      wlp = &fake_wlp;
      continue;
    } else if (i != COLLIDE_MISS) {
      world->ray_polygon_colls++;

      smash->what = COLLIDE_WALL + i;
      smash->target = (void *)wlp->this_wall;
      smash->next = shead;
      shead = smash;
      smash = (struct collision *)CHECKED_MEM_GET(sv->local_storage->coll,
                                                  "collision structure");
    }
  }

  double dx, dy, dz;
  dx = dy = dz = 0.0;
  int i = -10;
  if (v->x < 0.0) {
    dx = world->x_fineparts[sv->llf.x] - init_pos->x;
    i = 0;
  } else if (v->x > 0.0) {
    dx = world->x_fineparts[sv->urb.x] - init_pos->x;
    i = 1;
  }

  int j = -10;
  if (v->y < 0.0) {
    dy = world->y_fineparts[sv->llf.y] - init_pos->y;
    j = 0;
  } else if (v->y > 0.0) {
    dy = world->y_fineparts[sv->urb.y] - init_pos->y;
    j = 1;
  }

  int k = -10;
  if (v->z < 0.0) {
    dz = world->z_fineparts[sv->llf.z] - init_pos->z;
    k = 0;
  } else if (v->z > 0.0) {
    dz = world->z_fineparts[sv->urb.z] - init_pos->z;
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

  smash->loc.x = init_pos->x + smash->t * v->x;
  smash->loc.y = init_pos->y + smash->t * v->y;
  smash->loc.z = init_pos->z + smash->t * v->z;

  smash->target = sv;
  smash->next = shead;
  shead = smash;

  // Check molecule collisions
  for (; c != NULL; c = c->next) {
    struct abstract_molecule *a = (struct abstract_molecule *)c->target;
    if (a->properties == NULL)
      continue;

    i = collide_mol(init_pos, v, a, &(c->t), &(c->loc), world->rx_radius_3d);
    if (i != COLLIDE_MISS) {
      smash = (struct collision *)CHECKED_MEM_GET(sv->local_storage->coll,
                                                  "collision structure");
      memcpy(smash, c, sizeof(struct collision));

      smash->what = COLLIDE_VOL + i;

      smash->next = shead;
      shead = smash;
    }
  }

  return shead;
}

/******************************/
/** exact_disk stuff follows **/
/******************************/

/* This is the same as the normal vector3, but the coordinates have different
 * names */
struct exd_vector3 {
  double m, u, v;
};

/****************************************************************
exd_zetize:
In: y coordinate (as in atan2)
    x coordinate (as in atan2)
Out: Zeta value corresponding to (y,x), in the range 0 to 4.
     Zeta is a substitute for the angle theta, and this function
     is a substitute for atan2(y,x) which returns theta. Like
     theta, zeta increases throughout the unit circle, but it
     only has 8-fold symmetry instead of perfect symmetry.  Zeta
     values 0-1 are the first quadrant, 1-2 the second, and so
     on.  Zeta is a monotonically increasing function of theta,
     but requires ~9x less computation time--valuable for when
     you need to sort by angle but don't need the angle itself.
Note: This is a utility function in 'exact_disk()'.
****************************************************************/
/* Speed: 9ns (compare with 84ns for atan2) */
/* Added extra computations--speed not retested yet */
static double exd_zetize(double y, double x) {
  if (y >= 0.0) {
    if (x >= 0) {
      if (x < y)
        return 1.0 - 0.5 * x / y;
      else
        return 0.5 * y / x;
    } else {
      if (-x < y)
        return 1.0 - 0.5 * x / y;
      else
        return 2.0 + 0.5 * y / x;
    }
  } else {
    if (x <= 0) {
      if (y < x)
        return 3.0 - 0.5 * x / y;
      else
        return 2.0 + 0.5 * y / x;
    } else {
      if (x < -y)
        return 3.0 - 0.5 * x / y;
      else
        return 4.0 + 0.5 * y / x;
    }
  }
}

/*********************************************************************
exd_coordize:
In: movement vector
    place to store unit movement vector (first basis vector)
    place to store second basis vector
    place to store third basis vector
Out: No return value.  Unit vectors m,u,v are set such that vector m
     is in the direction of vector mv, and vectors u and v are
     orthogonal to m.  The vectors m,u,v, form a right-handed
     coordinate system.
Note: This is a utility function for 'exact_disk()'.
*********************************************************************/
/* Speed: 86ns on azzuri (as marked + 6ns function call overhead) */
static void exd_coordize(struct vector3 *mv, struct vector3 *m,
                         struct vector3 *u, struct vector3 *v) {
  double a;

  /* Normalize input vector -- 27ns */
  a = 1.0 / sqrt(mv->x * mv->x + mv->y * mv->y + mv->z * mv->z);
  m->x = a * mv->x;
  m->y = a * mv->y;
  m->z = a * mv->z;

  /* Find orthogonal vectors -- 21ns */
  if (m->x * m->x > m->y * m->y) {
    if (m->x * m->x > m->z * m->z) {
      if (m->y * m->y > m->z * m->z) {
        u->x = m->y;
        u->y = -m->x;
        u->z = 0.0;
        a = 1.0 - m->z * m->z;
        v->x = m->z * m->x;
        v->y = m->z * m->y;
        v->z = -a;
      } else {
        u->x = m->z;
        u->y = 0.0;
        u->z = -m->x;
        a = 1.0 - m->y * m->y;
        v->x = -m->y * m->x;
        v->y = a;
        v->z = -m->y * m->z;
      }
    } else {
      u->x = -m->z;
      u->y = 0.0;
      u->z = m->x;
      a = 1.0 - m->y * m->y;
      v->x = m->y * m->x;
      v->y = -a;
      v->z = m->y * m->z;
    }
  } else {
    if (m->y * m->y > m->z * m->z) {
      if (m->x * m->x > m->z * m->z) {
        u->x = -m->y;
        u->y = m->x;
        u->z = 0.0;
        a = 1.0 - m->z * m->z;
        v->x = -m->z * m->x;
        v->y = -m->z * m->y;
        v->z = a;
      } else {
        u->x = 0.0;
        u->y = m->z;
        u->z = -m->y;
        a = 1.0 - m->x * m->x;
        v->x = -a;
        v->y = m->x * m->y;
        v->z = m->x * m->z;
      }
    } else {
      u->x = 0.0;
      u->y = -m->z;
      u->z = m->y;
      a = 1.0 - m->x * m->x;
      v->x = a;
      v->y = -m->x * m->y;
      v->z = -m->x * m->z;
    }
  }

  /* Normalize orthogonal vectors -- 32ns */
  a = 1 / sqrt(a);
  u->x *= a;
  u->y *= a;
  u->z *= a;
  v->x *= a;
  v->y *= a;
  v->z *= a;
}

/* Exact Disk Flags */
/* Flags for the exact disk computation */
enum {
  EXD_HEAD,
  EXD_TAIL,
  EXD_CROSS,
  EXD_SPAN,
  EXD_OTHER
};

/* Negative numbers used as flags for reaction disks */
/* Note: TARGET_OCCLUDED is assumed for any negative number not defined here */
#define TARGET_OCCLUDED -1

/*************************************************************************
exact_disk:
  In: world: simulation state
      loc: location of moving molecule at time of collision
      mv: movement vector for moving molecule
      R: interaction radius
      sv:  subvolume the moving molecule is in
      moving: the moving molecule
      target: the target molecule at time of collision
      use_expanded_list:
      x_fineparts:
      y_fineparts:
      z_fineparts:
  Out: The fraction of a full interaction disk that is actually
       accessible to the moving molecule, computed exactly from the
       geometry, or TARGET_OCCLUDED if the path to the target molecule is
       blocked.
*************************************************************************/
double exact_disk(struct volume *world, struct vector3 *loc, struct vector3 *mv,
                  double R, struct subvolume *sv,
                  struct volume_molecule *moving,
                  struct volume_molecule *target, int use_expanded_list,
                  double *x_fineparts, double *y_fineparts,
                  double *z_fineparts) {
#define EXD_SPAN_CALC(v1, v2, p)                                               \
  ((v1)->u - (p)->u) * ((v2)->v - (p)->v) -                                    \
      ((v2)->u - (p)->u) * ((v1)->v - (p)->v)
#define EXD_TIME_CALC(v1, v2, p)                                               \
  ((p)->u *(v1)->v - (p)->v *(v1)->u) /                                        \
      ((p)->v *((v2)->u - (v1)->u) - (p)->u *((v2)->v - (v1)->v))
  struct wall_list *wl;
  struct wall *w;
  struct vector3 llf, urb;

  struct exd_vector3 v0muv, v1muv, v2muv;
  struct exd_vertex pa, pb;
  struct exd_vertex *ppa, *ppb, *pqa, *pqb, *vertex_head, *vp, *vq, *vr, *vs;
  double pa_pb;
  int n_verts, n_edges;
  int p_flags;

  double R2;
  struct vector3 m, u, v;
  struct exd_vector3 Lmuv;
  struct exd_vertex sm;
  double m2_i;
  double l_n, m_n;
  double a, b, c, d, r, s, t, A, zeta, last_zeta;
  int i;
  int num_matching_rxns = 0;
  struct rxn *matching_rxns[MAX_MATCHING_RXNS];

  /* Initialize */
  vertex_head = NULL;
  n_verts = 0;
  n_edges = 0;

  /* Partially set up coordinate systems for first pass */
  R2 = R * R;
  m2_i = 1.0 / (mv->x * mv->x + mv->y * mv->y + mv->z * mv->z);
  Lmuv.m = Lmuv.u = Lmuv.v = 0.0;      /* Keep compiler happy */
  sm.u = sm.v = sm.r2 = sm.zeta = 0.0; /* More compiler happiness */

  /* Set up coordinate system and convert vertices */
  exd_coordize(mv, &m, &u, &v);

  Lmuv.m = loc->x * m.x + loc->y * m.y + loc->z * m.z;
  Lmuv.u = loc->x * u.x + loc->y * u.y + loc->z * u.z;
  Lmuv.v = loc->x * v.x + loc->y * v.y + loc->z * v.z;

  if (!distinguishable_vec3(loc, &(target->pos), EPS_C)) { /* Hit target exactly! */
    sm.u = sm.v = sm.r2 = sm.zeta = 0.0;
  } else { /* Find location of target in moving-molecule-centric coords */
    sm.u = (target->pos.x - loc->x) * u.x + (target->pos.y - loc->y) * u.y +
           (target->pos.z - loc->z) * u.z;
    sm.v = (target->pos.x - loc->x) * v.x + (target->pos.y - loc->y) * v.y +
           (target->pos.z - loc->z) * v.z;
    sm.r2 = sm.u * sm.u + sm.v * sm.v;
    sm.zeta = exd_zetize(sm.v, sm.u);
  }

  /* Find walls that occlude the interaction disk (or block the reaction) */
  for (wl = sv->wall_head; wl != NULL; wl = wl->next) {
    w = wl->this_wall;

    /* Ignore this wall if it is too far away! */

    /* Find distance from plane of wall to molecule */
    l_n = loc->x * w->normal.x + loc->y * w->normal.y + loc->z * w->normal.z;
    d = w->d - l_n;

    /* See if we're within interaction distance of wall */
    m_n = mv->x * w->normal.x + mv->y * w->normal.y + mv->z * w->normal.z;

    if (d * d >= R2 * (1 - m2_i * m_n * m_n))
      continue;

    /* Ignore this wall if no overlap between wall & disk bounding boxes */

    /* Find wall bounding box */
    urb.x = llf.x = w->vert[0]->x;
    if (w->vert[1]->x < llf.x)
      llf.x = w->vert[1]->x;
    else
      urb.x = w->vert[1]->x;
    if (w->vert[2]->x < llf.x)
      llf.x = w->vert[2]->x;
    else if (w->vert[2]->x > urb.x)
      urb.x = w->vert[2]->x;

    urb.y = llf.y = w->vert[0]->y;
    if (w->vert[1]->y < llf.y)
      llf.y = w->vert[1]->y;
    else
      urb.y = w->vert[1]->y;
    if (w->vert[2]->y < llf.y)
      llf.y = w->vert[2]->y;
    else if (w->vert[2]->y > urb.y)
      urb.y = w->vert[2]->y;

    urb.z = llf.z = w->vert[0]->z;
    if (w->vert[1]->z < llf.z)
      llf.z = w->vert[1]->z;
    else
      urb.z = w->vert[1]->z;
    if (w->vert[2]->z < llf.z)
      llf.z = w->vert[2]->z;
    else if (w->vert[2]->z > urb.z)
      urb.z = w->vert[2]->z;

    /* Reject those without overlapping bounding boxes */
    b = R2 * (1.0 - mv->x * mv->x * m2_i);
    a = llf.x - loc->x;
    if (a > 0 && a * a >= b)
      continue;
    a = loc->x - urb.x;
    if (a > 0 && a * a >= b)
      continue;

    b = R2 * (1.0 - mv->y * mv->y * m2_i);
    a = llf.y - loc->y;
    if (a > 0 && a * a >= b)
      continue;
    a = loc->y - urb.y;
    if (a > 0 && a * a >= b)
      continue;

    b = R2 * (1.0 - mv->z * mv->z * m2_i);
    a = llf.z - loc->z;
    if (a > 0 && a * a >= b)
      continue;
    a = loc->z - urb.z;
    if (a > 0 && a * a >= b)
      continue;

    /* Ignore this wall if moving molecule can travel through it */

    /* Reject those that the moving particle can travel through */
    if ((moving->properties->flags & CAN_VOLWALL) != 0) {
      num_matching_rxns = trigger_intersect(
          world->reaction_hash, world->rx_hashsize, world->all_mols,
          world->all_volume_mols, world->all_surface_mols,
          moving->properties->hashval, (struct abstract_molecule *)moving, 0, w,
          matching_rxns, 1, 1, 0);
      if (num_matching_rxns == 0)
        continue;
      int blocked = 0;
      for (i = 0; i < num_matching_rxns; i++) {
        if (matching_rxns[i]->n_pathways == RX_REFLEC) {
          blocked = 1;
        }
      }
      if (!blocked) {
        continue;
      }
    }

    /* Find line of intersection between wall and disk */

#if 0
    /* Set up coordinate system and convert vertices */
    if (uncoordinated) {
      exd_coordize(mv, &m, &u, &v);

      Lmuv.m = loc->x * m.x + loc->y * m.y + loc->z * m.z;
      Lmuv.u = loc->x * u.x + loc->y * u.y + loc->z * u.z;
      Lmuv.v = loc->x * v.x + loc->y * v.y + loc->z * v.z;

      if (!distinguishable_vec3(loc, &(target->pos),
                                EPS_C)) /* Hit target exactly! */
      {
        sm.u = sm.v = sm.r2 = sm.zeta = 0.0;
      } else /* Find location of target in moving-molecule-centric coords */
      {
        sm.u = (target->pos.x - loc->x) * u.x + (target->pos.y - loc->y) * u.y +
               (target->pos.z - loc->z) * u.z;
        sm.v = (target->pos.x - loc->x) * v.x + (target->pos.y - loc->y) * v.y +
               (target->pos.z - loc->z) * v.z;
        sm.r2 = sm.u * sm.u + sm.v * sm.v;
        sm.zeta = exd_zetize(sm.v, sm.u);
      }
      uncoordinated = 0;
    }
#endif
    v0muv.m = w->vert[0]->x * m.x + w->vert[0]->y * m.y + w->vert[0]->z * m.z -
              Lmuv.m;
    v0muv.u = w->vert[0]->x * u.x + w->vert[0]->y * u.y + w->vert[0]->z * u.z -
              Lmuv.u;
    v0muv.v = w->vert[0]->x * v.x + w->vert[0]->y * v.y + w->vert[0]->z * v.z -
              Lmuv.v;

    v1muv.m = w->vert[1]->x * m.x + w->vert[1]->y * m.y + w->vert[1]->z * m.z -
              Lmuv.m;
    v1muv.u = w->vert[1]->x * u.x + w->vert[1]->y * u.y + w->vert[1]->z * u.z -
              Lmuv.u;
    v1muv.v = w->vert[1]->x * v.x + w->vert[1]->y * v.y + w->vert[1]->z * v.z -
              Lmuv.v;

    v2muv.m = w->vert[2]->x * m.x + w->vert[2]->y * m.y + w->vert[2]->z * m.z -
              Lmuv.m;
    v2muv.u = w->vert[2]->x * u.x + w->vert[2]->y * u.y + w->vert[2]->z * u.z -
              Lmuv.u;
    v2muv.v = w->vert[2]->x * v.x + w->vert[2]->y * v.y + w->vert[2]->z * v.z -
              Lmuv.v;

    /* Draw lines between points and pick intersections with plane of m=0 */
    if ((v0muv.m < 0) == (v1muv.m < 0)) /* v0,v1 on same side */
    {
      if ((v2muv.m < 0) == (v1muv.m < 0))
        continue;
      t = v0muv.m / (v0muv.m - v2muv.m);
      pa.u = v0muv.u + (v2muv.u - v0muv.u) * t;
      pa.v = v0muv.v + (v2muv.v - v0muv.v) * t;
      t = v1muv.m / (v1muv.m - v2muv.m);
      pb.u = v1muv.u + (v2muv.u - v1muv.u) * t;
      pb.v = v1muv.v + (v2muv.v - v1muv.v) * t;
    } else if ((v0muv.m < 0) == (v2muv.m < 0)) /* v0,v2 on same side */
    {
      t = v0muv.m / (v0muv.m - v1muv.m);
      pa.u = v0muv.u + (v1muv.u - v0muv.u) * t;
      pa.v = v0muv.v + (v1muv.v - v0muv.v) * t;
      t = v2muv.m / (v2muv.m - v1muv.m);
      pb.u = v2muv.u + (v1muv.u - v2muv.u) * t;
      pb.v = v2muv.v + (v1muv.v - v2muv.v) * t;
    } else /* v1, v2 on same side */
    {
      t = v1muv.m / (v1muv.m - v0muv.m);
      pa.u = v1muv.u + (v0muv.u - v1muv.u) * t;
      pa.v = v1muv.v + (v0muv.v - v1muv.v) * t;
      t = v2muv.m / (v2muv.m - v0muv.m);
      pb.u = v2muv.u + (v0muv.u - v2muv.u) * t;
      pb.v = v2muv.v + (v0muv.v - v2muv.v) * t;
    }

    /* Check to make sure endpoints are sensible */
    pa.r2 = pa.u * pa.u + pa.v * pa.v;
    pb.r2 = pb.u * pb.u + pb.v * pb.v;
    if (pa.r2 < EPS_C * R2 ||
        pb.r2 < EPS_C * R2) /* Can't tell where origin is relative to wall
                               endpoints */
    {
      if (vertex_head != NULL)
        mem_put_list(sv->local_storage->exdv, vertex_head);
      return TARGET_OCCLUDED;
    }
    if (!distinguishable(pa.u * pb.v, pb.u * pa.v, EPS_C) &&
        pa.u * pb.u + pa.v * pb.v <
            0) /* Antiparallel, can't tell which side of wall origin is on */
    {
      if (vertex_head != NULL)
        mem_put_list(sv->local_storage->exdv, vertex_head);
      return TARGET_OCCLUDED;
    }

    /* Intersect line with circle; skip this wall if no intersection */
    t = 0;
    s = 1;
    if (pa.r2 > R2 || pb.r2 > R2) {
      pa_pb = pa.u * pb.u + pa.v * pb.v;
      if (!distinguishable(pa.r2 + pb.r2, 2 * pa_pb,
                           EPS_C)) /* Wall endpoints are basically on top of
                                      each other */
      {
        /* Might this tiny bit of wall block the target?  If not, continue,
         * otherwise return TARGET_OCCLUDED */
        /* Safe if we're clearly closer; in danger if we're even remotely
         * parallel, otherwise surely safe */
        /* Note: use SQRT_EPS_C for cross products since previous test vs. EPS_C
         * was on squared values (linear difference term cancels) */
        if (sm.r2 < pa.r2 && sm.r2 < pb.r2 &&
            distinguishable(sm.r2, pa.r2, EPS_C) &&
            distinguishable(sm.r2, pa.r2, EPS_C))
          continue;
        if (!distinguishable(sm.u * pa.v, sm.v * pa.u, SQRT_EPS_C) ||
            !distinguishable(sm.u * pb.v, sm.v * pb.u, SQRT_EPS_C)) {
          if (vertex_head != NULL)
            mem_put_list(sv->local_storage->exdv, vertex_head);
          return TARGET_OCCLUDED;
        }
        continue;
      }
      a = 1.0 / (pa.r2 + pb.r2 - 2 * pa_pb);
      b = (pa_pb - pa.r2) * a;
      c = (R2 - pa.r2) * a;
      d = b * b + c;
      if (d <= 0)
        continue;
      d = sqrt(d);
      t = -b - d;
      if (t >= 1)
        continue;
      if (t < 0)
        t = 0;
      s = -b + d;
      if (s <= 0)
        continue;
      if (s > 1)
        s = 1;
    }

    /* Add this edge to the growing list, or return -1 if edge blocks target */

    /* Construct final endpoints and prepare to store them */
    ppa = (struct exd_vertex *)CHECKED_MEM_GET(sv->local_storage->exdv,
                                               "exact disk vertex");
    ppb = (struct exd_vertex *)CHECKED_MEM_GET(sv->local_storage->exdv,
                                               "exact disk vertex");
    if (t > 0) {
      ppa->u = pa.u + t * (pb.u - pa.u);
      ppa->v = pa.v + t * (pb.v - pa.v);
      ppa->r2 = ppa->u * ppa->u + ppa->v * ppa->v;
      ppa->zeta = exd_zetize(ppa->v, ppa->u);
    } else {
      ppa->u = pa.u;
      ppa->v = pa.v;
      ppa->r2 = pa.r2;
      ppa->zeta = exd_zetize(pa.v, pa.u);
    }
    if (s < 1) {
      ppb->u = pa.u + s * (pb.u - pa.u);
      ppb->v = pa.v + s * (pb.v - pa.v);
      ppb->r2 = ppb->u * ppb->u + ppb->v * ppb->v;
      ppb->zeta = exd_zetize(ppb->v, ppb->u);
    } else {
      ppb->u = pb.u;
      ppb->v = pb.v;
      ppb->r2 = pb.r2;
      ppb->zeta = exd_zetize(pb.v, pb.u);
    }

    /* It's convenient if ppa is earlier, ccw, than ppb */
    a = (ppb->zeta - ppa->zeta);
    if (a < 0)
      a += 4.0;
    if (a >= 2.0) {
      vp = ppb;
      ppb = ppa;
      ppa = vp;
      a = 4.0 - a;
    }

    /* Detect a blocked reaction: line is between origin and target */
    b = (sm.zeta - ppa->zeta);
    if (b < 0)
      b += 4.0;

    if (b < a) {
      c = (ppa->u - sm.u) * (ppb->v - sm.v) - (ppa->v - sm.v) * (ppb->u - sm.u);
      if (c < 0 || !distinguishable((ppa->u - sm.u) * (ppb->v - sm.v),
                                    (ppa->v - sm.v) * (ppb->u - sm.u),
                                    EPS_C)) /* Blocked! */
      {
        ppa->next = ppb;
        ppb->next = vertex_head;
        mem_put_list(sv->local_storage->exdv, ppa);
        return TARGET_OCCLUDED;
      }
    }

    ppa->role = EXD_HEAD;
    ppb->role = EXD_TAIL;
    ppa->e = ppb;
    ppb->e = NULL;

    ppb->next = vertex_head;
    ppa->next = ppb;
    vertex_head = ppa;
    n_verts += 2;
    n_edges++;
  }

  /* Find partition boundaries that occlude the interaction disk */
  if (!use_expanded_list) /* We'll hit partitions */
  {
    /* First see if any overlap */
    p_flags = 0;

    d = loc->x - x_fineparts[sv->llf.x];
    if (d < R) {
      c = R2 * (mv->y * mv->y + mv->z * mv->z) * m2_i;
      if (d * d < c)
        p_flags |= X_NEG_BIT;
      d = x_fineparts[sv->urb.x] - loc->x;
      if (d * d < c)
        p_flags |= X_POS_BIT;
    } else {
      d = x_fineparts[sv->urb.x] - loc->x;
      if (d < R && d * d < R2 * (mv->y * mv->y + mv->z * mv->z) * m2_i)
        p_flags |= X_POS_BIT;
    }

    d = loc->y - y_fineparts[sv->llf.y];
    if (d < R) {
      c = R2 * (mv->x * mv->x + mv->z * mv->z) * m2_i;
      if (d * d < c)
        p_flags |= Y_NEG_BIT;
      d = y_fineparts[sv->urb.y] - loc->y;
      if (d * d < c)
        p_flags |= Y_POS_BIT;
    } else {
      d = y_fineparts[sv->urb.y] - loc->y;
      if (d < R && d * d < R2 * (mv->x * mv->x + mv->z * mv->z) * m2_i)
        p_flags |= Y_POS_BIT;
    }

    d = loc->z - z_fineparts[sv->llf.z];
    if (d < R) {
      c = R2 * (mv->y * mv->y + mv->x * mv->x) * m2_i;
      if (d * d < c)
        p_flags |= Z_NEG_BIT;
      d = z_fineparts[sv->urb.z] - loc->z;
      if (d * d < c)
        p_flags |= Z_POS_BIT;
    } else {
      d = z_fineparts[sv->urb.z] - loc->z;
      if (d < R && d * d < R2 * (mv->y * mv->y + mv->x * mv->x) * m2_i)
        p_flags |= Z_POS_BIT;
    }

    /* Now find the lines created by any that do overlap */
    if (p_flags) {
      for (i = 1; i <= p_flags; i *= 2) {
        if ((i & p_flags) != 0) {
          /* Load up the relevant variables */
          switch (i) {
          case X_NEG_BIT:
            d = x_fineparts[sv->llf.x] - loc->x;
            a = u.x;
            b = v.x;
            break;
          case X_POS_BIT:
            d = x_fineparts[sv->urb.x] - loc->x;
            a = u.x;
            b = v.x;
            break;
          case Y_NEG_BIT:
            d = y_fineparts[sv->llf.y] - loc->y;
            a = u.y;
            b = v.y;
            break;
          case Y_POS_BIT:
            d = y_fineparts[sv->urb.y] - loc->y;
            a = u.y;
            b = v.y;
            break;
          case Z_NEG_BIT:
            d = z_fineparts[sv->llf.z] - loc->z;
            a = u.z;
            b = v.z;
            break;
          case Z_POS_BIT:
            d = z_fineparts[sv->urb.z] - loc->z;
            a = u.z;
            b = v.z;
            break;
          default:
            continue;
          }

          if (!distinguishable(a, 0, EPS_C)) {
            s = d / b;
            if (s * s > R2) {
              mcell_internal_error(
                  "Unexpected results in exact disk: s=%.2f s^2=%.2f R2=%.2f\n",
                  s, s * s, R2);
              /*continue;*/
            }
            t = sqrt(R2 - s * s);
            pa.u = t;
            pa.v = s;
            pb.u = -t;
            pb.v = s;
          } else if (!distinguishable(b, 0, EPS_C)) {
            t = d / a;
            if (t * t > R2) {
              mcell_internal_error(
                  "Unexpected results in exact disk: t=%.2f t^2=%.2f R2=%.2f\n",
                  t, t * t, R2);
              /*continue;*/
            }
            s = sqrt(R2 - t * t);
            pa.u = t;
            pa.v = s;
            pb.u = t;
            pb.v = -s;
          } else {
            c = a * a + b * b;
            s = d * b;
            if (d * d > R2 * c) {
              mcell_internal_error("Unexpected results in exact disk: d=%.2f "
                                   "d^2=%.2f R2=%.2f c=%.2f R2*c=%.2f\n",
                                   d, d * d, R2, c, R2 * c);
              /*continue;*/
            }
            t = sqrt(R2 * c - d * d);
            c = 1.0 / c;
            r = 1.0 / a;
            pa.v = c * (s + t * a);
            pa.u = (d - b * pa.v) * r;
            pb.v = c * (s - t * a);
            pb.u = (d - b * pb.v) * r;
          }

          /* Create memory for the pair of vertices */
          ppa = (struct exd_vertex *)CHECKED_MEM_GET(sv->local_storage->exdv,
                                                     "exact disk vertex");
          ppb = (struct exd_vertex *)CHECKED_MEM_GET(sv->local_storage->exdv,
                                                     "exact disk vertex");

          a = exd_zetize(pa.v, pa.u);
          b = exd_zetize(pb.v, pb.u);
          c = b - a;
          if (c < 0)
            c += 4;
          if (c < 2) {
            ppa->u = pa.u;
            ppa->v = pa.v;
            ppa->r2 = pa.u * pa.u + pa.v * pa.v;
            ppa->zeta = a;
            ppb->u = pb.u;
            ppb->v = pb.v;
            ppb->r2 = pb.u * pb.u + pb.v * pb.v;
            ppb->zeta = b;
          } else {
            ppb->u = pa.u;
            ppb->v = pa.v;
            ppb->r2 = pa.u * pa.u + pa.v * pa.v;
            ppb->zeta = a;
            ppa->u = pb.u;
            ppa->v = pb.v;
            ppa->r2 = pb.u * pb.u + pb.v * pb.v;
            ppa->zeta = b;
          }

          ppa->role = EXD_HEAD;
          ppb->role = EXD_TAIL;
          ppa->e = ppb;
          ppb->e = NULL;

          ppb->next = vertex_head;
          ppa->next = ppb;
          vertex_head = ppa;
          n_verts += 2;
          n_edges++;
        }
      }
    }
  }

  /* Now that we have everything, see if we can perform simple calculations */

  /* Did we even find anything?  If not, return full area */
  if (n_edges == 0) {
    return 1.0;
  }
  /* If there is only one edge, just calculate it */
  else if (n_edges == 1) {
    ppa = vertex_head;
    ppb = ppa->e;

    a = ppa->u * ppb->u + ppa->v * ppb->v;
    b = ppa->u * ppb->v - ppa->v * ppb->u;
    if (a <= 0) /* Angle > pi/2 */
    {
      s = atan(-a / b) + 0.5 * MY_PI;
    } else {
      s = atan(b / a);
    }
    A = (0.5 * b + R2 * (MY_PI - 0.5 * s)) / (MY_PI * R2);

    mem_put_list(sv->local_storage->exdv, vertex_head);
    return A;
  }

  /* If there are multiple edges, calculating area is more complex. */

  /* Insertion sort the multiple edges */
  vp = vertex_head->next;
  ppa = ppb = vertex_head;
  ppa->next = NULL;
  ppa->span = NULL;
  while (vp != NULL) {
    /* Snip off one item from old list to add */
    vp->span = NULL;
    vq = vp->next;

    /* Add it to list with ppa as head and ppb as tail */
    if (vp->zeta < ppa->zeta) {
      vp->next = ppa;
      ppa = vp;
    } else {
      for (pqa = ppa; pqa->next != NULL; pqa = pqa->next) {
        if (vp->zeta < pqa->next->zeta)
          break;
      }
      vp->next = pqa->next;
      pqa->next = vp;
      if (vp->next == NULL)
        ppb = vp;
    }

    /* Repeat for remainder of old list */
    vp = vq;
  }

  /* Close circular list */
  vertex_head = ppa;
  ppb->next = ppa;

  /* Walk around the circle, inserting points where lines cross */
  ppb = NULL;
  for (ppa = vertex_head; ppa != vertex_head || ppb == NULL; ppa = ppa->next) {
    if (ppa->role != EXD_HEAD)
      continue;
    ppb = ppa->e;

    for (pqa = ppa->next; pqa != ppb; pqa = pqa->next) {
      if (pqa->role != EXD_HEAD)
        continue;
      pqb = pqa->e;

      /* Create displacement vectors */
      pa.u = ppb->u - ppa->u;
      pa.v = ppb->v - ppa->v;
      pb.u = pqb->u - pqa->u;
      pb.v = pqb->v - pqa->v;
      r = pb.u * pa.v - pa.u * pb.v;

      /* Check if lines are parallel--combine if so */
      if (r * r <
          EPS_C * (pa.u * pa.u + pa.v * pa.v) * (pb.u * pb.u + pb.v * pb.v)) {
        pqa->e = NULL;
        pqa->role = EXD_OTHER;

        a = pqb->zeta - ppb->zeta;
        if (a < 0)
          a += 4.0;

        if (a > 2) /* Other line is completely contained inside us */
        {
          pqb->role = EXD_OTHER;
        } else /* We have a new endpoint, so we need to check crosses again */
        {
          ppa->e = pqb;
          ppb->role = EXD_OTHER;
          ppb = pqb;
          pqa = ppa;
        }
        continue;
      }

      /* Check if these lines cross and find times at which they do */
      s = (ppa->u - pqa->u) * pa.v - (ppa->v - pqa->v) * pa.u;
      if (s * r <= EPS_C * R2 * R2)
        continue;
      t = s / r;
      if (t >= 1 - EPS_C)
        continue;
      if (pa.u * pa.u > pa.v * pa.v) {
        s = (pqa->u - ppa->u + t * pb.u) * pa.u;
        if (s <= EPS_C * R2 || s >= pa.u * pa.u * (1.0 - EPS_C))
          continue;
      } else {
        s = (pqa->v - ppa->v + t * pb.v) * pa.v;
        if (s <= EPS_C * R2 || s >= pa.v * pa.v * (1.0 - EPS_C))
          continue;
      }

      /* Create intersection point */
      vq = (struct exd_vertex *)CHECKED_MEM_GET(sv->local_storage->exdv,
                                                "exact disk vertex");
      vq->u = pqa->u + t * pb.u;
      vq->v = pqa->v + t * pb.v;
      vq->r2 = vq->u * vq->u + vq->v * vq->v;
      vq->zeta = exd_zetize(vq->v, vq->u);
      vq->e = ppb;
      vq->span = NULL;
      vq->role = EXD_CROSS;

      /* Insert new point into the list */
      for (vp = ppa; vp != ppb; vp = vp->next) {
        a = vq->zeta - vp->next->zeta;
        if (a > 2.0)
          a -= 4.0;
        else if (a < -2.0)
          a += 4.0;

        if (a < 0)
          break;
      }

      vq->next = vp->next;
      vp->next = vq;
      if (vq->zeta < vertex_head->zeta)
        vertex_head = vq;
    }
  }

  /* Collapse nearby points in zeta and R */
  for (vp = vertex_head, vq = NULL; vq != vertex_head; vp = vq) {
    for (vq = vp->next; vq != vertex_head; vq = vq->next) {
      if (vq->zeta - vp->zeta < EPS_C) {
        vq->zeta = vp->zeta;
        if (-EPS_C < vq->r2 - vp->r2 && EPS_C > vq->r2 - vp->r2) {
          vq->r2 = vp->r2;
          /* Mark crosses that occur multiple times--only need one */
          //          if (vq->role==EXD_CROSS && vp->role != EXD_OTHER) vq->role
          // = EXD_OTHER;
          //          else if (vp->role==EXD_CROSS && vq->role != EXD_OTHER)
          // vp->role = EXD_OTHER;
        }
      } else
        break;
    }
  }

  /* Register all spanning line segments */
  vq = NULL;
  for (vp = vertex_head; vp != vertex_head || vq == NULL; vp = vp->next) {
    if (vp->role != EXD_HEAD)
      continue;

    for (vq = vp->next; vq != vp->e; vq = vq->next) {
      if (!distinguishable(vq->zeta, vp->zeta, EPS_C))
        continue;
      if (!distinguishable(vq->zeta, vp->e->zeta, EPS_C))
        break;
      if (vq->role == EXD_OTHER)
        continue;

      vr = (struct exd_vertex *)CHECKED_MEM_GET(sv->local_storage->exdv,
                                                "exact disk vertex");
      vr->next = vq->span;
      vq->span = vr;
      vr->e = vp;
      vr->zeta = vq->zeta;
      vr->role = EXD_SPAN;
    }
  }

  /* Now we finally walk through and calculate the area */
  A = 0.0;
  zeta = 0.0;
  last_zeta = -1;
  vs = NULL;
  for (vp = vertex_head; zeta < 4.0 - EPS_C; vp = vp->next) {
    if (vp->role == EXD_OTHER)
      continue;
    if (!distinguishable(vp->zeta, last_zeta, EPS_C))
      continue;
    last_zeta = vp->zeta;

    /* Store data for the next tentatively approved point */
    if (vs == &pa)
      vr = &pb;
    else
      vr = &pa;
    vr->u = vp->u;
    vr->v = vp->v;
    vr->zeta = vp->zeta;
    if (vp->role == EXD_TAIL) {
      vr->r2 = R2 * (1.0 + EPS_C);
      vr->e = NULL;
    } else {
      vr->r2 = vp->r2;
      vr->e = vp->e;
    }

    /* Check head points at same place to see if they're closer */
    for (vq = vp->next; (!distinguishable(vq->zeta, last_zeta, EPS_C));
         vq = vq->next) {
      if (vq->role == EXD_HEAD) {
        if (vq->r2 < vp->r2 || vr->e == NULL) {
          vr->u = vq->u;
          vr->v = vq->v;
          vr->r2 = vq->r2;
          vr->e = vq->e;
        } else if (!distinguishable(vq->r2, vr->r2, EPS_C)) {
          b = EXD_SPAN_CALC(vr, vr->e, vq->e);
          if (b > 0)
            vr->e = vq->e;
        }
      }
    }

    /* Check each span to see if anything is closer than our approval point */
    for (vq = vp->span; vq != NULL; vq = vq->next) {
      ppa = vq->e;
      ppb = ppa->e;
      b = EXD_SPAN_CALC(ppa, ppb, vr);
      c = b * b;
      if (c < R2 * R2 * EPS_C) /* Span crosses the point */
      {
        if (vr->e == NULL) {
          vr->r2 = vr->u * vr->u + vr->v * vr->v;
          vr->e = ppb;
        } else {
          b = EXD_SPAN_CALC(vr, vr->e, ppb);
          if (b > 0)
            vr->e = ppb;
        }
      } else if (b < 0 ||
                 vr->e == NULL) /* Span is inside the point or spans tail */
      {
        t = EXD_TIME_CALC(ppa, ppb, vp);
        vr->u = ppa->u + t * (ppb->u - ppa->u);
        vr->v = ppa->v + t * (ppb->v - ppa->v);
        vr->r2 = vr->u * vr->u + vr->v * vr->v;
        vr->e = ppb;
      }
    }

    /* Should have an approved point in vr */
    if (vs == NULL) /* No angle traversed yet */
    {
      vs = vr;
    } else {
      c = vr->zeta - vs->zeta;
      if (c < 0)
        c += 4.0;
      if (/*vr->e != vs->e || c+zeta >= 4.0-EPS_C*/ c > EPS_C) {
        zeta += c;
        if (vs->e == NULL ||
            (vs->e->zeta - vs->zeta) * (vs->e->zeta - vs->zeta) <
                EPS_C * EPS_C) {
          if (c >= 2.0) /* More than pi */
          {
            vs->u = -vs->u;
            vs->v = -vs->v;
            A += 0.5 * MY_PI * R2;
          }
          a = vs->u * vr->u + vs->v * vr->v;
          b = vs->u * vr->v - vs->v * vr->u;
          if (a <= 0) /* More than pi/2 */
          {
            s = atan(-a / b) + 0.5 * MY_PI;
          } else {
            s = atan(b / a);
          }
          A += 0.5 * s * R2;
        } else {
          if (!distinguishable(vs->e->zeta, vr->zeta, EPS_C)) {
            A += 0.5 * (vs->u * vs->e->v - vs->v * vs->e->u);
          } else {
            t = EXD_TIME_CALC(vs, vs->e, vr);
            b = vs->u + (vs->e->u - vs->u) * t;
            c = vs->v + (vs->e->v - vs->v) * t;
            A += 0.5 * (vs->u * c - vs->v * b);
          }
        }
        vs = vr;
      } else {
        if (vr->e != NULL)
          vs = vr;
      }
    }
  }

  /* Finally, let's clean up the mess we made! */

  /* Deallocate lists */
  /* Note: vertex_head points to a circular list at this point. */
  /*       We delete starting with vertex_head->next, and nil   */
  /*       that pointer to break the cycle in the list.         */
  ppa = vertex_head->next;
  vertex_head->next = NULL;

  /* Flatten out lists so that "span" elements are included... */
  for (ppb = ppa; ppb != NULL; ppb = ppb->next) {
    if (ppb->span != NULL) {
      struct exd_vertex *next = ppb->next;
      ppb->next = ppb->span;
      ppb->span = NULL;
      while (ppb->next != NULL)
        ppb = ppb->next;
      ppb->next = next;
    }
  }
  mem_put_list(sv->local_storage->exdv, ppa);

  /* Return fractional area */

  return A / (MY_PI * R2);

#undef EXD_TIME_CALC
#undef EXD_SPAN_CALC
}

/**************************/
/** done with exact_disk **/
/**************************/

/****************************************************************************
safe_diffusion_step:
  In: vm: molecule that is moving
      shead: linked list of potential collisions with molecules from the
      radial_subdivisions:
      r_step:
      x_fineparts:
      y_fineparts:
      z_fineparts:
  Out: The estimated number of diffusion steps this molecule can take before
       something interesting might happen to it, or 1.0 if something might
       happen within one timestep.
  Note: Each molecule uses its own timestep.  Only molecules that the moving
        molecule can react with directly are counted (secondary reaction
        products are ignored, so might be skipped).  "Might happen" is to
        the 99% confidence level (i.e. the distance you'd have to go before
        1% of the molecules will have gotten far enough to have a chance of
        reacting, although those 1% will probably not go in the right
        direction).  This doesn't take into account the diffusion of other
        target molecules, so it may introduce errors for clouds of molecules
        diffusing into each other from a distance.
        *FIXME*: Add a flag to make this be very conservative or to turn
        this off entirely, aside from the TIME_STEP_MAX= directive.
****************************************************************************/
double safe_diffusion_step(struct volume_molecule *vm, struct collision *shead,
                           u_int radial_subdivisions, double *r_step,
                           double *x_fineparts, double *y_fineparts,
                           double *z_fineparts) {
  double d2;
  double d2_nearmax;
  double d2min = GIGANTIC;
  struct subvolume *sv = vm->subvol;
  struct wall *w;
  struct wall_list *wl;
  struct collision *smash;
  double steps;
  struct volume_molecule *mp;

  d2_nearmax = vm->properties->space_step *
               r_step[(int)(radial_subdivisions * MULTISTEP_PERCENTILE)];
  d2_nearmax *= d2_nearmax;

  if ((vm->properties->flags & (CAN_VOLVOL | CANT_INITIATE)) == CAN_VOLVOL) {
    for (smash = shead; smash != NULL; smash = smash->next) {
      mp = (struct volume_molecule *)smash->target;
      d2 = (vm->pos.x - mp->pos.x) * (vm->pos.x - mp->pos.x) +
           (vm->pos.y - mp->pos.y) * (vm->pos.y - mp->pos.y) +
           (vm->pos.z - mp->pos.z) * (vm->pos.z - mp->pos.z);
      if (d2 < d2min)
        d2min = d2;
    }
  }
  for (wl = sv->wall_head; wl != NULL; wl = wl->next) {
    w = wl->this_wall;
    d2 = (w->normal.x * vm->pos.x + w->normal.y * vm->pos.y +
          w->normal.z * vm->pos.z) -
         w->d;
    d2 *= d2;
    if (d2 < d2min)
      d2min = d2;
  }

  d2 = (vm->pos.x - x_fineparts[sv->llf.x]);
  d2 *= d2;
  if (d2 < d2min)
    d2min = d2;

  d2 = (vm->pos.x - x_fineparts[sv->urb.x]);
  d2 *= d2;
  if (d2 < d2min)
    d2min = d2;

  d2 = (vm->pos.y - y_fineparts[sv->llf.y]);
  d2 *= d2;
  if (d2 < d2min)
    d2min = d2;

  d2 = (vm->pos.y - y_fineparts[sv->urb.y]);
  d2 *= d2;
  if (d2 < d2min)
    d2min = d2;

  d2 = (vm->pos.z - z_fineparts[sv->llf.z]);
  d2 *= d2;
  if (d2 < d2min)
    d2min = d2;

  d2 = (vm->pos.z - z_fineparts[sv->urb.z]);
  d2 *= d2;
  if (d2 < d2min)
    d2min = d2;

  if (d2min < d2_nearmax)
    steps = 1.0;
  else {
    double steps_sq = d2min / d2_nearmax;
    if (steps_sq < MULTISTEP_WORTHWHILE * MULTISTEP_WORTHWHILE)
      steps = 1.0;
    else
      steps = sqrt(steps_sq);
  }

  return steps;
}

/****************************************************************************
expand_collision_list_for_neighbor:
  This is a helper function to reduce duplicated code in expand_collision_list.
  This code will be called once for each adjacent subvolume.  The clipping
  indicators below (trim_[xyz]) are interpreted as follows:

    trim < 0: We went in the negative direction for this axis.  Search within
              -trim of the maximal partition boundary in the adjacent subvol
    trim > 0: We went in the positive direction for this axis.  Search within
              trim of the minimal partition boundary in the adjacent subvol
    trim = 0: The subvolume is adjacent along this axis.  Search the entire
              width of this axis of the subvolume.

  In: sv: the "current" subvolume
      vm: the current molecule
      new_sv: adjacent subvolume to search
      path_llf: path bounding box lower left front
      path_urb: path bounding box upper right back
      shead1: current list head
      trim_x: X clipping indicator
      trim_y: Y clipping indicator
      trim_z: Z clipping indicator
      x_fineparts:
      y_fineparts:
      z_fineparts:
      rx_hashsize:
      reaction_hash:
  Out: Returns linked list of molecules from neighbor subvolumes
       that are located within "interaction_radius" from the the subvolume
       border.
       The molecules are added only when the molecule displacement
       bounding box intersects with the subvolume bounding box.
****************************************************************************/
static struct collision *expand_collision_list_for_neighbor(
    struct subvolume *sv, struct volume_molecule *vm, struct subvolume *new_sv,
    struct vector3 *path_llf, struct vector3 *path_urb,
    struct collision *shead1, double trim_x, double trim_y, double trim_z,
    double *x_fineparts, double *y_fineparts, double *z_fineparts,
    int rx_hashsize, struct rxn **reaction_hash) {
  int num_matching_rxns = 0;
  struct rxn *matching_rxns[MAX_MATCHING_RXNS];

  /* Grab the subvolume boundaries */
  struct vector3 new_sv_llf, new_sv_urb;
  new_sv_llf.x = x_fineparts[new_sv->llf.x];
  new_sv_llf.y = y_fineparts[new_sv->llf.y];
  new_sv_llf.z = z_fineparts[new_sv->llf.z];
  new_sv_urb.x = x_fineparts[new_sv->urb.x];
  new_sv_urb.y = y_fineparts[new_sv->urb.y];
  new_sv_urb.z = z_fineparts[new_sv->urb.z];

  /* Quickly check if the subvolume bounds and the path bounds intersect */
  if (!test_bounding_boxes(path_llf, path_urb, &new_sv_llf, &new_sv_urb))
    return shead1;

  /* Find the bounds for molecules to check */
  double x_min, x_max;
  double y_min, y_max;
  double z_min, z_max;
  if (trim_x < 0.0) {
    x_min = new_sv_urb.x + trim_x;
    x_max = new_sv_urb.x + EPS_C;
  } else if (trim_x > 0.0) {
    x_min = new_sv_llf.x - EPS_C;
    x_max = new_sv_llf.x + trim_x;
  } else {
    x_min = new_sv_llf.x - EPS_C;
    x_max = new_sv_urb.x + EPS_C;
  }
  if (trim_y < 0.0) {
    y_min = new_sv_urb.y + trim_y;
    y_max = new_sv_urb.y + EPS_C;
  } else if (trim_y > 0.0) {
    y_min = new_sv_llf.y - EPS_C;
    y_max = new_sv_llf.y + trim_y;
  } else {
    y_min = new_sv_llf.y - EPS_C;
    y_max = new_sv_urb.y + EPS_C;
  }
  if (trim_z < 0.0) {
    z_min = new_sv_urb.z + trim_z;
    z_max = new_sv_urb.z + EPS_C;
  } else if (trim_z > 0.0) {
    z_min = new_sv_llf.z - EPS_C;
    z_max = new_sv_llf.z + trim_z;
  } else {
    z_min = new_sv_llf.z - EPS_C;
    z_max = new_sv_urb.z + EPS_C;
  }

  /* scan molecules from this SV */
  struct per_species_list *psl_next, *psl, **psl_head = &new_sv->species_head;
  for (psl = new_sv->species_head; psl != NULL; psl = psl_next) {
    psl_next = psl->next;
    if (psl->properties == NULL) {
      psl_head = &psl->next;
      continue;
    }

    /* Garbage collection of empty per-species lists */
    if (psl->head == NULL) {
      *psl_head = psl->next;
      ht_remove(&new_sv->mol_by_species, psl);
      mem_put(new_sv->local_storage->pslv, psl);
      continue;
    } else
      psl_head = &psl->next;

    /* no possible reactions. skip it. */
    if (!trigger_bimolecular_preliminary(
             reaction_hash, rx_hashsize, vm->properties->hashval,
             psl->properties->hashval, vm->properties, psl->properties))
      continue;

    for (struct volume_molecule *mp = psl->head; mp != NULL; mp = mp->next_v) {
      /* Skip defunct molecules */
      if (mp->properties == NULL)
        continue;

      /* skip molecules outside the region of interest */
      if (mp->pos.x < x_min || mp->pos.x > x_max)
        continue;
      if (mp->pos.y < y_min || mp->pos.y > y_max)
        continue;
      if (mp->pos.z < z_min || mp->pos.z > z_max)
        continue;
      // count only in the relevant periodic box
      if (!periodic_boxes_are_identical(vm->periodic_box, mp->periodic_box)) {
        continue;
      }

      /* check for possible reactions */
      num_matching_rxns = trigger_bimolecular(
          reaction_hash, rx_hashsize, vm->properties->hashval,
          mp->properties->hashval, (struct abstract_molecule *)vm,
          (struct abstract_molecule *)mp, 0, 0, matching_rxns);
      if (num_matching_rxns <= 0)
        continue;

      /* Add a collision for each matching reaction */
      for (int i = 0; i < num_matching_rxns; i++) {
        struct collision *smash = (struct collision *)CHECKED_MEM_GET(
            sv->local_storage->coll, "collision data");
        smash->target = (void *)mp;
        smash->intermediate = matching_rxns[i];
        smash->next = shead1;
        smash->what = 0;
        smash->what |= COLLIDE_VOL;
        shead1 = smash;
      }
    }
  }

  return shead1;
}

/****************************************************************************
expand_collision_list:
  In: vm: molecule that is moving
      mv: displacement to the new location
      sv: subvolume that we start in
      rx_radius_3d:
      ny_parts:
      nz_parts:
      x_fineparts:
      y_fineparts:
      z_fineparts:
  Out: Returns list of collisions with molecules from neighbor subvolumes
       that are located within "interaction_radius" from the the subvolume
       border.  The molecules are added only when the molecule displacement
       bounding box intersects with the subvolume bounding box.
****************************************************************************/
static struct collision *
expand_collision_list(struct volume_molecule *vm, struct vector3 *mv,
                      struct subvolume *sv, double rx_radius_3d,
                      int ny_parts, int nz_parts, double *x_fineparts,
                      double *y_fineparts, double *z_fineparts, int rx_hashsize,
                      struct rxn **reaction_hash) {
  struct collision *shead1 = NULL;
  /* neighbors of the current subvolume */
  struct vector3 path_llf, path_urb;
  double R = (rx_radius_3d);

  /* find the molecule path bounding box. */
  path_bounding_box(&vm->pos, mv, &path_llf, &path_urb, rx_radius_3d);

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

  /* go in the direction X_POS */
  if (x_pos) {
    struct subvolume *new_sv = sv + (nz_parts - 1) * (ny_parts - 1);
    shead1 = expand_collision_list_for_neighbor(
        sv, vm, new_sv, &path_llf, &path_urb, shead1, R, 0.0, 0.0, x_fineparts,
        y_fineparts, z_fineparts, rx_hashsize, reaction_hash);

    /* go +X, +Y) */
    if (y_pos) {
      struct subvolume *new_sv_y = new_sv + (nz_parts - 1);
      shead1 = expand_collision_list_for_neighbor(
          sv, vm, new_sv_y, &path_llf, &path_urb, shead1, R, R, 0.0, x_fineparts,
          y_fineparts, z_fineparts, rx_hashsize, reaction_hash);

      /* go +X, +Y, +Z) */
      if (z_pos)
        shead1 = expand_collision_list_for_neighbor(
            sv, vm, new_sv_y + 1, &path_llf, &path_urb, shead1, R, R, R,
            x_fineparts, y_fineparts, z_fineparts, rx_hashsize, reaction_hash);

      /* go +X, +Y, -Z */
      if (z_neg)
        shead1 = expand_collision_list_for_neighbor(
            sv, vm, new_sv_y - 1, &path_llf, &path_urb, shead1, R, R, -R,
            x_fineparts, y_fineparts, z_fineparts, rx_hashsize, reaction_hash);
    }

    /* go +X, -Y) */
    if (y_neg) {
      struct subvolume *new_sv_y = new_sv - (nz_parts - 1);
      shead1 = expand_collision_list_for_neighbor(
          sv, vm, new_sv_y, &path_llf, &path_urb, shead1, R, -R, 0.0,
          x_fineparts, y_fineparts, z_fineparts, rx_hashsize, reaction_hash);

      /* go +X, -Y, +Z) */
      if (z_pos)
        shead1 = expand_collision_list_for_neighbor(
            sv, vm, new_sv_y + 1, &path_llf, &path_urb, shead1, R, -R, R,
            x_fineparts, y_fineparts, z_fineparts, rx_hashsize, reaction_hash);

      /* go +X, -Y, -Z */
      if (z_neg)
        shead1 = expand_collision_list_for_neighbor(
            sv, vm, new_sv_y - 1, &path_llf, &path_urb, shead1, R, -R, -R,
            x_fineparts, y_fineparts, z_fineparts, rx_hashsize, reaction_hash);
    }

    /* go +X, +Z) */
    if (z_pos)
      shead1 = expand_collision_list_for_neighbor(
          sv, vm, new_sv + 1, &path_llf, &path_urb, shead1, R, 0.0, R,
          x_fineparts, y_fineparts, z_fineparts, rx_hashsize, reaction_hash);

    /* go +X, -Z */
    if (z_neg)
      shead1 = expand_collision_list_for_neighbor(
          sv, vm, new_sv - 1, &path_llf, &path_urb, shead1, R, 0.0, -R,
          x_fineparts, y_fineparts, z_fineparts, rx_hashsize, reaction_hash);
  }

  /* go in the direction X_NEG */
  if (x_neg) {
    struct subvolume *new_sv = sv - (nz_parts - 1) * (ny_parts - 1);
    shead1 = expand_collision_list_for_neighbor(
        sv, vm, new_sv, &path_llf, &path_urb, shead1, -R, 0.0, 0.0, x_fineparts,
        y_fineparts, z_fineparts, rx_hashsize, reaction_hash);

    /* go -X, +Y) */
    if (y_pos) {
      struct subvolume *new_sv_y = new_sv + (nz_parts - 1);
      shead1 = expand_collision_list_for_neighbor(
          sv, vm, new_sv_y, &path_llf, &path_urb, shead1, -R, R, 0.0,
          x_fineparts, y_fineparts, z_fineparts, rx_hashsize, reaction_hash);

      /* go -X, +Y, +Z) */
      if (z_pos)
        shead1 = expand_collision_list_for_neighbor(
            sv, vm, new_sv_y + 1, &path_llf, &path_urb, shead1, -R, R, R,
            x_fineparts, y_fineparts, z_fineparts, rx_hashsize, reaction_hash);

      /* go -X, +Y, -Z */
      if (z_neg)
        shead1 = expand_collision_list_for_neighbor(
            sv, vm, new_sv_y - 1, &path_llf, &path_urb, shead1, -R, R, -R,
            x_fineparts, y_fineparts, z_fineparts, rx_hashsize, reaction_hash);
    }

    /* go -X, -Y) */
    if (y_neg) {
      struct subvolume *new_sv_y = new_sv - (nz_parts - 1);
      shead1 = expand_collision_list_for_neighbor(
          sv, vm, new_sv_y, &path_llf, &path_urb, shead1, -R, -R, 0.0,
          x_fineparts, y_fineparts, z_fineparts, rx_hashsize, reaction_hash);

      /* go -X, -Y, +Z) */
      if (z_pos)
        shead1 = expand_collision_list_for_neighbor(
            sv, vm, new_sv_y + 1, &path_llf, &path_urb, shead1, -R, -R, R,
            x_fineparts, y_fineparts, z_fineparts, rx_hashsize, reaction_hash);

      /* go -X, -Y, -Z */
      if (z_neg)
        shead1 = expand_collision_list_for_neighbor(
            sv, vm, new_sv_y - 1, &path_llf, &path_urb, shead1, -R, -R, -R,
            x_fineparts, y_fineparts, z_fineparts, rx_hashsize, reaction_hash);
    }

    /* go -X, +Z) */
    if (z_pos)
      shead1 = expand_collision_list_for_neighbor(
          sv, vm, new_sv + 1, &path_llf, &path_urb, shead1, -R, 0.0, R,
          x_fineparts, y_fineparts, z_fineparts, rx_hashsize, reaction_hash);

    /* go -X, -Z */
    if (z_neg)
      shead1 = expand_collision_list_for_neighbor(
          sv, vm, new_sv - 1, &path_llf, &path_urb, shead1, -R, 0.0, -R,
          x_fineparts, y_fineparts, z_fineparts, rx_hashsize, reaction_hash);
  }

  /* go in the direction Y_POS */
  if (y_pos) {
    struct subvolume *new_sv = sv + (nz_parts - 1);
    shead1 = expand_collision_list_for_neighbor(
        sv, vm, new_sv, &path_llf, &path_urb, shead1, 0.0, R, 0.0, x_fineparts,
        y_fineparts, z_fineparts, rx_hashsize, reaction_hash);

    /* go +Y, +Z) */
    if (z_pos)
      shead1 = expand_collision_list_for_neighbor(
          sv, vm, new_sv + 1, &path_llf, &path_urb, shead1, 0.0, R, R,
          x_fineparts, y_fineparts, z_fineparts, rx_hashsize, reaction_hash);

    /* go +Y, -Z */
    if (z_neg)
      shead1 = expand_collision_list_for_neighbor(
          sv, vm, new_sv - 1, &path_llf, &path_urb, shead1, 0.0, R, -R,
          x_fineparts, y_fineparts, z_fineparts, rx_hashsize, reaction_hash);
  }

  /* go in the direction Y_NEG */
  if (y_neg) {
    struct subvolume *new_sv = sv - (nz_parts - 1);
    shead1 = expand_collision_list_for_neighbor(
        sv, vm, new_sv, &path_llf, &path_urb, shead1, 0.0, -R, 0.0, x_fineparts,
        y_fineparts, z_fineparts, rx_hashsize, reaction_hash);

    /* go -Y, +Z) */
    if (z_pos)
      shead1 = expand_collision_list_for_neighbor(
          sv, vm, new_sv + 1, &path_llf, &path_urb, shead1, 0.0, -R, R,
          x_fineparts, y_fineparts, z_fineparts, rx_hashsize, reaction_hash);

    /* go -Y, -Z */
    if (z_neg)
      shead1 = expand_collision_list_for_neighbor(
          sv, vm, new_sv - 1, &path_llf, &path_urb, shead1, 0.0, -R, -R,
          x_fineparts, y_fineparts, z_fineparts, rx_hashsize, reaction_hash);
  }

  /* go in the direction Z_POS */
  if (z_pos)
    shead1 = expand_collision_list_for_neighbor(
        sv, vm, sv + 1, &path_llf, &path_urb, shead1, 0.0, 0.0, R, x_fineparts,
        y_fineparts, z_fineparts, rx_hashsize, reaction_hash);

  /* go in the direction Z_NEG */
  if (z_neg)
    shead1 = expand_collision_list_for_neighbor(
        sv, vm, sv - 1, &path_llf, &path_urb, shead1, 0.0, 0.0, -R, x_fineparts,
        y_fineparts, z_fineparts, rx_hashsize, reaction_hash);

  return shead1;
}

/****************************************************************************
expand_collision_partner_list_for_neighbor:
  This is a helper function to reduce duplicated code in
  expand_collision_partner_list.  This code will be called once for each
  adjacent subvolume.  The clipping indicators below (trim_[xyz]) are
  interpreted as follows:

    trim < 0: We went in the negative direction for this axis.  Search within
              -trim of the maximal partition boundary in the adjacent subvol
    trim > 0: We went in the positive direction for this axis.  Search within
              trim of the minimal partition boundary in the adjacent subvol
    trim = 0: The subvolume is adjacent along this axis.  Search the entire
              width of this axis of the subvolume.

  In: struct subvolume *sv - the "current" subvolume
      struct volume_molecule *vm - the current molecule
      struct vector3 *mv - displacement to the new location
      struct subvolume *new_sv - adjacent subvolume to search
      struct vector3 *path_llf - path bounding box lower left front
      struct vector3 *path_urb - path bounding box upper right back
      struct sp_collision *shead1 - current list head
      double trim_x - X clipping indicator
      double trim_y - Y clipping indicator
      double trim_z - Z clipping indicator
  Out: Returns linked list of molecules from neighbor subvolumes
       that are located within "interaction_radius" from the the subvolume
       border.
       The molecules are added only when the molecule displacement
       bounding box intersects with the subvolume bounding box.
****************************************************************************/
struct sp_collision *expand_collision_partner_list_for_neighbor(
    struct subvolume *sv, struct volume_molecule *vm, struct vector3 *mv,
    struct subvolume *new_sv, struct vector3 *path_llf,
    struct vector3 *path_urb, struct sp_collision *shead1, double trim_x,
    double trim_y, double trim_z, double *x_fineparts, double *y_fineparts,
    double *z_fineparts, int rx_hashsize, struct rxn **reaction_hash) {
  struct species *spec = vm->properties;
  struct sp_collision *smash;

  /* Grab the subvolume boundaries */
  struct vector3 new_sv_llf, new_sv_urb;
  new_sv_llf.x = x_fineparts[new_sv->llf.x];
  new_sv_llf.y = y_fineparts[new_sv->llf.y];
  new_sv_llf.z = z_fineparts[new_sv->llf.z];
  new_sv_urb.x = x_fineparts[new_sv->urb.x];
  new_sv_urb.y = y_fineparts[new_sv->urb.y];
  new_sv_urb.z = z_fineparts[new_sv->urb.z];

  /* Quickly check if the subvolume bounds and the path bounds intersect */
  if (!test_bounding_boxes(path_llf, path_urb, &new_sv_llf, &new_sv_urb))
    return shead1;

  int moving_tri_molecular_flag = 0, moving_bi_molecular_flag = 0,
      moving_mol_mol_grid_flag = 0;
  /* collision flags */
  int col_tri_molecular_flag = 0, col_bi_molecular_flag = 0,
      col_mol_mol_grid_flag = 0;

  moving_tri_molecular_flag =
      ((spec->flags & (CAN_VOLVOLVOL | CANT_INITIATE)) == CAN_VOLVOLVOL);
  moving_bi_molecular_flag =
      ((spec->flags & (CAN_VOLVOL | CANT_INITIATE)) == CAN_VOLVOL);
  moving_mol_mol_grid_flag =
      ((spec->flags & (CAN_VOLVOLSURF | CANT_INITIATE)) == CAN_VOLVOLSURF);

  /* Find the bounds for molecules to check */
  double x_min, x_max;
  double y_min, y_max;
  double z_min, z_max;
  if (trim_x < 0.0) {
    x_min = new_sv_urb.x + trim_x;
    x_max = new_sv_urb.x + EPS_C;
  } else if (trim_x > 0.0) {
    x_min = new_sv_llf.x - EPS_C;
    x_max = new_sv_llf.x + trim_x;
  } else {
    x_min = new_sv_llf.x - EPS_C;
    x_max = new_sv_urb.x + EPS_C;
  }
  if (trim_y < 0.0) {
    y_min = new_sv_urb.y + trim_y;
    y_max = new_sv_urb.y + EPS_C;
  } else if (trim_y > 0.0) {
    y_min = new_sv_llf.y - EPS_C;
    y_max = new_sv_llf.y + trim_y;
  } else {
    y_min = new_sv_llf.y - EPS_C;
    y_max = new_sv_urb.y + EPS_C;
  }
  if (trim_z < 0.0) {
    z_min = new_sv_urb.z + trim_z;
    z_max = new_sv_urb.z + EPS_C;
  } else if (trim_z > 0.0) {
    z_min = new_sv_llf.z - EPS_C;
    z_max = new_sv_llf.z + trim_z;
  } else {
    z_min = new_sv_llf.z - EPS_C;
    z_max = new_sv_urb.z + EPS_C;
  }

  /* scan molecules from this SV */
  struct per_species_list *psl_next, *psl, **psl_head = &new_sv->species_head;
  for (psl = new_sv->species_head; psl != NULL; psl = psl_next) {
    psl_next = psl->next;
    if (psl->properties == NULL) {
      psl_head = &psl->next;
      continue;
    }

    /* Garbage collection of empty per-species lists */
    if (psl->head == NULL) {
      *psl_head = psl->next;
      ht_remove(&new_sv->mol_by_species, psl);
      mem_put(new_sv->local_storage->pslv, psl);
      continue;
    } else
      psl_head = &psl->next;

    col_tri_molecular_flag =
        moving_tri_molecular_flag &&
        ((psl->properties->flags & CAN_VOLVOLVOL) == CAN_VOLVOLVOL);
    col_bi_molecular_flag =
        moving_bi_molecular_flag &&
        ((psl->properties->flags & CAN_VOLVOL) == CAN_VOLVOL) &&
        trigger_bimolecular_preliminary(reaction_hash, rx_hashsize,
                                        spec->hashval, psl->properties->hashval,
                                        spec, psl->properties);
    col_mol_mol_grid_flag =
        moving_mol_mol_grid_flag &&
        ((psl->properties->flags & CAN_VOLVOLSURF) == CAN_VOLVOLSURF);
    if (col_bi_molecular_flag || col_tri_molecular_flag ||
        col_mol_mol_grid_flag) {
      struct volume_molecule *mp;
      for (mp = psl->head; mp != NULL; mp = mp->next_v) {
        /* Skip defunct molecules */
        if (mp->properties == NULL)
          continue;

        /* skip molecules outside the region of interest */
        if (mp->pos.x < x_min || mp->pos.x > x_max)
          continue;
        if (mp->pos.y < y_min || mp->pos.y > y_max)
          continue;
        if (mp->pos.z < z_min || mp->pos.z > z_max)
          continue;

        smash = (struct sp_collision *)CHECKED_MEM_GET(
            sv->local_storage->sp_coll, "collision data");
        smash->t = 0.0;
        smash->t_start = 0.0;
        smash->pos_start.x = vm->pos.x;
        smash->pos_start.y = vm->pos.y;
        smash->pos_start.z = vm->pos.z;
        smash->sv_start = sv;
        smash->disp.x = mv->x;
        smash->disp.y = mv->y;
        smash->disp.z = mv->z;
        smash->loc.x = 0.0;
        smash->loc.y = 0.0;
        smash->loc.z = 0.0;
        smash->moving = spec;
        smash->target = (void *)mp;
        smash->what = 0;
        if (col_bi_molecular_flag) {
          smash->what |= COLLIDE_VOL;
        }
        if (col_tri_molecular_flag) {
          smash->what |= COLLIDE_VOL_VOL;
        }
        if (col_mol_mol_grid_flag) {
          smash->what |= COLLIDE_VOL_SURF;
        }
        smash->next = shead1;
        shead1 = smash;
      }
    }
  }
  return shead1;
}

/*************************************************************************
diffuse_3D:
  In: world: simulation state
      vm: molecule that is moving
      max_time: maximum time we can spend diffusing
  Out: Pointer to the molecule if it still exists (may have been
       reallocated), NULL otherwise.
       Position and time are updated, but molecule is not rescheduled.
  Note: This version takes into account only 2-way reactions and 3-way
        reactions of type MOL_GRID_GRID
*************************************************************************/
struct volume_molecule *diffuse_3D(
    struct volume *world,
    struct volume_molecule *vm,
    double max_time) {

  struct species* spec = vm->properties;
  if (spec == NULL) {
    mcell_internal_error(
        "Attempted to take a diffusion step for a defunct molecule.");
  }

  /* flags related to the possible reaction between volume molecule
     and one or two surface molecules */
  int mol_grid_flag = ((spec->flags & CAN_VOLSURF) == CAN_VOLSURF);
  int mol_grid_grid_flag = ((spec->flags & CAN_VOLSURFSURF) == CAN_VOLSURFSURF);

  if (spec->space_step <= 0.0) {
    vm->t += max_time;
    return vm;
  }

  int inertness = 0;
  set_inertness_and_maxtime(world, vm, &max_time, &inertness);

  /* Done housekeeping, now let's do something fun! */
  int calculate_displacement = 1;

  /* this flag is set to 1 only after reflection from a wall and only with
   * expanded lists. */
  int redo_expand_collision_list_flag = 0;

  /* initialize before fake recursion */
  double steps = 1.0;
  double t_steps = 1.0;
  double rate_factor = 1.0;
  double r_rate_factor = 1.0;
  struct vector3 displacement;  /* Molecule moves along this vector */
  struct vector3 displacement2; /* Used for 3D mol-mol unbinding */

pretend_to_call_diffuse_3D: ; /* Label to allow fake recursion */

  struct subvolume *sv = vm->subvol;
  struct collision *shead = NULL; /* Things we might hit (can interact with) */
  struct collision *stail = NULL; /* Things we might hit (can interact with -
                                     tail of the collision linked list) */
  struct collision *shead_exp = NULL; /* Things we might hit (can interact with)
                                         from neighbor subvolumes */
  /* scan subvolume for possible mol-mol reactions with vm */
  if ((spec->flags & (CAN_VOLVOL | CANT_INITIATE)) == CAN_VOLVOL &&
      inertness < inert_to_all) {
    determine_mol_mol_reactions(world, vm, &shead, &stail, inertness);
  }

  if (calculate_displacement) {
    compute_displacement(world, shead, vm, &displacement, &displacement2,
      &rate_factor, &r_rate_factor, &steps, &t_steps, max_time);
  }

  if (world->use_expanded_list &&
      ((vm->properties->flags & (CAN_VOLVOL | CANT_INITIATE)) == CAN_VOLVOL) &&
      !inertness) {
    shead_exp = expand_collision_list(
      vm, &displacement, sv, world->rx_radius_3d, world->ny_parts,
      world->nz_parts, world->x_fineparts, world->y_fineparts,
      world->z_fineparts, world->rx_hashsize, world->reaction_hash);
    if (stail != NULL)
      stail->next = shead_exp;
    else {
      if (shead != NULL)
        mcell_internal_error("Collision lists corrupted.  While expanding the "
                             "collision lists, expected shead to be NULL, but "
                             "it wasn't.");
      shead = shead_exp;
    }
  }

  struct wall* reflectee = NULL;
  struct collision *smash;      /* Thing we've hit that's under consideration */
  do {
    /* due to redo_expand_collision_list_flag this only happens after reflection */
    if (world->use_expanded_list && redo_expand_collision_list_flag) {
      redo_collision_list(world, &shead, &stail, &shead_exp, vm, &displacement, sv);
    }

    struct collision* shead2 = ray_trace(world, &(vm->pos), shead, sv, &displacement, reflectee);
    if (shead2 == NULL) {
      mcell_internal_error("ray_trace returned NULL.");
    }

    if (shead2->next != NULL) {
      shead2 =
          (struct collision *)ae_list_sort((struct abstract_element *)shead2);
    }

    struct vector3* loc_certain = NULL;
    struct collision *tentative = shead2;

    for (smash = shead2; smash != NULL; smash = smash->next) {
      if (world->notify->molecule_collision_report == NOTIFY_FULL) {
        if (((smash->what & COLLIDE_VOL) != 0) &&
            (world->rxn_flags.vol_vol_reaction_flag)) {
          world->vol_vol_colls++;
        }
      }

      if (smash->t >= 1.0 || smash->t < 0.0) {
        if ((smash->what & COLLIDE_VOL) != 0) {
          mcell_internal_error(
              "Detected a mol-mol collision outside of the 0.0...1.0 time "
              "window.  Iteration %lld, time of collision %.8e, mol1=%s, "
              "mol2=%s",
              world->current_iterations, smash->t, vm->properties->sym->name,
              ((struct volume_molecule *)smash->target)->properties->sym->name);
        }
        smash = NULL;
        break;
      }

      if ((smash->what & COLLIDE_VOL) != 0) {
        if (smash->t < EPS_C) {
          continue;
        }
        if (collide_and_react_with_vol_mol(world, smash, vm, &tentative,
          &displacement, loc_certain, t_steps, r_rate_factor) == 1) {
          FREE_COLLISION_LISTS();
          return NULL;
        } else {
          continue;
        }
      } else if ((smash->what & COLLIDE_WALL) != 0) {

        struct wall* w = (struct wall *)smash->target;
        if (w->grid != NULL && (mol_grid_flag || mol_grid_grid_flag) &&
          inertness < inert_to_all) {
          int destroyed = collide_and_react_with_surf_mol(world, smash, vm,
            &tentative, &loc_certain, t_steps, mol_grid_flag, mol_grid_grid_flag,
            r_rate_factor);
          // if destroyed = -1 we didn't react with any molecules and keep going
          // to check for wall collisions
          if (destroyed == 1) {
            FREE_COLLISION_LISTS();
            return NULL;
          } else if (destroyed == 0) {
            continue;
          }
        }

        if ((spec->flags & CAN_VOLWALL) != 0) {
          int destroyed = collide_and_react_with_walls(world, smash, vm,
            &tentative, &loc_certain, t_steps, inertness, r_rate_factor);
          // if destroyed = -1 we didn't react with any walls and keep going to
          // either reflect or encounter periodic bc
          if (destroyed == 1) {
            FREE_COLLISION_LISTS();
            return NULL;
          } else if (destroyed == 0) {
            continue;
          }
        }

        if (reflect_or_periodic_bc(world, smash, &displacement, &vm, &reflectee,
            &tentative, &t_steps) == 1) {
          FREE_COLLISION_LISTS();
          calculate_displacement = 0;
          if (vm->properties == NULL) {
            mcell_internal_error("A defunct molecule is diffusing.");
          }
          // molecule was reflected into periodic box. Need to start over
          // figuring out targets based on the current displacement
          goto pretend_to_call_diffuse_3D;
        }
        // Only useful if using expanded lists, but easier to always set it
        redo_expand_collision_list_flag = 1; 
                                             
        break;
      } else if ((smash->what & COLLIDE_SUBVOL) != 0) {
        collide_and_react_with_subvol(
          world, smash, &displacement, &vm, &tentative, &t_steps);
        FREE_COLLISION_LISTS();
        calculate_displacement = 0;

        if (vm->properties == NULL) {
          mcell_internal_error("A defunct molecule is diffusing.");
        }
        // molecule entered a new subvolume. Continue motion from the beginning.
        goto pretend_to_call_diffuse_3D;
      }
    }

    if (shead2 != NULL) {
      mem_put_list(sv->local_storage->coll, shead2);
    }
  } while (smash != NULL);

  vm->pos.x += displacement.x;
  vm->pos.y += displacement.y;
  vm->pos.z += displacement.z;
  vm->t += t_steps;

  /* Done with traversing disk, now do real motion */
  if (inertness == inert_to_all) 
  {
    inertness = inert_to_mol;
    t_steps = spec->time_step;
    displacement = displacement2;
    calculate_displacement = 0;
    goto pretend_to_call_diffuse_3D;
  }

  vm->index = -1;
  vm->previous_wall = NULL;

  if (shead != NULL)
    mem_put_list(sv->local_storage->coll, shead);

  return vm;
}

/*************************************************************************
move_sm_on_same_triangle:

  This is a helper function for diffuse_2D.

  In: world: simulation state
      sm: molecule that is moving
      new_loc: this is the location we are moving to.
      previous_box: this is the periodic box we were in previously.
      new_wall: this is the new wall we ended up on
      hd_info:
  Out: The grid is created on a new triangle and we place the molecule if
       possible. Counts are updated.
*************************************************************************/
int move_sm_on_same_triangle(
    struct volume *state,
    struct surface_molecule *sm,
    struct vector2 *new_loc,
    struct periodic_image *previous_box,
    struct wall *new_wall,
    struct hit_data *hd_info) {
  unsigned int new_idx = uv2grid(new_loc, new_wall->grid);
  if (new_idx >= sm->grid->n_tiles) {
    mcell_internal_error("After ray_trace_2D, selected u, v coordinates "
                         "map to an out-of-bounds grid cell.  uv=(%.2f, "
                         "%.2f) sm=%d/%d",
                         new_loc->u, new_loc->v, new_idx, sm->grid->n_tiles);
  }
  // We're on a new part of the grid
  struct surface_molecule_list *sm_list = sm->grid->sm_list[new_idx];
  if (new_idx != sm->grid_index) {
    if ((state->periodic_box_obj && periodicbox_in_surfmol_list(sm->periodic_box, sm_list)) ||
        (!state->periodic_box_obj && sm_list && sm_list->sm)) {
      if (hd_info != NULL) {
        delete_void_list((struct void_list *)hd_info);
        hd_info = NULL;
      }
      return 1; /* Pick again--full here */
    }

    remove_surfmol_from_list(&sm->grid->sm_list[sm->grid_index], sm);
    sm->grid_index = new_idx;
    sm->grid->sm_list[new_idx] = add_surfmol_with_unique_pb_to_list(
      sm->grid->sm_list[new_idx], sm);
    assert(sm->grid->sm_list[new_idx] != NULL);
    count_moved_surface_mol(
      state, sm, sm->grid, new_loc, state->count_hashmask,
      state->count_hash, &state->ray_polygon_colls, previous_box);
  // We ended up on the same exact grid element! 
  // XXX: do we even need to update counts??
  } else {
    count_moved_surface_mol(
      state, sm, sm->grid, new_loc, state->count_hashmask,
      state->count_hash, &state->ray_polygon_colls, previous_box);
  }

  sm->s_pos.u = new_loc->u;
  sm->s_pos.v = new_loc->v;
  return 0;
}

/*************************************************************************
move_sm_to_new_triangle: 

  This is a helper function for diffuse_2D.

  In: world: simulation state
      sm: molecule that is moving
      new_loc: this is the location we are moving to.
      previous_box: this is the periodic box we were in previously.
      new_wall: this is the new wall we ended up on
      hd_info:
  Out: The grid is created on a new triangle and we place the molecule if
       possible. Counts are updated.
*************************************************************************/
int move_sm_to_new_triangle(
    struct volume *state,
    struct surface_molecule *sm,
    struct vector2 *new_loc,
    struct periodic_image *previous_box,
    struct wall *new_wall,
    struct hit_data *hd_info) {
  // No SM has been here before, so we need to make a grid on this wall.
  if (new_wall->grid == NULL) {
    if (create_grid(state, new_wall, NULL))
      mcell_allocfailed("Failed to create a grid for a wall.");
  }

  /* Move to new tile */
  unsigned int new_idx = uv2grid(new_loc, new_wall->grid);
  if (new_idx >= new_wall->grid->n_tiles) {
    mcell_internal_error(
        "After ray_trace_2D to a new wall, selected u, v coordinates map "
        "to an out-of-bounds grid cell.  uv=(%.2f, %.2f) sm=%d/%d",
        new_loc->u, new_loc->v, new_idx, new_wall->grid->n_tiles);
  }

  struct surface_molecule_list *sm_list = new_wall->grid->sm_list[new_idx];
  if ((state->periodic_box_obj && periodicbox_in_surfmol_list(sm->periodic_box, sm_list)) ||
      (!state->periodic_box_obj && sm_list && sm_list->sm)) {
    if (hd_info != NULL) {
      delete_void_list((struct void_list *)hd_info);
      hd_info = NULL;
    }
    return 1; /* Pick again--full here */
  }

  count_moved_surface_mol(
    state, sm, new_wall->grid, new_loc, state->count_hashmask,
    state->count_hash, &state->ray_polygon_colls, previous_box);

  remove_surfmol_from_list(&sm->grid->sm_list[sm->grid_index], sm);
  sm->grid->n_occupied--;
  sm->grid = new_wall->grid;
  sm->grid_index = new_idx;
  sm_list = add_surfmol_with_unique_pb_to_list(sm->grid->sm_list[new_idx], sm);
  assert(sm_list != NULL);
  sm->grid->sm_list[sm->grid_index] = sm_list;
  sm->grid->n_occupied++;

  sm->s_pos.u = new_loc->u;
  sm->s_pos.v = new_loc->v;

  return 0;
}

/*************************************************************************
diffuse_2D:
  In: world: simulation state
      sm: molecule that is moving
      max_time: maximum time we can spend diffusing
      advance_time: how much to advance molecule internal time (return value)
  Out: Pointer to the molecule, or NULL if there was an error (right now
       there is no reallocation)
       Position and time are updated, but molecule is not rescheduled,
       nor does it react
  To-do: This doesn't work with triggers.  Change style of counting code
         so that it can update as we go, like with 3D diffusion.
*************************************************************************/
struct surface_molecule *diffuse_2D(
    struct volume *world,
    struct surface_molecule *sm,
    double max_time,
    double *advance_time) {

  struct species *spec = sm->properties;
  if (spec == NULL) {
    mcell_internal_error(
        "Attempted to take a 2-D diffusion step for a defunct molecule.");
  }

  // XXX: When would this ever happen, and shouldn't it just be an error?
  if (spec->space_step <= 0.0) {
    sm->t += max_time;
    return sm;
  }

  // Using global SPACE_STEP or per species CUSTOM_SPACE_STEP/CUSTOM_TIME_STEP
  if (spec->time_step > 1.0) {
    double sched_time = convert_iterations_to_seconds(
        world->start_iterations, world->time_unit,
        world->simulation_start_seconds, sm->t);
    /* Newly created particles with long custom time steps gradually increase
     * their timestep to the full value... Why??? */
    double f = 1 + 0.2 * ((sched_time - sm->birthday)/world->time_unit);
    if (f < 1)
      mcell_internal_error("A %s molecule is scheduled to move before it was "
                           "born [birthday=%.15g, t=%.15g]",
                           spec->sym->name, sm->birthday, sched_time);
    if (max_time > f)
      max_time = f;
  }

  double steps = 0.0;
  double t_steps = 0.0;
  double space_factor = 0.0;
  /* Where are we going? */
  if (spec->time_step > max_time) {
    t_steps = max_time;
    steps = max_time / spec->time_step;
  } else {
    t_steps = spec->time_step;
    steps = 1.0;
  }
  if (steps < EPS_C) {
    steps = EPS_C;
    t_steps = EPS_C * spec->time_step;
  }

  if (steps == 1.0) {
    space_factor = spec->space_step;
  } else {
    space_factor = spec->space_step * sqrt(steps);
  }

  world->diffusion_number++;
  world->diffusion_cumtime += steps;

  struct periodic_image previous_box = { .x = sm->periodic_box->x,
                                         .y = sm->periodic_box->y,
                                         .z = sm->periodic_box->z
                                       };
  struct hit_data *hd_info = NULL;
  for (int find_new_position = (SURFACE_DIFFUSION_RETRIES + 1);
       find_new_position > 0; find_new_position--) {
    hd_info = NULL;

    struct vector2 displacement;
    pick_2D_displacement(&displacement, space_factor, world->rng);

    if (sm->properties->flags & SET_MAX_STEP_LENGTH) {
      double disp_length = sqrt(displacement.u * displacement.u +
                         displacement.v * displacement.v);
      if (disp_length > sm->properties->max_step_length) {
        /* rescale displacement to the level of MAXIMUM_STEP_LENGTH */
        displacement.u *= (sm->properties->max_step_length / disp_length);
        displacement.v *= (sm->properties->max_step_length / disp_length);
      }
    }

    struct vector2 new_loc;
    struct rxn *rxp = NULL;
    int kill_me = 0;
    struct wall *new_wall = ray_trace_2D(world, sm, &displacement, &new_loc,
      &kill_me, &rxp, &hd_info);
    // Either something ambiguous happened or we hit absorptive border
    if (new_wall == NULL) {
      if (kill_me == 1) {
        /* molecule hit ABSORPTIVE region border */
        if (rxp == NULL) {
          mcell_internal_error("Error in 'ray_trace_2D()' after hitting "
                               "ABSORPTIVE region border.");
        }
        if (hd_info != NULL) {
          count_region_border_update(world, sm->properties, hd_info);
        }
        int result = outcome_unimolecular(world, rxp, 0,
                                      (struct abstract_molecule *)sm, sm->t);
        if (result != RX_DESTROY) {
          mcell_internal_error("Molecule should disappear after hitting "
                               "ABSORPTIVE region border.");
        }
        delete_void_list((struct void_list *)hd_info);
        hd_info = NULL;
        return NULL;
      }

      if (hd_info != NULL) {
        delete_void_list((struct void_list *)hd_info);
        hd_info = NULL;
      }
      continue; /* Something went wrong--try again */
    }

    // After diffusing, we are still on the SAME triangle.
    if (new_wall == sm->grid->surface) {
      if (move_sm_on_same_triangle(world, sm, &new_loc, &previous_box, new_wall, hd_info)) {
        continue; 
      }
    }
    // After diffusing, we ended up on a NEW triangle.
    else {
      if (move_sm_to_new_triangle(world, sm, &new_loc, &previous_box, new_wall, hd_info)) {
        continue;
      }
    }
    find_new_position = 0;
  }

  if (hd_info != NULL) {
    count_region_border_update(world, sm->properties, hd_info);
    delete_void_list((struct void_list *)hd_info);
    hd_info = NULL;
  }

  *advance_time = t_steps;
  return sm;
}

/***************************************************************************
react_2D_all_neighbors:
  In: world: simulation state
      sm: molecule that may react
      t: maximum duration we have to react
      molecule_collision_report:
      grid_grid_reaction_flag:
      surf_surf_colls:
  Out: Pointer to the molecule if it still exists (may have been
       destroyed), NULL otherwise.
  Note: Time is not updated--assume that's already taken care of
        elsewhere.
        This function takes into account variable number of neighbors.
  Note: If surface molecule (reaction initiator) or potential reaction
        partner are located on the different regions and any of them is
        behind the restrictive region boundary - we do not even
        test for reaction.
****************************************************************************/
struct surface_molecule *
react_2D_all_neighbors(struct volume *world, struct surface_molecule *sm,
                       double t, enum notify_level_t molecule_collision_report,
                       int grid_grid_reaction_flag,
                       long long *surf_surf_colls) {

  int i;     /* points to the pathway of the reaction */
  int j;     /* points to the the reaction */
  int n = 0; /* total number of possible reactions for a given molecules
                with all its neighbors */

  int l = 0;
  int num_matching_rxns = 0;
  struct rxn *matching_rxns[MAX_MATCHING_RXNS];

  /* linked list of the tile neighbors */
  struct tile_neighbor *tile_nbr_head = NULL, *curr;
  int list_length = 0; /* length of the linked lists above */

  if ((u_int)sm->grid_index >= sm->grid->n_tiles) {
    mcell_internal_error("tile index %u is greater or equal number_of_tiles %u",
                         (u_int)sm->grid_index, sm->grid->n_tiles);
  }

  find_neighbor_tiles(world, sm, sm->grid, sm->grid_index, 0, 1, &tile_nbr_head,
                      &list_length);

  if (tile_nbr_head == NULL)
    return sm; /* no reaction may happen */

  const int num_nbrs = list_length;
  int max_size = num_nbrs * MAX_MATCHING_RXNS;
  struct rxn *rxn_array[max_size]; /* array of reaction objects with neighbor
                                     molecules */
  double local_prob_factor; /* local probability factor for the
                                 reactions */
  double cf[max_size]; /* Correction factors for area for those molecules */

  struct surface_molecule *smol[max_size]; /* points to neighbor molecules */

  /* Calculate local_prob_factor for the reaction probability.
     Here we convert from 3 neighbor tiles (upper probability
     limit) to the real "num_nbrs" neighbor tiles. */

  local_prob_factor = 3.0 / num_nbrs;

  for (int kk = 0; kk < max_size; kk++) {
    rxn_array[kk] = NULL;
    smol[kk] = NULL;
    cf[kk] = 0;
  }

  /* step through the neighbors */
  for (curr = tile_nbr_head; curr != NULL; curr = curr->next) {
    /* Neighboring molecule */
    struct surface_molecule_list *sm_list = curr->grid->sm_list[curr->idx]; 
    if (sm_list == NULL || sm_list->sm == NULL)
      continue;
    struct surface_molecule *smp = curr->grid->sm_list[curr->idx]->sm;

    /* check whether the neighbor molecule is behind
       the restrictive region boundary   */
    if ((sm->properties->flags & CAN_REGION_BORDER) ||
        (smp->properties->flags & CAN_REGION_BORDER)) {
      if (sm->grid->surface != smp->grid->surface) {
        /* INSIDE-OUT check */
        if (walls_belong_to_at_least_one_different_restricted_region(
                world, sm->grid->surface, sm, smp->grid->surface, smp))
          continue;

        /* OUTSIDE-IN check */
        if (walls_belong_to_at_least_one_different_restricted_region(
                world, sm->grid->surface, smp, smp->grid->surface, sm))
          continue;
      }
    }

    num_matching_rxns = trigger_bimolecular(
        world->reaction_hash, world->rx_hashsize, sm->properties->hashval,
        smp->properties->hashval, (struct abstract_molecule *)sm,
        (struct abstract_molecule *)smp, sm->orient, smp->orient,
        matching_rxns);

    if (num_matching_rxns > 0) {
      if (molecule_collision_report == NOTIFY_FULL) {
        if (grid_grid_reaction_flag)
          surf_surf_colls++;
      }

      for (int jj = 0; jj < num_matching_rxns; jj++) {
        if (matching_rxns[jj] != NULL) {
          if (matching_rxns[jj]->prob_t != NULL)
            update_probs(world, matching_rxns[jj], sm->t);
          rxn_array[l] = matching_rxns[jj];
          cf[l] = t / (curr->grid->binding_factor);
          smol[l] = smp;
          l++;
        }
      }

      n += num_matching_rxns;
    }
  }

  delete_tile_neighbor_list(tile_nbr_head);

  if (n == 0) {
    return sm; /* Nobody to react with */
  } else if (n == 1) {
    i = test_bimolecular(rxn_array[0], cf[0], local_prob_factor, NULL, NULL,
                         world->rng);
    j = 0;
  } else {
    // previously "test_many_bimolecular_all_neighbors"
    int all_neighbors_flag = 1;
    j = test_many_bimolecular(rxn_array, cf, local_prob_factor, n, &(i), 
                              world->rng, all_neighbors_flag);
  }

  if ((j == RX_NO_RX) || (i < RX_LEAST_VALID_PATHWAY)) {
    return sm; /* No reaction */
  }

  /* run the reaction */
  int outcome_bimol_result = outcome_bimolecular(
      world, rxn_array[j], i, (struct abstract_molecule *)sm,
      (struct abstract_molecule *)smol[j], sm->orient, smol[j]->orient, sm->t,
      NULL, NULL);

  if (outcome_bimol_result == RX_DESTROY) {
    mem_put(sm->birthplace, sm);
    return NULL;
  }

  return sm;
}

/*************************************************************************
clean_up_old_molecules:

 This function just removes defunct molecules from the scheduler.
*************************************************************************/
void clean_up_old_molecules(struct storage *local) {
  if (local->timer->defunct_count > MIN_DEFUNCT_FOR_GC &&
      MAX_DEFUNCT_FRAC * (local->timer->count) < local->timer->defunct_count) {
    struct abstract_molecule *am;
    am = (struct abstract_molecule *)schedule_cleanup(local->timer,
                                                      *is_defunct_molecule);
    while (am != NULL) {
      struct abstract_molecule *temp = am;
      am = am->next;
      if ((temp->flags & IN_MASK) == IN_SCHEDULE) {
        temp->next = NULL;
        mem_put(temp->birthplace, temp);
      } else {
        temp->flags &= ~IN_SCHEDULE;
      }
    }
  }
}

/*************************************************************************
reschedule_surface_molecules:

 If a surface molecules moves across a memory subdivision boundary, it might
 need to be reallocated and moved to a new scheduler.
*************************************************************************/
void reschedule_surface_molecules(
    struct volume *state, struct storage *local,
    struct abstract_molecule *am) {
  struct vector3 pos3d;
  struct surface_molecule *sm = (struct surface_molecule *)(void *)am;
  uv2xyz(&sm->s_pos, sm->grid->surface, &pos3d);

  struct subvolume *sv = find_subvolume(state, &pos3d, sm->grid->subvol);
  if (sv->local_storage != local) {
    struct surface_molecule *sm_new =
        (struct surface_molecule *)CHECKED_MEM_GET(sv->local_storage->smol,
                                                   "surface molecule");
    memcpy(sm_new, sm, sizeof(struct surface_molecule));
    sm_new->next = NULL;
    sm_new->birthplace = sv->local_storage->smol;
    if (sm->grid->sm_list[sm->grid_index] && 
        (sm->grid->sm_list[sm->grid_index]->sm == sm)) {
      sm->grid->sm_list[sm->grid_index]->sm = sm_new;
      sm->grid = NULL;
      sm->grid_index = 0;
    }

    mem_put(sm->birthplace, sm);
    if (schedule_add(sv->local_storage->timer, sm_new))
      mcell_allocfailed("Failed to add a '%s' surface molecule to scheduler "
                        "after migrating to a new memory store.",
                        am->properties->sym->name);
  } else {
    if (schedule_add(local->timer, am))
      mcell_allocfailed("Failed to add a '%s' surface molecule to scheduler "
                        "after taking a diffusion step.",
                        am->properties->sym->name);
  }
}

/*************************************************************************
run_timestep:
  In: state: simulation state
      local: local storage area to use
      release_time: time of the next release event
      checkpt_time: time of the next checkpoint
  Out: No return value.  Every molecule in the subvolume is updated in
       position and rescheduled at least one timestep ahead.
  Note: This also occasionally does garbage collection on the scheduling
        queue.
*************************************************************************/
void run_timestep(struct volume *state, struct storage *local,
                  double release_time, double checkpt_time) {
  struct abstract_molecule *am;

  // Check for garbage collection first
  clean_up_old_molecules(local);

  // Now run the timestep

  /* Do not trigger the scheduler to advance!  This will be done
   * by the main loop. */
  while (local->timer->current != NULL) {
    am = (struct abstract_molecule *)schedule_next(local->timer);
    if (am->properties == NULL) /* Defunct!  Remove molecule. */
    {
      if ((am->flags & IN_MASK) == IN_SCHEDULE) {
        am->next = NULL;
        mem_put(am->birthplace, am);
      } else
        am->flags &= ~IN_SCHEDULE;
      if (local->timer->defunct_count > 0)
        local->timer->defunct_count--;

      continue;
    }

    am->flags &= ~IN_SCHEDULE;

    // Check for unimolecular reactions
    // If molec is new or need rescheduled, this just computes a new lifetime
    if (am->t2 < EPS_C || am->t2 < EPS_C * am->t) {
      if (!check_for_unimolecular_reaction(state, am)) {
        continue;
      }
    }

    // How to advance surface molecule scheduling time
    double surface_mol_advance_time = 0;

    struct wall *current_wall = NULL;
    // The maximum time we can spend diffusing or looking for reactions
    double max_time;
    int can_diffuse = ((am->flags & ACT_DIFFUSE) != 0);
    if (can_diffuse) {
      max_time = checkpt_time - am->t;
      if (local->max_timestep < max_time)
        max_time = local->max_timestep;
      if ((am->flags & (ACT_REACT)) != 0 && am->t2 < max_time)
        max_time = am->t2;

      if ((am->flags & TYPE_VOL) != 0) {
        double save_sched_time = am->t;
        if (max_time > release_time - am->t)
          max_time = release_time - am->t;
        if (am->properties->flags & (CAN_VOLVOLVOL | CAN_VOLVOLSURF))
          am = (struct abstract_molecule *)diffuse_3D_big_list(
              state, (struct volume_molecule *)am, max_time);
        else
          am = (struct abstract_molecule *)diffuse_3D(
              state, (struct volume_molecule *)am, max_time);
        if (am != NULL) /* We still exist */
        {
          // Perform only for unimolecular reactions
          if ((am->flags & ACT_REACT) != 0) {
            am->t2 -= am->t - save_sched_time;
            if (am->t2 < 0)
              am->t2 = 0;
          }
        } else
          continue;
      } else {
        if (max_time > release_time - am->t) {
          max_time = release_time - am->t;
        }

        // Remember current wall
        current_wall = ((struct surface_molecule *)am)->grid->surface;

        am = (struct abstract_molecule *)diffuse_2D(
            state, (struct surface_molecule *)am, max_time,
            &surface_mol_advance_time);
        if (am == NULL) {
          continue;
        }
      }
    }

    int can_surface_mol_react =
        (am->properties->flags & (CAN_SURFSURFSURF | CAN_SURFSURF));
    if (((am->flags & TYPE_SURF) != 0) && can_surface_mol_react) {
      // Didn't move, so we need to figure out how long to react for
      if (!can_diffuse) 
      {
        max_time = checkpt_time - am->t;
        if (am->t2 < max_time && (am->flags & (ACT_REACT)) != 0)
          max_time = am->t2;
        if (max_time > release_time - am->t)
          max_time = release_time - am->t;
        if (am->properties->time_step < max_time)
          max_time = am->properties->time_step;
        surface_mol_advance_time = max_time;
      } else
        max_time = surface_mol_advance_time;

      if (can_surface_mol_react) {
        if ((am->properties->flags & (CANT_INITIATE | CAN_SURFSURF)) ==
            CAN_SURFSURF) {
          am = (struct abstract_molecule *)react_2D_all_neighbors(
              state, (struct surface_molecule *)am, max_time,
              state->notify->molecule_collision_report,
              state->rxn_flags.surf_surf_reaction_flag,
              &(state->surf_surf_colls));
          if (am == NULL)
            continue;
        }
        if ((am->properties->flags & (CANT_INITIATE | CAN_SURFSURFSURF)) ==
            CAN_SURFSURFSURF) {
          am = (struct abstract_molecule *)react_2D_trimol_all_neighbors(
              state, (struct surface_molecule *)am, max_time,
              state->notify->molecule_collision_report,
              state->notify->final_summary,
              state->rxn_flags.surf_surf_surf_reaction_flag,
              &(state->surf_surf_surf_colls));
          if (am == NULL)
            continue;
        }
      }
    }

    // Advance surface molecule scheduling time
    if ((am->flags & TYPE_SURF) != 0 && (can_diffuse || can_surface_mol_react)) {
      am->t += surface_mol_advance_time;

      // Perform only for unimolecular reactions
      if ((am->flags & ACT_REACT) != 0) {
        /* This case takes care of newly created surface products A which
         * only have a unimolecular surface reaction defined (A @surf) and
         * are thus scheduled am->t2 = FOREVER */
        int can_surf_react = ((am->properties->flags & CAN_SURFWALL) != 0);
        if (can_surf_react && !distinguishable(am->t2, (double)FOREVER, EPS_C)) {
          am->t2 = 0;
          am->flags |= ACT_CHANGE; /* Reschedule reaction time */
        }
        else {
          am->t2 -= surface_mol_advance_time;
          if (am->t2 < 0) {
            am->t2 = 0;
          }
          /* If the molecule didn't leave its wall AND its lifetime hasn't run
           * out, then we don't need to reschedule.
           * NOTE: We really only have to make sure that the new wall has the
           * same collection of surface classes. Doing so could be
           * significantly more efficient. */
          if ((current_wall !=
              ((struct surface_molecule *)am)->grid->surface) &&
              (am->t2 > EPS_C || am->t2 > EPS_C * am->t)) {
            am->t2 = 0;
            am->flags |= ACT_CHANGE; /* Reschedule reaction time */
          }
        }
      }
    } else if (!can_diffuse) {
      // NOTE: t2 should only be 0 at this point if "am" is inert. This is
      // basically just a clunky way to ignore it, although it's not clear why
      // we can't just set t to the total number of iterations.
      if (am->t2 == 0)
        am->t += MAX_UNI_TIMESKIP;
      else {
        am->t += am->t2;
        am->t2 = 0;
      }
    }

    am->flags |= IN_SCHEDULE;

    /* If we're near an integer boundary, advance to the next integer */
    double t = ceil(am->t) * (1.0 + 0.1 * EPS_C);
    if (!distinguishable(t, am->t, EPS_C))
      am->t = t;

    if (am->flags & TYPE_SURF) {
      reschedule_surface_molecules(state, local, am);
    } else {
      if (schedule_add(
              ((struct volume_molecule *)am)->subvol->local_storage->timer, am))
        mcell_allocfailed("Failed to add a '%s' volume molecule to scheduler "
                          "after taking a diffusion step.",
                          am->properties->sym->name);
    }
  }
  if (local->timer->error)
    mcell_internal_error("Scheduler reported an out-of-memory error while "
                         "retrieving molecules, but this should never happen.");
}


/*************************************************************************
run_concentration_clamp:
  In: world: simulation state
      t_now: the current time.
  Out: No return value.  Molecules are released at concentration-clamped
       surfaces to maintain the desired concentation.
*************************************************************************/
void run_concentration_clamp(struct volume *world, double t_now) {
  int this_count = 0;
  static int total_count = 0;
  for (struct ccn_clamp_data *ccd = world->clamp_list; ccd != NULL; ccd = ccd->next) {
    if (ccd->objp == NULL) {
      continue;
    }
    for (struct ccn_clamp_data *ccdo = ccd; ccdo != NULL; ccdo = ccdo->next_obj) {
      for (struct ccn_clamp_data *ccdm = ccdo; ccdm != NULL; ccdm = ccdm->next_mol) {
        double n_collisions = ccdo->scaling_factor * ccdm->mol->space_step *
                       ccdm->concentration / ccdm->mol->time_step;
        if (ccdm->orient != 0) {
          n_collisions *= 0.5;
        }
        int n_emitted = poisson_dist(n_collisions, rng_dbl(world->rng));

        if (n_emitted == 0)
          continue;

        struct volume_molecule vm;
        vm.t = t_now + 0.5;
        vm.t2 = 0;
        vm.flags = IN_SCHEDULE | ACT_NEWBIE | TYPE_VOL | IN_VOLUME |
                  ACT_CLAMPED | ACT_DIFFUSE;
        vm.properties = ccdm->mol;
        vm.mesh_name = NULL;
        vm.birthplace = NULL;
        vm.birthday = convert_iterations_to_seconds(
            world->start_iterations, world->time_unit,
            world->simulation_start_seconds, t_now);
        vm.subvol = NULL;
        vm.previous_wall = NULL;
        vm.index = 0;
        struct volume_molecule *vmp = NULL;

        this_count += n_emitted;
        while (n_emitted > 0) {
          int idx = bisect_high(ccdo->cum_area, ccdo->n_sides,
                            rng_dbl(world->rng) *
                                ccdo->cum_area[ccd->n_sides - 1]);
          struct wall *w = ccdo->objp->wall_p[ccdo->side_idx[idx]];

          double s1 = sqrt(rng_dbl(world->rng));
          double s2 = rng_dbl(world->rng) * s1;

          struct vector3 v;
          v.x = w->vert[0]->x + s1 * (w->vert[1]->x - w->vert[0]->x) +
                s2 * (w->vert[2]->x - w->vert[1]->x);
          v.y = w->vert[0]->y + s1 * (w->vert[1]->y - w->vert[0]->y) +
                s2 * (w->vert[2]->y - w->vert[1]->y);
          v.z = w->vert[0]->z + s1 * (w->vert[1]->z - w->vert[0]->z) +
                s2 * (w->vert[2]->z - w->vert[1]->z);

          if (ccdm->orient == 1) {
            vm.index = 1;
          }
          else if (ccdm->orient == -1) {
            vm.index = -1;
          }
          else {
            vm.index = (rng_uint(world->rng) & 2) - 1;
          }

          double eps = EPS_C * vm.index;

          s1 = fabs(v.x);
          s2 = fabs(v.y);
          if (s1 < s2) {
            s1 = s2;
          }
          s2 = fabs(v.z);
          if (s1 < s2) {
            s1 = s2;
          }
          if (s1 > 1.0){
            eps *= s1;
          }

          vm.pos.x = v.x + w->normal.x * eps;
          vm.pos.y = v.y + w->normal.y * eps;
          vm.pos.z = v.z + w->normal.z * eps;
          vm.previous_wall = w;
          // TODO: This isn't right. We need to figure out what PB these should
          // really be created in.
          struct periodic_image periodic_box = {.x = 0,
                                                .y = 0,
                                                .z = 0
                                               };
          
          vm.periodic_box = &periodic_box;

          if (vmp == NULL) {
            vmp = insert_volume_molecule(world, &vm, vmp);
            if (vmp == NULL)
              mcell_allocfailed("Failed to insert a '%s' volume molecule while "
                                "concentration clamping.",
                                vm.properties->sym->name);
            if (trigger_unimolecular(world->reaction_hash, world->rx_hashsize,
                                     ccdm->mol->hashval,
                                     (struct abstract_molecule *)vmp) != NULL) {
              vm.flags |= ACT_REACT;
              vmp->flags |= ACT_REACT;
            }
          } else {
            vmp = insert_volume_molecule(world, &vm, vmp);
            if (vmp == NULL)
              mcell_allocfailed("Failed to insert a '%s' volume molecule while "
                                "concentration clamping.",
                                vm.properties->sym->name);
          }

          n_emitted--;
        }
      }
    }
  }

  total_count += this_count;
}


/******************************************************************************
 *
 * redo_collision list is a helper function used in diffuse_3D to compute the
 * list of possible collisions in neighboring subvolumes.
 *
 ******************************************************************************/
void redo_collision_list(struct volume* world, struct collision** shead,
  struct collision** stail, struct collision** shead_exp, struct volume_molecule* m,
  struct vector3* displacement, struct subvolume* sv) {

  struct collision* st = *stail;
  struct collision* sh = *shead_exp;
  if (st != NULL) {
    st->next = NULL;
    if (sh != NULL) {
      mem_put_list(sv->local_storage->coll, sh);
      sh = NULL;
    }
  } else if (sh != NULL) {
    mem_put_list(sv->local_storage->coll, sh);
    sh = NULL;
    *shead = NULL;
  }
  if ((m->properties->flags & (CAN_VOLVOL | CANT_INITIATE)) == CAN_VOLVOL) {
    sh = expand_collision_list(m, displacement, sv, world->rx_radius_3d,
      world->ny_parts, world->nz_parts, world->x_fineparts,
      world->y_fineparts, world->z_fineparts, world->rx_hashsize,
      world->reaction_hash);
    if (st != NULL)
      st->next = sh;
    else {
      if (*shead != NULL)
        mcell_internal_error("Collision lists corrupted.  While expanding "
                             "the collision lists, expected shead to be "
                             "NULL, but it wasn't.");
      *shead = sh;
    }
  }
  *stail = st;
  *shead_exp = sh;
}

/******************************************************************************
 *
 * collide_and_react_with_vol_mol is a helper function used in diffuse_3D to
 * handle collision of a diffusing molecule with a molecular target.
 *
 * Returns 1 if reaction does happen and 0 otherwise.
 *
 ******************************************************************************/
static int collide_and_react_with_vol_mol(struct volume* world,
  struct collision* smash, struct volume_molecule* m, struct collision** tentative,
  struct vector3* displacement, struct vector3* loc_certain, double t_steps,
  double r_rate_factor) {

  struct abstract_molecule* am = (struct abstract_molecule *)smash->target;
  double factor = exact_disk(
      world, &(smash->loc), displacement, world->rx_radius_3d, m->subvol,
      m, (struct volume_molecule *)am, world->use_expanded_list,
      world->x_fineparts, world->y_fineparts, world->z_fineparts);

  if (factor < 0) { /* Probably hit a wall, might have run out of memory */
    return 0; /* Reaction blocked by a wall */
  }

  double scaling = factor * r_rate_factor;
  struct rxn* rx = smash->intermediate;
  if ((rx != NULL) && (rx->prob_t != NULL)) {
    update_probs(world, rx, m->t);
  }

  struct species *spec = m->properties;
  struct periodic_image *periodic_box = m->periodic_box;
  int i = test_bimolecular(
    rx, scaling, 0, am, (struct abstract_molecule *)m, world->rng);

  if (i < RX_LEAST_VALID_PATHWAY) {
    return 0;
  }

  int j = outcome_bimolecular(world, rx, i, (struct abstract_molecule *)m, am,
    0, 0, m->t + t_steps * smash->t, &(smash->loc), loc_certain);

  if (j != RX_DESTROY) {
    return 0;
  } else {
    /* Count the hits up until we were destroyed */
    struct collision* ttv = *tentative;
    for (; ttv != NULL && ttv->t <= smash->t; ttv = ttv->next) {
      if (!(ttv->what & COLLIDE_WALL)) {
        continue;
      }
      if (m->properties == NULL) {
        continue;
      }
      if (!(m->properties->flags & ((struct wall *)ttv->target)->flags &
            COUNT_SOME_MASK)) {
        continue;
      }
      count_region_update(world, spec, periodic_box,
        ((struct wall *)ttv->target)->counting_regions,
        ((ttv->what & COLLIDE_MASK) == COLLIDE_FRONT) ? 1 : -1, 0, &(ttv->loc), ttv->t);
      if (ttv == smash) {
        break;
      }
    }
    *tentative = ttv;
  }
  return 1;
}

/******************************************************************************
 *
 * collide_and_react_with_surf_mol is a helper function used in diffuse_3D to
 * handle collision of a diffusing 3D molecule with a surface molecule
 *
 * Return values:
 *
 * -1 : nothing happened - continue on with next smash targets
 *  0 : reaction happened and we still exist but are done with the current smash
 *  1 : reaction happened and we are destroyed
 *      target
 *
 ******************************************************************************/
int collide_and_react_with_surf_mol(struct volume* world, struct collision* smash,
  struct volume_molecule* m, struct collision** tentative,
  struct vector3** loc_certain, double t_steps, int mol_grid_flag,
  int mol_grid_grid_flag, double r_rate_factor) {

  struct collision* ttv = *tentative;
  struct vector3* loc = *loc_certain;
  struct wall* w = (struct wall *)smash->target;

  double t_confident = 0.0;
  if (smash->next == NULL) {
    t_confident = smash->t;
  } else if (smash->next->t * (1.0 - EPS_C) > smash->t) {
    t_confident = smash->t;
  } else {
    t_confident = smash->t * (1.0 - EPS_C);
  }

  int k = -1;
  if ((smash->what & COLLIDE_MASK) == COLLIDE_FRONT) {
    k = 1;
  }

  int j = xyz2grid(&(smash->loc), w->grid);
  struct surface_molecule_list *sm_list = w->grid->sm_list[j]; 
  if (sm_list == NULL || sm_list->sm == NULL) {
    return -1;
  }
  struct surface_molecule* sm = w->grid->sm_list[j]->sm;
  if (m->index == j && m->previous_wall == w) {
    m->index = -1; // Avoided rebinding, but next time it's OK
    return -1;
  }

  int num_matching_rxns = 0;
  struct rxn *matching_rxns[MAX_MATCHING_RXNS];
  double scaling_coef[MAX_MATCHING_RXNS];
  struct species* spec = m->properties;
  struct periodic_image *periodic_box = m->periodic_box;
  int ii = 0, jj = 0;
  if (mol_grid_flag) {
    num_matching_rxns = trigger_bimolecular(
      world->reaction_hash, world->rx_hashsize, spec->hashval,
      sm->properties->hashval, (struct abstract_molecule *)m,
      (struct abstract_molecule *)sm, k, sm->orient, matching_rxns);
    if (num_matching_rxns > 0) {
      if (world->notify->molecule_collision_report == NOTIFY_FULL) {
        if (world->rxn_flags.vol_surf_reaction_flag)
          world->vol_surf_colls++;
      }

      for (int l = 0; l < num_matching_rxns; l++) {
        if (matching_rxns[l]->prob_t != NULL)
          update_probs(world, matching_rxns[l], m->t);
        scaling_coef[l] = r_rate_factor / w->grid->binding_factor;
      }

      if (num_matching_rxns == 1) {
        ii = test_bimolecular(matching_rxns[0], scaling_coef[0], 0,
          (struct abstract_molecule *)m, (struct abstract_molecule *)sm,
          world->rng);
        jj = 0;
      } else {
        jj = test_many_bimolecular(matching_rxns, scaling_coef, 0,
          num_matching_rxns, &(ii), world->rng, 0);
      }
      if ((jj > RX_NO_RX) && (ii >= RX_LEAST_VALID_PATHWAY)) {
        /* Save m flags in case m gets collected in outcome_bimolecular */
        short mflags = m->flags;
        int l = outcome_bimolecular(world, matching_rxns[jj], ii,
          (struct abstract_molecule *)m, (struct abstract_molecule *)sm,
          k, sm->orient, m->t + t_steps * smash->t, &(smash->loc), loc);

        if (l == RX_FLIP) {
          if ((m->flags & COUNT_ME) != 0 && (spec->flags & COUNT_SOME_MASK) != 0) {
            /* Count as far up as we can unambiguously */
            int destroy_flag = 0;
            count_tentative_collisions(
              world, &ttv, smash, spec, t_confident, destroy_flag,
              periodic_box);
          }
          *tentative = ttv;
          *loc_certain = &(ttv->loc);
          return 0; /* pass through */
        } else if (l == RX_DESTROY) {
          if ((mflags & COUNT_ME) != 0 && (spec->flags & COUNT_HITS) != 0) {
            /* Count the hits up until we were destroyed */
            int destroy_flag = 0;
            count_tentative_collisions(
              world, &ttv, smash, spec, t_confident, destroy_flag,
              periodic_box);
          }
          *tentative = ttv;
          return 1;
        }
      }
    }
  }

  /* test for the trimolecular reactions of the type MOL_GRID_GRID */
  if (mol_grid_grid_flag) {
    struct surface_molecule *smp; /* Neighboring molecules */
    struct tile_neighbor *tile_nbr_head = NULL, *curr;
    int list_length = 0;
    int n = 0; /* total number of possible reactions for a given
                   molecule with all its neighbors */

    /* find neighbor molecules to react with */
    find_neighbor_tiles(world, sm, sm->grid, sm->grid_index, 0, 1,
      &tile_nbr_head, &list_length);
    if (tile_nbr_head != NULL) {
      const int num_nbrs = (const int)list_length;
      double local_prob_factor; /*local probability factor for the
                                   reaction */
      int max_size = num_nbrs * MAX_MATCHING_RXNS;
      /* array of reaction objects with neighbor mols */
      struct rxn *rxn_array[max_size];
      /* correction factors for areas for these mols */
      double cf[max_size];
      struct surface_molecule *smol[max_size]; /* points to neighbor mols */

      local_prob_factor = 3.0 / num_nbrs;
      jj = RX_NO_RX;
      ii = RX_LEAST_VALID_PATHWAY - 1;
      for (int kk = 0; kk < max_size; kk++) {
        smol[kk] = NULL;
        rxn_array[kk] = NULL;
        cf[kk] = 0;
      }

      /* step through the neighbors */
      int ll = 0;
      for (curr = tile_nbr_head; curr != NULL; curr = curr->next) {
        sm_list = curr->grid->sm_list[curr->idx];
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
                    world, sm->grid->surface, sm, smp->grid->surface, smp)) {
              continue;
            }

            /* OUTSIDE-IN check */
            if (walls_belong_to_at_least_one_different_restricted_region(
                    world, sm->grid->surface, smp, smp->grid->surface, sm)) {
              continue;
            }
          }
        }

        num_matching_rxns = trigger_trimolecular(world->reaction_hash,
          world->rx_hashsize, spec->hashval, sm->properties->hashval,
          smp->properties->hashval, spec, sm->properties, smp->properties,
          k, sm->orient, smp->orient, matching_rxns);

        if (num_matching_rxns > 0) {
          if (world->notify->molecule_collision_report == NOTIFY_FULL &&
              world->rxn_flags.vol_surf_surf_reaction_flag) {
              world->vol_surf_surf_colls++;
          }
          for (j = 0; j < num_matching_rxns; j++) {
            if (matching_rxns[j]->prob_t != NULL) {
              update_probs(world, matching_rxns[j], m->t);
            }
            rxn_array[ll] = matching_rxns[j];
            cf[ll] = r_rate_factor / (w->grid->binding_factor *
                                      curr->grid->binding_factor);
            smol[ll] = smp;
            ll++;
          }
          n += num_matching_rxns;
        }
      }
      delete_tile_neighbor_list(tile_nbr_head);

      if (n == 1) {
        ii = test_bimolecular(rxn_array[0], cf[0], local_prob_factor,
          NULL, NULL, world->rng);
        jj = 0;
      } else if (n > 1) {
        // previously "test_many_bimolecular_all_neighbors"
        int all_neighbors_flag = 1;
        jj = test_many_bimolecular(rxn_array, cf, local_prob_factor,
          n, &(ii), world->rng, all_neighbors_flag);
      }

      if (n > max_size)
        mcell_internal_error(
            "The size of the reactions array is not sufficient.");

      if ((n > 0) && (ii >= RX_LEAST_VALID_PATHWAY) && (jj > RX_NO_RX)) {
        /* Save m flags in case it gets collected in outcome_trimolecular */
        int mflags = m->flags;
        int l = outcome_trimolecular(world, rxn_array[jj], ii,
          (struct abstract_molecule *)m, (struct abstract_molecule *)sm,
          (struct abstract_molecule *)smol[jj], k, sm->orient,
          smol[jj]->orient, m->t + t_steps * smash->t, &smash->loc, &m->pos);

        if (l == RX_FLIP) {
          if ((m->flags & COUNT_ME) != 0 && (spec->flags & COUNT_SOME_MASK) != 0) {
            /* Count as far up as we can unambiguously */
            int destroy_flag = 0;
            count_tentative_collisions(
              world, &ttv, smash, spec, t_confident, destroy_flag,
              periodic_box);
          }
          *loc_certain = &(ttv->loc);
          *tentative = ttv;
          return 0; /* pass through */
        } else if (l == RX_DESTROY) {
          if ((mflags & COUNT_ME) != 0 && (spec->flags & COUNT_HITS) != 0) {
            /* Count the hits up until we were destroyed */
            int destroy_flag = 0;
            count_tentative_collisions(
              world, &ttv, smash, spec, t_confident, destroy_flag,
              periodic_box);
          }
          *tentative = ttv;
          return 1;
        }
      }
    }
  }
  return -1;
}


/******************************************************************************
 *
 * check_collisions_with_walls is a helper function used in diffuse_3D to handle
 * collision of a diffusing molecule with a wall
 *
 * Return values:
 *
 * -1 : nothing happened - continue on with next smash targets
 *  0 : reaction happened and we still exist but are done with the current smash
 *      target
 *  1 : reaction happened and we are destroyed
 *
 ******************************************************************************/
int collide_and_react_with_walls(struct volume* world, struct collision* smash,
  struct volume_molecule* m, struct collision** tentative,
  struct vector3** loc_certain, double t_steps, int inertness,
  double r_rate_factor) {

  struct collision *ttv = *tentative;
  struct vector3 *loc = *loc_certain;

  double t_confident = 0.0;
  if (smash->next == NULL) {
    t_confident = smash->t;
  } else if (smash->next->t * (1.0 - EPS_C) > smash->t) {
    t_confident = smash->t;
  } else {
    t_confident = smash->t * (1.0 - EPS_C);
  }

  int k = -1;
  if ((smash->what & COLLIDE_MASK) == COLLIDE_FRONT) {
    k = 1;
  }

  // check for reactions with walls
  m->index = -1;
  int num_matching_rxns = 0;
  struct rxn* rx = NULL;
  struct species* spec = m->properties;
  struct wall* w = (struct wall *)smash->target;
  struct rxn *matching_rxns[MAX_MATCHING_RXNS];
  num_matching_rxns = trigger_intersect(world->reaction_hash,
    world->rx_hashsize, world->all_mols, world->all_volume_mols,
    world->all_surface_mols, spec->hashval, (struct abstract_molecule *)m, k, w,
    matching_rxns, 1, 0, 0);
  if (num_matching_rxns == 0) {
    return -1;
  }

  int is_transp_flag = 0;
  struct rxn *transp_rx = NULL;
  for (int ii = 0; ii < num_matching_rxns; ii++) {
    rx = matching_rxns[ii];
    if (rx->n_pathways == RX_TRANSP) {
      is_transp_flag = 1;
      transp_rx = matching_rxns[ii];
      break;
    }
  }

  if ((!is_transp_flag) && (world->notify->molecule_collision_report == NOTIFY_FULL) &&
       world->rxn_flags.vol_wall_reaction_flag) {
    world->vol_wall_colls++;
  }

  struct periodic_image *periodic_box = m->periodic_box;
  if (is_transp_flag) {
    transp_rx->n_occurred++;
    if ((m->flags & COUNT_ME) != 0 && (spec->flags & COUNT_SOME_MASK) != 0) {
      /* Count as far up as we can unambiguously */
      int destroy_flag = 0;
      count_tentative_collisions(
        world, tentative, smash, spec, smash->t, destroy_flag, periodic_box);
      // XXX: RX_FLIP case below should probably be handled this way too.
      for (; ttv != NULL && ttv->t <= t_confident; ttv = ttv->next) {
        *loc_certain = &(ttv->loc);
      }
    }
    return 0; /* Ignore this wall and keep going */
  } else if (inertness < inert_to_all) {
    /* Collisions with the surfaces declared REFLECTIVE are treated similar to
     * the default surfaces after this loop. */
    for (int l = 0; l < num_matching_rxns; l++) {
      if (matching_rxns[l]->prob_t != NULL) {
        update_probs(world, matching_rxns[l], m->t);
      }
    }
    int jj = 0;
    int i = 0;
    if (num_matching_rxns == 1) {
      i = test_intersect(matching_rxns[0], r_rate_factor, world->rng);
      jj = 0;
    } else {
      jj = test_many_intersect(matching_rxns, r_rate_factor,
                               num_matching_rxns, &(i), world->rng);
    }

    if ((i >= RX_LEAST_VALID_PATHWAY) && (jj > RX_NO_RX)) {
      /* Save m flags in case it gets collected in outcome_intersect */
      rx = matching_rxns[jj];
      int mflags = m->flags;
      int j = outcome_intersect(world, rx, i, w, (struct abstract_molecule *)m, k,
        m->t + t_steps * smash->t, &(smash->loc), loc);

      if (j == RX_FLIP) {
        if ((m->flags & COUNT_ME) != 0 && (spec->flags & COUNT_SOME_MASK) != 0) {
          /* Count as far up as we can unambiguously */
          int destroy_flag = 0;
          count_tentative_collisions(
            world, &ttv, smash, spec, t_confident, destroy_flag, periodic_box);
        }
        *loc_certain = &(ttv->loc);
        *tentative = ttv;
        return 0; /* pass through */
      } else if (j == RX_DESTROY) {
        if ((mflags & COUNT_ME) != 0 && (spec->flags & COUNT_HITS) != 0) {
          /* Count the hits up until we were destroyed */
          int destroy_flag = 1;
          count_tentative_collisions(
            world, tentative, smash, spec, smash->t, destroy_flag,
            periodic_box);
        }
        return 1;
      }
    }
  }
  return -1;
}

/******************************************************************************
 *
 * the reflect_or_periodic_bc helper function is used in diffuse_3D to handle
 * either reflections or periodic boundary conditions for a diffusing molecule
 * encountering a wall
 *
 * Return values:
 *
 *  0 : indicates that the molecule reflected off a wall
 *  1 : indicates that the molecule hit a periodic box and was moved to a
 *      position in the neighboring image
 *
 ******************************************************************************/
int reflect_or_periodic_bc(
    struct volume* world, struct collision* smash,
    struct vector3* displacement, struct volume_molecule** mol,
    struct wall** reflectee, struct collision** tentative, double* t_steps) {

  struct wall* w = (struct wall*)smash->target;
  struct wall *reflect_w = w;
  double reflect_t = smash->t;
  struct volume_molecule* vm = *mol;
  bool periodic_traditional = world->periodic_traditional;
  register_hits(world, vm, tentative, &reflect_w, &reflect_t, displacement,
    smash, t_steps);
  struct vector3 orig_pos = {vm->pos.x, vm->pos.y, vm->pos.z};
  (*reflectee) = reflect_w;

  int k = -1;
  if ((smash->what & COLLIDE_MASK) == COLLIDE_FRONT) {
    k = 1;
  }

  bool periodic_x = w->parent_object->periodic_x && (k == -1);
  bool periodic_y = w->parent_object->periodic_y && (k == -1);
  bool periodic_z = w->parent_object->periodic_z && (k == -1);

  // Lower left and upper right corners of the periodic box
  double llx = 0.0;
  double urx = 0.0;
  double lly = 0.0;
  double ury = 0.0;
  double llz = 0.0;
  double urz = 0.0;
  // in the presence of periodic boundary conditions we retrieve the box size
  if (periodic_x || periodic_y || periodic_z) {
    struct object* o = w->parent_object;
    assert(o->object_type == BOX_OBJ);
    struct polygon_object* p = (struct polygon_object*)(o->contents);
    struct subdivided_box* sb = p->sb;
    llx = sb->x[0];
    urx = sb->x[1];
    lly = sb->y[0];
    ury = sb->y[1];
    llz = sb->z[0];
    urz = sb->z[1];
  }

  double reflectFactor = -2.0 * (displacement->x * reflect_w->normal.x +
    displacement->y * reflect_w->normal.y + displacement->z *
    reflect_w->normal.z);

  int box_inc_x = 0;
  int box_inc_y = 0;
  int box_inc_z = 0;
  double x_pos = 0;
  double y_pos = 0;
  double z_pos = 0;

  // X direction: reflect or periodic BC
  if (periodic_x) {
    int x_inc = (vm->periodic_box->x % 2 == 0) ? 1 : -1;
    if (!distinguishable(vm->pos.x, llx, EPS_C)) {
      x_pos = urx - EPS_C;
      box_inc_x = -x_inc;
    } else if (!distinguishable(vm->pos.x, urx, EPS_C)) {
      x_pos = llx + EPS_C;
      box_inc_x = x_inc;
    }
    // Wrap molecule around to other side of box
    if (periodic_traditional && x_pos) {
      vm->pos.x = x_pos;
    }
  }
  // Set displacement for remainder of step length
  //
  // No PBCs or non-traditional PBCs
  if ((!periodic_x) || (periodic_x && !periodic_traditional)) {
    displacement->x = (displacement->x + reflectFactor * reflect_w->normal.x) *
      (1.0 - reflect_t);
  }
  // Traditional PBCs
  else {
    displacement->x *= (1.0 - reflect_t);
  }

  // Y direction: reflect or periodic BC
  if (periodic_y) {
    int y_inc = (vm->periodic_box->y % 2 == 0) ? 1 : -1;
    if (!distinguishable(vm->pos.y, lly, EPS_C)) {
      y_pos = ury - EPS_C;
      box_inc_y = -y_inc;
    } else if (!distinguishable(vm->pos.y, ury, EPS_C)) {
      y_pos = lly + EPS_C;
      box_inc_y = y_inc;
    }
    // Wrap molecule around to other side of box
    if (periodic_traditional && y_pos) {
      vm->pos.y = y_pos;
    }
  }

  // Set displacement for remainder of step length
  //
  // No PBCs or non-traditional PBCs
  if ((!periodic_y) || (periodic_y && !periodic_traditional)) {
    displacement->y = (displacement->y + reflectFactor * reflect_w->normal.y) *
      (1.0 - reflect_t);
  }
  // Traditional PBCs
  else {
    displacement->y *= (1.0 - reflect_t);
  }

  // Z direction: reflect or periodic BC
  if (periodic_z) {
    int z_inc = (vm->periodic_box->z % 2 == 0) ? 1 : -1;
    if (!distinguishable(vm->pos.z, llz, EPS_C)) {
      z_pos = urz - EPS_C;
      box_inc_z = -z_inc;
    } else if (!distinguishable(vm->pos.z, urz, EPS_C)) {
      z_pos = llz + EPS_C;
      box_inc_z = z_inc;
    }
    // Wrap molecule around to other side of box
    if (periodic_traditional && z_pos) {
      vm->pos.z =  z_pos;
    }
  }

  // Set displacement for remainder of step length
  //
  // No PBCs or non-traditional PBCs
  if ((!periodic_z) || (periodic_z && !periodic_traditional)) {
    displacement->z = (displacement->z + reflectFactor * reflect_w->normal.z) *
      (1.0 - reflect_t);
  }
  // Traditional PBCs
  else {
    displacement->z *= (1.0 - reflect_t);
  }

  // if we changed our position by periodic BC in either x, y or z we need
  // to check for migration to the proper subvolume.
  if ((periodic_traditional) && (periodic_x || periodic_y || periodic_z)) {
    (*reflectee) = NULL;
    struct subvolume *nsv = find_subvolume(world, &vm->pos, NULL);
    if (nsv == NULL) {
      struct species* spec = vm->properties;
      mcell_internal_error(
          "A %s molecule escaped the periodic box at [%.2f, %.2f, %.2f]",
          spec->sym->name, vm->pos.x * world->length_unit,
          vm->pos.y * world->length_unit, vm->pos.z * world->length_unit);
    } else {
      // decrement counts of regions we are leaving
      if (vm->properties->flags & (COUNT_CONTENTS | COUNT_ENCLOSED)) {
        count_region_from_scratch(world, (struct abstract_molecule *)vm, NULL,
                                  -1, &(orig_pos), NULL, reflect_t, NULL);
      }
      struct volume_molecule *new_m = migrate_volume_molecule(vm, nsv);
      // increment counts of regions we are entering
      if (new_m->properties->flags & (COUNT_CONTENTS | COUNT_ENCLOSED)) {
      count_region_from_scratch(world, (struct abstract_molecule *)new_m, NULL, 1,
                                &(new_m->pos), NULL, reflect_t, NULL);
      }
      *mol = new_m;
    }
    return 1;
  }

  if (!(periodic_traditional) && (box_inc_x || box_inc_y || box_inc_z)) {
    (*reflectee) = NULL;
    struct subvolume *nsv = find_subvolume(world, &vm->pos, NULL);
    if (nsv == NULL) {
      struct species* spec = vm->properties;
      mcell_internal_error(
          "A %s molecule escaped the periodic box at [%.2f, %.2f, %.2f]",
          spec->sym->name, vm->pos.x * world->length_unit,
          vm->pos.y * world->length_unit, vm->pos.z * world->length_unit);
    } else {
      // decrement counts of regions we are leaving
      if (vm->properties->flags & (COUNT_CONTENTS | COUNT_ENCLOSED)) {
        count_region_from_scratch(world, (struct abstract_molecule *)vm, NULL,
                                  -1, &(orig_pos), NULL, reflect_t,
                                  vm->periodic_box);
      }
      struct volume_molecule *new_m = migrate_volume_molecule(vm, nsv);
      vm->periodic_box->x += box_inc_x;
      vm->periodic_box->y += box_inc_y;
      vm->periodic_box->z += box_inc_z;
      // increment counts of regions we are entering
      if (new_m->properties->flags & (COUNT_CONTENTS | COUNT_ENCLOSED)) {
        count_region_from_scratch(world, (struct abstract_molecule *)new_m,
                                  NULL, 1, &(new_m->pos), NULL, reflect_t,
                                  new_m->periodic_box);
      }
      *mol = new_m;
    }
    return 1;
  }

  *mol = vm;
  return 0;
}


/******************************************************************************
 *
 * register_hits during 3D diffusion when colliding with a wall before
 * reflecting or encountering periodic boundary conditions.
 *
 * Return values:
 *
 * No return value.
 *
 * By default, we will reflect from the point of collision on the last
 * wall we hit; however, if there were one or more transparent walls we
 * hit at the same time or slightly before, we did not count them as
 * crossings, so we'd better be sure we don't cross them now.  Due to
 * round-off error, if we are counting, we need to make sure we don't
 * go back through these "tentative" surfaces again.  This involves
 * finding the first "tentative" surface, and traveling back a tiny
 * bit from that.
 * If we're doing counting, register hits for all "tentative" surfaces,
 * and update the point of reflection as explained above.
 *
 ******************************************************************************/
 void register_hits(struct volume* world, struct volume_molecule* m,
  struct collision** tentative, struct wall** reflect_w, double* reflect_t,
  struct vector3* displacement, struct collision* smash, double* t_steps) {

  struct collision* ttv = *tentative;
  struct species* spec = m->properties;
  struct vector3 reflect_pt = smash->loc;
  if ((m->flags & COUNT_ME) != 0 && (spec->flags & COUNT_SOME_MASK) != 0) {
    /* Find the first wall among the tentative collisions. */
    while (ttv != NULL && ttv->t <= smash->t &&
           !(ttv->what & COLLIDE_WALL)) {
      ttv = ttv->next;
    }
    assert(ttv != NULL);

    /* Grab out the relevant details. */
    (*reflect_w) = ((struct wall *)ttv->target);
    reflect_pt = ttv->loc;
    (*reflect_t) = ttv->t * (1 - EPS_C);

    /* Now, since we're reflecting before passing through these surfaces,
     * register them as hits, but not as crossings. */
    for (; ttv != NULL && ttv->t <= smash->t; ttv = ttv->next) {
      if (!(ttv->what & COLLIDE_WALL)) {
        continue;
      }
      if (!(spec->flags & ((struct wall *)ttv->target)->flags &
            COUNT_SOME_MASK)) {
        continue;
      }
      count_region_update(world, m->properties, m->periodic_box,
        ((struct wall *)ttv->target)->counting_regions,
        ((ttv->what & COLLIDE_MASK) == COLLIDE_FRONT) ? 1 : -1, 0, &(ttv->loc), ttv->t);
      if (ttv == smash)
        break;
    }
  }
  *tentative = ttv;

  /* Update molecule location to the point of reflection */
  m->pos = reflect_pt;
  m->t += *t_steps * (*reflect_t);

  /* Reduce our remaining available time. */
  *t_steps *= (1.0 - (*reflect_t));
}


/******************************************************************************
 *
 * the collide_and_react_with_subvol helper function is used in diffuse_3D to
 * collisions of diffusing molecule with a subvolume
 *
 * Return values:
 *
 * No return value. We just have to ensure that we update the counts properly
 * and then migrate the molecule to the proper subvolume.
 *
 ******************************************************************************/
void collide_and_react_with_subvol(struct volume* world, struct collision *smash,
  struct vector3* displacement, struct volume_molecule** mol,
  struct collision** tentative, double* t_steps) {

  struct collision* ttv = *tentative;
  struct volume_molecule* m = *mol;
  struct species* spec = m->properties;
  if ((m->flags & COUNT_ME) != 0 && (spec->flags & COUNT_SOME_MASK) != 0) {
    /* We're leaving the SV so we actually crossed everything we thought
     * we might have crossed */
    for (; ttv != NULL && ttv != smash; ttv = ttv->next) {
      if (!(ttv->what & COLLIDE_WALL)) {
        continue;
      }
      if (!(spec->flags & ((struct wall *)ttv->target)->flags & COUNT_SOME_MASK)) {
        continue;
      }
      count_region_update(world, spec, m->periodic_box,
          ((struct wall *)ttv->target)->counting_regions,
          ((ttv->what & COLLIDE_MASK) == COLLIDE_FRONT) ? 1 : -1, 1, &(ttv->loc), ttv->t);
    }
  }

  m->pos.x = smash->loc.x;
  m->pos.y = smash->loc.y;
  m->pos.z = smash->loc.z;

  displacement->x *= (1.0 - smash->t);
  displacement->y *= (1.0 - smash->t);
  displacement->z *= (1.0 - smash->t);

  m->t += (*t_steps) * smash->t;
  (*t_steps) *= (1.0 - smash->t);
  if (*t_steps < EPS_C) {
    *t_steps = EPS_C;
  }

  struct subvolume *nsv = traverse_subvol(
    m->subvol, smash->what - COLLIDE_SV_NX - COLLIDE_SUBVOL, world->ny_parts,
    world->nz_parts);
  if (nsv == NULL) {
    mcell_internal_error(
        "A %s molecule escaped the world at [%.2f, %.2f, %.2f]",
        spec->sym->name, m->pos.x * world->length_unit,
        m->pos.y * world->length_unit, m->pos.z * world->length_unit);
  } else {
    m = migrate_volume_molecule(m, nsv);
  }

  *mol = m;
  *tentative = ttv;
}


/******************************************************************************
 *
 * the compute_displacement helper function is used in diffuse_3D to compute the
 * displacement for the currently diffusing volume molecule
 *
 * Return values:
 *
 * this function does not return anything
 *
 ******************************************************************************/
void compute_displacement(struct volume* world, struct collision* shead,
  struct volume_molecule* m, struct vector3* displacement,
  struct vector3* displacement2, double* rate_factor, double* r_rate_factor,
  double* steps, double* t_steps, double max_time) {

  struct species* spec = m->properties;
  if (m->flags & ACT_CLAMPED) { /* Surface clamping and microscopic reversibility */
    if (m->index <= DISSOCIATION_MAX) { /* Volume microscopic reversibility */
      pick_release_displacement(displacement, displacement2, spec->space_step,
        world->r_step_release, world->d_step, world->radial_subdivisions,
        world->directions_mask, world->num_directions, world->rx_radius_3d,
        world->rng);
      *t_steps = 0;
    } else { /* Clamping or surface microscopic reversibility */
      pick_clamped_displacement(displacement, m, world->r_step_surface,
        world->rng, world->radial_subdivisions);
      *t_steps = spec->time_step;
      m->previous_wall = NULL;
      m->index = -1;
    }

    m->flags -= ACT_CLAMPED;
    *r_rate_factor = *rate_factor = 1.0;
    *steps = 1.0;
  } else {
    if (max_time > MULTISTEP_WORTHWHILE) {
      *steps = safe_diffusion_step(m, shead, world->radial_subdivisions,
        world->r_step, world->x_fineparts, world->y_fineparts, world->z_fineparts);
    } else {
      *steps = 1.0;
    }

    *t_steps = *steps * spec->time_step;
    if (*t_steps > max_time) {
      *t_steps = max_time;
      *steps = max_time / spec->time_step;
    }
    if (*steps < EPS_C) {
      *steps = EPS_C;
      *t_steps = EPS_C * spec->time_step;
    }

    if (*steps == 1.0) {
      pick_displacement(displacement, spec->space_step, world->rng);
      *r_rate_factor = *rate_factor = 1.0;
    } else {
      *rate_factor = sqrt(*steps);
      *r_rate_factor = 1.0 / *rate_factor;
      pick_displacement(displacement, *rate_factor * spec->space_step, world->rng);
    }
  }

  if (spec->flags & SET_MAX_STEP_LENGTH) {
    double disp_length = vect_length(displacement);
    if (disp_length > spec->max_step_length) {
      /* rescale displacement to the level of MAXIMUM_STEP_LENGTH */
      displacement->x *= (spec->max_step_length / disp_length);
      displacement->y *= (spec->max_step_length / disp_length);
      displacement->z *= (spec->max_step_length / disp_length);
    }
  }
  world->diffusion_number++;
  world->diffusion_cumtime += *steps;
}


/******************************************************************************
 *
 * the determine_mol_mol_reactions helper function is used in diffuse_3D to
 * compute all possible molecule molecule reactions between the diffusing
 * molecule m and all other volume molecules in the subvolume.
 *
 * Return values:
 *
 * this function does not return anything
 *
 ******************************************************************************/
void determine_mol_mol_reactions(struct volume* world, struct volume_molecule* m,
  struct collision** shead, struct collision** stail, int inertness) {

  struct subvolume* sv = m->subvol;
  struct rxn *matching_rxns[MAX_MATCHING_RXNS];
  int num_matching_rxns = 0;
  struct species* spec = m->properties;
  struct per_species_list *psl_next, *psl, **psl_head = &sv->species_head;
  for (psl = sv->species_head; psl != NULL; psl = psl_next) {
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

    /* no possible reactions. skip it. */
    if (!trigger_bimolecular_preliminary(world->reaction_hash, world->rx_hashsize,
      m->properties->hashval, psl->properties->hashval, m->properties, psl->properties)) {
      continue;
    }

    for (struct volume_molecule* mp = psl->head; mp != NULL; mp = mp->next_v) {
      if (mp == m) {
        continue;
      }

      if (inertness == inert_to_mol && m->index == mp->index) {
        continue;
      }

      // count only in the relevant periodic box
      if (!periodic_boxes_are_identical(m->periodic_box, mp->periodic_box)) {
        continue;
      }

      num_matching_rxns = trigger_bimolecular(world->reaction_hash,
        world->rx_hashsize, spec->hashval, psl->properties->hashval,
        (struct abstract_molecule *)m, (struct abstract_molecule *)mp, 0, 0,
        matching_rxns);

      if (num_matching_rxns > 0) {
        for (int i = 0; i < num_matching_rxns; i++) {
          struct collision* smash =
           (struct collision *)CHECKED_MEM_GET(sv->local_storage->coll,
            "collision data");
          smash->target = (void *)mp;
          smash->what = COLLIDE_VOL;
          smash->intermediate = matching_rxns[i];
          smash->next = *shead;
          *shead = smash;
          if (*stail == NULL)
            *stail = *shead;
        }
      }
    }
  }
}


/******************************************************************************
 *
 * the set_inertness_and_maxtime helper function is used in diffuse_3D to
 * set the maxtime used for diffusion as well as the inertness for clamped
 * molecules.
 *
 * Return values:
 *
 * this function does not return anything
 *
 ******************************************************************************/
void set_inertness_and_maxtime(
    struct volume* world, struct volume_molecule* m, double* max_time,
    int* inertness) {

  struct species* spec = m->properties;
  if (world->volume_reversibility || world->surface_reversibility) {
    if (world->volume_reversibility &&
        m->index <= DISSOCIATION_MAX) { /* Only set if volume_reversibility is */
      if ((m->flags & ACT_CLAMPED) != 0) {
        *inertness = inert_to_all;
      } else {
        m->index = -1;
      }
    } else if (!world->surface_reversibility) {
      if (m->flags & ACT_CLAMPED) { /* Pretend we were already moving */
        m->birthday -= 5 * spec->time_step; /* Pretend to be old */
      }
    }
  } else {
    if (m->flags & ACT_CLAMPED) { /* Pretend we were already moving */
      m->birthday -= 5 * spec->time_step; /* Pretend to be old */
    } else if ((m->flags & MATURE_MOLECULE) == 0) {
      /* Newly created particles that have long time steps gradually increase */
      /* their timestep to the full value */
      if (spec->time_step > 1.0) {
        double f = 1.0 + 0.2 * (m->t - m->birthday);
        if (f < 1)
          mcell_internal_error("A %s molecule is scheduled to move before it "
                               "was born [birthday=%.15g, t=%.15g]",
                               spec->sym->name, m->birthday * world->time_unit,
                               m->t * world->time_unit);
        if (*max_time > f) {
          *max_time = f;
        }
        if (f > m->subvol->local_storage->max_timestep) {
          m->flags |= MATURE_MOLECULE;
        }
      }
    }
  }
}

/*******************************************************************************
 *
 * count_tentative_collisions is a helper function counting hits of a diffusing
 * volume molecules with time-ordered collision targets along a ray up to time
 * t_confident.
 *
 ******************************************************************************/
void count_tentative_collisions(
  struct volume *world, struct collision **tc, struct collision *smash,
  struct species *spec, double t_confident, int destroy_flag,
  struct periodic_image *box) {

  int crossed_flag = 1;
  if (destroy_flag == 1) {
    // If we are destroyed, then we haven't crossed.
    crossed_flag = 0; 
  }
  struct collision *ttv = *tc;
  for (; ttv != NULL && ttv->t <= t_confident; ttv = ttv->next) {
    if (!(ttv->what & COLLIDE_WALL)) {
      continue;
    }
    if (!(spec->flags & ((struct wall *)ttv->target)->flags & COUNT_SOME_MASK)) {
      continue;
    }
    count_region_update(
      world, spec, box, ((struct wall *)ttv->target)->counting_regions,
      ((ttv->what & COLLIDE_MASK) == COLLIDE_FRONT) ? 1 : -1, crossed_flag,
      &(ttv->loc), ttv->t);
    if ((destroy_flag) && (ttv == smash)) {
      break;
    }
  }
  *tc = ttv;
}

/*************************************************************************
periodicbox_in_surfmol_list:
  Is the periodic box of one surface molecule (e.g. [0,0,0]) the same as any of
  the periodic boxes in a surface molecule list?
  In: periodic_box: the periodic box of a surface molecule
      sml: a list of surface molecules, each belonging to their own PB
  Out: true if the PB appears in the SM list, false otherwise.
*************************************************************************/
bool periodicbox_in_surfmol_list(
    struct periodic_image *periodic_box,
    struct surface_molecule_list *sml) {
  for (struct surface_molecule_list *sml_curr = sml;
       sml_curr != NULL;
       sml_curr = sml_curr->next) {
    struct surface_molecule *sm = sml_curr->sm;
    if (sm && periodic_boxes_are_identical(periodic_box, sm->periodic_box)) {
      return true;
    }
  }
  return false;
}
