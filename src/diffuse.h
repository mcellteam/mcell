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

#pragma once

#include "mcell_structs.h"

#define MULTISTEP_WORTHWHILE 2
#define MULTISTEP_PERCENTILE 0.99
#define MULTISTEP_FRACTION 0.9
#define MAX_UNI_TIMESKIP 100000

struct vector3* reflect_periodic_2D(
    struct volume *state,
    int index_edge_was_hit,
    struct vector2 *origin_uv,
    struct wall *curr_wall,
    struct vector2 *disp_uv,
    struct vector2 *boundary_uv,
    struct vector3 *origin_xyz);

void pick_displacement(struct vector3 *v, double scale, struct rng_state *rng);

void pick_2D_displacement(struct vector2 *v, double scale,
                          struct rng_state *rng);

void pick_release_displacement(struct vector3 *in_disk, struct vector3 *away,
                               double scale, double *r_step_release,
                               double *d_step, u_int radial_subdivisions,
                               int directions_mask, u_int num_directions,
                               double rx_radius_3d, struct rng_state *rng);

void pick_clamped_displacement(struct vector3 *v, struct volume_molecule *m,
                               double *r_step_surfce, struct rng_state *rng,
                               u_int radial_subdivision);

struct wall *ray_trace_2D(struct volume *world, struct surface_molecule *sm,
                          struct vector2 *disp, struct vector2 *loc,
                          int *kill_me, struct rxn **rxp,
                          struct hit_data **hd_info);

struct collision *ray_trace(struct volume *world, struct vector3 *init_pos,
                            struct collision *c, struct subvolume *sv,
                            struct vector3 *v, struct wall *reflectee);

struct sp_collision *ray_trace_trimol(struct volume *world,
                                      struct volume_molecule *m,
                                      struct sp_collision *c,
                                      struct subvolume *sv, struct vector3 *v,
                                      struct wall *reflectee,
                                      double walk_start_time);

struct volume_molecule *diffuse_3D(struct volume *world,
                                   struct volume_molecule *m, double max_time);

struct volume_molecule *diffuse_3D_big_list(struct volume *world,
                                            struct volume_molecule *m,
                                            double max_time);

struct surface_molecule *diffuse_2D(struct volume *world,
                                    struct surface_molecule *sm,
                                    double max_time, double *advance_time);

struct surface_molecule *react_2D(struct volume *world,
                                  struct surface_molecule *sm, double t,
                                  enum notify_level_t molecule_collision_report,
                                  int grid_grid_reaction_flag,
                                  long long *surf_surf_colls);

struct surface_molecule *
react_2D_all_neighbors(struct volume *world, struct surface_molecule *sm,
                       double t, enum notify_level_t molecule_collision_report,
                       int grid_grid_reaction_flag, long long *surf_surf_colls);

struct surface_molecule *react_2D_trimol_all_neighbors(
    struct volume *world, struct surface_molecule *sm, double t,
    enum notify_level_t molecule_collision_report,
    enum notify_level_t final_summary, int grid_grid_reaction_flag,
    long long *surf_surf_colls);

void clean_up_old_molecules(struct storage *local);

void reschedule_surface_molecules(
    struct volume *state, struct storage *local, struct abstract_molecule *am);

void run_timestep(struct volume *world, struct storage *local,
                  double release_time, double checkpt_time);

void run_concentration_clamp(struct volume *world, double t_now);

struct sp_collision *expand_collision_partner_list_for_neighbor(
    struct subvolume *sv, struct volume_molecule *m, struct vector3 *mv,
    struct subvolume *new_sv, struct vector3 *path_llf,
    struct vector3 *path_urb, struct sp_collision *shead1, double trim_x,
    double trim_y, double trim_z, double *x_fineparts, double *y_fineparts,
    double *z_fineparts, int rx_hashsize, struct rxn **reaction_hash);

double safe_diffusion_step(struct volume_molecule *m, struct collision *shead,
                           u_int radial_subdivisions, double *r_step,
                           double *x_fineparts, double *y_fineparts,
                           double *z_fineparts);

double exact_disk(struct volume *world, struct vector3 *loc, struct vector3 *mv,
                  double R, struct subvolume *sv,
                  struct volume_molecule *moving,
                  struct volume_molecule *target, int use_expanded_list,
                  double *x_fineparts, double *y_fineparts,
                  double *z_fineparts);

bool periodicbox_in_surfmol_list(
    struct periodic_image *periodic_box,
    struct surface_molecule_list *sml);
