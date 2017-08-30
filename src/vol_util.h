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

int inside_subvolume(struct vector3 *point, struct subvolume *subvol,
                     double *x_fineparts, double *y_fineparts,
                     double *z_fineparts);

struct subvolume *find_coarse_subvol(struct volume *world, struct vector3 *loc);

struct subvolume *traverse_subvol(struct subvolume *here,
                                  int which, int ny_parts, int nz_parts);

struct subvolume *next_subvol(struct vector3 *here, struct vector3 *move,
                              struct subvolume *sv, double *x_fineparts,
                              double *y_fineparts, double *z_fineparts,
                              int ny_parts, int nz_parts);

struct subvolume *find_subvolume(struct volume *world, struct vector3 *loc,
                                 struct subvolume *guess);

double collide_sv_time(struct vector3 *point, struct vector3 *move,
                       struct subvolume *sv, double *x_fineparts,
                       double *y_fineparts, double *z_fineparts);

int is_defunct_molecule(struct abstract_element *e);

struct wall* find_closest_wall(
    struct volume *state, struct vector3 *loc, double search_diam,
    struct vector2 *best_uv, int *grid_index, struct species *s, char *mesh_name,
    struct string_buffer *reg_names, struct string_buffer *regions_to_ignore);

struct surface_molecule *
place_surface_molecule(struct volume *state, struct species *s,
                       struct vector3 *loc, short orient, double search_diam,
                       double t, struct subvolume **psv, char *mesh_name,
                       struct string_buffer *reg_names,
                       struct string_buffer *regions_to_ignore,
                       struct periodic_image *periodic_box);

struct surface_molecule *
insert_surface_molecule(struct volume *state, struct species *s,
                        struct vector3 *loc, short orient, double search_diam,
                        double t, char *mesh_name,
                        struct string_buffer *reg_names,
                        struct string_buffer *regions_to_ignore,
                        struct periodic_image *periodic_box);

struct volume_molecule *insert_volume_molecule(struct volume *world,
                                               struct volume_molecule *vm,
                                               struct volume_molecule *guess);

struct volume_molecule *migrate_volume_molecule(struct volume_molecule *vm,
                                                struct subvolume *new_sv);

int eval_rel_region_3d(struct release_evaluator *expr, struct waypoint *wp,
                       struct region_list *in_regions,
                       struct region_list *out_regions);

int release_molecules(struct volume *world, struct release_event_queue *req);

int release_by_list(struct volume *state, struct release_event_queue *req,
                    struct volume_molecule *vm);

int release_ellipsoid_or_rectcuboid(struct volume *state,
                                    struct release_event_queue *req,
                                    struct volume_molecule *vm, int number);

int set_partitions(struct volume *world);

double increase_fine_partition_size(struct volume *state, double *fineparts,
                                    double *f_min, double *f_max,
                                    double smallest_spacing);

void set_fineparts(double min, double max, double *partitions,
                   double *fineparts, int n_parts, int in, int start);

void set_auto_partitions(struct volume *state, double steps_min,
                         double steps_max, struct vector3 *part_min,
                         struct vector3 *part_max, double f_max,
                         double smallest_spacing);

void set_user_partitions(struct volume *state, double dfx, double dfy,
                         double dfz);

void find_closest_fine_part(double *partitions, double *fineparts,
                            int n_fineparts, int n_parts);

double *add_extra_outer_partitions(double *partitions, double bb_llf_val,
                                   double bb_urb_val, double df_val,
                                   int *n_parts);

void path_bounding_box(struct vector3 *loc, struct vector3 *displacement,
                       struct vector3 *llf, struct vector3 *urb,
                       double rx_radius_3d);

void ht_add_molecule_to_list(struct pointer_hash *h, struct volume_molecule *vm);
void ht_remove(struct pointer_hash *h, struct per_species_list *psl);

void collect_molecule(struct volume_molecule *vm);

bool periodic_boxes_are_identical(const struct periodic_image *b1,
  const struct periodic_image *b2);

int convert_relative_to_abs_PBC_coords(
    struct object *periodic_box_obj,
    struct periodic_image *periodic_box,
    bool periodic_traditional,
    struct vector3 *pos,
    struct vector3 *pos_output);

struct surface_molecule_list* add_surfmol_with_unique_pb_to_list(
    struct surface_molecule_list *sm_list,
    struct surface_molecule *sm);

void remove_surfmol_from_list(
    struct surface_molecule_list **sm_head,
    struct surface_molecule *sm);
