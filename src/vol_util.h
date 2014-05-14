/***********************************************************************************
 *                                                                                 *
 * Copyright (C) 2006-2014 by *
 * The Salk Institute for Biological Studies and *
 * Pittsburgh Supercomputing Center, Carnegie Mellon University *
 *                                                                                 *
 * This program is free software; you can redistribute it and/or *
 * modify it under the terms of the GNU General Public License *
 * as published by the Free Software Foundation; either version 2 *
 * of the License, or (at your option) any later version. *
 *                                                                                 *
 * This program is distributed in the hope that it will be useful, *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the *
 * GNU General Public License for more details. *
 *                                                                                 *
 * You should have received a copy of the GNU General Public License *
 * along with this program; if not, write to the Free Software *
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301,
 *USA. *
 *                                                                                 *
 ***********************************************************************************/

#ifndef MCELL_VOL_UTIL
#define MCELL_VOL_UTIL

#include "mcell_structs.h"

int inside_subvolume(struct vector3 *point, struct subvolume *subvol,
                     double *x_fineparts, double *y_fineparts,
                     double *z_fineparts);

struct subvolume *find_coarse_subvol(struct volume *world, struct vector3 *loc);

struct subvolume *traverse_subvol(struct subvolume *here, struct vector3 *point,
                                  int which, int nx_parts, int ny_parts,
                                  int nz_parts);

struct subvolume *next_subvol(struct vector3 *here, struct vector3 *move,
                              struct subvolume *sv, double *x_fineparts,
                              double *y_fineparts, double *z_fineparts,
                              int nx_parts, int ny_parts, int nz_parts);

struct subvolume *find_subvolume(struct volume *world, struct vector3 *loc,
                                 struct subvolume *guess);

double collide_sv_time(struct vector3 *point, struct vector3 *move,
                       struct subvolume *sv, double *x_fineparts,
                       double *y_fineparts, double *z_fineparts);

int is_defunct_molecule(struct abstract_element *e);

struct grid_molecule *place_grid_molecule(struct volume *world,
                                          struct species *s,
                                          struct vector3 *loc, short orient,
                                          double search_diam, double t,
                                          struct subvolume **psv,
                                          struct grid_molecule **cmplx);

struct grid_molecule *insert_grid_molecule(struct volume *world,
                                           struct species *s,
                                           struct vector3 *loc, short orient,
                                           double search_diam, double t,
                                           struct grid_molecule **cmplx);

struct volume_molecule *insert_volume_molecule(struct volume *world,
                                               struct volume_molecule *m,
                                               struct volume_molecule *guess);

void exsert_volume_molecule(struct volume *world, struct volume_molecule *m);

int insert_volume_molecule_list(struct volume *world,
                                struct volume_molecule *m);

struct volume_molecule *migrate_volume_molecule(struct volume_molecule *m,
                                                struct subvolume *new_sv);

int eval_rel_region_3d(struct release_evaluator *expr, struct waypoint *wp,
                       struct region_list *in_regions,
                       struct region_list *out_regions);

int release_molecules(struct volume *world, struct release_event_queue *req);

void randomize_vol_mol_position(struct volume *world,
                                struct volume_molecule *mp,
                                struct vector3 *low_end, double size_x,
                                double size_y, double size_z);

int set_partitions(struct volume *world);

void path_bounding_box(struct vector3 *loc, struct vector3 *displacement,
                       struct vector3 *llf, struct vector3 *urb,
                       double rx_radius_3d);

void ht_add_molecule_to_list(struct pointer_hash *h, struct volume_molecule *m);
void ht_remove(struct pointer_hash *h, struct per_species_list *psl);

void collect_molecule(struct volume_molecule *m);

#endif
