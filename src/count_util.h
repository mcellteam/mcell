/******************************************************************************
 *
 * Copyright (C) 2006-2015 by
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

#ifndef MCELL_COUNT_UTIL
#define MCELL_COUNT_UTIL

#include "mcell_structs.h"

int region_listed(struct region_list *rl, struct region *r);

void count_region_update(struct volume *world, struct species *sp,
                         struct region_list *rl, int direction, int crossed,
                         struct vector3 *loc, double t);

void count_region_border_update(struct volume *world, struct species *sp,
                                struct hit_data *hd_info);

void count_region_from_scratch(struct volume *world,
                               struct abstract_molecule *am,
                               struct rxn_pathname *rxpn, int n,
                               struct vector3 *loc, struct wall *my_wall,
                               double t);

void count_moved_surface_mol(struct volume *world, struct surface_molecule *sm,
                             struct surface_grid *sg, struct vector2 *loc,
                             int count_hashmask, struct counter **count_hash);

void fire_count_event(struct volume *world, struct counter *event, int n,
                      struct vector3 *where, byte what);

int place_waypoints(struct volume *world);

int prepare_counters(struct volume *world);

int check_counter_geometry(int count_hashmask, struct counter **count_hash,
                           byte *place_waypoints_flag);

int expand_object_output(struct output_request *request, struct object *obj);
int object_has_geometry(struct object *obj);

/* hit data for region borders */
void update_hit_data(struct hit_data **hd_head, struct wall *current,
                     struct wall *target, struct surface_molecule *sm,
                     struct vector2 boundary_pos, int direction, int crossed);

/************************************************************
 * Complex counting
 ************************************************************/
int count_complex(struct volume *world, struct volume_molecule *cmplex,
                  struct volume_molecule *replaced_subunit,
                  int replaced_subunit_idx);
int count_complex_surface(struct surface_molecule *cmplex,
                          struct surface_molecule *replaced_subunit,
                          int replaced_subunit_idx);
int count_complex_surface_new(struct surface_molecule *cmplex);

#endif
