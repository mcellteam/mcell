/******************************************************************************
 *
 * Copyright (C) 2006-2017 by
 * The Salk Institute for Biological Studies and
 * Pittsburgh Supercomputing Center, Carnegie Mellon University
 *
 * Use of this source code is governed by an MIT-style
 * license that can be found in the LICENSE file or at
 * https://opensource.org/licenses/MIT.
 *
******************************************************************************/

#pragma once

#include "mcell_structs.h"
#include "dyngeom_parse_extras.h"
#include "vol_util.h"

int region_listed(struct region_list *rl, struct region *r);

void count_region_update(
    struct volume *world,
    struct volume_molecule *vm,
    struct species *sp,
    u_long id,
    struct periodic_image *img,
    struct region_list *rl,
    int direction,
    int crossed,
    struct vector3 *loc,
    double t);

void count_region_border_update(struct volume *world, struct species *sp,
                                struct hit_data *hd_info, u_long id);

void count_region_from_scratch(struct volume *world,
                               struct abstract_molecule *am,
                               struct rxn_pathname *rxpn, int n,
                               struct vector3 *loc, struct wall *my_wall,
                               double t, struct periodic_image *periodic_box);

void count_moved_surface_mol(struct volume *world, struct surface_molecule *sm,
  struct surface_grid *sg, struct vector2 *loc, int count_hashmask,
  struct counter **count_hash, long long *ray_polygon_colls,
  struct periodic_image *previous_box);

void fire_count_event(struct volume *world, struct counter *event, int n,
                      struct vector3 *where, byte what, u_long id);

int place_waypoints(struct volume *world);

int prepare_counters(struct volume *world);

int check_counter_geometry(int count_hashmask, struct counter **count_hash,
                           byte *place_waypoints_flag);

int expand_object_output(struct output_request *request,
                         struct geom_object *obj,
                         struct sym_table_head *reg_sym_table);
int object_has_geometry(struct geom_object *obj);

/* hit data for region borders */
void update_hit_data(struct hit_data **hd_head, struct wall *current,
                     struct wall *target, struct surface_molecule *sm,
                     struct vector2 boundary_pos, int direction, int crossed);
int is_object_instantiated(struct sym_entry *entry,
                           struct geom_object *root_instance);
