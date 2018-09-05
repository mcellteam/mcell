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
#include "mdlparse_aux.h"

int init_notifications(struct volume *world);
int init_variables(struct volume *world);
int init_data_structures(struct volume *world);
int init_species(struct volume *world);
int init_bounding_box(struct volume *world);
int init_partitions(struct volume *world);
int init_vertices_walls(struct volume *world);
int init_regions(struct volume *world);
int init_checkpoint_state(struct volume *world, long long *exec_iterations);
int init_viz_data(struct volume *world);
int init_reaction_data(struct volume *world);
int init_timers(struct volume *world);
int init_counter_name_hash(struct sym_table_head **counter_by_name,
                           struct output_block *output_block_head);

int parse_input(struct volume *world);

int load_checkpoint(struct volume *world);

int instance_obj(struct volume *world, struct object *objp, double (*im)[4]);

int instance_release_site(struct mem_helper *magic_mem,
                          struct schedule_helper *releaser, struct object *objp,
                          double (*im)[4]);

int instance_polygon_object(enum warn_level_t degenerate_polys,
                            struct object *objp);

void init_clamp_lists(struct ccn_clamp_data *clamp_list);

int instance_obj_regions(struct volume *world, struct object *objp);

int init_wall_regions(double length_unit, struct ccn_clamp_data *clamp_list,
                      struct species **species_list, int n_species,
                      struct object *objp);

int init_surf_mols(struct volume *world);
int instance_obj_surf_mols(struct volume *world, struct object *objp);
int init_wall_surf_mols(struct volume *world, struct object *objp);

int init_surf_mols_by_density(struct volume *world, struct wall *w,
                              struct sm_dat *sm_dat_head);

int init_surf_mols_by_number(struct volume *world, struct object *objp,
                             struct region_list *rlp);

void cube_corners(struct vector3 *p1, struct vector3 *p2,
                  struct vector3 *corner);

void cube_face(struct vector3 *corner, struct vector3 **face, int i);

void cube_faces(struct vector3 *corner, struct vector3 *(*face)[4]);

int init_releases(struct schedule_helper *releaser);

int init_dynamic_geometry(struct volume *state);

int schedule_dynamic_geometry(struct mdlparse_vars *parse_state);

int reschedule_release_events(struct volume *world);

void publish_special_reactions_report(struct species *sp,
                                      struct name_list *vol_species_name_list,
                                      struct name_list *surf_species_name_list,
                                      int n_species,
                                      struct species **species_list);

int accumulate_vertex_counts_per_storage(struct volume *world,
                                         struct object *objp,
                                         int *num_vertices_this_storage,
                                         double (*im)[4]);

int accumulate_vertex_counts_per_storage_polygon_object(
    struct volume *world, struct object *objp, int *num_vertices_this_storage,
    double (*im)[4]);

int which_storage_contains_vertex(struct volume *world, struct vector3 *v);

int fill_world_vertices_array(struct volume *world, struct object *objp,
                              int *num_vertices_this_storage, double (*im)[4]);

int fill_world_vertices_array_polygon_object(struct volume *world,
                                             struct object *objp,
                                             int *num_vertices_this_storage,
                                             double (*im)[4]);
void check_for_conflicting_surface_classes(struct wall *w, int n_species,
                                           struct species **species_list);
void check_for_conflicts_in_surface_class(struct volume *world,
                                          struct species *sp);
struct species *get_species_by_name(char *name, int n_species,
                                    struct species **species_list);
void create_name_lists_of_volume_and_surface_mols(
    struct volume *world, struct name_list **vol_species_name_list,
    struct name_list **surf_species_name_list);
void remove_molecules_name_list(struct name_list **nlist);
int check_for_overlapped_walls(
    struct rng_state *rng, int n_subvols, struct subvolume *subvol);
struct vector3 *create_region_bbox(struct region *r);
