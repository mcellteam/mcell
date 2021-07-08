/******************************************************************************
 *
 * Copyright (C) 2006-2014 by
 * The Salk Institute for Biological Studies and
 * Pittsburgh Supercomputing Center, Carnegie Mellon University
 *
 * Use of this source code is governed by an MIT-style
 * license that can be found in the LICENSE file or at
 * https://opensource.org/licenses/MIT.
 *
 *****************************************************************************/

#ifndef DYNGEOM_H
#define DYNGEOM_H

#include "mdlparse_aux.h"

#define MAX_NUM_REGIONS 100
#define MAX_NUM_OBJECTS 100

struct mesh_transparency {
  struct mesh_transparency *next;
  char *name;
  int in_to_out;
  int out_to_in;
  int transp_top_front;
  int transp_top_back;
};

struct name_hits {
  struct name_hits *next;
  char *name; /* molecule name */
  int hits; /* molecule orientation */
};

struct n_parts {
  int nx_parts;
  int ny_parts;
  int nz_parts;
};

struct molecule_info ** save_all_molecules(
    struct volume *state, struct storage_list *storage_head);

void save_common_molecule_properties(struct molecule_info *mol_info,
                                     struct abstract_molecule *am_ptr,
                                     struct string_buffer *reg_names,
                                     struct string_buffer *mesh_names,
                                     char *mesh_name);

void save_volume_molecule(struct volume *state, struct molecule_info *mol_info,
                          struct abstract_molecule *am_ptr,
                          struct string_buffer **mesh_names);

int save_surface_molecule(struct molecule_info *mol_info,
                          struct abstract_molecule *am_ptr,
                          struct string_buffer **reg_names,
                          char **mesh_name);

void cleanup_names_molecs(
    int num_all_molecules, struct molecule_info **all_molecules);

int place_all_molecules(
    struct volume *state,
    struct string_buffer *names_to_ignore,
    struct string_buffer *regions_to_ignore);

void check_for_large_molecular_displacement(
    struct vector3 *old_pos,
    struct vector3 *new_pos,
    struct volume_molecule *vm,
    double *time_unit,
    enum warn_level_t large_molecular_displacement_warning);

const char *compare_molecule_nesting(int *move_molecule,
                               int *out_to_in, 
                               struct string_buffer *mesh_names_old,
                               struct string_buffer *mesh_names_new,
                               struct mesh_transparency *mesh_transp);

const char *check_overlapping_meshes(
    int *move_molecule, int *out_to_in, int difference,
    struct string_buffer *compare_this, const char *best_mesh,
    struct mesh_transparency *mesh_transp);

const char *check_nonoverlapping_meshes(int *move_molecule,
                                  int *out_to_in,
                                  struct string_buffer *mesh_names_old,
                                  struct string_buffer *mesh_names_new,
                                  const char *best_mesh,
                                  struct mesh_transparency *mesh_transp);

const char *check_outin_or_inout(
    int start, int increment, int end, int *move_molecule,
    int *out_to_in, const char *best_mesh, struct string_buffer *mesh_names,
    struct mesh_transparency *mesh_transp);

struct volume_molecule *insert_volume_molecule_encl_mesh(
    struct volume *state,
    struct volume_molecule *vm,
    struct volume_molecule *vm_guess,
    struct string_buffer *mesh_names_old,
    struct string_buffer *names_to_ignore);

int hit_wall(
    struct wall *w, struct name_hits **name_head,
    struct name_hits **name_tail, struct vector3 *rand_vector);

void hit_subvol(
    struct n_parts *np, struct string_buffer *mesh_names,
    struct collision *smash, struct collision *shead,
    struct name_hits *name_head, struct subvolume *sv,
    struct volume_molecule *virt_mol);

struct string_buffer *find_enclosing_meshes(
    struct volume *state,
    struct volume_molecule *vm,
    struct string_buffer *names_to_ignore);

void place_mol_relative_to_mesh(
    struct volume *state, struct vector3 *loc, struct subvolume *sv,
    const char *mesh_name, struct vector3 *new_pos, int out_to_in);

void destroy_mesh_transp_data(
    struct sym_table_head *mol_sym_table,
    struct pointer_hash *species_mesh_transp);
int destroy_everything(struct volume *state);
void destroy_walls(struct volume *state);
void destroy_partitions(struct volume *state);
int destroy_objects(struct geom_object *obj_ptr, int free_poly_flag);
int destroy_poly_object(struct geom_object *obj_ptr, int free_poly_flag);

int reset_current_counts(struct sym_table_head *mol_sym_table,
                         int count_hashmask,
                         struct counter **count_hash);

void check_count_validity(struct output_request *output_request_head,
                          struct string_buffer *regions_to_ignore,
                          struct string_buffer *new_region_names,
                          struct string_buffer *meshes_to_ignore,
                          struct string_buffer *new_mesh_names);

void reset_count_type(const char *name,
                      struct output_request *request,
                      struct string_buffer *names_to_ignore,
                      struct string_buffer *new_names);

int init_species_mesh_transp(struct volume *world);

int find_sm_region_transp(struct geom_object *obj_ptr,
                          struct mesh_transparency **mesh_transp_head,
                          struct mesh_transparency **mesh_transp_tail,
                          const char *species_name);

void check_surf_class_properties(
    const char *species_name, struct mesh_transparency *mesh_transp,
  struct name_orient *surf_class_props);

int find_vm_obj_region_transp(struct geom_object *obj_ptr,
                              struct mesh_transparency **mesh_transp_head,
                              struct mesh_transparency **mesh_transp_tail,
                              const char *species_name);

int find_all_obj_region_transp(struct geom_object *obj_ptr,
                               struct mesh_transparency **mesh_transp_head,
                               struct mesh_transparency **mesh_transp_tail,
                               const char *species_name, int sm_flag);

int add_dynamic_geometry_events(
    struct mdlparse_vars *parse_state,
    const char *dynamic_geometry_filepath,
    double timestep,
    struct mem_helper *dynamic_geometry_events_mem,
    struct dg_time_filename **dg_time_fname_head);

const char *get_mesh_instantiation_names(struct geom_object *obj_ptr,
                                   struct string_buffer *mesh_names);

void diff_string_buffers(
    struct string_buffer *diff_names,
    struct string_buffer *names_a,
    struct string_buffer *names_b);

void sym_diff_string_buffers(
    struct string_buffer *diff_names,
    struct string_buffer *names_a,
    struct string_buffer *names_b,
    enum warn_level_t add_remove_mesh_warning);

int get_reg_names_all_objects(
    struct geom_object *obj_ptr, struct string_buffer *regions_to_ignore);

int get_reg_names_this_object(
    struct geom_object *obj_ptr, struct string_buffer *regions_to_ignore);

void update_geometry(struct volume *state, struct dg_time_filename *dyn_geom);

#endif
