/******************************************************************************
 *
 * Copyright (C) 2006-2014 by
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
 *****************************************************************************/

#ifndef DYNGEOM_H
#define DYNGEOM_H

struct object_transparency {
  struct object_transparency *next;
  char *obj_name;
  int in_to_out;
  int out_to_in;
};

struct molecule_info ** save_all_molecules(
    struct volume *state, struct storage_list *storage_head);

void save_common_molecule_properties(struct molecule_info *mol_info,
                                     struct abstract_molecule *am_ptr,
                                     struct string_buffer *reg_names,
                                     struct string_buffer *obj_names,
                                     char *mesh_name);

void save_volume_molecule(struct volume *state, struct molecule_info *mol_info,
                          struct abstract_molecule *am_ptr,
                          struct string_buffer **obj_names);

int save_surface_molecule(struct molecule_info *mol_info,
                          struct abstract_molecule *am_ptr,
                          struct string_buffer **reg_names,
                          char **mesh_name);

int place_all_molecules(struct volume *state);

void check_for_large_molecular_displacement(
    struct vector3 *old_pos,
    struct vector3 *new_pos,
    struct volume_molecule *vm,
    double *time_unit,
    enum warn_level_t large_molecular_displacement_warning);

char *compare_molecule_nesting(int *move_molecule,
                               int *out_to_in, 
                               struct string_buffer *obj_names_old,
                               struct string_buffer *obj_names_new,
                               struct object_transparency *obj_transp);

char *check_overlapping_meshes(
    int *move_molecule, int *out_to_in, int difference,
    struct string_buffer *compare_this, char *best_location,
    struct object_transparency *obj_transp);

struct volume_molecule *insert_volume_molecule_encl_mesh(
    struct volume *state, struct volume_molecule *vm,
    struct volume_molecule *vm_guess, struct string_buffer *obj_names_old);

struct string_buffer *find_enclosing_mesh_name(struct volume *state,
                                               struct volume_molecule *vm);

void place_mol_relative_to_mesh(
    struct volume *state, struct vector3 *loc, struct subvolume *sv,
    char *mesh_name, struct vector3 *new_pos, int out_to_in);

int destroy_everything(struct volume *state);
int destroy_objects(struct object *obj_ptr, int free_poly_flag);
int destroy_poly_object(struct object *obj_ptr, int free_poly_flag);

int reset_current_counts(struct volume *state);
int enable_counting_for_all_objects(struct object *obj_ptr);
int enable_counting_for_object(struct object *obj_ptr);

int init_mol_obj_transp(struct volume *world);
int find_obj_region_transp(struct object *obj_ptr,
                           struct object_transparency **obj_transp_head,
                           struct object_transparency **obj_transp_tail,
                           char *species_name);
int find_all_obj_region_transp(struct object *obj_ptr,
                               struct object_transparency **obj_transp_head,
                               struct object_transparency **obj_transp_tail,
                               char *species_name);

#endif
