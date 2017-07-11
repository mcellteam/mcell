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
* ****************************************************************************/

#include <stdlib.h>
#include <string.h>

#include "mcell_misc.h"
#include "init.h"
#include "mcell_dyngeom.h"
#include "mcell_objects.h"
#include "count_util.h"
#include "logging.h"
#include "dyngeom.h"

/* simple wrapper for executing the supplied function call. In case
 * of an error returns with MCELL_FAIL and prints out error_message */
#define CHECKED_CALL(function, error_message)                                  \
  {                                                                            \
    if (function) {                                                            \
      mcell_log(error_message);                                                \
      return MCELL_FAIL;                                                       \
    }                                                                          \
  }

int mcell_add_dynamic_geometry_file(char *dynamic_geometry_filepath,
                                    struct mdlparse_vars *parse_state) {
  struct volume *state = parse_state->vol;
  char *dynamic_geometry_filename =
      mcell_find_include_file(dynamic_geometry_filepath, state->curr_file);
  state->dynamic_geometry_filename = dynamic_geometry_filename;
#ifdef NOSWIG
  schedule_dynamic_geometry(parse_state);
#endif
  free(dynamic_geometry_filepath);
  return 0;
}

int mcell_change_geometry(struct volume *state, struct poly_object_list *pobj_list) {
  state->all_molecules = save_all_molecules(state, state->storage_head);

  // Turn off progress reports to avoid spamming mostly useless info to stdout
  state->notify->progress_report = NOTIFY_NONE;
  if (state->dynamic_geometry_flag != 1) {
    free(state->mdl_infile_name);
  }

  // Make list of already existing regions with fully qualified names.
  struct string_buffer *old_region_names =
      CHECKED_MALLOC_STRUCT(struct string_buffer, "string buffer");
  initialize_string_buffer(old_region_names, MAX_NUM_OBJECTS);
  get_reg_names_all_objects(state->root_instance, old_region_names);

  // Make list of already existing meshes with fully qualified names.
  struct string_buffer *old_inst_mesh_names = 
      CHECKED_MALLOC_STRUCT(struct string_buffer, "string buffer");
  initialize_string_buffer(old_inst_mesh_names, MAX_NUM_OBJECTS);
  get_mesh_instantiation_names(state->root_instance, old_inst_mesh_names);

  // We set this mainly to take care of some issues with counting, triggers,
  // memory cleanup.
  state->dynamic_geometry_flag = 1;

  CHECKED_CALL(reset_current_counts(
    state->mol_sym_table,
    state->count_hashmask,
    state->count_hash),
    "Error when reseting counters.");

  struct vector3 llf;
  struct vector3 urb;
  if (state->periodic_box_obj) {
    struct polygon_object* p = (struct polygon_object*)(state->periodic_box_obj->contents);
    struct subdivided_box* sb = p->sb;
    llf.x = sb->x[0];
    llf.y = sb->y[0];
    llf.z = sb->z[0];
       
    urb.x = sb->x[1];
    urb.y = sb->y[1];
    urb.z = sb->z[1];
  }

  CHECKED_CALL(destroy_everything(state), "Error when freeing memory.");
  // We need to reenable the ability to parse geometry
  state->disable_polygon_objects = 0;
  // Reparse the geometry and instantiations. Nothing else should be included
  // in these other MDLs.

  struct object *world_object = NULL;
  mcell_create_instance_object(state, "Scene", &world_object);

  while (pobj_list != NULL) {
    struct poly_object polygon = {
      pobj_list->obj_name,
      pobj_list->vertices,
      pobj_list->num_vert,
      pobj_list->connections,
      pobj_list->num_conn
    };
    struct object *new_mesh = NULL;
    mcell_create_poly_object(state, world_object, &polygon, &new_mesh);

    struct region *test_region = mcell_create_region(state, new_mesh, pobj_list->reg_name);
    mcell_set_region_elements(test_region, pobj_list->surf_reg_faces, 1);
    pobj_list = pobj_list->next;
  }

  CHECKED_CALL(init_bounding_box(state), "Error initializing bounding box.");
  free(state->subvol);
  if (state->periodic_box_obj) {
    mcell_create_periodic_box(state, "PERIODIC_BOX_INST", &llf, &urb);
  }
  CHECKED_CALL(init_partitions(state), "Error initializing partitions.");
  CHECKED_CALL(init_vertices_walls(state),
               "Error initializing vertices and walls.");
  CHECKED_CALL(init_regions(state), "Error initializing regions.");

  if (state->place_waypoints_flag) {
    CHECKED_CALL(place_waypoints(state), "Error while placing waypoints.");
  }

  if (state->with_checks_flag) {
    CHECKED_CALL(check_for_overlapped_walls(
        state->rng, state->n_subvols, state->subvol),
        "Error while checking for overlapped walls.");
  }
  CHECKED_CALL(init_species_mesh_transp(state),
               "Error while initializing species-mesh transparency list.");

  // Make NEW list of fully qualified region names.
  struct string_buffer *new_region_names =
      CHECKED_MALLOC_STRUCT(struct string_buffer, "string buffer");
  initialize_string_buffer(new_region_names, MAX_NUM_OBJECTS);
  get_reg_names_all_objects(state->root_instance, new_region_names);

  // Make NEW list of instantiated, fully qualified mesh names.
  struct string_buffer *new_inst_mesh_names =
      CHECKED_MALLOC_STRUCT(struct string_buffer, "string buffer");
  initialize_string_buffer(new_inst_mesh_names, MAX_NUM_OBJECTS);
  get_mesh_instantiation_names(state->root_instance, new_inst_mesh_names);

  struct string_buffer *meshes_to_ignore =
      CHECKED_MALLOC_STRUCT(struct string_buffer, "string buffer");
  initialize_string_buffer(meshes_to_ignore, MAX_NUM_OBJECTS);
  // Compare old list of mesh names with new list. 
  sym_diff_string_buffers(
    meshes_to_ignore, old_inst_mesh_names, new_inst_mesh_names,
    state->notify->add_remove_mesh_warning);

  struct string_buffer *regions_to_ignore =
      CHECKED_MALLOC_STRUCT(struct string_buffer, "string buffer");
  initialize_string_buffer(regions_to_ignore, MAX_NUM_OBJECTS);
  // Compare old list of region names with new list. 
  sym_diff_string_buffers(
    regions_to_ignore, old_region_names, new_region_names,
    state->notify->add_remove_mesh_warning);

  check_count_validity(
      state->output_request_head,
      regions_to_ignore,
      new_region_names,
      meshes_to_ignore,
      new_inst_mesh_names);
  place_all_molecules(state, meshes_to_ignore, regions_to_ignore);

  destroy_string_buffer(old_region_names);
  destroy_string_buffer(new_region_names);
  destroy_string_buffer(old_inst_mesh_names);
  destroy_string_buffer(new_inst_mesh_names);
  destroy_string_buffer(meshes_to_ignore);
  destroy_string_buffer(regions_to_ignore);
  free(old_region_names);
  free(new_region_names);
  free(old_inst_mesh_names);
  free(new_inst_mesh_names);
  free(meshes_to_ignore);
  free(regions_to_ignore);

  return 0;
}
