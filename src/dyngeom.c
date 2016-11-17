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

#include "config.h"

#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "chkpt.h"
#include "vol_util.h"
#include "grid_util.h"
#include "wall_util.h"
#include "init.h"
#include "logging.h"
#include "count_util.h"
#include "diffuse.h"
#include "dyngeom.h"
#include "mcell_misc.h"
#include "dyngeom_parse_extras.h"
#include "mdlparse_aux.h"
#include "react.h"

#define NO_MESH "\0"

/***************************************************************************
 save_all_molecules: Save all the molecules currently in the scheduler.

 In:  state: MCell state
      storage_head: we will pull all the molecules out of the scheduler from
        this
 Out: An array of all the molecules to be saved
***************************************************************************/
struct molecule_info **save_all_molecules(struct volume *state,
                                          struct storage_list *storage_head) {

  // Find total number of molecules in the scheduler.
  unsigned long long num_all_molecules = count_items_in_scheduler(storage_head);
  int ctr = 0;
  struct molecule_info **all_molecules = CHECKED_MALLOC_ARRAY(
      struct molecule_info *, num_all_molecules, "all molecules");

  // Iterate over all the molecules in every scheduler of every storage.
  for (struct storage_list *sl_ptr = storage_head; sl_ptr != NULL;
       sl_ptr = sl_ptr->next) {
    for (struct schedule_helper *sh_ptr = sl_ptr->store->timer; sh_ptr != NULL;
         sh_ptr = sh_ptr->next_scale) {
      for (int i = -1; i < sh_ptr->buf_len; i++) {
        for (struct abstract_element *ae_ptr =
                 (i < 0) ? sh_ptr->current : sh_ptr->circ_buf_head[i];
             ae_ptr != NULL; ae_ptr = ae_ptr->next) {
          struct abstract_molecule *am_ptr = (struct abstract_molecule *)ae_ptr;
          if (am_ptr->properties == NULL)
            continue;

          struct molecule_info *mol_info =
              CHECKED_MALLOC_STRUCT(struct molecule_info, "molecule info");
          all_molecules[ctr] = mol_info;
          mol_info->molecule = CHECKED_MALLOC_STRUCT(struct abstract_molecule,
                                                     "abstract molecule");

          // Mesh names needed for VOLUME molecules
          struct string_buffer *nested_mesh_names = NULL;

          // Region names and mesh name needed for SURFACE molecules
          struct string_buffer *reg_names =
              CHECKED_MALLOC_STRUCT(struct string_buffer, "string buffer");
          if (initialize_string_buffer(reg_names, MAX_NUM_REGIONS)) {
            return NULL;
          }
          char *mesh_name = NULL;

          if ((am_ptr->properties->flags & NOT_FREE) == 0) {
            save_volume_molecule(state, mol_info, am_ptr, &nested_mesh_names);
          } else if ((am_ptr->properties->flags & ON_GRID) != 0) {
            if (save_surface_molecule(mol_info, am_ptr, &reg_names, &mesh_name))
              return NULL;
          } else {
            continue;
          }

          save_common_molecule_properties(
              mol_info, am_ptr, reg_names, nested_mesh_names, mesh_name);
          ctr += 1;
        }
      }
    }
  }

  state->num_all_molecules = ctr;

  return all_molecules;
}

/***************************************************************************
 save_common_molecule_properties:

 In:  mol_info: holds all the information for recreating and placing a molecule
      am_ptr: abstract molecule pointer
      reg_names: the region names the molecule is on (surface molecules)
      nested_mesh_names: the meshes the molecule is nested inside of (volume molecs)
      mesh_name: mesh name that molecule is in (surface molecs)
 Out: Nothing. The common properties of surface and volume molecules are saved
      in mol_info.
***************************************************************************/
void save_common_molecule_properties(struct molecule_info *mol_info,
                                     struct abstract_molecule *am_ptr,
                                     struct string_buffer *reg_names,
                                     struct string_buffer *nested_mesh_names,
                                     char *mesh_name) {
  mol_info->molecule->t = am_ptr->t;
  mol_info->molecule->t2 = am_ptr->t2;
  mol_info->molecule->flags = am_ptr->flags;
  mol_info->molecule->properties = am_ptr->properties;
  mol_info->molecule->birthday = am_ptr->birthday;
  mol_info->molecule->id = am_ptr->id;
  mol_info->molecule->periodic_box = am_ptr->periodic_box;
  mol_info->molecule->mesh_name = CHECKED_STRDUP(mesh_name, "mesh name");
  // Only free temporary object names we just allocated above.
  // Don't want to accidentally free symbol names of objects.
  if (mesh_name && (strcmp(mesh_name, NO_MESH) != 0) &&
      ((am_ptr->properties->flags & NOT_FREE) == 0)) {
    free(mesh_name);
  }
  mol_info->reg_names = reg_names;
  mol_info->mesh_names = nested_mesh_names;
}

/***************************************************************************
 save_volume_molecule:

 In:  state: MCell state
      mol_info: holds all the information for recreating and placing a molecule
      am_ptr: abstract molecule pointer
      nested_mesh_names: mesh names that molecule is inside of
 Out: Nothing. Molecule info and mesh name are updated
***************************************************************************/
void save_volume_molecule(struct volume *state,
                          struct molecule_info *mol_info,
                          struct abstract_molecule *am_ptr,
                          struct string_buffer **nested_mesh_names) {
  struct volume_molecule *vm_ptr = (struct volume_molecule *)am_ptr;

  *nested_mesh_names = find_enclosing_meshes(state, vm_ptr, NULL);
  mol_info->pos.x = vm_ptr->pos.x;
  mol_info->pos.y = vm_ptr->pos.y;
  mol_info->pos.z = vm_ptr->pos.z;
  mol_info->orient = 0;
}

/***************************************************************************
 save_surface_molecule:

 In:  mol_info: holds all the information for recreating and placing a molecule
      am_ptr: abstract molecule pointer
      reg_names: surface region names that the molecule is on get stored here
      mesh_name: mesh name that molecule is on gets stored here
 Out: Zero on success. One otherwise. Save relevant surface molecule data in
      mol_info. Set the mesh_name. Populate reg_names with the region names the
      sm is on.
***************************************************************************/
int save_surface_molecule(struct molecule_info *mol_info,
                          struct abstract_molecule *am_ptr,
                          struct string_buffer **reg_names,
                          char **mesh_name) {
  struct vector3 where;
  struct surface_molecule *sm_ptr = (struct surface_molecule *)am_ptr;
  uv2xyz(&sm_ptr->s_pos, sm_ptr->grid->surface, &where);
  mol_info->pos.x = where.x;
  mol_info->pos.y = where.y;
  mol_info->pos.z = where.z;
  mol_info->orient = sm_ptr->orient;
  *mesh_name = sm_ptr->grid->surface->parent_object->sym->name;
  struct name_list *reg_name_list_head, *reg_name_list;
  reg_name_list_head = find_regions_names_by_wall(sm_ptr->grid->surface, NULL);
  // Add the names from reg_name_list_head to reg_names
  for (reg_name_list = reg_name_list_head; reg_name_list != NULL;
       reg_name_list = reg_name_list->next) {
    char *str = CHECKED_STRDUP(reg_name_list->name, "region name");
    if (add_string_to_buffer(*reg_names, str)) {
      free(str);
      destroy_string_buffer(*reg_names);
      return 1;
    }
  }
  if (reg_name_list_head != NULL) {
    remove_molecules_name_list(&reg_name_list_head);
  }
  remove_surfmol_from_list(&sm_ptr->grid->sm_list[sm_ptr->grid_index], sm_ptr);
  return 0;
}

/***************************************************************************
 cleanup_names_molecs: Cleanup molecule data and string buffers for mesh and
                       region names

 In:  num_all_molecules: the number of all the molecules we just placed
      all_molecules: array of info about all the molecules we just placed
 Out: Nothing
***************************************************************************/
void cleanup_names_molecs(
    int num_all_molecules,
    struct molecule_info **all_molecules) {
  for (int i = 0; i < num_all_molecules; i++) {
    char *mesh_name = all_molecules[i]->molecule->mesh_name;
    if (mesh_name && (strcmp(mesh_name, NO_MESH) != 0)) {
      free(mesh_name);
    }
    free(all_molecules[i]->molecule);
    // XXX: Could consolidate "reg_names" and "mesh_names", but we might want
    // both for surface molecules eventually.
    destroy_string_buffer(all_molecules[i]->reg_names);
    free(all_molecules[i]->reg_names);
    // This is currently unused for surface molecules.
    if (all_molecules[i]->mesh_names != NULL) {
      destroy_string_buffer(all_molecules[i]->mesh_names);
      free(all_molecules[i]->mesh_names);
    }
    free(all_molecules[i]);
  }
  free(all_molecules);
}

/***************************************************************************
 place_all_molecules: Place all molecules currently in the scheduler into the
                      world.

 In:  state: MCell state
      meshes_to_ignore: don't place molecules on these meshes
      regions_to_ignore: don't place molecules on these regions
 Out: Zero on success. One otherwise.
***************************************************************************/
int place_all_molecules(
    struct volume *state,
    struct string_buffer *meshes_to_ignore,
    struct string_buffer *regions_to_ignore) {

  struct volume_molecule vm;
  memset(&vm, 0, sizeof(struct volume_molecule));
  struct volume_molecule *vm_ptr = &vm;
  struct volume_molecule *vm_guess = NULL;

  int num_all_molecules = state->num_all_molecules;

  for (int n_mol = 0; n_mol < num_all_molecules; n_mol++) {

    struct molecule_info *mol_info = state->all_molecules[n_mol];
    struct abstract_molecule *am_ptr = mol_info->molecule;
    // Insert volume molecule into world.
    if ((am_ptr->properties->flags & NOT_FREE) == 0) {
      vm_ptr->t = am_ptr->t;
      vm_ptr->t2 = am_ptr->t2;
      vm_ptr->flags = am_ptr->flags;
      vm_ptr->properties = am_ptr->properties;
      vm_ptr->birthday = am_ptr->birthday;
      vm_ptr->id = am_ptr->id;
      vm_ptr->pos.x = mol_info->pos.x;
      vm_ptr->pos.y = mol_info->pos.y;
      vm_ptr->pos.z = mol_info->pos.z;
      vm_ptr->periodic_box = am_ptr->periodic_box;

      vm_guess = insert_volume_molecule_encl_mesh(
          state, vm_ptr, vm_guess, mol_info->mesh_names, meshes_to_ignore);

      if (vm_guess == NULL) {
        mcell_error("Cannot insert copy of molecule of species '%s' into "
                    "world.\nThis may be caused by a shortage of memory.",
                    vm_ptr->properties->sym->name);
      }
    }
    // Insert surface molecule into world.
    else if ((am_ptr->properties->flags & ON_GRID) != 0) {
      char *mesh_name = am_ptr->mesh_name;
      struct surface_molecule *sm = insert_surface_molecule(
          state, am_ptr->properties, &mol_info->pos, mol_info->orient,
          state->vacancy_search_dist2, am_ptr->t, mesh_name,
          mol_info->reg_names, regions_to_ignore, am_ptr->periodic_box);
      free(am_ptr->periodic_box);
      if (sm == NULL) {
        mcell_warn("Unable to find surface upon which to place molecule %s.",
                   am_ptr->properties->sym->name);
      }
    }
  }

  cleanup_names_molecs(state->num_all_molecules, state->all_molecules);

  return 0;
}

/***************************************************************************
 compare_molecule_nesting:

  Compare where the molecule was (prior to the geometry change) with where
  it is now relative to the meshes. This also takes into account meshes nested
  inside other meshes.

  First, a little discussion on notation for the sake of being clear and
  concise. When we say something like A->B->C->null, this means that mesh A is
  inside mesh B which is inside mesh C. Lastly, C is an outermost mesh. Now,
  onto the algorithm itself...
  
  First, assume there is overlap, which we will check later.
  For example:
    mesh_names_old: A->B->C->D->null
    mesh_names_new:       C->D->null
  A was the closest enclosing mesh and then C became the closest enclosing
  mesh. That means the molecule moved from A to C.
  
  The order could be reversed like this:
    mesh_names_old:       C->D->null
    mesh_names_new: A->B->C->D->null
  This means the molecule moved from C to A.
  
  Check the last overlapping entry for both. If they actually overlap, these
  strings will be the same. Otherwise, this is nonoverlapping (aside for null)
  like this:
    mesh_names_old: C->D->null
    mesh_names_new: E->F->null
  This means the molecule moved from C to E.

  Next, we see if movement is possible from the starting position to ending
  position.

 In: move_molecule: if set, we need to move the molecule
     out_to_in: if set, the molecule moved from outside to inside
     mesh_names_old: a list of names that the molecule was nested in
     mesh_names_new: a list of names that the molecule is nested in
     mesh_transp: the object transparency rules for this species
 Out: The name of the mesh that we are either immediately inside or outside of.
      Also move_molecule and out_to_in are set.
***************************************************************************/
char *compare_molecule_nesting(int *move_molecule,
                               int *out_to_in, 
                               struct string_buffer *mesh_names_old,
                               struct string_buffer *mesh_names_new,
                               struct mesh_transparency *mesh_transp) {

  int old_n_strings = mesh_names_old->n_strings;
  int new_n_strings = mesh_names_new->n_strings;
  int difference;
  char *old_mesh_name;
  char *new_mesh_name;
  char *best_mesh = mesh_names_old->strings[0];
  struct string_buffer *compare_this;

  // mesh_names_old example:       C->D->null
  // mesh_names_new example: A->B->C->D->null
  if (old_n_strings < new_n_strings) {
    difference = new_n_strings - old_n_strings;
    old_mesh_name = mesh_names_old->strings[0];
    new_mesh_name = mesh_names_new->strings[difference];
    compare_this = mesh_names_new;
    *out_to_in = 1;
  }
  // mesh_names_old example: A->B->C->D->null
  // mesh_names_new example:       C->D->null
  else if (new_n_strings < old_n_strings) {
    difference = old_n_strings - new_n_strings;
    old_mesh_name = mesh_names_old->strings[difference];
    new_mesh_name = mesh_names_new->strings[0];
    compare_this = mesh_names_old;
    *out_to_in = 0;
  }
  // Same amount of nesting
  else {
    difference = 0;
    old_mesh_name = mesh_names_old->strings[0];
    new_mesh_name = mesh_names_new->strings[0];
    // Doesn't really matter if we use old or new one
    compare_this = mesh_names_old; 
  }

  if (old_mesh_name == NULL) {
    old_mesh_name = NO_MESH;
  }
  if (new_mesh_name == NULL) {
    new_mesh_name = NO_MESH;
  }
  if (best_mesh == NULL) {
    best_mesh = NO_MESH;
  }
  int names_match = (strcmp(old_mesh_name, new_mesh_name) == 0);
  if (names_match && (difference != 0)) {
    best_mesh = check_overlapping_meshes(
        move_molecule, out_to_in, difference, compare_this, best_mesh,
        mesh_transp);
  }
  else if (!names_match) {
    best_mesh = check_nonoverlapping_meshes(
        move_molecule, out_to_in, mesh_names_old, mesh_names_new, best_mesh,
        mesh_transp);
  }

  return best_mesh;
}

/***************************************************************************
 check_overlapping_meshes:

 In: move_molecule: if set, we need to move the molecule
     out_to_in: if set, the molecule moved from outside to inside
     difference: number of different meshes between old and new
     compare_this: the mesh names to check
     best_mesh: the current best mesh name
     mesh_transp: the object transparency rules for this species
 Out: The name of the mesh that we are either immediately inside or outside of.
      Also move_molecule and out_to_in are set.
***************************************************************************/
char *check_overlapping_meshes(
    int *move_molecule,
    int *out_to_in,
    int difference,
    struct string_buffer *compare_this,
    char *best_mesh,
    struct mesh_transparency *mesh_transp) {
  int start;
  int increment;
  int end;
  // The molecule moved from *outside* a mesh to *inside* a mesh
  if (*out_to_in) {
    // We offset the starting value by one, because, even though the molecule
    // started at the mesh corresponding to the "difference" index value, we
    // only care if it was stopped by the next mesh inward.
    start = difference-1;
    increment = -1;
    end = 0;
  }
  // The molecule moved from *inside* a mesh to *outside* a mesh
  else {
    start = 0;
    increment = 1;
    end = difference;
  }
  best_mesh = check_outin_or_inout(
      start, increment, end, move_molecule, out_to_in, best_mesh,
      compare_this, mesh_transp);
  if (strcmp(best_mesh, NO_MESH) == 0) {
    return NULL;
  }
  else {
    return best_mesh;
  }
}

/***************************************************************************
 check_nonoverlapping_meshes:

 By nonoverlapping, we mean that a molecule moved from one mesh (or set of
 nested meshes) to another mesh (or nested set) which it neither contains nor
 is inside of. For example, a molecule in A->B->null that moved to C->D->null
 would be considered nonoverlapping. A is not in C and C is not in A.

 In: move_molecule: if set, we need to move the molecule
     out_to_in: if set, the molecule moved from outside to inside
     mesh_names_old: the nested mesh names prior to this dyngeom event
     mesh_names_new: the nested mesh names during to this dyngeom event
     best_mesh: the mesh the molecule should actully be inside
     mesh_transp: the object transparency rules for this species
 Out: The name of the mesh that we are either immediately inside or outside of.
      Also move_molecule and out_to_in are set.
***************************************************************************/
char *check_nonoverlapping_meshes(int *move_molecule,
                                  int *out_to_in,
                                  struct string_buffer *mesh_names_old,
                                  struct string_buffer *mesh_names_new,
                                  char *best_mesh,
                                  struct mesh_transparency *mesh_transp) {

  *out_to_in = 0;
  int start = 0;
  int increment = 1;
  int end = mesh_names_old->n_strings;

  // Moving in to out
  best_mesh = check_outin_or_inout(
      start, increment, end, move_molecule, out_to_in, best_mesh,
      mesh_names_old, mesh_transp); 

  // Moving out to in
  if (!(*move_molecule)) {
    *out_to_in = 1;
    start = mesh_names_new->n_strings;
    end = 0;
    increment = -1;
    best_mesh = check_outin_or_inout(
        start, increment, end, move_molecule, out_to_in, best_mesh,
        mesh_names_new, mesh_transp); 
  }

  return best_mesh;
}

/***************************************************************************
 check_outin_or_inout:
 
 See if a molecule can move from through the meshes (mesh_names) in the
 direction specified (out_to_in). If it has to stop, return the name of the
 mesh that blocks it.

 In: start: start checking at this index
     increment: move forward or backward through the string buffer
     end: stop checking at this index
     move_molecule: if set, we need to move the molecule
     out_to_in: if set, the molecule moved from outside to inside
     best_mesh: the current best mesh name
     mesh_names: mesh names that molecule is inside of
     mesh_transp: the object transparency rules for this species
 Out: The name of the mesh that we are either immediately inside or outside of.
      Also move_molecule and out_to_in are set.
***************************************************************************/
char *check_outin_or_inout(
    int start,
    int increment,
    int end,
    int *move_molecule,
    int *out_to_in,
    char *best_mesh,
    struct string_buffer *mesh_names,
    struct mesh_transparency *mesh_transp) {
  int done = 0;
  int mesh_idx = start;
  while (!done) {
    char *mesh_name = mesh_names->strings[mesh_idx];
    if (mesh_name == NULL) {
      mesh_name = NO_MESH;
    }
    struct mesh_transparency *mt = mesh_transp;
    for (; mt != NULL; mt = mt->next) {
      if (strcmp(mesh_name, mt->name) == 0) {
        if (((*out_to_in) && !mt->out_to_in) ||
            (!(*out_to_in) && !mt->in_to_out)) {
          best_mesh = mesh_name;
          *move_molecule = 1;
          done = 1;
        }
        break;
      }
    }
    if (mesh_idx == end) {
      done = 1; 
    }
    mesh_idx = mesh_idx + increment;
  }
  return best_mesh;
}

/*************************************************************************
insert_volume_molecule_encl_mesh:
  In: state: MCell state
      vm: pointer to volume_molecule that we're going to place in local storage
      vm_guess: pointer to a volume_molecule that may be nearby
      nested_mesh_names_old: the meshes this molecule was inside of previously
      meshes_to_ignore: the meshes we should ignore when placing this molecule
  Out: pointer to the new volume_molecule (copies data from volume molecule
       passed in), or NULL if out of memory.  Molecule is placed in scheduler
       also.
*************************************************************************/
struct volume_molecule *insert_volume_molecule_encl_mesh(
    struct volume *state,
    struct volume_molecule *vm,
    struct volume_molecule *vm_guess,
    struct string_buffer *nested_mesh_names_old,
    struct string_buffer *meshes_to_ignore) {
  struct subvolume *sv;

  // We should only have to do this the first time this function gets called
  if (vm_guess == NULL)
    sv = find_subvolume(state, &(vm->pos), NULL);
  // This should speed things up if the last molecule was close to this one
  else if (inside_subvolume(&(vm->pos), vm_guess->subvol, state->x_fineparts,
                            state->y_fineparts, state->z_fineparts)) {
    sv = vm_guess->subvol;
  } else
    sv = find_subvolume(state, &(vm->pos), vm_guess->subvol);

  struct volume_molecule *new_vm = CHECKED_MEM_GET(
    sv->local_storage->mol, "volume molecule");
  memcpy(new_vm, vm, sizeof(struct volume_molecule));
  new_vm->mesh_name = NULL;
  new_vm->prev_v = NULL;
  new_vm->next_v = NULL;
  new_vm->next = NULL;
  new_vm->subvol = sv;
  new_vm->periodic_box = vm->periodic_box;

  struct string_buffer *nested_mesh_names_new = find_enclosing_meshes(
      state, new_vm, meshes_to_ignore);

  // Make a new string buffer without all the meshes we don't care about (i.e.
  // the ones we *removed* in this dyn_geom_event). We are already ingoring the
  // ones just *added* in this dyn_geom_event (in find_enclosing_meshes). Maybe
  // we should do that here to be consistent and keep the logic decoupled.
  struct string_buffer *nested_mesh_names_old_filtered =
      CHECKED_MALLOC_STRUCT(struct string_buffer, "string buffer");
  initialize_string_buffer(nested_mesh_names_old_filtered, MAX_NUM_OBJECTS);
  diff_string_buffers(
    nested_mesh_names_old_filtered, nested_mesh_names_old, meshes_to_ignore);

  char *species_name = new_vm->properties->sym->name;
  unsigned int keyhash = (unsigned int)(intptr_t)(species_name);
  void *key = (void *)(species_name);
  struct mesh_transparency *mesh_transp = (
      struct mesh_transparency *)pointer_hash_lookup(state->species_mesh_transp,
                                                     key, keyhash);

  int move_molecule = 0;
  int out_to_in = 0;
  char *mesh_name = compare_molecule_nesting(
    &move_molecule,
    &out_to_in,
    nested_mesh_names_old_filtered,
    nested_mesh_names_new,
    mesh_transp);

  struct vector3 new_pos;
  if (move_molecule) {
    /* move molecule to another location so that it is directly inside or
     * outside of "mesh_name" */
    place_mol_relative_to_mesh(
        state, &(vm->pos), sv, mesh_name, &new_pos, out_to_in);
    check_for_large_molecular_displacement(
        &(vm->pos), &new_pos, vm, &(state->time_unit),
        state->notify->large_molecular_displacement);
    new_vm->pos = new_pos;
    struct subvolume *new_sv = find_subvolume(state, &(new_vm->pos), NULL);
    new_vm->subvol = new_sv;
    state->dyngeom_molec_displacements++;
  }

  destroy_string_buffer(nested_mesh_names_old_filtered);
  free(nested_mesh_names_old_filtered);
  destroy_string_buffer(nested_mesh_names_new);
  free(nested_mesh_names_new);

  new_vm->birthplace = new_vm->subvol->local_storage->mol;
  ht_add_molecule_to_list(&(new_vm->subvol->mol_by_species), new_vm);
  new_vm->subvol->mol_count++;
  new_vm->properties->population++;

  if ((new_vm->properties->flags & COUNT_SOME_MASK) != 0) {
    new_vm->flags |= COUNT_ME;
  }
  // XXX: need to set periodic box properly
  if (new_vm->properties->flags & (COUNT_CONTENTS | COUNT_ENCLOSED)) {
    count_region_from_scratch(state, (struct abstract_molecule *)new_vm, NULL,
                              1, &(new_vm->pos), NULL, new_vm->t,
                              new_vm->periodic_box);
  }

  if (schedule_add(new_vm->subvol->local_storage->timer, new_vm))
    mcell_allocfailed("Failed to add volume molecule to scheduler.");

  return new_vm;
}

/*************************************************************************
check_for_large_molecular_displacement:
  In:  old_pos: The current position of the molecule
       new_pos: The position we are trying to move the molecule to
       vm: volume molecule
       timestep: global timestep in seconds
       large_molecular_displacement_warning: the warning value 
          (ignore, warn, error) set for molecular displacement
  Out: 0 on success, 1 otherwise.
************************************************************************/
void check_for_large_molecular_displacement(
    struct vector3 *old_pos,
    struct vector3 *new_pos,
    struct volume_molecule *vm,
    double *timestep,
    enum warn_level_t large_molecular_displacement_warning) {

  double displacement = distance_vec3(old_pos, new_pos) / 100.0;
  double l_perp_bar = sqrt(4 * 1.0e8 * vm->properties->D * *timestep / MY_PI);
  double l_r_bar = 2 * l_perp_bar;
  if (displacement >= l_r_bar) {
    switch (large_molecular_displacement_warning) {
    case WARN_COPE:
      break;

    case WARN_WARN:
      mcell_warn("Displacement of '%s' is greater than l_r_bar.\n"
                 "\tdisplacement = %.9g microns\n"
                 "\tl_r_bar = %.9g microns\n",
                 vm->properties->sym->name, displacement, l_r_bar);
      break;

    case WARN_ERROR:
      mcell_error("Displacement of '%s' is greater than l_r_bar.\n"
                  "\tdisplacement = %.9g microns\n"
                  "\tl_r_bar = %.9g microns\n",
                  vm->properties->sym->name, displacement, l_r_bar);

    }
  }
}

/*************************************************************************
hit_wall:
  In:  w: wall
       name_hits_head: the head of a list that tracks the meshes we've hit and
         how many times they've been hit
       name_tail_head: the tail of the same list
       displace_vector: a large displacement vector that spans the whole
         simulation space
  Out: Head and tail are updated. In other words, we track meshes we've hit
************************************************************************/
int hit_wall(
    struct wall *w,
    struct name_hits **name_hits_head,
    struct name_hits **name_tail_head,
    struct vector3 *displace_vector) {

  // Discard open-type meshes, like planes, etc.
  if (w->parent_object->is_closed <= 0)
    return 1;

  /* Discard the cases when the random vector just grazes the mesh at the
   * encounter point */
  double d_prod = dot_prod(displace_vector, &(w->normal));
  if (!distinguishable(d_prod, 0, EPS_C))
    return 1;

  // First time hitting *any* object
  if (*name_hits_head == NULL) {
    // name and count of hits
    struct name_hits *nh = CHECKED_MALLOC_STRUCT(
        struct name_hits, "struct name_hits");
    nh->name = CHECKED_STRDUP(w->parent_object->sym->name,
                              "w->parent_object->sym->name");
    nh->hits = 1;
    nh->next = NULL;
    *name_hits_head = nh;
    *name_tail_head = *name_hits_head;
  // We've hit at least one object already
  } else {
    int found = 0; // flag
    for (struct name_hits *nhl = *name_hits_head; nhl != NULL; nhl = nhl->next) {
      // Keep track of how many times we hit *this* object
      if (strcmp(nhl->name, w->parent_object->sym->name) == 0) {
        nhl->hits++;
        found = 1;
        break;
      }
    }
    // First time hitting *this* object
    if (!found) {
      // Add to the end of list
      struct name_hits *nh =
          CHECKED_MALLOC_STRUCT(struct name_hits, "struct name_hits");
      nh->name = CHECKED_STRDUP(w->parent_object->sym->name,
                                "w->parent_object->sym->name");
      nh->hits = 1;
      nh->next = (*name_tail_head)->next;
      (*name_tail_head)->next = nh;
      *name_tail_head = (*name_tail_head)->next;
    }
  }
  return 0;
}

/*************************************************************************
hit_subvol:
  In:  state: MCell state
       mesh_names: meshes that the molecule is inside of
       smash: the thing that the current molecule has collided with
       shead: the head of a list of what the current molecule has collided with
       name_hits_head: the head of a list that tracks the meshes we've hit and
         how many times they've been hit
       sv: subvolume
       virt_mol:
  Out: Compile list of meshes we are inside of (hit odd number of times) or
       update next subvolume
************************************************************************/
void hit_subvol(
    struct n_parts *np,
    struct string_buffer *mesh_names,
    struct collision *smash,
    struct collision *shead,
    struct name_hits *name_hits_head,
    struct subvolume *sv,
    struct volume_molecule *virt_mol) {

  virt_mol->pos.x = smash->loc.x;
  virt_mol->pos.y = smash->loc.y;
  virt_mol->pos.z = smash->loc.z;
  virt_mol->subvol = NULL;

  struct subvolume *new_sv = traverse_subvol(
      sv, smash->what - COLLIDE_SV_NX - COLLIDE_SUBVOL,
      np->ny_parts, np->nz_parts);
  // Hit the edge of the world
  if (new_sv == NULL) {
    if (shead != NULL)
      mem_put_list(sv->local_storage->coll, shead);

    // Compile the final list of meshes (names and counts) that we are inside
    for (struct name_hits *nhl = name_hits_head; nhl != NULL; nhl = nhl->next) {
      if (nhl->hits % 2 != 0) {
        char *mesh_name = CHECKED_STRDUP(nhl->name, "mesh name");
        if (add_string_to_buffer(mesh_names, mesh_name)) {
          free(mesh_name);
          destroy_string_buffer(mesh_names);
        }
      }
    }
    
    // Clean up
    while (name_hits_head != NULL) {
      struct name_hits *nnext = name_hits_head->next;
      free(name_hits_head->name);
      free(name_hits_head);
      name_hits_head = nnext;
    }
    return;
  }

  if (shead != NULL)
    mem_put_list(sv->local_storage->coll, shead);
  virt_mol->subvol = new_sv;
}

/*************************************************************************
find_enclosing_meshes:
  In:  state: MCell state
       vm: volume molecule
       meshes_to_ignore: ignore these meshes when checking what this molecule
         is inside of
  Out: String buffer of enclosing meshes if they exist. NULL otherwise.
************************************************************************/
struct string_buffer *find_enclosing_meshes(
    struct volume *state,
    struct volume_molecule *vm,
    struct string_buffer *meshes_to_ignore) {

  // We create a virtual molecule, so that we don't displace the real one (vm).
  struct volume_molecule virt_mol;
  memcpy(&virt_mol, vm, sizeof(struct volume_molecule));
  virt_mol.prev_v = NULL;
  virt_mol.next_v = NULL;
  virt_mol.next = NULL;

  // This is where we will store the names of the meshes we are nested in.
  struct string_buffer *mesh_names =
      CHECKED_MALLOC_STRUCT(struct string_buffer, "string buffer");
  if (initialize_string_buffer(mesh_names, MAX_NUM_OBJECTS)) {
    return NULL;
  }

  // Displacement vector along arbitray cardinal axis.
  struct vector3 displace_vector = {0.0, 0.0, 1.0};
  // Find the diagonal of the world's bounding box
  struct vector3 world_diag;
  vectorize(&state->bb_urb, &state->bb_llf, &world_diag);
  // Length of the world's bounding box diagonal
  double world_diag_length = vect_length(&world_diag);
  // Set world diagonal to nonzero value so we don't fail in traverse_subvol
  if (!distinguishable(world_diag_length, 0, EPS_C)) {
    world_diag_length = 1.0;
  }

  // Scale displacement vector by the length of the world
  displace_vector.x *= world_diag_length;
  displace_vector.y *= world_diag_length;
  displace_vector.z *= world_diag_length;

  struct collision *smash; /* Thing we've hit that's under consideration */
  struct collision *shead = NULL; // Head of the linked list of collisions
  struct subvolume *sv = virt_mol.subvol;
  struct name_hits *nh_head = NULL, *nh_tail = NULL;
  do {
    // Get collision list for walls and a subvolume. We don't care about
    // colliding with other molecules like we do with reactions
    shead = ray_trace(state, &(virt_mol.pos), NULL, sv, &displace_vector, NULL);
    if (shead == NULL)
      mcell_internal_error("ray_trace() returned NULL.");

    if (shead->next != NULL) {
      shead =
          (struct collision *)ae_list_sort((struct abstract_element *)shead);
    }

    for (smash = shead; smash != NULL; smash = smash->next) {
      // We hit a wall
      if ((smash->what & COLLIDE_WALL) != 0) {
        // Only check this when we are placing molecules, not when we are
        // saving them.
        struct wall *w = (struct wall *)smash->target;
        if ((meshes_to_ignore) && (is_string_present_in_string_array(
              w->parent_object->sym->name,
              meshes_to_ignore->strings,
              meshes_to_ignore->n_strings))) {
          continue; 
        }
        if (hit_wall(w, &nh_head, &nh_tail, &displace_vector)) {
          continue;
        }

      // We hit a subvolume
      } else if ((smash->what & COLLIDE_SUBVOL) != 0) {

        // Numbers of coarse partitions
        struct n_parts np = {
          state->nx_parts, state->ny_parts, state->nz_parts};
        hit_subvol(&np, mesh_names, smash, shead, nh_head, sv, &virt_mol);
        // We hit the edge of the world
        if (virt_mol.subvol == NULL) {
          return mesh_names; 
        }
        sv = virt_mol.subvol;
        break;
      }
    }

  } while (smash != NULL);

  return NULL;
}

/**********************************************************************
place_mol_relative_to_mesh:
  In: state: MCell state
      loc: 3D location of molecule
      sv: start subvolume
      mesh_name: name of closest enclosing mesh (NULL means that molecule
        is outside of all meshes)
      new_pos: new position of molecule (return value)
      out_to_in: if set, the molecule moved from outside to inside
  Note: new position of molecule that is just behind the closest
       wall that belongs to object called "mesh_name" (if "mesh_name != NULL")
       or just outside the closest wall that belongs to the farthest
       enclosing mesh object (if "mesh_name == NULL").
  Note: we call this function when geometry changes after checkpoint so that:
        1) before checkpoint molecule was inside the mesh, and after
           checkpoint it is outside the mesh,
        2) before checkpoint molecule was outside all meshes, and after
           checkpoint it is inside the mesh,
        3) before checkpoint molecule was inside one mesh, and after
           checkpoint it is inside another mesh.
        Here by mesh we mean the closest enclosing mesh.
**********************************************************************/
void place_mol_relative_to_mesh(struct volume *state,
                                struct vector3 *loc,
                                struct subvolume *sv,
                                char *mesh_name,
                                struct vector3 *new_pos,
                                int out_to_in) {
  struct vector2 s_loc;
  struct vector2 best_s_loc;
  struct wall *best_w = NULL;
  double d2;
  double best_d2 = GIGANTIC + 1;

  for (struct wall_list *wl = sv->wall_head; wl != NULL; wl = wl->next) {
    if (strcmp(wl->this_wall->parent_object->sym->name, mesh_name) != 0) {
      continue;
    }

    d2 = closest_interior_point(loc, wl->this_wall, &s_loc, GIGANTIC);
    if (d2 < best_d2) {
      best_d2 = d2;
      best_w = wl->this_wall;
      best_s_loc = s_loc;
    }
  }

  // Look into neighbor subvolumes
  const int sv_index = sv - state->subvol;
  int sv_remain = sv_index;

  // Turn linear sv_index into part_x, part_y, part_z triple.
  const int part_x =
      sv_remain / ((state->ny_parts - 1) * (state->nz_parts - 1));
  sv_remain -= part_x * ((state->ny_parts - 1) * (state->nz_parts - 1));
  const int part_y = sv_remain / (state->nz_parts - 1);
  sv_remain -= part_y * (state->nz_parts - 1);
  const int part_z = sv_remain;

  // Find min x partition.
  int x_min;
  for (x_min = part_x; x_min > 0; x_min--) {
    d2 = loc->x - state->x_partitions[x_min];
    d2 *= d2;
    if (d2 >= best_d2)
      break;
  }

  // Find max x partition.
  int x_max;
  for (x_max = part_x; x_max < state->nx_parts - 1; x_max++) {
    d2 = loc->x - state->x_partitions[x_max + 1];
    d2 *= d2;
    if (d2 >= best_d2)
      break;
  }

  // Find min y partition.
  int y_min;
  for (y_min = part_y; y_min > 0; y_min--) {
    d2 = loc->y - state->y_partitions[y_min];
    d2 *= d2;
    if (d2 >= best_d2)
      break;
  }

  // Find max y partition.
  int y_max;
  for (y_max = part_y; y_max < state->ny_parts - 1; y_max++) {
    d2 = loc->y - state->y_partitions[y_max + 1];
    d2 *= d2;
    if (d2 >= best_d2)
      break;
  }

  // Find min z partition.
  int z_min;
  for (z_min = part_z; z_min > 0; z_min--) {
    d2 = loc->z - state->z_partitions[z_min];
    d2 *= d2;
    if (d2 >= best_d2)
      break;
  }

  // Find max z partition.
  int z_max;
  for (z_max = part_z; z_max < state->nz_parts - 1; z_max++) {
    d2 = loc->z - state->z_partitions[z_max + 1];
    d2 *= d2;
    if (d2 >= best_d2)
      break;
  }

  if (x_min < part_x || x_max > part_x || y_min < part_y || y_max > part_y ||
      z_min < part_z || z_max > part_z) {
    for (int px = x_min; px < x_max; px++) {
      for (int py = y_min; py < y_max; py++) {
        for (int pz = z_min; pz < z_max; pz++) {
          const int this_sv =
              pz + (state->nz_parts - 1) * (py + (state->ny_parts - 1) * px);
          if (this_sv == sv_index)
            continue;

          for (struct wall_list *wl = state->subvol[this_sv].wall_head;
               wl != NULL; wl = wl->next) {
            if (strcmp(wl->this_wall->parent_object->sym->name, mesh_name) != 0) {
              continue;
            }

            d2 = closest_interior_point(loc, wl->this_wall, &s_loc, GIGANTIC);
            if (d2 < best_d2) {
              best_d2 = d2;
              best_w = wl->this_wall;
              best_s_loc = s_loc;
            }
          }
        }
      }
    }
  }

  if (best_w == NULL) {
    mcell_internal_error("Error in function 'place_mol_relative_to_mesh()'.");
  }

  // We will return the point just behind or in front of the closest enclosing
  // mesh

  /* the parametric equation of ray is L(t) = A + t(B - A),
     where A - start vector, B - some point on the ray, and parameter t >= 0 */

  struct vector3 v;
  if (state->dynamic_geometry_molecule_placement == 1) {
    double s1 = sqrt(rng_dbl(state->rng));
    double s2 = rng_dbl(state->rng) * s1;

    v.x = best_w->vert[0]->x + s1 * (best_w->vert[1]->x - best_w->vert[0]->x) + s2 * (best_w->vert[2]->x - best_w->vert[1]->x);
    v.y = best_w->vert[0]->y + s1 * (best_w->vert[1]->y - best_w->vert[0]->y) + s2 * (best_w->vert[2]->y - best_w->vert[1]->y);
    v.z = best_w->vert[0]->z + s1 * (best_w->vert[1]->z - best_w->vert[0]->z) + s2 * (best_w->vert[2]->z - best_w->vert[1]->z);
  }
  else if (state->dynamic_geometry_molecule_placement == 0) {
    uv2xyz(&best_s_loc, best_w, &v);
  }

  double bump = (out_to_in > 0) ? EPS_C : -EPS_C;
  struct vector3 displacement = { .x = 2 * bump * best_w->normal.x,
                                  .y = 2 * bump * best_w->normal.y,
                                  .z = 2 * bump * best_w->normal.z,
                                };
  struct subvolume *new_sv = find_subvolume(state, &v, NULL);
  tiny_diffuse_3D(state, new_sv, &displacement, &v, best_w);

  // Make sure we didn't end up on a neighbor's wall, which is kind of easy to
  // do with, for example, a shrinking box/cuboid.
  for (int edge = 0; edge < 3; edge++) {
    struct wall *neighbor = NULL;
    if (best_w != best_w->edges[edge]->forward) {
      neighbor = best_w->edges[edge]->forward;
    }
    else if (best_w != best_w->edges[edge]->backward) {
      neighbor = best_w->edges[edge]->backward;
    }
    else {
      continue;
    }
    d2 = closest_interior_point(&v, neighbor, &best_s_loc, GIGANTIC);
    if (!distinguishable(d2, 0, EPS_C)) {
      new_sv = find_subvolume(state, &v, NULL);
      bump = (out_to_in > 0) ? EPS_C : -EPS_C;
      displacement.x = 2 * bump * neighbor->normal.x;
      displacement.y = 2 * bump * neighbor->normal.y;
      displacement.z = 2 * bump * neighbor->normal.z;
      new_sv = find_subvolume(state, &v, NULL);
      tiny_diffuse_3D(state, new_sv, &displacement, &v, neighbor);
      break;
    }
  }

  new_pos->x = v.x;
  new_pos->y = v.y;
  new_pos->z = v.z;
}


/***************************************************************************
destroy_mesh_transp_data:
  In:  mol_sym_table:
       species_mesh_transp:
  Out: Destroy mesh-species transparency data structure 
***************************************************************************/
void destroy_mesh_transp_data(
    struct sym_table_head *mol_sym_table,
    struct pointer_hash *species_mesh_transp) {
  for (int n_mol_bin = 0; n_mol_bin < mol_sym_table->n_bins; n_mol_bin++) {
    for (struct sym_entry *sym_ptr = mol_sym_table->entries[n_mol_bin];
         sym_ptr != NULL; sym_ptr = sym_ptr->next) {
      char *species_name = sym_ptr->name;
      if (strcmp(species_name, "ALL_MOLECULES") == 0)
        continue;
      else if (strcmp(species_name, "ALL_VOLUME_MOLECULES") == 0)
        continue;
      else if (strcmp(species_name, "ALL_SURFACE_MOLECULES") == 0)
        continue;
      unsigned int keyhash = (unsigned int)(intptr_t)(species_name);
      void *key = (void *)(species_name);
      struct mesh_transparency *mesh_transp = (
          struct mesh_transparency *)pointer_hash_lookup(
              species_mesh_transp, key, keyhash);
      delete_void_list((struct void_list *)mesh_transp);
    }
  }
  pointer_hash_destroy(species_mesh_transp);
  free(species_mesh_transp);
}

/***************************************************************************
destroy_everything:
  In:  state: MCell state
  Out: Zero on success. One otherwise. This wipes out almost everything in the
       simulation except for things like the symbol tables, reactions, etc.
       Currently, this is necessary for dynamic geometries. In principle, it
       would make more sense to only trash the meshes changing, but there are
       so many tightly coupled dependencies that it's difficult to do it at
       this time.
***************************************************************************/
int destroy_everything(struct volume *state) {
  destroy_objects(state->root_instance, 1);
  destroy_objects(state->root_object, 0);
  state->root_instance->first_child = NULL; 
  state->root_instance->last_child = NULL; 
  state->root_instance->n_walls = 0;
  state->root_instance->n_walls_actual = 0;
  state->root_instance->n_verts = 0;
  state->root_object->first_child = NULL; 
  state->root_object->last_child = NULL; 
  state->root_object->n_walls = 0;
  state->root_object->n_walls_actual = 0;
  state->root_object->n_verts = 0;

  if (state->clamp_list) {
    free(state->clamp_list->side_idx);
    free(state->clamp_list->cum_area);
  }

  destroy_walls(state);

  // Destroy memory helpers
  delete_mem(state->coll_mem);
  delete_mem(state->exdv_mem);

  struct storage_list *mem;
  for (mem = state->storage_head; mem != NULL; mem = mem->next) {
    delete_mem(mem->store->list);
    delete_mem(mem->store->mol);
    delete_mem(mem->store->smol);
    delete_mem(mem->store->face);
    delete_mem(mem->store->join);
    delete_mem(mem->store->grids);
    delete_mem(mem->store->regl);
    delete_mem(mem->store->pslv);
  }

  // Destroy subvolumes
  for (int i = 0; i < state->n_subvols; i++) {
    struct subvolume *sv = &state->subvol[i];
    pointer_hash_destroy(&sv->mol_by_species);
    sv->local_storage->wall_head = NULL;
    sv->local_storage->wall_count = 0;
    sv->local_storage->vert_count = 0;
    sv->wall_head = NULL;
  }

  for (mem = state->storage_head; mem != NULL; mem = mem->next) {
    delete_scheduler(mem->store->timer);
    free(mem->store);
  }
  state->storage_head->store = NULL;
  state->storage_head = NULL;

  delete_mem(state->storage_allocator);
  delete_mem(state->sp_coll_mem);
  delete_mem(state->tri_coll_mem);

  destroy_partitions(state);

  free(state->waypoints);

  // Destroy mesh-species transparency data structure
  destroy_mesh_transp_data(state->mol_sym_table, state->species_mesh_transp);

  for (struct ccn_clamp_data *clamp_list = state->clamp_list;
       clamp_list != NULL;
       clamp_list = clamp_list->next) {
    clamp_list->n_sides = 0; 
  }

  return 0;
}

/***************************************************************************
destroy_walls:
  In:  state: MCell state
  Out: Free up the memory for the walls and vertices
***************************************************************************/
void destroy_walls(struct volume *state) {
  if (state->walls_using_vertex) {
    for (int i = 0; i<state->n_verts; i++) {
      delete_void_list((struct void_list *)state->walls_using_vertex[i]);
    }
  }

  free(state->walls_using_vertex);
  free(state->all_vertices);
  state->n_walls = 0;
  state->n_verts = 0;
}

/***************************************************************************
destroy_partitions:
  In:  state: MCell state
  Out: Free up the memory for the fine and coarse partitions
***************************************************************************/
void destroy_partitions(struct volume *state) {
  state->n_fineparts = 0;
  free(state->x_fineparts);
  free(state->y_fineparts);
  free(state->z_fineparts);
  state->x_fineparts = NULL;
  state->y_fineparts = NULL;
  state->z_fineparts = NULL;

  free(state->x_partitions);
  free(state->y_partitions);
  free(state->z_partitions);
  state->x_partitions = NULL;
  state->y_partitions = NULL;
  state->z_partitions = NULL;
}

/***************************************************************************
destroy_objects:
  In: obj_ptr: object to be destroyed
      free_poly_flag: see explanation in destroy_poly_object
  Out: Zero on success. One otherwise. Recursively destroys objects.
  Note: Currently, this ultimately only destroys polygon objects. I don't know
        if there's a need to trash release objects that use release patterns.
***************************************************************************/
int destroy_objects(struct object *obj_ptr, int free_poly_flag) {
  obj_ptr->sym->count = 0;
  switch (obj_ptr->object_type) {
  case META_OBJ:
    for (struct object *child_obj_ptr = obj_ptr->first_child;
         child_obj_ptr != NULL; child_obj_ptr = child_obj_ptr->next) {
      destroy_objects(child_obj_ptr, free_poly_flag);
      child_obj_ptr->n_walls = 0;
      child_obj_ptr->n_walls_actual = 0;
      child_obj_ptr->n_verts = 0;
      child_obj_ptr->first_child = NULL;
      child_obj_ptr->last_child = NULL;
    }
    break;
  case BOX_OBJ:
  case POLY_OBJ:
    destroy_poly_object(obj_ptr, free_poly_flag);
    break;

  // do nothing
  case REL_SITE_OBJ:
  case VOXEL_OBJ:
    break;
  }

  return 0;
}

/***************************************************************************
destroy_poly_object:
  In: obj_ptr: object
      free_poly_flag: Destroy polygon_object if set. There's a pointer to this
      in the object definition (child of root_object) AND the instantiated
      object (child of root_insance), so we only want to free it once.
  Out: Zero on success. One otherwise. Polygon object is destroyed
***************************************************************************/
int destroy_poly_object(struct object *obj_ptr, int free_poly_flag) {
  if (free_poly_flag) {
    for (int wall_num = 0; wall_num < obj_ptr->n_walls; wall_num++) {
      struct wall *w = obj_ptr->wall_p[wall_num];
      if (w->grid) {
        /*free(w->grid->mol);*/
        delete_void_list((struct void_list *)w->grid->sm_list);
      } 
      delete_void_list((struct void_list *)w->surf_class_head);
    }
    struct polygon_object *poly_obj_ptr = obj_ptr->contents;
    free(poly_obj_ptr->side_removed);
    poly_obj_ptr->side_removed = NULL;
    free(poly_obj_ptr->element);
    poly_obj_ptr->element = NULL;
    poly_obj_ptr->references--;
    // Clean up when there are no other instances of this object
    if (poly_obj_ptr->references == 0) {
      free(obj_ptr->contents);
      obj_ptr->contents = NULL;
    }
    free(obj_ptr->last_name);
    obj_ptr->last_name = NULL;
  }
  obj_ptr->sym->count = 0;
  free(obj_ptr->walls);
  free(obj_ptr->wall_p);
  obj_ptr->wall_p = NULL;
  free(obj_ptr->vertices);
  obj_ptr->vertices = NULL;

  obj_ptr->wall_p = NULL;
  obj_ptr->n_walls = 0;
  obj_ptr->n_walls_actual = 0;
  obj_ptr->n_verts = 0;
  struct region_list *regs, *next_regs;
  for (regs = obj_ptr->regions; regs != NULL;) {
    if (free_poly_flag && (strcmp(regs->reg->region_last_name, "ALL") != 0)) {
      //XXX: this doesn't work with dynamic geometries and pymcell...
      /*free(regs->reg->region_last_name);*/
      regs->reg->region_last_name = NULL;
    }
    delete_void_list((struct void_list *)regs->reg->sm_dat_head);
    regs->reg->sm_dat_head = NULL;
    free(regs->reg->membership);
    regs->reg->membership = NULL;
    free(regs->reg->bbox);
    regs->reg->bbox = NULL;
    if (regs->reg->boundaries) {
      pointer_hash_destroy(regs->reg->boundaries);
    }
    free(regs->reg->boundaries); 
    regs->reg->boundaries = NULL;
    regs->reg->sym->count = 0;
    next_regs = regs->next;
    free(regs);
    regs = next_regs;
  }
  obj_ptr->regions = NULL;
  obj_ptr->num_regions = 0;
  obj_ptr->total_area = 0;
  init_matrix(obj_ptr->t_matrix);

  return 0;
}

/***************************************************************************
reset_current_counts:
  In: mol_sym_table:
      count_hashmask:
      count_hash:
  Out: Zero on success. Species populations are set to zero. Counts on/in
       regions are also set to zero.
***************************************************************************/
int reset_current_counts(struct sym_table_head *mol_sym_table,
                         int count_hashmask,
                         struct counter **count_hash) {
  // Set global populations of species back to zero, since they will get set to
  // the proper values when we insert the molecules into the world
  for (int n_mol_bin = 0; n_mol_bin < mol_sym_table->n_bins; n_mol_bin++) {
    for (struct sym_entry *sym_ptr = mol_sym_table->entries[n_mol_bin];
         sym_ptr != NULL; sym_ptr = sym_ptr->next) {
      struct species *mol = (struct species *)sym_ptr->value;
      mol->population = 0;
    }
  }

  // Set counts on/in regions back to zero for the same reasons listed above.
  for (int i = 0; i <= count_hashmask; i++) {
    if (count_hash[i] != NULL) {
      struct counter *c;
      for (c = count_hash[i]; c != NULL; c = c->next) {
        if ((c->counter_type & MOL_COUNTER) != 0) {
          c->data.move.n_enclosed = 0; 
          c->data.move.n_at = 0; 
        }
      }
    }
  }
  return 0;
}

/***************************************************************************
check_count_validity:
  Check if all the regions that we are counting in/on are
  valid/defined/existent or invalid/undefined/non-existent. This is necessary
  since objects can appear or disappear with dynamic geometries.

  In: output_request_head:
      regions_to_ignore:
      new_region_names:
      meshes_to_ignore:
      new_mesh_names:
  Out: If a count type needs changed, it is done by reset_count_type.
***************************************************************************/
void check_count_validity(struct output_request *output_request_head,
                          struct string_buffer *regions_to_ignore,
                          struct string_buffer *new_region_names,
                          struct string_buffer *meshes_to_ignore,
                          struct string_buffer *new_mesh_names) {

  for (struct output_request *request = output_request_head;
       request != NULL; request = request->next) {
    if (request->count_location != NULL) {
      if (request->count_location->sym_type != REG) {
        mcell_internal_error(
            "Non-region location symbol (type=%d) in count request.",
            request->count_location->sym_type);
      }
      char *reg_name = request->count_location->name;
      // Counting in/on an object
      if (is_reverse_abbrev(",ALL", reg_name)) {
        struct region *reg_of_count = (
            struct region *)request->count_location->value;
        char *obj_name = reg_of_count->parent->sym->name;
        reset_count_type(
            obj_name,
            request,
            meshes_to_ignore,
            new_mesh_names);
      }
      // Counting in/on a region
      else {
        reset_count_type(
            reg_name,
            request,
            regions_to_ignore,
            new_region_names);
      }
    }
  }
}

/***************************************************************************
reset_count_type:
  We may need to either change the count type from valid->invalid or
  invalid->valid. For example, if we were counting all the molecules in an
  object that disappeared, we would need to change the type from COUNT_INT to
  UNSET.

  In: name: the name of the mesh/region
      request: information about count statement we want to reset
      names_to_ignore: a string buffer of meshes/regions to ignore
      new_names: a string buffer of meshes/regions that were just added
  Out: none
***************************************************************************/
void reset_count_type(char *name,
                      struct output_request *request,
                      struct string_buffer *names_to_ignore,
                      struct string_buffer *new_names) {

  int num_ignore = names_to_ignore->n_strings;
  int num_add = new_names->n_strings;
  struct output_buffer *buffer = request->requester->column->buffer;
  struct output_block *block = request->requester->column->set->block;
  int buf_index = block->buf_index;
  int buffersize = block->buffersize;
  int trig_bufsize = block->trig_bufsize;
  // Reset count type that were potentially unset. This should be more
  // efficient.
  if (is_string_present_in_string_array(
      name, new_names->strings, num_add)) {
    // XXX: We can't assume that a count which was originally unset should
    // now be an int.
    if (buffer[0].data_type == COUNT_UNSET) {
      for (int idx = buf_index; idx < buffersize; idx++) {
        buffer[idx].data_type = COUNT_INT;
      }
    }
    else if (buffer[0].data_type == COUNT_TRIG_STRUCT) {
      for (int idx = buf_index; idx < trig_bufsize; idx++) {
        buffer[idx].data_type = COUNT_TRIG_STRUCT;
      }
    }
    else if (buffer[0].data_type == COUNT_INT) {
      for (int idx = buf_index; idx < buffersize; idx++) {
        buffer[idx].data_type = COUNT_INT;
      }
    }
    else if (buffer[0].data_type == COUNT_DBL) {
      for (int idx = buf_index; idx < buffersize; idx++) {
        buffer[idx].data_type = COUNT_DBL;
      }
    }
  }
  // Unset count types for meshes that were removed.
  else if (is_string_present_in_string_array(
      name, names_to_ignore->strings, num_ignore)) {
    if (buffer[0].data_type == COUNT_TRIG_STRUCT) {
      for (int idx = buf_index; idx < trig_bufsize; idx++) {
        buffer[idx].val.tval->name = NULL;
      }
    }
    else {
      for (int idx = buf_index; idx < buffersize; idx++) {
        buffer[idx].data_type = COUNT_UNSET;
      }
    }
  }
}

/***************************************************************************
init_species_mesh_transp:
  In:  state: MCell state
  Out: Zero on success. Create a data structure so we can quickly check if a
       molecule species can move in or out of any given surface region
***************************************************************************/
int init_species_mesh_transp(struct volume *state) {
  struct pointer_hash *species_mesh_transp;
  if ((species_mesh_transp = CHECKED_MALLOC_STRUCT(
      struct pointer_hash, "pointer_hash")) == NULL) {
    mcell_internal_error("Out of memory while creating molecule-object "
                         "transparency pointer hash");
  }
  if (pointer_hash_init(species_mesh_transp, state->n_species)) {
    mcell_error(
      "Failed to initialize data structure for molecule-object transparency.");
  }
  // Initialize pointer hash with species names as keys and values are a linked
  // list of pointers of mesh_transparency. The mesh_transparency struct
  // contains the mesh name and whether the species can go in_to_out and/or
  // out_to_in. These are initially set to 0 (not transparent), and can be set
  // to 1 (transparent) in find_vm_obj_region_transp.
  state->species_mesh_transp = species_mesh_transp;
  for (int i = 0; i < state->n_species; i++) {
    struct species *spec = state->species_list[i];
    char *species_name = spec->sym->name;
    if (spec->flags & IS_SURFACE)
      continue;
    if (strcmp(species_name, "ALL_MOLECULES") == 0)
      continue;
    else if (strcmp(species_name, "ALL_VOLUME_MOLECULES") == 0)
      continue;
    else if (strcmp(species_name, "ALL_SURFACE_MOLECULES") == 0)
      continue;
    unsigned int keyhash = (unsigned int)(intptr_t)(species_name);
    void *key = (void *)(species_name);
    struct mesh_transparency *mesh_transp_head = NULL;
    struct mesh_transparency *mesh_transp_tail = NULL;
    int sm_flag = 0;
    if (spec->flags & ON_GRID) {
      sm_flag = 1;
    }
    find_all_obj_region_transp(state->root_instance, &mesh_transp_head,
                               &mesh_transp_tail, species_name, sm_flag);
    if (pointer_hash_add(
        state->species_mesh_transp, key, keyhash, (void *)mesh_transp_head)) {
      mcell_allocfailed("Failed to store species-mesh transparency in"
                        " pointer_hash table.");
    }
  }
  return 0;
}

/***************************************************************************
find_sm_region_transp:
  In: obj_ptr: the mesh object
      mesh_transp_head: Head of the mesh transparency list
      mesh_transp_tail: Tail of the mesh transparency list
      species_name: the species/molecule name
  Out: Zero on success. Create a data structure so we can quickly check if a
  surface molecule species can move in or out of any given surface region
***************************************************************************/
int find_sm_region_transp(struct object *obj_ptr,
                          struct mesh_transparency **mesh_transp_head,
                          struct mesh_transparency **mesh_transp_tail,
                          char *species_name) {

  // Check every region on the object
  for (struct region_list *reg_list_ptr = obj_ptr->regions;
       reg_list_ptr != NULL;
       reg_list_ptr = reg_list_ptr->next) {
    struct region *reg_ptr = reg_list_ptr->reg;
    if (reg_ptr->surf_class != NULL) {
      // Unlike volume molecules, we need to create a mesh transparency entry
      // for every region instead of every object.
      struct mesh_transparency *mesh_transp;
      mesh_transp = CHECKED_MALLOC_STRUCT(
          struct mesh_transparency, "object transparency");
      mesh_transp->next = NULL;
      mesh_transp->name = reg_ptr->sym->name;
      // Ignore this until I merge in Markus' experimental changes
      mesh_transp->in_to_out = 0;
      mesh_transp->out_to_in = 0;
      // Default state for surface molecules is transparent, so we just
      // need to check reflective and absorptive regions
      mesh_transp->transp_top_front = 1;
      mesh_transp->transp_top_back = 1;
      if (*mesh_transp_tail == NULL) {
        *mesh_transp_head = mesh_transp; 
        *mesh_transp_tail = mesh_transp; 
        (*mesh_transp_head)->next = NULL;
        (*mesh_transp_tail)->next = NULL;
      }
      else {
        (*mesh_transp_tail)->next = mesh_transp;
        *mesh_transp_tail = mesh_transp;
      }
      // Check if species_name is in the absorptive list for this region
      check_surf_class_properties(
          species_name, mesh_transp, reg_ptr->surf_class->absorb_mols);
      // Check if species_name is in the relective list for this region
      check_surf_class_properties(
          species_name, mesh_transp, reg_ptr->surf_class->refl_mols);
    }
  }
  return 0;
}

/***************************************************************************
check_surf_class_properties:
  In:  species_name: the name of the molecule/species that we are checking
       mesh_transp: contains info about whether the species is transparent to
                    the regions on this mesh
       surf_class_props: absorptive/reflective surface class properties
  Out: None. Check if species_name is in the absorptive/reflective list for
       a given region
***************************************************************************/
void check_surf_class_properties(
  char *species_name,
  struct mesh_transparency *mesh_transp,
  struct name_orient *surf_class_props) {

  struct name_orient *no;
  for (no = surf_class_props; no != NULL; no = no->next) {
    if (strcmp(no->name, species_name) == 0) {
      // Absorptive/Reflective to top front molecules
      if (no->orient == 1) {
        mesh_transp->transp_top_front = 0;
        break;
      }
      // Absorptive/Reflective to top back molecules
      else if (no->orient == -1) {
        mesh_transp->transp_top_back = 0;
        break;
      }
      // Absorptive/Reflective from top front or top back (e.g. ABSORPTIVE = A;)
      else if (no->orient == 0) {
        mesh_transp->transp_top_front = 0;
        mesh_transp->transp_top_back = 0;
        break;
      }
    }
  }
}

/***************************************************************************
find_vm_obj_region_transp:
  In:  obj_ptr: The object we are currently checking for transparency
       mesh_transp_head: Head of the mesh transparency list
       mesh_transp_tail: Tail of the mesh transparency list
       species_name: The name of the molecule/species we are checking
  Out: Zero on success. Check every region on obj_ptr to see if any of them are
       transparent to the volume molecules with species_name.
***************************************************************************/
int find_vm_obj_region_transp(struct object *obj_ptr,
                              struct mesh_transparency **mesh_transp_head,
                              struct mesh_transparency **mesh_transp_tail,
                              char *species_name) {

  struct mesh_transparency *mesh_transp;
  mesh_transp =
      CHECKED_MALLOC_STRUCT(struct mesh_transparency, "object transparency");
  mesh_transp->next = NULL;
  mesh_transp->name = obj_ptr->sym->name;
  mesh_transp->in_to_out = 0;
  mesh_transp->out_to_in = 0;
  if (*mesh_transp_tail == NULL) {
    *mesh_transp_head = mesh_transp; 
    *mesh_transp_tail = mesh_transp; 
    (*mesh_transp_head)->next = NULL;
    (*mesh_transp_tail)->next = NULL;
  }
  else {
    (*mesh_transp_tail)->next = mesh_transp;
    *mesh_transp_tail = mesh_transp;
  }
  // Assuming the first region in the region list will always be ALL.
  // Could be unsafe.
  double volume = obj_ptr->regions->reg->volume;
  // Set volume for the ALL reg if it hasn't been done already. Bit clunky...
  if (!distinguishable(volume, 0.0, EPS_C)) {
    int count_regions_flag = 0;
    is_manifold(obj_ptr->regions->reg, count_regions_flag);
    volume = obj_ptr->regions->reg->volume;
  }
  // Check every region on the object
  for (struct region_list *reg_list_ptr = obj_ptr->regions;
       reg_list_ptr != NULL;
       reg_list_ptr = reg_list_ptr->next) {
    struct region *reg_ptr = reg_list_ptr->reg;
    if (reg_ptr->surf_class != NULL) {
      struct name_orient *no;
      // Check if species_name is in the transparency list for this region
      for (no = reg_ptr->surf_class->transp_mols; no != NULL; no = no->next) {
        if (strcmp(no->name, species_name) == 0) {
          // Side note about reg_ptr->volume stuff below: if volume is
          // positive, then we have outward facing normals. this is the typical
          // case. if volume is negative, then you have inward facing normals.

          // Transparent from outside to inside
          if ((no->orient == 1 && volume > 0) ||
              (no->orient == -1 && volume < 0)) {
            mesh_transp->out_to_in = 1;
            break;
          }
          // Transparent from inside to outside
          else if ((no->orient == -1 && volume > 0) ||
              (no->orient == 1 && volume < 0)) {
            mesh_transp->in_to_out = 1;
            break;
          }
          // Transparent from either direction (e.g. TRANSPARENT = A;)
          else if (no->orient == 0) {
            mesh_transp->in_to_out = 1;
            mesh_transp->out_to_in = 1;
            break;
          }
        }
      }
    }
  }
  return 0;
}

/***************************************************************************
find_all_obj_region_transp:
  In: obj_ptr: The mesh object
      mesh_transp_head: Head of the object transparency list
      mesh_transp_tail: Tail of the object transparency list
      species_name: The name of the molecule/species we are checking
      sm_flag: surface molecule flag
  Out: Zero on success. Check every polygon object to see if it is transparent
       to species_name.
***************************************************************************/
int find_all_obj_region_transp(struct object *obj_ptr,
                               struct mesh_transparency **mesh_transp_head,
                               struct mesh_transparency **mesh_transp_tail,
                               char *species_name,
                               int sm_flag) {

  switch (obj_ptr->object_type) {
  case META_OBJ:
    for (struct object *child_obj_ptr = obj_ptr->first_child;
         child_obj_ptr != NULL; child_obj_ptr = child_obj_ptr->next) {
      if (find_all_obj_region_transp(
          child_obj_ptr, mesh_transp_head, mesh_transp_tail, species_name,
          sm_flag))
        return 1;
    }
    break;

  case BOX_OBJ:
  case POLY_OBJ:
    if (sm_flag) {
      if (find_sm_region_transp(
          obj_ptr, mesh_transp_head, mesh_transp_tail, species_name)) {
        return 1;
      }
    }
    else {
      if (find_vm_obj_region_transp(
          obj_ptr, mesh_transp_head, mesh_transp_tail, species_name)) {
        return 1;
      }
    }
    break;

  case VOXEL_OBJ:
  case REL_SITE_OBJ:
    break;
  }

  return 0;
}

/************************************************************************
 add_dynamic_geometry_events:
 In:  dynamic_geometry_filename: filename and path for dyngeom file
      dynamic_geometry_filepath: XXX: identical to above? remove?
      timestep: global timestep in seconds
      dynamic_geometry_events_mem: memory to store time and MDL names for
                                   dynamic geometry
      dg_time_fname_head: the head of the dynamic geometry event list
 Out: 0 on success, 1 on failure. dynamic geometry events are added to
      dg_time_fname_head from which they will eventually be added to a
      scheduler.
 ***********************************************************************/
int add_dynamic_geometry_events(
    struct mdlparse_vars *parse_state,
    char *dynamic_geometry_filepath,
    double timestep,
    struct mem_helper *dynamic_geometry_events_mem,
    struct dg_time_filename **dg_time_fname_head) {

  struct volume *state = parse_state->vol;
  struct dyngeom_parse_vars *dg_parse = create_dg_parse(state);
  state->dg_parse = dg_parse;
  FILE *f = fopen(dynamic_geometry_filepath, "r");

  if (!f) {
    return 1;
  } else {
    const char *SEPARATORS = "\f\n\r\t\v ,;";
    const char *FIRST_DIGIT = "+-0123456789";
    struct dg_time_filename *dg_time_fname_tail = NULL;
    char buf[2048];
    char *char_ptr;
    char *zero_file_name = NULL;
    int linecount = 0;
    int i;

    while (fgets(buf, 2048, f)) {
      linecount++;
      // Ignore leading whitespace
      for (i = 0; i < 2048; i++) {
        if (!strchr(SEPARATORS, buf[i]))
          break;
      }

      if (i < 2048 && strchr(FIRST_DIGIT, buf[i])) {
        // Grab time
        double time = strtod((buf + i), &char_ptr);
        if (char_ptr == (buf + i))
          continue; /* Conversion error. */

        // Skip over whitespace between time and filename
        for (i = char_ptr - buf; i < 2048; i++) {
          if (!strchr(SEPARATORS, buf[i]))
            break;
        }

        // Grab mdl filename. This could probably be cleaned up
        char *line_ending = strchr(buf + i, '\n');
        int line_ending_idx = line_ending - buf;
        int file_name_length = line_ending_idx - i + 1;
        char file_name[file_name_length];
        strncpy(file_name, buf + i, file_name_length);
        file_name[file_name_length - 1] = '\0';
        // Expand path name if needed
        char *full_file_name = mcell_find_include_file(
          file_name, dynamic_geometry_filepath);
        // Treat time 0 as if it were an include file.
        if (!distinguishable(time, 0, EPS_C)) {
          zero_file_name = full_file_name;
          continue;
        }
        // Do the normal DG parsing on every other time
        else {
          parse_dg_init(dg_parse, full_file_name, state);
          destroy_objects(state->root_instance, 0);
          destroy_objects(state->root_object, 0);

          struct dg_time_filename *dyn_geom;
          dyn_geom = CHECKED_MEM_GET(dynamic_geometry_events_mem,
                                     "time-varying dynamic geometry");
          if (dyn_geom == NULL)
            return 1;

          dyn_geom->event_time = round(time / timestep);
          dyn_geom->mdl_file_path = full_file_name;
          dyn_geom->next = NULL;

          // Append each entry to end of dg_time_fname_head list
          if (*dg_time_fname_head == NULL) {
            *dg_time_fname_head = dyn_geom;
            dg_time_fname_tail = dyn_geom;
          } else {
            dg_time_fname_tail->next = dyn_geom;
            dg_time_fname_tail = dyn_geom;
          }
        }
      }
    }

    fclose(f);
    parse_state->current_object = parse_state->vol->root_object;
#ifdef NOSWIG
    if (zero_file_name && mdlparse_file(parse_state, zero_file_name))
    {
      free(zero_file_name);
      return 1;
    }
#endif
    free(zero_file_name);
  }

  free(dynamic_geometry_filepath);
  // Disable parsing of geometry for the rest of the MDL. It should only happen
  // via files referenced in the DG file.
  state->disable_polygon_objects = 1;
  return 0;
}

/************************************************************************
 get_mesh_instantiation_names:
 In:  obj_ptr: the root instance meta object
      mesh_names: an initialized but empty string buffer
 Out: mesh_names is updated so that it contains a list of all the mesh objects
      with their fully qualified names.
 ***********************************************************************/
char *get_mesh_instantiation_names(struct object *obj_ptr,
                                   struct string_buffer *mesh_names) {
  switch (obj_ptr->object_type) {
  case META_OBJ:
    for (struct object *child_obj_ptr = obj_ptr->first_child;
         child_obj_ptr != NULL; child_obj_ptr = child_obj_ptr->next) {
      char *mesh_name = get_mesh_instantiation_names(
          child_obj_ptr, mesh_names);
      if ((mesh_name != NULL) &&
          (add_string_to_buffer(mesh_names, mesh_name))) {
        free(mesh_name);
        destroy_string_buffer(mesh_names);
        return NULL;
      }
    }
    break;
  case BOX_OBJ:
  case POLY_OBJ:
    return CHECKED_STRDUP(obj_ptr->sym->name, "mesh name");

  // do nothing
  case REL_SITE_OBJ:
  case VOXEL_OBJ:
    break;
  }
  return NULL;
}

/************************************************************************
 diff_string_buffers:
 In:  diff_names: The new names are stored here
      names_a: One set of names
      names_b: Another set of names
 Out: Assign difference of names_a and names_b to diff_names.
      Example: {A,B,C} - {B,C,D} = {A}. There might not be any if the the old
      and new list are identical. I'm sure this could be much more efficient,
      but it is sufficient for now. Note: This is very similar to
      sym_diff_string_buffers. Consolidate these.
 ***********************************************************************/
void diff_string_buffers(
    struct string_buffer *diff_names,
    struct string_buffer *names_a,
    struct string_buffer *names_b) {

  int n_strings_old = names_a->n_strings;
  int n_strings_new = names_b->n_strings;
  for (int i = 0; i < n_strings_old; i++) {
    if (!(is_string_present_in_string_array(
        names_a->strings[i], names_b->strings, n_strings_new))) {
      char *diff_name = CHECKED_STRDUP(names_a->strings[i], "name");
      if (add_string_to_buffer(diff_names, diff_name)) {
        destroy_string_buffer(diff_names);
      }
    }
  }
}

/************************************************************************
 sym_diff_string_buffers:
 In:  diff_names: The new names are stored here
      names_a: One set of names
      names_b: Another set of names
 Out: Assign symmetric difference of names_a and names_b to diff_names.
      Example: {A,B,C} + {B,C,D} = {A,D}. There might not be any if the the old
      and new list are identical. I'm sure this could be much more efficient,
      but it is sufficient for now.
 ***********************************************************************/
void sym_diff_string_buffers(
    struct string_buffer *diff_names,
    struct string_buffer *names_a,
    struct string_buffer *names_b,
    enum warn_level_t add_remove_mesh_warning) {

  int n_strings_old = names_a->n_strings;
  int n_strings_new = names_b->n_strings;

  // Track objects which have been removed
  for (int i = 0; i < n_strings_old; i++) {
    if (!(is_string_present_in_string_array(
        names_a->strings[i], names_b->strings, n_strings_new))) {
      char *diff_name = CHECKED_STRDUP(names_a->strings[i], "name");

      switch (add_remove_mesh_warning) {
        case WARN_COPE:
          break;
        case WARN_WARN:
          mcell_warn("\"%s\" removed.", diff_name);
          break;
        case WARN_ERROR:
          mcell_error("\"%s\" removed.", diff_name);
      }

      if (add_string_to_buffer(diff_names, diff_name)) {
        destroy_string_buffer(diff_names);
      }
    }
  }

  // Track objects which have been added
  for (int i = 0; i < n_strings_new; i++) {
    if (!(is_string_present_in_string_array(
        names_b->strings[i], names_a->strings, n_strings_old))) {
      char *diff_name = CHECKED_STRDUP(names_b->strings[i], "name");

      switch (add_remove_mesh_warning) {
        case WARN_COPE:
          break;
        case WARN_WARN:
          mcell_warn("\"%s\" added.", diff_name);
          break;
        case WARN_ERROR:
          mcell_error("\"%s\" added.", diff_name);
      }

      if (add_string_to_buffer(diff_names, diff_name)) {
        destroy_string_buffer(diff_names);
      }
    }
  }

}

/***************************************************************************
get_reg_names_all_objects:
  In: obj_ptr: grab the region names of this object (if it's a poly object) or
               its children (if it's a meta object) 
      region_names: we will store all of the region names in this
  Out: 0 on success, 1 on failure.
***************************************************************************/
int get_reg_names_all_objects(
    struct object *obj_ptr,
    struct string_buffer *region_names) {

  switch (obj_ptr->object_type) {
  case META_OBJ:
    for (struct object *child_obj_ptr = obj_ptr->first_child;
         child_obj_ptr != NULL; child_obj_ptr = child_obj_ptr->next) {
      if (get_reg_names_all_objects(child_obj_ptr, region_names))
        return 1;
    }
    break;

  case BOX_OBJ:
  case POLY_OBJ:
    if (get_reg_names_this_object(obj_ptr, region_names))
      return 1;
    break;

  case VOXEL_OBJ:
  case REL_SITE_OBJ:
    break;
  }

  return 0;
}

/***************************************************************************
get_reg_names_this_object:
  In: obj_ptr: we will get region names from this object
      region_names: we will store all of the region names in this
  Out: 0 on success, 1 on failure.
***************************************************************************/
int get_reg_names_this_object(
    struct object *obj_ptr,
    struct string_buffer *region_names) {
  for (struct region_list *reg_list_ptr = obj_ptr->regions;
       reg_list_ptr != NULL;
       reg_list_ptr = reg_list_ptr->next) {

    struct region *reg_ptr = reg_list_ptr->reg;

    // We only care about regions with boundaries so ignore these
    if ((strcmp(reg_ptr->region_last_name, "ALL") == 0) ||
        (reg_ptr->region_has_all_elements))
      continue;

    char *region_name = CHECKED_STRDUP(reg_ptr->sym->name, "region name");
    if ((region_name != NULL) &&
        (add_string_to_buffer(region_names, region_name))) {
      free(region_name);
      destroy_string_buffer(region_names);
      return 1;
    }
  }
  return 0;
}

/***************************************************************************
update_geometry:
  In:  state: MCell state
       dyn_geom: info about next dyngeom event (time and geom filename)
  Out: None. Molecule positions are saved. Old geometry is trashed. New
       geometry is created. Molecules are placed (and moved if necessary).
***************************************************************************/
void update_geometry(struct volume *state,
                     struct dg_time_filename *dyn_geom) {
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

  state->mdl_infile_name = dyn_geom->mdl_file_path;
  if (mcell_redo_geom(state)) {
    mcell_error("An error occurred while processing geometry changes.");
  }

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

}
