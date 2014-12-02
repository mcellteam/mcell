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

#define NO_MESH "\0"

/***************************************************************************
 save_all_molecules: Save all the molecules currently in the scheduler.

 In:  state: MCell state
      storage_head:
 Out: An array of all the molecules to be saved
***************************************************************************/
struct molecule_info **save_all_molecules(struct volume *state,
                                          struct storage_list *storage_head) {

  // Find total number of molecules in the scheduler.
  unsigned long long total_items = count_items_in_scheduler(storage_head);
  int ctr = 0;
  struct molecule_info **all_molecules;
  all_molecules = CHECKED_MALLOC_ARRAY(struct molecule_info *, total_items,
                                       "all molecules");

  // Iterate over all molecules in the scheduler.
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

          const int MAX_NUM_REGIONS = 100;
          struct string_buffer *reg_names =
              CHECKED_MALLOC_STRUCT(struct string_buffer, "string buffer");
          if (initialize_string_buffer(reg_names, MAX_NUM_REGIONS)) {
            return NULL;
          }

          char *mesh_name;
          if ((am_ptr->properties->flags & NOT_FREE) == 0) {
            save_volume_molecule(state, mol_info, am_ptr, &mesh_name);
          } else if ((am_ptr->properties->flags & ON_GRID) != 0) {
            if (save_surface_molecule(mol_info, am_ptr, &reg_names, &mesh_name))
              return NULL;
          } else {
            continue;
          }

          save_common_molecule_properties(mol_info, am_ptr, reg_names,
                                          mesh_name);
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

 In:  mol_info:
      am_ptr: abstract molecule pointer
      reg_names: region names
      mesh_name: mesh name that molecule is on or in
 Out: Nothing. The common properties of surface and volume molecules are saved
      in mol_info.
***************************************************************************/
void save_common_molecule_properties(struct molecule_info *mol_info,
                                     struct abstract_molecule *am_ptr,
                                     struct string_buffer *reg_names,
                                     char *mesh_name) {
  mol_info->molecule->t = am_ptr->t;
  mol_info->molecule->t2 = am_ptr->t2;
  mol_info->molecule->flags = am_ptr->flags;
  mol_info->molecule->properties = am_ptr->properties;
  mol_info->molecule->birthday = am_ptr->birthday;
  mol_info->molecule->mesh_name = CHECKED_STRDUP(mesh_name, "mesh name");
  // Only free temporary object names we just allocated above.
  // Don't want to accidentally free symbol names of objects.
  if ((mesh_name != NO_MESH) && ((am_ptr->properties->flags & NOT_FREE) == 0)) {
    free(mesh_name);
  }
  mol_info->reg_names = reg_names;
}

/***************************************************************************
 save_volume_molecule:

 In:  state: MCell state
      mol_info:
      am_ptr: abstract molecule pointer
      mesh_name: mesh name that molecule is in
 Out: Nothing. Molecule info and mesh name are updated
***************************************************************************/
void save_volume_molecule(struct volume *state, struct molecule_info *mol_info,
                          struct abstract_molecule *am_ptr, char **mesh_name) {
  struct volume_molecule *vm_ptr = (struct volume_molecule *)am_ptr;

  int farthest_flag = 0;
  *mesh_name = find_enclosing_mesh_name(state, vm_ptr, farthest_flag);
  if (*mesh_name == NULL) {
    *mesh_name = NO_MESH;
  }
  mol_info->pos.x = vm_ptr->pos.x;
  mol_info->pos.y = vm_ptr->pos.y;
  mol_info->pos.z = vm_ptr->pos.z;
  mol_info->orient = 0;
}

/***************************************************************************
 save_surface_molecule:

 In:  mol_info:
      am_ptr: abstract molecule pointer
      reg_names: region names
      mesh_name: mesh name that molecule is on
 Out: Zero on success. One otherwise.
***************************************************************************/
int save_surface_molecule(struct molecule_info *mol_info,
                          struct abstract_molecule *am_ptr,
                          struct string_buffer **reg_names, char **mesh_name) {
  struct vector3 where;
  struct surface_molecule *sm_ptr = (struct surface_molecule *)am_ptr;
  uv2xyz(&sm_ptr->s_pos, sm_ptr->grid->surface, &where);
  mol_info->pos.x = where.x;
  mol_info->pos.y = where.y;
  mol_info->pos.z = where.z;
  mol_info->orient = sm_ptr->orient;
  *mesh_name = sm_ptr->grid->surface->parent_object->sym->name;
  int num_regions;
  struct name_list *reg_name_list_head, *reg_name_list;
  reg_name_list_head =
      find_regions_names_by_wall(sm_ptr->grid->surface, &num_regions);
  int k;
  for (reg_name_list = reg_name_list_head, k = 0; reg_name_list != NULL;
       reg_name_list = reg_name_list->next, k++) {
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
  return 0;
}

/***************************************************************************
 place_all_molecules: Place all molecules currently in the scheduler into the
                      world.

 In:  state: MCell state
 Out: Zero on success. One otherwise.
***************************************************************************/
int place_all_molecules(struct volume *state) {

  struct volume_molecule vm;
  memset(&vm, 0, sizeof(struct volume_molecule));
  struct volume_molecule *vm_ptr = &vm;
  struct volume_molecule *vm_guess = NULL;

  unsigned long long total_items = state->num_all_molecules;

  for (unsigned long long n_mol = 0; n_mol < total_items; n_mol++) {

    struct molecule_info *mol_info = state->all_molecules[n_mol];
    struct abstract_molecule *am_ptr = mol_info->molecule;
    char *mesh_name = am_ptr->mesh_name;
    // Insert volume molecule into world.
    if ((am_ptr->properties->flags & NOT_FREE) == 0) {
      vm_ptr->t = am_ptr->t;
      vm_ptr->t2 = am_ptr->t2;
      vm_ptr->flags = am_ptr->flags;
      vm_ptr->properties = am_ptr->properties;
      vm_ptr->birthday = am_ptr->birthday;
      vm_ptr->pos.x = mol_info->pos.x;
      vm_ptr->pos.y = mol_info->pos.y;
      vm_ptr->pos.z = mol_info->pos.z;

      if (strlen(mesh_name) == 0) {
        vm_guess =
            insert_volume_molecule_encl_mesh(state, vm_ptr, vm_guess, NULL);
      } else {
        vm_guess = insert_volume_molecule_encl_mesh(state, vm_ptr, vm_guess,
                                                    mesh_name);
      }

      if (vm_guess == NULL) {
        mcell_error("Cannot insert copy of molecule of species '%s' into "
                    "world.\nThis may be caused by a shortage of memory.",
                    vm_ptr->properties->sym->name);
      }
    }
        // Insert surface molecule into world.
        else if ((am_ptr->properties->flags & ON_GRID) != 0) {
      insert_surface_molecule(state, am_ptr->properties, &mol_info->pos,
                              mol_info->orient, CHKPT_GRID_TOLERANCE, am_ptr->t,
                              NULL, mesh_name, mol_info->reg_names);
    }
  }

  // Do some cleanup.
  for (int i = 0; i < state->num_all_molecules; i++) {
    if (state->all_molecules[i]->molecule->mesh_name != NO_MESH) {
      free(state->all_molecules[i]->molecule->mesh_name);
    }
    free(state->all_molecules[i]->molecule);
    destroy_string_buffer(state->all_molecules[i]->reg_names);
    free(state->all_molecules[i]->reg_names);
    free(state->all_molecules[i]);
  }
  free(state->all_molecules);

  return 0;
}

/*************************************************************************
insert_volume_molecule_encl_mesh:
  In: state: MCell state
      vm: pointer to volume_molecule that we're going to place in local storage
      vm_guess: pointer to a volume_molecule that may be nearby
      mesh_name: closest enclosing mesh name
  Out: pointer to the new volume_molecule (copies data from volume molecule
       passed in), or NULL if out of memory.  Molecule is placed in scheduler
       also.
*************************************************************************/

struct volume_molecule *insert_volume_molecule_encl_mesh(
    struct volume *state, struct volume_molecule *vm,
    struct volume_molecule *vm_guess, char *mesh_name) {
  struct volume_molecule *new_vm;
  struct subvolume *sv, *new_sv;
  char *mesh_name_try;
  int move_molecule = 0;
  struct vector3 new_pos;

  if (vm_guess == NULL)
    sv = find_subvolume(state, &(vm->pos), NULL);
  else if (inside_subvolume(&(vm->pos), vm_guess->subvol, state->x_fineparts,
                            state->y_fineparts, state->z_fineparts)) {
    sv = vm_guess->subvol;
  } else
    sv = find_subvolume(state, &(vm->pos), vm_guess->subvol);

  new_vm = CHECKED_MEM_GET(sv->local_storage->mol, "volume molecule");
  memcpy(new_vm, vm, sizeof(struct volume_molecule));
  new_vm->mesh_name = NULL;
  new_vm->prev_v = NULL;
  new_vm->next_v = NULL;
  new_vm->next = NULL;
  new_vm->subvol = sv;

  int farthest_flag = 0;
  mesh_name_try = find_enclosing_mesh_name(state, new_vm, farthest_flag);

  //char *species_name = new_vm->properties->sym->name;
  //unsigned int keyhash = (unsigned int)(intptr_t)(species_name);
  //void *key = (void *)(species_name);
  //struct object_transparency *obj_transp = (
  //    struct object_transparency *)pointer_hash_lookup(state->mol_obj_transp,
  //                                                     key, keyhash);

  /* mol was inside mesh, now it is outside mesh */
  if ((mesh_name_try == NULL) && (mesh_name != NULL)) {
    move_molecule = 1;
  }
  /* mol was outside mesh, now it is inside mesh */
  if ((mesh_name_try != NULL) && (mesh_name == NULL)) {
    move_molecule = 1;
    free(mesh_name_try);
  }
  if ((mesh_name_try != NULL) && (mesh_name != NULL)) {
    /* mol was inside one mesh, now it is inside another mesh */
    if (strcmp(mesh_name_try, mesh_name) != 0)
      move_molecule = 1;
    free(mesh_name_try);
  }

  if (move_molecule) {
    /* move molecule to another location so that closest
       enclosing mesh name is "mesh_name" */

    place_mol_relative_to_mesh(state, &(vm->pos), sv, mesh_name, &new_pos);
    check_for_large_molecular_displacement(
        &(vm->pos), &new_pos, vm, &(state->time_unit),
        state->notify->large_molecular_displacement);
    new_vm->pos = new_pos;
    new_sv = find_subvolume(state, &(new_vm->pos), NULL);
    new_vm->subvol = new_sv;
  }

  new_vm->birthplace = new_vm->subvol->local_storage->mol;
  ht_add_molecule_to_list(&(new_vm->subvol->mol_by_species), new_vm);
  new_vm->subvol->mol_count++;
  new_vm->properties->population++;

  if ((new_vm->properties->flags & COUNT_SOME_MASK) != 0) {
    new_vm->flags |= COUNT_ME;
  }
  if (new_vm->properties->flags & (COUNT_CONTENTS | COUNT_ENCLOSED)) {
    count_region_from_scratch(state, (struct abstract_molecule *)new_vm, NULL,
                              1, &(new_vm->pos), NULL, new_vm->t);
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
       time_unit:
       large_molecular_displacement_warning: the warning value 
          (ignore, warn, error) set for molecular displacement
  Out: 0 on success, 1 otherwise.
************************************************************************/
void check_for_large_molecular_displacement(
    struct vector3 *old_pos,
    struct vector3 *new_pos,
    struct volume_molecule *vm,
    double *time_unit,
    enum warn_level_t large_molecular_displacement_warning) {

  double displacement = distance_vec3(old_pos, new_pos) / 100.0;
  double l_perp_bar = sqrt(4 * 1.0e8 * vm->properties->D * *time_unit / MY_PI);
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

    default:
      UNHANDLED_CASE(large_molecular_displacement_warning);
    }
  }
}

/*************************************************************************
find_enclosing_mesh_name:
  In:  state: MCell state
       vm: volume molecule
  Out: Name of the closest enclosing mesh if such exists,
       NULL otherwise.
************************************************************************/
char *find_enclosing_mesh_name(struct volume *state,
                               struct volume_molecule *vm,
                               int farthest_flag) {
  struct collision *smash; /* Thing we've hit that's under consideration */
  struct collision *shead; /* Head of the linked list of collisions */

  struct subvolume *sv;
  struct wall *w;
  struct vector3 rand_vector, world_diag;
  struct volume_molecule virt_mol; /* volume_molecule template */
  double world_diag_length; /* length of the world bounding box diagonal */

  /* we will reuse this struct in the different context of
     registering the meshes names through "no->name" and
     the number of hits with the mesh through "no->orient" */
  struct name_orient *no, *nol, *nnext, *no_head = NULL, *tail = NULL;
  char *return_name = NULL, *farthest_name = NULL;

  memcpy(&virt_mol, vm, sizeof(struct volume_molecule));
  virt_mol.prev_v = NULL;
  virt_mol.next_v = NULL;
  virt_mol.next = NULL;

  int calculate_random_vector = 1; /* flag */
  int found;                       /* flag */

pretend_to_call_find_enclosing_mesh: /* Label to allow fake recursion */

  sv = virt_mol.subvol;
  shead = NULL;

  /* pick up a random vector */
  if (calculate_random_vector) {
    srand((unsigned int)time(NULL));
    rand_vector.x = (double)rand() / (double)RAND_MAX;
    rand_vector.y = (double)rand() / (double)RAND_MAX;
    rand_vector.z = (double)rand() / (double)RAND_MAX;

    /* find the diagonal of the world */
    vectorize(&state->bb_urb, &state->bb_llf, &world_diag);
    world_diag_length = vect_length(&world_diag);

    /* scale random vector by the size of the world */
    rand_vector.x *= world_diag_length;
    rand_vector.y *= world_diag_length;
    rand_vector.z *= world_diag_length;
  }

  do {
    shead = ray_trace(state, &virt_mol, NULL, sv, &rand_vector, NULL);
    if (shead == NULL)
      mcell_internal_error("ray_trace() returned NULL.");

    if (shead->next != NULL) {
      shead =
          (struct collision *)ae_list_sort((struct abstract_element *)shead);
    }

    for (smash = shead; smash != NULL; smash = smash->next) {
      if ((smash->what & COLLIDE_WALL) != 0) {
        w = (struct wall *)smash->target;

        /* discard open-type meshes, like planes, etc. */
        if (w->parent_object->is_closed <= 0)
          continue;

        /* discard the cases when the random vector just grazes
           the mesh at the encounter point */
        double d_prod = dot_prod(&rand_vector, &(w->normal));
        if (!distinguishable(d_prod, 0, EPS_C))
          continue;

        if (no_head == NULL) {
          no = CHECKED_MALLOC_STRUCT(struct name_orient, "struct name_orient");
          no->name = CHECKED_STRDUP(w->parent_object->sym->name,
                                    "w->parent_object->sym->name");
          no->orient = 1;
          no->next = NULL;
          no_head = no;
          tail = no_head;
        } else {
          found = 0;
          for (nol = no_head; nol != NULL; nol = nol->next) {
            if (strcmp(nol->name, w->parent_object->sym->name) == 0) {
              nol->orient++;
              found = 1;
              break;
            }
          }
          if (!found) {
            /* add to the end of list */
            no =
                CHECKED_MALLOC_STRUCT(struct name_orient, "struct name_orient");
            no->name = CHECKED_STRDUP(w->parent_object->sym->name,
                                      "w->parent_object->sym->name");
            no->orient = 1;
            no->next = tail->next;
            tail->next = no;
            tail = tail->next;
          }
        }

      } else if ((smash->what & COLLIDE_SUBVOL) != 0) {

        struct subvolume *nsv;

        virt_mol.pos.x = smash->loc.x;
        virt_mol.pos.y = smash->loc.y;
        virt_mol.pos.z = smash->loc.z;

        nsv = traverse_subvol(
            sv, &(virt_mol.pos), smash->what - COLLIDE_SV_NX - COLLIDE_SUBVOL,
            state->nx_parts, state->ny_parts, state->nz_parts);
        // Hit the edge of the world
        if (nsv == NULL) {
          if (shead != NULL)
            mem_put_list(sv->local_storage->coll, shead);

          for (nol = no_head; nol != NULL; nol = nol->next) {
            if (nol->orient % 2 != 0) {
              if (farthest_flag) {
                farthest_name = nol->name;
              }
              else {
                return_name = CHECKED_STRDUP(nol->name, "nol->name");
                break;
              }
            }
          }
          
          if (farthest_name != NULL) {
            return_name = CHECKED_STRDUP(farthest_name, "nol->name");
          }

          while (no_head != NULL) {

            nnext = no_head->next;
            free(no_head->name);
            free(no_head);
            no_head = nnext;
          }
          
          return return_name;
        }

        if (shead != NULL)
          mem_put_list(sv->local_storage->coll, shead);
        calculate_random_vector = 0;
        virt_mol.subvol = nsv;

        // Jump to beginning of function
        goto pretend_to_call_find_enclosing_mesh;
      }
    }

  } while (smash != NULL);

  // I'm not sure when we would ever get to this point since our random vector
  // is big enough to always hit the edge of the world and return that way.
  // This might be dead code, but I'll leave it in for now.
  if (shead != NULL)
    mem_put_list(sv->local_storage->coll, shead);

  for (nol = no_head; nol != NULL; nol = nol->next) {

    if (nol->orient % 2 != 0) {
      if (farthest_flag) {
        farthest_name = nol->name;
      }
      else {
        return_name = CHECKED_STRDUP(nol->name, "nol->name");
        break;
      }
    }
  }
  
  if (farthest_name != NULL) {
    return_name = CHECKED_STRDUP(farthest_name, "nol->name");
  }

  while (no_head != NULL) {
    nnext = no_head->next;
    free(no_head->name);
    free(no_head);
    no_head = nnext;
  }

  if (return_name != NULL)
    return return_name;
  else
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
void place_mol_relative_to_mesh(struct volume *state, struct vector3 *loc,
                                struct subvolume *sv, char *mesh_name,
                                struct vector3 *new_pos) {
  struct wall *best_w = NULL;
  struct wall_list *wl;
  double d2, best_d2;
  struct vector2 s_loc;
  struct vector3 best_xyz;
  char *mesh_name_try = NULL; /* farthest enclosing mesh name */
  struct volume_molecule virt_mol;

  best_w = NULL;
  best_d2 = GIGANTIC + 1;
  best_xyz.x = 0;
  best_xyz.y = 0;
  best_xyz.z = 0;

  if (mesh_name == NULL) {
    /* we have to move molecule that now appeared to be inside the mesh outside
     * of all enclosing meshes */
    virt_mol.prev_v = NULL;
    virt_mol.next_v = NULL;
    virt_mol.next = NULL;
    virt_mol.pos.x = loc->x;
    virt_mol.pos.y = loc->y;
    virt_mol.pos.z = loc->z;
    virt_mol.subvol = find_subvolume(state, loc, NULL);

    int farthest_flag = 1;
    mesh_name_try = find_enclosing_mesh_name(state, &virt_mol, farthest_flag);
    if (mesh_name_try == NULL) {
      mcell_internal_error("Cannot find the farthest enclosing mesh.");
    }
  }

  for (wl = sv->wall_head; wl != NULL; wl = wl->next) {
    if (mesh_name != NULL) {
      if (strcmp(wl->this_wall->parent_object->sym->name, mesh_name) != 0) {
        continue;
      }
    } else {
      if (strcmp(wl->this_wall->parent_object->sym->name, mesh_name_try) != 0) {
        continue;
      }
    }

    d2 = closest_interior_point(loc, wl->this_wall, &s_loc, GIGANTIC);
    if (d2 < best_d2) {
      best_d2 = d2;
      best_w = wl->this_wall;
    }
  }

  /* look into neighbor subvolumes */
  const int sv_index = sv - state->subvol;
  int sv_remain = sv_index;

  /* Turn linear sv_index into part_x, part_y, part_z triple. */
  const int part_x =
      sv_remain / ((state->ny_parts - 1) * (state->nz_parts - 1));
  sv_remain -= part_x * ((state->ny_parts - 1) * (state->nz_parts - 1));
  const int part_y = sv_remain / (state->nz_parts - 1);
  sv_remain -= part_y * (state->nz_parts - 1);
  const int part_z = sv_remain;

  /* Find min x partition. */
  int x_min;
  for (x_min = part_x; x_min > 0; x_min--) {
    d2 = loc->x - state->x_partitions[x_min];
    d2 *= d2;
    if (d2 >= best_d2)
      break;
  }

  /* Find max x partition. */
  int x_max;
  for (x_max = part_x; x_max < state->nx_parts - 1; x_max++) {
    d2 = loc->x - state->x_partitions[x_max + 1];
    d2 *= d2;
    if (d2 >= best_d2)
      break;
  }

  /* Find min y partition. */
  int y_min;
  for (y_min = part_y; y_min > 0; y_min--) {
    d2 = loc->y - state->y_partitions[y_min];
    d2 *= d2;
    if (d2 >= best_d2)
      break;
  }

  /* Find max y partition. */
  int y_max;
  for (y_max = part_y; y_max < state->ny_parts - 1; y_max++) {
    d2 = loc->y - state->y_partitions[y_max + 1];
    d2 *= d2;
    if (d2 >= best_d2)
      break;
  }

  /* Find min z partition. */
  int z_min;
  for (z_min = part_z; z_min > 0; z_min--) {
    d2 = loc->z - state->z_partitions[z_min];
    d2 *= d2;
    if (d2 >= best_d2)
      break;
  }

  /* Find max z partition. */
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

          for (wl = state->subvol[this_sv].wall_head; wl != NULL;
               wl = wl->next) {
            if (mesh_name != NULL) {
              if (strcmp(wl->this_wall->parent_object->sym->name, mesh_name) !=
                  0) {
                continue;
              }
            } else {
              if (strcmp(wl->this_wall->parent_object->sym->name,
                         mesh_name_try) !=
                  0) {
                continue;
              }
            }

            d2 = closest_interior_point(loc, wl->this_wall, &s_loc, GIGANTIC);
            if (d2 < best_d2) {
              best_d2 = d2;
              best_w = wl->this_wall;
            }
          }
        }
      }
    }
  }

  if (mesh_name_try != NULL) {
    free(mesh_name_try);
  }

  if (best_w != NULL) {
    find_wall_center(best_w, &best_xyz);
  } else
    mcell_internal_error(
        "Error in function 'place_mol_relative_to_mesh()'.");

  /* We will return the point just behind the closest (or farthest)
     enclosing mesh */

  /* the parametric equation of ray is L(t) = A + t(B - A),
     where A - start vector, B - some point on the ray,
     and parameter t >= 0 */

  /* If we need to place molecule inside the closest enclosing mesh
     we select a point that lies on the inward directed wall normal that
     starts at the point "best_xyz" and is located at the distance of 2*EPC_C.
     If we need to place molecule outside of the farthest enclosing mesh
     we select a point that lies on the outward directed wall normal that
     starts at the point "best_xyz" and is located at the distance of 2*EPC_C.
   */

  if (mesh_name != NULL) {
    new_pos->x = best_xyz.x - 2 * MESH_DISTINCTIVE * (best_w->normal.x);
    new_pos->y = best_xyz.y - 2 * MESH_DISTINCTIVE * (best_w->normal.y);
    new_pos->z = best_xyz.z - 2 * MESH_DISTINCTIVE * (best_w->normal.z);

  } else {
    new_pos->x = best_xyz.x + 2 * MESH_DISTINCTIVE * (best_w->normal.x);
    new_pos->y = best_xyz.y + 2 * MESH_DISTINCTIVE * (best_w->normal.y);
    new_pos->z = best_xyz.z + 2 * MESH_DISTINCTIVE * (best_w->normal.z);
  }
}

/***************************************************************************
destroy_everything:
  In: state
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
  state->root_object->first_child = NULL; 
  state->root_object->last_child = NULL; 

  free(state->all_vertices);
  state->n_walls = 0;
  state->n_verts = 0;

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

  // Destroy partitions and boundaries
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

  free(state->waypoints);

  return 0;
}

/***************************************************************************
destroy_objects:
  In: obj_ptr - object to be destroyed
      free_poly_flag - see explanation in destroy_poly_object
  Out: Zero on success. One otherwise. Recursively destroys objects.
  Note: Currently, this ultimately only destroys polygon objects. I don't know
        if there's a need to trash release objects that use release patterns.
***************************************************************************/
int destroy_objects(struct object *obj_ptr, int free_poly_flag) {
  switch (obj_ptr->object_type) {
  case META_OBJ:
    for (struct object *child_obj_ptr = obj_ptr->first_child;
         child_obj_ptr != NULL; child_obj_ptr = child_obj_ptr->next) {
      destroy_objects(child_obj_ptr, free_poly_flag);
    }
    break;
  case BOX_OBJ:
  case POLY_OBJ:
    destroy_poly_object(obj_ptr, free_poly_flag);
    break;

  // do nothing
  case REL_SITE_OBJ:
  case VOXEL_OBJ:
  default:
    break;
  }

  return 0;
}

/***************************************************************************
destroy_poly_object:
  In: obj_ptr - object
      freep_poly_flag - Destroy polygon_object if set. There's a link to this
      in the object definition (child of root_object) AND the instantiated
      object (child of root_insance), so we only want to free it once.
  Out: Zero on success. One otherwise. Polygon object is destroyed
***************************************************************************/
int destroy_poly_object(struct object *obj_ptr, int free_poly_flag) {
  if (free_poly_flag) {
    for (int wall_num = 0; wall_num < obj_ptr->n_walls; wall_num++) {
      struct wall *w = obj_ptr->wall_p[wall_num];
      if (w->grid) {
        free(w->grid->mol);
      } 
    }
    struct polygon_object *poly_obj_ptr = obj_ptr->contents;
    free(poly_obj_ptr->side_removed);
    free(poly_obj_ptr->element);
    free(obj_ptr->contents);
    obj_ptr->contents = NULL;
  }
  free(obj_ptr->walls);
  free(obj_ptr->wall_p);
  free(obj_ptr->vertices);

  obj_ptr->wall_p = NULL;
  obj_ptr->n_walls = 0;
  obj_ptr->n_walls_actual = 0;
  obj_ptr->n_verts = 0;
  struct region_list *regs, *next_regs;
  for (regs = obj_ptr->regions; regs != NULL;) {
    free(regs->reg->membership);
    regs->reg->membership = NULL;
    free(regs->reg->bbox);
    regs->reg->bbox = NULL;
    if (regs->reg->boundaries) {
      pointer_hash_destroy(regs->reg->boundaries);
    }
    free(regs->reg->boundaries); 
    regs->reg->boundaries = NULL;
    next_regs = regs->next;
    free(regs);
    regs = next_regs;
  }
  obj_ptr->regions = NULL;
  obj_ptr->num_regions = 0;

  return 0;
}

/***************************************************************************
reset_current_counts:
  In: state
  Out: Zero on success. Species populations are set to zero. Counts on/in
       regions are also set to zero.
***************************************************************************/
int reset_current_counts(struct volume *state) {
  // Set global populations of species back to zero, since they will get set to
  // the proper values when we insert the molecules into the world
  struct sym_table_head *mol_sym_table = state->mol_sym_table;
  for (int n_mol_bin = 0; n_mol_bin < mol_sym_table->n_bins; n_mol_bin++) {
    for (struct sym_table *sym_ptr = mol_sym_table->entries[n_mol_bin];
         sym_ptr != NULL; sym_ptr = sym_ptr->next) {
      struct species *mol = (struct species *)sym_ptr->value;
      mol->population = 0;
    }
  }

  // Set counts on/in regions back to zero for the same reasons listed above.
  for (int i = 0; i <= state->count_hashmask; i++)
    if (state->count_hash[i] != NULL) {
      struct counter *c;
      for (c = state->count_hash[i]; c != NULL; c = c->next) {
        if ((c->counter_type & MOL_COUNTER) != 0) {
          c->data.move.n_enclosed = 0; 
          c->data.move.n_at = 0; 
        }
      }
    }
  return 0;
}

/***************************************************************************
enable_counting_for_all_objects:
  In: obj_ptr
  Out: Zero on success. Enable counting for every polygon object.
***************************************************************************/
int enable_counting_for_all_objects(struct object *obj_ptr) {
  switch (obj_ptr->object_type) {
  case META_OBJ:
    for (struct object *child_obj_ptr = obj_ptr->first_child;
         child_obj_ptr != NULL; child_obj_ptr = child_obj_ptr->next) {
      enable_counting_for_all_objects(child_obj_ptr);
    }
    break;
  case BOX_OBJ:
  case POLY_OBJ:
    enable_counting_for_object(obj_ptr);
    break;
  // do nothing
  case REL_SITE_OBJ:
  case VOXEL_OBJ:
  default:
    break;
  }
  return 0;
}

/***************************************************************************
enable_counting_for_object:
  In: obj_ptr
  Out: Zero on success. Enable counting for every region on an object.
***************************************************************************/
int enable_counting_for_object(struct object *obj_ptr) {
  struct region_list *regs;
  for (regs = obj_ptr->regions; regs != NULL; regs=regs->next) {
    regs->reg->flags |= COUNT_CONTENTS;
    //regs->reg->flags |= COUNT_ENCLOSED;
  }
  return 0;
}

/***************************************************************************
init_mol_obj_transp:
  In: state: simulation state
  Out: Zero on success. Create a data structure so we can quickly check if a
  molecule species can move in or out of any given surface region
***************************************************************************/
int init_mol_obj_transp(struct volume *state) {
  struct pointer_hash *mol_obj_transp;
  if ((mol_obj_transp = CHECKED_MALLOC_STRUCT(struct pointer_hash,
                                              "pointer_hash")) == NULL) {
    mcell_internal_error("Out of memory while creating molecule-object "
                         "transparency pointer hash");
  }
  if (pointer_hash_init(mol_obj_transp, state->n_species)) {
    mcell_error(
      "Failed to initialize data structure for molecule-object transparency.");
    return 1;
  }
  // Initialize pointer hash with molecule/species names as keys and values are
  // a linked list of pointers of object_transparency. The object_transparency
  // struct contains the object name and whether the species can go in_to_out
  // and/or out_to_in. These are initially set to 0 (not transparent), and can
  // be set to 1 (transparent) in find_obj_region_transp.
  state->mol_obj_transp = mol_obj_transp;
  for (int i = 0; i < state->n_species; i++) {
    struct species *spec = state->species_list[i];
    char *species_name = spec->sym->name;
    // Only check volume molecules (for now).
    if (spec->flags & ON_GRID)
      continue;
    else if (spec->flags & IS_SURFACE)
      continue;
    if (strcmp(species_name, "ALL_MOLECULES") == 0)
      continue;
    else if (strcmp(species_name, "ALL_VOLUME_MOLECULES") == 0)
      continue;
    else if (strcmp(species_name, "ALL_SURFACE_MOLECULES") == 0)
      continue;
    unsigned int keyhash = (unsigned int)(intptr_t)(species_name);
    void *key = (void *)(species_name);
    struct object_transparency *obj_transp_head = NULL;
    struct object_transparency *obj_transp_tail = NULL;
    find_all_obj_region_transp(state->root_instance, &obj_transp_head,
                               &obj_transp_tail, species_name);
    if (pointer_hash_add(
        state->mol_obj_transp, key, keyhash, (void *)obj_transp_head)) {
      mcell_allocfailed("Failed to store molecule-object transparency in"
                        " pointer_hash table.");
    }
  }
  return 0;
}

/***************************************************************************
find_obj_region_transp:
  In:  obj_ptr: The object we are currently checking for transparency
       obj_transp_head: Head of the object transparency list
       obj_transp_tail: Tail of the object transparency list
       species_name: The name of the molecule/species we are checking
  Out: Zero on success. Check every region on obj_ptr to see if any of them are
       transparent to species_name.
***************************************************************************/
int find_obj_region_transp(struct object *obj_ptr,
                           struct object_transparency **obj_transp_head,
                           struct object_transparency **obj_transp_tail,
                           char *species_name) {
  struct object_transparency *obj_transp;
  obj_transp =
      CHECKED_MALLOC_STRUCT(struct object_transparency, "object transparency");
  obj_transp->next = NULL;
  obj_transp->obj_name = obj_ptr->sym->name;
  obj_transp->in_to_out = 0;
  obj_transp->out_to_in = 0;
  if (*obj_transp_tail == NULL) {
    *obj_transp_head = obj_transp; 
    *obj_transp_tail = obj_transp; 
    (*obj_transp_head)->next = NULL;
    (*obj_transp_tail)->next = NULL;
  }
  else {
    (*obj_transp_tail)->next = obj_transp;
    *obj_transp_tail = obj_transp;
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
             obj_transp->out_to_in = 1;
             break;
          }
          // Transparent from inside to outside
          else if ((no->orient == -1 && volume > 0) ||
              (no->orient == 1 && volume < 0)) {
             obj_transp->in_to_out = 1;
             break;
          }
          // Transparent from either direction (e.g. TRANSPARENT = A;)
          else if (no->orient == 0) {
             obj_transp->in_to_out = 1;
             obj_transp->out_to_in = 1;
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
  In: obj_ptr:
      obj_transp_head: Head of the object transparency list
      obj_transp_tail: Tail of the object transparency list
      species_name: The name of the molecule/species we are checking
  Out: Zero on success. Check every polygon object to see if it is transparent
       to species_name.
***************************************************************************/
int find_all_obj_region_transp(struct object *obj_ptr,
                               struct object_transparency **obj_transp_head,
                               struct object_transparency **obj_transp_tail,
                               char *species_name) {
  switch (obj_ptr->object_type) {
  case META_OBJ:
    for (struct object *child_obj_ptr = obj_ptr->first_child;
         child_obj_ptr != NULL; child_obj_ptr = child_obj_ptr->next) {
      if (find_all_obj_region_transp(
          child_obj_ptr, obj_transp_head, obj_transp_tail, species_name))
        return 1;
    }
    break;

  case REL_SITE_OBJ:
    break;

  case BOX_OBJ:
  case POLY_OBJ:
    if (find_obj_region_transp(obj_ptr, obj_transp_head, obj_transp_tail,
                               species_name))
      return 1;
    break;

  default:
    UNHANDLED_CASE(obj_ptr->object_type);
  }

  return 0;
}
