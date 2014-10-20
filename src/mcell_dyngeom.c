#include <stdlib.h>
#include <string.h>

#include "mcell_misc.h"
#include "mcell_dyngeom.h"
#include "chkpt.h"
#include "vol_util.h"
#include "grid_util.h"
#include "wall_util.h"
#include "init.h"
#include "logging.h"

#define NO_MESH "\0"

/************************************************************************
 mcell_add_dynamic_geometry:
 In:  dynamic_geometry_filepath:
      curr_mdl_filepath:
      time_unit:
      dynamic_geometry_mem:
      dynamic_geometry_head:
 Out: 0 on success, 1 on failure. dynamic geometry events are added to
      dynamic_geometry_head from which they will eventually be added to a
      scheduler.
 ***********************************************************************/
int mcell_add_dynamic_geometry(
    char const *dynamic_geometry_filepath, char const *curr_mdl_filepath,
    double time_unit, struct mem_helper *dynamic_geometry_mem, 
    struct dynamic_geometry **dynamic_geometry_head) {

  char *dyn_geom_filename = mcell_find_include_file(
      dynamic_geometry_filepath, curr_mdl_filepath);
  FILE *f = fopen(dyn_geom_filename, "r");

  if (!f) {
    free(dyn_geom_filename);
    return 1;
  }
  else {
    const char *RATE_SEPARATORS = "\f\n\r\t\v ,;";
    const char *FIRST_DIGIT = "+-0123456789";
    struct dynamic_geometry *dyn_geom_tail = NULL;
    char buf[2048];
    char *char_ptr;
    int linecount = 0;
    int i;

    while (fgets(buf, 2048, f)) {
      linecount++;
      // ignore leading whitespace
      for (i = 0; i < 2048; i++) {
        if (!strchr(RATE_SEPARATORS, buf[i]))
          break;
      }

      if (i < 2048 && strchr(FIRST_DIGIT, buf[i])) {
        // Grab time
        double time = strtod((buf + i), &char_ptr);
        if (char_ptr == (buf + i))
          continue; /* Conversion error. */

        // Skip over whitespace between time and filename
        for (i = char_ptr - buf; i < 2048; i++) {
          if (!strchr(RATE_SEPARATORS, buf[i]))
            break;
        }

        // Grab mdl filename. This could probably be cleaned up
        char *line_ending = strchr(buf+i,'\n');
        int line_ending_idx = line_ending-buf;
        int file_name_length = line_ending_idx-i+1;
        char file_name[file_name_length];
        strncpy(file_name, buf+i, file_name_length);
        file_name[file_name_length-1] = '\0';

        struct dynamic_geometry *dyn_geom;
        dyn_geom = CHECKED_MEM_GET(dynamic_geometry_mem,
                                   "time-varying dynamic geometry");
        if (dyn_geom == NULL)
          return 1;
         
        // I think we should probably wait until sim init to do these kinds of
        // conversions, but it's here for consistency
        dyn_geom->event_time = time / time_unit;
        dyn_geom->mdl_file_path = strdup(file_name);
        dyn_geom->next = NULL;

        // Append each entry to end of dynamic_geometry_head list
        if (*dynamic_geometry_head == NULL) {
          *dynamic_geometry_head = dyn_geom;
          dyn_geom_tail = dyn_geom;
        }
        else {
          dyn_geom_tail->next = dyn_geom;
          dyn_geom_tail = dyn_geom;
        }
      }
    }

    fclose(f);
  }

  free(dyn_geom_filename);
  return 0;

}

/***************************************************************************
 save_all_molecules: Save all the molecules currently in the scheduler.

 In:  state: MCell state
      storage_head:
 Out: An array of all the molecules to be saved
 Note: This is meant to be used with dynamic geometries. Currently only saves
       volume molecules.
***************************************************************************/
struct molecule_info ** save_all_molecules(
    struct volume *state, struct storage_list *storage_head) {
  
  // Find total number of molecules in the scheduler.
  unsigned long long total_items = count_items_in_scheduler(storage_head);
  int ctr = 0;
  struct molecule_info **all_molecules;
  all_molecules = CHECKED_MALLOC_ARRAY(
      struct molecule_info *, total_items, "all molecules");

  // Iterate over all molecules in the scheduler.
  for (struct storage_list *sl_ptr = storage_head; sl_ptr != NULL;
       sl_ptr = sl_ptr->next) {
    for (struct schedule_helper *sh_ptr = sl_ptr->store->timer; sh_ptr != NULL;
         sh_ptr = sh_ptr->next_scale) {
      for (int i = -1; i < sh_ptr->buf_len; i++) {
        for (struct abstract_element *ae_ptr = (i < 0) ? sh_ptr->current
                                                    : sh_ptr->circ_buf_head[i];
             ae_ptr != NULL; ae_ptr = ae_ptr->next) {
          struct abstract_molecule *am_ptr = (
              struct abstract_molecule *)ae_ptr;
          if (am_ptr->properties == NULL)
            continue;

          struct molecule_info *mol_info = CHECKED_MALLOC_STRUCT(
              struct molecule_info, "molecule info");
          all_molecules[ctr] = mol_info;
          mol_info->molecule = CHECKED_MALLOC_STRUCT(
              struct abstract_molecule, "abstract molecule");
  
          const int MAX_NUM_REGIONS = 100;
          struct string_buffer *reg_names = CHECKED_MALLOC_STRUCT(
              struct string_buffer, "string buffer");
          if (initialize_string_buffer(reg_names, MAX_NUM_REGIONS)) {
            return NULL; 
          }

          char *mesh_name;
          if ((am_ptr->properties->flags & NOT_FREE) == 0) {
            save_volume_molecule(state, mol_info, am_ptr, &mesh_name);
          }
          else if ((am_ptr->properties->flags & ON_GRID) != 0) { 
            if (save_surface_molecule(mol_info, am_ptr, &reg_names, &mesh_name))
              return NULL;
          }
          else {
            continue; 
          }

          save_common_molecule_properties(
              mol_info, am_ptr, reg_names, mesh_name);
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
      mesh_name:
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
  mol_info->molecule->mesh_name = strdup(mesh_name);
  // Only free temporary object names we just allocated above.
  // Don't want to accidentally free symbol names of objects.
  if ((mesh_name != NO_MESH) &&
      ((am_ptr->properties->flags & NOT_FREE) == 0)) {
    free(mesh_name); 
  }
  mol_info->reg_names = reg_names;
}

/***************************************************************************
 save_volume_molecule: 

 In:  state:
      mol_info:
      am_ptr: abstract molecule pointer
      mesh_name:
 Out: Nothing. Molecule info and mesh name are updated
***************************************************************************/
void save_volume_molecule(struct volume *state, struct molecule_info *mol_info,
                          struct abstract_molecule *am_ptr, char **mesh_name) {
  struct volume_molecule *vm_ptr = (struct volume_molecule *)am_ptr;

  *mesh_name = find_closest_enclosing_mesh_name(state, vm_ptr); 
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
      mesh_name:
 Out: Zero on success. One otherwise.
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
  int num_regions;
  struct name_list *reg_name_list_head, *reg_name_list;
  reg_name_list_head = find_regions_names_by_wall(
      sm_ptr->grid->surface, &num_regions);
  int k;
  for (reg_name_list = reg_name_list_head, k = 0; reg_name_list != NULL;
       reg_name_list = reg_name_list->next, k++) {
    char *str = strdup(reg_name_list->name);
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
 Note: This is meant to be used with dynamic geometries. Currently only places
       volume molecules.
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

      if(strlen(mesh_name) == 0) {
          vm_guess = insert_volume_molecule_encl_mesh(
              state, vm_ptr, vm_guess, NULL);  
      }
      else {
          vm_guess = insert_volume_molecule_encl_mesh(
              state, vm_ptr, vm_guess, mesh_name);  
      }

      if (vm_guess == NULL) {
        mcell_error("Cannot insert copy of molecule of species '%s' into "
                    "world.\nThis may be caused by a shortage of memory.",
                    vm_ptr->properties->sym->name);
      }
    }
    // Insert surface molecule into world. 
    else if ((am_ptr->properties->flags & ON_GRID) != 0) { 
      insert_surface_molecule(
          state, am_ptr->properties, &mol_info->pos, mol_info->orient,
          CHKPT_GRID_TOLERANCE, am_ptr->t, NULL, mesh_name,
          mol_info->reg_names);

    }
  }

  // Do some cleanup.
  for (int i=0; i<state->num_all_molecules; i++) {
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
