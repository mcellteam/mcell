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
#include "count_util.h"
#include "diffuse.h"

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

/*************************************************************************
insert_volume_molecule_encl_mesh:
  In: pointer to a volume_molecule that we're going to place in local storage
      pointer to a volume_molecule that may be nearby
      closest enclosing mesh name
  Out: pointer to the new volume_molecule (copies data from volume molecule
       passed in), or NULL if out of memory.  Molecule is placed in scheduler
       also.
*************************************************************************/

struct volume_molecule* insert_volume_molecule_encl_mesh(
    struct volume *state, struct volume_molecule *m,
    struct volume_molecule *vm_guess, char *mesh_name)
{
  struct volume_molecule *new_m;
  struct subvolume *sv, *new_sv;;
  char *mesh_name_try;
  int move_molecule = 0;
  struct vector3 new_pos;

  if (vm_guess == NULL) sv = find_subvolume(state, &(m->pos), NULL);
  else if (inside_subvolume(
      &(m->pos), vm_guess->subvol, state->x_fineparts, state->y_fineparts,
      state->z_fineparts)) {
    sv = vm_guess->subvol;
  }
  else sv = find_subvolume(state, &(m->pos), vm_guess->subvol);
  
  new_m = CHECKED_MEM_GET(sv->local_storage->mol, "volume molecule");
  memcpy(new_m, m, sizeof(struct volume_molecule));
  new_m->mesh_name = NULL;
  new_m->prev_v = NULL;
  new_m->next_v = NULL;
  new_m->next = NULL;
  new_m->subvol = sv;

  mesh_name_try = find_closest_enclosing_mesh_name(state, new_m);

  /* mol was inside mesh, now it is outside mesh */
  if ((mesh_name_try == NULL) && (mesh_name != NULL)) {
    move_molecule = 1;    
  }
  /* mol was outside mesh, now it is inside mesh */
  if ((mesh_name_try != NULL) && (mesh_name == NULL)) {
    move_molecule = 1;
    free(mesh_name_try);
  }
  if ((mesh_name_try != NULL) && (mesh_name != NULL))
  {
    /* mol was inside one mesh, now it is inside another mesh */
    if (strcmp(mesh_name_try, mesh_name) != 0) move_molecule = 1;
    free(mesh_name_try);
  }

  if (move_molecule)
  {
     /* move molecule to another location so that closest
        enclosing mesh name is "mesh_name" */

     place_mol_relative_to_mesh(state, &(m->pos), sv, mesh_name, &new_pos);
     new_m->pos = new_pos;
     new_sv = find_subvolume(state, &(new_m->pos), NULL);
     new_m->subvol = new_sv;
       
  }  

  new_m->birthplace = new_m->subvol->local_storage->mol;
  ht_add_molecule_to_list(&(new_m->subvol->mol_by_species), new_m);
  new_m->subvol->mol_count++;
  new_m->properties->population++;

  if ((new_m->properties->flags&COUNT_SOME_MASK) != 0) new_m->flags |= COUNT_ME;
  if (new_m->properties->flags & (COUNT_CONTENTS|COUNT_ENCLOSED))
  {
    count_region_from_scratch(state, (struct abstract_molecule*)new_m, NULL, 1,
                              &(new_m->pos), NULL, new_m->t);
  }
  
  if ( schedule_add(new_m->subvol->local_storage->timer,new_m) )
    mcell_allocfailed("Failed to add volume molecule to scheduler.");

  return new_m;
}

/*************************************************************************
find_closest_enclosing_mesh_name:
  In: volume molecule 
  Out: Name of the closest enclosing mesh if such exists,
       NULL otherwise.
************************************************************************/
char * find_closest_enclosing_mesh_name(
    struct volume *state, struct volume_molecule *m)
{
  struct collision *smash;     /* Thing we've hit that's under consideration */
  struct collision *shead;     /* Head of the linked list of collisions */

  struct subvolume *sv;
  struct wall *w;
  struct vector3 rand_vector, world_diag;
  struct volume_molecule virt_mol;  /* volume_molecule template */
  double world_diag_length;  /* length of the world bounding box diagonal */

  /* we will reuse this struct in the different context of 
     registering the meshes names through "no->name" and 
     the number of hits with the mesh through "no->orient" */
  struct name_orient *no, *nol, *nnext, *no_head = NULL, *tail = NULL;
  char *return_name = NULL;
     
  memcpy(&virt_mol, m, sizeof(struct volume_molecule));
  virt_mol.prev_v = NULL;
  virt_mol.next_v = NULL;
  virt_mol.next = NULL;

  int calculate_random_vector = 1;   /* flag */
  int found; /* flag */

pretend_to_call_find_enclosing_mesh:   /* Label to allow fake recursion */

  sv = virt_mol.subvol;

  shead = NULL;

  /* pick up a random vector */
  if(calculate_random_vector)
  {
    srand((unsigned int)time(NULL));
    rand_vector.x = (double)rand()/(double)RAND_MAX;
    rand_vector.y = (double)rand()/(double)RAND_MAX;
    rand_vector.z = (double)rand()/(double)RAND_MAX;

    /* find the diagonal of the world */
    vectorize(&state->bb_urb, &state->bb_llf, &world_diag);
    world_diag_length = vect_length(&world_diag);

    /* scale random vector by the size of the world */
    rand_vector.x *= world_diag_length;
    rand_vector.y *= world_diag_length;
    rand_vector.z *= world_diag_length;

  }

  do
  {
     shead = ray_trace(state, &virt_mol,NULL,sv,&rand_vector,NULL);
     if (shead==NULL) mcell_internal_error("ray_trace() returned NULL.");

     if (shead->next!=NULL)
     {
        shead = (struct collision*)ae_list_sort((
            struct abstract_element*)shead);
     }

     for (smash = shead; smash != NULL; smash = smash->next)
     {
         if ( (smash->what & COLLIDE_WALL) != 0 )
         {
            w = (struct wall*) smash->target;
    
            /* discard open-type meshes, like planes, etc. */
            if(w->parent_object->is_closed <= 0) continue; 
         
            /* discard the cases when the random vector just grazes 
               the mesh at the encounter point */
            double d_prod = dot_prod(&rand_vector, &(w->normal));
            if(!distinguishable(d_prod, 0, EPS_C)) continue;
  
            if(no_head == NULL)
            {
              no = CHECKED_MALLOC_STRUCT(
                  struct name_orient, "struct name_orient");
              no->name = CHECKED_STRDUP(
                  w->parent_object->sym->name, "w->parent_object->sym->name");
              no->orient = 1;
              no->next = NULL; 
              no_head = no;
              tail = no_head;
            }else{
              found = 0;
              for(nol = no_head; nol != NULL; nol = nol->next)
              {
                if(strcmp(nol->name, w->parent_object->sym->name) == 0) 
                {
                  nol->orient++;
                  found = 1;
                  break;
                }
              }
              if(!found)
              {
                /* add to the end of list */
                no = CHECKED_MALLOC_STRUCT(
                    struct name_orient, "struct name_orient");
                no->name = CHECKED_STRDUP(
                    w->parent_object->sym->name, "w->parent_object->sym->name");
                no->orient = 1;
                no->next = tail->next;
                tail->next = no;
                tail = tail->next;
              }

            }

         }else if ((smash->what & COLLIDE_SUBVOL) != 0){

              struct subvolume *nsv;

              virt_mol.pos.x = smash->loc.x;
              virt_mol.pos.y = smash->loc.y;
              virt_mol.pos.z = smash->loc.z;

              nsv = traverse_subvol(
                  sv, &(virt_mol.pos),
                  smash->what - COLLIDE_SV_NX - COLLIDE_SUBVOL,
                  state->nx_parts, state->ny_parts, state->nz_parts); 
              if (nsv==NULL)
              {

                if (shead != NULL) mem_put_list(sv->local_storage->coll,shead);

                for(nol = no_head; nol != NULL; nol = nol->next)
                {
                  if(nol->orient%2 != 0) 
                  {
                     return_name = CHECKED_STRDUP(nol->name, "nol->name");
                     break;
                  }
                }
  
                while(no_head != NULL)
                {

                  nnext = no_head->next;
                  free(no_head->name);
                  free(no_head);
                  no_head = nnext;
                }

                return return_name;

              }
        
              if (shead != NULL) mem_put_list(sv->local_storage->coll,shead);
              calculate_random_vector = 0;
              virt_mol.subvol = nsv;

              // Jump to beginning of function
              goto pretend_to_call_find_enclosing_mesh;          
         }
     }
  
  }while(smash != NULL);

  if (shead != NULL) mem_put_list(sv->local_storage->coll,shead);

  for(nol = no_head; nol != NULL; nol = nol->next)
  {

    if(nol->orient%2 != 0) 
    {
      return_name = CHECKED_STRDUP(nol->name, "nol->name");
      break;
    }
  }

  while(no_head != NULL)
  {
    nnext = no_head->next;
    free(no_head->name);
    free(no_head);
    no_head = nnext;
  }

  return return_name;

}

/*************************************************************************
find_farthest_enclosing_mesh_name:
  In: volume molecule 
  Out: Name of the farthest enclosing mesh if such exists,
       NULL otherwise.
************************************************************************/
char * find_farthest_enclosing_mesh_name(
    struct volume *state, struct volume_molecule *m)
{
  struct collision *smash;     /* Thing we've hit that's under consideration */
  struct collision *shead;     /* Head of the linked list of collisions */

  struct subvolume *sv;
  struct wall *w;
  struct vector3 rand_vector, world_diag;
  struct volume_molecule virt_mol;  /* volume_molecule template */
  double world_diag_length;  /* length of the world bounding box diagonal */
 
  memcpy(&virt_mol, m, sizeof(struct volume_molecule));
  virt_mol.prev_v = NULL;
  virt_mol.next_v = NULL;
  virt_mol.next = NULL;
  
  /* we will reuse this struct in the different context of 
     registering the meshes names through "no->name" and 
     the number of hits with the mesh through "no->orient" */
  struct name_orient *no, *nol, *nnext, *no_head = NULL, *tail = NULL;
  char *return_name = NULL, *best_name = NULL;   

  int calculate_random_vector = 1;   /* flag */
  int found; /* flag */

pretend_to_call_find_enclosing_mesh:   /* Label to allow fake recursion */

  sv = virt_mol.subvol;
  shead = NULL;

  /* pick up a random vector */
  if(calculate_random_vector)
  {
    srand((unsigned int)time(NULL));
    rand_vector.x = (double)rand()/(double)RAND_MAX;
    rand_vector.y = (double)rand()/(double)RAND_MAX;
    rand_vector.z = (double)rand()/(double)RAND_MAX;

    /* find the diagonal of the world */
    vectorize(&state->bb_urb, &state->bb_llf, &world_diag);
    world_diag_length = vect_length(&world_diag);

    /* scale random vector by the size of the world */
    rand_vector.x *= world_diag_length;
    rand_vector.y *= world_diag_length;
    rand_vector.z *= world_diag_length;

  }

  do
  {

     shead = ray_trace(state, &virt_mol,NULL,sv,&rand_vector,NULL);
     if (shead==NULL) mcell_internal_error("ray_trace() returned NULL.");

     if (shead->next!=NULL)
     {
        shead = (struct collision*)ae_list_sort((
            struct abstract_element*)shead);
     }

     for (smash = shead; smash != NULL; smash = smash->next)
     {
         if ( (smash->what & COLLIDE_WALL) != 0 )
         {
            w = (struct wall*) smash->target;
            
            /* discard open-type meshes, like planes, etc. */
            if(w->parent_object->is_closed <= 0) continue; 
         
            /* discard the cases when the random vector just grazes 
               the mesh at the encounter point */
            double d_prod = dot_prod(&rand_vector, &(w->normal));
            if(!distinguishable(d_prod, 0, EPS_C)) continue;
  
            if(no_head == NULL)
            {
              no = CHECKED_MALLOC_STRUCT(
                  struct name_orient, "struct name_orient");
              no->name = CHECKED_STRDUP(
                  w->parent_object->sym->name, "w->parent_object->sym->name");
              no->orient = 1;
              no->next = NULL; 
              no_head = no;
              tail = no_head;
            }else{
              found = 0;
              for(nol = no_head; nol != NULL; nol = nol->next)
              {
                if(strcmp(nol->name, w->parent_object->sym->name) == 0) 
                {
                  nol->orient++;
                  found = 1;
                  break;
                }
              }
              if(!found)
              {
                /* add to the end of list */
                no = CHECKED_MALLOC_STRUCT(
                    struct name_orient, "struct name_orient");
                no->name = CHECKED_STRDUP(
                    w->parent_object->sym->name, "w->parent_object->sym->name");
                no->orient = 1;
                no->next = tail->next;
                tail->next = no;
                tail = tail->next;
              }

            }
                  
         }else if ((smash->what & COLLIDE_SUBVOL) != 0){

              struct subvolume *nsv;

              virt_mol.pos.x = smash->loc.x;
              virt_mol.pos.y = smash->loc.y;
              virt_mol.pos.z = smash->loc.z;

              nsv = traverse_subvol(
                  sv, &(virt_mol.pos),
                  smash->what - COLLIDE_SV_NX - COLLIDE_SUBVOL, state->nx_parts,
                  state->ny_parts, state->nz_parts); 
              if (nsv==NULL)
              {
                if (shead != NULL) mem_put_list(sv->local_storage->coll,shead);

                for(nol = no_head; nol != NULL; nol = nol->next)
                {
                  if(nol->orient%2 != 0) 
                  {
                    best_name = nol->name;
                  }
                }
  
                if(best_name != NULL) {
                   return_name = CHECKED_STRDUP(best_name, "nol->name"); 
                }

                while(no_head != NULL)
                {
                  nnext = no_head->next;
                  free(no_head->name);
                  free(no_head);
                  no_head = nnext;
                }

                if(return_name != NULL) return return_name;
                else return NULL;

              }
        
              if (shead != NULL) mem_put_list(sv->local_storage->coll,shead);
              calculate_random_vector = 0;
              virt_mol.subvol = nsv;

              // Jump to beginning of function
              goto pretend_to_call_find_enclosing_mesh;          
         }
     }
  
  }while(smash != NULL);


  if (shead != NULL) mem_put_list(sv->local_storage->coll,shead);

  for(nol = no_head; nol != NULL; nol = nol->next)
  {
    if(nol->orient%2 != 0) 
    {
       best_name = nol->name;
    }
  }
  
  if(best_name != NULL) {
      return_name = CHECKED_STRDUP(best_name, "nol->name"); 
  }

  while(no_head != NULL)
  {
    nnext = no_head->next;
    free(no_head->name);
    free(no_head);
    no_head = nnext;
  }

  if(return_name != NULL) return return_name;
  else return NULL;

}

/**********************************************************************
place_mol_relative_to_mesh:
  In: 3D location of molecule
      start subvolume
      name of closest enclosing mesh (NULL means that molecule
        is outside of all meshes)
      new position of molecule (return value)
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
void place_mol_relative_to_mesh(
    struct volume *state, struct vector3 *loc, struct subvolume *sv,
    char *mesh_name, struct vector3 *new_pos)
{
  struct wall *best_w = NULL;
  struct wall_list *wl;
  double d2, best_d2;
  struct vector2 s_loc;
  /* struct vector2 best_uv;  */
  struct vector3 best_xyz;
  char *mesh_name_try = NULL;  /* farthest enclosing mesh name */
  struct volume_molecule virt_mol;

  best_w = NULL;
  best_d2 = GIGANTIC  + 1;
  best_xyz.x = 0;
  best_xyz.y = 0;
  best_xyz.z = 0;

  if(mesh_name == NULL)
  {

     /* we have to move molecule that now appeared to be 
        inside the mesh outside of all enclosing meshes */
    virt_mol.prev_v = NULL;
    virt_mol.next_v = NULL;
    virt_mol.next = NULL;
    virt_mol.pos.x = loc->x;
    virt_mol.pos.y = loc->y;
    virt_mol.pos.z = loc->z;
    virt_mol.subvol = find_subvolume(state, loc, NULL);

    mesh_name_try = find_farthest_enclosing_mesh_name(state, &virt_mol);
    if(mesh_name_try == NULL) {
      mcell_internal_error("Cannot find the farthest enclosing mesh.");
    }
  }

  for(wl = sv->wall_head; wl != NULL; wl = wl->next)
  {
    if(mesh_name != NULL)
    {   
      if (strcmp(wl->this_wall->parent_object->sym->name, mesh_name) != 0) {
        continue;
      }
    }else{
      if (strcmp(wl->this_wall->parent_object->sym->name,
                 mesh_name_try) != 0) { continue; }
    }

    d2 = closest_interior_point(loc, wl->this_wall, &s_loc, GIGANTIC);
    if (d2 < best_d2)
    {
       best_d2 = d2;
       best_w = wl->this_wall;
       /* best_uv.u = s_loc.u;
       best_uv.v = s_loc.v; */
    }
  }


  /* look into neighbor subvolumes */
  const int sv_index = sv - state->subvol;
  int sv_remain = sv_index;

  /* Turn linear sv_index into part_x, part_y, part_z triple. */
  const int part_x = sv_remain / ((state->ny_parts-1)*(state->nz_parts-1));
  sv_remain -= part_x * ((state->ny_parts-1)*(state->nz_parts-1));
  const int part_y = sv_remain / (state->nz_parts-1);
  sv_remain -= part_y * (state->nz_parts-1);
  const int part_z = sv_remain;

  /* Find min x partition. */
  int x_min;
  for (x_min=part_x; x_min>0; x_min--)
  {
    d2 = loc->x - state->x_partitions[x_min]; d2 *= d2;
    if (d2 >= best_d2) break;
  }

  /* Find max x partition. */
  int x_max;
  for (x_max=part_x; x_max<state->nx_parts-1 ; x_max++)
  {
    d2 = loc->x - state->x_partitions[x_max + 1]; d2 *= d2;
    if (d2 >= best_d2) break;
  }

  /* Find min y partition. */
  int y_min;
  for (y_min=part_y; y_min>0; y_min--)
  {
    d2 = loc->y - state->y_partitions[y_min]; d2 *= d2;
    if (d2 >= best_d2) break;
  }

  /* Find max y partition. */
  int y_max;
  for (y_max=part_y; y_max<state->ny_parts-1; y_max++)
  {
    d2 = loc->y - state->y_partitions[y_max+1]; d2 *= d2;
    if (d2 >= best_d2) break;
  }

  /* Find min z partition. */
  int z_min;
  for (z_min=part_z; z_min>0; z_min--)
  {
    d2 = loc->z - state->z_partitions[z_min]; d2 *= d2;
    if (d2 >= best_d2) break;
  }

  /* Find max z partition. */
  int z_max;
  for (z_max=part_z; z_max<state->nz_parts-1; z_max++)
  {
    d2 = loc->z - state->z_partitions[z_max+1]; d2 *= d2;
    if (d2 >= best_d2) break;
  }
   

  if (x_min<part_x || x_max>part_x || y_min<part_y || y_max>part_y ||
      z_min<part_z || z_max>part_z)
  {
    for (int px=x_min; px<x_max; px++)
    {
      for (int py=y_min; py<y_max; py++)
      {
        for (int pz=z_min; pz<z_max; pz++)
        {
          const int this_sv = \
              pz + (state->nz_parts-1)*(py + (state->ny_parts-1)*px);
          if (this_sv == sv_index) continue;

          for (wl=state->subvol[this_sv].wall_head; wl!=NULL; wl=wl->next)
          {
            if(mesh_name != NULL)
            {
              if(strcmp(wl->this_wall->parent_object->sym->name,
                        mesh_name) != 0) { continue; }
            } else {
                if(strcmp(wl->this_wall->parent_object->sym->name,
                          mesh_name_try) != 0) { continue; }
            }              
                
            d2 = closest_interior_point(loc,wl->this_wall,&s_loc,GIGANTIC);
            if (d2 < best_d2)
            {
              best_d2 = d2;
              best_w = wl->this_wall;
              /* best_uv.u = s_loc.u;
              best_uv.v = s_loc.v; */
            }
          }
        }
      }
    }
  }

  if (mesh_name_try != NULL) {
    free(mesh_name_try); 
  }
    
  if (best_w!=NULL)
  {
    /* uv2xyz(&best_uv,best_w,&best_xyz);  */
     find_wall_center(best_w, &best_xyz);

  }else mcell_internal_error(
      "Error in function 'place_behind_closest_mesh_position()'.");

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

    if (mesh_name != NULL)
    {
      new_pos->x = best_xyz.x - 2*MESH_DISTINCTIVE*(best_w->normal.x);
      new_pos->y = best_xyz.y - 2*MESH_DISTINCTIVE*(best_w->normal.y);
      new_pos->z = best_xyz.z - 2*MESH_DISTINCTIVE*(best_w->normal.z);

    } else {
      new_pos->x = best_xyz.x + 2*MESH_DISTINCTIVE*(best_w->normal.x);
      new_pos->y = best_xyz.y + 2*MESH_DISTINCTIVE*(best_w->normal.y);
      new_pos->z = best_xyz.z + 2*MESH_DISTINCTIVE*(best_w->normal.z);
    }

}
