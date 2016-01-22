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

#include "config.h"

#include <assert.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <float.h>
#include <time.h>
#include <sys/types.h>
#ifndef _WIN32
#include <sys/resource.h>
#endif

#include "version_info.h"
#include "logging.h"
#include "rng.h"
#include "mcell_structs.h"
#include "sym_table.h"
#include "count_util.h"
#include "vol_util.h"
#include "wall_util.h"
#include "grid_util.h"
#include "viz_output.h"
#include "react.h"
#include "react_output.h"
#include "chkpt.h"
#include "init.h"
#include "mdlparse_aux.h"
#include "mcell_objects.h"
#include "triangle_overlap.h"

#define MESH_DISTINCTIVE EPS_C

struct reschedule_helper {
  struct reschedule_helper *next;
  struct release_event_queue *req;
};

/* Initialize the visualization output (frame_data_lists). */
static int init_viz_output(struct volume *world);

static int compute_bb(struct volume *world, struct object *objp,
                      double (*im)[4]);
static int compute_bb_release_site(struct volume *world, struct object *objp,
                                   double (*im)[4]);
static int compute_bb_polygon_object(struct volume *world, struct object *objp,
                                     double (*im)[4]);

static int init_species_defaults(struct volume *world);
static int init_regions_helper(struct volume *world);

static struct ccn_clamp_data* find_clamped_object_in_list(struct ccn_clamp_data *ccd,
  struct object *obj);

#define MICROSEC_PER_YEAR 365.25 * 86400.0 * 1e6

/* Sets default notification values */
int init_notifications(struct volume *world) {
  if (!(world->notify = CHECKED_MALLOC_STRUCT_NODIE(struct notifications,
                                                    "notification states"))) {
    return 1;
  }

  /* Notifications */
  if (world->quiet_flag) {
    world->notify->progress_report = NOTIFY_NONE;
    world->notify->diffusion_constants = NOTIFY_NONE;
    world->notify->reaction_probabilities = NOTIFY_NONE;
    world->notify->time_varying_reactions = NOTIFY_NONE;
    world->notify->reaction_prob_notify = 0.0;
    world->notify->partition_location = NOTIFY_NONE;
    world->notify->box_triangulation = NOTIFY_NONE;
    world->notify->iteration_report = NOTIFY_NONE;
    world->notify->custom_iteration_value = 0;
    world->notify->release_events = NOTIFY_NONE;
    world->notify->file_writes = NOTIFY_NONE;
    world->notify->final_summary = NOTIFY_NONE;
    world->notify->throughput_report = NOTIFY_NONE;
    world->notify->checkpoint_report = NOTIFY_NONE;
    world->notify->reaction_output_report = NOTIFY_NONE;
    world->notify->volume_output_report = NOTIFY_NONE;
    world->notify->viz_output_report = NOTIFY_NONE;
    world->notify->molecule_collision_report = NOTIFY_NONE;
  } else {
    world->notify->progress_report = NOTIFY_FULL;
    world->notify->diffusion_constants = NOTIFY_BRIEF;
    world->notify->reaction_probabilities = NOTIFY_FULL;
    world->notify->time_varying_reactions = NOTIFY_FULL;
    world->notify->reaction_prob_notify = 0.0;
    world->notify->partition_location = NOTIFY_NONE;
    world->notify->box_triangulation = NOTIFY_NONE;
    world->notify->iteration_report = NOTIFY_FULL;
    world->notify->custom_iteration_value = 0;
    world->notify->release_events = NOTIFY_FULL;
    world->notify->file_writes = NOTIFY_NONE;
    world->notify->final_summary = NOTIFY_FULL;
    world->notify->throughput_report = NOTIFY_FULL;
    world->notify->checkpoint_report = NOTIFY_FULL;
    world->notify->reaction_output_report = NOTIFY_NONE;
    world->notify->volume_output_report = NOTIFY_NONE;
    world->notify->viz_output_report = NOTIFY_NONE;
    world->notify->molecule_collision_report = NOTIFY_NONE;
  }
  /* Warnings */
  world->notify->neg_diffusion = WARN_WARN;
  world->notify->neg_reaction = WARN_WARN;
  world->notify->high_reaction_prob = WARN_COPE;
  world->notify->reaction_prob_warn = 1.0;
  world->notify->close_partitions = WARN_WARN;
  world->notify->degenerate_polys = WARN_WARN;
  world->notify->overwritten_file = WARN_COPE;
  world->notify->short_lifetime = WARN_WARN;
  world->notify->short_lifetime_value = 50;
  world->notify->missed_reactions = WARN_WARN;
  world->notify->missed_reaction_value = 0.001;
  world->notify->missed_surf_orient = WARN_ERROR;
  world->notify->useless_vol_orient = WARN_WARN;
  world->notify->complex_placement_failure_threshold = 25;
  world->notify->complex_placement_failure = WARN_WARN;
  world->notify->mol_placement_failure = WARN_WARN;
  world->notify->invalid_output_step_time = WARN_WARN;

  if (world->log_freq != 0 && world->log_freq != ULONG_MAX) /* User set this */
  {
    world->notify->custom_iteration_value = world->log_freq;
  }

  return 0;
}

/*
 * Initialize the volume data output.  This will create a scheduler for the
 * volume data, and add volume data output items to the scheduler.
 */
static void init_volume_data_output(struct volume *wrld) {
  struct volume_output_item *vo, *vonext;

  wrld->volume_output_scheduler = create_scheduler(
      1.0, 100.0, 100, wrld->simulation_start_seconds / wrld->time_unit);
  if (wrld->volume_output_scheduler == NULL)
    mcell_allocfailed("Failed to create scheduler for volume output data.");

  double r_time_unit = 1.0 / wrld->time_unit;
  for (vo = wrld->volume_output_head; vo != NULL; vo = vonext) {
    vonext = vo->next; /* schedule_add overwrites 'next' */

    if (vo->timer_type == OUTPUT_BY_STEP) {
      if (wrld->chkpt_seq_num == 1)
        vo->t = 0.0;
      else {
        /* Get step time in internal units, find next scheduled output time */
        double f = vo->step_time * r_time_unit;
        vo->t = f * ceil(wrld->volume_output_scheduler->now / f);
      }
    } else if (vo->num_times > 0) {
      /* Set time scaling factor depending on output type */
      double time_scale = 0.0;
      if (vo->timer_type == OUTPUT_BY_ITERATION_LIST)
        time_scale = 1.0;
      else
        time_scale = r_time_unit;

      /* Find the time of next output */
      if (wrld->chkpt_seq_num == 1) {
        vo->next_time = vo->times;
        vo->t = time_scale * *vo->next_time;
      } else /* Scan forward to find first output after checkpoint time */
      {
        int idx = bisect_high(vo->times, vo->num_times,
                              wrld->volume_output_scheduler->now / time_scale);

        /* If we've already passed the last time for this one, skip it! */
        if (idx < 0 || idx >= vo->num_times)
          continue;
        if (wrld->volume_output_scheduler->now / time_scale > vo->times[idx])
          continue;

        vo->t = vo->times[idx] * time_scale;
        vo->next_time = vo->times + idx;
      }

      /* Advance the next_time pointer */
      ++vo->next_time;
    }

    if (schedule_add(wrld->volume_output_scheduler, vo))
      mcell_allocfailed("Failed to add item to schedule for volume output.");
  }
}

/***********************************************************************
 *
 * initialize the state of variables in world
 *
 ***********************************************************************/
int init_variables(struct volume *world) {
  world->t_start = time(NULL);

  if (world->notify->progress_report != NOTIFY_NONE)
    mcell_log("MCell initializing simulation...");

  // XXX: This is in the wrong place here and should be moved
  //      to a separate function perhaps
  install_emergency_output_hooks(world);
  emergency_output_hook_enabled = 0;

  world->curr_file = world->mdl_infile_name;
  world->chkpt_iterations = 0;
  world->last_checkpoint_iteration = 0;
  world->chkpt_seq_num = 0;
  world->keep_chkpts = 0;

  world->last_timing_time = (struct timeval) { 0, 0 };
  world->last_timing_iteration = 0;

  world->chkpt_flag = 0;
  world->viz_blocks = NULL;
  world->ray_voxel_tests = 0;
  world->ray_polygon_tests = 0;
  world->ray_polygon_colls = 0;
  world->vol_vol_colls = 0;
  world->vol_surf_colls = 0;
  world->surf_surf_colls = 0;
  world->vol_wall_colls = 0;
  world->vol_vol_vol_colls = 0;
  world->vol_vol_surf_colls = 0;
  world->vol_surf_surf_colls = 0;
  world->surf_surf_surf_colls = 0;
  world->chkpt_start_time_seconds = 0;
  world->chkpt_byte_order_mismatch = 0;
  world->diffusion_number = 0;
  world->diffusion_cumtime = 0.0;
  world->current_iterations = 0;
  world->elapsed_time = 0;
  world->time_unit = 0;
  world->time_step_max = 0;
  world->start_iterations = 0;
  world->current_time_seconds = 0;
  world->simulation_start_seconds = 0;
  world->grid_density = 10000;
  world->r_length_unit = sqrt(world->grid_density);
  world->length_unit = 1.0 / world->r_length_unit;
  world->rx_radius_3d = 0;
  world->radial_directions = 16384;
  world->radial_subdivisions = 1024;
  world->fully_random = 0;
  world->num_directions = world->radial_directions;
  world->r_step = NULL;
  world->r_step_surface = NULL;
  world->r_step_release = NULL;
  world->d_step = NULL;
  world->dissociation_index = DISSOCIATION_MAX;
  world->place_waypoints_flag = 0;
  world->count_scheduler = NULL;
  world->volume_output_scheduler = NULL;
  world->storage_head = NULL;
  world->storage_allocator = NULL;
  world->x_partitions = NULL;
  world->y_partitions = NULL;
  world->z_partitions = NULL;
  world->x_fineparts = NULL;
  world->y_fineparts = NULL;
  world->z_fineparts = NULL;
  world->n_fineparts = 0;
  world->mem_part_x = 14;
  world->mem_part_y = 14;
  world->mem_part_z = 14;
  world->mem_part_pool = 0;
  world->complex_placement_attempts = 100;
  world->all_vertices = NULL;
  world->walls_using_vertex = NULL;

  world->use_expanded_list = 1;
  world->randomize_smol_pos = 1;
  world->vacancy_search_dist2 = 0.1;
  world->surface_reversibility = 0;
  world->volume_reversibility = 0;
  world->n_reactions = 0;
  world->current_mol_id = 0;

  world->rxn_flags.vol_vol_reaction_flag = 0;
  world->rxn_flags.vol_surf_reaction_flag = 0;
  world->rxn_flags.surf_surf_reaction_flag = 0;
  world->rxn_flags.vol_wall_reaction_flag = 0;
  world->rxn_flags.vol_vol_vol_reaction_flag = 0;
  world->rxn_flags.vol_vol_surf_reaction_flag = 0;
  world->rxn_flags.vol_surf_surf_reaction_flag = 0;
  world->rxn_flags.surf_surf_surf_reaction_flag = 0;

  world->create_shared_walls_info_flag = 0;
  world->reaction_prob_limit_flag = 0;

  world->mcell_version = mcell_version;
  world->clamp_list = NULL;

  return 0;
}

/***********************************************************************
 *
 * initialize the main data structures in world
 *
 ***********************************************************************/
int init_data_structures(struct volume *world) {
  int i;

  if (!(world->rng = CHECKED_MALLOC_STRUCT_NODIE(
            struct rng_state, "random number generator state"))) {
    return 1;
  }

  if (world->seed_seq < 1 || world->seed_seq > INT_MAX)
    mcell_error(
        "Random sequence number must be in the range 1 to 2^31-1 [2147483647]");
  rng_init(world->rng, world->seed_seq);
  if (world->notify->progress_report != NOTIFY_NONE)
    mcell_log("MCell[%d]: random sequence %d", world->procnum, world->seed_seq);

  world->count_hashmask = COUNT_HASHMASK;
  if (!(world->count_hash =
            CHECKED_MALLOC_ARRAY(struct counter *, (world->count_hashmask + 1),
                                 "counter hash table"))) {
    return 1;
  }

  for (i = 0; i <= world->count_hashmask; i++)
    world->count_hash[i] = NULL;

  world->oexpr_mem = create_mem_named(sizeof(struct output_expression), 128,
                                      "output expression");
  if (world->oexpr_mem == NULL) {
    mcell_allocfailed_nodie("Failed to create memory pool for reaction data "
                            "output expressions.");
    return 1;
  }

  world->outp_request_mem =
      create_mem_named(sizeof(struct output_request), 64, "output request");
  if (world->outp_request_mem == NULL) {
    mcell_allocfailed_nodie("Failed to create memory pool for reaction data "
                            "output commands.");
    return 1;
  }

  world->counter_mem = create_mem_named(sizeof(struct counter), 32, "counter");
  if (world->counter_mem == NULL) {
    mcell_allocfailed_nodie("Failed to create memory pool for reaction and "
                            "molecule counts.");
    return 1;
  }

  world->trig_request_mem =
      create_mem_named(sizeof(struct trigger_request), 32, "trigger request");
  if (world->trig_request_mem == NULL) {
    mcell_allocfailed_nodie("Failed to create memory pool for reaction and "
                            "molecule output triggers.");
    return 1;
  }

  world->magic_mem = create_mem_named(sizeof(struct magic_list), 1024,
                                      "reaction-triggered release");
  if (world->magic_mem == NULL) {
    mcell_allocfailed_nodie("Failed to create memory pool for "
                            "reaction-triggered release lists.");
    return 1;
  }

  if ((world->fstream_sym_table = init_symtab(1024)) == NULL) {
    mcell_allocfailed_nodie(
        "Failed to initialize symbol table for file streams.");
    return 1;
  }

  if ((world->mol_sym_table = init_symtab(1024)) == NULL) {
    mcell_allocfailed_nodie("Failed to initialize symbol table for molecules.");
    return 1;
  }

  if ((world->rxn_sym_table = init_symtab(1024)) == NULL) {
    mcell_allocfailed_nodie("Failed to initialize symbol table for reactions.");
    return 1;
  }

  if ((world->obj_sym_table = init_symtab(1024)) == NULL) {
    mcell_allocfailed_nodie("Failed to initialize symbol table for objects.");
    return 1;
  }

  if ((world->reg_sym_table = init_symtab(1024)) == NULL) {
    mcell_allocfailed_nodie("Failed to initialize symbol table for regions.");
    return 1;
  }

  if ((world->rpat_sym_table = init_symtab(1024)) == NULL) {
    mcell_allocfailed_nodie(
        "Failed to initialize symbol table for release patterns.");
    return 1;
  }

  if ((world->rxpn_sym_table = init_symtab(1024)) == NULL) {
    mcell_allocfailed_nodie(
        "Failed to initialize symbol table for reaction pathways.");
    return 1;
  }

  struct sym_entry *sym;
  if ((sym = store_sym("WORLD_OBJ", OBJ, world->obj_sym_table, NULL)) == NULL) {
    mcell_allocfailed_nodie(
        "Failed to store the world root object in the object symbol table.");
    return 1;
  }

  world->root_object = (struct object *)sym->value;
  world->root_object->object_type = META_OBJ;
  if (!(world->root_object->last_name = CHECKED_STRDUP_NODIE("", NULL))) {
    return 1;
  }

  if ((sym = store_sym("WORLD_INSTANCE", OBJ, world->obj_sym_table, NULL)) ==
      NULL) {
    mcell_allocfailed_nodie(
        "Failed to store the world root instance in the object symbol table.");
    return 1;
  }

  world->root_instance = (struct object *)sym->value;
  world->root_instance->object_type = META_OBJ;
  if (!(world->root_instance->last_name = CHECKED_STRDUP("", NULL))) {
    return 1;
  }

  if ((sym = store_sym("DEFAULT_RELEASE_PATTERN", RPAT, world->rpat_sym_table,
                      NULL)) == NULL) {
    mcell_allocfailed_nodie("Failed to store the default release pattern in "
                            "the release patterns symbol table.");
    return 1;
  }
  world->default_release_pattern = (struct release_pattern *)sym->value;
  world->default_release_pattern->delay = 0;
  world->default_release_pattern->release_interval = FOREVER;
  world->default_release_pattern->train_interval = FOREVER;
  world->default_release_pattern->train_duration = FOREVER;
  world->default_release_pattern->number_of_trains = 1;

  if ((sym = store_sym("ALL_VOLUME_MOLECULES", MOL, world->mol_sym_table,
                      NULL)) == NULL) {
    mcell_allocfailed_nodie(
        "Failed to store all volume molecules in the molecule symbol table.");
    return 1;
  }
  world->all_volume_mols = (struct species *)sym->value;

  if ((sym = store_sym("ALL_SURFACE_MOLECULES", MOL, world->mol_sym_table,
                      NULL)) == NULL) {
    mcell_allocfailed_nodie(
        "Failed to store the surface molecules in the molecule symbol table.");
    return 1;
  }
  world->all_surface_mols = (struct species *)sym->value;

  if ((sym = store_sym("ALL_MOLECULES", MOL, world->mol_sym_table, NULL)) ==
      NULL) {
    mcell_allocfailed_nodie(
        "Failed to store all molecules in the molecule symbol table.");
    return 1;
  }
  world->all_mols = (struct species *)sym->value;

  world->volume_output_head = NULL;
  world->output_block_head = NULL;
  world->output_request_head = NULL;

  world->releaser = create_scheduler(1.0, 100.0, 100, 0.0);
  if (world->releaser == NULL) {
    mcell_allocfailed_nodie("Failed to create release scheduler.");
    return 1;
  }

  return 0;
}

/***********************************************************************
 *
 * parse the model's mdl files and update our global state
 *
 ***********************************************************************/
int parse_input(struct volume *world) {
  /* Parse the MDL file: */
  no_printf("Node %d parsing MDL file %s\n", world->procnum,
            world->mdl_infile_name);
  if (mdlparse_init(world)) {
    return (1);
  }
  no_printf("Done parsing MDL file: %s\n", world->mdl_infile_name);

  /* we do not want to count collisions if the policy is not to print */
  if (world->notify->final_summary == NOTIFY_NONE)
    world->notify->molecule_collision_report = NOTIFY_NONE;

  if (world->iterations == INT_MIN) {
    mcell_error_nodie(
        "Total number of iterations is not specified either "
        "through the ITERATIONS keyword or through the command line option "
        "'-iterations'.");
    return 1;
  }

  return 0;
}

/***********************************************************************
 *
 * initialize the models' species table
 *
 ***********************************************************************/
int init_species(struct volume *world) {

  int reactants_3D_present = 0; /* flag to check whether there are 3D reactants
                             (participants in the reactions
                              between 3D molecules) in the simulation */

  /* Set up the array of species */
  if (init_species_defaults(world)) {
    mcell_error_nodie("Unknown error while initializing species table.");
    return 0;
  }
  no_printf("Done setting up species.\n");

  /* Create linked lists of volume and surface molecules names */
  struct name_list *vol_species_name_list = NULL;
  struct name_list *surf_species_name_list = NULL;
  if (world->notify->reaction_probabilities == NOTIFY_FULL) {
    create_name_lists_of_volume_and_surface_mols(world, &vol_species_name_list,
                                                 &surf_species_name_list);
  }

  for (int i = 0; i < world->n_species; i++) {
    struct species *sp = world->species_list[i];

    if (sp->flags & IS_SURFACE) {
      check_for_conflicts_in_surface_class(world, sp);
      if (world->notify->reaction_probabilities == NOTIFY_FULL) {
        publish_special_reactions_report(sp, vol_species_name_list,
                                         surf_species_name_list,
                                         world->n_species, world->species_list);
      }
    }
  }

  /* Memory deallocate linked lists of volume molecules names and surface
     molecules names */
  if (vol_species_name_list != NULL)
    remove_molecules_name_list(&vol_species_name_list);
  if (surf_species_name_list != NULL)
    remove_molecules_name_list(&surf_species_name_list);

  /* If there are no 3D molecules-reactants in the simulation
     set up the "use_expanded_list" flag to zero. */
  for (int i = 0; i < world->n_species; i++) {
    struct species *sp = world->species_list[i];
    if (sp == world->all_mols)
      continue;
    if (sp == world->all_volume_mols)
      continue;
    if (sp == world->all_surface_mols)
      continue;
    if (sp->flags & ON_GRID)
      continue;
    if (sp->flags & IS_SURFACE)
      continue;

    if ((sp->flags & (CAN_VOLVOL | CAN_VOLVOLVOL)) != 0) {
      reactants_3D_present = 1;
      break;
    }
  }

  if (reactants_3D_present == 0) {
    world->use_expanded_list = 0;
  }

  return 0;
}

/***********************************************************************
 *
 * initialize the models' vertices and walls
 *
 ***********************************************************************/
int init_vertices_walls(struct volume *world) {
  struct storage_list *sl;
  int num_storages = 0;           /* number of storages in the world */
  int *num_vertices_this_storage; /* array of vertex counts belonging to each
                                     storage */

  for (sl = world->storage_head; sl != NULL; sl = sl->next) {
    num_storages++;
  }

  /* Allocate array of counts (note: the ordering of this array follows the
   * ordering of the linked list "world->storage_list") */
  if (!(num_vertices_this_storage = CHECKED_MALLOC_ARRAY_NODIE(
            int, num_storages, "array of vertex counts per storage"))) {
    return 1;
  }
  memset(num_vertices_this_storage, 0, sizeof(int) * num_storages);

  double tm[4][4];
  init_matrix(tm);

  /* Accumulate vertex counts and rescale each vertex coordinates if needed */
  if (accumulate_vertex_counts_per_storage(world, world->root_instance,
                                           num_vertices_this_storage, tm)) {
    mcell_error_nodie("Error while accumulating vertex counts per storage.");
    return 1;
  }

  /* Cumulate the vertex count */
  for (int kk = 1; kk < num_storages; ++kk) {
    num_vertices_this_storage[kk] += num_vertices_this_storage[kk - 1];
  }

  /* Allocate the  global "all_vertices" array */
  if (!(world->all_vertices = CHECKED_MALLOC_ARRAY_NODIE(
            struct vector3, num_vertices_this_storage[num_storages - 1],
            "array of all vertices in the world"))) {
    return 1;
  }

  /* Allocate the global "walls_using_vertex" array */
  if (world->create_shared_walls_info_flag) {
    if (!(world->walls_using_vertex = CHECKED_MALLOC_ARRAY_NODIE(
              struct wall_list *, num_vertices_this_storage[num_storages - 1],
              "wall list pointers"))) {
      return 1;
    }

    for (int kk = 0; kk < num_vertices_this_storage[num_storages - 1]; kk++) {
      world->walls_using_vertex[kk] = NULL;
    }
  }
  init_matrix(tm);

  /* Copy vertices into the global array "world->all_vertices"
    and fill "objp->vertices"  for each object in the world */
  if (fill_world_vertices_array(world, world->root_instance,
                                num_vertices_this_storage, tm))
    return 1;
  /* free memory */
  free(num_vertices_this_storage);

  init_matrix(tm);
  /* Instantiate all objects */
  if (world->notify->progress_report != NOTIFY_NONE)
    mcell_log("Instantiating objects...");
  if (instance_obj(world, world->root_instance, tm))
    return 1;

  if (world->notify->progress_report != NOTIFY_NONE)
    mcell_log("Creating walls...");
  if (distribute_world(world)) {
    mcell_error_nodie("Unknown error while distributing geometry "
                      "among partitions.");
    return 1;
  }

  if (world->notify->progress_report != NOTIFY_NONE)
    mcell_log("Creating edges...");
  if (sharpen_world(world)) {
    mcell_error_nodie("Unknown error while adding edges to geometry.");
    return 1;
  }

  return 0;
}

/***********************************************************************
 *
 * initialize the models' regions
 *
 ***********************************************************************/
int init_regions(struct volume *world) {
  if (prepare_counters(world)) {
    mcell_error_nodie(
        "Unknown error while preparing counters for reaction data output.");
    return 1;
  }

  if (init_regions_helper(world)) {
    mcell_error_nodie("Unknown error while initializing object regions.");
    return 1;
  }

  if (check_counter_geometry(world->count_hashmask, world->count_hash,
                             &world->place_waypoints_flag)) {
    mcell_error_nodie(
        "Unknown error while validating geometry of counting regions.");
    return 1;
  }

  /* flags that tell whether there are regions set with surface classes
     that contain ALL_MOLECULES or ALL_SURFACE_MOLECULES keywords.*/
  int all_mols_region_present = 0, all_surf_mols_region_present = 0;
  struct species *all_mols_sp = get_species_by_name(
      "ALL_MOLECULES", world->n_species, world->species_list);
  struct species *all_surf_mols_sp = get_species_by_name(
      "ALL_SURFACE_MOLECULES", world->n_species, world->species_list);

  if ((all_mols_sp != NULL) && (all_mols_sp->flags & REGION_PRESENT)) {
    all_mols_region_present = 1;
  }
  if ((all_surf_mols_sp != NULL) &&
      (all_surf_mols_sp->flags & REGION_PRESENT)) {
    all_surf_mols_region_present = 1;
  }

  /* if surface molecules are defined as part of SURFACE_CLASS
     definitions but there are no regions with this SURFACE_CLASS
     in the model, then to speed up model simulation remove
     the flag CAN_REGION_BORDER from the surface molecule */
  if ((!all_mols_region_present) && (!all_surf_mols_region_present)) {
    for (int i = 0; i < world->n_species; i++) {
      struct species *sp = world->species_list[i];
      if (sp->flags & ON_GRID) {
        if ((sp->flags & CAN_REGION_BORDER) &&
            ((sp->flags & REGION_PRESENT) == 0)) {
          sp->flags &= ~CAN_REGION_BORDER;
        }
      }
    }
  }

  return 0;
}

/**********************************************************************
 *
 * initialize the hash which provides a mapping from the name of
 * count statments (currently identical to the name of the output
 * files) to the underlying output_block data structure containing
 * the counts.
 *
 **********************************************************************/
int init_counter_name_hash(struct sym_table_head **counter_by_name,
                           struct output_block *output_block_head) {
  *counter_by_name = init_symtab(2048);
  if (*counter_by_name == NULL) {
    mcell_log("error creating count symbol table");
    return 1;
  }

  // insert count data
  for (struct output_block *out_block = output_block_head;
       out_block != NULL; out_block = out_block->next) {
    for (struct output_set *set = out_block->data_set_head; set != NULL;
         set = set->next) {
      store_sym(set->outfile_name, COUNT_OBJ_PTR, *counter_by_name, set);
    }
  }

  return 0;
}

/***********************************************************************
 *
 * load the model checkpoint data
 *
 ***********************************************************************/
int load_checkpoint(struct volume *world) {
  FILE *chkpt_infs = NULL;
  if ((chkpt_infs = fopen(world->chkpt_infile, "rb")) == NULL) {
    world->chkpt_seq_num = 1;
  } else {
    mcell_log("Reading from checkpoint file '%s'.", world->chkpt_infile);
    if (read_chkpt(world, chkpt_infs)) {
      mcell_error_nodie("Failed to read checkpoint file '%s'.",
                        world->chkpt_infile);
      return 1;
    }
    fclose(chkpt_infs);
  }

  return 0;
}

/***********************************************************************
 *
 * initialize the model's viz data output
 *
 ***********************************************************************/
int init_viz_data(struct volume *world) {
  /* Initialize the frame data for the visualization and reaction output. */
  if (init_viz_output(world)) {
    mcell_error_nodie("Unknown error while initializing VIZ output.");
    return 1;
  }

  /* Initialize the volume output */
  init_volume_data_output(world);

  return 0;
}

/***********************************************************************
 *
 * initialize the model's reaction data output
 *
 ***********************************************************************/
int init_reaction_data(struct volume *world) {
  struct output_block *obp, *obpn;
  struct output_set *set;
  double f;

  world->count_scheduler = create_scheduler(1.0, 100.0, 100, world->start_iterations);
  if (world->count_scheduler == NULL) {
    mcell_allocfailed_nodie(
        "Failed to create scheduler for reaction data output.");
    return 1;
  }

  /* Schedule the reaction data output events */
  obp = world->output_block_head;
  while (obp != NULL) {
    obpn = obp->next; /* Save this--will be lost when we schedule obp */

    if (obp->timer_type == OUTPUT_BY_STEP) {
      if (world->chkpt_seq_num == 1)
        obp->t = 0.0;
      else {
        double stepInt = obp->step_time/world->time_unit; /* Step time (internal units) */
        double start = world->count_scheduler->now;
        start += stepInt - fmod(start, stepInt);
        obp->t = start;
      }
    } else if (obp->time_now == NULL) /* When would this be non-NULL?? */
    {
      /* Set time scaling factor depending on output type */
      if (obp->timer_type == OUTPUT_BY_ITERATION_LIST)
        f = 1.0;
      else
        f = 1.0 / world->time_unit;

      /* Find the time of next output */
      if (world->chkpt_seq_num == 1) {
        obp->time_now = obp->time_list_head;
        obp->t = f * obp->time_now->value;
      } else /* Scan forward to find first output after checkpoint time */
      {
        for (obp->time_now = obp->time_list_head; obp->time_now != NULL;
             obp->time_now = obp->time_now->next) {
          if (obp->timer_type == OUTPUT_BY_ITERATION_LIST) {
            obp->t = f * obp->time_now->value;
            if (!(obp->t < world->iterations + 1 &&
                  obp->t <= world->count_scheduler->now))
              break;
          } else if (obp->timer_type == OUTPUT_BY_TIME_LIST) {
            if (obp->time_now->value > world->simulation_start_seconds) {
              obp->t = world->count_scheduler->now +
                       (obp->time_now->value - world->simulation_start_seconds) /
                           world->time_unit;
              break;
            }
          }
        }
      }
    }

    for (set = obp->data_set_head; set != NULL; set = set->next) {
      if (set->file_flags == FILE_SUBSTITUTE) {
        if (world->chkpt_seq_num == 1) {
          FILE *file = fopen(set->outfile_name, "w");
          if (file == NULL) {
            mcell_perror_nodie(errno, "Failed to open reaction data output "
                                      "file '%s' for writing",
                               set->outfile_name);
            return 1;
          }

          fclose(file);
        } else if (obp->timer_type == OUTPUT_BY_ITERATION_LIST) {
          if (obp->time_now == NULL)
            continue;
          if (truncate_output_file(set->outfile_name, obp->t)) {
            mcell_error_nodie("Failed to prepare reaction data output file "
                              "'%s' to receive output.",
                              set->outfile_name);
            return 1;
          }
        } else if (obp->timer_type == OUTPUT_BY_TIME_LIST) {
          if (obp->time_now == NULL)
            continue;
          if (truncate_output_file(set->outfile_name,
                                   obp->t * world->time_unit)) {
            mcell_error_nodie("Failed to prepare reaction data output file "
                              "'%s' to receive output.",
                              set->outfile_name);
            return 1;
          }
        } else {
          /* we need to truncate up until the start of the new checkpoint
            * simulation plus a single TIMESTEP */
          double startTime =
              world->chkpt_start_time_seconds + world->time_unit;
          if (truncate_output_file(set->outfile_name, startTime)) {
            mcell_error_nodie("Failed to prepare reaction data output file "
                              "'%s' to receive output.",
                              set->outfile_name);
            return 1;
          }
        }
      }
    }

    if (schedule_add(world->count_scheduler, obp)) {
      mcell_allocfailed_nodie(
          "Failed to add reaction data output item to scheduler.");
      return 1;
    }
    obp = obpn;
  }

  return 0;
}

/***********************************************************************
 *
 * initialize the model's timers
 *
 ***********************************************************************/
int init_timers(struct volume *world) {
  struct rusage init_time = { .ru_utime = { 0, 0 }, .ru_stime = { 0, 0 } };
  getrusage(RUSAGE_SELF, &init_time);

  world->u_init_time.tv_sec = init_time.ru_utime.tv_sec;
  world->u_init_time.tv_usec = init_time.ru_utime.tv_usec;
  world->s_init_time.tv_sec = init_time.ru_stime.tv_sec;
  world->s_init_time.tv_usec = init_time.ru_stime.tv_usec;

  no_printf("Done initializing simulation\n");
  return 0;
}

/***********************************************************************
 *
 * initialize the model's checkpoint state
 *
 ***********************************************************************/
int init_checkpoint_state(struct volume *world, long long *exec_iterations) {
  if (world->notify->checkpoint_report != NOTIFY_NONE)
    mcell_log("MCell: checkpoint sequence number %d begins at elapsed "
              "time %1.15g seconds",
              world->chkpt_seq_num, world->chkpt_start_time_seconds);

  if (world->iterations < world->start_iterations) {
    mcell_error_nodie("Start time after checkpoint %lld is greater than "
                      "total number of iterations specified %lld.",
                      world->start_iterations, world->iterations);
    return 1;
  }

  if (world->chkpt_iterations &&
      (world->iterations - world->start_iterations) < world->chkpt_iterations) {
    world->chkpt_iterations = world->iterations - world->start_iterations;
  }

  if (world->chkpt_iterations)
    *exec_iterations = world->chkpt_iterations;
  else if (world->chkpt_infile)
    *exec_iterations = world->iterations - world->start_iterations;
  else
    *exec_iterations = world->iterations;
  if (*exec_iterations < 0) {
    mcell_error_nodie(
        "Number of iterations to execute is zero or negative. "
        "Please verify ITERATIONS and/or CHECKPOINT_ITERATIONS commands.");
    return 1;
  }
  if (world->notify->progress_report != NOTIFY_NONE)
    mcell_log(
        "MCell: executing %lld iterations starting at iteration number %lld.",
        *exec_iterations, world->start_iterations);

  return 0;
}

/*****************************************************************************
 *
 * reschedule_release_events reschedules release events during restarts
 * from a checkpoint.
 *
 * During parse time and initialization, release events (e.g. via release
 * patterns combined with a delay) are scheduled based on the assumption that
 * the timestep throughout the simulation was consistent and identical to the
 * one in the parsed mdl file. This assumption may be wrong for restarts from
 * a checkpoint in which the timestep was changed with respect to previous runs
 * thus resulting in release events being scheduled at a wrong internal time.
 * The current function computes the proper internal release time based on
 * the start time of the checkpoint and the current iteration number and then
 * reschedules all future events accordingly.
 *
 * This function returns 0 on success and 1 otherwise.
 ******************************************************************************/
int reschedule_release_events(struct volume *world) {

  struct reschedule_helper *helper = NULL;
  for (struct schedule_helper *sh = world->releaser; sh != NULL; sh = sh->next_scale) {
    for (int i = -1; i < sh->buf_len; i++) {
      for (struct abstract_element *ae = (i == -1) ? sh->current : sh->circ_buf_head[i];
        ae != NULL; ae = ae->next) {
        struct release_event_queue *req = (struct release_event_queue *)ae;
        struct reschedule_helper *tmp =
          CHECKED_MALLOC_STRUCT(struct reschedule_helper,
            "Error creating reschedule helper");
        if (helper == NULL) {
          tmp->next = NULL;
        } else {
          tmp->next = helper;
        }
        tmp->req = req;
        helper = tmp;
      }
    }
  }

  struct reschedule_helper *rh = helper;
  while (rh != NULL) {
    struct release_event_queue *req = rh->req;
    rh = rh->next;

    // adjust event time
    double sched_time = req->event_time * world->time_unit;
    double real_sched_time = convert_seconds_to_iterations(
        world->start_iterations, world->time_unit,
        world->chkpt_start_time_seconds, sched_time);

    // adjust time of start of train
    double train_time = req->train_high_time * world->time_unit;
    req->train_high_time = convert_seconds_to_iterations(
        world->start_iterations, world->time_unit,
        world->chkpt_start_time_seconds, train_time);

    schedule_reschedule(world->releaser, req, real_sched_time);
  }

  while (helper != NULL) {
    rh = helper->next;
    free(helper);
    helper = rh;
  }

  return 0;
}

/*************************************************************************
 Mark an object and all of its children for inclusion in a particular viz
 output block.

 In: vizblk: the viz output block in which to include the object
     objp: the object to include
     viz_state: the desired viz state
 Out: No return value.  vizblk is updated.
*************************************************************************/
static void set_viz_state_include(struct viz_output_block *vizblk,
                                  struct object *objp, int viz_state) {
  struct sym_entry *symp;
  struct viz_child *vcp;
  switch (objp->object_type) {
  case META_OBJ:
    for (struct object *child_objp = objp->first_child; child_objp != NULL;
         child_objp = child_objp->next)
      set_viz_state_include(vizblk, child_objp, viz_state);
    break;

  case BOX_OBJ:
  case POLY_OBJ:
    symp = retrieve_sym(objp->sym->name, vizblk->viz_children);
    if (symp == NULL) {
      vcp =
          CHECKED_MALLOC_STRUCT(struct viz_child, "visualization child object");
      vcp->obj = objp;
      vcp->viz_state = NULL;
      vcp->next = NULL;
      vcp->parent = NULL;
      vcp->children = NULL;
      if (store_sym(objp->sym->name, VIZ_CHILD, vizblk->viz_children, vcp) ==
          NULL)
        mcell_allocfailed("Failed to store visualization child object in the "
                          "visualization objects table.");
    } else
      vcp = (struct viz_child *)symp->value;

    if (vcp->viz_state == NULL) {
      vcp->viz_state =
          CHECKED_MALLOC_ARRAY(int, objp->n_walls, "visualization state array");
      for (int i = 0; i < objp->n_walls; i++)
        vcp->viz_state[i] = viz_state;
    } else {
      /* Do not override any specific viz states already set. */
      for (int i = 0; i < objp->n_walls; i++)
        if (vcp->viz_state[i] == EXCLUDE_OBJ)
          vcp->viz_state[i] = viz_state;
    }
    break;

  case REL_SITE_OBJ:
    /* just do nothing */
    break;

  case VOXEL_OBJ:
  default:
    mcell_internal_error("Invalid object type (%d) while setting viz state.",
                         objp->object_type);
  }
}

/*************************************************************************
 Mark all mesh objects for inclusion in the specified viz output block.

 In: vizblk: the viz output block in which to include the object
     viz_state: the desired viz state
 Out: No return value.  vizblk is updated.
*************************************************************************/
static void set_viz_all_meshes(struct object *root_instance,
                               struct viz_output_block *vizblk, int viz_state) {
  set_viz_state_include(vizblk, root_instance, viz_state);
}

/*************************************************************************
 Mark all molecule objects for inclusion in the specified viz output block.

 In: vizblk: the viz output block in which to include the object
     viz_state: the visualization state desired
 Out: No return value.  vizblk is updated.
*************************************************************************/
static void set_viz_all_molecules(struct volume *world,
                                  struct viz_output_block *vizblk,
                                  int viz_state) {
  for (int i = 0; i < world->n_species; i++) {
    struct species *sp = world->species_list[i];
    if (sp->flags & IS_SURFACE)
      continue;
    if (sp == world->all_mols)
      continue;
    if (sp == world->all_volume_mols)
      continue;
    if (sp == world->all_surface_mols)
      continue;
    if (vizblk->species_viz_states[i] != EXCLUDE_OBJ)
      continue;

    /* set viz_state to INCLUDE_OBJ for the molecule we want to visualize
       but will not assign state value */
    vizblk->species_viz_states[i] = viz_state;
  }
}

/*************************************************************************
 Free all viz_child objects which represent either meta objects, or unrendered
 mesh objects.

 In: vizblk: the viz output block whose children should be trimmed
 Out: No return value; parse-time viz output data structures are freed, as are
      any excess viz_child objects.
*************************************************************************/
static void free_extra_viz_children(struct viz_output_block *vizblk) {
  for (int i = 0; i < vizblk->viz_children->n_bins; ++i) {
    for (struct sym_entry *sym = vizblk->viz_children->entries[i]; sym != NULL;
         sym = sym->next) {
      free(sym->name);
      struct viz_child *vcp = (struct viz_child *)sym->value;
      vcp->next = NULL;
      vcp->parent = NULL;
      vcp->children = NULL;
      if (vcp->viz_state == NULL)
        free(vcp);
    }
  }

  destroy_symtab(vizblk->viz_children);
  vizblk->viz_children = NULL;
}

/*************************************************************************
 Comparison function for viz_children, suitable for use with qsort.  Used to
 order viz_child objects in ascending order of their 'obj' pointers, which
 allows us to use binary search to find if an object is included.  (DREAMM
 modes only).

 In: vc1: first child for comparison
     vc2: second child for comparison
 Out: -1, 0, or 1 if first child is lt, eq, or gt the second child, resp.
*************************************************************************/
static int viz_child_compare(void const *vc1, void const *vc2) {
  struct viz_child const *c1 = *(struct viz_child const **)vc1;
  struct viz_child const *c2 = *(struct viz_child const **)vc2;
  if (c1->obj < c2->obj)
    return -1;
  else if (c1->obj > c2->obj)
    return 1;
  else
    return 0;
}

/*************************************************************************
 Convert the viz_child objects in the given viz_output block into an array,
 sorted in ascending order of their 'obj' pointers, excluding any viz_child
 objects which are not marked for inclusion.

 In: vizblk: the viz output block whose children to copy to an array.
 Out: vizblk is updated (dreamm_objects and n_dreamm_objects are set).
*************************************************************************/
static void convert_viz_objects_to_array(struct viz_output_block *vizblk) {
  int count = 0;
  for (int i = 0; i < vizblk->viz_children->n_bins; ++i) {
    for (struct sym_entry *sym = vizblk->viz_children->entries[i]; sym != NULL;
         sym = sym->next) {
      struct viz_child *vcp = (struct viz_child *)sym->value;
      if (vcp->viz_state)
        ++count;
    }
  }

  /* Stash info in the visualization block. */
  vizblk->n_dreamm_objects = count;
  vizblk->dreamm_object_info =
      CHECKED_MALLOC_ARRAY(struct viz_child *, count, "DREAMM mesh objects");
  vizblk->dreamm_objects =
      CHECKED_MALLOC_ARRAY(struct object *, count, "DREAMM mesh objects");

  /* Now copy data in, unsorted. */
  count = 0;
  for (int i = 0; i < vizblk->viz_children->n_bins; ++i) {
    for (struct sym_entry *sym = vizblk->viz_children->entries[i]; sym != NULL;
         sym = sym->next) {
      struct viz_child *vcp = (struct viz_child *)sym->value;
      if (vcp->viz_state)
        vizblk->dreamm_object_info[count++] = vcp;
    }
  }
  assert(count == vizblk->n_dreamm_objects);

  /* Sort the data. */
  qsort(vizblk->dreamm_object_info, vizblk->n_dreamm_objects,
        sizeof(struct viz_child *), viz_child_compare);

  /* Copy out just the objects. */
  for (int i = 0; i < count; ++i)
    vizblk->dreamm_objects[i] = vizblk->dreamm_object_info[i]->obj;
}

/*************************************************************************
 Expand the mesh info for all viz children on the given viz output block.
 (DREAMM/DX modes only).

 In: vizblk: the viz output block whose children to expand
 Out: vizblk is updated
*************************************************************************/
static void expand_viz_children(struct viz_output_block *vizblk) {
  switch (vizblk->viz_mode) {
  case DREAMM_V3_GROUPED_MODE:
  case DREAMM_V3_MODE:
    convert_viz_objects_to_array(vizblk);
    free_extra_viz_children(vizblk);
    break;

  case NO_VIZ_MODE:
  case ASCII_MODE:
  case CELLBLENDER_MODE:
  default:
    /* Do nothing. */
    break;
  }
}

/*************************************************************************
 Initialize the species state array for a given viz output block.

 In: vizblk: the viz output block whose species table to update
 Out: vizblk is updated
*************************************************************************/
static int init_viz_species_states(int n_species,
                                   struct viz_output_block *vizblk) {
  vizblk->species_viz_states =
      CHECKED_MALLOC_ARRAY(int, n_species, "species viz states array");
  if (vizblk->species_viz_states == NULL)
    return 1;
  for (int i = 0; i < n_species; ++i)
    vizblk->species_viz_states[i] = EXCLUDE_OBJ;

  int n_entries = vizblk->parser_species_viz_states.num_items;
  int n_bins = vizblk->parser_species_viz_states.table_size;
  for (int i = 0; n_entries > 0 && i < n_bins; ++i) {
    struct species *specp =
        (struct species *)(vizblk->parser_species_viz_states.keys[i]);
    if (specp != NULL) {
      int viz_state =
          (int)(intptr_t)vizblk->parser_species_viz_states.values[i];

      vizblk->species_viz_states[specp->species_id] = viz_state;
      --n_entries;
    }
  }

  pointer_hash_destroy(&vizblk->parser_species_viz_states);
  return 0;
}

/*************************************************************************
 Initialize all viz output blocks for this simulation.

 In: None.
 Out: 0 on success, 1 if an error occurs
*************************************************************************/
static int init_viz_output(struct volume *world) {
  for (struct viz_output_block *vizblk = world->viz_blocks; vizblk != NULL;
       vizblk = vizblk->next) {
    /* Copy species states into an array. */
    if (init_viz_species_states(world->n_species, vizblk))
      return 1;

    /* If ALL_MESHES or ALL_MOLECULES were requested, mark them all for
     * inclusion. */
    if (vizblk->viz_output_flag & VIZ_ALL_MESHES)
      set_viz_all_meshes(world->root_instance, vizblk, vizblk->default_mesh_state);
    if (vizblk->viz_output_flag & VIZ_ALL_MOLECULES)
      set_viz_all_molecules(world, vizblk, vizblk->default_mol_state);

    /* Copy viz children to the appropriate array. */
    expand_viz_children(vizblk);

    /* Initialize each data frame in this block. */
    if (init_frame_data_list(world, vizblk)) {
      mcell_internal_error("Unknown error while initializing VIZ output.");
      /*return 1;*/
    }
  }

  return 0;
}

/********************************************************************
init_species:
   Initializes array of molecules types to the default properties values.

*********************************************************************/
static int init_species_defaults(struct volume *world) {

  world->speed_limit = 0;

  world->n_species = world->mol_sym_table->n_entries;
  world->species_list =
      CHECKED_MALLOC_ARRAY(struct species *, world->n_species, "species table");
  unsigned int count = 0;
  for (int i = 0; i < world->mol_sym_table->n_bins; i++) {
    for (struct sym_entry *sym = world->mol_sym_table->entries[i];
         sym != NULL; sym = sym->next) {
      if (sym->sym_type == MOL) {
        struct species *s = (struct species *)sym->value;
        world->species_list[count] = s;
        world->species_list[count]->species_id = count;
        world->species_list[count]->chkpt_species_id = UINT_MAX;
        world->species_list[count]->population = 0;
        world->species_list[count]->n_deceased = 0;
        world->species_list[count]->cum_lifetime_seconds = 0;

        if (!(world->species_list[count]->flags & SET_MAX_STEP_LENGTH)) {
          world->species_list[count]->max_step_length = DBL_MAX;
        }

        // If volume molecule, set max speed per time step.
        if ((s->flags & NOT_FREE) == 0) {
          double speed = 6.0 * s->space_step / sqrt(MY_PI);
          if (speed > world->speed_limit)
            world->speed_limit = speed;
        }
        count++;
      }
    }
  }

  return 0;
}

/********************************************************************
 create_storage:

    This is a helper function to create the storage associated with a
    particular subdivision (i.e. an IxJxK box of subvolumes with a common
    scheduler and common memory allocation).

    When we create a storage, we need to decide how big the memory pools should
    be.

    If the memory pools are too small, you lose some of the benefits of memory
    pooling, since your memory blocks are more scattered.  You also take on
    extra overhead because of the extra arenas you end up allocating when you
    need more objects than the original pool could provide.

    If the memory pools are too large, you waste a lot of memory, as the memory
    pool will allocate far more objects than you end up using.  This can be a
    serious issue.

    In the best of all possible worlds, we'd tune the memory pools to coincide
    exactly with the peak usage of each type of object in that subdivision.
    Unfortunately, this would require uncanny prescience.

    In the absence of such foresight, this code uses a fairly dumb heuristic
    for deciding how large the memory pools should be in a subdivision.  It
    uses the number of subvolumes in the subdivision to set the size.  This
    seems, in practice, to be adequate.  It's almost certain that we can
    improve upon it.

    A fairly simple improvement here might multiply the arena sizes by
    different factors depending upon the type of object.  There will likely be
    an approximate relation between, say, the number of edges and the number of
    walls (roughly a factor of 1.5 for triangular walls on a manifold surface).

    Another possible improvement here might be to adaptively size the memory
    arenas in the pooled allocators as we create objects.  The first arena may
    be small, but each subsequent arena would be larger.

    Obviously, there are many more complicated things we could do, but it's not
    clear that they would gain us much.

    In:  int nsubvols - how many subvolumes will share this storage
    Out: A freshly allocated storage with initialized memory pools.
 *******************************************************************/
static struct storage *create_storage(struct volume *world, int nsubvols) {
  struct storage *shared_mem = NULL;
  shared_mem =
      CHECKED_MALLOC_STRUCT(struct storage, "memory storage partition");
  memset(shared_mem, 0, sizeof(struct storage));

  if (world->mem_part_pool != 0)
    nsubvols = world->mem_part_pool;
  if (nsubvols < 8)
    nsubvols = 8;
  if (nsubvols > 4096)
    nsubvols = 4096;
  /* We should tune the algorithm for selecting allocation block sizes.  */
  /* XXX: Round up to power of 2?  Shouldn't matter, I think. */
  if ((shared_mem->list = create_mem_named(sizeof(struct wall_list), nsubvols,
                                           "wall list")) == NULL)
    mcell_allocfailed("Failed to create memory pool for wall list.");
  if ((shared_mem->mol = create_mem_named(sizeof(struct volume_molecule),
                                          nsubvols, "vol mol")) == NULL)
    mcell_allocfailed("Failed to create memory pool for volume molecules.");
  if ((shared_mem->smol = create_mem_named(sizeof(struct surface_molecule),
                                           nsubvols, "surface mol")) == NULL)
    mcell_allocfailed("Failed to create memory pool for surface molecules.");
  if ((shared_mem->face =
           create_mem_named(sizeof(struct wall), nsubvols, "wall")) == NULL)
    mcell_allocfailed("Failed to create memory pool for walls.");
  if ((shared_mem->join =
           create_mem_named(sizeof(struct edge), nsubvols, "edge")) == NULL)
    mcell_allocfailed("Failed to create memory pool for edges.");
  if ((shared_mem->grids = create_mem_named(sizeof(struct surface_grid),
                                            nsubvols, "surface grid")) == NULL)
    mcell_allocfailed("Failed to create memory pool for surface grids.");
  if ((shared_mem->regl = create_mem_named(sizeof(struct region_list), nsubvols,
                                           "region list")) == NULL)
    mcell_allocfailed("Failed to create memory pool for region lists.");
  if ((shared_mem->pslv = create_mem_named(sizeof(struct per_species_list), 32,
                                           "per species list")) == NULL)
    mcell_allocfailed(
        "Failed to create memory pool for per-species molecule lists.");
  shared_mem->coll = world->coll_mem;
  shared_mem->sp_coll = world->sp_coll_mem;
  shared_mem->tri_coll = world->tri_coll_mem;
  shared_mem->exdv = world->exdv_mem;

  if (world->chkpt_init) {
    if ((shared_mem->timer = create_scheduler(1.0, 100.0, 100, 0.0)) == NULL)
      mcell_allocfailed("Failed to create molecule scheduler.");
    shared_mem->current_time = 0.0;
  }

  if (world->time_step_max == 0.0)
    shared_mem->max_timestep = MICROSEC_PER_YEAR;
  else {
    if (world->time_step_max < world->time_unit)
      shared_mem->max_timestep = 1.0;
    else
      shared_mem->max_timestep = world->time_step_max / world->time_unit;
  }

  return shared_mem;
}

static void sanity_check_memory_subdivision(struct volume *world) {
  if (world->mem_part_x <= 0) {
    if (world->mem_part_x < 0) {
      mcell_warn("X-axis memory partition bin size set to a negative value. "
                 "Setting to default value of 14.");
      world->mem_part_x = 14;
    } else
      world->mem_part_x = 10000000;
  }
  if (world->mem_part_y <= 0) {
    if (world->mem_part_y < 0) {
      mcell_warn("Y-axis memory partition bin size set to a negative value. "
                 "Setting to default value of 14.");
      world->mem_part_y = 14;
    } else
      world->mem_part_y = 10000000;
  }
  if (world->mem_part_z <= 0) {
    if (world->mem_part_z < 0) {
      mcell_warn("Z-axis memory partition bin size set to a negative value. "
                 "Setting to default value of 14.");
      world->mem_part_z = 14;
    } else
      world->mem_part_z = 10000000;
  }
}

/********************************************************************
 init_partitions:

    Initialize the partitions for the simulation.

    When we create a storage, we need to decide how big the memory pools should
    be.

    If the memory pools are too small, you lose some of the benefits of memory
    pooling, since your memory blocks are more scattered.  You also take on
    extra overhead because of the extra arenas you end up allocating when you
    need more objects than the original pool could provide.

    If the memory pools are too large, you waste a lot of memory, as the memory
    pool will allocate far more objects than you end up using.  This can be a
    serious issue.

    In the best of all possible worlds, we'd tune the memory pools to coincide
    exactly with the peak usage of each type of object in that subdivision.
    Unfortunately, this would require uncanny prescience.

    In the absence of such foresight, this code uses a fairly dumb heuristic
    for deciding how large the memory pools should be in a subdivision.  It
    uses the number of subvolumes in the subdivision to set the size.  This
    seems, in practice, to be adequate.  It's almost certain that we can
    improve upon it.

    A fairly simple improvement here might multiply the arena sizes by
    different factors depending upon the type of object.  There will likely be
    an approximate relation between, say, the number of edges and the number of
    walls (roughly a factor of 1.5 for triangular walls on a manifold surface).

    Another possible improvement here might be to adaptively size the memory
    arenas in the pooled allocators as we create objects.  The first arena may
    be small, but each subsequent arena would be larger.

    Obviously, there are many more complicated things we could do, but it's not
    clear that they would gain us much.

    In:  int nsubvols - how many subvolumes will share this storage
    Out: A freshly allocated storage with initialized memory pools, or NULL if
         memory allocation fails.  Program state remains valid upon failure of
         this function.
 *******************************************************************/
int init_partitions(struct volume *world) {

  /* Initialize the partitions, themselves */
  if (set_partitions(world))
    return 1;

  /* Initialize dummy waypoints (why do we do this?) */
  world->n_waypoints = 1;
  world->waypoints = CHECKED_MALLOC_ARRAY(struct waypoint, world->n_waypoints,
                                          "dummy waypoint");

  /* Allocate the subvolumes */
  world->n_subvols =
      (world->nz_parts - 1) * (world->ny_parts - 1) * (world->nx_parts - 1);
  if (world->notify->progress_report != NOTIFY_NONE)
    mcell_log("Creating %d subvolumes (%d,%d,%d per axis).", world->n_subvols,
              world->nx_parts - 1, world->ny_parts - 1, world->nz_parts - 1);
  world->subvol = CHECKED_MALLOC_ARRAY(struct subvolume, world->n_subvols,
                                       "spatial subvolumes");

  /* Decide how fine-grained to make the memory subdivisions */
  sanity_check_memory_subdivision(world);

  /* Allocate the data structures which are shared between storages */
  if ((world->coll_mem = create_mem_named(sizeof(struct collision), 128,
                                          "collision")) == NULL)
    mcell_allocfailed("Failed to create memory pool for collisions.");
  if ((world->sp_coll_mem = create_mem_named(sizeof(struct sp_collision), 128,
                                             "sp collision")) == NULL)
    mcell_allocfailed(
        "Failed to create memory pool for trimolecular-pathway collisions.");
  if ((world->tri_coll_mem = create_mem_named(sizeof(struct tri_collision), 128,
                                              "tri collision")) == NULL)
    mcell_allocfailed(
        "Failed to create memory pool for trimolecular collisions.");
  if ((world->exdv_mem = create_mem_named(sizeof(struct exd_vertex), 64,
                                          "exact disk vertex")) == NULL)
    mcell_allocfailed(
        "Failed to create memory pool for exact disk calculation vertices.");

  /* How many storage subdivisions along each axis? */
  int nx = (world->nx_parts + (world->mem_part_x) - 2) / (world->mem_part_x);
  int ny = (world->ny_parts + (world->mem_part_y) - 2) / (world->mem_part_y);
  int nz = (world->nz_parts + (world->mem_part_z) - 2) / (world->mem_part_z);
  if (world->notify->progress_report != NOTIFY_NONE)
    mcell_log("Creating %d memory partitions (%d,%d,%d per axis).",
              nx * ny * nz, nx, ny, nz);

  /* Create memory pool for storages */
  if ((world->storage_allocator =
           create_mem_named(sizeof(struct storage_list), nx * ny * nz,
                            "storage allocator")) == NULL)
    mcell_allocfailed("Failed to create memory pool for storage list.");

  /* Allocate the storages */
  struct storage *shared_mem[nx * ny * nz];
  int cx = 0, cy = 0, cz = 0;
  for (int i = 0; i < nx * ny * nz; ++i) {
    /* Determine the number of subvolumes included in this subdivision */
    int xd = world->mem_part_x, yd = world->mem_part_y, zd = world->mem_part_z;
    if (cx == nx - 1)
      xd = (world->nx_parts - 1) % world->mem_part_x;
    if (cy == ny - 1)
      yd = (world->ny_parts - 1) % world->mem_part_y;
    if (cz == nz - 1)
      zd = (world->nz_parts - 1) % world->mem_part_z;
    if (++cx == nx) {
      cx = 0;
      if (++cy == ny) {
        cy = 0;
        ++cz;
      }
    }

    /* Allocate this storage */
    if ((shared_mem[i] = create_storage(world, xd * yd * zd)) == NULL)
      mcell_internal_error("Unknown error while creating a storage.");

    /* Add to the storage list */
    struct storage_list *l = (struct storage_list *)CHECKED_MEM_GET(
        world->storage_allocator, "storage list item");
    l->next = world->storage_head;
    l->store = shared_mem[i];
    world->storage_head = l;
  }

  /* Initialize each subvolume */
  for (int i = 0; i < world->nx_parts - 1; i++)
    for (int j = 0; j < world->ny_parts - 1; j++)
      for (int k = 0; k < world->nz_parts - 1; k++) {
        int h = k + (world->nz_parts - 1) * (j + (world->ny_parts - 1) * i);
        struct subvolume *sv = &(world->subvol[h]);
        sv->wall_head = NULL;
        memset(&sv->mol_by_species, 0, sizeof(struct pointer_hash));
        sv->species_head = NULL;
        sv->mol_count = 0;

        sv->llf.x = bisect_near(world->x_fineparts, world->n_fineparts,
                                world->x_partitions[i]);
        sv->llf.y = bisect_near(world->y_fineparts, world->n_fineparts,
                                world->y_partitions[j]);
        sv->llf.z = bisect_near(world->z_fineparts, world->n_fineparts,
                                world->z_partitions[k]);
        sv->urb.x = bisect_near(world->x_fineparts, world->n_fineparts,
                                world->x_partitions[i + 1]);
        sv->urb.y = bisect_near(world->y_fineparts, world->n_fineparts,
                                world->y_partitions[j + 1]);
        sv->urb.z = bisect_near(world->z_fineparts, world->n_fineparts,
                                world->z_partitions[k + 1]);

        /* Set flags so we know which directions to not go (we will fall off the
         * world!) */
        sv->world_edge =
            0; /* Assume we're not at the edge of the world in any direction */
        if (i == 0)
          sv->world_edge |= X_NEG_BIT;
        if (i == world->nx_parts - 2)
          sv->world_edge |= X_POS_BIT;
        if (j == 0)
          sv->world_edge |= Y_NEG_BIT;
        if (j == world->ny_parts - 2)
          sv->world_edge |= Y_POS_BIT;
        if (k == 0)
          sv->world_edge |= Z_NEG_BIT;
        if (k == world->nz_parts - 2)
          sv->world_edge |= Z_POS_BIT;

        /* Bind this subvolume to the appropriate storage */
        int shidx =
            (i / (world->mem_part_x)) +
            nx * (j / (world->mem_part_y) + ny * (k / (world->mem_part_z)));
        sv->local_storage = shared_mem[shidx];
      }
  return 0;
}

/**
 * Initializes the bounding boxes of the world.
 */
int init_bounding_box(struct volume *world) {
  double tm[4][4];
  double vol_infinity;

  no_printf("Initializing physical objects\n");
  vol_infinity = sqrt(DBL_MAX) / 4;
  world->bb_llf.x = vol_infinity;
  world->bb_llf.y = vol_infinity;
  world->bb_llf.z = vol_infinity;
  world->bb_urb.x = -vol_infinity;
  world->bb_urb.y = -vol_infinity;
  world->bb_urb.z = -vol_infinity;
  init_matrix(tm);

  if (compute_bb(world, world->root_instance, tm))
    return 1;

  if ((!distinguishable(world->bb_llf.x, vol_infinity, EPS_C)) &&
      (!distinguishable(world->bb_llf.y, vol_infinity, EPS_C)) &&
      (!distinguishable(world->bb_llf.z, vol_infinity, EPS_C)) &&
      (!distinguishable(world->bb_urb.x, -vol_infinity, EPS_C)) &&
      (!distinguishable(world->bb_urb.y, -vol_infinity, EPS_C)) &&
      (!distinguishable(world->bb_urb.z, -vol_infinity, EPS_C))) {
    world->bb_llf.x = 0;
    world->bb_llf.y = 0;
    world->bb_llf.z = 0;
    world->bb_urb.x = 0;
    world->bb_urb.y = 0;
    world->bb_urb.z = 0;
  }
  if (world->procnum == 0) {
    if (world->notify->progress_report) {
      mcell_log("MCell: world bounding box in microns =");
      mcell_log("         [ %.9g %.9g %.9g ] [ %.9g %.9g %.9g ]",
                world->bb_llf.x * world->length_unit,
                world->bb_llf.y * world->length_unit,
                world->bb_llf.z * world->length_unit,
                world->bb_urb.x * world->length_unit,
                world->bb_urb.y * world->length_unit,
                world->bb_urb.z * world->length_unit);
    }
  }

  world->n_walls = world->root_instance->n_walls;
  world->n_verts = world->root_instance->n_verts;
  no_printf("World object contains %d walls and %d vertices\n", world->n_walls,
            world->n_verts);

  return 0;
}

/**
 * Instantiates all physical objects.
 * This function is recursively called on the tree object objp until
 * all the objects in the data structure have been instantiated.
 * <br>
 * This function actually calls instance_release_site() and
 * instance_polygon_object() to handle the actual instantiation of
 * those objects.
 */
int instance_obj(struct volume *world, struct object *objp, double (*im)[4]) {
  double tm[4][4];
  mult_matrix(objp->t_matrix, im, tm, 4, 4, 4);

  switch (objp->object_type) {
  case META_OBJ:
    for (struct object *child_objp = objp->first_child; child_objp != NULL;
         child_objp = child_objp->next) {
      if (instance_obj(world, child_objp, tm))
        return 1;
    }
    break;

  case REL_SITE_OBJ:
    if (instance_release_site(world->magic_mem, world->releaser, objp, tm))
      return 1;
    break;

  case BOX_OBJ:
  case POLY_OBJ:
    if (instance_polygon_object(world->notify->degenerate_polys, objp))
      return 1;
    break;

  case VOXEL_OBJ:
  default:
    UNHANDLED_CASE(objp->object_type);
  }

  return 0;
}

/************************************************************************
accumulate_vertex_counts_per_storage:
    Calculates total number of vertices that belong to each storage.
    This function is recursively called on the tree object objp until
    all the vertices in the object and its children have been counted.

        In: object
        array of vertex counts per storage
            transformation matrix
        Out: 0 - on success, and 1 - on failure.
             Array of vertex counts per storage is updated for each
             object vertex and recursively for object's children
************************************************************************/
int accumulate_vertex_counts_per_storage(struct volume *world,
                                         struct object *objp,
                                         int *num_vertices_this_storage,
                                         double (*im)[4]) {
  double tm[4][4];
  mult_matrix(objp->t_matrix, im, tm, 4, 4, 4);

  switch (objp->object_type) {
  case META_OBJ:
    for (struct object *child_objp = objp->first_child; child_objp != NULL;
         child_objp = child_objp->next) {
      if (accumulate_vertex_counts_per_storage(world, child_objp,
                                               num_vertices_this_storage, tm))
        return 1;
    }
    break;

  case BOX_OBJ:
  case POLY_OBJ:
    if (accumulate_vertex_counts_per_storage_polygon_object(
            world, objp, num_vertices_this_storage, tm))
      return 1;
    break;

  case REL_SITE_OBJ:
  case VOXEL_OBJ:
  default:
    break;
  }

  return 0;
}

/*************************************************************************
accumulate_vertex_counts_per_storage_polygon_object:
        Array of vertex counts per storage is updated for each
             polygon object vertex.

    In: polygon object
        array of vertex counts per storage
            transformation matrix
        Out: 0 - on success, and 1 - on failure.
             Array of vertex counts per storage is updated for each
             polygon object vertex
**************************************************************************/
int accumulate_vertex_counts_per_storage_polygon_object(
    struct volume *world, struct object *objp, int *num_vertices_this_storage,
    double (*im)[4]) {
  struct vertex_list *vl;
  struct vector3 v;
  struct polygon_object *pop;
  /* index in the "simulated" array of storages that follows
     the linked list "world->storage_list" */
  int idx;

  pop = (struct polygon_object *)objp->contents;

  for (vl = pop->parsed_vertices; vl != NULL; vl = vl->next) {
    double p[4][4];
    p[0][0] = vl->vertex->x;
    p[0][1] = vl->vertex->y;
    p[0][2] = vl->vertex->z;
    p[0][3] = 1.0;
    mult_matrix(p, im, p, 1, 4, 4);

    v.x = p[0][0];
    v.y = p[0][1];
    v.z = p[0][2];

    idx = which_storage_contains_vertex(world, &v);
    if (idx < 0)
      return 1;

    ++num_vertices_this_storage[idx];
  }

  return 0;
}

/*************************************************************************
which_storage_contains_vertex:
    Returns the index of storage in the list of storages where
           the vertex resides (through the subvolume to which it belongs).
    In:  vertex

    Out: index of the storage in "world->storage_head" list
             or (-1) when not found
**************************************************************************/
int which_storage_contains_vertex(struct volume *world, struct vector3 *v) {
  struct subvolume *sv;
  struct storage_list *sl;
  int kk;

  sv = find_subvolume(world, v, NULL);

  for (sl = world->storage_head, kk = 0; sl != NULL; sl = sl->next, kk++) {
    if (sl->store == sv->local_storage)
      return kk;
  }

  /* if we came here, the vertex was not found in any of the storages */
  return -1;
}

/***********************************************************************
fill_world_vertices_array:
    Fills the array "world->all_vertices" with information
    about the vertices in the world by going recursively
        through all children objects.

        In: object
            array of number of vertices per storage
            transformation matrix
        Out: 0 if successful (the array "world->all_vertices" is filled),
             1 - on failure
************************************************************************/
int fill_world_vertices_array(struct volume *world, struct object *objp,
                              int *num_vertices_this_storage, double (*im)[4]) {
  double tm[4][4];
  mult_matrix(objp->t_matrix, im, tm, 4, 4, 4);

  switch (objp->object_type) {
  case META_OBJ:
    for (struct object *child_objp = objp->first_child; child_objp != NULL;
         child_objp = child_objp->next) {
      if (fill_world_vertices_array(world, child_objp,
                                    num_vertices_this_storage, tm))
        return 1;
    }
    break;

  case BOX_OBJ:
  case POLY_OBJ:
    if (fill_world_vertices_array_polygon_object(world, objp,
                                                 num_vertices_this_storage, tm))
      return 1;
    break;

  case REL_SITE_OBJ:
  case VOXEL_OBJ:
  default:
    break;
  }

  return 0;
}

/***********************************************************************
fill_world_vertices_array_polygon_object:
    Fills the array "world->all_vertices" with information
    about the vertices in the polygon object.
        Also creates and fills "objp->vertices" array

        In: object
            array of number of vertices per storage
            transformation matrix
        Out: 0 if successful (the array "world->all_vertices" is filled with
             info about vertices in the polygon object)
             1 - on failure
************************************************************************/
int fill_world_vertices_array_polygon_object(struct volume *world,
                                             struct object *objp,
                                             int *num_vertices_this_storage,
                                             double (*im)[4]) {

  struct polygon_object *pop;
  struct vertex_list *vl;
  int cur_vtx = 0; /* index */
  int which_storage, where_in_array;
  struct vector3 *v, vv;

  pop = (struct polygon_object *)objp->contents;

  objp->vertices =
      CHECKED_MALLOC_ARRAY(struct vector3 *, objp->n_verts, "polygon vertices");

  for (vl = pop->parsed_vertices; vl != NULL; vl = vl->next) {
    double p[4][4];
    p[0][0] = vl->vertex->x;
    p[0][1] = vl->vertex->y;
    p[0][2] = vl->vertex->z;
    p[0][3] = 1.0;
    mult_matrix(p, im, p, 1, 4, 4);

    vv.x = p[0][0];
    vv.y = p[0][1];
    vv.z = p[0][2];

    which_storage = which_storage_contains_vertex(world, &vv);
    where_in_array = --num_vertices_this_storage[which_storage];
    v = world->all_vertices + where_in_array;
    *v = vv;
    objp->vertices[cur_vtx++] = v;
  }

  return 0;
}

/**
 * Instantiates a release site.
 * Creates a new release site from a template release site
 * as defined in the MDL file after applying the necessary
 * geometric transformations (rotation and translation).
 * Adds the rel
 */
int instance_release_site(struct mem_helper *magic_mem,
                          struct schedule_helper *releaser, struct object *objp,
                          double (*im)[4]) {

  struct release_event_queue *reqp;
  struct release_site_obj *rsop = (struct release_site_obj *)objp->contents;

  no_printf("Instancing release site object %s\n", objp->sym->name);
  if (!distinguishable(rsop->release_prob, MAGIC_PATTERN_PROBABILITY, EPS_C)) {

    struct magic_list *ml = (struct magic_list *)CHECKED_MEM_GET(
        magic_mem, "rxn-triggered release descriptor");
    ml->data = rsop;
    ml->type = magic_release;

    struct rxn_pathname *rxpn = (struct rxn_pathname *)rsop->pattern;
    ml->next = rxpn->magic;
    rxpn->magic = ml;

    /* Region releases need to be in release queue to get initialized */
    /* Release code itself is smart enough to ignore MAGIC_PATTERNs */
    if (rsop->release_shape == SHAPE_REGION) {
      reqp = CHECKED_MALLOC_STRUCT(struct release_event_queue, "release site");
      reqp->release_site = rsop;
      reqp->event_time = 0;
      reqp->train_counter = 0;
      reqp->train_high_time = 0;
      if (schedule_add(releaser, reqp))
        mcell_allocfailed("Failed to schedule molecule release.");
    }
  } else {
    reqp = CHECKED_MALLOC_STRUCT(struct release_event_queue, "release site");
    reqp->release_site = rsop;
    reqp->event_time = rsop->pattern->delay;
    reqp->train_counter = 0;
    reqp->train_high_time = rsop->pattern->delay;
    for (int i = 0; i < 4; i++)
      for (int j = 0; j < 4; j++)
        reqp->t_matrix[i][j] = im[i][j];

    /* Schedule the release event */
    if (schedule_add(releaser, reqp))
      mcell_allocfailed("Failed to schedule molecule release.");

    if (rsop->pattern->train_duration > rsop->pattern->train_interval)
      mcell_error(
          "Release pattern train duration is greater than train interval.");
  }

  no_printf("Done instancing release site object %s\n", objp->sym->name);

  return 0;
}

/**
 * Computes the bounding box for the entire simulation world.
 * Does things recursively in a manner similar to instance_obj().
 */
static int compute_bb(struct volume *world, struct object *objp,
                      double (*im)[4]) {
  double tm[4][4];
  mult_matrix(objp->t_matrix, im, tm, 4, 4, 4);

  switch (objp->object_type) {
  case META_OBJ:
    for (struct object *child_objp = objp->first_child; child_objp != NULL;
         child_objp = child_objp->next) {
      if (compute_bb(world, child_objp, tm))
        return 1;
    }
    break;

  case REL_SITE_OBJ:
    if (compute_bb_release_site(world, objp, tm))
      return 1;
    break;

  case BOX_OBJ:
  case POLY_OBJ:
    if (compute_bb_polygon_object(world, objp, tm))
      return 1;
    break;

  case VOXEL_OBJ:
  default:
    UNHANDLED_CASE(objp->object_type);
  }

  return 0;
}

/**
 * Updates the bounding box of the world based on the size
 * and location of a release site.
 * Used by compute_bb().
 */
static int compute_bb_release_site(struct volume *world, struct object *objp,
                                   double (*im)[4]) {

  struct release_site_obj *rsop = (struct release_site_obj *)objp->contents;

  if (rsop->release_shape == SHAPE_REGION)
    return 0;

  if (rsop->location == NULL)
    mcell_error("Location is not specified for the geometrical shape release "
                "site '%s'.",
                objp->sym->name);

  double location[1][4];
  location[0][0] = rsop->location->x;
  location[0][1] = rsop->location->y;
  location[0][2] = rsop->location->z;
  location[0][3] = 1.0;
  mult_matrix(location, im, location, 1, 4, 4);

  double diam_x, diam_y, diam_z; /* diameters of the release_site */
  if (rsop->diameter == NULL) {
    diam_x = diam_y = diam_z = 0;
  } else {
    diam_x = rsop->diameter->x;
    diam_y = rsop->diameter->y;
    diam_z = rsop->diameter->z;
  }

  if (location[0][0] - diam_x < world->bb_llf.x) {
    world->bb_llf.x = location[0][0] - diam_x;
  }
  if (location[0][1] - diam_y < world->bb_llf.y) {
    world->bb_llf.y = location[0][1] - diam_y;
  }
  if (location[0][2] - diam_z < world->bb_llf.z) {
    world->bb_llf.z = location[0][2] - diam_z;
  }
  if (location[0][0] + diam_x > world->bb_urb.x) {
    world->bb_urb.x = location[0][0] + diam_x;
  }
  if (location[0][1] + diam_y > world->bb_urb.y) {
    world->bb_urb.y = location[0][1] + diam_y;
  }
  if (location[0][2] + diam_z > world->bb_urb.z) {
    world->bb_urb.z = location[0][2] + diam_z;
  }

  return 0;
}

/**
 * Updates the bounding box of the world based on the size
 * and location of a polygon_object.  Also updates the vertices in
   "pop->parsed_vertices" array.
 * Used by compute_bb().
 */
static int compute_bb_polygon_object(struct volume *world, struct object *objp,
                                     double (*im)[4]) {

  struct polygon_object *pop = (struct polygon_object *)objp->contents;

  for (struct vertex_list *vl = pop->parsed_vertices;
       vl != NULL;
       vl = vl->next) {
    double p[1][4];
    p[0][0] = vl->vertex->x;
    p[0][1] = vl->vertex->y;
    p[0][2] = vl->vertex->z;
    p[0][3] = 1.0;
    mult_matrix(p, im, p, 1, 4, 4);
    if (p[0][0] < world->bb_llf.x)
      world->bb_llf.x = p[0][0];
    if (p[0][1] < world->bb_llf.y)
      world->bb_llf.y = p[0][1];
    if (p[0][2] < world->bb_llf.z)
      world->bb_llf.z = p[0][2];
    if (p[0][0] > world->bb_urb.x)
      world->bb_urb.x = p[0][0];
    if (p[0][1] > world->bb_urb.y)
      world->bb_urb.y = p[0][1];
    if (p[0][2] > world->bb_urb.z)
      world->bb_urb.z = p[0][2];
  }

  return 0;
}

/**
 * Instantiates a polygon_object.
 * Creates walls from a template polygon_object or box object
 * as defined in the MDL file after applying the necessary geometric
 * transformations (scaling, rotation and translation).
 * <br>
 */
int instance_polygon_object(enum warn_level_t degenerate_polys,
                            struct object *objp) {

  int index_0, index_1, index_2;
  unsigned int degenerate_count;

  struct polygon_object *pop = (struct polygon_object *)objp->contents;
  int n_walls = pop->n_walls;
  double total_area = 0;

  /* Allocate and initialize walls and vertices */
  struct wall *w = CHECKED_MALLOC_ARRAY(struct wall, n_walls, "polygon walls");
  struct wall **wp = CHECKED_MALLOC_ARRAY(struct wall *, n_walls, "polygon wall pointers");

  objp->walls = w;
  objp->wall_p = wp;

  /* we do not need "parsed_vertices" info */
  if (pop->parsed_vertices != NULL) {
    free_vertex_list(pop->parsed_vertices);
    pop->parsed_vertices = NULL;
  }

  degenerate_count = 0;
  for (int n_wall = 0; n_wall < n_walls; ++n_wall) {
    if (!get_bit(pop->side_removed, n_wall)) {
      wp[n_wall] = &w[n_wall];
      index_0 = pop->element[n_wall].vertex_index[0];
      index_1 = pop->element[n_wall].vertex_index[1];
      index_2 = pop->element[n_wall].vertex_index[2];

      /* sanity check that the vertex indices are in range */
      if ((index_0 > pop->n_verts) || (index_1 > pop->n_verts) ||
          (index_2 > pop->n_verts)) {
        mcell_error("object %s has elements with out of bounds vertex indices",
                    objp->sym->name);
      }

      init_tri_wall(objp, n_wall, objp->vertices[index_0],
                    objp->vertices[index_1], objp->vertices[index_2]);
      total_area += wp[n_wall]->area;

      if (!distinguishable(wp[n_wall]->area, 0, EPS_C)) {
        if (degenerate_polys != WARN_COPE) {
          if (degenerate_polys == WARN_ERROR) {
            mcell_error("Degenerate polygon found: %s %d\n"
                        "  Vertex 0: %.5e %.5e %.5e\n"
                        "  Vertex 1: %.5e %.5e %.5e\n"
                        "  Vertex 2: %.5e %.5e %.5e",
                        objp->sym->name, n_wall, objp->vertices[index_0]->x,
                        objp->vertices[index_0]->y, objp->vertices[index_0]->z,
                        objp->vertices[index_1]->x, objp->vertices[index_1]->y,
                        objp->vertices[index_1]->z, objp->vertices[index_2]->x,
                        objp->vertices[index_2]->y, objp->vertices[index_2]->z);
          } else
            mcell_warn(
                "Degenerate polygon found and automatically removed: %s %d\n"
                "  Vertex 0: %.5e %.5e %.5e\n"
                "  Vertex 1: %.5e %.5e %.5e\n"
                "  Vertex 2: %.5e %.5e %.5e",
                objp->sym->name, n_wall, objp->vertices[index_0]->x,
                objp->vertices[index_0]->y, objp->vertices[index_0]->z,
                objp->vertices[index_1]->x, objp->vertices[index_1]->y,
                objp->vertices[index_1]->z, objp->vertices[index_2]->x,
                objp->vertices[index_2]->y, objp->vertices[index_2]->z);
        }
        set_bit(pop->side_removed, n_wall, 1);
        objp->n_walls_actual--;
        degenerate_count++;
        wp[n_wall] = NULL;
      }
    } else {
      wp[n_wall] = NULL;
    }
  }
  if (degenerate_count)
    remove_gaps_from_regions(objp);

  objp->total_area = total_area;

#ifdef DEBUG
  printf("n_walls = %d\n", n_walls);
  printf("n_walls_actual = %d\n", objp->n_walls_actual);
#endif

  return 0;
}

/********************************************************************
 init_regions:

    Traverse the world initializing regions on each object.

    In:  none
    Out: 0 on success, 1 on failure
 *******************************************************************/
static int init_regions_helper(struct volume *world) {
  if (world->clamp_list != NULL)
    init_clamp_lists(world->clamp_list);

  return instance_obj_regions(world, world->root_instance);
}

/* First part of concentration clamp initialization. */
/* After this, list is grouped by surface class. */
/* Second part (list of objects) happens with regions. */
void init_clamp_lists(struct ccn_clamp_data *clamp_list) {

  /* Sort by memory order of surface_class pointer--handy way to collect like
   * classes */
  clamp_list = (struct ccn_clamp_data *)void_list_sort(
      (struct void_list *)clamp_list);

  /* Toss other molecules in same surface class into next_mol lists */
  for (struct ccn_clamp_data *ccd = clamp_list; ccd != NULL; ccd = ccd->next) {
    while (ccd->next != NULL && ccd->surf_class == ccd->next->surf_class) {
      ccd->next->next_mol = ccd->next_mol;
      ccd->next_mol = ccd->next;
      ccd->next = ccd->next->next;
    }
    for (struct ccn_clamp_data *temp = ccd->next_mol;
         temp != NULL;
         temp = temp->next_mol) {
      temp->next = ccd->next;
    }
  }
}

/**
 * Traverse through metaobjects, placing regions on real objects as we find
 * them.
 */
int instance_obj_regions(struct volume *world, struct object *objp) {
  switch (objp->object_type) {
  case META_OBJ:
    for (struct object *child_objp = objp->first_child; child_objp != NULL;
         child_objp = child_objp->next) {
      if (instance_obj_regions(world, child_objp))
        return 1;
    }
    break;

  case REL_SITE_OBJ:
    break;

  case BOX_OBJ:
  case POLY_OBJ:
    if (init_wall_regions(world->length_unit, world->clamp_list,
                          world->species_list, world->n_species, objp))
      return 1;
    break;

  case VOXEL_OBJ:
  default:
    UNHANDLED_CASE(objp->object_type);
  }

  return 0;
}

/**
 * Initialize data associated with wall regions.
 * This function is called during wall instantiation Pass #3
 * after walls have been copied to sub-volume local memory.
 * Sets wall surf_class by region.
 * Creates surface grids.
 * Populates surface molecule tiles by region.
 * Creates virtual regions on which to clamp concentration
 */
int init_wall_regions(double length_unit, struct ccn_clamp_data *clamp_list,
                      struct species **species_list, int n_species,
                      struct object *objp) {
  struct wall *w;
  struct region *rp;
  struct region_list *rlp, *wrlp;
  struct surf_class_list *scl;
  struct edge_list *el;
  struct void_list *temp_list;
  int num_boundaries;
  struct pointer_hash *borders;
  struct edge_list *rp_borders_head;
  int surf_class_present;

  struct species *sp;
  struct name_orient *no;

  const struct polygon_object *pop = (struct polygon_object *)objp->contents;
  int n_walls = pop->n_walls;

  no_printf("Processing %d regions in polygon list object: %s\n",
            objp->num_regions, objp->sym->name);

  /* prepend a copy of sm_dat for each element referenced in each region
     of this object to the sm_prop list for the referenced element */
  for (rlp = objp->regions; rlp != NULL; rlp = rlp->next) {
    rp = rlp->reg;

    if (rp->membership == NULL)
      mcell_internal_error("Missing region information for '%s'.",
                           rp->sym->name);

    /* This code is used in the description of "restrictive regions"
       for surface molecules. The flag REGION_SET indicates that for
       the surface molecule that has CAN_REGION_BORDER flag set through
       the SURFACE_CLASS definition there are regions defined with
       this surface_class assigned. */
    if (rp->surf_class != NULL) {
      for (no = rp->surf_class->refl_mols; no != NULL; no = no->next) {
        sp = get_species_by_name(no->name, n_species, species_list);
        if (sp != NULL) {
          if ((sp->flags & REGION_PRESENT) == 0) {
            sp->flags |= REGION_PRESENT;
          }
        }
      }

      for (no = rp->surf_class->absorb_mols; no != NULL; no = no->next) {
        sp = get_species_by_name(no->name, n_species, species_list);
        if (sp != NULL) {
          if ((sp->flags & REGION_PRESENT) == 0) {
            sp->flags |= REGION_PRESENT;
          }
        }
      }

      for (no = rp->surf_class->transp_mols; no != NULL; no = no->next) {
        sp = get_species_by_name(no->name, n_species, species_list);
        if (sp != NULL) {
          if ((sp->flags & REGION_PRESENT) == 0) {
            sp->flags |= REGION_PRESENT;
          }
        }
      }
    }

    rp_borders_head = NULL;
    int count = 0;

    for (int n_wall = 0; n_wall < rp->membership->nbits; ++n_wall) {
      if (get_bit(rp->membership, n_wall)) {
        count++;
      }
    }
    if (count == n_walls)
      rp->region_has_all_elements = 1;

    for (int n_wall = 0; n_wall < rp->membership->nbits; ++n_wall) {
      if (get_bit(rp->membership, n_wall)) {
        /* prepend this region to wall region list of i_th wall only if the
         * region is used in counting */
        w = objp->wall_p[n_wall];

        rp->area += w->area;
        if (rp->surf_class != NULL) {
          /* check whether this region's surface class is already
             assigned to the wall's surface class list */
          surf_class_present = 0;
          for (scl = w->surf_class_head; scl != NULL; scl = scl->next) {
            if (scl->surf_class == rp->surf_class) {
              surf_class_present = 1;
              break;
            }
          }
          if (!surf_class_present) {
            scl = CHECKED_MALLOC_STRUCT(struct surf_class_list,
                                        "surf_class_list");
            scl->surf_class = rp->surf_class;
            if (w->surf_class_head == NULL) {
              scl->next = NULL;
              w->surf_class_head = scl;
            } else {
              scl->next = w->surf_class_head;
              w->surf_class_head = scl;
            }
            w->num_surf_classes++;
          }
        }

        if ((rp->flags & COUNT_SOME_MASK) != 0) {
          wrlp = (struct region_list *)CHECKED_MEM_GET(w->birthplace->regl,
                                                       "wall region list");
          wrlp->reg = rp;
          wrlp->next = w->counting_regions;
          w->counting_regions = wrlp;
          w->flags |= rp->flags;
        }

        /* add edges of this wall to the region's edge list */
        if ((strcmp(rp->region_last_name, "ALL") != 0) &&
            (!(rp->region_has_all_elements))) {
          for (int ii = 0; ii < 3; ii++) {
            if ((el = CHECKED_MALLOC_STRUCT(struct edge_list, "edge_list")) ==
                NULL) {
              mcell_internal_error(
                  "Out of memory while creating edge list for the region '%s'",
                  rp->sym->name);
            }
            el->ed = w->edges[ii];
            el->next = rp_borders_head;
            rp_borders_head = el;
          }
        }
      }
    } /* end for */

    /* from all edges in the region collect ones that
       constitute external borders of the region into "rp->boundaries" */

    if ((strcmp(rp->region_last_name, "ALL") != 0) &&
        (!(rp->region_has_all_elements))) {
      /* sort the linked list */
      temp_list = void_list_sort((struct void_list *)rp_borders_head);
      /* remove all internal edges */
      num_boundaries = remove_both_duplicates(&temp_list);
      rp_borders_head = (struct edge_list *)temp_list;

      if ((borders = CHECKED_MALLOC_STRUCT(struct pointer_hash,
                                           "pointer_hash")) == NULL) {
        mcell_internal_error("Out of memory while creating boundary pointer "
                             "hash for the region %s",
                             rp->sym->name);
      }

      if (pointer_hash_init(borders, 2 * num_boundaries)) {
        mcell_error(
            "Failed to initialize data structure for region boundaries.");
        /*return 1;*/
      }
      rp->boundaries = borders;

      for (el = rp_borders_head; el != NULL; el = el->next) {
        unsigned int keyhash = (unsigned int)(intptr_t)(el->ed);
        void *key = (void *)(el->ed);
        if (pointer_hash_add(rp->boundaries, key, keyhash, (void *)(el->ed))) {
          mcell_allocfailed(
              "Failed to store edge in the region pointer_hash table.");
        }
      }

      delete_void_list((struct void_list *)rp_borders_head);
      rp_borders_head = NULL;
    }

  } /*end loop over all regions in object */

  for (int n_wall = 0; n_wall < n_walls; n_wall++) {
    if (get_bit(pop->side_removed, n_wall))
      continue;

    w = objp->wall_p[n_wall];
    if (w->counting_regions != NULL) {
      w->counting_regions = (struct region_list *)void_list_sort((
          struct void_list *)w->counting_regions); /* Helpful for comparisons */
    }
    if (w->num_surf_classes > 1)
      check_for_conflicting_surface_classes(w, n_species, species_list);
  }

  /* Check to see if we need to generate virtual regions for */
  /* concentration clamps on this object */
  if (clamp_list != NULL) {
    struct ccn_clamp_data *ccd;
    struct ccn_clamp_data *temp;
    int j;
    int found_something = 0;

    for (int n_wall = 0; n_wall < n_walls; n_wall++) {
      if (get_bit(pop->side_removed, n_wall))
        continue;
      if (objp->wall_p[n_wall]->surf_class_head != NULL) {

        for (scl = objp->wall_p[n_wall]->surf_class_head; scl != NULL;
             scl = scl->next) {
          for (ccd = clamp_list; ccd != NULL; ccd = ccd->next) {
            if (scl->surf_class == ccd->surf_class) {
              if (ccd->objp != objp) {
                if (ccd->objp == NULL)
                  ccd->objp = objp;
                else if ((temp = find_clamped_object_in_list(ccd, objp)) != NULL)
                {
                  ccd = temp;
                }
                else {
                  temp = CHECKED_MALLOC_STRUCT(struct ccn_clamp_data,
                                               "concentration clamp data");
                  memcpy(temp, ccd, sizeof(struct ccn_clamp_data));
                  temp->objp = objp;
                  temp->sides = NULL;
                  temp->n_sides = 0;
                  temp->side_idx = NULL;
                  temp->cum_area = NULL;
                  ccd->next_obj = temp;
                  ccd = temp;
                }
              }
              if (ccd->sides == NULL) {
                ccd->sides = new_bit_array(n_walls);
                if (ccd->sides == NULL)
                  mcell_allocfailed("Failed to allocate membership bit array "
                                    "for concentration clamp data.");
                set_all_bits(ccd->sides, 0);
              }
              set_bit(ccd->sides, n_wall, 1);
              ccd->n_sides++;
              found_something = 1;
            }
          }
        }
      }
    }

    if (found_something) {
      for (ccd = clamp_list; ccd != NULL; ccd = ccd->next) {
        if (ccd->objp != objp) {
          if (ccd->next_obj != NULL && ccd->next_obj->objp == objp)
            ccd = ccd->next_obj;
          else
            continue;
        }

        ccd->side_idx = CHECKED_MALLOC_ARRAY(
            int, ccd->n_sides, "concentration clamp polygon side index");
        ccd->cum_area = CHECKED_MALLOC_ARRAY(
            double, ccd->n_sides,
            "concentration clamp polygon side cumulative area");

        j = 0;
        for (int n_wall = 0; n_wall < n_walls; n_wall++) {
          if (get_bit(ccd->sides, n_wall)) {
            ccd->side_idx[j] = n_wall;
            ccd->cum_area[j] = objp->wall_p[n_wall]->area;
            j++;
          }
        }
        if (j != ccd->n_sides)
          mcell_internal_error("Miscounted the number of walls for "
                               "concentration clamp.  object=%s  surface "
                               "class=%s",
                               objp->sym->name, ccd->surf_class->sym->name);

        for (j = 1; j < ccd->n_sides; j++)
          ccd->cum_area[j] += ccd->cum_area[j - 1];

        ccd->scaling_factor =
            ccd->cum_area[ccd->n_sides - 1] * length_unit *
            length_unit * length_unit /
            2.9432976599069717358e-9; /* sqrt(MY_PI)/(1e-15*N_AV) */
      }
    }
  }
  return 0;
}

/********************************************************************
 init_surf_mols:

    Traverse the world placing surface molecules.

    In:  none
    Out: 0 on success, 1 on failure
 *******************************************************************/
int init_surf_mols(struct volume *world) {
  return instance_obj_surf_mols(world, world->root_instance);
}

/********************************************************************
 instance_obj_surf_mols:

    Place any appropriate surface molecules on this object and/or its children.

    In:  struct object *objp - the object upon which to instantiate molecules
    Out: 0 on success, 1 on failure
 *******************************************************************/
int instance_obj_surf_mols(struct volume *world, struct object *objp) {
  struct object *child_objp;

  switch (objp->object_type) {
  case META_OBJ:
    for (child_objp = objp->first_child; child_objp != NULL;
         child_objp = child_objp->next) {
      if (instance_obj_surf_mols(world, child_objp))
        return 1;
    }
    break;
  case REL_SITE_OBJ:
    break;
  case BOX_OBJ:
  case POLY_OBJ:
    if (init_wall_surf_mols(world, objp))
      return 1;
    break;
  case VOXEL_OBJ:
  default:
    break;
  }

  return 0;
}

/********************************************************************
 init_wall_surf_mols:

    Place any appropriate surface molecules on this wall.  The object passed in
    must be a box or a polygon.

    In:  struct object *objp - the object upon which to instantiate molecules
    Out: 0 on success, 1 on failure
 *******************************************************************/
int init_wall_surf_mols(struct volume *world, struct object *objp) {
  struct sm_dat *smdp, *dup_smdp, **sm_prop;
  struct region_list *rlp, *rlp2, *reg_sm_num_head;
  /* byte all_region; */ /* flag that points to the region called ALL */
  struct surf_class_list *scl;

  const struct polygon_object *pop = (struct polygon_object *)objp->contents;
  int n_walls = pop->n_walls;

  /* allocate scratch storage to hold surface molecule info for each wall */
  sm_prop = CHECKED_MALLOC_ARRAY(struct sm_dat *, n_walls,
                                 "surface molecule data scratch space");

  for (int n_wall = 0; n_wall < n_walls; ++n_wall)
    sm_prop[n_wall] = NULL;

  /* prepend a copy of sm_dat for each element referenced in each region
     of this object to the sm_prop list for the referenced element */
  reg_sm_num_head = NULL;

  for (rlp = objp->regions; rlp != NULL; rlp = rlp->next) {
    struct region *rp = rlp->reg;
    byte reg_sm_num = 0;

    /* all_region = (strcmp(rp->region_last_name, "ALL") == 0); */

    /* Place molecules defined through DEFINE_SURFACE_REGIONS */
    for (int n_wall = 0; n_wall < rp->membership->nbits; n_wall++) {
      if (get_bit(rp->membership, n_wall)) {
        /* prepend region sm data for this region to sm_prop for i_th wall */
        for (smdp = rp->sm_dat_head; smdp != NULL; smdp = smdp->next) {
          if (smdp->quantity_type == SURFMOLDENS) {
            dup_smdp =
                CHECKED_MALLOC_STRUCT(struct sm_dat, "surface molecule data");
            dup_smdp->sm = smdp->sm;
            dup_smdp->quantity_type = smdp->quantity_type;
            dup_smdp->quantity = smdp->quantity;
            dup_smdp->orientation = smdp->orientation;
            dup_smdp->next = sm_prop[n_wall];
            sm_prop[n_wall] = dup_smdp;
          } else
            reg_sm_num = 1;
        }
      }
    } /* done checking each wall */

    if (rp->surf_class != NULL) {
      for (smdp = rp->surf_class->sm_dat_head; smdp != NULL;
           smdp = smdp->next) {

        if (smdp->quantity_type == SURFMOLNUM) {
          reg_sm_num = 1;
          break;
        }
      }
    }

    if (reg_sm_num) {
      rlp2 = CHECKED_MALLOC_STRUCT(struct region_list,
                                   "surface molecule placement region list");
      rlp2->reg = rp;
      rlp2->next = reg_sm_num_head;
      reg_sm_num_head = rlp2;
    }
  } /*end for (... ; rlp != NULL ; ...) */

  /* Place molecules defined through DEFINE_SURFACE_CLASSES */
  for (int n_wall = 0; n_wall < n_walls; n_wall++) {
    struct wall *w = objp->wall_p[n_wall];
    if (w == NULL)
      continue;

    for (scl = w->surf_class_head; scl != NULL; scl = scl->next) {
      for (smdp = scl->surf_class->sm_dat_head; smdp != NULL;
           smdp = smdp->next) {
        if (smdp->quantity_type == SURFMOLDENS) {
          dup_smdp =
              CHECKED_MALLOC_STRUCT(struct sm_dat, "surface molecule data");
          dup_smdp->sm = smdp->sm;
          dup_smdp->quantity_type = smdp->quantity_type;
          dup_smdp->quantity = smdp->quantity;
          dup_smdp->orientation = smdp->orientation;
          dup_smdp->next = sm_prop[n_wall];
          sm_prop[n_wall] = dup_smdp;
        }
      }
    }
  }

  /* Place regular (non-macro) molecules by density */
  for (int n_wall = 0; n_wall < n_walls; n_wall++) {
    if (!get_bit(pop->side_removed, n_wall)) {
      if (sm_prop[n_wall] != NULL) {
        if (init_surf_mols_by_density(world, objp->wall_p[n_wall],
                                      sm_prop[n_wall]))
          return 1;
      }
    }
  }

  /* Place regular (non-macro) molecules by number */
  if (reg_sm_num_head != NULL) {
    if (init_surf_mols_by_number(world, objp, reg_sm_num_head))
      return 1;

    /* free region list created to hold regions populated by number */
    rlp = reg_sm_num_head;
    while (rlp != NULL) {
      rlp2 = rlp;
      rlp = rlp->next;
      free(rlp2);
    }
  }

  /* free sm_prop array and contents */
  for (int n_wall = 0; n_wall < n_walls; n_wall++) {
    if (sm_prop[n_wall] != NULL) {
      smdp = sm_prop[n_wall];
      while (smdp != NULL) {
        dup_smdp = smdp;
        smdp = smdp->next;
        free(dup_smdp);
      }
    }
  }
  free(sm_prop);

  return 0;
}


/********************************************************************
 init_surf_mols_by_density:

    Place surface molecules on the specified wall.  This occurs after placing
    surface macromolecules, but before placing surface molecules by number.
    This is done by computing a per-tile probability, and releasing a molecule
    onto each tile with the appropriate probability.

    In:  struct wall *w - wall upon which to place
         struct sm_dat *smdp - description of what to release
    Out: 0 on success, 1 on failure
 *******************************************************************/
int init_surf_mols_by_density(struct volume *world, struct wall *w,
                              struct sm_dat *smdp_head) {

  no_printf("Initializing surface molecules by density...\n");

  if (create_grid(world, w, NULL))
    mcell_allocfailed("Failed to create grid for wall.");
  struct object *objp = w->parent_object;

  int num_sm_dat = 0;
  for (struct sm_dat *smdp = smdp_head; smdp != NULL; smdp = smdp->next)
    ++num_sm_dat;

  struct species **sm =
      CHECKED_MALLOC_ARRAY(struct species *, num_sm_dat,
                           "surface-molecule-by-density placement array");
  memset(sm, 0, num_sm_dat * sizeof(struct species *));

  double *prob = CHECKED_MALLOC_ARRAY(
      double, num_sm_dat, "surface-molecule-by-density placement array");
  memset(prob, 0, num_sm_dat * sizeof(double));

  short *orientation = CHECKED_MALLOC_ARRAY(
      short, num_sm_dat, "surface-molecule-by-density placement array");
  memset(orientation, 0, num_sm_dat * sizeof(short));

  struct surface_grid *sg = w->grid;
  unsigned int n_tiles = sg->n_tiles;
  double area = w->area;
  objp->n_tiles += n_tiles;
  no_printf("Initializing %d surf_mols...\n", n_tiles);
  no_printf("  Area = %.9g\n", area);
  no_printf("  Grid_size = %d\n", sg->n);
  no_printf("  Number of surface molecule types in wall = %d\n", num_sm_dat);

  unsigned int n_sm_entry = 0;
  double tot_prob = 0;
  double tot_density = 0;
  for (struct sm_dat *smdp = smdp_head; smdp != NULL; smdp = smdp->next) {
    no_printf("  Adding surface molecule %s to wall at density %.9g\n",
              smdp->sm->sym->name, smdp->quantity);
    tot_prob += (area * smdp->quantity) / (n_tiles * world->grid_density);
    prob[n_sm_entry] = tot_prob;
    if (smdp->orientation > 0)
      orientation[n_sm_entry] = 1;
    else if (smdp->orientation < 0)
      orientation[n_sm_entry] = -1;
    else
      orientation[n_sm_entry] = 0;
    sm[n_sm_entry++] = smdp->sm;
    tot_density += smdp->quantity;
  }

  if (tot_density > world->grid_density)
    mcell_warn(
        "Total surface molecule density too high: %f.  Filling all available "
        "surface molecule sites.",
        tot_density);

  if (world->chkpt_init) {
    for (unsigned int n_tile = 0; n_tile < n_tiles; ++n_tile) {
      if (sg->mol[n_tile] != NULL)
        continue;

      int p_index = -1;
      double rnd = rng_dbl(world->rng);
      for (int n_sm = 0; n_sm < num_sm_dat; ++n_sm) {
        if (rnd <= prob[n_sm]) {
          p_index = n_sm;
          break;
        }
      }

      if (p_index == -1)
        continue;

      short flags = TYPE_SURF | ACT_NEWBIE | IN_SCHEDULE | IN_SURFACE;
      struct surface_molecule *new_sm = place_single_molecule(
          world, w, n_tile, sm[p_index], flags, orientation[p_index], 0, 0, 0);
      if (trigger_unimolecular(world->reaction_hash, world->rx_hashsize,
                               sm[p_index]->hashval,
                               (struct abstract_molecule *)new_sm) != NULL ||
          (sm[p_index]->flags & CAN_SURFWALL) != 0) {
        new_sm->flags |= ACT_REACT;
      }
    }
  }

  unsigned int n_occupied = w->grid->n_occupied;
  sg->n_occupied = n_occupied;
  objp->n_occupied_tiles += n_occupied;

#ifdef DEBUG
  for (int n_sm = 0; n_sm < num_sm_dat; ++n_sm)
    no_printf("Total number of surface molecules %s = %d\n",
              sm[n_sm]->sym->name, sm[n_sm]->population);
#endif

  free(sm);
  free(prob);
  free(orientation);

  no_printf("Done initializing %u surface molecules by density\n", n_occupied);

  return 0;
}

/********************************************************************
 init_surf_mols_by_number:

    Place surface molecules on the specified object.  This occurs after placing
    surface macromolecules, and after placing surface molecules by density.

    In:  struct object *objp - object upon which to place
         struct region_list *reg_sm_num_head - list of what to place
    Out: 0 on success, 1 on failure
 *******************************************************************/
int init_surf_mols_by_number(struct volume *world, struct object *objp,
                             struct region_list *reg_sm_num_head) {
  static struct surface_molecule DUMMY_MOLECULE;
  static struct surface_molecule *bread_crumb = &DUMMY_MOLECULE;

  short flags = TYPE_SURF | ACT_NEWBIE | IN_SCHEDULE | IN_SURFACE;
  unsigned int n_free_sm;
  // struct subvolume *gsv = NULL;

  no_printf("Initializing surface molecules by number...\n");
  /* traverse region list and add surface molecule sites by number to whole
    regions as appropriate */
  for (struct region_list *rlp = reg_sm_num_head; rlp != NULL;
       rlp = rlp->next) {
    struct region *rp = rlp->reg;
    /* initialize surface molecule grids in region as needed and */
    /* count total number of free surface molecule sites in region */
    n_free_sm = 0;
    for (int n_wall = 0; n_wall < rp->membership->nbits; n_wall++) {
      if (get_bit(rp->membership, n_wall)) {
        struct wall *w = objp->wall_p[n_wall];
        if (create_grid(world, w, NULL))
          mcell_allocfailed("Failed to allocate grid for wall.");
        struct surface_grid *sg = w->grid;
        n_free_sm = n_free_sm + (sg->n_tiles - sg->n_occupied);
      }
    }
    no_printf("Number of free surface molecule tiles in region %s = %d\n",
              rp->sym->name, n_free_sm);

    if (n_free_sm == 0) {
      mcell_warn("Number of free surface molecule tiles in region %s = %d",
                 rp->sym->name, n_free_sm);
      continue;
    }

    if (world->chkpt_init) { /* only needed for denovo initiliazation */

      /* allocate memory to hold array of pointers to all free tiles */
      struct surface_molecule ***tiles =
          CHECKED_MALLOC_ARRAY(struct surface_molecule **, n_free_sm,
                               "surface molecule placement tiles array");

      unsigned int *idx = CHECKED_MALLOC_ARRAY(
          unsigned int, n_free_sm, "surface molecule placement indices array");

      struct wall **walls = CHECKED_MALLOC_ARRAY(
          struct wall *, n_free_sm, "surface molecule placement walls array");

      /* initialize array of pointers to all free tiles */
      int n_slot = 0;
      for (int n_wall = 0; n_wall < rp->membership->nbits; n_wall++) {
        if (get_bit(rp->membership, n_wall)) {
          struct wall *w = objp->wall_p[n_wall];
          struct surface_grid *sg = w->grid;
          if (sg != NULL) {
            for (unsigned int n_tile = 0; n_tile < sg->n_tiles; n_tile++) {
              if (sg->mol[n_tile] == NULL) {
                tiles[n_slot] = &(sg->mol[n_tile]);
                idx[n_slot] = n_tile;
                walls[n_slot++] = w;
              }
            }
          }
        }
      }

      /* distribute desired number of surface molecule sites */
      /* for each surface molecule type to add */
      /* place molecules BY NUMBER when it is defined through
       * DEFINE_SURFACE_REGION */
      for (struct sm_dat *smdp = rp->sm_dat_head; smdp != NULL;
           smdp = smdp->next) {
        if (smdp->quantity_type == SURFMOLNUM) {
          struct species *sm = smdp->sm;
          short orientation;
          unsigned int n_set = (unsigned int)smdp->quantity;
          unsigned int n_clear = n_free_sm - n_set;

          /* Compute orientation */
          if (smdp->orientation > 0)
            orientation = 1;
          else if (smdp->orientation < 0)
            orientation = -1;
          else
            orientation = 0;

          /* Clamp n_set to number of available slots (w/ warning). */
          if (n_set > n_free_sm) {
            mcell_warn(
                "Number of %s surface molecules to place (%d) exceeds number "
                "of free surface molecule tiles (%d) in region %s[%s].\n"
                "Surface molecule %s placed on all available surface molecule "
                "sites.",
                sm->sym->name, n_set, n_free_sm, rp->parent->sym->name,
                rp->region_last_name, sm->sym->name);
            n_set = n_free_sm;
            n_clear = 0;
          }

          no_printf("distribute %d of surface molecule %s\n", n_set,
                    sm->sym->name);
          no_printf("n_set = %d  n_clear = %d  n_free_sm = %d\n", n_set,
                    n_clear, n_free_sm);

          /* if filling more than half the free tiles
             init all with bread_crumbs
             choose which tiles to free again
             and then convert remaining bread_crumbs to actual molecules */
          if (n_set > n_free_sm / 2) {
            no_printf("filling more than half the free tiles: init all with "
                      "bread_crumb\n");
            for (unsigned int j = 0; j < n_free_sm; j++) {
              *tiles[j] = bread_crumb;
            }

            no_printf("choose which tiles to free again\n");
            for (unsigned int j = 0; j < n_clear; j++) {

              /* Loop until we find a vacant tile. */
              while (1) {
                int slot_num = (int)(rng_dbl(world->rng) * n_free_sm);
                if (*tiles[slot_num] == bread_crumb) {
                  *tiles[slot_num] = NULL;
                  break;
                }
              }
            }

            no_printf("convert remaining bread_crumbs to actual molecules\n");
            for (unsigned int j = 0; j < n_free_sm; j++) {
              if (*tiles[j] == bread_crumb) {
                struct surface_molecule *new_sm = place_single_molecule(
                    world, walls[j], idx[j], sm, flags, orientation, 0, 0, 0);
                if (trigger_unimolecular(
                        world->reaction_hash, world->rx_hashsize, sm->hashval,
                        (struct abstract_molecule *)new_sm) != NULL ||
                    (sm->flags & CAN_SURFWALL) != 0) {
                  new_sm->flags |= ACT_REACT;
                }
              }
            }
          } else { /* just fill only the tiles we need */
            no_printf("fill only the tiles we need\n");
            for (unsigned int j = 0; j < n_set; j++) {

              /* Loop until we find a vacant tile. */
              while (1) {
                int slot_num = (int)(rng_dbl(world->rng) * n_free_sm);
                if (*tiles[slot_num] == NULL) {
                  struct surface_molecule *new_sm = place_single_molecule(
                      world, walls[slot_num], idx[slot_num], sm, flags,
                      orientation, 0, 0, 0);
                  if (trigger_unimolecular(
                          world->reaction_hash, world->rx_hashsize, sm->hashval,
                          (struct abstract_molecule *)new_sm) != NULL ||
                      (sm->flags & CAN_SURFWALL) != 0) {
                    new_sm->flags |= ACT_REACT;
                  }
                  break;
                }
              }
            }
          }

          if (n_clear > 0) {
            struct surface_molecule ***tiles_tmp;
            unsigned int *idx_tmp;
            struct wall **walls_tmp;

            /* allocate memory to hold array of pointers to remaining free tiles
             */
            tiles_tmp =
                CHECKED_MALLOC_ARRAY(struct surface_molecule **, n_clear,
                                     "surface molecule placement tiles array");
            idx_tmp = CHECKED_MALLOC_ARRAY(
                unsigned int, n_clear,
                "surface molecule placement indices array");
            walls_tmp =
                CHECKED_MALLOC_ARRAY(struct wall *, n_clear,
                                     "surface molecule placement walls array");

            n_slot = 0;
            for (unsigned int n_sm = 0; n_sm < n_free_sm; n_sm++) {
              if (*tiles[n_sm] == NULL) {
                tiles_tmp[n_slot] = tiles[n_sm];
                idx_tmp[n_slot] = idx[n_sm];
                walls_tmp[n_slot++] = walls[n_sm];
              }
            }
            /* free original array of pointers to all free tiles */
            free(tiles);
            free(idx);
            free(walls);
            tiles = tiles_tmp;
            idx = idx_tmp;
            walls = walls_tmp;
            n_free_sm = n_free_sm - n_set;
          }

          /* update n_occupied for each surface molecule grid */
          for (int n_wall = 0; n_wall < rp->membership->nbits; n_wall++) {
            if (get_bit(rp->membership, n_wall)) {
              struct surface_grid *sg = objp->wall_p[n_wall]->grid;
              if (sg != NULL) {
                sg->n_occupied = 0;
                for (unsigned int n_tile = 0; n_tile < sg->n_tiles; ++n_tile) {
                  if (sg->mol[n_tile] != NULL)
                    sg->n_occupied++;
                }
              }
            }
          }
        }
      }

      /* place molecules BY NUMBER when it is defined through
       * DEFINE_SURFACE_CLASS */
      if (rp->surf_class != NULL) {
        for (struct sm_dat *smdp = rp->surf_class->sm_dat_head; smdp != NULL;
             smdp = smdp->next) {
          if (smdp->quantity_type == SURFMOLNUM) {
            struct species *sm = smdp->sm;
            short orientation;
            unsigned int n_set = (unsigned int)smdp->quantity;
            unsigned int n_clear = n_free_sm - n_set;

            /* Compute orientation */
            if (smdp->orientation > 0)
              orientation = 1;
            else if (smdp->orientation < 0)
              orientation = -1;
            else
              orientation = 0;

            /* Clamp n_set to number of available slots (w/ warning). */
            if (n_set > n_free_sm) {
              mcell_warn(
                  "Number of %s surface molecules to place (%d) exceeds number "
                  "of free surface molecule tiles (%d) in region %s[%s].\n"
                  "Surface molecules %s placed on all available surface "
                  "molecule sites.",
                  sm->sym->name, n_set, n_free_sm, rp->parent->sym->name,
                  rp->region_last_name, sm->sym->name);
              n_set = n_free_sm;
              n_clear = 0;
            }

            no_printf("distribute %d of surface molecule %s\n", n_set,
                      sm->sym->name);
            no_printf("n_set = %d  n_clear = %d  n_free_sm = %d\n", n_set,
                      n_clear, n_free_sm);

            /* if filling more than half the free tiles
               init all with bread_crumbs
               choose which tiles to free again
               and then convert remaining bread_crumbs to actual molecules */
            if (n_set > n_free_sm / 2) {
              no_printf("filling more than half the free tiles: init all with "
                        "bread_crumb\n");
              for (unsigned int j = 0; j < n_free_sm; j++) {
                *tiles[j] = bread_crumb;
              }

              no_printf("choose which tiles to free again\n");
              for (unsigned int j = 0; j < n_clear; j++) {

                /* Loop until we find a vacant tile. */
                while (1) {
                  int slot_num = (int)(rng_dbl(world->rng) * n_free_sm);
                  if (*tiles[slot_num] == bread_crumb) {
                    *tiles[slot_num] = NULL;
                    break;
                  }
                }
              }

              no_printf("convert remaining bread_crumbs to actual molecules\n");
              for (unsigned int j = 0; j < n_free_sm; j++) {
                if (*tiles[j] == bread_crumb) {
                  struct surface_molecule *new_sm = place_single_molecule(
                      world, walls[j], idx[j], sm, flags, orientation, 0, 0, 0);
                  if (trigger_unimolecular(
                          world->reaction_hash, world->rx_hashsize, sm->hashval,
                          (struct abstract_molecule *)new_sm) != NULL ||
                      (sm->flags & CAN_SURFWALL) != 0) {
                    new_sm->flags |= ACT_REACT;
                  }
                }
              }
            } else { /* just fill only the tiles we need */
              no_printf("fill only the tiles we need\n");
              for (unsigned int j = 0; j < n_set; j++) {

                /* Loop until we find a vacant tile. */
                while (1) {
                  int slot_num = (int)(rng_dbl(world->rng) * n_free_sm);
                  if (*tiles[slot_num] == NULL) {
                    struct surface_molecule *new_sm = place_single_molecule(
                        world, walls[slot_num], idx[slot_num], sm, flags,
                        orientation, 0, 0, 0);
                    if (trigger_unimolecular(
                            world->reaction_hash, world->rx_hashsize,
                            sm->hashval,
                            (struct abstract_molecule *)new_sm) != NULL ||
                        (sm->flags & CAN_SURFWALL) != 0) {
                      new_sm->flags |= ACT_REACT;
                    }
                    break;
                  }
                }
              }
            }

            if (n_clear > 0) {
              struct surface_molecule ***tiles_tmp;
              unsigned int *idx_tmp;
              struct wall **walls_tmp;

              /* allocate memory to hold array of pointers to remaining free
               * tiles */
              tiles_tmp = CHECKED_MALLOC_ARRAY(
                  struct surface_molecule **, n_clear,
                  "surface molecule placement tiles array");
              idx_tmp = CHECKED_MALLOC_ARRAY(
                  unsigned int, n_clear,
                  "surface molecule placement indices array");
              walls_tmp = CHECKED_MALLOC_ARRAY(
                  struct wall *, n_clear,
                  "surface molecule placement walls array");

              n_slot = 0;
              for (unsigned int n_sm = 0; n_sm < n_free_sm; n_sm++) {
                if (*tiles[n_sm] == NULL) {
                  tiles_tmp[n_slot] = tiles[n_sm];
                  idx_tmp[n_slot] = idx[n_sm];
                  walls_tmp[n_slot++] = walls[n_sm];
                }
              }
              /* free original array of pointers to all free tiles */
              free(tiles);
              free(idx);
              free(walls);
              tiles = tiles_tmp;
              idx = idx_tmp;
              walls = walls_tmp;
              n_free_sm = n_free_sm - n_set;
            }

            /* update n_occupied for each surface molecule grid */
            for (int n_wall = 0; n_wall < rp->membership->nbits; n_wall++) {
              if (get_bit(rp->membership, n_wall)) {
                struct surface_grid *sg = objp->wall_p[n_wall]->grid;
                if (sg != NULL) {
                  sg->n_occupied = 0;
                  for (unsigned int n_tile = 0; n_tile < sg->n_tiles;
                       ++n_tile) {
                    if (sg->mol[n_tile] != NULL)
                      sg->n_occupied++;
                  }
                }
              }
            }
          }
        }
      } /* end of if (rp->surf_clas != NULL) */

      /* free array of pointers to all free tiles */
      if (tiles != NULL) {
        free(tiles);
      }
      if (idx != NULL) {
        free(idx);
      }
      if (walls != NULL) {
        free(walls);
      }
    } /* end if (world->chkpt_init) */
  }
  no_printf("Done initializing surface molecules by number.\n");
  return 0;
}

/***************************************************************************
rel_expr_grab_obj:
  In: release expression
      place to allocate memory for temporary void_list
  Out: a linked list containing all the objects referred to in the
       release expression (including duplcates), or NULL if there are
       no such objects.
***************************************************************************/

/* Not the most efficient due to slow merging, but it works. */
static struct void_list *rel_expr_grab_obj(struct release_evaluator *root,
                                           struct mem_helper *voidmem) {
  struct void_list *vl = NULL;
  struct void_list *vr = NULL;

  if (root->left != NULL) {
    if (root->op & REXP_LEFT_REGION) {
      vl = CHECKED_MEM_GET(voidmem, "temporary list for region release");
      if (vl == NULL)
        return NULL;
      vl->data = ((struct region *)(root->left))->parent;
      vl->next = NULL;
    } else
      vl = rel_expr_grab_obj(root->left, voidmem);
  }
  if (root->right != NULL) {
    if (root->op & REXP_RIGHT_REGION) {
      vr = CHECKED_MEM_GET(voidmem, "temporary list for region release");
      if (vr == NULL)
        return NULL;
      vr->data = ((struct region *)(root->right))->parent;
      vr->next = NULL;
    } else
      vr = rel_expr_grab_obj(root->right, voidmem);
  }

  if (vl == NULL) {
    if (vr == NULL)
      return NULL;
    return vr;
  } else if (vr == NULL) {
    return vl;
  } else {
    struct void_list *vp;

    for (vp = vl; vp->next != NULL; vp = vp->next) {
    }

    vp->next = vr;

    return vl;
  }
  return NULL;
}

/***************************************************************************
find_unique_rev_objects:
  In: release expression
      place to store the number of unique objects we find
  Out: an array of pointers to each object listed in the release
       expression (no duplicates), or NULL if out of memory.  The
       second argument is set to the length of the array.
***************************************************************************/
static struct object **find_unique_rev_objects(struct release_evaluator *root,
                                               int *n) {
  struct object **o_array;
  struct void_list *vp, *vq;
  struct mem_helper *voidmem;
  int n_unique;

  voidmem = create_mem(sizeof(struct void_list), 1024);
  if (voidmem == NULL)
    mcell_allocfailed("Failed to create temporary list memory pool.");

  vp = rel_expr_grab_obj(root, voidmem);
  if (vp == NULL)
    return NULL;

  vp = void_list_sort(vp);

  for (n_unique = 1, vq = vp; vq != NULL && vq->next != NULL;
       vq = vq->next, n_unique++) {
    while (vq->data == vq->next->data) {
      vq->next = vq->next->next;
      if (vq->next == NULL)
        break;
    }
  }

  if (vq == NULL)
    n_unique--;
  *n = n_unique;

  o_array = CHECKED_MALLOC_ARRAY(struct object *, n_unique,
                                 "object array for region release");
  vq = vp;
  for (unsigned int n_obj = 0; vq != NULL; vq = vq->next, ++n_obj)
    o_array[n_obj] = (struct object *)vq->data;

  delete_mem(voidmem);

  return o_array;
}

/***************************************************************************
eval_rel_region_expr:
  In: release expression for a 2D region release
      the number of distinct objects in the world listed in the expression
      array of pointers to each of those objects
      array of pointers to bit arrays specifying which walls of each
        object are included in this release
  Out: 0 on success, 1 on failure.  On success, the bit arrays are set
       so that they indicate which walls of each object are included in
       this release site.
***************************************************************************/
static int eval_rel_region_expr(struct release_evaluator *expr, int n,
                                struct object **objs,
                                struct bit_array **result) {
  char bit_op;

  if (expr->left != NULL) {
    if (expr->op & REXP_LEFT_REGION) {
      int pos = void_array_search((void **)objs, n,
                                  ((struct region *)(expr->left))->parent);
      result[pos] =
          duplicate_bit_array(((struct region *)(expr->left))->membership);
      if (result[pos] == NULL)
        return 1;
    } else {
      if (eval_rel_region_expr(expr->left, n, objs, result))
        return 1;
    }

    if (expr->right == NULL) {
      if (expr->op & REXP_NO_OP)
        return 0;
      else
        return 1;
    }

    if (expr->op & REXP_RIGHT_REGION) {
      int pos = void_array_search((void **)objs, n,
                                  ((struct region *)(expr->right))->parent);
      if (result[pos] == NULL) {
        result[pos] =
            duplicate_bit_array(((struct region *)(expr->right))->membership);
        if (result[pos] == NULL)
          return 1;
      } else {
        if (expr->op & REXP_UNION)
          bit_op = '|';
        else if (expr->op & REXP_SUBTRACTION)
          bit_op = '-';
        else if (expr->op & REXP_INTERSECTION)
          bit_op = '&';
        else
          return 1;

        bit_operation(result[pos], ((struct region *)(expr->right))->membership,
                      bit_op);
      }
    } else {
      struct bit_array *res2[n];
      for (int i = 0; i < n; i++)
        res2[i] = NULL;

      if (eval_rel_region_expr(expr->right, n, objs, res2))
        return 1;

      for (int i = 0; i < n; i++) {
        if (res2[i] == NULL)
          continue;
        if (result[i] == NULL)
          result[i] = res2[i];
        else {
          if (expr->op & REXP_UNION)
            bit_op = '|';
          else if (expr->op & REXP_SUBTRACTION)
            bit_op = '-';
          else if (expr->op & REXP_INTERSECTION)
            bit_op = '&';
          else
            return 1;

          bit_operation(result[i], res2[i], bit_op);
          free_bit_array(res2[i]);
        }
      }
    }
  } else
    return 1; /* Left should always have something! */

  return 0;
}

/***************************************************************************
init_rel_region_data_2d:
  In: release data for a release of 2D molecules onto a region
  Out: 0 on success, 1 on failure.  A summary of all potentially available
       space over all objects contained in the region expression is
       generated and stored in arrays (typically of length equal to the
       number of walls in the region expression).
***************************************************************************/
static int init_rel_region_data_2d(struct release_site_obj *rsop,
                                   struct release_region_data *rrd) {
  rrd->owners = find_unique_rev_objects(rrd->expression, &(rrd->n_objects));
  if (rrd->owners == NULL)
    mcell_error("No objects were found matching the 2-D region release request "
                "for release site '%s'.",
                rsop->name);

  rrd->in_release =
      CHECKED_MALLOC_ARRAY(struct bit_array *, rrd->n_objects,
                           "region membership array for 2D region release");
  for (int n_object = 0; n_object < rrd->n_objects; ++n_object)
    rrd->in_release[n_object] = NULL;

  if (eval_rel_region_expr(rrd->expression, rrd->n_objects, rrd->owners,
                           rrd->in_release))
    mcell_error("Could not evaluate region expression for release site '%s'.",
                rsop->name);

  for (int n_object = 0; n_object < rrd->n_objects; n_object++) {
    if (rrd->owners[n_object] == NULL)
      mcell_internal_error("Object %d of %d in region expression for release "
                           "site '%s' was not found!",
                           n_object + 1, rrd->n_objects, rsop->name);
  }

  rrd->walls_per_obj = CHECKED_MALLOC_ARRAY(
      int, rrd->n_objects, "wall counts for 2D region release");

  rrd->n_walls_included = 0;
  for (int n_object = 0; n_object < rrd->n_objects; ++n_object) {
    if (rrd->in_release[n_object] == NULL)
      rrd->walls_per_obj[n_object] = 0;
    else
      rrd->walls_per_obj[n_object] = count_bits(rrd->in_release[n_object]);
    rrd->n_walls_included += rrd->walls_per_obj[n_object];
  }

  rrd->cum_area_list =
      CHECKED_MALLOC_ARRAY(double, rrd->n_walls_included,
                           "cumulative area list for 2D region release");
  rrd->wall_index = CHECKED_MALLOC_ARRAY(int, rrd->n_walls_included,
                                         "wall indices for 2D region release");
  rrd->obj_index = CHECKED_MALLOC_ARRAY(int, rrd->n_walls_included,
                                        "object indices for 2D region release");

  unsigned int n_wall_overall = 0;
  for (int n_object = 0; n_object < rrd->n_objects; ++n_object) {
    if (rrd->walls_per_obj[n_object] == 0)
      continue;
    int owner_type = rrd->owners[n_object]->object_type;
    if (owner_type != POLY_OBJ && owner_type != BOX_OBJ)
      mcell_internal_error("Found a region on an object which is neither a box "
                           "nor a polygon (type=%d).",
                           owner_type);

    struct polygon_object *po =
        (struct polygon_object *)(rrd->owners[n_object]->contents);
    int n_walls = po->n_walls;
    for (int n_wall = 0; n_wall < n_walls; ++n_wall) {
      if (get_bit(rrd->in_release[n_object], n_wall)) {
        rrd->cum_area_list[n_wall_overall] =
            rrd->owners[n_object]->wall_p[n_wall]->area;
        rrd->obj_index[n_wall_overall] = n_object;
        rrd->wall_index[n_wall_overall] = n_wall;
        ++n_wall_overall;
      }
    }
  }

  for (int n_wall = 1; n_wall < rrd->n_walls_included; n_wall++) {
    rrd->cum_area_list[n_wall] += rrd->cum_area_list[n_wall - 1];
  }

  return 0;
}

/***************************************************************************
create_region_bbox:
  In: a region
  Out: pointer to a 2-element array contining the LLF and URB corners of
       a bounding box around the region, or NULL if out of memory.
***************************************************************************/
struct vector3 *create_region_bbox(struct region *r) {
  struct vector3 *bbox =
      CHECKED_MALLOC_ARRAY(struct vector3, 2, "region bounding box");

  int found_first_wall = 0;
  for (int n_wall = 0; n_wall < r->membership->nbits; ++n_wall) {
    if (get_bit(r->membership, n_wall)) {
      if (!found_first_wall) {
        bbox[0].x = bbox[1].x = r->parent->wall_p[n_wall]->vert[0]->x;
        bbox[0].y = bbox[1].y = r->parent->wall_p[n_wall]->vert[0]->y;
        bbox[0].z = bbox[1].z = r->parent->wall_p[n_wall]->vert[0]->z;
        found_first_wall = 1;
      }
      for (unsigned int n_vert = 0; n_vert < 3; ++n_vert) {
        struct vector3 *v = r->parent->wall_p[n_wall]->vert[n_vert];
        if (bbox[0].x > v->x)
          bbox[0].x = v->x;
        else if (bbox[1].x < v->x)
          bbox[1].x = v->x;
        if (bbox[0].y > v->y)
          bbox[0].y = v->y;
        else if (bbox[1].y < v->y)
          bbox[1].y = v->y;
        if (bbox[0].z > v->z)
          bbox[0].z = v->z;
        else if (bbox[1].z < v->z)
          bbox[1].z = v->z;
      }
    }
  }

  return bbox;
}

/***************************************************************************
eval_rel_region_bbox:
  In: release expression for a 3D region release
      place to store LLF corner of the bounding box for the release
      place to store URB corner
  Out: 0 on success, 1 on failure.  Bounding box is set based on release
       expression (based boolean intersection of bounding boxes for each
       region).  The function reports failure if any region is unclosed.
***************************************************************************/
static int eval_rel_region_bbox(struct release_evaluator *expr,
                                struct vector3 *llf, struct vector3 *urb) {
  struct region *r;

  if (expr->left != NULL) {
    if (expr->op & REXP_LEFT_REGION) {
      r = (struct region *)(expr->left);

      if (r->bbox == NULL)
        r->bbox = create_region_bbox(r);

      llf->x = r->bbox[0].x;
      llf->y = r->bbox[0].y;
      llf->z = r->bbox[0].z;
      urb->x = r->bbox[1].x;
      urb->y = r->bbox[1].y;
      urb->z = r->bbox[1].z;

      if (r->manifold_flag == MANIFOLD_UNCHECKED) {
        if (is_manifold(r))
          r->manifold_flag = IS_MANIFOLD;
        else
          mcell_error(
              "Cannot release a 3D molecule inside the unclosed region '%s'.",
              r->sym->name);
      }

    } else {
      if (eval_rel_region_bbox(expr->left, llf, urb))
        return 1;
    }

    if (expr->right == NULL) {
      if (expr->op & REXP_NO_OP)
        return 0;
      else
        mcell_internal_error(
            "Right subtree of release expression is unexpectedly NULL.");
    }

    if (expr->op & REXP_SUBTRACTION)
      return 0;
    else {
      struct vector3 llf2;
      struct vector3 urb2;

      if (expr->op & REXP_RIGHT_REGION) {
        r = (struct region *)(expr->right);
        if (r->manifold_flag == MANIFOLD_UNCHECKED) {
          if (is_manifold(r))
            r->manifold_flag = IS_MANIFOLD;
          else
            mcell_error(
                "Cannot release a 3D molecule inside the unclosed region '%s'.",
                r->sym->name);
        }

        if (r->bbox == NULL)
          r->bbox = create_region_bbox(r);

        llf2.x = r->bbox[0].x;
        llf2.y = r->bbox[0].y;
        llf2.z = r->bbox[0].z;
        urb2.x = r->bbox[1].x;
        urb2.y = r->bbox[1].y;
        urb2.z = r->bbox[1].z;
      } else {
        if (eval_rel_region_bbox(expr->right, &llf2, &urb2))
          return 1;
      }

      if (expr->op & REXP_UNION) {
        if (llf->x > llf2.x)
          llf->x = llf2.x;
        if (llf->y > llf2.y)
          llf->y = llf2.y;
        if (llf->z > llf2.z)
          llf->z = llf2.z;
        if (urb->x < urb2.x)
          urb->x = urb2.x;
        if (urb->y < urb2.y)
          urb->y = urb2.y;
        if (urb->z < urb2.z)
          urb->z = urb2.z;
      } else if (expr->op & (REXP_INTERSECTION)) {
        if (llf->x < llf2.x)
          llf->x = llf2.x;
        if (llf->y < llf2.y)
          llf->y = llf2.y;
        if (llf->z < llf2.z)
          llf->z = llf2.z;
        if (urb->x > urb2.x)
          urb->x = urb2.x;
        if (urb->y > urb2.y)
          urb->y = urb2.y;
        if (urb->z > urb2.z)
          urb->z = urb2.z;
      } else
        mcell_internal_error("Release expression contains an unknown or "
                             "unexpected operator: (%d).",
                             expr->op);
    }
  } else
    mcell_internal_error(
        "Left subtree of release expression is unexpectedly NULL.");

  return 0;
}

/***************************************************************************
init_rel_region_data_3d:
  In: release region data structure
  Out: 0 on success, 1 on failure, -1 if there is no volume contained
       in the release (user error).  The release must be for a 3D release
       of molecules.  eval_rel_region_bbox is called to perform the
       initialization.
***************************************************************************/
static int init_rel_region_data_3d(struct release_region_data *rrd) {
  rrd->n_walls_included = 0;

  if (eval_rel_region_bbox(rrd->expression, &(rrd->llf), &(rrd->urb)))
    return 1;

  if (rrd->llf.x >= rrd->urb.x || rrd->llf.y >= rrd->urb.y ||
      rrd->llf.z >= rrd->urb.z) {
    return -1; /* Special signal to print out "nothing in here" error msg */
  }

  return 0;
}

/***************************************************************************
output_regrel_eval_tree:
  In: file to print tree on
      prefix string to put before the current branch of the tree
      prefix character for left half of current branch
      prefix character for right half of current branch
      release expression
  Out: no return value.  The tree is printed to the file.
***************************************************************************/
static void output_relreg_eval_tree(FILE *f, char *prefix, char cA, char cB,
                                    struct release_evaluator *expr) {
  size_t l = strlen(prefix);
  char my_op;

  if (expr->op & REXP_NO_OP) {
    fprintf(f, "%s >%s\n", prefix, ((struct region *)(expr->left))->sym->name);
  } else {
    char prefixA[l + 3];
    char prefixB[l + 3];
    strncpy(prefixA, prefix, l);
    strncpy(prefixB, prefix, l);
    prefixA[l] = cA;
    prefixB[l] = cB;
    prefixA[l + 1] = prefixB[l + 1] = ' ';
    prefixA[l + 2] = prefixB[l + 2] = 0;

    if (expr->op & REXP_LEFT_REGION) {
      fprintf(f, "%s >%s\n", prefix,
              ((struct region *)(expr->left))->sym->name);
    } else {
      output_relreg_eval_tree(f, prefixA, ' ', '|', expr->left);
    }

    my_op = '?';
    if (expr->op & REXP_UNION)
      my_op = '+';
    else if (expr->op & REXP_INTERSECTION)
      my_op = '*';
    else if (expr->op & REXP_SUBTRACTION)
      my_op = '-';

    fprintf(f, "%s%c\n", prefix, my_op);

    if (expr->op & REXP_RIGHT_REGION) {
      fprintf(f, "%s >%s\n", prefix,
              ((struct region *)(expr->left))->sym->name);
    } else {
      output_relreg_eval_tree(f, prefixA, '|', ' ', expr->right);
    }
  }
}

/***************************************************************************
init_releases:
  In: nothing
  Out: 0 on success, 1 on failure.  All release sites are initialized.
       Right now, the only release sites that need to be initialized are
       releases on regions.
***************************************************************************/
int init_releases(struct schedule_helper *releaser) {
  struct release_event_queue *req;
  struct abstract_element *ae;
  struct schedule_helper *sh;
  int i;

  for (sh = releaser; sh != NULL; sh = sh->next_scale) {
    for (i = -1; i < sh->buf_len; i++) {
      for (ae = (i == -1) ? sh->current : sh->circ_buf_head[i]; ae != NULL;
           ae = ae->next) {
        req = (struct release_event_queue *)ae;
        switch ((int)req->release_site->release_shape) {
        case SHAPE_REGION:
          if (req->release_site->mol_type == NULL)
            mcell_error("Molecule type was not specified for the region "
                        "release site '%s'.",
                        req->release_site->name);
          if ((req->release_site->mol_type->flags & NOT_FREE) == 0) {
            switch (init_rel_region_data_3d(req->release_site->region_data)) {
            case 0:
              break;

            case -1:
              mcell_warn("Region release site '%s' is empty!  Ignoring!  "
                         "Evaluation tree:\n",
                         req->release_site->name);
              output_relreg_eval_tree(
                  mcell_get_error_file(), " ", ' ', ' ',
                  req->release_site->region_data->expression);
              req->release_site->release_number_method = CONSTNUM;
              req->release_site->release_number = 0;
              break;

            default:
              mcell_error("Unexpected error while initializing 3-D region "
                          "releases for release site '%s'.",
                          req->release_site->name);
              /*break;*/
            }
          } else {
            if (init_rel_region_data_2d(req->release_site,
                                        req->release_site->region_data))
              mcell_error("Unexpected error while initializing 2-D region "
                          "releases for release site '%s'.",
                          req->release_site->name);
          }
          break;

        case SHAPE_LIST:
          if (req->release_site->mol_list == NULL)
            mcell_error("Molecule positions for the LIST release site '%s' are "
                        "not specified.",
                        req->release_site->name);
          break;

        case SHAPE_SPHERICAL:
        case SHAPE_CUBIC:
        case SHAPE_ELLIPTIC:
        case SHAPE_RECTANGULAR:
        case SHAPE_SPHERICAL_SHELL:
          /* geometrical release sites */
          if (req->release_site->mol_type == NULL)
            mcell_error(
                "Molecule type for the release site '%s' is not specified.",
                req->release_site->name);
          if (req->release_site->diameter == NULL)
            mcell_error("Diameter for the geometrical shape release site '%s' "
                        "is not specified.",
                        req->release_site->name);
          break;

        case SHAPE_UNDEFINED:
        default:
          UNHANDLED_CASE(req->release_site->release_shape);
        }
      }
    }
  }

  return 0;
}

/***************************************************************************
publish_special_reactions_report:
  In: species that is a SURFACE_CLASS
      lists of names of volume and surface molecules
  Out: None.  If species is a surface (surface class) and special reactions
       like TRANSPARENT, REFLECTIVE or ABSORPTIVE are defined for it,
       the reactions report is printed out.
***************************************************************************/
void publish_special_reactions_report(struct species *sp,
                                      struct name_list *vol_species_name_list,
                                      struct name_list *surf_species_name_list,
                                      int n_species,
                                      struct species **species_list) {
  struct name_orient *no;
  FILE *log_file;
  struct species *spec;
  /* orientation of ALL_MOLECULES, ALL_VOLUME_MOLECULES, ALL_SURFACE_MOLECULES
   */
  int all_mol_orient = INT_MIN;
  int all_volume_mol_orient = INT_MIN;
  int all_surface_mol_orient = INT_MIN;
  /* flags */
  int refl_mols_all_volume_mol = 0;
  int refl_mols_all_surface_mol = 0;
  int refl_mols_all_mol = 0;
  int transp_mols_all_volume_mol = 0;
  int transp_mols_all_surface_mol = 0;
  int transp_mols_all_mol = 0;
  int absorb_mols_all_volume_mol = 0;
  int absorb_mols_all_surface_mol = 0;
  int absorb_mols_all_mol = 0;

  struct name_list *nl;
  int surf_refl_title_printed;      /*flag*/
  int borders_refl_title_printed;   /*flag*/
  int surf_transp_title_printed;    /*flag*/
  int borders_transp_title_printed; /*flag*/
  int surf_absorb_title_printed;    /*flag*/
  int borders_absorb_title_printed; /*flag*/

  /* Below I will employ following set of rules for printing out
     the relative orientations of surface classes and molecules;
     1. The orientation of surface class is always printed out as {1}.
     2. The orientation of molecule is printed out
        as {1} when it is positive, {-1} when it is negative,
        and {0} when it is zero, or absent.
  */

  log_file = mcell_get_log_file();

  if (sp->refl_mols != NULL) {
    /* search for ALL_MOLECULES, ALL_VOLUME_MOLECULES, ALL_SURFACE_MOLECULES */
    for (no = sp->refl_mols; no != NULL; no = no->next) {
      if (strcmp(no->name, "ALL_VOLUME_MOLECULES") == 0) {
        all_volume_mol_orient = no->orient;
        refl_mols_all_volume_mol = 1;
      }
      if (strcmp(no->name, "ALL_SURFACE_MOLECULES") == 0) {
        all_surface_mol_orient = no->orient;
        refl_mols_all_surface_mol = 1;
      }
      if (strcmp(no->name, "ALL_MOLECULES") == 0) {
        all_mol_orient = no->orient;
        refl_mols_all_mol = 1;
      }
    }

    if ((refl_mols_all_volume_mol || refl_mols_all_mol) &&
        (vol_species_name_list != NULL)) {
      fprintf(log_file, "Surfaces with surface class \"%s{1}\" are REFLECTIVE "
                        "for volume molecules  ",
              sp->sym->name);
      for (nl = vol_species_name_list; nl != NULL; nl = nl->next) {
        if (refl_mols_all_volume_mol) {
          fprintf(log_file, "%s{%d}", nl->name, all_volume_mol_orient);
        } else if (refl_mols_all_mol) {
          fprintf(log_file, "%s{%d}", nl->name, all_mol_orient);
        }
        if (nl->next != NULL)
          fprintf(log_file, " ");
        else
          fprintf(log_file, ".");
      }
    }
    fprintf(log_file, "\n");

    if ((refl_mols_all_surface_mol || refl_mols_all_mol) &&
        (surf_species_name_list != NULL)) {
      fprintf(log_file, "Borders of regions with surface class \"%s{1}\" are "
                        "REFLECTIVE for surface molecules  ",
              sp->sym->name);
      for (nl = surf_species_name_list; nl != NULL; nl = nl->next) {
        if (refl_mols_all_surface_mol) {
          fprintf(log_file, "%s{%d}", nl->name, all_surface_mol_orient);
        } else if (refl_mols_all_mol) {
          fprintf(log_file, "%s{%d}", nl->name, all_mol_orient);
        }
        if (nl->next != NULL)
          fprintf(log_file, " ");
        else
          fprintf(log_file, ".");
      }
    }
    fprintf(log_file, "\n");

    surf_refl_title_printed = 0;
    borders_refl_title_printed = 0;
    /* First go through all volume molecules */
    if ((!refl_mols_all_mol) && (!refl_mols_all_volume_mol)) {
      for (no = sp->refl_mols; no != NULL; no = no->next) {
        spec = get_species_by_name(no->name, n_species, species_list);
        if (spec == NULL)
          mcell_internal_error("Cannot find molecule name %s", no->name);
        if (spec->flags & ON_GRID)
          continue;
        if (strcmp(no->name, "ALL_MOLECULES") == 0)
          continue;
        if (strcmp(no->name, "ALL_VOLUME_MOLECULES") == 0)
          continue;
        if (strcmp(no->name, "ALL_SURFACE_MOLECULES") == 0)
          continue;

        if (!surf_refl_title_printed) {
          surf_refl_title_printed = 1;
          fprintf(log_file, "Surfaces with surface class \"%s{1}\" are "
                            "REFLECTIVE for volume molecules  ",
                  sp->sym->name);
        }
        fprintf(log_file, "%s{%d}", no->name, no->orient);
        if (no->next != NULL)
          fprintf(log_file, " ");
      }
      if (surf_refl_title_printed)
        fprintf(log_file, ".\n");
    }

    /* Now go through all surface molecules */
    if ((!refl_mols_all_mol) && (!refl_mols_all_surface_mol)) {
      for (no = sp->refl_mols; no != NULL; no = no->next) {
        spec = get_species_by_name(no->name, n_species, species_list);
        if (spec == NULL)
          mcell_internal_error("Cannot find molecule name %s", no->name);
        if ((spec->flags & NOT_FREE) == 0)
          continue;
        if (strcmp(no->name, "ALL_MOLECULES") == 0)
          continue;
        if (strcmp(no->name, "ALL_VOLUME_MOLECULES") == 0)
          continue;
        if (strcmp(no->name, "ALL_SURFACE_MOLECULES") == 0)
          continue;
        if (!borders_refl_title_printed) {
          borders_refl_title_printed = 1;
          fprintf(log_file, "Borders of regions with surface class \"%s{1}\" "
                            "are REFLECTIVE for surface molecules  ",
                  sp->sym->name);
        }
        fprintf(log_file, "%s{%d}", no->name, no->orient);
        if (no->next != NULL)
          fprintf(log_file, " ");
      }
      if (borders_refl_title_printed)
        fprintf(log_file, ".\n");
    }
  }

  if (sp->transp_mols != NULL) {
    /* search for ALL_MOLECULES, ALL_VOLUME_MOLECULES, ALL_SURFACE_MOLECULES */
    for (no = sp->transp_mols; no != NULL; no = no->next) {
      if (strcmp(no->name, "ALL_VOLUME_MOLECULES") == 0) {
        all_volume_mol_orient = no->orient;
        transp_mols_all_volume_mol = 1;
      }
      if (strcmp(no->name, "ALL_SURFACE_MOLECULES") == 0) {
        all_surface_mol_orient = no->orient;
        transp_mols_all_surface_mol = 1;
      }
      if (strcmp(no->name, "ALL_MOLECULES") == 0) {
        all_mol_orient = no->orient;
        transp_mols_all_mol = 1;
      }
    }

    if ((transp_mols_all_volume_mol || transp_mols_all_mol) &&
        (vol_species_name_list != NULL)) {
      fprintf(log_file, "Surfaces with surface class \"%s{1}\" are TRANSPARENT "
                        "for volume molecules  ",
              sp->sym->name);
      for (nl = vol_species_name_list; nl != NULL; nl = nl->next) {
        if (transp_mols_all_volume_mol) {
          fprintf(log_file, "%s{%d}", nl->name, all_volume_mol_orient);
        } else if (transp_mols_all_mol) {
          fprintf(log_file, "%s{%d}", nl->name, all_mol_orient);
        }
        if (nl->next != NULL)
          fprintf(log_file, " ");
        else
          fprintf(log_file, ".");
      }
    }
    fprintf(log_file, "\n");

    if ((transp_mols_all_surface_mol || transp_mols_all_mol) &&
        (surf_species_name_list != NULL)) {
      fprintf(log_file, "Borders of regions with surface class \"%s{1}\" are "
                        "TRANSPARENT for surface molecules  ",
              sp->sym->name);
      for (nl = surf_species_name_list; nl != NULL; nl = nl->next) {
        if (transp_mols_all_surface_mol) {
          fprintf(log_file, "%s{%d}", nl->name, all_surface_mol_orient);
        } else if (transp_mols_all_mol) {
          fprintf(log_file, "%s{%d}", nl->name, all_mol_orient);
        }
        if (nl->next != NULL)
          fprintf(log_file, " ");
        else
          fprintf(log_file, ".");
      }
    }
    fprintf(log_file, "\n");

    surf_transp_title_printed = 0;
    borders_transp_title_printed = 0;
    /* First go through all volume molecules */
    if ((!transp_mols_all_mol) && (!transp_mols_all_volume_mol)) {
      for (no = sp->transp_mols; no != NULL; no = no->next) {
        spec = get_species_by_name(no->name, n_species, species_list);
        if (spec == NULL)
          mcell_internal_error("Cannot find molecule name %s", no->name);
        if (spec->flags & ON_GRID)
          continue;
        if (strcmp(no->name, "ALL_MOLECULES") == 0)
          continue;
        if (strcmp(no->name, "ALL_VOLUME_MOLECULES") == 0)
          continue;
        if (strcmp(no->name, "ALL_SURFACE_MOLECULES") == 0)
          continue;

        if (!surf_transp_title_printed) {
          surf_transp_title_printed = 1;
          fprintf(log_file, "Surfaces with surface class \"%s{1}\" are "
                            "TRANSPARENT for volume molecules  ",
                  sp->sym->name);
        }
        fprintf(log_file, "%s{%d}", no->name, no->orient);
        if (no->next != NULL)
          fprintf(log_file, " ");
      }
      if (surf_transp_title_printed)
        fprintf(log_file, ".\n");
    }

    /* Now go through all surface molecules */
    if ((!transp_mols_all_mol) && (!transp_mols_all_surface_mol)) {
      for (no = sp->transp_mols; no != NULL; no = no->next) {
        spec = get_species_by_name(no->name, n_species, species_list);
        if (spec == NULL)
          mcell_internal_error("Cannot find molecule name %s", no->name);
        if ((spec->flags & NOT_FREE) == 0)
          continue;
        if (strcmp(no->name, "ALL_MOLECULES") == 0)
          continue;
        if (strcmp(no->name, "ALL_VOLUME_MOLECULES") == 0)
          continue;
        if (strcmp(no->name, "ALL_SURFACE_MOLECULES") == 0)
          continue;

        if (!borders_transp_title_printed) {
          borders_transp_title_printed = 1;
          fprintf(log_file, "Borders of regions with surface class \"%s{1}\" "
                            "are TRANSPARENT for surface molecules  ",
                  sp->sym->name);
        }
        fprintf(log_file, "%s{%d}", no->name, no->orient);
        if (no->next != NULL)
          fprintf(log_file, " ");
      }
      if (borders_transp_title_printed)
        fprintf(log_file, ".\n");
    }
  }

  if (sp->absorb_mols != NULL) {
    /* search for ALL_MOLECULES, ALL_VOLUME_MOLECULES, ALL_SURFACE_MOLECULES */
    for (no = sp->absorb_mols; no != NULL; no = no->next) {
      if (strcmp(no->name, "ALL_VOLUME_MOLECULES") == 0) {
        all_volume_mol_orient = no->orient;
        absorb_mols_all_volume_mol = 1;
      }
      if (strcmp(no->name, "ALL_SURFACE_MOLECULES") == 0) {
        all_surface_mol_orient = no->orient;
        absorb_mols_all_surface_mol = 1;
      }
      if (strcmp(no->name, "ALL_MOLECULES") == 0) {
        all_mol_orient = no->orient;
        absorb_mols_all_mol = 1;
      }
    }

    if ((absorb_mols_all_volume_mol || absorb_mols_all_mol) &&
        (vol_species_name_list != NULL)) {
      fprintf(log_file, "Surfaces with surface class \"%s{1}\" are ABSORPTIVE "
                        "for volume molecules  ",
              sp->sym->name);
      for (nl = vol_species_name_list; nl != NULL; nl = nl->next) {
        if (absorb_mols_all_volume_mol) {
          fprintf(log_file, "%s{%d}", nl->name, all_volume_mol_orient);
        } else if (absorb_mols_all_mol) {
          fprintf(log_file, "%s{%d}", nl->name, all_mol_orient);
        }
        if (nl->next != NULL)
          fprintf(log_file, " ");
        else
          fprintf(log_file, ".");
      }
    }
    fprintf(log_file, "\n");

    if ((absorb_mols_all_surface_mol || absorb_mols_all_mol) &&
        (surf_species_name_list != NULL)) {
      fprintf(log_file, "Borders of regions with surface class \"%s{1}\" are "
                        "ABSORPTIVE for surface molecules  ",
              sp->sym->name);
      for (nl = surf_species_name_list; nl != NULL; nl = nl->next) {
        if (absorb_mols_all_surface_mol) {
          fprintf(log_file, "%s{%d}", nl->name, all_surface_mol_orient);
        } else if (absorb_mols_all_mol) {
          fprintf(log_file, "%s{%d}", nl->name, all_mol_orient);
        }
        if (nl->next != NULL)
          fprintf(log_file, " ");
        else
          fprintf(log_file, ".");
      }
    }
    fprintf(log_file, "\n");

    surf_absorb_title_printed = 0;
    borders_absorb_title_printed = 0;
    /* First go through all volume molecules */
    if ((!absorb_mols_all_mol) && (!absorb_mols_all_volume_mol)) {
      for (no = sp->absorb_mols; no != NULL; no = no->next) {
        spec = get_species_by_name(no->name, n_species, species_list);
        if (spec == NULL)
          mcell_internal_error("Cannot find molecule name %s", no->name);
        if (spec->flags & ON_GRID)
          continue;
        if (strcmp(no->name, "ALL_MOLECULES") == 0)
          continue;
        if (strcmp(no->name, "ALL_VOLUME_MOLECULES") == 0)
          continue;
        if (strcmp(no->name, "ALL_SURFACE_MOLECULES") == 0)
          continue;

        if (!surf_absorb_title_printed) {
          surf_absorb_title_printed = 1;
          fprintf(log_file, "Surfaces with surface class \"%s{1}\" are "
                            "ABSORPTIVE for volume molecules  ",
                  sp->sym->name);
        }
        fprintf(log_file, "%s{%d}", no->name, no->orient);
        if (no->next != NULL)
          fprintf(log_file, " ");
      }
      if (surf_absorb_title_printed)
        fprintf(log_file, ".\n");
    }

    /* Now go through all surface molecules */
    if ((!absorb_mols_all_mol) && (!absorb_mols_all_surface_mol)) {
      for (no = sp->absorb_mols; no != NULL; no = no->next) {
        spec = get_species_by_name(no->name, n_species, species_list);
        if (spec == NULL)
          mcell_internal_error("Cannot find molecule name %s", no->name);
        if ((spec->flags & NOT_FREE) == 0)
          continue;
        if (strcmp(no->name, "ALL_MOLECULES") == 0)
          continue;
        if (strcmp(no->name, "ALL_VOLUME_MOLECULES") == 0)
          continue;
        if (strcmp(no->name, "ALL_SURFACE_MOLECULES") == 0)
          continue;

        if (!borders_absorb_title_printed) {
          borders_absorb_title_printed = 1;
          fprintf(log_file, "Borders of regions with surface class \"%s{1}\" "
                            "are ABSORPTIVE for surface molecules  ",
                  sp->sym->name);
        }
        fprintf(log_file, "%s{%d}", no->name, no->orient);
        if (no->next != NULL)
          fprintf(log_file, " ");
      }
      if (borders_absorb_title_printed)
        fprintf(log_file, ".\n");
    }
  }

  if ((sp->refl_mols != NULL) || (sp->transp_mols != NULL) ||
      (sp->absorb_mols != NULL)) {
    fprintf(log_file, "\n");
  }
}

/***************************************************************************
check_for_conflicts_in_surface_class:
  In: species that is a SURFACE_CLASS
  Out: None.  If species is a SURFACE_CLASS we check for conflicting
       properties, like ABSORPTIVE/TRANSPARENT for the same molecule.
       Also we search for conflicts like combination of
       TRANSPARENT=A and REFLECTIVE=ALL_MOLECULES.
       The combination of regular reaction declared through
       DEFINE_REACTIONS and special reaction (TRANSPARENT/ABSORPTIVE)
       is not allowed for volume molecule.
       If we find the conflicts, the fatal error is generated.
***************************************************************************/
void check_for_conflicts_in_surface_class(struct volume *world,
                                          struct species *sp) {
  struct name_orient *no, *no1, *no2;
  /* orientation of ALL_MOLECULES, ALL_VOLUME_MOLECULES, ALL_SURFACE_MOLECULES
   */
  int all_mol_orient = INT_MIN;
  int all_volume_mol_orient = INT_MIN;
  int all_surface_mol_orient = INT_MIN;
  /* flags */
  int refl_mols_all_volume_mol = 0;
  int refl_mols_all_surface_mol = 0;
  int refl_mols_all_mol = 0;
  int transp_mols_all_volume_mol = 0;
  int transp_mols_all_surface_mol = 0;
  int transp_mols_all_mol = 0;
  int absorb_mols_all_volume_mol = 0;
  int absorb_mols_all_surface_mol = 0;
  int absorb_mols_all_mol = 0;
  struct species *spec;

  if ((sp->flags & IS_SURFACE) == 0)
    return;

  /* search for ALL_MOLECULES, ALL_VOLUME_MOLECULES, ALL_SURFACE_MOLECULES */
  if (sp->refl_mols != NULL) {
    for (no = sp->refl_mols; no != NULL; no = no->next) {
      if (strcmp(no->name, "ALL_VOLUME_MOLECULES") == 0) {
        all_volume_mol_orient = no->orient;
        refl_mols_all_volume_mol = 1;
      }
      if (strcmp(no->name, "ALL_SURFACE_MOLECULES") == 0) {
        all_surface_mol_orient = no->orient;
        refl_mols_all_surface_mol = 1;
      }
      if (strcmp(no->name, "ALL_MOLECULES") == 0) {
        all_mol_orient = no->orient;
        refl_mols_all_mol = 1;
      }
    }
  }

  if (refl_mols_all_mol && refl_mols_all_volume_mol) {
    mcell_error(
        "REFLECTIVE properties are specified through the use of both"
        " ALL_MOLECULES and ALL_VOLUME_MOLECULES within the surface class "
        "'%s'.",
        sp->sym->name);
  }

  if (refl_mols_all_mol && refl_mols_all_surface_mol) {
    mcell_error(
        "REFLECTIVE properties are specified through the use of both"
        " ALL_MOLECULES and ALL_SURFACE_MOLECULES within the surface class "
        "'%s'.",
        sp->sym->name);
  }

  if (sp->transp_mols != NULL) {
    for (no = sp->transp_mols; no != NULL; no = no->next) {
      if (strcmp(no->name, "ALL_VOLUME_MOLECULES") == 0) {
        all_volume_mol_orient = no->orient;
        transp_mols_all_volume_mol = 1;
      }
      if (strcmp(no->name, "ALL_SURFACE_MOLECULES") == 0) {
        all_surface_mol_orient = no->orient;
        transp_mols_all_surface_mol = 1;
      }
      if (strcmp(no->name, "ALL_MOLECULES") == 0) {
        all_mol_orient = no->orient;
        transp_mols_all_mol = 1;
      }
    }
  }

  if (transp_mols_all_mol && transp_mols_all_volume_mol) {
    mcell_error(
        "TRANSPARENT properties are specified through the use of both "
        "ALL_MOLECULES and ALL_VOLUME_MOLECULES within the surface class "
        "'%s'.",
        sp->sym->name);
  }

  if (transp_mols_all_mol && transp_mols_all_surface_mol) {
    mcell_error(
        "TRANSPARENT properties are specified through the use of both "
        "ALL_MOLECULES and ALL_SURFACE_MOLECULES within the surface class "
        "'%s'.",
        sp->sym->name);
  }

  if (sp->absorb_mols != NULL) {
    for (no = sp->absorb_mols; no != NULL; no = no->next) {
      if (strcmp(no->name, "ALL_VOLUME_MOLECULES") == 0) {
        all_volume_mol_orient = no->orient;
        absorb_mols_all_volume_mol = 1;
      }
      if (strcmp(no->name, "ALL_SURFACE_MOLECULES") == 0) {
        all_surface_mol_orient = no->orient;
        absorb_mols_all_surface_mol = 1;
      }
      if (strcmp(no->name, "ALL_MOLECULES") == 0) {
        all_mol_orient = no->orient;
        absorb_mols_all_mol = 1;
      }
    }
  }

  if (absorb_mols_all_mol && absorb_mols_all_volume_mol) {
    mcell_error("ABSORPTIVE properties are specified through the use of both"
                " ALL_MOLECULES and ALL_VOLUME_MOLECULES within the surface "
                "class '%s'.",
                sp->sym->name);
  }

  if (absorb_mols_all_mol && absorb_mols_all_surface_mol) {
    mcell_error("ABSORPTIVE properties are specified through the use of both"
                " ALL_MOLECULES and ALL_SURFACE_MOLECULES within the surface"
                " class '%s'.",
                sp->sym->name);
  }

  for (no = sp->absorb_mols; no != NULL; no = no->next) {
    if (transp_mols_all_mol) {
      if ((all_mol_orient == no->orient) || (all_mol_orient == 0) ||
          (no->orient == 0)) {
        mcell_error("TRANSPARENT and ABSORPTIVE properties are "
                    "simultaneously specified for molecule %s through the "
                    "use of ALL_MOLECULES within the surface class '%s'.",
                    no->name, sp->sym->name);
      }
    }
    if (transp_mols_all_volume_mol) {
      spec =
          get_species_by_name(no->name, world->n_species, world->species_list);
      if (spec == NULL)
        mcell_internal_error("Cannot find molecule name %s", no->name);
      if (spec->flags & ON_GRID)
        continue;

      if ((all_volume_mol_orient == no->orient) ||
          (all_volume_mol_orient == 0) || (no->orient == 0)) {
        mcell_error("TRANSPARENT and ABSORPTIVE properties are "
                    "simultaneously specified for molecule %s through the "
                    "use of ALL_VOLUME_MOLECULES within the surface class "
                    "'%s'.",
                    no->name, sp->sym->name);
      }
    }
    if (transp_mols_all_surface_mol) {
      spec =
          get_species_by_name(no->name, world->n_species, world->species_list);
      if (spec == NULL)
        mcell_internal_error("Cannot find molecule name %s", no->name);
      if ((spec->flags & NOT_FREE) == 0)
        continue;

      if ((all_surface_mol_orient == no->orient) ||
          (all_surface_mol_orient == 0) || (no->orient == 0)) {
        mcell_error("TRANSPARENT and ABSORPTIVE properties are "
                    "simultaneously specified for molecule %s through the "
                    "use of ALL_SURFACE_MOLECULES within the surface class"
                    " '%s'.",
                    no->name, sp->sym->name);
      }
    }

    if (refl_mols_all_mol) {
      if ((all_mol_orient == no->orient) || (all_mol_orient == 0) ||
          (no->orient == 0)) {
        mcell_error("REFLECTIVE and ABSORPTIVE properties are "
                    "simultaneously specified for molecule %s through the "
                    "use of ALL_MOLECULES within the surface class '%s'.",
                    no->name, sp->sym->name);
      }
    }
    if (refl_mols_all_volume_mol) {
      spec =
          get_species_by_name(no->name, world->n_species, world->species_list);
      if (spec == NULL)
        mcell_internal_error("Cannot find molecule name %s", no->name);
      if (spec->flags & ON_GRID)
        continue;

      if ((all_volume_mol_orient == no->orient) ||
          (all_volume_mol_orient == 0) || (no->orient == 0)) {
        mcell_error("REFLECTIVE and ABSORPTIVE properties are "
                    "simultaneously specified for molecule %s through the "
                    "use of ALL_VOLUME_MOLECULES within the surface class "
                    "'%s'.",
                    no->name, sp->sym->name);
      }
    }
    if (refl_mols_all_surface_mol) {
      spec =
          get_species_by_name(no->name, world->n_species, world->species_list);
      if (spec == NULL)
        mcell_internal_error("Cannot find molecule name %s", no->name);
      if ((spec->flags & NOT_FREE) == 0)
        continue;

      if ((all_surface_mol_orient == no->orient) ||
          (all_surface_mol_orient == 0) || (no->orient == 0)) {
        mcell_error("REFLECTIVE and ABSORPTIVE properties are "
                    "simultaneously specified for molecule %s through the "
                    "use of ALL_SURFACE_MOLECULES within the surface class"
                    " '%s'.",
                    no->name, sp->sym->name);
      }
    }
  }

  for (no = sp->transp_mols; no != NULL; no = no->next) {
    if (absorb_mols_all_mol) {
      if ((all_mol_orient == no->orient) || (all_mol_orient == 0) ||
          (no->orient == 0)) {
        mcell_error("TRANSPARENT and ABSORPTIVE properties are "
                    "simultaneously specified for molecule %s through the "
                    "use of ALL_MOLECULES within the surface class '%s'.",
                    no->name, sp->sym->name);
      }
    }
    if (absorb_mols_all_volume_mol) {
      spec =
          get_species_by_name(no->name, world->n_species, world->species_list);
      if (spec == NULL)
        mcell_internal_error("Cannot find molecule name %s", no->name);
      if (spec->flags & ON_GRID)
        continue;

      if ((all_volume_mol_orient == no->orient) ||
          (all_volume_mol_orient == 0) || (no->orient == 0)) {
        mcell_error("TRANSPARENT and ABSORPTIVE properties are "
                    "simultaneously specified for molecule %s through the "
                    "use of ALL_VOLUME_MOLECULES within the surface class "
                    "'%s'.",
                    no->name, sp->sym->name);
      }
    }
    if (absorb_mols_all_surface_mol) {
      spec =
          get_species_by_name(no->name, world->n_species, world->species_list);
      if (spec == NULL)
        mcell_internal_error("Cannot find molecule name %s", no->name);
      if ((spec->flags & NOT_FREE) == 0)
        continue;

      if ((all_surface_mol_orient == no->orient) ||
          (all_surface_mol_orient == 0) || (no->orient == 0)) {
        mcell_error("TRANSPARENT and ABSORPTIVE properties are "
                    "simultaneously specified for molecule %s through the "
                    "use of ALL_SURFACE_MOLECULES within the surface class"
                    " '%s'.",
                    no->name, sp->sym->name);
      }
    }

    if (refl_mols_all_mol) {
      if ((all_mol_orient == no->orient) || (all_mol_orient == 0) ||
          (no->orient == 0)) {
        mcell_error("REFLECTIVE and TRANSPARENT properties are "
                    "simultaneously specified for molecule %s through the "
                    "use of ALL_MOLECULES within the surface class '%s'.",
                    no->name, sp->sym->name);
      }
    }
    if (refl_mols_all_volume_mol) {
      spec =
          get_species_by_name(no->name, world->n_species, world->species_list);
      if (spec == NULL)
        mcell_internal_error("Cannot find molecule name %s", no->name);
      if (spec->flags & ON_GRID)
        continue;

      if ((all_volume_mol_orient == no->orient) ||
          (all_volume_mol_orient == 0) || (no->orient == 0)) {
        mcell_error("REFLECTIVE and TRANSPARENT properties are "
                    "simultaneously specified for molecule %s through the "
                    "use of ALL_VOLUME_MOLECULES within the surface class "
                    "'%s'.",
                    no->name, sp->sym->name);
      }
    }
    if (refl_mols_all_surface_mol) {
      spec =
          get_species_by_name(no->name, world->n_species, world->species_list);
      if (spec == NULL)
        mcell_internal_error("Cannot find molecule name %s", no->name);
      if ((spec->flags & NOT_FREE) == 0)
        continue;

      if ((all_surface_mol_orient == no->orient) ||
          (all_surface_mol_orient == 0) || (no->orient == 0)) {
        mcell_error("REFLECTIVE and TRANSPARENT properties are "
                    "simultaneously specified for molecule %s through the "
                    "use of ALL_SURFACE_MOLECULES within the surface class"
                    " '%s'.",
                    no->name, sp->sym->name);
      }
    }
  }

  for (no = sp->refl_mols; no != NULL; no = no->next) {
    if (absorb_mols_all_mol) {
      if ((all_mol_orient == no->orient) || (all_mol_orient == 0) ||
          (no->orient == 0)) {
        mcell_error("REFLECTIVE and ABSORPTIVE properties are "
                    "simultaneously specified for molecule %s through the "
                    "use of ALL_MOLECULES within the surface class '%s'.",
                    no->name, sp->sym->name);
      }
    }
    if (absorb_mols_all_volume_mol) {
      spec =
          get_species_by_name(no->name, world->n_species, world->species_list);
      if (spec == NULL)
        mcell_internal_error("Cannot find molecule name %s", no->name);
      if (spec->flags & ON_GRID)
        continue;

      if ((all_volume_mol_orient == no->orient) ||
          (all_volume_mol_orient == 0) || (no->orient == 0)) {
        mcell_error("REFLECTIVE and ABSORPTIVE properties are "
                    "simultaneously specified for molecule %s through the "
                    "use of ALL_VOLUME_MOLECULES within the surface class "
                    "'%s'.",
                    no->name, sp->sym->name);
      }
    }
    if (absorb_mols_all_surface_mol) {
      spec =
          get_species_by_name(no->name, world->n_species, world->species_list);
      if (spec == NULL)
        mcell_internal_error("Cannot find molecule name %s", no->name);
      if ((spec->flags & NOT_FREE) == 0)
        continue;

      if ((all_surface_mol_orient == no->orient) ||
          (all_surface_mol_orient == 0) || (no->orient == 0)) {
        mcell_error("REFLECTIVE and ABSORPTIVE properties are "
                    "simultaneously specified for molecule %s through the "
                    "use of ALL_SURFACE_MOLECULES within the surface class"
                    " '%s'.",
                    no->name, sp->sym->name);
      }
    }

    if (transp_mols_all_mol) {
      if ((all_mol_orient == no->orient) || (all_mol_orient == 0) ||
          (no->orient == 0)) {
        mcell_error("REFLECTIVE and TRANSPARENT properties are "
                    "simultaneously specified for molecule %s through the "
                    "use of ALL_MOLECULES within the surface class '%s'.",
                    no->name, sp->sym->name);
      }
    }
    if (transp_mols_all_volume_mol) {
      spec =
          get_species_by_name(no->name, world->n_species, world->species_list);
      if (spec == NULL)
        mcell_internal_error("Cannot find molecule name %s", no->name);
      if (spec->flags & ON_GRID)
        continue;

      if ((all_volume_mol_orient == no->orient) ||
          (all_volume_mol_orient == 0) || (no->orient == 0)) {
        mcell_error("REFLECTIVE and TRANSPARENT properties are "
                    "simultaneously specified for molecule %s through the "
                    "use of ALL_VOLUME_MOLECULES within the surface class "
                    "'%s'.",
                    no->name, sp->sym->name);
      }
    }
    if (transp_mols_all_surface_mol) {
      spec =
          get_species_by_name(no->name, world->n_species, world->species_list);
      if (spec == NULL)
        mcell_internal_error("Cannot find molecule name %s", no->name);
      if ((spec->flags & NOT_FREE) == 0)
        continue;

      if ((all_surface_mol_orient == no->orient) ||
          (all_surface_mol_orient == 0) || (no->orient == 0)) {
        mcell_error("REFLECTIVE and TRANSPARENT properties are "
                    "simultaneously specified for molecule %s through "
                    "the use of ALL_SURFACE_MOLECULES within the surface "
                    "class '%s'.",
                    no->name, sp->sym->name);
      }
    }
  }

  for (no = sp->clamp_conc_mols; no != NULL; no = no->next) {
    if (absorb_mols_all_mol) {
      if ((all_mol_orient == no->orient) || (all_mol_orient == 0) ||
          (no->orient == 0)) {
        mcell_error("CLAMP_CONCENTRATION and ABSORPTIVE properties are "
                    "simultaneously specified for molecule %s through the "
                    "use of ALL_MOLECULES within the surface class '%s'.",
                    no->name, sp->sym->name);
      }
    }
    if (absorb_mols_all_volume_mol) {
      spec =
          get_species_by_name(no->name, world->n_species, world->species_list);
      if (spec == NULL)
        mcell_internal_error("Cannot find molecule name %s", no->name);
      if (spec->flags & ON_GRID)
        continue;

      if ((all_volume_mol_orient == no->orient) ||
          (all_volume_mol_orient == 0) || (no->orient == 0)) {
        mcell_error("CLAMP_CONCENTRATION and ABSORPTIVE properties are "
                    "simultaneously specified for molecule %s through the "
                    "use of ALL_VOLUME_MOLECULES within the surface class "
                    "'%s'.",
                    no->name, sp->sym->name);
      }
    }

    if (transp_mols_all_mol) {
      if ((all_mol_orient == no->orient) || (all_mol_orient == 0) ||
          (no->orient == 0)) {
        mcell_error("CLAMP_CONCENTRATION and TRANSPARENT properties "
                    "are simultaneously specified for molecule %s through "
                    "the use of ALL_MOLECULES within the surface class "
                    "'%s'.",
                    no->name, sp->sym->name);
      }
    }
    if (transp_mols_all_volume_mol) {
      spec =
          get_species_by_name(no->name, world->n_species, world->species_list);
      if (spec == NULL)
        mcell_internal_error("Cannot find molecule name %s", no->name);
      if (spec->flags & ON_GRID)
        continue;

      if ((all_volume_mol_orient == no->orient) ||
          (all_volume_mol_orient == 0) || (no->orient == 0)) {
        mcell_error("CLAMP_CONCENTRATION and TRANSPARENT properties are "
                    "simultaneously specified for molecule %s through the "
                    "use of ALL_VOLUME_MOLECULES within the surface class "
                    "'%s'.",
                    no->name, sp->sym->name);
      }
    }

    if (refl_mols_all_mol) {
      if ((all_mol_orient == no->orient) || (all_mol_orient == 0) ||
          (no->orient == 0)) {
        mcell_error("CLAMP_CONCENTRATION and REFLECTIVE properties are "
                    "simultaneously specified for molecule %s through the "
                    "use of ALL_MOLECULES within the surface class '%s'.",
                    no->name, sp->sym->name);
      }
    }
    if (refl_mols_all_volume_mol) {
      spec =
          get_species_by_name(no->name, world->n_species, world->species_list);
      if (spec == NULL)
        mcell_internal_error("Cannot find molecule name %s", no->name);
      if (spec->flags & ON_GRID)
        continue;

      if ((all_volume_mol_orient == no->orient) ||
          (all_volume_mol_orient == 0) || (no->orient == 0)) {
        mcell_error("CLAMP_CONCENTRATION and REFLECTIVE properties are "
                    "simultaneously specified for molecule %s through the "
                    "use of ALL_VOLUME_MOLECULES within the surface class "
                    "'%s'.",
                    no->name, sp->sym->name);
      }
    }
  }

  for (no1 = sp->transp_mols; no1 != NULL; no1 = no1->next) {
    for (no2 = sp->absorb_mols; no2 != NULL; no2 = no2->next) {
      if (strcmp(no1->name, no2->name) == 0) {
        if ((no1->orient == no2->orient) || (no1->orient == 0) ||
            (no2->orient == 0)) {
          mcell_error("TRANSPARENT and ABSORPTIVE properties are "
                      "simultaneously specified for the same molecule %s "
                      "within the surface class '%s'.",
                      no1->name, sp->sym->name);
        }
      }
    }
    for (no2 = sp->refl_mols; no2 != NULL; no2 = no2->next) {
      if (strcmp(no1->name, no2->name) == 0) {
        if ((no1->orient == no2->orient) || (no1->orient == 0) ||
            (no2->orient == 0)) {
          mcell_error("TRANSPARENT and REFLECTIVE properties are "
                      "simultaneously specified for the same molecule %s "
                      "within the surface class '%s'.",
                      no1->name, sp->sym->name);
        }
      }
    }
    for (no2 = sp->clamp_conc_mols; no2 != NULL; no2 = no2->next) {
      if (strcmp(no1->name, no2->name) == 0) {
        if ((no1->orient == no2->orient) || (no1->orient == 0) ||
            (no2->orient == 0)) {
          mcell_error("TRANSPARENT and CLAMP_CONCENTRATION properties are "
                      "simultaneously specified for the same molecule %s "
                      "within the surface class '%s'.",
                      no1->name, sp->sym->name);
        }
      }
    }
  }

  for (no1 = sp->refl_mols; no1 != NULL; no1 = no1->next) {
    for (no2 = sp->absorb_mols; no2 != NULL; no2 = no2->next) {
      if (strcmp(no1->name, no2->name) == 0) {
        if ((no1->orient == no2->orient) || (no1->orient == 0) ||
            (no2->orient == 0)) {
          mcell_error(
              "REFLECTIVE and ABSORPTIVE properties are "
              "simultaneously specified for the same molecule %s within the "
              "surface class '%s'.",
              no1->name, sp->sym->name);
        }
      }
    }
    for (no2 = sp->clamp_conc_mols; no2 != NULL; no2 = no2->next) {
      if (strcmp(no1->name, no2->name) == 0) {
        if ((no1->orient == no2->orient) || (no1->orient == 0) ||
            (no2->orient == 0)) {
          mcell_error(
              "REFLECTIVE and CLAMP_CONCENTRATION properties are "
              "simultaneously specified for the same molecule %s within the "
              "surface class '%s'.",
              no1->name, sp->sym->name);
        }
      }
    }
  }

  for (no1 = sp->absorb_mols; no1 != NULL; no1 = no1->next) {
    for (no2 = sp->clamp_conc_mols; no2 != NULL; no2 = no2->next) {
      if (strcmp(no1->name, no2->name) == 0) {
        if ((no1->orient == no2->orient) || (no1->orient == 0) ||
            (no2->orient == 0)) {
          mcell_error(
              "ABSORPTIVE and CLAMP_CONCENTRATION properties are "
              "simultaneously specified for the same molecule %s within the "
              "surface class '%s'.",
              no1->name, sp->sym->name);
        }
      }
    }
  }

  /* Check for combinations of ALL_MOLECULES or ALL_VOLUME_MOLECULES
    with regular reactions that involve volume molecules. */
  /* We will check for volume molecules only
    since special reactions with surface molecules
    have a meaning as reactions on the region border
    and should be allowed. */

  u_int hash_value, hashW, hashM;
  struct species *mol_sp;
  struct rxn *inter;
  /* flags */
  int special_rx_same_orient, regular_rx_same_orient;
  int count_sim_orient_rxns, i0, i1;

  hashW = sp->hashval;

  for (int i = 0; i < world->n_species; i++) {
    if (world->species_list[i]->flags & IS_SURFACE)
      continue;
    mol_sp = world->species_list[i];
    if (strcmp(mol_sp->sym->name, "ALL_MOLECULES") == 0)
      continue;
    if (strcmp(mol_sp->sym->name, "ALL_VOLUME_MOLECULES") == 0)
      continue;
    if (strcmp(mol_sp->sym->name, "ALL_SURFACE_MOLECULES") == 0)
      continue;
    if (mol_sp->flags & ON_GRID)
      continue;

    hashM = mol_sp->hashval;
    hash_value = (hashM + hashW) & (world->rx_hashsize - 1);

    /* Are there regular reactions with volume molecules? */
    for (inter = world->reaction_hash[hash_value]; inter != NULL;
         inter = inter->next) {
      if ((inter->players[0]->flags & ON_GRID) != 0)
        continue;
      if (strcmp(inter->players[0]->sym->name, mol_sp->sym->name) != 0)
        continue;

      if (transp_mols_all_mol) {
        if ((all_mol_orient == inter->geometries[0]) || (all_mol_orient == 0) ||
            (inter->geometries[0] == 0)) {
          mcell_error(
              "Combination of similar oriented TRANSPARENT reaction "
              "using ALL_MOLECULES and regular reaction for molecule '%s' "
              "for the same surface class '%s' is not allowed.",
              inter->players[0]->sym->name, sp->sym->name);
        }
      }
      if (absorb_mols_all_mol) {
        if ((all_mol_orient == inter->geometries[0]) || (all_mol_orient == 0) ||
            (inter->geometries[0] == 0)) {
          mcell_error(
              "Combination of similar oriented ABSORPTIVE reaction "
              "using ALL_MOLECULES and regular reaction for molecule '%s' "
              "for the same surface class '%s' is not allowed.",
              inter->players[0]->sym->name, sp->sym->name);
        }
      }
      if (transp_mols_all_volume_mol) {
        if ((all_volume_mol_orient == inter->geometries[0]) ||
            (all_volume_mol_orient == 0) || (inter->geometries[0] == 0)) {
          mcell_error(
              "Combination of similar oriented TRANSPARENT reaction "
              "using ALL_VOLUME_MOLECULES and regular reaction for molecule "
              "'%s' for the same surface class '%s' is not allowed.",
              inter->players[0]->sym->name, sp->sym->name);
        }
      }
      if (absorb_mols_all_volume_mol) {
        if ((all_volume_mol_orient == inter->geometries[0]) ||
            (all_volume_mol_orient == 0) || (inter->geometries[0] == 0)) {
          mcell_error(
              "Combination of similar oriented ABSORPTIVE reaction "
              "using ALL_VOLUME_MOLECULES and regular reaction for molecule "
              "'%s' for the same surface class '%s' is not allowed.",
              inter->players[0]->sym->name, sp->sym->name);
        }
      }
    }
  }

  /* Below we will check for the conflicting reactions, like
    A; @ surf_class; -> TRANSPARENT and
    A; @ surf_class -> B[some_rate]
  */
  /* We will check for volume molecules only
    since special reactions with surface molecules
    have a meaning as reactions on the region border
    and should be allowed. */

  for (no = sp->transp_mols; no != NULL; no = no->next) {
    /* surface_class always has orientation {1} */
    if (no->orient == 1)
      special_rx_same_orient = 1;
    else
      special_rx_same_orient = 0;

    char *mol_name = no->name;
    if (strcmp(mol_name, "ALL_MOLECULES") == 0)
      continue;
    if (strcmp(mol_name, "ALL_VOLUME_MOLECULES") == 0)
      continue;
    if (strcmp(mol_name, "ALL_SURFACE_MOLECULES") == 0)
      continue;

    mol_sp =
        get_species_by_name(mol_name, world->n_species, world->species_list);
    if (mol_sp == NULL)
      mcell_error("Cannot find '%s' among molecules.", mol_name);
    if ((mol_sp->flags & ON_GRID) != 0)
      continue;
    hashM = mol_sp->hashval;

    hash_value = (hashM + hashW) & (world->rx_hashsize - 1);

    /* below we will compare TRANSP reaction only with
       regular reaction for the same molecule. */
    for (inter = world->reaction_hash[hash_value]; inter != NULL;
         inter = inter->next) {
      if (inter->n_pathways <= RX_SPECIAL) {
        continue;
      }

      if ((inter->players[0] == mol_sp) && (inter->players[1] == sp)) {
        if (inter->geometries[0] == inter->geometries[1]) {
          regular_rx_same_orient = 1;
        } else {
          regular_rx_same_orient = 0;
        }
        if ((no->orient == 0) || (inter->geometries[0] == 0)) {
          mcell_error("Combination of similar oriented TRANSPARENT and regular "
                      "reactions for molecule '%s' on the same surface class "
                      "'%s' is not allowed.",
                      mol_sp->sym->name, sp->sym->name);
        }
        if (special_rx_same_orient && regular_rx_same_orient) {
          mcell_error("Combination of similar oriented TRANSPARENT and regular "
                      "reactions for molecule '%s' on the same surface class "
                      "'%s' is not allowed.",
                      mol_sp->sym->name, sp->sym->name);
        }
        if (!special_rx_same_orient && !regular_rx_same_orient) {
          mcell_error("Combination of similar oriented TRANSPARENT and regular "
                      "reactions for molecule '%s' on the same surface class "
                      "'%s' is not allowed.",
                      mol_sp->sym->name, sp->sym->name);
        }
      }
    }
  }

  /* Below we will compare ABSORPTIVE reaction only with
     regular reaction.  The difficulty is that internally
     ABSORPTIVE reaction is implemented as regular reaction
     with no products. */
  for (no = sp->absorb_mols; no != NULL; no = no->next) {

    char *mol_name = no->name;
    if (strcmp(mol_name, "ALL_MOLECULES") == 0)
      continue;
    if (strcmp(mol_name, "ALL_VOLUME_MOLECULES") == 0)
      continue;
    if (strcmp(mol_name, "ALL_SURFACE_MOLECULES") == 0)
      continue;
    mol_sp =
        get_species_by_name(mol_name, world->n_species, world->species_list);
    if (mol_sp == NULL)
      mcell_error("Cannot find '%s' among molecules.", mol_name);
    /* we will check for volume molecules only
       since special reactions with surface molecules
       have a meaning as reactions on the region border
       and should be allowed. */
    if ((mol_sp->flags & ON_GRID) != 0)
      continue;

    if (no->orient == 1)
      special_rx_same_orient = 1;
    else
      special_rx_same_orient = 0;

    hashM = mol_sp->hashval;
    hash_value = (hashM + hashW) & (world->rx_hashsize - 1);

    /* count similar oriented ABSORPTIVE reactions */
    count_sim_orient_rxns = 0;
    for (inter = world->reaction_hash[hash_value]; inter != NULL;
         inter = inter->next) {
      if (inter->n_pathways <= RX_SPECIAL) {
        continue;
      }

      if ((inter->players[0] == mol_sp) && (inter->players[1] == sp)) {
        if ((no->orient == inter->geometries[0]) &&
            (inter->geometries[1] == 1) && (inter->n_pathways == 1)) {
          i0 = inter->product_idx[0];
          i1 = inter->product_idx[1];
          if (((i1 - i0) == 2) && (inter->players[2] == NULL) &&
              (inter->players[3] == NULL)) {
            /* this is an original ABSORPTIVE reaction */
            continue;
          }
        }

        if (inter->geometries[0] == inter->geometries[1]) {
          regular_rx_same_orient = 1;
        } else {
          regular_rx_same_orient = 0;
        }

        if ((no->orient == 0) || (inter->geometries[0] == 0)) {
          count_sim_orient_rxns++;
        } else if (special_rx_same_orient && regular_rx_same_orient) {
          count_sim_orient_rxns++;
        } else if (!special_rx_same_orient && !regular_rx_same_orient) {
          count_sim_orient_rxns++;
        }
      }
      if (count_sim_orient_rxns > 0) {
        mcell_error("Combination of similar oriented ABSORPTIVE reaction and "
                    "regular reaction for molecule '%s' on the same surface "
                    "class '%s' is not allowed.",
                    mol_name, sp->sym->name);
      }
    }
  }

  /* Internally CLAMP_CONC reaction is implemented as regular reaction */
  for (no = sp->clamp_conc_mols; no != NULL; no = no->next) {
    char *mol_name = no->name;
    if (strcmp(mol_name, "ALL_MOLECULES") == 0)
      continue;
    if (strcmp(mol_name, "ALL_VOLUME_MOLECULES") == 0)
      continue;
    if (strcmp(mol_name, "ALL_SURFACE_MOLECULES") == 0)
      continue;
    mol_sp =
        get_species_by_name(mol_name, world->n_species, world->species_list);
    if (mol_sp == NULL)
      mcell_error("Cannot find '%s' among molecules.", mol_name);
    /* we will check for volume molecules only
       since special reactions with surface molecules
       have a meaning as reactions on the region border
       and should be allowed. */
    if ((mol_sp->flags & ON_GRID) != 0)
      continue;

    if (no->orient == 1)
      special_rx_same_orient = 1;
    else
      special_rx_same_orient = 0;

    hashM = mol_sp->hashval;
    hash_value = (hashM + hashW) & (world->rx_hashsize - 1);

    /* count similar oriented CLAMP_CONC and regular reactions */
    count_sim_orient_rxns = 0;
    for (inter = world->reaction_hash[hash_value]; inter != NULL;
         inter = inter->next) {
      if (inter->n_pathways <= RX_SPECIAL) {
        continue;
      }

      if ((inter->players[0] == mol_sp) && (inter->players[1] == sp)) {
        if ((no->orient == inter->geometries[0]) &&
            (inter->geometries[1] == 1)) {
          i0 = inter->product_idx[0];
          i1 = inter->product_idx[1];
          if (((i1 - i0) == 2) && (inter->players[2] == NULL) &&
              (inter->players[3] == NULL)) {
            /* this is an original CLAMP_CONC reaction */
            continue;
          }
        }

        if (inter->geometries[0] == inter->geometries[1]) {
          regular_rx_same_orient = 1;
        } else {
          regular_rx_same_orient = 0;
        }

        if ((no->orient == 0) || (inter->geometries[0] == 0)) {
          count_sim_orient_rxns++;
        } else if (special_rx_same_orient && regular_rx_same_orient) {
          count_sim_orient_rxns++;
        } else if (!special_rx_same_orient && !regular_rx_same_orient) {
          count_sim_orient_rxns++;
        }
      }
      if (count_sim_orient_rxns > 0) {
        mcell_error("Combination of similar oriented CLAMP_CONCENTRATION "
                    "reaction and regular reaction for molecule '%s' on the "
                    "same surface class '%s' is not allowed.",
                    mol_name, sp->sym->name);
      }
    }
  }
}

/***************************************************************************
check_for_conflicting_surface_classes:
  In: wall
  Out: The wall has a linked list of surface classes. Here we check for
       conflicts, like ABSORPTIVE/TRANSPARENT for the same molecule
       between different surface classes. Possible application -
       overlapping regions.
       If we find the conflicts, the fatal error is generated.
***************************************************************************/
void check_for_conflicting_surface_classes(struct wall *w, int n_species,
                                           struct species **species_list) {

  struct surf_class_list *scl, *scl2;
  struct species *sp, *sp2;
  struct name_orient *no, *no2;
  /* orientation of ALL_MOLECULES, ALL_VOLUME_MOLECULES, ALL_SURFACE_MOLECULES
   */
  int sp_refl_all_mols_orient, sp_transp_all_mols_orient,
      sp_absorb_all_mols_orient;
  int sp_refl_all_volume_mols_orient, sp_transp_all_volume_mols_orient,
      sp_absorb_all_volume_mols_orient;
  int sp_refl_all_surface_mols_orient, sp_transp_all_surface_mols_orient,
      sp_absorb_all_surface_mols_orient;
  int sp2_refl_all_mols_orient, sp2_transp_all_mols_orient,
      sp2_absorb_all_mols_orient;
  int sp2_refl_all_volume_mols_orient, sp2_transp_all_volume_mols_orient,
      sp2_absorb_all_volume_mols_orient;
  int sp2_refl_all_surface_mols_orient, sp2_transp_all_surface_mols_orient,
      sp2_absorb_all_surface_mols_orient;
  /* flags */
  int sp_refl_all_mols, sp_transp_all_mols, sp_absorb_all_mols;
  int sp_refl_all_volume_mols, sp_transp_all_volume_mols,
      sp_absorb_all_volume_mols;
  int sp_refl_all_surface_mols, sp_transp_all_surface_mols,
      sp_absorb_all_surface_mols;
  int sp2_refl_all_mols, sp2_transp_all_mols, sp2_absorb_all_mols;
  int sp2_refl_all_volume_mols, sp2_transp_all_volume_mols,
      sp2_absorb_all_volume_mols;
  int sp2_refl_all_surface_mols, sp2_transp_all_surface_mols,
      sp2_absorb_all_surface_mols;
  struct species *spec;

  /* check for conflicts like TRANSPARENT/ABSORPTIVE type for the same molecule
   */
  for (scl = w->surf_class_head; scl != NULL; scl = scl->next) {
    sp = scl->surf_class;
    for (scl2 = scl->next; scl2 != NULL; scl2 = scl2->next) {
      sp2 = scl2->surf_class;

      for (no = sp->transp_mols; no != NULL; no = no->next) {
        for (no2 = sp2->absorb_mols; no2 != NULL; no2 = no2->next) {
          if (strcmp(no->name, no2->name) == 0) {
            if ((no->orient == no2->orient) || (no->orient == 0) ||
                (no2->orient == 0)) {
              mcell_error("Conflicting TRANSPARENT and ABSORPTIVE properties "
                          "are simultaneously specified for the same molecule "
                          "'%s' on the same wall through the surface classes "
                          "'%s' and '%s'.",
                          no->name, sp->sym->name, sp2->sym->name);
            }
          }
        }
        for (no2 = sp2->refl_mols; no2 != NULL; no2 = no2->next) {
          if (strcmp(no->name, no2->name) == 0) {
            if ((no->orient == no2->orient) || (no->orient == 0) ||
                (no2->orient == 0)) {
              mcell_error("Conflicting TRANSPARENT and REFLECTIVE properties "
                          "are simultaneously specified for the same molecule "
                          "'%s' on the same wall through the surface classes "
                          "'%s' and '%s'.",
                          no->name, sp->sym->name, sp2->sym->name);
            }
          }
        }
        for (no2 = sp2->clamp_conc_mols; no2 != NULL; no2 = no2->next) {
          if (strcmp(no->name, no2->name) == 0) {
            if ((no->orient == no2->orient) || (no->orient == 0) ||
                (no2->orient == 0)) {
              mcell_error("Conflicting TRANSPARENT and CLAMP_CONCENTRATION "
                          "properties are simultaneously specified for the "
                          "same molecule '%s' on the same wall through the "
                          "surface classes '%s' and '%s'.",
                          no->name, sp->sym->name, sp2->sym->name);
            }
          }
        }
      }

      for (no = sp->refl_mols; no != NULL; no = no->next) {
        for (no2 = sp2->absorb_mols; no2 != NULL; no2 = no2->next) {
          if (strcmp(no->name, no2->name) == 0) {
            if ((no->orient == no2->orient) || (no->orient == 0) ||
                (no2->orient == 0)) {
              mcell_error("Conflicting REFLECTIVE and ABSORPTIVE properties "
                          "are simultaneously specified for the same molecule "
                          "'%s' on the same wall through the surface classes "
                          "'%s' and '%s'.",
                          no->name, sp->sym->name, sp2->sym->name);
            }
          }
        }
        for (no2 = sp2->transp_mols; no2 != NULL; no2 = no2->next) {
          if (strcmp(no->name, no2->name) == 0) {
            if ((no->orient == no2->orient) || (no->orient == 0) ||
                (no2->orient == 0)) {
              mcell_error("Conflicting REFLECTIVE and TRANSPARENT properties "
                          "are simultaneously specified for the same molecule "
                          "'%s' on the same wall through the surface classes "
                          "'%s' and '%s'.",
                          no->name, sp->sym->name, sp2->sym->name);
            }
          }
        }
        for (no2 = sp2->clamp_conc_mols; no2 != NULL; no2 = no2->next) {
          if (strcmp(no->name, no2->name) == 0) {
            if ((no->orient == no2->orient) || (no->orient == 0) ||
                (no2->orient == 0)) {
              mcell_error("Conflicting REFLECTIVE and CLAMP_CONCENTRATION "
                          "properties are simultaneously specified for the "
                          "same molecule '%s' on the same wall through the "
                          "surface classes '%s' and '%s'.",
                          no->name, sp->sym->name, sp2->sym->name);
            }
          }
        }
      }

      for (no = sp->absorb_mols; no != NULL; no = no->next) {
        for (no2 = sp2->transp_mols; no2 != NULL; no2 = no2->next) {
          if (strcmp(no->name, no2->name) == 0) {
            if ((no->orient == no2->orient) || (no->orient == 0) ||
                (no2->orient == 0)) {
              mcell_error("Conflicting ABSORPTIVE and TRANSPARENT properties "
                          "are simultaneously specified for the same molecule "
                          "'%s' on the same wall through the surface classes "
                          "'%s' and '%s'.",
                          no->name, sp->sym->name, sp2->sym->name);
            }
          }
        }
        for (no2 = sp2->refl_mols; no2 != NULL; no2 = no2->next) {
          if (strcmp(no->name, no2->name) == 0) {
            if ((no->orient == no2->orient) || (no->orient == 0) ||
                (no2->orient == 0)) {
              mcell_error("Conflicting ABSORPTIVE and REFLECTIVE properties "
                          "are simultaneously specified for the same molecule "
                          "'%s' on the same wall through the surface classes "
                          "'%s' and '%s'.",
                          no->name, sp->sym->name, sp2->sym->name);
            }
          }
        }
        for (no2 = sp2->clamp_conc_mols; no2 != NULL; no2 = no2->next) {
          if (strcmp(no->name, no2->name) == 0) {
            if ((no->orient == no2->orient) || (no->orient == 0) ||
                (no2->orient == 0)) {
              mcell_error("Conflicting ABSORPTIVE and CLAMP_CONCENTRATION "
                          "properties are simultaneously specified for the "
                          "same molecule '%s' on the same wall through the "
                          "surface classes '%s' and '%s'.",
                          no->name, sp->sym->name, sp2->sym->name);
            }
          }
        }
      }

      for (no = sp->clamp_conc_mols; no != NULL; no = no->next) {
        for (no2 = sp2->transp_mols; no2 != NULL; no2 = no2->next) {
          if (strcmp(no->name, no2->name) == 0) {
            if ((no->orient == no2->orient) || (no->orient == 0) ||
                (no2->orient == 0)) {
              mcell_error("Conflicting CLAMP_CONCENTRATION and TRANSPARENT "
                          "properties are simultaneously specified for the "
                          "same molecule '%s' on the same wall through the "
                          "surface classes '%s' and '%s'.",
                          no->name, sp->sym->name, sp2->sym->name);
            }
          }
        }
        for (no2 = sp2->refl_mols; no2 != NULL; no2 = no2->next) {
          if (strcmp(no->name, no2->name) == 0) {
            if ((no->orient == no2->orient) || (no->orient == 0) ||
                (no2->orient == 0)) {
              mcell_error("Conflicting CLAMP_CONCENTRATION and REFLECTIVE "
                          "properties are simultaneously specified for the "
                          "same molecule '%s' on the same wall through the "
                          "surface classes '%s' and '%s'.",
                          no->name, sp->sym->name, sp2->sym->name);
            }
          }
        }
        for (no2 = sp2->absorb_mols; no2 != NULL; no2 = no2->next) {
          if (strcmp(no->name, no2->name) == 0) {
            if ((no->orient == no2->orient) || (no->orient == 0) ||
                (no2->orient == 0)) {
              mcell_error("Conflicting CLAMP_CONCENTRATION and ABSORPTIVE "
                          "properties are simultaneously specified for the "
                          "same molecule '%s' on the same wall through the "
                          "surface classes '%s' and '%s'.",
                          no->name, sp->sym->name, sp2->sym->name);
            }
          }
        }
      }
    }
  }

  /* check for conflicts resulting from ALL_MOLECULES,
     ALL_VOLUME_MOLECULES, ALL_SURFACE_MOLECULES keywords */
  for (scl = w->surf_class_head; scl != NULL; scl = scl->next) {
    sp = scl->surf_class;

    sp_transp_all_mols = 0;
    sp_absorb_all_mols = 0;
    sp_refl_all_mols = 0;
    sp_transp_all_volume_mols = 0;
    sp_absorb_all_volume_mols = 0;
    sp_refl_all_volume_mols = 0;
    sp_transp_all_surface_mols = 0;
    sp_absorb_all_surface_mols = 0;
    sp_refl_all_surface_mols = 0;

    sp_transp_all_mols_orient = INT_MIN;
    sp_absorb_all_mols_orient = INT_MIN;
    sp_refl_all_mols_orient = INT_MIN;
    sp_transp_all_volume_mols_orient = INT_MIN;
    sp_absorb_all_volume_mols_orient = INT_MIN;
    sp_refl_all_volume_mols_orient = INT_MIN;
    sp_transp_all_surface_mols_orient = INT_MIN;
    sp_absorb_all_surface_mols_orient = INT_MIN;
    sp_refl_all_surface_mols_orient = INT_MIN;

    if (sp->refl_mols != NULL) {
      for (no = sp->refl_mols; no != NULL; no = no->next) {
        if (strcmp(no->name, "ALL_VOLUME_MOLECULES") == 0) {
          sp_refl_all_volume_mols_orient = no->orient;
          sp_refl_all_volume_mols = 1;
        }
        if (strcmp(no->name, "ALL_SURFACE_MOLECULES") == 0) {
          sp_refl_all_surface_mols_orient = no->orient;
          sp_refl_all_surface_mols = 1;
        }
        if (strcmp(no->name, "ALL_MOLECULES") == 0) {
          sp_refl_all_mols_orient = no->orient;
          sp_refl_all_mols = 1;
        }
      }
    }
    if (sp->transp_mols != NULL) {
      for (no = sp->transp_mols; no != NULL; no = no->next) {
        if (strcmp(no->name, "ALL_VOLUME_MOLECULES") == 0) {
          sp_transp_all_volume_mols_orient = no->orient;
          sp_transp_all_volume_mols = 1;
        }
        if (strcmp(no->name, "ALL_SURFACE_MOLECULES") == 0) {
          sp_transp_all_surface_mols_orient = no->orient;
          sp_transp_all_surface_mols = 1;
        }
        if (strcmp(no->name, "ALL_MOLECULES") == 0) {
          sp_transp_all_mols_orient = no->orient;
          sp_transp_all_mols = 1;
        }
      }
    }

    if (sp->absorb_mols != NULL) {
      for (no = sp->absorb_mols; no != NULL; no = no->next) {
        if (strcmp(no->name, "ALL_VOLUME_MOLECULES") == 0) {
          sp_absorb_all_volume_mols_orient = no->orient;
          sp_absorb_all_volume_mols = 1;
        }
        if (strcmp(no->name, "ALL_SURFACE_MOLECULES") == 0) {
          sp_absorb_all_surface_mols_orient = no->orient;
          sp_absorb_all_surface_mols = 1;
        }
        if (strcmp(no->name, "ALL_MOLECULES") == 0) {
          sp_absorb_all_mols_orient = no->orient;
          sp_absorb_all_mols = 1;
        }
      }
    }

    for (scl2 = scl->next; scl2 != NULL; scl2 = scl2->next) {
      sp2 = scl2->surf_class;

      sp2_transp_all_mols = 0;
      sp2_absorb_all_mols = 0;
      sp2_refl_all_mols = 0;
      sp2_transp_all_volume_mols = 0;
      sp2_absorb_all_volume_mols = 0;
      sp2_refl_all_volume_mols = 0;
      sp2_transp_all_surface_mols = 0;
      sp2_absorb_all_surface_mols = 0;
      sp2_refl_all_surface_mols = 0;

      sp2_transp_all_mols_orient = INT_MIN;
      sp2_absorb_all_mols_orient = INT_MIN;
      sp2_refl_all_mols_orient = INT_MIN;
      sp2_transp_all_volume_mols_orient = INT_MIN;
      sp2_absorb_all_volume_mols_orient = INT_MIN;
      sp2_refl_all_volume_mols_orient = INT_MIN;
      sp2_transp_all_surface_mols_orient = INT_MIN;
      sp2_absorb_all_surface_mols_orient = INT_MIN;
      sp2_refl_all_surface_mols_orient = INT_MIN;

      if (sp2->refl_mols != NULL) {
        for (no = sp2->refl_mols; no != NULL; no = no->next) {
          if (strcmp(no->name, "ALL_VOLUME_MOLECULES") == 0) {
            sp2_refl_all_volume_mols_orient = no->orient;
            sp2_refl_all_volume_mols = 1;
          }
          if (strcmp(no->name, "ALL_SURFACE_MOLECULES") == 0) {
            sp2_refl_all_surface_mols_orient = no->orient;
            sp2_refl_all_surface_mols = 1;
          }
          if (strcmp(no->name, "ALL_MOLECULES") == 0) {
            sp2_refl_all_mols_orient = no->orient;
            sp2_refl_all_mols = 1;
          }
        }
      }
      if (sp2->transp_mols != NULL) {
        for (no = sp2->transp_mols; no != NULL; no = no->next) {
          if (strcmp(no->name, "ALL_VOLUME_MOLECULES") == 0) {
            sp2_transp_all_volume_mols_orient = no->orient;
            sp2_transp_all_volume_mols = 1;
          }
          if (strcmp(no->name, "ALL_SURFACE_MOLECULES") == 0) {
            sp2_transp_all_surface_mols_orient = no->orient;
            sp2_transp_all_surface_mols = 1;
          }
          if (strcmp(no->name, "ALL_MOLECULES") == 0) {
            sp2_transp_all_mols_orient = no->orient;
            sp2_transp_all_mols = 1;
          }
        }
      }

      if (sp2->absorb_mols != NULL) {
        for (no = sp2->absorb_mols; no != NULL; no = no->next) {
          if (strcmp(no->name, "ALL_VOLUME_MOLECULES") == 0) {
            sp2_absorb_all_volume_mols_orient = no->orient;
            sp2_absorb_all_volume_mols = 1;
          }
          if (strcmp(no->name, "ALL_SURFACE_MOLECULES") == 0) {
            sp2_absorb_all_surface_mols_orient = no->orient;
            sp2_absorb_all_surface_mols = 1;
          }
          if (strcmp(no->name, "ALL_MOLECULES") == 0) {
            sp2_absorb_all_mols_orient = no->orient;
            sp2_absorb_all_mols = 1;
          }
        }
      }

      /* Below we will check for conflicts when two surface classes
         both have specified ALL_MOLECULES, or ALL_VOLUME_MOLECULES,
         or ALL_SURFACE_MOLECULES keywords in different combinations */
      if (sp_transp_all_mols) {
        if (sp2_absorb_all_mols) {
          if ((sp_transp_all_mols_orient == 0) ||
              (sp2_absorb_all_mols_orient == 0) ||
              (sp_transp_all_mols_orient == sp2_absorb_all_mols_orient)) {
            mcell_error("Conflicting TRANSPARENT and ABSORPTIVE properties are "
                        "simultaneously specified on the same wall using "
                        "ALL_MOLECULES through the surface classes '%s' and "
                        "'%s'.",
                        sp->sym->name, sp2->sym->name);
          }
        }
        if (sp2_absorb_all_volume_mols) {
          if ((sp_transp_all_mols_orient == 0) ||
              (sp2_absorb_all_volume_mols_orient == 0) ||
              (sp_transp_all_mols_orient ==
               sp2_absorb_all_volume_mols_orient)) {
            mcell_error("Conflicting TRANSPARENT and ABSORPTIVE properties are "
                        "simultaneously specified on the same wall using "
                        "ALL_MOLECULES and ALL_VOLUME_MOLECULES through the "
                        "surface classes '%s' and '%s'.",
                        sp->sym->name, sp2->sym->name);
          }
        }
        if (sp2_absorb_all_surface_mols) {
          if ((sp_transp_all_mols_orient == 0) ||
              (sp2_absorb_all_surface_mols_orient == 0) ||
              (sp_transp_all_mols_orient ==
               sp2_absorb_all_surface_mols_orient)) {
            mcell_error("Conflicting TRANSPARENT and ABSORPTIVE properties are "
                        "simultaneously specified on the same wall using "
                        "ALL_MOLECULES and ALL_SURFACE_MOLECULES through the "
                        "surface classes '%s' and '%s'.",
                        sp->sym->name, sp2->sym->name);
          }
        }
        if (sp2_refl_all_mols) {
          if ((sp_transp_all_mols_orient == 0) ||
              (sp2_refl_all_mols_orient == 0) ||
              (sp_transp_all_mols_orient == sp2_refl_all_mols_orient)) {
            mcell_error("Conflicting TRANSPARENT and REFLECTIVE properties are "
                        "simultaneously specified on the same wall using "
                        "ALL_MOLECULES through the surface classes '%s' and "
                        "'%s'.",
                        sp->sym->name, sp2->sym->name);
          }
        }
        if (sp2_refl_all_volume_mols) {
          if ((sp_transp_all_mols_orient == 0) ||
              (sp2_refl_all_volume_mols_orient == 0) ||
              (sp_transp_all_mols_orient == sp2_refl_all_volume_mols_orient)) {
            mcell_error("Conflicting TRANSPARENT and REFLECTIVE properties are "
                        "simultaneously specified on the same wall using "
                        "ALL_MOLECULES and ALL_VOLUME_MOLECULES through the "
                        "surface classes '%s' and '%s'.",
                        sp->sym->name, sp2->sym->name);
          }
        }
        if (sp2_refl_all_surface_mols) {
          if ((sp_transp_all_mols_orient == 0) ||
              (sp2_refl_all_surface_mols_orient == 0) ||
              (sp_transp_all_mols_orient == sp2_refl_all_surface_mols_orient)) {
            mcell_error("Conflicting TRANSPARENT and REFLECTIVE properties are "
                        "simultaneously specified on the same wall using "
                        "ALL_MOLECULES and ALL_SURFACE_MOLECULES through the "
                        "surface classes '%s' and '%s'.",
                        sp->sym->name, sp2->sym->name);
          }
        }
      }
      if (sp_absorb_all_mols) {
        if (sp2_transp_all_mols) {
          if ((sp2_transp_all_mols_orient == 0) ||
              (sp_absorb_all_mols_orient == 0) ||
              (sp2_transp_all_mols_orient == sp_absorb_all_mols_orient)) {
            mcell_error("Conflicting ABSORPTIVE and TRANSPARENT properties are "
                        "simultaneously specified on the same wall using "
                        "ALL_MOLECULES through the surface classes '%s' and "
                        "'%s'.",
                        sp->sym->name, sp2->sym->name);
          }
        }
        if (sp2_transp_all_volume_mols) {
          if ((sp2_transp_all_volume_mols_orient == 0) ||
              (sp_absorb_all_mols_orient == 0) ||
              (sp2_transp_all_volume_mols_orient ==
               sp_absorb_all_mols_orient)) {
            mcell_error("Conflicting ABSORPTIVE and TRANSPARENT properties are "
                        "simultaneously specified on the same wall using "
                        "ALL_MOLECULES and ALL_VOLUME_MOLECULES through the "
                        "surface classes '%s' and '%s'.",
                        sp->sym->name, sp2->sym->name);
          }
        }
        if (sp2_transp_all_surface_mols) {
          if ((sp2_transp_all_surface_mols_orient == 0) ||
              (sp_absorb_all_mols_orient == 0) ||
              (sp2_transp_all_surface_mols_orient ==
               sp_absorb_all_mols_orient)) {
            mcell_error("Conflicting ABSORPTIVE and TRANSPARENT properties are "
                        "simultaneously specified on the same wall using "
                        "ALL_MOLECULES and ALL_SURFACE_MOLECULES through the "
                        "surface classes '%s' and '%s'.",
                        sp->sym->name, sp2->sym->name);
          }
        }
        if (sp2_refl_all_mols) {
          if ((sp_absorb_all_mols_orient == 0) ||
              (sp2_refl_all_mols_orient == 0) ||
              (sp_absorb_all_mols_orient == sp2_refl_all_mols_orient)) {
            mcell_error("Conflicting ABSORPTIVE and REFLECTIVE properties are "
                        "simultaneously specified on the same wall using "
                        "ALL_MOLECULES through the surface classes '%s' and "
                        "'%s'.",
                        sp->sym->name, sp2->sym->name);
          }
        }
        if (sp2_refl_all_volume_mols) {
          if ((sp_absorb_all_mols_orient == 0) ||
              (sp2_refl_all_volume_mols_orient == 0) ||
              (sp_absorb_all_mols_orient == sp2_refl_all_volume_mols_orient)) {
            mcell_error("Conflicting ABSORPTIVE and REFLECTIVE properties are "
                        "simultaneously specified on the same wall using "
                        "ALL_MOLECULES and ALL_VOLUME_MOLECULES through the "
                        "surface classes '%s' and '%s'.",
                        sp->sym->name, sp2->sym->name);
          }
        }
        if (sp2_refl_all_surface_mols) {
          if ((sp_absorb_all_mols_orient == 0) ||
              (sp2_refl_all_surface_mols_orient == 0) ||
              (sp_absorb_all_mols_orient == sp2_refl_all_surface_mols_orient)) {
            mcell_error("Conflicting ABSORPTIVE and REFLECTIVE properties are "
                        "simultaneously specified on the same wall using "
                        "ALL_MOLECULES and ALL_SURFACE_MOLECULES through the "
                        "surface classes '%s' and '%s'.",
                        sp->sym->name, sp2->sym->name);
          }
        }
      }
      if (sp_refl_all_mols) {
        if (sp2_transp_all_mols) {
          if ((sp2_transp_all_mols_orient == 0) ||
              (sp_refl_all_mols_orient == 0) ||
              (sp2_transp_all_mols_orient == sp_refl_all_mols_orient)) {
            mcell_error("Conflicting REFLECTIVE and TRANSPARENT properties are "
                        "simultaneously specified on the same wall using "
                        "ALL_MOLECULES through the surface classes '%s' and "
                        "'%s'.",
                        sp->sym->name, sp2->sym->name);
          }
        }
        if (sp2_transp_all_volume_mols) {
          if ((sp2_transp_all_volume_mols_orient == 0) ||
              (sp_refl_all_mols_orient == 0) ||
              (sp2_transp_all_volume_mols_orient == sp_refl_all_mols_orient)) {
            mcell_error("Conflicting REFLECTIVE and TRANSPARENT properties are "
                        "simultaneously specified on the same wall using "
                        "ALL_MOLECULES and ALL_VOLUME_MOLECULES through the "
                        "surface classes '%s' and '%s'.",
                        sp->sym->name, sp2->sym->name);
          }
        }
        if (sp2_transp_all_surface_mols) {
          if ((sp2_transp_all_surface_mols_orient == 0) ||
              (sp_refl_all_mols_orient == 0) ||
              (sp2_transp_all_surface_mols_orient == sp_refl_all_mols_orient)) {
            mcell_error("Conflicting REFLECTIVE and TRANSPARENT properties are "
                        "simultaneously specified on the same wall using "
                        "ALL_MOLECULES and ALL_SURFACE_MOLECULES through the "
                        "surface classes '%s' and '%s'.",
                        sp->sym->name, sp2->sym->name);
          }
        }
        if (sp2_absorb_all_mols) {
          if ((sp2_absorb_all_mols_orient == 0) ||
              (sp_refl_all_mols_orient == 0) ||
              (sp2_absorb_all_mols_orient == sp_refl_all_mols_orient)) {
            mcell_error("Conflicting ABSORPTIVE and REFLECTIVE properties are "
                        "simultaneously specified on the same wall using "
                        "ALL_MOLECULES through the surface classes '%s' and "
                        "'%s'.",
                        sp->sym->name, sp2->sym->name);
          }
        }
        if (sp2_absorb_all_volume_mols) {
          if ((sp2_absorb_all_volume_mols_orient == 0) ||
              (sp_refl_all_mols_orient == 0) ||
              (sp2_absorb_all_volume_mols_orient == sp_refl_all_mols_orient)) {
            mcell_error("Conflicting ABSORPTIVE and REFLECTIVE properties are "
                        "simultaneously specified on the same wall using "
                        "ALL_MOLECULES and ALL_VOLUME_MOLECULES through the "
                        "surface classes '%s' and '%s'.",
                        sp->sym->name, sp2->sym->name);
          }
        }
        if (sp2_absorb_all_surface_mols) {
          if ((sp2_absorb_all_surface_mols_orient == 0) ||
              (sp_refl_all_mols_orient == 0) ||
              (sp2_absorb_all_surface_mols_orient == sp_refl_all_mols_orient)) {
            mcell_error("Conflicting ABSORPTIVE and REFLECTIVE properties are "
                        "simultaneously specified on the same wall using "
                        "ALL_MOLECULES and ALL_SURFACE_MOLECULES through the "
                        "surface classes '%s' and '%s'.",
                        sp->sym->name, sp2->sym->name);
          }
        }
      }

      if (sp_transp_all_volume_mols) {
        if (sp2_absorb_all_volume_mols) {
          if ((sp_transp_all_volume_mols_orient == 0) ||
              (sp2_absorb_all_volume_mols_orient == 0) ||
              (sp_transp_all_volume_mols_orient ==
               sp2_absorb_all_volume_mols_orient)) {
            mcell_error("Conflicting TRANSPARENT and ABSORPTIVE properties are "
                        "simultaneously specified on the same wall using "
                        "ALL_VOLUME_MOLECULES through the surface classes '%s' "
                        "and '%s'.",
                        sp->sym->name, sp2->sym->name);
          }
        }
        if (sp2_refl_all_volume_mols) {
          if ((sp_transp_all_volume_mols_orient == 0) ||
              (sp2_refl_all_volume_mols_orient == 0) ||
              (sp_transp_all_volume_mols_orient ==
               sp2_refl_all_volume_mols_orient)) {
            mcell_error("Conflicting TRANSPARENT and REFLECTIVE properties are "
                        "simultaneously specified on the same wall using "
                        "ALL_VOLUME_MOLECULES through the surface classes '%s' "
                        "and '%s'.",
                        sp->sym->name, sp2->sym->name);
          }
        }
      }

      if (sp_absorb_all_volume_mols) {
        if (sp2_transp_all_volume_mols) {
          if ((sp2_transp_all_volume_mols_orient == 0) ||
              (sp_absorb_all_volume_mols_orient == 0) ||
              (sp2_transp_all_volume_mols_orient ==
               sp_absorb_all_volume_mols_orient)) {
            mcell_error("Conflicting ABSORPTIVE and TRANSPARENT properties are "
                        "simultaneously specified on the same wall using "
                        "ALL_VOLUME_MOLECULES through the surface classes '%s' "
                        "and '%s'.",
                        sp->sym->name, sp2->sym->name);
          }
        }
        if (sp2_refl_all_volume_mols) {
          if ((sp_absorb_all_volume_mols_orient == 0) ||
              (sp2_refl_all_volume_mols_orient == 0) ||
              (sp_absorb_all_volume_mols_orient ==
               sp2_refl_all_volume_mols_orient)) {
            mcell_error("Conflicting ABSORPTIVE and REFLECTIVE properties are "
                        "simultaneously specified on the same wall using "
                        "ALL_VOLUME_MOLECULES through the surface classes '%s' "
                        "and '%s'.",
                        sp->sym->name, sp2->sym->name);
          }
        }
      }

      if (sp_refl_all_volume_mols) {
        if (sp2_transp_all_volume_mols) {
          if ((sp2_transp_all_volume_mols_orient == 0) ||
              (sp_refl_all_volume_mols_orient == 0) ||
              (sp2_transp_all_volume_mols_orient ==
               sp_refl_all_volume_mols_orient)) {
            mcell_error("Conflicting REFLECTIVE and TRANSPARENT properties are "
                        "simultaneously specified on the same wall using "
                        "ALL_VOLUME_MOLECULES through the surface classes '%s' "
                        "and '%s'.",
                        sp->sym->name, sp2->sym->name);
          }
        }
        if (sp2_absorb_all_volume_mols) {
          if ((sp2_absorb_all_volume_mols_orient == 0) ||
              (sp_refl_all_volume_mols_orient == 0) ||
              (sp2_absorb_all_volume_mols_orient ==
               sp_refl_all_volume_mols_orient)) {
            mcell_error("Conflicting ABSORPTIVE and REFLECTIVE properties are "
                        "simultaneously specified on the same wall using "
                        "ALL_VOLUME_MOLECULES through the surface classes '%s' "
                        "and '%s'.",
                        sp->sym->name, sp2->sym->name);
          }
        }
      }

      if (sp_transp_all_surface_mols) {
        if (sp2_absorb_all_surface_mols) {
          if ((sp_transp_all_surface_mols_orient == 0) ||
              (sp2_absorb_all_surface_mols_orient == 0) ||
              (sp_transp_all_surface_mols_orient ==
               sp2_absorb_all_surface_mols_orient)) {
            mcell_error("Conflicting TRANSPARENT and ABSORPTIVE properties are "
                        "simultaneously specified on the same wall using "
                        "ALL_SURFACE_MOLECULES through the surface classes "
                        "'%s' and '%s'.",
                        sp->sym->name, sp2->sym->name);
          }
        }
        if (sp2_refl_all_surface_mols) {
          if ((sp_transp_all_surface_mols_orient == 0) ||
              (sp2_refl_all_surface_mols_orient == 0) ||
              (sp_transp_all_surface_mols_orient ==
               sp2_refl_all_surface_mols_orient)) {
            mcell_error("Conflicting TRANSPARENT and REFLECTIVE properties are "
                        "simultaneously specified on the same wall using "
                        "ALL_SURFACE_MOLECULES through the surface classes "
                        "'%s' and '%s'.",
                        sp->sym->name, sp2->sym->name);
          }
        }
      }

      if (sp_absorb_all_surface_mols) {
        if (sp2_transp_all_surface_mols) {
          if ((sp2_transp_all_surface_mols_orient == 0) ||
              (sp_absorb_all_surface_mols_orient == 0) ||
              (sp2_transp_all_surface_mols_orient ==
               sp_absorb_all_surface_mols_orient)) {
            mcell_error("Conflicting ABSORPTIVE and TRANSPARENT properties are "
                        "simultaneously specified on the same wall using "
                        "ALL_SURFACE_MOLECULES through the surface classes "
                        "'%s' and '%s'.",
                        sp->sym->name, sp2->sym->name);
          }
        }
        if (sp2_refl_all_surface_mols) {
          if ((sp_absorb_all_surface_mols_orient == 0) ||
              (sp2_refl_all_surface_mols_orient == 0) ||
              (sp_absorb_all_surface_mols_orient ==
               sp2_refl_all_surface_mols_orient)) {
            mcell_error("Conflicting ABSORPTIVE and REFLECTIVE properties are "
                        "simultaneously specified on the same wall using "
                        "ALL_SURFACE_MOLECULES through the surface classes "
                        "'%s' and '%s'.",
                        sp->sym->name, sp2->sym->name);
          }
        }
      }

      if (sp_refl_all_surface_mols) {
        if (sp2_transp_all_surface_mols) {
          if ((sp2_transp_all_surface_mols_orient == 0) ||
              (sp_refl_all_surface_mols_orient == 0) ||
              (sp2_transp_all_surface_mols_orient ==
               sp_refl_all_surface_mols_orient)) {
            mcell_error("Conflicting REFLECTIVE and TRANSPARENT properties are "
                        "simultaneously specified on the same wall using "
                        "ALL_SURFACE_MOLECULES through the surface classes "
                        "'%s' and '%s'.",
                        sp->sym->name, sp2->sym->name);
          }
        }
        if (sp2_absorb_all_surface_mols) {
          if ((sp2_absorb_all_surface_mols_orient == 0) ||
              (sp_refl_all_surface_mols_orient == 0) ||
              (sp2_absorb_all_surface_mols_orient ==
               sp_refl_all_surface_mols_orient)) {
            mcell_error("Conflicting ABSORPTIVE and REFLECTIVE properties are "
                        "simultaneously specified on the same wall using "
                        "ALL_SURFACE_MOLECULES through the surface classes "
                        "'%s' and '%s'.",
                        sp->sym->name, sp2->sym->name);
          }
        }
      }

      /* Below we will check for the conflicts in combination of
         ALL_MOLECULES (or ALL_VOLUME_MOLECULES, or ALL_SURFACE_MOLECULES)
         with regular molecule */
      if (sp_transp_all_mols) {
        for (no = sp2->absorb_mols; no != NULL; no = no->next) {
          if ((sp_transp_all_mols_orient == no->orient) ||
              (sp_transp_all_mols_orient == 0) || (no->orient == 0)) {
            mcell_error("Conflicting TRANSPARENT and ABSORPTIVE properties are "
                        "simultaneously specified on the same wall using "
                        "ALL_MOLECULES and molecule %s through the surface "
                        "classes '%s' and '%s'.",
                        no->name, sp->sym->name, sp2->sym->name);
          }
        }
        for (no = sp2->refl_mols; no != NULL; no = no->next) {
          if ((sp_transp_all_mols_orient == no->orient) ||
              (sp_transp_all_mols_orient == 0) || (no->orient == 0)) {
            mcell_error("Conflicting TRANSPARENT and REFLECTIVE properties are "
                        "simultaneously specified on the same wall using "
                        "ALL_MOLECULES and molecule %s through the surface "
                        "classes '%s' and '%s'.",
                        no->name, sp->sym->name, sp2->sym->name);
          }
        }
        for (no = sp2->clamp_conc_mols; no != NULL; no = no->next) {
          if ((sp_transp_all_mols_orient == no->orient) ||
              (sp_transp_all_mols_orient == 0) || (no->orient == 0)) {
            mcell_error("Conflicting TRANSPARENT and CLAMP_CONCENTRATION "
                        "properties are simultaneously specified using "
                        "ALL_MOLECULES and molecule %s through the surface "
                        "classes '%s' and '%s'.",
                        no->name, sp->sym->name, sp2->sym->name);
          }
        }
      }

      if (sp_transp_all_volume_mols) {
        for (no = sp2->absorb_mols; no != NULL; no = no->next) {
          spec = get_species_by_name(no->name, n_species, species_list);
          if (spec == NULL) {
            mcell_error("Cannot find species %s in simulation", no->name);
          }
          if (spec->flags & ON_GRID)
            continue;

          if ((sp_transp_all_volume_mols_orient == no->orient) ||
              (sp_transp_all_volume_mols_orient == 0) || (no->orient == 0)) {
            mcell_error("Conflicting TRANSPARENT and ABSORPTIVE properties are "
                        "simultaneously specified on the same wall using "
                        "ALL_VOLUME_MOLECULES and molecule %s through the "
                        "surface classes '%s' and '%s'.",
                        no->name, sp->sym->name, sp2->sym->name);
          }
        }
        for (no = sp2->refl_mols; no != NULL; no = no->next) {
          spec = get_species_by_name(no->name, n_species, species_list);
          if (spec == NULL) {
            mcell_error("Cannot find species %s in simulation", no->name);
          }
          if (spec->flags & ON_GRID)
            continue;

          if ((sp_transp_all_volume_mols_orient == no->orient) ||
              (sp_transp_all_volume_mols_orient == 0) || (no->orient == 0)) {
            mcell_error("Conflicting TRANSPARENT and REFLECTIVE properties "
                        "are simultaneously specified on the same wall using"
                        " ALL_VOLUME_MOLECULES and molecule %s through the "
                        " surface classes '%s' and '%s'.",
                        no->name, sp->sym->name, sp2->sym->name);
          }
        }
        for (no = sp2->clamp_conc_mols; no != NULL; no = no->next) {
          if ((sp_transp_all_volume_mols_orient == no->orient) ||
              (sp_transp_all_volume_mols_orient == 0) || (no->orient == 0)) {
            mcell_error("Conflicting TRANSPARENT and CLAMP_CONCENTRATION "
                        "properties are simultaneously specified using "
                        "ALL_VOLUME_MOLECULES and molecule %s through the "
                        "surface classes '%s' and '%s'.",
                        no->name, sp->sym->name, sp2->sym->name);
          }
        }
      }

      if (sp_transp_all_surface_mols) {
        for (no = sp2->absorb_mols; no != NULL; no = no->next) {
          spec = get_species_by_name(no->name, n_species, species_list);
          if (spec == NULL) {
            mcell_error("Cannot find species %s in simulation", no->name);
          }
          if ((spec->flags & NOT_FREE) == 0)
            continue;

          if ((sp_transp_all_surface_mols_orient == no->orient) ||
              (sp_transp_all_surface_mols_orient == 0) || (no->orient == 0)) {
            mcell_error("Conflicting TRANSPARENT and ABSORPTIVE properties "
                        "are simultaneously specified on the same wall using"
                        "ALL_SURFACE_MOLECULES and molecule %s through the "
                        "surface classes '%s' and '%s'.",
                        no->name, sp->sym->name, sp2->sym->name);
          }
        }
        for (no = sp2->refl_mols; no != NULL; no = no->next) {
          spec = get_species_by_name(no->name, n_species, species_list);
          if (spec == NULL) {
            mcell_error("Cannot find species %s in simulation", no->name);
          }
          if ((spec->flags & NOT_FREE) == 0)
            continue;

          if ((sp_transp_all_surface_mols_orient == no->orient) ||
              (sp_transp_all_surface_mols_orient == 0) || (no->orient == 0)) {
            mcell_error("Conflicting TRANSPARENT and REFLECTIVE properties "
                        "are simultaneously specified on the same wall using"
                        "ALL_SURFACE_MOLECULES and molecule %s through the "
                        "surface classes '%s' and '%s'.",
                        no->name, sp->sym->name, sp2->sym->name);
          }
        }
      }

      if (sp2_transp_all_mols) {
        for (no = sp->absorb_mols; no != NULL; no = no->next) {
          if ((sp2_transp_all_mols_orient == no->orient) ||
              (sp2_transp_all_mols_orient == 0) || (no->orient == 0)) {
            mcell_error("Conflicting TRANSPARENT and ABSORPTIVE properties "
                        "are simultaneously specified on the same wall using"
                        "ALL_MOLECULES and molecule %s through the surface "
                        "classes '%s' and '%s'.",
                        no->name, sp->sym->name, sp2->sym->name);
          }
        }
        for (no = sp->refl_mols; no != NULL; no = no->next) {
          if ((sp2_transp_all_mols_orient == no->orient) ||
              (sp2_transp_all_mols_orient == 0) || (no->orient == 0)) {
            mcell_error("Conflicting TRANSPARENT and REFLECTIVE properties "
                        "are simultaneously specified on the same wall using"
                        "ALL_MOLECULES and molecule %s through the surface "
                        "classes '%s' and '%s'.",
                        no->name, sp->sym->name, sp2->sym->name);
          }
        }
        for (no = sp->clamp_conc_mols; no != NULL; no = no->next) {
          if ((sp2_transp_all_mols_orient == no->orient) ||
              (sp2_transp_all_mols_orient == 0) || (no->orient == 0)) {
            mcell_error("Conflicting TRANSPARENT and CLAMP_CONCENTRATION "
                        "properties are simultaneously specified using "
                        "ALL_MOLECULES and molecule %s through the surface "
                        "classes '%s' and '%s'.",
                        no->name, sp->sym->name, sp2->sym->name);
          }
        }
      }

      if (sp2_transp_all_volume_mols) {
        for (no = sp->absorb_mols; no != NULL; no = no->next) {
          spec = get_species_by_name(no->name, n_species, species_list);
          if (spec == NULL) {
            mcell_error("Cannot find species %s in simulation", no->name);
          }
          if (spec->flags & ON_GRID)
            continue;

          if ((sp2_transp_all_volume_mols_orient == no->orient) ||
              (sp2_transp_all_volume_mols_orient == 0) || (no->orient == 0)) {
            mcell_error("Conflicting TRANSPARENT and ABSORPTIVE properties "
                        "are simultaneously specified on the same wall using"
                        "ALL_VOLUME_MOLECULES and molecule %s through the "
                        "surface classes '%s' and '%s'.",
                        no->name, sp->sym->name, sp2->sym->name);
          }
        }
        for (no = sp->refl_mols; no != NULL; no = no->next) {
          spec = get_species_by_name(no->name, n_species, species_list);
          if (spec == NULL) {
            mcell_error("Cannot find species %s in simulation", no->name);
          }
          if (spec->flags & ON_GRID)
            continue;

          if ((sp2_transp_all_volume_mols_orient == no->orient) ||
              (sp2_transp_all_volume_mols_orient == 0) || (no->orient == 0)) {
            mcell_error("Conflicting TRANSPARENT and REFLECTIVE properties "
                        "are simultaneously specified on the same wall using"
                        "ALL_VOLUME_MOLECULES and molecule %s through the "
                        "surface classes '%s' and '%s'.",
                        no->name, sp->sym->name, sp2->sym->name);
          }
        }
        for (no = sp->clamp_conc_mols; no != NULL; no = no->next) {
          if ((sp2_transp_all_volume_mols_orient == no->orient) ||
              (sp2_transp_all_volume_mols_orient == 0) || (no->orient == 0)) {
            mcell_error("Conflicting TRANSPARENT and CLAMP_CONCENTRATION "
                        "properties are simultaneously specified using "
                        "ALL_VOLUME_MOLECULES and molecule %s through the "
                        "surface classes '%s' and '%s'.",
                        no->name, sp->sym->name, sp2->sym->name);
          }
        }
      }

      if (sp2_transp_all_surface_mols) {
        for (no = sp->absorb_mols; no != NULL; no = no->next) {
          spec = get_species_by_name(no->name, n_species, species_list);
          if (spec == NULL) {
            mcell_error("Cannot find species %s in simulation", no->name);
          }
          if ((spec->flags & NOT_FREE) == 0)
            continue;

          if ((sp2_transp_all_surface_mols_orient == no->orient) ||
              (sp2_transp_all_surface_mols_orient == 0) || (no->orient == 0)) {
            mcell_error("Conflicting TRANSPARENT and ABSORPTIVE properties "
                        "are simultaneously specified on the same wall using"
                        "ALL_SURFACE_MOLECULES and molecule %s through the "
                        "surface classes '%s' and '%s'.",
                        no->name, sp->sym->name, sp2->sym->name);
          }
        }
        for (no = sp->refl_mols; no != NULL; no = no->next) {
          spec = get_species_by_name(no->name, n_species, species_list);
          if (spec == NULL) {
            mcell_error("Cannot find species %s in simulation", no->name);
          }
          if ((spec->flags & NOT_FREE) == 0)
            continue;

          if ((sp2_transp_all_surface_mols_orient == no->orient) ||
              (sp2_transp_all_surface_mols_orient == 0) || (no->orient == 0)) {
            mcell_error("Conflicting TRANSPARENT and REFLECTIVE properties "
                        "are simultaneously specified on the same wall using"
                        " ALL_SURFACE_MOLECULES and molecule %s through the "
                        " surface classes '%s' and '%s'.",
                        no->name, sp->sym->name, sp2->sym->name);
          }
        }
      }

      if (sp_absorb_all_mols) {
        for (no = sp2->transp_mols; no != NULL; no = no->next) {
          if ((sp_absorb_all_mols_orient == no->orient) ||
              (sp_absorb_all_mols_orient == 0) || (no->orient == 0)) {
            mcell_error("Conflicting ABSORPTIVE and TRANSPARENT properties "
                        "are simultaneously specified on the same wall using"
                        "ALL_MOLECULES through the surface classes '%s' and "
                        "'%s'.",
                        sp->sym->name, sp2->sym->name);
          }
        }
        for (no = sp2->refl_mols; no != NULL; no = no->next) {
          if ((sp_absorb_all_mols_orient == no->orient) ||
              (sp_absorb_all_mols_orient == 0) || (no->orient == 0)) {
            mcell_error("Conflicting ABSORPTIVE and REFLECTIVE properties "
                        "are simultaneously specified on the same wall using"
                        " ALL_MOLECULES through the surface classes '%s' and"
                        " '%s'.",
                        sp->sym->name, sp2->sym->name);
          }
        }
        for (no = sp2->clamp_conc_mols; no != NULL; no = no->next) {
          if ((sp_absorb_all_mols_orient == no->orient) ||
              (sp_absorb_all_mols_orient == 0) || (no->orient == 0)) {
            mcell_error("Conflicting ABSORPTIVE and CLAMP_CONCENTRATION "
                        "properties are simultaneously specified using "
                        "ALL_MOLECULES through the surface classes '%s' "
                        "and '%s'.",
                        sp->sym->name, sp2->sym->name);
          }
        }
      }

      if (sp_absorb_all_volume_mols) {
        for (no = sp2->transp_mols; no != NULL; no = no->next) {
          spec = get_species_by_name(no->name, n_species, species_list);
          if (spec == NULL) {
            mcell_error("Cannot find species %s in simulation", no->name);
          }
          if (spec->flags & ON_GRID)
            continue;

          if ((sp_absorb_all_volume_mols_orient == no->orient) ||
              (sp_absorb_all_volume_mols_orient == 0) || (no->orient == 0)) {
            mcell_error("Conflicting ABSORPTIVE and TRANSPARENT properties "
                        "are simultaneously specified on the same wall using"
                        " ALL_VOLUME_MOLECULES through the surface classes "
                        " '%s' and '%s'.",
                        sp->sym->name, sp2->sym->name);
          }
        }
        for (no = sp2->refl_mols; no != NULL; no = no->next) {
          spec = get_species_by_name(no->name, n_species, species_list);
          if (spec == NULL) {
            mcell_error("Cannot find species %s in simulation", no->name);
          }
          if (spec->flags & ON_GRID)
            continue;

          if ((sp_absorb_all_volume_mols_orient == no->orient) ||
              (sp_absorb_all_volume_mols_orient == 0) || (no->orient == 0)) {
            mcell_error("Conflicting ABSORPTIVE and REFLECTIVE properties "
                        "are simultaneously specified on the same wall using"
                        "ALL_VOLUME_MOLECULES through the surface classes "
                        "'%s' and '%s'.",
                        sp->sym->name, sp2->sym->name);
          }
        }
        for (no = sp2->clamp_conc_mols; no != NULL; no = no->next) {
          if ((sp_absorb_all_volume_mols_orient == no->orient) ||
              (sp_absorb_all_volume_mols_orient == 0) || (no->orient == 0)) {
            mcell_error("Conflicting ABSORPTIVE and CLAMP_CONCENTRATION "
                        "properties are simultaneously specified using "
                        "ALL_VOLUME_MOLECULES through the surface classes "
                        "'%s' and '%s'.",
                        sp->sym->name, sp2->sym->name);
          }
        }
      }

      if (sp_absorb_all_surface_mols) {
        for (no = sp2->transp_mols; no != NULL; no = no->next) {
          spec = get_species_by_name(no->name, n_species, species_list);
          if (spec == NULL) {
            mcell_error("Cannot find species %s in simulation", no->name);
          }
          if ((spec->flags & NOT_FREE) == 0)
            continue;

          if ((sp_absorb_all_surface_mols_orient == no->orient) ||
              (sp_absorb_all_surface_mols_orient == 0) || (no->orient == 0)) {
            mcell_error("Conflicting ABSORPTIVE and TRANSPARENT properties "
                        "are simultaneously specified on the same wall using"
                        "ALL_SURFACE_MOLECULES and molecule %s through the "
                        "surface classes '%s' and '%s'.",
                        no->name, sp->sym->name, sp2->sym->name);
          }
        }
        for (no = sp2->refl_mols; no != NULL; no = no->next) {
          spec = get_species_by_name(no->name, n_species, species_list);
          if (spec == NULL) {
            mcell_error("Cannot find species %s in simulation", no->name);
          }
          if ((spec->flags & NOT_FREE) == 0)
            continue;

          if ((sp_absorb_all_surface_mols_orient == no->orient) ||
              (sp_absorb_all_surface_mols_orient == 0) || (no->orient == 0)) {
            mcell_error("Conflicting ABSORPTIVE and REFLECTIVE properties "
                        "are simultaneously specified on the same wall using"
                        "ALL_SURFACE_MOLECULES and molecule %s through the "
                        "surface classes '%s' and '%s'.",
                        no->name, sp->sym->name, sp2->sym->name);
          }
        }
      }

      if (sp2_absorb_all_mols) {
        for (no = sp->transp_mols; no != NULL; no = no->next) {
          if ((sp2_absorb_all_mols_orient == no->orient) ||
              (sp2_absorb_all_mols_orient == 0) || (no->orient == 0)) {
            mcell_error("Conflicting ABSORPTIVE and TRANSPARENT properties "
                        "are simultaneously specified on the same wall using "
                        "ALL_MOLECULES and molecule %s through the surface "
                        "classes '%s' and '%s'.",
                        no->name, sp->sym->name, sp2->sym->name);
          }
        }
        for (no = sp->refl_mols; no != NULL; no = no->next) {
          if ((sp2_absorb_all_mols_orient == no->orient) ||
              (sp2_absorb_all_mols_orient == 0) || (no->orient == 0)) {
            mcell_error("Conflicting ABSORPTIVE and REFLECTIVE properties are"
                        " simultaneously specified on the same wall using "
                        "ALL_MOLECULES and molecule %s through the surface "
                        "classes '%s' and '%s'.",
                        no->name, sp->sym->name, sp2->sym->name);
          }
        }
        for (no = sp->clamp_conc_mols; no != NULL; no = no->next) {
          if ((sp2_absorb_all_mols_orient == no->orient) ||
              (sp2_absorb_all_mols_orient == 0) || (no->orient == 0)) {
            mcell_error("Conflicting ABSORPTIVE and CLAMP_CONCENTRATION "
                        "properties are simultaneously specified using "
                        "ALL_MOLECULES and molecule %s through the surface "
                        "classes '%s' and '%s'.",
                        no->name, sp->sym->name, sp2->sym->name);
          }
        }
      }

      if (sp2_absorb_all_volume_mols) {
        for (no = sp->transp_mols; no != NULL; no = no->next) {
          spec = get_species_by_name(no->name, n_species, species_list);
          if (spec == NULL)
            mcell_error("Cannot find species %s in simulation", no->name);
          if (spec->flags & ON_GRID)
            continue;

          if ((sp2_absorb_all_volume_mols_orient == no->orient) ||
              (sp2_absorb_all_volume_mols_orient == 0) || (no->orient == 0)) {
            mcell_error("Conflicting ABSORPTIVE and TRANSPARENT properties "
                        "are simultaneously specified on the same wall using "
                        "ALL_VOLUME_MOLECULES and molecule %s through the "
                        "surface classes '%s' and '%s'.",
                        no->name, sp->sym->name, sp2->sym->name);
          }
        }
        for (no = sp->refl_mols; no != NULL; no = no->next) {
          spec = get_species_by_name(no->name, n_species, species_list);
          if (spec == NULL)
            mcell_error("Cannot find species %s in simulation", no->name);
          if (spec->flags & ON_GRID)
            continue;

          if ((sp2_absorb_all_volume_mols_orient == no->orient) ||
              (sp2_absorb_all_volume_mols_orient == 0) || (no->orient == 0)) {
            mcell_error("Conflicting ABSORPTIVE and REFLECTIVE properties are"
                        " simultaneously specified on the same wall using "
                        "ALL_VOLUME_MOLECULES and molecule %s through the "
                        "surface classes '%s' and '%s'.",
                        no->name, sp->sym->name, sp2->sym->name);
          }
        }
        for (no = sp->clamp_conc_mols; no != NULL; no = no->next) {
          if ((sp2_absorb_all_volume_mols_orient == no->orient) ||
              (sp2_absorb_all_volume_mols_orient == 0) || (no->orient == 0)) {
            mcell_error("Conflicting ABSORPTIVE and CLAMP_CONCENTRATION "
                        "properties are simultaneously specified using "
                        "ALL_VOLUME_MOLECULES and molecule %s through the "
                        "surface classes '%s' and '%s'.",
                        no->name, sp->sym->name, sp2->sym->name);
          }
        }
      }

      if (sp2_absorb_all_surface_mols) {
        for (no = sp->transp_mols; no != NULL; no = no->next) {
          spec = get_species_by_name(no->name, n_species, species_list);
          if (spec == NULL)
            mcell_error("Cannot find species %s in simulation", no->name);
          if ((spec->flags & NOT_FREE) == 0)
            continue;

          if ((sp2_absorb_all_surface_mols_orient == no->orient) ||
              (sp2_absorb_all_surface_mols_orient == 0) || (no->orient == 0)) {
            mcell_error("Conflicting ABSORPTIVE and TRANSPARENT properties "
                        "are simultaneously specified on the same wall using "
                        "ALL_SURFACE_MOLECULES and molecule %s through the "
                        "surface classes '%s' and '%s'.",
                        no->name, sp->sym->name, sp2->sym->name);
          }
        }
        for (no = sp->refl_mols; no != NULL; no = no->next) {
          spec = get_species_by_name(no->name, n_species, species_list);
          if (spec == NULL)
            mcell_error("Cannot find species %s in simulation", no->name);
          if ((spec->flags & NOT_FREE) == 0)
            continue;

          if ((sp2_absorb_all_surface_mols_orient == no->orient) ||
              (sp2_absorb_all_surface_mols_orient == 0) || (no->orient == 0)) {
            mcell_error("Conflicting ABSORPTIVE and REFLECTIVE properties are"
                        "simultaneously specified on the same wall using "
                        "ALL_SURFACE_MOLECULES and molecule %s through the "
                        "surface classes '%s' and '%s'.",
                        no->name, sp->sym->name, sp2->sym->name);
          }
        }
      }

      if (sp_refl_all_mols) {
        for (no = sp2->transp_mols; no != NULL; no = no->next) {
          if ((sp_refl_all_mols_orient == no->orient) ||
              (sp_refl_all_mols_orient == 0) || (no->orient == 0)) {
            mcell_error("Conflicting REFLECTIVE and TRANSPARENT properties "
                        "are simultaneously specified on the same wall using "
                        "ALL_MOLECULES and molecule %s through the surface "
                        "classes '%s' and '%s'.",
                        no->name, sp->sym->name, sp2->sym->name);
          }
        }
        for (no = sp2->absorb_mols; no != NULL; no = no->next) {
          if ((sp_refl_all_mols_orient == no->orient) ||
              (sp_refl_all_mols_orient == 0) || (no->orient == 0)) {
            mcell_error("Conflicting REFLECTIVE and ABSORPTIVE properties are"
                        "simultaneously specified on the same wall using "
                        "ALL_MOLECULES and molecule %s through the surface "
                        "classes '%s' and '%s'.",
                        no->name, sp->sym->name, sp2->sym->name);
          }
        }
        for (no = sp2->clamp_conc_mols; no != NULL; no = no->next) {
          if ((sp_refl_all_mols_orient == no->orient) ||
              (sp_refl_all_mols_orient == 0) || (no->orient == 0)) {
            mcell_error("Conflicting REFLECTIVE and CLAMP_CONCENTRATION "
                        "properties are simultaneously specified using "
                        "ALL_MOLECULES and molecule %s through the surface "
                        "classes '%s' and '%s'.",
                        no->name, sp->sym->name, sp2->sym->name);
          }
        }
      }

      if (sp_refl_all_volume_mols) {
        for (no = sp2->transp_mols; no != NULL; no = no->next) {
          spec = get_species_by_name(no->name, n_species, species_list);
          if (spec == NULL)
            mcell_error("Cannot find species %s in simulation", no->name);
          if (spec->flags & ON_GRID)
            continue;

          if ((sp_refl_all_volume_mols_orient == no->orient) ||
              (sp_refl_all_volume_mols_orient == 0) || (no->orient == 0)) {
            mcell_error("Conflicting REFLECTIVE and TRANSPARENT properties "
                        "are simultaneously specified on the same wall using "
                        "ALL_VOLUME_MOLECULES and molecule %s through the "
                        "surface classes '%s' and '%s'.",
                        no->name, sp->sym->name, sp2->sym->name);
          }
        }
        for (no = sp2->absorb_mols; no != NULL; no = no->next) {
          spec = get_species_by_name(no->name, n_species, species_list);
          if (spec == NULL)
            mcell_error("Cannot find species %s in simulation", no->name);
          if (spec->flags & ON_GRID)
            continue;

          if ((sp_refl_all_volume_mols_orient == no->orient) ||
              (sp_refl_all_volume_mols_orient == 0) || (no->orient == 0)) {
            mcell_error("Conflicting REFLECTIVE and ABSORPTIVE properties are"
                        " simultaneously specified on the same wall using "
                        "ALL_VOLUME_MOLECULES and molecule %s through the "
                        "surface classes '%s' and '%s'.",
                        no->name, sp->sym->name, sp2->sym->name);
          }
        }
        for (no = sp2->clamp_conc_mols; no != NULL; no = no->next) {
          if ((sp_refl_all_volume_mols_orient == no->orient) ||
              (sp_refl_all_volume_mols_orient == 0) || (no->orient == 0)) {
            mcell_error("Conflicting REFLECTIVE and CLAMP_CONCENTRATION "
                        "properties are simultaneously specified using "
                        "ALL_VOLUME_MOLECULES and molecule %s through the "
                        "surface classes '%s' and '%s'.",
                        no->name, sp->sym->name, sp2->sym->name);
          }
        }
      }

      if (sp_refl_all_surface_mols) {
        for (no = sp2->transp_mols; no != NULL; no = no->next) {
          spec = get_species_by_name(no->name, n_species, species_list);
          if (spec == NULL)
            mcell_error("Cannot find species %s in simulation", no->name);
          if ((spec->flags & NOT_FREE) == 0)
            continue;

          if ((sp_refl_all_surface_mols_orient == no->orient) ||
              (sp_refl_all_surface_mols_orient == 0) || (no->orient == 0)) {
            mcell_error("Conflicting REFLECTIVE and TRANSPARENT properties "
                        "are simultaneously specified on the same wall using "
                        "ALL_SURFACE_MOLECULES and molecule %s through the "
                        "surface classes '%s' and '%s'.",
                        no->name, sp->sym->name, sp2->sym->name);
          }
        }
        for (no = sp2->absorb_mols; no != NULL; no = no->next) {
          spec = get_species_by_name(no->name, n_species, species_list);
          if (spec == NULL)
            mcell_error("Cannot find species %s in simulation", no->name);
          if ((spec->flags & NOT_FREE) == 0)
            continue;

          if ((sp_refl_all_surface_mols_orient == no->orient) ||
              (sp_refl_all_surface_mols_orient == 0) || (no->orient == 0)) {
            mcell_error("Conflicting REFLECTIVE and ABSORPTIVE properties are"
                        "simultaneously specified on the same wall using "
                        "ALL_SURFACE_MOLECULES and molecule %s through the "
                        "surface classes '%s' and '%s'.",
                        no->name, sp->sym->name, sp2->sym->name);
          }
        }
      }

      if (sp2_refl_all_mols) {
        for (no = sp->transp_mols; no != NULL; no = no->next) {
          if ((sp2_refl_all_mols_orient == no->orient) ||
              (sp2_refl_all_mols_orient == 0) || (no->orient == 0)) {
            mcell_error("Conflicting REFLECTIVE and TRANSPARENT properties "
                        "are simultaneously specified on the same wall using "
                        "ALL_MOLECULES and molecule %s through the surface "
                        "classes '%s' and '%s'.",
                        no->name, sp->sym->name, sp2->sym->name);
          }
        }
        for (no = sp->absorb_mols; no != NULL; no = no->next) {
          if ((sp2_refl_all_mols_orient == no->orient) ||
              (sp2_refl_all_mols_orient == 0) || (no->orient == 0)) {
            mcell_error("Conflicting REFLECTIVE and ABSORPTIVE properties are"
                        " simultaneously specified on the same wall using "
                        "ALL_MOLECULES and molecule %s through the surface "
                        "classes '%s' and '%s'.",
                        no->name, sp->sym->name, sp2->sym->name);
          }
        }
        for (no = sp->clamp_conc_mols; no != NULL; no = no->next) {
          if ((sp2_refl_all_mols_orient == no->orient) ||
              (sp2_refl_all_mols_orient == 0) || (no->orient == 0)) {
            mcell_error("Conflicting REFLECTIVE and CLAMP_CONCENTRATION "
                        "properties are simultaneously specified using "
                        "ALL_MOLECULES and molecule %s through the surface "
                        "classes '%s' and '%s'.",
                        no->name, sp->sym->name, sp2->sym->name);
          }
        }
      }

      if (sp2_refl_all_volume_mols) {
        for (no = sp->transp_mols; no != NULL; no = no->next) {
          spec = get_species_by_name(no->name, n_species, species_list);
          if (spec == NULL)
            mcell_error("Cannot find species %s in simulation", no->name);
          if (spec->flags & ON_GRID)
            continue;

          if ((sp2_refl_all_volume_mols_orient == no->orient) ||
              (sp2_refl_all_volume_mols_orient == 0) || (no->orient == 0)) {
            mcell_error("Conflicting REFLECTIVE and TRANSPARENT properties "
                        "are simultaneously specified on the same wall using"
                        " ALL_VOLUME_MOLECULES and molecule %s through the "
                        "surface classes '%s' and '%s'.",
                        no->name, sp->sym->name, sp2->sym->name);
          }
        }
        for (no = sp->absorb_mols; no != NULL; no = no->next) {
          spec = get_species_by_name(no->name, n_species, species_list);
          if (spec == NULL)
            mcell_error("Cannot find species %s in simulation", no->name);
          if (spec->flags & ON_GRID)
            continue;

          if ((sp2_refl_all_volume_mols_orient == no->orient) ||
              (sp2_refl_all_volume_mols_orient == 0) || (no->orient == 0)) {
            mcell_error("Conflicting REFLECTIVE and ABSORPTIVE properties are"
                        "simultaneously specified on the same wall using "
                        "ALL_VOLUME_MOLECULES and molecule %s through the "
                        "surface classes '%s' and '%s'.",
                        no->name, sp->sym->name, sp2->sym->name);
          }
        }
        for (no = sp->clamp_conc_mols; no != NULL; no = no->next) {
          if ((sp2_refl_all_volume_mols_orient == no->orient) ||
              (sp2_refl_all_volume_mols_orient == 0) || (no->orient == 0)) {
            mcell_error("Conflicting REFLECTIVE and CLAMP_CONCENTRATION "
                        "properties are simultaneously specified using "
                        "ALL_VOLUME_MOLECULES and molecule %s through the "
                        "surface classes '%s' and '%s'.",
                        no->name, sp->sym->name, sp2->sym->name);
          }
        }
      }

      if (sp2_refl_all_surface_mols) {
        for (no = sp->transp_mols; no != NULL; no = no->next) {
          spec = get_species_by_name(no->name, n_species, species_list);
          if (spec == NULL)
            mcell_error("Cannot find species %s in simulation", no->name);
          if ((spec->flags & NOT_FREE) == 0)
            continue;

          if ((sp2_refl_all_surface_mols_orient == no->orient) ||
              (sp2_refl_all_surface_mols_orient == 0) || (no->orient == 0)) {
            mcell_error("Conflicting REFLECTIVE and TRANSPARENT properties "
                        "are simultaneously specified on the same wall using "
                        "ALL_SURFACE_MOLECULES and molecule %s through the "
                        "surface classes '%s' and '%s'.",
                        no->name, sp->sym->name, sp2->sym->name);
          }
        }
        for (no = sp->absorb_mols; no != NULL; no = no->next) {
          spec = get_species_by_name(no->name, n_species, species_list);
          if (spec == NULL)
            mcell_error("Cannot find species %s in simulation", no->name);
          if ((spec->flags & NOT_FREE) == 0)
            continue;

          if ((sp2_refl_all_surface_mols_orient == no->orient) ||
              (sp2_refl_all_surface_mols_orient == 0) || (no->orient == 0)) {
            mcell_error("Conflicting REFLECTIVE and ABSORPTIVE properties are"
                        " simultaneously specified on the same wall using "
                        "ALL_SURFACE_MOLECULES and molecule %s through the "
                        "surface classes '%s' and '%s'.",
                        no->name, sp->sym->name, sp2->sym->name);
          }
        }
      }
    }
  }
}

/***************************************************************************
get_species_by_name:
  In: name of species
  Out: Species object or NULL if we can't find one.
***************************************************************************/
struct species *get_species_by_name(char *name, int n_species,
                                    struct species **species_list) {
  for (int i = 0; i < n_species; i++) {
    struct species *sp = species_list[i];

    if (strcmp(sp->sym->name, name) == 0)
      return sp;
  }

  return NULL;
}

/***************************************************************************
create_name_lists_of_volume_and_surface_mols:
  In: pointer to the empty linked lists.
  Out: none. Two linked lists of volume molecules names and surface molecules
       names are created.
***************************************************************************/
void create_name_lists_of_volume_and_surface_mols(
    struct volume *world, struct name_list **vol_species_name_list,
    struct name_list **surf_species_name_list) {
  struct species *spec;
  struct name_list *nl, *nl_vol_head = NULL, *nl_surf_head = NULL;
  int i;

  /* create name lists of volume and surface species */
  for (i = 0; i < world->n_species; i++) {
    spec = world->species_list[i];
    if (spec == NULL)
      mcell_internal_error("Cannot find molecule name");
    if (spec == world->all_mols)
      continue;
    if ((spec == world->all_volume_mols) || (spec == world->all_surface_mols))
      continue;
    if (spec->flags & IS_SURFACE)
      continue;

    if (spec->flags & ON_GRID) {
      nl = CHECKED_MALLOC_STRUCT(struct name_list, "name_list");
      nl->name = CHECKED_STRDUP(spec->sym->name, "species name");
      nl->prev = NULL; /* we will use only FORWARD feature */

      if (nl_surf_head == NULL) {
        nl->next = NULL;
        nl_surf_head = nl;
      } else {
        nl->next = nl_surf_head;
        nl_surf_head = nl;
      }
    } else {
      nl = CHECKED_MALLOC_STRUCT(struct name_list, "name_list");
      nl->name = CHECKED_STRDUP(spec->sym->name, "species name");
      nl->prev = NULL; /* we will use only FORWARD feature */

      if (nl_vol_head == NULL) {
        nl->next = NULL;
        nl_vol_head = nl;
      } else {
        nl->next = nl_vol_head;
        nl_vol_head = nl;
      }
    }
  }

  if (nl_vol_head != NULL)
    *vol_species_name_list = nl_vol_head;
  if (nl_surf_head != NULL)
    *surf_species_name_list = nl_surf_head;
}

/***************************************************************************
remove_molecules_name_list:
  In: none
  Out: none. Global linked list of molecules names is memory
       deallocated.
***************************************************************************/
void remove_molecules_name_list(struct name_list **nlist) {
  struct name_list *nnext, *nl_head;

  nl_head = *nlist;

  /* remove name list */
  while (nl_head != NULL) {
    nnext = nl_head->next;
    if (nl_head->name != NULL)
      free(nl_head->name);
    free(nl_head);
    nl_head = nnext;
  }
  nl_head = NULL;

  *nlist = NULL;
}

/*****************************************************************
check_for_overlapped_walls:
  In: rng: random number generator
      n_subvols: number of subvolumes
      subvol: a subvolume
  Out: 0 if no errors, the world geometry is successfully checked for
       overlapped walls.
       1 if there are any overlapped walls.
******************************************************************/
int check_for_overlapped_walls(
    struct rng_state *rng, int n_subvols, struct subvolume *subvol) {

  /* pick up a random vector */
  struct vector3 rand_vector;
  rand_vector.x = rng_dbl(rng);
  rand_vector.y = rng_dbl(rng);
  rand_vector.z = rng_dbl(rng);

  for (int i = 0; i < n_subvols; i++) {
    struct subvolume *sv = &(subvol[i]);
    struct wall_aux_list *head = NULL;

    for (struct wall_list *wlp = sv->wall_head; wlp != NULL; wlp = wlp->next) {
      double d_prod = dot_prod(&rand_vector, &(wlp->this_wall->normal));
      /* we want to place walls with opposite normals into
         neighboring positions in the sorted linked list */
      if (d_prod < 0)
        d_prod = -d_prod;

      struct wall_aux_list *newNode;
      newNode = CHECKED_MALLOC_STRUCT(struct wall_aux_list, "wall_aux_list");
      newNode->this_wall = wlp->this_wall;
      newNode->d_prod = d_prod;

      sorted_insert_wall_aux_list(&head, newNode);
    }

    for (struct wall_aux_list *curr = head; curr != NULL; curr = curr->next) {
      struct wall *w1 = curr->this_wall;

      struct wall_aux_list *next_curr = curr->next;
      while ((next_curr != NULL) &&
             (!distinguishable(curr->d_prod, next_curr->d_prod, EPS_C))) {
        /* there may be several walls with the same (or mirror)
           oriented normals */
        struct wall *w2 = next_curr->this_wall;

        if (are_walls_coplanar(w1, w2, MESH_DISTINCTIVE)) {
          if ((are_walls_coincident(w1, w2, MESH_DISTINCTIVE) ||
               coplanar_tri_overlap(w1, w2))) {
            mcell_error(
                "walls are overlapped: wall %d from '%s' and wall "
                "%d from '%s'.",
                w1->side, w1->parent_object->sym->name, w2->side,
                w2->parent_object->sym->name);
          }
        }
        next_curr = next_curr->next;
      }
    }
    /* free memory */
    if (head != NULL)
      delete_wall_aux_list(head);
  }

  return 0;
}


/* this function checks if obj is contained in the linked list of objects
 * which have the passed-in concentration clamp */
struct ccn_clamp_data* find_clamped_object_in_list(struct ccn_clamp_data *ccd,
  struct object *obj) {
  struct ccn_clamp_data *c = ccd;
  while (c->next_obj != NULL) {
    if (c->next_obj->objp == obj) {
      return c->next_obj;
    }
    c = c->next_obj;
  }
  return NULL;
}
