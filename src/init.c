#include <assert.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <float.h>
#include <time.h>
#include <netdb.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <sys/time.h>
#include <sys/resource.h>

#ifdef KELP
#include <kelp.h>
#endif

#include "version_info.h"
#include "logging.h"
#include "rng.h"
#include "mcell_structs.h"
#include "strfunc.h"
#include "vector.h"
#include "sym_table.h"
#include "count_util.h"
#include "vol_util.h"
#include "wall_util.h"
#include "grid_util.h"
#include "viz_output.h"
#include "react.h"
#include "react_output.h"
#include "util.h"
#include "chkpt.h"
#include "mdlparse_util.h"
#include "init.h"
#include "mdlparse_aux.h"

extern struct volume *world;

/* Initialize the surface macromolecules on a given object */
static int init_complex_effectors(struct object *objp, struct region_list *head);

/* Initialize the visualization output (frame_data_lists). */
static int init_viz_output(void);

static int compute_bb(struct object *objp, double (*im)[4]);
static int compute_bb_release_site(struct object *objp, double (*im)[4]);
static int compute_bb_polygon_object(struct object *objp, double (*im)[4]);

#define MICROSEC_PER_YEAR 365.25*86400.0*1e6

/* Sets default notification values */
int init_notifications(void)
{
  world->notify = CHECKED_MALLOC_STRUCT(struct notifications,
                                        "notification states");
 
  /* Notifications */
  if (world->quiet_flag)
  {
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
  }
  else
  {
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
  
  if (world->log_freq != 0  &&  world->log_freq != ULONG_MAX) /* User set this */
  {
    world->notify->custom_iteration_value = world->log_freq;
  }

  return 0;
}

/*
 * Initialize the volume data output.  This will create a scheduler for the
 * volume data, and add volume data output items to the scheduler.
 */
static void init_volume_data_output(struct volume *wrld)
{
  struct volume_output_item *vo, *vonext;

  wrld->volume_output_scheduler = create_scheduler(1.0, 100.0, 100, wrld->current_start_real_time / wrld->time_unit);
  if (wrld->volume_output_scheduler == NULL)
    mcell_allocfailed("Failed to create scheduler for volume output data.");

  double r_time_unit = 1.0 / world->time_unit;
  for (vo = wrld->volume_output_head; vo != NULL; vo = vonext)
  {
    vonext = vo->next;  /* schedule_add overwrites 'next' */

    if (vo->timer_type==OUTPUT_BY_STEP)
    {
      if (world->chkpt_seq_num == 1) vo->t=0.0;
      else
      {
        /* Get step time in internal units, find next scheduled output time */
        double f = vo->step_time * r_time_unit;
        vo->t = f * ceil(wrld->volume_output_scheduler->now / f);
      }
    }
    else if (vo->num_times > 0)
    {
      /* Set time scaling factor depending on output type */
      double time_scale = 0.0;
      if (vo->timer_type==OUTPUT_BY_ITERATION_LIST) time_scale = 1.0;
      else time_scale = r_time_unit;

      /* Find the time of next output */
      if (world->chkpt_seq_num == 1) {
        vo->next_time = vo->times;
        vo->t = time_scale * *vo->next_time;
      }
      else /* Scan forward to find first output after checkpoint time */
      {
        int idx = bisect_high(vo->times, vo->num_times, world->volume_output_scheduler->now / time_scale);

        /* If we've already passed the last time for this one, skip it! */
        if (idx < 0 || idx >= vo->num_times)
          continue;
        if (world->volume_output_scheduler->now / time_scale > vo->times[idx])
          continue;

        vo->t = vo->times[idx] * time_scale;
        vo->next_time = vo->times + idx;
      }

      /* Advance the next_time pointer */
      ++ vo->next_time;
    }

    if (schedule_add(world->volume_output_scheduler, vo))
      mcell_allocfailed("Failed to add item to schedule for volume output.");
  }
}

/**
 * Initializes the parameters required for the simulation (duh!).
 * This driver function sets up everything needed to get the simulation
 * up and running. It does the following (in that approximate order):
 * 	- Initializes variables to nice values
 * 	- Initialize all the compute nodes (init_nodes())
 * 	- Creates all the data structures
 * 	- Parse the MDL input file (mdlparse_init())
 * 	- Setup the geometry (init_geom())
 * \todo Need more info here.
 */
int init_sim(void)
{
  struct sym_table *gp;
  struct output_block *obp,*obpn;
  struct output_set *set;
  int i;
  double f;
  int reactants_3D_present = 0; /* flag to check whether there are 3D reactants
                             (participants in the reactions
                              between 3D molecules) in the simulation */

#ifdef KELP
  if (world->procnum == 0) {
#endif
    if (world->notify->progress_report != NOTIFY_NONE)
      mcell_log("MCell initializing simulation...");
#ifdef KELP
  }
#endif


  /* Initialize variables to reasonably safe values */

  world->curr_file=world->mdl_infile_name;

  /* by Erhan Gokcay 5/3/2002 ========================== */
  /* We can not initialize chkpt_infile anymore. It is set in mcell.c */
  /*  chkpt_infile=NULL; */
  /* =================================================== */

  world->chkpt_iterations=0;
  world->chkpt_seq_num=0;

  /*world->chkpt_init=1; */  /* set in the main() */
  world->chkpt_flag=0;
  world->viz_blocks=NULL;
  world->ray_voxel_tests=0;
  world->ray_polygon_tests=0;
  world->ray_polygon_colls=0;
  world->mol_mol_colls=0;
  world->mol_grid_colls=0;
  world->grid_grid_colls=0;
  world->mol_wall_colls=0;
  world->mol_mol_mol_colls=0;
  world->mol_mol_grid_colls=0;
  world->mol_grid_grid_colls=0;
  world->grid_grid_grid_colls=0;
  world->chkpt_elapsed_real_time=0;
  world->chkpt_elapsed_real_time_start=0;
  world->chkpt_byte_order_mismatch = 0;
  world->it_time=0;
  world->elapsed_time=0;
  world->time_unit=0;
  world->time_step_max=0;
  world->start_time=0;
  world->current_real_time=0;
  world->current_start_real_time=0;
  world->grid_density=10000;
  world->r_length_unit=sqrt(world->grid_density);
  world->length_unit=1.0/world->r_length_unit;
  world->rx_radius_3d = 0;
  world->radial_directions=16384;
  world->radial_subdivisions=1024;
  world->fully_random=0;
  world->num_directions=world->radial_directions;
  world->r_step=NULL;
  world->r_step_surface=NULL;
  world->r_step_release=NULL;
  world->d_step=NULL;
  world->dissociation_index = DISSOCIATION_MAX;
  world->place_waypoints_flag=0;
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
  
  world->use_expanded_list=1;
  world->randomize_gmol_pos=1;
  world->vacancy_search_dist2=0;
  world->surface_reversibility=0;
  world->volume_reversibility=0;
  world->n_reactions = 0;

  world->mol_mol_reaction_flag = 0;
  world->mol_grid_reaction_flag = 0;
  world->grid_grid_reaction_flag = 0;
  world->mol_wall_reaction_flag = 0;
  world->mol_mol_mol_reaction_flag = 0;
  world->mol_mol_grid_reaction_flag = 0;
  world->mol_grid_grid_reaction_flag = 0;
  world->grid_grid_grid_reaction_flag = 0;

  world->mcell_version = mcell_version;
  
  world->clamp_list = NULL;

  world->rng = CHECKED_MALLOC_STRUCT(struct rng_state,
                                     "random number generator state");
  if (world->seed_seq < 1 || world->seed_seq > INT_MAX)
    mcell_error("Random sequence number must be in the range 1 to 2^31-1 [2147483647]");
  rng_init(world->rng,world->seed_seq);
  if (world->notify->progress_report != NOTIFY_NONE)
    mcell_log("MCell[%d]: random sequence %d",world->procnum,world->seed_seq);

  world->count_hashmask = COUNT_HASHMASK;
  world->count_hash = CHECKED_MALLOC_ARRAY(struct counter*,
                                           (world->count_hashmask+1),
                                           "counter hash table");
  for (i=0;i<=world->count_hashmask;i++) world->count_hash[i] = NULL;
  
  world->oexpr_mem = create_mem_named(sizeof(struct output_expression),128,"output expression");
  if (world->oexpr_mem==NULL)
    mcell_allocfailed("Failed to create memory pool for reaction data output expressions.");
  world->outp_request_mem = create_mem_named(sizeof(struct output_request),64,"output request");
  if (world->outp_request_mem==NULL)
    mcell_allocfailed("Failed to create memory pool for reaction data output commands.");
  world->counter_mem = create_mem_named(sizeof(struct counter),32,"counter");
  if (world->counter_mem==NULL)
    mcell_allocfailed("Failed to create memory pool for reaction and molecule counts.");
  world->trig_request_mem = create_mem_named(sizeof(struct trigger_request),32,"trigger request");
  if (world->trig_request_mem==NULL)
    mcell_allocfailed("Failed to create memory pool for reaction and molecule output triggers.");
  world->magic_mem = create_mem_named(sizeof(struct magic_list),1024,"reaction-triggered release");
  if (world->magic_mem==NULL)
    mcell_allocfailed("Failed to create memory pool for reaction-triggered release lists.");

  if ((world->fstream_sym_table = init_symtab(1024)) == NULL)
    mcell_allocfailed("Failed to initialize symbol table for file streams.");
  if ((world->mol_sym_table = init_symtab(1024)) == NULL)
    mcell_allocfailed("Failed to initialize symbol table for molecules.");
  if ((world->rxn_sym_table = init_symtab(1024)) == NULL)
    mcell_allocfailed("Failed to initialize symbol table for reactions.");
  if ((world->obj_sym_table = init_symtab(1024)) == NULL)
    mcell_allocfailed("Failed to initialize symbol table for objects.");
  if ((world->reg_sym_table = init_symtab(1024)) == NULL)
    mcell_allocfailed("Failed to initialize symbol table for regions.");
  if ((world->rpat_sym_table = init_symtab(1024)) == NULL)
    mcell_allocfailed("Failed to initialize symbol table for release patterns.");
  if ((world->rxpn_sym_table = init_symtab(1024)) == NULL)
    mcell_allocfailed("Failed to initialize symbol table for reaction pathways.");

  if ((gp = store_sym("WORLD_OBJ", OBJ, world->obj_sym_table, NULL)) == NULL)
    mcell_allocfailed("Failed to store the world root object in the object symbol table.");
  world->root_object = (struct object *) gp->value;
  world->root_object->object_type = META_OBJ;
  world->root_object->last_name = CHECKED_STRDUP("", NULL);

  if ((gp = store_sym("WORLD_INSTANCE", OBJ, world->obj_sym_table, NULL)) == NULL)
    mcell_allocfailed("Failed to store the world root instance in the object symbol table.");
  world->root_instance = (struct object *)gp->value;
  world->root_instance->object_type = META_OBJ;
  world->root_instance->last_name = CHECKED_STRDUP("", NULL);

  if ((gp = store_sym("DEFAULT_RELEASE_PATTERN", RPAT, world->rpat_sym_table, NULL)) == NULL)
    mcell_allocfailed("Failed to store the default release pattern in the release patterns symbol table.");
  world->default_release_pattern = (struct release_pattern *) gp->value;
  world->default_release_pattern->delay = 0;
  world->default_release_pattern->release_interval = FOREVER;
  world->default_release_pattern->train_interval = FOREVER;
  world->default_release_pattern->train_duration = FOREVER;
  world->default_release_pattern->number_of_trains = 1;
   
  if ((gp = store_sym("GENERIC_MOLECULE", MOL, world->mol_sym_table, NULL)) == NULL)
    mcell_allocfailed("Failed to store the generic molecule in the molecule symbol table.");
  world->g_mol = (struct species *) gp->value;

  if ((gp = store_sym("GENERIC_SURFACE", MOL, world->mol_sym_table, NULL)) == NULL)
    mcell_allocfailed("Failed to store the generic surface class in the molecule symbol table.");
  world->g_surf = (struct species *) gp->value;
  world->g_surf->flags = IS_SURFACE;

  world->volume_output_head = NULL;
  world->output_block_head=NULL;
  world->output_request_head=NULL;

  world->releaser = create_scheduler(1.0, 100.0, 100, 0.0);
  if (world->releaser == NULL)
    mcell_allocfailed("Failed to create release scheduler.");

  /* Parse the MDL file: */
  no_printf("Node %d parsing MDL file %s\n",world->procnum,world->mdl_infile_name);
  if (mdlparse_init(world)) {
    return(1);
  }
  no_printf("Done parsing MDL file: %s\n",world->mdl_infile_name);
  install_emergency_output_hooks();
  emergency_output_hook_enabled = 0;

  /* we do not want to count collisions if the policy is not to print */
  if(world->notify->final_summary == NOTIFY_NONE) world->notify->molecule_collision_report = NOTIFY_NONE;

  if (world->iterations == INT_MIN)
    mcell_error("Total number of iterations is not specified either through the ITERATIONS keyword or through the command line option '-iterations'.");

  /* Set up the array of species */
  if (init_species())
    mcell_error("Unknown error while initializing species table.");
  no_printf("Done setting up species.\n");

  for(i = 0; i < world->n_species; i++)
  {
    struct species *sp = world->species_list[i];

    if((sp->flags & IS_SURFACE) && (world->notify->reaction_probabilities==NOTIFY_FULL)){
       publish_special_reactions_report(sp);
    }
  }
  

 /* If there are no 3D molecules-reactants in the simulation
    set up the"use_expanded_list" flag to zero. */
  for(i = 0; i < world->n_species; i++)
  {
    struct species *sp = world->species_list[i];
    if ((sp == world->g_mol)  ||  (sp == world->g_surf)) continue;
    if(sp->flags & ON_GRID) continue;
    if(sp->flags & IS_SURFACE) continue;
    
    if ((sp->flags & (CAN_MOLMOL|CAN_MOLMOLMOL)) != 0)
    {
      reactants_3D_present = 1;
      break;
    }
  }

  if(reactants_3D_present == 0){
	world->use_expanded_list = 0;
  }

/* Instantiation Pass #1: Initialize the geometry */
  if (init_geom())
    mcell_internal_error("Unknown error while initializing world geometry.");
  no_printf("Done setting up geometry.\n");
  
/* Instantiation Pass #2: Partition geometry */
  if (init_partitions())
    mcell_internal_error("Unknown error while initializing partitions.");
  if (distribute_world())
    mcell_internal_error("Unknown error while distributing geometry among partitions.");
  if (sharpen_world())
    mcell_internal_error("Unknown error while adding edges to geometry.");


/* Instantiation Pass #3: Initialize regions */
  if (prepare_counters())
    mcell_internal_error("Unknown error while preparing counters for reaction data output.");
  if (init_regions())
    mcell_internal_error("Unknown error while initializing object regions.");
  if (check_counter_geometry())
    mcell_internal_error("Unknown error while validating geometry of counting regions.");
  
  if (world->place_waypoints_flag)
  {
    if (place_waypoints())
      mcell_internal_error("Unknown error while placing waypoints.");
  }
  
  if (init_effectors())
    mcell_internal_error("Unknown error while placing effectors on regions.");
  
  if (init_releases())
    mcell_internal_error("Unknown error while initializing release sites.");

  if (world->chkpt_infile)
  {
    FILE *chkpt_infs = NULL;
    if ((chkpt_infs = fopen(world->chkpt_infile,"rb")) == NULL)
      world->chkpt_seq_num = 1;
    else
    {
      mcell_log("Reading from checkpoint file '%s'.", world->chkpt_infile);
      if (read_chkpt(chkpt_infs))
        mcell_error("Failed to read checkpoint file '%s'.", world->chkpt_infile);
      fclose(chkpt_infs);
    }
  }
  else {
    world->chkpt_seq_num=1;
  }

  /* Initialize the frame data for the visualization and reaction output. */
  if (init_viz_output())
    mcell_internal_error("Unknown error while initializing VIZ output.");

  /* Initialize the volume output */
  init_volume_data_output(world);

  world->count_scheduler = create_scheduler(1.0,100.0,100,world->start_time);
  if (world->count_scheduler == NULL)
    mcell_allocfailed("Failed to create scheduler for reaction data output.");

  /* Schedule the reaction data output events */
  obp = world->output_block_head;
  while(obp != NULL)
  {
    obpn = obp->next; /* Save this--will be lost when we schedule obp */
    
    if (obp->timer_type==OUTPUT_BY_STEP)
    {
      if (world->chkpt_seq_num==1) obp->t=0.0;
      else
      {
        f = obp->step_time/world->time_unit; /* Step time (internal units) */
        obp->t = f*ceil(world->count_scheduler->now / f) + f; /* Round up */
      }      
    }
    else if (obp->time_now==NULL) /* When would this be non-NULL?? */
    {
      /* Set time scaling factor depending on output type */
      if (obp->timer_type==OUTPUT_BY_ITERATION_LIST) f=1.0;
      else f=1.0/world->time_unit;
      
      /* Find the time of next output */
      if (world->chkpt_seq_num == 1)
      {
        obp->time_now = obp->time_list_head;
        obp->t = f*obp->time_now->value;
      }
      else /* Scan forward to find first output after checkpoint time */
      {
        for (obp->time_now=obp->time_list_head ; obp->time_now!=NULL ; obp->time_now=obp->time_now->next)
        {
          if(obp->timer_type == OUTPUT_BY_ITERATION_LIST){
          obp->t=f*obp->time_now->value;
          if (!(obp->t < world->iterations+1 && obp->t <= world->count_scheduler->now)) break;
          }else if(obp->timer_type == OUTPUT_BY_TIME_LIST){
            if (obp->time_now->value > world->current_start_real_time){
               obp->t = world->count_scheduler->now + (obp->time_now->value - world->current_start_real_time)/world->time_unit;
               break;  
            }
          }
        }
      }
    }

      for (set=obp->data_set_head ; set!=NULL ; set=set->next)
      {
        if (set->file_flags==FILE_SUBSTITUTE)
        {
          if (world->chkpt_seq_num==1)
          {
            FILE *file = fopen(set->outfile_name,"w");
            if (file==NULL)
              mcell_perror(errno, "Failed to open reaction data output file '%s' for writing", set->outfile_name);
            fclose(file);
          }
          else if (obp->timer_type==OUTPUT_BY_ITERATION_LIST)
          {
            if(obp->time_now == NULL) continue;
            if (truncate_output_file(set->outfile_name,obp->t))
              mcell_error("Failed to prepare reaction data output file '%s' to receive output.", set->outfile_name);
          }
          else if (obp->timer_type==OUTPUT_BY_TIME_LIST)
          {
            if(obp->time_now == NULL) continue;
            if (truncate_output_file(set->outfile_name,obp->t*world->time_unit))
              mcell_error("Failed to prepare reaction data output file '%s' to receive output.", set->outfile_name);
          }
          else
          {
            /* we need to truncate up until the start of the new checkpoint 
             * simulation plus a single TIMESTEP */
            double startTime = world->chkpt_elapsed_real_time_start +  world->time_unit;
            if (truncate_output_file(set->outfile_name, startTime))
              mcell_error("Failed to prepare reaction data output file '%s' to receive output.", set->outfile_name);
          }
        }
      }
    
    if (schedule_add(world->count_scheduler , obp))
      mcell_allocfailed("Failed to add reaction data output item to scheduler.");
    obp = obpn;
  }

  no_printf("Done initializing simulation\n");
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
                                  struct object *objp,
                                  int viz_state)
{
  struct sym_table *symp;
  struct viz_child *vcp;
  switch (objp->object_type)
  {
    case META_OBJ:
      for (struct object *child_objp = objp->first_child;
           child_objp != NULL;
           child_objp=child_objp->next)
        set_viz_state_include(vizblk, child_objp, viz_state);
      break;

    case BOX_OBJ:
    case POLY_OBJ:
      symp = retrieve_sym(objp->sym->name, vizblk->viz_children);
      if (symp == NULL)
      {
        vcp = CHECKED_MALLOC_STRUCT(struct viz_child, "visualization child object");
        vcp->obj = objp;
        vcp->viz_state = NULL;
        vcp->next = NULL;
        vcp->parent = NULL;
        vcp->children = NULL;
        if (store_sym(objp->sym->name, VIZ_CHILD, vizblk->viz_children, vcp) == NULL)
          mcell_allocfailed("Failed to store visualization child object in the visualization objects table.");
      }
      else
        vcp = (struct viz_child *) symp->value;

      if (vcp->viz_state == NULL)
      {
        vcp->viz_state = CHECKED_MALLOC_ARRAY(int, objp->n_walls,
                                              "visualization state array");
        for (int i=0;i<objp->n_walls;i++)
          vcp->viz_state[i] = viz_state;
      }
      else
      {
        /* Do not override any specific viz states already set. */
        for (int i=0;i<objp->n_walls;i++)
          if (vcp->viz_state[i] == EXCLUDE_OBJ)
            vcp->viz_state[i] = viz_state;
      }
      break;

    case REL_SITE_OBJ:
      /* just do nothing */
      break;

    default:
      mcell_internal_error("Invalid object type (%d) while setting viz state.", objp->object_type);
  }
}

/*************************************************************************
 Mark all mesh objects for inclusion in the specified viz output block.

 In: vizblk: the viz output block in which to include the object
     viz_state: the desired viz state
 Out: No return value.  vizblk is updated.
*************************************************************************/
static void set_viz_all_meshes(struct viz_output_block *vizblk,
                               int viz_state)
{
  set_viz_state_include(vizblk, world->root_instance, viz_state);
}

/*************************************************************************
 Mark all molecule objects for inclusion in the specified viz output block.

 In: vizblk: the viz output block in which to include the object
     viz_state: the visualization state desired
 Out: No return value.  vizblk is updated.
*************************************************************************/
static void set_viz_all_molecules(struct viz_output_block *vizblk,
                                  int viz_state)
{
  for (int i = 0; i < world->n_species; i++)
  {
    struct species *sp = world->species_list[i];
    if (sp == world->g_mol) continue;
    if (sp->flags & IS_SURFACE) continue;
    if (vizblk->species_viz_states[i] != EXCLUDE_OBJ) continue;

    /* set viz_state to INCLUDE_OBJ for the molecule we want to visualize
       but will not assign state value */
    vizblk->species_viz_states[i] = viz_state;
  }
}

/*************************************************************************
 Count the number of selected children of a given viz_child.  (For DX mode
 only.)

 In: vcp: the viz child root to count
 Out: The total number of viz_child objects marked for inclusion
*************************************************************************/
static int count_viz_children(struct viz_child *vcp)
{
  int count = 0;

  for (; vcp != NULL; vcp = vcp->next)
  {
    if (vcp->children)
      count += count_viz_children(vcp->children);
    if (vcp->viz_state != NULL)
      ++ count;
  }

  return count;
}

/*************************************************************************
 Copy viz_child objects into an array (for DX mode only).

 In: viz: the viz object whose array should be populated
     vcp: the viz child root to copy into the array
     pos: pointer to an index into the "actual_objects" array.
 Out: No return value; viz->actual_objects is filled and *pos is updated.
*************************************************************************/
static void populate_viz_children_array(struct viz_dx_obj *viz,
                                        struct viz_child *vcp,
                                        int *pos)
{
  for (; vcp != NULL; vcp = vcp->next)
  {
    if (vcp->children)
      populate_viz_children_array(viz, vcp->children, pos);
    if (vcp->viz_state != NULL)
      viz->actual_objects[(*pos) ++] = vcp;
  }
}

/*************************************************************************
 Convert a viz_dx_obj's tree of viz_child objects into an array (for DX mode
 only).

 In: viz: the viz object whose children should be moved to an array
 Out: No return value; viz->actual_objects/viz->n_actual_objects are updated.
*************************************************************************/
static void convert_viz_children_to_array(struct viz_dx_obj *viz)
{
  int count = count_viz_children(viz->viz_child_head);
  if (count == 0)
    viz->actual_objects = NULL;
  else
    viz->actual_objects = CHECKED_MALLOC_ARRAY(struct viz_child *,
                                               count,
                                               "visualization children array");
  viz->n_actual_objects = count;

  /* Copy items into the array. */
  int pos = 0;
  populate_viz_children_array(viz, viz->viz_child_head, &pos);
  assert(pos == count);

  /* Clear the tree. */
  viz->viz_child_head = NULL;
}

/*************************************************************************
 Free all viz_child objects which represent either meta objects, or unrendered
 mesh objects.

 In: vizblk: the viz output block whose children should be trimmed
 Out: No return value; parse-time viz output data structures are freed, as are
      any excess viz_child objects.
*************************************************************************/
static void free_extra_viz_children(struct viz_output_block *vizblk)
{
  for (int i=0; i<vizblk->viz_children->n_bins; ++ i)
  {
    for (struct sym_table *sym = vizblk->viz_children->entries[i];
         sym != NULL;
         sym = sym->next)
    {
      struct viz_child *vcp = (struct viz_child *) sym->value;
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
static int viz_child_compare(void const *vc1, void const *vc2)
{
  struct viz_child const *c1 = *(struct viz_child const **) vc1;
  struct viz_child const *c2 = *(struct viz_child const **) vc2;
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
static void convert_viz_objects_to_array(struct viz_output_block *vizblk)
{
  int count = 0;
  for (int i=0; i<vizblk->viz_children->n_bins; ++ i)
  {
    for (struct sym_table *sym = vizblk->viz_children->entries[i];
         sym != NULL;
         sym = sym->next)
    {
      struct viz_child *vcp = (struct viz_child *) sym->value;
      if (vcp->viz_state)
        ++ count;
    }
  }

  /* Stash info in the visualization block. */
  vizblk->n_dreamm_objects = count;
  vizblk->dreamm_object_info = CHECKED_MALLOC_ARRAY(struct viz_child *,
                                                    count,
                                                    "DREAMM mesh objects");
  vizblk->dreamm_objects = CHECKED_MALLOC_ARRAY(struct object *,
                                                count,
                                                "DREAMM mesh objects");

  /* Now copy data in, unsorted. */
  count = 0;
  for (int i=0; i<vizblk->viz_children->n_bins; ++ i)
  {
    for (struct sym_table *sym = vizblk->viz_children->entries[i];
         sym != NULL;
         sym = sym->next)
    {
      struct viz_child *vcp = (struct viz_child *) sym->value;
      if (vcp->viz_state)
        vizblk->dreamm_object_info[count ++] = vcp;
    }
  }
  assert(count == vizblk->n_dreamm_objects);

  /* Sort the data. */
  qsort(vizblk->dreamm_object_info,
        vizblk->n_dreamm_objects,
        sizeof(struct viz_child *),
        viz_child_compare);

  /* Copy out just the objects. */
  for (int i=0; i<count; ++i)
    vizblk->dreamm_objects[i] = vizblk->dreamm_object_info[i]->obj;
}

/*************************************************************************
 Expand the mesh info for all viz children on the given viz output block.
 (DREAMM/DX modes only).

 In: vizblk: the viz output block whose children to expand
 Out: vizblk is updated
*************************************************************************/
static void expand_viz_children(struct viz_output_block *vizblk)
{
  switch (vizblk->viz_mode)
  {
    case DX_MODE:
      /* Convert viz_child tables to viz_child pointer array on each viz_dx_obj. */
      for (struct viz_dx_obj *viz = vizblk->dx_obj_head; viz != NULL; viz = viz->next)
        convert_viz_children_to_array(viz);
      free_extra_viz_children(vizblk);
      break;

    case DREAMM_V3_GROUPED_MODE:
    case DREAMM_V3_MODE:
      convert_viz_objects_to_array(vizblk);
      free_extra_viz_children(vizblk);
      break;

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
static int init_viz_species_states(struct viz_output_block *vizblk)
{
  vizblk->species_viz_states = CHECKED_MALLOC_ARRAY(int,
                                                    world->n_species,
                                                    "species viz states array");
  if (vizblk->species_viz_states == NULL)
    return 1;
  for (int i=0; i<world->n_species; ++i)
    vizblk->species_viz_states[i] = EXCLUDE_OBJ;

  int n_entries = vizblk->parser_species_viz_states.num_items;
  int n_bins = vizblk->parser_species_viz_states.table_size;
  for (int i=0; n_entries > 0  &&  i < n_bins; ++i)
  {
    struct species *specp =
          (struct species *) (vizblk->parser_species_viz_states.keys[i]);
    if (specp != NULL)
    {
      int viz_state =
            (int) (intptr_t) vizblk->parser_species_viz_states.values[i];

      /* In ASCII and RK mode, fold INCLUDE_OBJ states to 0. */
      if (vizblk->viz_mode == ASCII_MODE  &&  viz_state == INCLUDE_OBJ)
          viz_state = 0;
      if (vizblk->viz_mode == RK_MODE  &&  viz_state == INCLUDE_OBJ)
          viz_state = 0;

      vizblk->species_viz_states[specp->species_id] = viz_state;
      -- n_entries;
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
static int init_viz_output(void)
{
  for (struct viz_output_block *vizblk = world->viz_blocks;
       vizblk != NULL;
       vizblk = vizblk->next)
  {
    /* Copy species states into an array. */
    if (init_viz_species_states(vizblk))
      return 1;

    /* If ALL_MESHES or ALL_MOLECULES were requested, mark them all for inclusion. */
    if (vizblk->viz_mode != DX_MODE  &&  (vizblk->viz_output_flag & VIZ_ALL_MESHES))
      set_viz_all_meshes(vizblk, vizblk->default_mesh_state);
    if (vizblk->viz_mode != DX_MODE  &&  (vizblk->viz_output_flag & VIZ_ALL_MOLECULES))
      set_viz_all_molecules(vizblk, vizblk->default_mol_state);

    /* Copy viz children to the appropriate array. */
    expand_viz_children(vizblk);

    /* Initialize each data frame in this block. */
    if (init_frame_data_list(vizblk))
    {
      mcell_internal_error("Unknown error while initializing VIZ output.");
      return 1;
    }
  }

  return 0;
}


/********************************************************************
init_species: 
   Initializes array of molecules types to the default properties values.

*********************************************************************/
int init_species(void)
{
  int i;
  int count = 0;
  struct sym_table *gp;
  struct species *s;
  double speed;
 
  world->speed_limit = 0;
  
  count = world->mol_sym_table->n_entries;
  world->n_species = count;
  world->species_list = CHECKED_MALLOC_ARRAY(struct species*,
                                             world->n_species,
                                             "species table");
  count = 0;
  for (i=0;i<world->mol_sym_table->n_bins;i++)
  {
    for (gp = world->mol_sym_table->entries[i] ; gp != NULL ; gp = gp->next)
    {    
      if (gp->sym_type==MOL)
      {
        s = (struct species*) gp->value;
        world->species_list[count] = s;
        world->species_list[count]->species_id = count;
        world->species_list[count]->chkpt_species_id = UINT_MAX;
        world->species_list[count]->population = 0;
	world->species_list[count]->n_deceased = 0;
	world->species_list[count]->cum_lifetime = 0;
	
    
        if(!(world->species_list[count]->flags & SET_MAX_STEP_LENGTH))
        {
	   world->species_list[count]->max_step_length = DBL_MAX;
        }

        if ((s->flags & NOT_FREE) == 0)
        {
          speed = 6.0*s->space_step/sqrt(MY_PI);
          if (speed > world->speed_limit) world->speed_limit = speed;
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
static struct storage *create_storage(int nsubvols)
{
  struct storage *shared_mem = NULL;
  shared_mem = CHECKED_MALLOC_STRUCT(struct storage,
                                     "memory storage partition");
  memset(shared_mem, 0, sizeof(struct storage));

  if (world->mem_part_pool != 0)
    nsubvols = world->mem_part_pool;
  if (nsubvols < 8) nsubvols = 8;
  if (nsubvols > 4096) nsubvols = 4096;
  /* We should tune the algorithm for selecting allocation block sizes.  */
  /* XXX: Round up to power of 2?  Shouldn't matter, I think. */
  if ((shared_mem->list  = create_mem_named(sizeof(struct wall_list),nsubvols,"wall list")) == NULL)
    mcell_allocfailed("Failed to create memory pool for wall list.");
  if ((shared_mem->mol   = create_mem_named(sizeof(struct volume_molecule),nsubvols,"vol mol")) == NULL)
    mcell_allocfailed("Failed to create memory pool for volume molecules.");
  if ((shared_mem->gmol  = create_mem_named(sizeof(struct grid_molecule),nsubvols,"grid mol")) == NULL)
    mcell_allocfailed("Failed to create memory pool for grid molecules.");
  if ((shared_mem->face  = create_mem_named(sizeof(struct wall),nsubvols,"wall")) == NULL)
    mcell_allocfailed("Failed to create memory pool for walls.");
  if ((shared_mem->join  = create_mem_named(sizeof(struct edge),nsubvols,"edge")) == NULL)
    mcell_allocfailed("Failed to create memory pool for edges.");
  if ((shared_mem->tree  = create_mem_named(sizeof(struct vertex_tree),nsubvols,"vertex tree")) == NULL)
    mcell_allocfailed("Failed to create memory pool for the vertex tree.");
  if ((shared_mem->grids = create_mem_named(sizeof(struct surface_grid),nsubvols,"surface grid")) == NULL)
    mcell_allocfailed("Failed to create memory pool for surface grids.");
  if ((shared_mem->regl  = create_mem_named(sizeof(struct region_list),nsubvols,"region list")) == NULL)
    mcell_allocfailed("Failed to create memory pool for region lists.");
  if ((shared_mem->pslv  = create_mem_named(sizeof(struct per_species_list),32,"per species list")) == NULL)
    mcell_allocfailed("Failed to create memory pool for per-species molecule lists.");
  shared_mem->coll = world->coll_mem;
  shared_mem->sp_coll = world->sp_coll_mem;
  shared_mem->tri_coll = world->tri_coll_mem;
  shared_mem->exdv = world->exdv_mem;

  if (world->chkpt_init)
  {
    if ((shared_mem->timer = create_scheduler(1.0,100.0,100,0.0)) == NULL)
      mcell_allocfailed("Failed to create molecule scheduler.");
    shared_mem->current_time = 0.0;
  }

  if (world->time_step_max==0.0) shared_mem->max_timestep = MICROSEC_PER_YEAR;
  else
  {
    if (world->time_step_max < world->time_unit) shared_mem->max_timestep = 1.0;
    else shared_mem->max_timestep = world->time_step_max/world->time_unit;
  }

  return shared_mem;
}

static void sanity_check_memory_subdivision(void)
{
  if (world->mem_part_x <= 0)
  {
    if (world->mem_part_x < 0)
    {
      mcell_warn("X-axis memory partition bin size set to a negative value.  Setting to default value of 14.");
      world->mem_part_x = 14;
    }
    else
      world->mem_part_x = 10000000;
  }
  if (world->mem_part_y <= 0)
  {
    if (world->mem_part_y < 0)
    {
      mcell_warn("Y-axis memory partition bin size set to a negative value.  Setting to default value of 14.");
      world->mem_part_y = 14;
    }
    else
      world->mem_part_y = 10000000;
  }
  if (world->mem_part_z <= 0)
  {
    if (world->mem_part_z < 0)
    {
      mcell_warn("Z-axis memory partition bin size set to a negative value.  Setting to default value of 14.");
      world->mem_part_z = 14;
    }
    else
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
int init_partitions(void)
{
  int i,j,k,h;
  struct subvolume *sv;

  /* Initialize the partitions, themselves */
  if (set_partitions()) return 1;
  
  /* Initialize dummy waypoints (why do we do this?) */
  world->n_waypoints = 1;
  world->waypoints = CHECKED_MALLOC_ARRAY(struct waypoint,
                                          world->n_waypoints,
                                          "dummy waypoint");

  /* Allocate the subvolumes */
  world->n_subvols = (world->nz_parts-1) * (world->ny_parts-1) * (world->nx_parts-1);
  if (world->notify->progress_report!=NOTIFY_NONE)
    mcell_log("Creating %d subvolumes (%d,%d,%d per axis).", world->n_subvols, world->nx_parts-1, world->ny_parts-1, world->nz_parts-1);
  world->subvol = CHECKED_MALLOC_ARRAY(struct subvolume, world->n_subvols, "spatial subvolumes");

  /* Decide how fine-grained to make the memory subdivisions */
  sanity_check_memory_subdivision();

  /* Allocate the data structures which are shared between storages */
  if ((world->coll_mem  = create_mem_named(sizeof(struct collision),128,"collision")) == NULL)
    mcell_allocfailed("Failed to create memory pool for collisions.");
  if ((world->sp_coll_mem  = create_mem_named(sizeof(struct sp_collision),128,"sp collision")) == NULL)
    mcell_allocfailed("Failed to create memory pool for trimolecular-pathway collisions.");
  if ((world->tri_coll_mem  = create_mem_named(sizeof(struct tri_collision),128,"tri collision")) == NULL)
    mcell_allocfailed("Failed to create memory pool for trimolecular collisions.");
  if ((world->exdv_mem  = create_mem_named(sizeof(struct exd_vertex),64,"exact disk vertex")) == NULL)
    mcell_allocfailed("Failed to create memory pool for exact disk calculation vertices.");

  /* How many storage subdivisions along each axis? */
  int nx = (world->nx_parts + (world->mem_part_x) - 2) / (world->mem_part_x);
  int ny = (world->ny_parts + (world->mem_part_y) - 2) / (world->mem_part_y);
  int nz = (world->nz_parts + (world->mem_part_z) - 2) / (world->mem_part_z);
  if (world->notify->progress_report!=NOTIFY_NONE)
    mcell_log("Creating %d memory partitions (%d,%d,%d per axis).", nx*ny*nz, nx, ny, nz);

  /* Create memory pool for storages */
  if ((world->storage_allocator = create_mem_named(sizeof(struct storage_list),nx*ny*nz,"storage allocator")) == NULL)
    mcell_allocfailed("Failed to create memory pool for storage list.");

  /* Allocate the storages */
  struct storage *shared_mem[nx*ny*nz];
  int cx = 0, cy = 0, cz = 0;
  for (i=0; i<nx*ny*nz; ++i)
  {
    /* Determine the number of subvolumes included in this subdivision */
    int xd = world->mem_part_x, yd = world->mem_part_y, zd = world->mem_part_z;
    if (cx == nx-1)
      xd = (world->nx_parts - 1) % world->mem_part_x;
    if (cy == ny-1)
      yd = (world->ny_parts - 1) % world->mem_part_y;
    if (cz == nz-1)
      zd = (world->nz_parts - 1) % world->mem_part_z;
    if (++ cx == nx)
    {
      cx = 0;
      if (++ cy == ny)
      {
        cy = 0;
        ++ cz;
      }
    }

    /* Allocate this storage */
    if ((shared_mem[i] = create_storage(xd*yd*zd)) == NULL)
      mcell_internal_error("Unknown error while creating a storage.");

    /* Add to the storage list */
    struct storage_list *l = (struct storage_list*) CHECKED_MEM_GET(world->storage_allocator, "storage list item");
    l->next = world->storage_head;
    l->store = shared_mem[i];
    world->storage_head = l;
  }

  /* Initialize each subvolume */
  for (i=0;i<world->nx_parts-1;i++)
  for (j=0;j<world->ny_parts-1;j++)
  for (k=0;k<world->nz_parts-1;k++)
  {
    h = k + (world->nz_parts-1)*(j + (world->ny_parts-1)*i);
    sv = & (world->subvol[ h ]);
    sv->wall_head = NULL;
    memset(&sv->mol_by_species, 0, sizeof(struct pointer_hash));
    sv->species_head = NULL;
    sv->mol_count = 0;
    
    sv->llf.x = bisect_near( world->x_fineparts , world->n_fineparts , world->x_partitions[i] );
    sv->llf.y = bisect_near( world->y_fineparts , world->n_fineparts , world->y_partitions[j] );
    sv->llf.z = bisect_near( world->z_fineparts , world->n_fineparts , world->z_partitions[k] );
    sv->urb.x = bisect_near( world->x_fineparts , world->n_fineparts , world->x_partitions[i+1] );
    sv->urb.y = bisect_near( world->y_fineparts , world->n_fineparts , world->y_partitions[j+1] );
    sv->urb.z = bisect_near( world->z_fineparts , world->n_fineparts , world->z_partitions[k+1] );
    
    /* Set flags so we know which directions to not go (we will fall off the world!) */
    sv->world_edge=0;  /* Assume we're not at the edge of the world in any direction */
    if (i==0) sv->world_edge |= X_NEG_BIT;
    if (i==world->nx_parts-2) sv->world_edge |= X_POS_BIT;
    if (j==0) sv->world_edge |= Y_NEG_BIT;
    if (j==world->ny_parts-2) sv->world_edge |= Y_POS_BIT;
    if (k==0) sv->world_edge |= Z_NEG_BIT;
    if (k==world->nz_parts-2) sv->world_edge |= Z_POS_BIT;
   
    /* This part is commented because we are not using
       bsp_trees at this time */ 
   /*
     sv->is_bsp = 0; 
         
    if (i==0) sv->neighbor[X_NEG] = NULL;
    else sv->neighbor[X_NEG] = &(world->subvol[ h - (world->nz_parts-1)*(world->ny_parts-1) ]);
    
    if (i==world->nx_parts-2) sv->neighbor[X_POS] = NULL;
    else sv->neighbor[X_POS] = &(world->subvol[ h + (world->nz_parts-1)*(world->ny_parts-1) ]);
    
    if (j==0) sv->neighbor[Y_NEG] = NULL;
    else sv->neighbor[Y_NEG] = &(world->subvol[ h - (world->nz_parts-1) ]);
    
    if (j==world->ny_parts-2) sv->neighbor[Y_POS] = NULL;
    else sv->neighbor[Y_POS] = &(world->subvol[ h + (world->nz_parts-1) ]);
    
    if (k==0) sv->neighbor[Z_NEG] = NULL;
    else sv->neighbor[Z_NEG] = &(world->subvol[ h - 1 ]);
    
    if (k==world->nz_parts-2) sv->neighbor[Z_POS] = NULL;
    else sv->neighbor[Z_POS] = &(world->subvol[ h + 1 ]);
       */

    /* Bind this subvolume to the appropriate storage */
    int shidx = (i / (world->mem_part_x)) + nx * (j / (world->mem_part_y) + ny * (k / (world->mem_part_z)));
    sv->local_storage = shared_mem[shidx];
  }
  
  return 0;
}

/**
 * Initializes the geometry of the world.
 * Calls instance_obj() to instantiate all physical objects.
 * (Meta objects, box objects, polygon objects and release sites)
 */
int init_geom(void)
{
  double tm[4][4];
  double vol_infinity;
  
  no_printf("Initializing physical objects\n");
  vol_infinity=sqrt(DBL_MAX)/4;
  world->bb_llf.x=vol_infinity;
  world->bb_llf.y=vol_infinity;
  world->bb_llf.z=vol_infinity;
  world->bb_urb.x=-vol_infinity;
  world->bb_urb.y=-vol_infinity;
  world->bb_urb.z=-vol_infinity;
  init_matrix(tm);
  
  if (compute_bb(world->root_instance,tm))
    return 1;

  if (world->bb_llf.x==vol_infinity 
      && world->bb_llf.y==vol_infinity
      && world->bb_llf.z==vol_infinity
      && world->bb_urb.x==-vol_infinity
      && world->bb_urb.y==-vol_infinity
      && world->bb_urb.z==-vol_infinity) {
    world->bb_llf.x=0;
    world->bb_llf.y=0;
    world->bb_llf.z=0;
    world->bb_urb.x=0;
    world->bb_urb.y=0;
    world->bb_urb.z=0;
  }
  if (world->procnum == 0) {
    if (world->notify->progress_report)
    {
      mcell_log("MCell: world bounding box in microns =");
      mcell_log("         [ %.9g %.9g %.9g ] [ %.9g %.9g %.9g ]",
        world->bb_llf.x*world->length_unit,world->bb_llf.y*world->length_unit,
        world->bb_llf.z*world->length_unit,world->bb_urb.x*world->length_unit,
        world->bb_urb.y*world->length_unit,world->bb_urb.z*world->length_unit);
    }
  }

  world->n_walls=world->root_instance->n_walls;
  world->n_verts=world->root_instance->n_verts;
  no_printf("World object contains %d walls and %d vertices\n",
    world->n_walls,world->n_verts);
  
  if (instance_obj(world->root_instance, tm))
    return 1;

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
int instance_obj(struct object *objp, double (*im)[4])
{
  double tm[4][4];
  mult_matrix(objp->t_matrix, im, tm, 4, 4, 4);

  switch (objp->object_type)
  {
    case META_OBJ:
      for (struct object *child_objp = objp->first_child;
           child_objp != NULL;
           child_objp = child_objp->next)
      {
        if (instance_obj(child_objp, tm))
          return 1;
      }
      break;

    case REL_SITE_OBJ:
      if (instance_release_site(objp, tm))
        return 1;
      break;

    case BOX_OBJ:
    case POLY_OBJ:
      if (instance_polygon_object(objp, tm))
        return 1;
      break;

    default:
      UNHANDLED_CASE(objp->object_type);
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
int instance_release_site(struct object *objp, double (*im)[4])
{
  struct release_site_obj *rsop;
  struct release_event_queue *reqp;
  int i,j;

  rsop=(struct release_site_obj *)objp->contents;
  
  no_printf("Instancing release site object %s\n",objp->sym->name);
  if (rsop->release_prob==MAGIC_PATTERN_PROBABILITY)
  {
    struct magic_list *ml;
    struct rxn_pathname *rxpn;
    
    ml = (struct magic_list*) CHECKED_MEM_GET(world->magic_mem, "rxn-triggered release descriptor");
    ml->data = rsop;
    ml->type = magic_release;
    
    rxpn = (struct rxn_pathname*)rsop->pattern;
    ml->next=rxpn->magic;
    rxpn->magic=ml;
    
    /* Region releases need to be in release queue to get initialized */
    /* Release code itself is smart enough to ignore MAGIC_PATTERNs */
    if (rsop->release_shape==SHAPE_REGION) 
    {
      reqp = CHECKED_MALLOC_STRUCT(struct release_event_queue, "release site");
      reqp->release_site=rsop;
      reqp->event_time=0;
      reqp->train_counter=0;
      reqp->train_high_time=0;
      if (schedule_add(world->releaser,reqp))
        mcell_allocfailed("Failed to schedule molecule release.");
    }
  }
  else
  {
    reqp = CHECKED_MALLOC_STRUCT(struct release_event_queue, "release site");
    reqp->release_site=rsop;
    reqp->event_time=rsop->pattern->delay;
    reqp->train_counter=0;
    reqp->train_high_time=rsop->pattern->delay;
    for (i=0;i<4;i++) for (j=0;j<4;j++) reqp->t_matrix[i][j]=im[i][j];
  
    /* Schedule the release event */
    if (schedule_add(world->releaser,reqp))
      mcell_allocfailed("Failed to schedule molecule release.");
  
    if (rsop->pattern->train_duration > rsop->pattern->train_interval)
      mcell_error("Release pattern train duration is greatter than train interval.");
  }
  
  no_printf("Done instancing release site object %s\n",objp->sym->name);

  return 0;
}



/**
 * Computes the bounding box for the entire simulation world.
 * Does things recursively in a manner similar to instance_obj().
 */
static int compute_bb(struct object *objp, double (*im)[4])
{
  double tm[4][4];
  mult_matrix(objp->t_matrix,im,tm,4,4,4);

  switch (objp->object_type)
  {
    case META_OBJ:
      for (struct object *child_objp = objp->first_child;
           child_objp != NULL;
           child_objp = child_objp->next)
      {
        if (compute_bb(child_objp, tm))
          return 1;
      }
      break;

    case REL_SITE_OBJ:
      if (compute_bb_release_site(objp, tm))
        return 1;
      break;

    case BOX_OBJ:
    case POLY_OBJ:
      if (compute_bb_polygon_object(objp, tm))
        return 1;
      break;

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
static int compute_bb_release_site(struct object *objp, double (*im)[4])
{
  struct release_site_obj *rsop;
  double location[1][4];
  double diam_x, diam_y, diam_z; /* diameters of the release_site */ 
 
  rsop=(struct release_site_obj *)objp->contents;
  
  if (rsop->release_shape == SHAPE_REGION) return 0;

  if (rsop->location == NULL)
    mcell_error("Location is not specified for the geometrical shape release site '%s'.", objp->sym->name);

  location[0][0]=rsop->location->x;
  location[0][1]=rsop->location->y;
  location[0][2]=rsop->location->z;
  location[0][3]=1.0;
  mult_matrix(location,im,location,1,4,4);
  
  if(rsop->diameter == NULL){
	diam_x = diam_y = diam_z = 0;
  }else{
        diam_x = rsop->diameter->x;
        diam_y = rsop->diameter->y;
        diam_z = rsop->diameter->z;
  }


  if (location[0][0]  - diam_x < world->bb_llf.x) {
    world->bb_llf.x=location[0][0] - diam_x;
  }
  if (location[0][1] - diam_y < world->bb_llf.y) {
    world->bb_llf.y=location[0][1] - diam_y;
  }
  if (location[0][2] - diam_z < world->bb_llf.z) {
    world->bb_llf.z=location[0][2] - diam_z;
  }
  if (location[0][0] + diam_x > world->bb_urb.x) {
    world->bb_urb.x=location[0][0] + diam_x;
  }
  if (location[0][1] + diam_y > world->bb_urb.y) {
    world->bb_urb.y=location[0][1] + diam_y;
  }
  if (location[0][2] + diam_z > world->bb_urb.z) {
    world->bb_urb.z=location[0][2] + diam_z;
  }

  return 0;
}

/**
 * Updates the bounding box of the world based on the size
 * and location of a polygon_object.
 * Used by compute_bb().
 */
static int compute_bb_polygon_object(struct object *objp, double (*im)[4])
{
  struct polygon_object *pop;

  pop=(struct polygon_object *)objp->contents;
  const unsigned int n_verts = pop->n_verts;

  for (unsigned int n_vert=0; n_vert<n_verts; ++ n_vert)
  {
    double p[1][4];
    p[0][0] = pop->vertex[n_vert].x;
    p[0][1] = pop->vertex[n_vert].y;
    p[0][2] = pop->vertex[n_vert].z;
    p[0][3] = 1.0;
    mult_matrix(p,im,p,1,4,4);
    if (p[0][0]<world->bb_llf.x) world->bb_llf.x = p[0][0];
    if (p[0][1]<world->bb_llf.y) world->bb_llf.y = p[0][1];
    if (p[0][2]<world->bb_llf.z) world->bb_llf.z = p[0][2];
    if (p[0][0]>world->bb_urb.x) world->bb_urb.x = p[0][0];
    if (p[0][1]>world->bb_urb.y) world->bb_urb.y = p[0][1];
    if (p[0][2]>world->bb_urb.z) world->bb_urb.z = p[0][2];
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
int instance_polygon_object(struct object *objp, double (*im)[4])
{
// #define INIT_VERTEX_NORMALS
// Uncomment to compute vertex normals
  struct polygon_object *pop;
  struct vector3 *v,**vp;
  struct wall *w,**wp;
  double p[1][4];
#ifdef INIT_VERTEX_NORMALS
  struct vector3 *vertex_normal;
  double origin[1][4];
  byte compute_vertex_normals = 0;
#endif
  double total_area;
  int index_0,index_1,index_2;
  unsigned int degenerate_count;

  pop=(struct polygon_object *)objp->contents;
  const unsigned int n_walls = pop->n_walls;
  const unsigned int n_verts = pop->n_verts;
  total_area=0;

/* Allocate and initialize walls and vertices */
    w  = CHECKED_MALLOC_ARRAY(struct wall,      n_walls, "polygon walls");
    wp = CHECKED_MALLOC_ARRAY(struct wall *,    n_walls, "polygon wall pointers");
    v  = CHECKED_MALLOC_ARRAY(struct vector3,   n_verts, "polygon vertices");
    vp = CHECKED_MALLOC_ARRAY(struct vector3 *, n_verts, "polygon vertex pointers");
    objp->shared_walls = CHECKED_MALLOC_ARRAY(struct wall_list *, n_verts, "wall list pointers");  
    objp->walls=w;
    objp->wall_p=wp;
    objp->verts=v;
           
    for(u_int i = 0; i < n_verts; i++)
    {
       objp->shared_walls[i] = NULL;
    }
              

/* If we want vertex normals we'll have to add a place to store them
   in struct object.
*/
#ifdef INIT_VERTEX_NORMALS
    if (pop->normal!=NULL)
      compute_vertex_normals = 1;
#endif

  for (unsigned int n_vert=0; n_vert<n_verts; ++ n_vert)
  {
    vp[n_vert] = &v[n_vert];
    p[0][0] = pop->vertex[n_vert].x;
    p[0][1] = pop->vertex[n_vert].y;
    p[0][2] = pop->vertex[n_vert].z;
    p[0][3] = 1.0;
    mult_matrix(p,im,p,1,4,4);
    v[n_vert].x = p[0][0];
    v[n_vert].y = p[0][1];
    v[n_vert].z = p[0][2];

#ifdef INIT_VERTEX_NORMALS
    if (compute_vertex_normals)
    {
      p[0][0] = pop->normal[n_vert].x;
      p[0][1] = pop->normal[n_vert].y;
      p[0][2] = pop->normal[n_vert].z;
      p[0][3] = 1.0;
      origin[0][0]=0;
      origin[0][1]=0;
      origin[0][2]=0;
      origin[0][3]=1.0;
      mult_matrix(p,im,p,1,4,4);
      mult_matrix(origin,im,origin,1,4,4);
      vertex_normal[n_vert].x = p[0][0]-origin[0][0];
      vertex_normal[n_vert].y = p[0][1]-origin[0][1];
      vertex_normal[n_vert].z = p[0][2]-origin[0][2];
      normalize(&vertex_normal[n_vert]);
    }
#endif
  }
  
  degenerate_count=0;
  for (unsigned int n_wall=0; n_wall<n_walls; ++ n_wall)
  {
    if (!get_bit(pop->side_removed, n_wall)) {
      wp[n_wall] = &w[n_wall];
      index_0 = pop->element[n_wall].vertex_index[0];
      index_1 = pop->element[n_wall].vertex_index[1];
      index_2 = pop->element[n_wall].vertex_index[2];

      init_tri_wall(objp,n_wall,vp[index_0],vp[index_1],vp[index_2], index_0, index_1, index_2);
      total_area+=wp[n_wall]->area;

      if (wp[n_wall]->area==0)
      {
        if (world->notify->degenerate_polys != WARN_COPE)
        {
          if (world->notify->degenerate_polys==WARN_ERROR)
          {
            mcell_error("Degenerate polygon found: %s %d\n"
                        "  Vertex 0: %.5e %.5e %.5e\n"
                        "  Vertex 1: %.5e %.5e %.5e\n"
                        "  Vertex 2: %.5e %.5e %.5e",
                        objp->sym->name, n_wall,
                        vp[index_0]->x, vp[index_0]->y, vp[index_0]->z,
                        vp[index_1]->x, vp[index_1]->y, vp[index_1]->z,
                        vp[index_2]->x, vp[index_2]->y, vp[index_2]->z);
          }
          else
            mcell_warn("Degenerate polygon found and automatically removed: %s %d\n"
                       "  Vertex 0: %.5e %.5e %.5e\n"
                       "  Vertex 1: %.5e %.5e %.5e\n"
                       "  Vertex 2: %.5e %.5e %.5e",
                       objp->sym->name, n_wall,
                       vp[index_0]->x, vp[index_0]->y, vp[index_0]->z,
                       vp[index_1]->x, vp[index_1]->y, vp[index_1]->z,
                       vp[index_2]->x, vp[index_2]->y, vp[index_2]->z);
        }
        set_bit(pop->side_removed,n_wall,1);
        objp->n_walls_actual--;
        degenerate_count++;
        wp[n_wall]=NULL;
      }
    }
    else {
      wp[n_wall]=NULL;
    }
  }
  if (degenerate_count) mdl_remove_gaps_from_regions(objp);
  
  objp->total_area=total_area;
  
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
int init_regions(void)
{
  if (world->clamp_list!=NULL) init_clamp_lists();

  return instance_obj_regions(world->root_instance);
}


/* First part of concentration clamp initialization. */
/* After this, list is grouped by surface class. */
/* Second part (list of objects) happens with regions. */
void init_clamp_lists(void)
{
  struct ccn_clamp_data *ccd,*temp;
  
  /* Sort by memory order of surface_class pointer--handy way to collect like classes */
  world->clamp_list = (struct ccn_clamp_data*)void_list_sort((struct void_list*)world->clamp_list);
  
  /* Toss other molecules in same surface class into next_mol lists */
  for (ccd = world->clamp_list ; ccd!=NULL ; ccd=ccd->next)
  {
    while (ccd->next != NULL && ccd->surf_class==ccd->next->surf_class)
    {
      ccd->next->next_mol = ccd->next_mol;
      ccd->next_mol = ccd->next;
      ccd->next = ccd->next->next;
    }
    for (temp=ccd->next_mol ; temp!=NULL ; temp=temp->next_mol)
    {
      temp->next = ccd->next;
    }
  }
}


/**
 * Traverse through metaobjects, placing regions on real objects as we find
 * them.
 */
int instance_obj_regions(struct object *objp)
{
  switch (objp->object_type)
  {
    case META_OBJ:
      for (struct object *child_objp = objp->first_child;
           child_objp != NULL;
           child_objp = child_objp->next)
      {
        if (instance_obj_regions(child_objp))
          return 1;
      }
      break;

    case REL_SITE_OBJ:
      break;

    case BOX_OBJ:
    case POLY_OBJ:
      if (init_wall_regions(objp))
        return 1;
      break;

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
 * Populates effector tiles by region.
 * Creates virtual regions on which to clamp concentration
 */
int init_wall_regions(struct object *objp)
{
  struct wall *w;
  struct region *rp;
  struct region_list *rlp,*wrlp;

  const struct polygon_object *pop = (struct polygon_object *) objp->contents;
  const unsigned int n_walls = pop->n_walls;

  no_printf("Processing %d regions in polygon list object: %s\n",objp->num_regions, objp->sym->name);  

  /* prepend a copy of eff_dat for each element referenced in each region
     of this object to the eff_prop list for the referenced element */
  for (rlp=objp->regions; rlp!=NULL; rlp=rlp->next)
  {
    rp = rlp->reg;
    if (rp->membership==NULL)
      mcell_internal_error("Missing region information for '%s'.", rp->sym->name);

    for (int n_wall=0; n_wall<rp->membership->nbits; ++ n_wall)
    {
      if (get_bit(rp->membership, n_wall))
      {
	/* prepend this region to wall region list of i_th wall only if the region is used in counting */
        w = objp->wall_p[n_wall];
	rp->area += w->area;
	if (rp->surf_class!=NULL) w->surf_class = rp->surf_class;

	if ((rp->flags & COUNT_SOME_MASK) != 0)
	{  
          wrlp = (struct region_list *) CHECKED_MEM_GET(w->birthplace->regl, "wall region list");
	  wrlp->reg=rp;
	  wrlp->next=w->counting_regions;
	  w->counting_regions=wrlp;
	  w->flags|=rp->flags;
	}
      }
    }

  } /*end loop over all regions in object */
  for (unsigned int n_wall=0; n_wall<n_walls; n_wall++)
  {
    if (get_bit(pop->side_removed, n_wall)) continue; 
     
    w = objp->wall_p[n_wall];
    if (w->counting_regions!=NULL)
    {
      w->counting_regions = (struct region_list*)void_list_sort((struct void_list*)w->counting_regions); /* Helpful for comparisons */
    }
  }
      
  no_printf("Total area of object %s = %.9g um^2\n",objp->sym->name,objp->total_area/world->grid_density);
  no_printf("  number of tiles = %u\n",objp->n_tiles);
  no_printf("  number of occupied tiles = %u\n",objp->n_occupied_tiles);
  no_printf("  grid molecule density = %.9g\n",objp->n_occupied_tiles*world->grid_density/objp->total_area);
    
  /* Check to see if we need to generate virtual regions for */
  /* concentration clamps on this object */
  if (world->clamp_list!=NULL)
  {
    struct ccn_clamp_data *ccd;
    struct ccn_clamp_data *temp;
    int j;
    int found_something = 0;
    
    for (unsigned int n_wall=0;n_wall<n_walls;n_wall++)
    {
      if (get_bit(pop->side_removed, n_wall)) continue;
      if (objp->wall_p[n_wall]->surf_class != world->g_surf)
      {
        for (ccd=world->clamp_list ; ccd!=NULL ; ccd=ccd->next)
        {
          if (objp->wall_p[n_wall]->surf_class == ccd->surf_class)
          {
            if (ccd->objp!=objp)
            {
              if (ccd->objp==NULL) ccd->objp=objp;
              else if (ccd->next_obj != NULL && ccd->next_obj->objp==objp) ccd=ccd->next_obj;
              else
              {
                temp = CHECKED_MALLOC_STRUCT(struct ccn_clamp_data, "concentration clamp data");
                memcpy(temp,ccd,sizeof(struct ccn_clamp_data));
                temp->objp = objp;
                temp->sides = NULL;
                temp->n_sides = 0;
                temp->side_idx = NULL;
                temp->cum_area = NULL;
                ccd->next_obj = temp;
                ccd = temp;
              }
            }
            if (ccd->sides==NULL)
            {
              ccd->sides = new_bit_array(n_walls);
              if (ccd->sides==NULL)
                mcell_allocfailed("Failed to allocate membership bit array for concentration clamp data.");
              set_all_bits(ccd->sides,0);
            }
            set_bit(ccd->sides, n_wall, 1);
            ccd->n_sides++;
            found_something=1;
          }
        }
      }
    }
    
    if (found_something)
    {
      for (ccd=world->clamp_list ; ccd!=NULL ; ccd=ccd->next)
      {
        if (ccd->objp!=objp)
        {
          if (ccd->next_obj!=NULL && ccd->next_obj->objp==objp) ccd=ccd->next_obj;
          else continue;
        }
        
        ccd->side_idx = CHECKED_MALLOC_ARRAY(int,
                                             ccd->n_sides,
                                             "concentration clamp polygon side index");
        ccd->cum_area = CHECKED_MALLOC_ARRAY(double,
                                             ccd->n_sides,
                                             "concentration clamp polygon side cumulative area");
        
        j=0;
        for (unsigned int n_wall=0; n_wall<n_walls; n_wall++)
        {
          if (get_bit(ccd->sides, n_wall))
          {
            ccd->side_idx[j] = n_wall;
            ccd->cum_area[j] = objp->wall_p[n_wall]->area;
            j++;
          }
        }
        if (j != ccd->n_sides)
          mcell_internal_error("Miscounted the number of walls for concentration clamp.  object=%s  surface class=%s",
                               objp->sym->name,
                               ccd->surf_class->sym->name);
        
        for (j=1;j<ccd->n_sides;j++) ccd->cum_area[j] += ccd->cum_area[j-1];
        
        ccd->scaling_factor = ccd->cum_area[ccd->n_sides-1] * world->length_unit * world->length_unit * world->length_unit /
                              2.9432976599069717358e-9;  /* sqrt(MY_PI)/(1e-15*N_AV) */
        if (ccd->orient!=0) ccd->scaling_factor *= 0.5;
      }
    }
    
  }


#ifdef KELP
  cdp->sym->ref_count--;
  if (!cdp->sym->ref_count) {	/* Done with the geometry information */
	destroy_sym_value(cdp->sym);	/* free up memory */
  }
#endif

  return 0;
}



/********************************************************************
 init_effectors:

    Traverse the world placing grid molecules.

    In:  none
    Out: 0 on success, 1 on failure
 *******************************************************************/
int init_effectors(void)
{
  if (instance_obj_effectors(world->root_instance)) return 1;
  return 0;
}


/********************************************************************
 instance_obj_effectors:

    Place any appropriate grid molecules on this object and/or its children.

    In:  struct object *objp - the object upon which to instantiate molecules
    Out: 0 on success, 1 on failure
 *******************************************************************/
int instance_obj_effectors(struct object *objp)
{
  struct object *child_objp;

  switch (objp->object_type)
  {
    case META_OBJ:
      for (child_objp=objp->first_child ; child_objp!=NULL ; child_objp=child_objp->next)
      {
	if (instance_obj_effectors(child_objp)) return 1;
      }
      break;
    case REL_SITE_OBJ:
      break;
    case BOX_OBJ:
    case POLY_OBJ:
      if (init_wall_effectors(objp)) return 1;
      break;
    default:
      break;
  }

  return 0;
}


/********************************************************************
 init_wall_effectors:

    Place any appropriate grid molecules on this wall.  The object passed in
    must be a box or a polygon.

    In:  struct object *objp - the object upon which to instantiate molecules
    Out: 0 on success, 1 on failure
 *******************************************************************/
int init_wall_effectors(struct object *objp)
{
  struct wall *w;
  struct eff_dat *effdp,*dup_effdp,**eff_prop;
  struct region *rp;
  struct region_list *rlp,*rlp2,*reg_eff_num_head,*complex_head;
  byte reg_eff_num;
  byte complex_eff;

  const struct polygon_object *pop = (struct polygon_object *) objp->contents;
  const unsigned int n_walls = pop->n_walls;

  /* allocate scratch storage to hold effector info for each wall */
  eff_prop = CHECKED_MALLOC_ARRAY(struct eff_dat *, n_walls, "effector data scratch space");

  for (unsigned int n_wall=0; n_wall<n_walls; ++ n_wall) eff_prop[n_wall]=NULL; 

  /* prepend a copy of eff_dat for each element referenced in each region
     of this object to the eff_prop list for the referenced element */
  reg_eff_num_head=NULL;

  /* List of regions which need macromol processing */
  complex_head=NULL;
  
  for (rlp=objp->regions ; rlp!=NULL ; rlp=rlp->next)
  {
    rp=rlp->reg;
    reg_eff_num=0;
    complex_eff=0;

    for (int n_wall=0; n_wall<rp->membership->nbits; n_wall++)
    {
      if (get_bit(rp->membership, n_wall))
      {
        w = objp->wall_p[n_wall];

	/* prepend region eff data for this region to eff_prop for i_th wall */
        for ( effdp=rp->eff_dat_head ; effdp!=NULL ; effdp=effdp->next )
	{
          if (effdp->eff->flags & IS_COMPLEX)
            complex_eff = 1;
          else if (effdp->quantity_type==EFFDENS)
	  {
            dup_effdp = CHECKED_MALLOC_STRUCT(struct eff_dat, "effector data");
	    dup_effdp->eff=effdp->eff;
	    dup_effdp->quantity_type=effdp->quantity_type;
	    dup_effdp->quantity=effdp->quantity;
	    dup_effdp->orientation=effdp->orientation;
            dup_effdp->next = eff_prop[n_wall];
            eff_prop[n_wall] = dup_effdp;
	  }
	  else reg_eff_num=1;
	}

	/* prepend surf_class eff data for this region to eff_prop for i_th wall on last region */
	if (w->surf_class != NULL && rlp->next==NULL)
	{
	  for ( effdp=w->surf_class->eff_dat_head ; effdp!=NULL ; effdp=effdp->next )
	  {
            if (effdp->eff->flags & IS_COMPLEX)
              complex_eff = 1;
            else if (effdp->quantity_type==EFFDENS)
	    {
              dup_effdp = CHECKED_MALLOC_STRUCT(struct eff_dat, "effector data");
	      dup_effdp->eff=effdp->eff;
	      dup_effdp->quantity_type=effdp->quantity_type;
	      dup_effdp->quantity=effdp->quantity;
	      dup_effdp->orientation=effdp->orientation;
              dup_effdp->next = eff_prop[n_wall];
              eff_prop[n_wall] = dup_effdp;
	    }
	    else reg_eff_num=1;
	  }
	}

      }
    } /* done checking each wall */
    
    if (complex_eff)
    {
      rlp2 = CHECKED_MALLOC_STRUCT(struct region_list, "complex effector placement region list");
      rlp2->reg=rp;
      rlp2->next=complex_head;
      complex_head=rlp2;
    }
    else if (reg_eff_num)
    {
      rlp2 = CHECKED_MALLOC_STRUCT(struct region_list, "effector placement region list");
      rlp2->reg=rp;
      rlp2->next=reg_eff_num_head;
      reg_eff_num_head=rlp2;
    }
  } /*end for (... ; rlp != NULL ; ...) */

  /* Place macromolecular complexes, if any */
  if (complex_head!=NULL)
  {
    if (init_complex_effectors(objp, complex_head))
      return 1;
    
    /* free list of regions with complex effectors */
    rlp=complex_head;
    while(rlp!=NULL)
    {
      rlp2=rlp;
      rlp=rlp->next;
      free(rlp2);
    }
  }

  /* Place regular (non-macro) molecules by density */
  for (unsigned int n_wall=0; n_wall<n_walls; n_wall++)
  {
    if (!get_bit(pop->side_removed, n_wall))
    {
      if (eff_prop[n_wall]!=NULL)
      {
        if (init_effectors_by_density(objp->wall_p[n_wall],eff_prop[n_wall]))
          return 1;
      }
    }
  }

  /* Place regular (non-macro) molecules by number */
  if (reg_eff_num_head!=NULL)
  {
    if (init_effectors_by_number(objp,reg_eff_num_head))
      return 1;
    
    /* free region list created to hold regions populated by number */
    rlp=reg_eff_num_head;
    while(rlp!=NULL)
    {
      rlp2=rlp;
      rlp=rlp->next;
      free(rlp2);
    }
  }

  /* free eff_prop array and contents */
  for (unsigned int n_wall=0; n_wall<n_walls; n_wall++)
  {
    if (eff_prop[n_wall] != NULL)
    {
      effdp = eff_prop[n_wall];
      while(effdp!=NULL)
      {
	dup_effdp=effdp;
	effdp=effdp->next;
	free(dup_effdp);
      }
    }
  }
  free(eff_prop);

    
  return 0;
}

/********************************************************************
 init_effectors_place_complex:

    Place a surface macromolecule on this wall.  Note that if a region is
    specified, the release will constrain all subunits to be within the region.
    Technically, it constrains the release to only occur if the path from the
    macromolecule's origin to the location of each subunit does not leave the
    region.  In practice, the distinction is unlikely to be important.

    In:  struct wall *w - the wall upon which to instantiate molecule
         struct region *rp - the region in which to instantiate molecule, or
                    NULL if the molecule is not restricted to a single region
         struct eff_dat const *effdp - description of what to release
    Out: 0 on success, 1 on failure

    Program state is not corrupted if this fails; either the entire
    macromolecule is released, or none of it is.
 *******************************************************************/
static int init_effectors_place_complex(struct wall *w,
                                        struct region *rp,
                                        struct eff_dat const *effdp)
{
  struct grid_molecule *gp;
  unsigned int grid_idx;
  double p;
  short orient = effdp->orientation;

  /* Pick orientation */
  if (orient == 0)
  {
    orient = (rng_uint(world->rng) & 1) ? 1 : -1;
  }

  /* Pick location on wall */
  p = rng_dbl(world->rng);
  grid_idx = p * (double) (w->grid->n * w->grid->n);
  if (grid_idx >= w->grid->n_tiles)
    grid_idx = w->grid->n_tiles - 1;

  gp = macro_insert_molecule_grid_2(effdp->eff, orient, w, grid_idx, 0.0, rp, NULL);
  return (gp != NULL) ? 0 : 1;
}

/********************************************************************
 init_effectors_place_complexes:

    Place any appropriate surface macromolecules on these walls.  Note that if
    a region is specified, the release will constrain all subunits to be within
    the region.  Technically, it constrains each molecule release to only occur
    if the path from the macromolecule's origin to the location of each subunit
    within its complex lies entirely within the region.  In practice, the
    distinction is unlikely to be important most of the time.

    In:  int n_to_place - how many to release
         int nwalls - upon how many walls do we release?
         double *weights - array of cumulative areas for the first n walls in
                  the walls array; used for weighted samping of the area
         struct wall * const *walls - the walls upon which to release
         struct region *rp - the region in which to instantiate molecules, or
                    NULL if the molecule is not restricted to a single region
         struct eff_dat const *effdp - description of what to release
    Out: 0 on success, 1 on failure

    Program state is not corrupted if this fails.  This may release fewer
    macromolecules than the user intended, if, for instance, we run out of
    space, but it will never leave partial complexes lying around.

    N.B. We may need to revisit the failure criteria and the retry behavior
         here.
 *******************************************************************/
static int init_effectors_place_complexes(int n_to_place,
                                          int nwalls,
                                          double *weights,
                                          struct wall * const *walls,
                                          struct region *rp,
                                          struct eff_dat const *effdp)
{
  if (nwalls <= 0)
    return 1;

  double max_weight = weights[nwalls - 1];
  long long n_failures = 0;
  int n_total = n_to_place;
  while (n_to_place > 0)
  {
    int num_tries = world->complex_placement_attempts;
    int chosen_wall = 0;
    double p = rng_dbl(world->rng) * max_weight;

    /* Pick a wall */
    chosen_wall = bisect_high(weights, nwalls, p);

    /* Try to find a spot for the release */
    while (-- num_tries >= 0)
    {
      if (! init_effectors_place_complex(walls[chosen_wall], rp, effdp))
        break;
    }

    if (num_tries >= 0)
    {
      -- n_to_place;
      n_failures = 0;
    }
    else
    {
      if (++ n_failures >= world->notify->complex_placement_failure_threshold)
      {
        switch (world->notify->complex_placement_failure)
        {
          case WARN_COPE:
            break;

          case WARN_WARN:
            mcell_warn("Unable to place some surface complexes of species '%s' (placed %d of %d).\n",
                       effdp->eff->sym->name,
                       n_total - n_to_place,
                       n_total);
            break;

          case WARN_ERROR:
            mcell_error("Error: Unable to place some surface complexes of species '%s' (placed %d of %d).\n",
                        effdp->eff->sym->name,
                        n_total - n_to_place,
                        n_total);


          default: UNHANDLED_CASE(world->notify->complex_placement_failure);
        }
        break;
      }
    }
  }

  return 0;
}

/********************************************************************
 init_complex_effectors:

    Place surface macromolecules on all walls on the specified object.

    In:  int n_to_place - how many to release
         int nwalls - upon how many walls do we release?
         double *weights - array of cumulative areas for the first n walls in
                  the walls array; used for weighted samping of the area
         struct wall * const *walls - the walls upon which to release
         struct region *rp - the region in which to instantiate molecules, or
                    NULL if the molecule is not restricted to a single region
         struct eff_dat const *effdp - description of what to release
    Out: 0 on success, 1 on failure

    Program state is not corrupted if this fails.  This may release fewer
    macromolecules than the user intended, if, for instance, we run out of
    space, but it will never leave partial complexes lying around.

    N.B. We may need to revisit the failure criteria and the retry behavior
         here.
 *******************************************************************/
static int init_complex_effectors(struct object *objp, struct region_list *head)
{
  /* Do not place molecules if we're restoring from a checkpoint */
  if (world->chkpt_init == 0)
    return 0;

  /* Now, handle each region release */
  for (; head != NULL; head = head->next)
  {
    struct wall *w = NULL;                      /* current wall */
    struct region *rp = head->reg;              /* current region */
    double total_area = 0.0;                    /* cumulative wall area so far */
    struct wall *walls[rp->membership->nbits];  /* walls in region */
    double weights[rp->membership->nbits];      /* cumulative area for walls array */
    int num_fill = 0;                           /* num walls in array so far */
    int avail_slots = 0;                        /* num free slots for molecules */

    /* Collect all walls in this region */
    for (int n_wall=0; n_wall<rp->membership->nbits; ++n_wall)
    {
      if (get_bit(rp->membership, n_wall))
      {
        w = objp->wall_p[n_wall];
        if (w->grid == NULL  &&  create_grid(w, NULL))
          mcell_allocfailed("Failed to create grid for wall.");

        walls[num_fill] = w;
        weights[num_fill] = total_area += w->area;
        avail_slots += w->grid->n_tiles - w->grid->n_occupied;
        ++ num_fill;
      }
    }

    /* Process each mol. type to release on this region */
    struct eff_dat *effdp;
    for (effdp = rp->eff_dat_head; effdp != NULL; effdp = effdp->next)
    {
      int n_to_place;

      /* Skip it if it isn't a macromol */
      if ((effdp->eff->flags & IS_COMPLEX) == 0)
        continue;

      /* Compute the total number to distribute in this region */
      if (effdp->quantity_type == EFFNUM)
        n_to_place = (int) (effdp->quantity + 0.5);
      else if (effdp->quantity_type == EFFDENS)
      {
        double dn_to_place = (effdp->quantity * total_area) / world->grid_density;
        n_to_place = (int) dn_to_place;
        if (rng_dbl(world->rng) < (dn_to_place - n_to_place))
          ++ n_to_place;
      }
      else
        mcell_internal_error("Unknown effector quantity type (%d).", effdp->quantity_type);

      /* Place them */
      if (init_effectors_place_complexes(n_to_place, num_fill, weights, walls, rp, effdp))
        return 1;
    }
  }
  return 0;
}

/********************************************************************
 init_effectors_by_density:

    Place surface molecules on the specified wall.  This occurs after placing
    surface macromolecules, but before placing surface molecules by number.
    This is done by computing a per-tile probability, and releasing a molecule
    onto each tile with the appropriate probability.

    In:  struct wall *w - wall upon which to place
         struct eff_dat *effdp - description of what to release
    Out: 0 on success, 1 on failure
 *******************************************************************/
int init_effectors_by_density(struct wall *w, struct eff_dat *effdp_head)
{
  struct object *objp;
  struct species **eff;
  struct surface_grid *sg;
  struct grid_molecule *mol;
  short *orientation;
  unsigned int n_eff_entry;
  unsigned int n_tiles;
  unsigned int n_occupied;
  int num_eff_dat;
  int p_index;
  double rnd,*prob,area,tot_prob,tot_density;
  struct subvolume *gsv = NULL;

  no_printf("Initializing effectors by density...\n");

  if (create_grid(w,NULL))
    mcell_allocfailed("Failed to create grid for wall.");
  sg=w->grid;
  objp=w->parent_object;

  num_eff_dat = 0;
  for (struct eff_dat *effdp = effdp_head; effdp != NULL; effdp = effdp->next)
    ++ num_eff_dat;

  eff         = CHECKED_MALLOC_ARRAY(struct species *, num_eff_dat, "effector-by-density placement array");
  prob        = CHECKED_MALLOC_ARRAY(double,           num_eff_dat, "effector-by-density placement array");
  orientation = CHECKED_MALLOC_ARRAY(short,            num_eff_dat, "effector-by-density placement array");
  memset(eff,         0, num_eff_dat*sizeof(struct species *));
  memset(prob,        0, num_eff_dat*sizeof(double));
  memset(orientation, 0, num_eff_dat*sizeof(short));

  n_tiles = sg->n_tiles;
  area = w->area;
  objp->n_tiles += n_tiles;
  no_printf("Initializing %d effectors...\n", n_tiles);
  no_printf("  Area = %.9g\n", area);
  no_printf("  Grid_size = %d\n", sg->n);
  no_printf("  Number of effector types in wall = %d\n", num_eff_dat);

  n_eff_entry = 0;
  tot_prob=0;
  tot_density=0;
  for (struct eff_dat *effdp = effdp_head; effdp != NULL; effdp = effdp->next)
  {
    no_printf("  Adding effector %s to wall at density %.9g\n",effdp->eff->sym->name,effdp->quantity);
    tot_prob+=(area*effdp->quantity)/(n_tiles*world->grid_density);
    prob[n_eff_entry]=tot_prob;
    if (effdp->orientation > 0) orientation[n_eff_entry] = 1;
    else if (effdp->orientation < 0) orientation[n_eff_entry] = -1;
    else orientation[n_eff_entry] = 0;
    eff[n_eff_entry ++] = effdp->eff;
    tot_density+=effdp->quantity;
  }

  if (tot_density>world->grid_density)
    mcell_warn("Total effector density too high: %f.  Filling all available effector sites.", tot_density);

  n_occupied=0;
  if (world->chkpt_init) {
    for (unsigned int n_tile = 0; n_tile<n_tiles; ++ n_tile)
    {
      if (sg->mol[n_tile] != NULL)
        continue;

      p_index=-1;
      rnd = rng_dbl(world->rng);
      for (int n_eff = 0; n_eff < num_eff_dat; ++ n_eff)
      {
        if (rnd <= prob[n_eff])
        {
          p_index = n_eff;
          break;
        }
      }

      if (p_index == -1)
        continue;

      struct vector2 s_pos;
      struct vector3 pos3d;
      n_occupied++;
      eff[p_index]->population++;

      if (world->randomize_gmol_pos) grid2uv_random(sg, n_tile, &s_pos);
      else grid2uv(sg, n_tile, &s_pos);
      uv2xyz(&s_pos, w, &pos3d);
      gsv = find_subvolume(&pos3d, gsv);
      mol = (struct grid_molecule *) CHECKED_MEM_GET(gsv->local_storage->gmol, "grid molecule");
      sg->mol[n_tile] = mol;
      mol->t=0;
      mol->t2=0;
      mol->birthday = 0;
      mol->cmplx = NULL;
      mol->s_pos.u = s_pos.u;
      mol->s_pos.v = s_pos.v;
      mol->properties=eff[p_index];
      mol->birthplace=w->birthplace->gmol;
      mol->grid_index = n_tile;
      mol->grid=sg;
      mol->orient=orientation[p_index];
      if (mol->orient == 0)
        mol->orient = (rng_uint(world->rng) & 1) ? 1 : -1;

      mol->flags=TYPE_GRID|ACT_NEWBIE|IN_SCHEDULE|IN_SURFACE;
      if (mol->properties->space_step>0) mol->flags |= ACT_DIFFUSE;
      if ( trigger_unimolecular(eff[p_index]->hashval,(struct abstract_molecule *)mol)!=NULL
           || (eff[p_index]->flags&CAN_GRIDWALL)!=0 ) {
        mol->flags|=ACT_REACT;
      }
      if ((mol->properties->flags&COUNT_ENCLOSED) != 0) mol->flags |= COUNT_ME;

      if ((mol->properties->flags & (COUNT_CONTENTS|COUNT_ENCLOSED)) != 0)
        count_region_from_scratch((struct abstract_molecule*)mol,NULL,1,NULL,NULL,mol->t);

      if (schedule_add(gsv->local_storage->timer,mol))
        mcell_allocfailed("Failed to add grid molecule '%s' to scheduler.", mol->properties->sym->name);
    }
  }

  sg->n_occupied=n_occupied;
  objp->n_occupied_tiles+=n_occupied;

#ifdef DEBUG
  for (int n_eff=0; n_eff < num_eff_dat; ++ n_eff)
    no_printf("Total number of effector %s = %d\n", eff[n_eff]->sym->name, eff[n_eff]->population);
#endif

  free(eff);
  free(prob);
  free(orientation);

  no_printf("Done initializing %u effectors by density\n",n_occupied);

  return 0;
}


/********************************************************************
 init_effectors_by_number:

    Place surface molecules on the specified object.  This occurs after placing
    surface macromolecules, and after placing surface molecules by density.

    In:  struct object *objp - object upon which to place
         struct region_list *reg_eff_num_head - list of what to place
    Out: 0 on success, 1 on failure
 *******************************************************************/
int init_effectors_by_number(struct object *objp, struct region_list *reg_eff_num_head)
{
  static struct grid_molecule DUMMY_MOLECULE;
  static struct grid_molecule *bread_crumb = &DUMMY_MOLECULE;

  unsigned int n_free_eff;
  struct subvolume *gsv=NULL;

    no_printf("Initializing effectors by number...\n");

    /* traverse region list and add effector sites by number to whole regions
       as appropriate */
    for (struct region_list *rlp = reg_eff_num_head;
         rlp != NULL;
         rlp = rlp->next)
    {
      struct region *rp=rlp->reg;
        /* initialize effector grids in region as needed and */
        /* count total number of free effector sites in region */
        n_free_eff=0;
        for (int n_wall=0; n_wall<rp->membership->nbits; n_wall++)
        {
          if (get_bit(rp->membership, n_wall))
          {
            struct wall *w=objp->wall_p[n_wall];
            if (create_grid(w,NULL))
              mcell_allocfailed("Failed to allocate grid for wall.");

	    struct surface_grid *sg=w->grid;
	    n_free_eff=n_free_eff+(sg->n_tiles-sg->n_occupied);
          }
        }
        no_printf("Number of free effector tiles in region %s = %d\n",rp->sym->name,n_free_eff);  

        if (n_free_eff == 0) {
          mcell_warn("Number of free effector tiles in region %s = %d", rp->sym->name, n_free_eff);
          continue;
        }
 
      if (world->chkpt_init) {  /* only needed for denovo initiliazation */
        struct grid_molecule ***tiles;
        unsigned int *idx;
        struct wall **walls;

        /* allocate memory to hold array of pointers to all free tiles */
        tiles = CHECKED_MALLOC_ARRAY(struct grid_molecule **, n_free_eff, "effector placement tiles array");
        idx   = CHECKED_MALLOC_ARRAY(unsigned int,            n_free_eff, "effector placement indices array");
        walls = CHECKED_MALLOC_ARRAY(struct wall *,           n_free_eff, "effector placement walls array");

        /* initialize array of pointers to all free tiles */
        int n_slot = 0;
        for (int n_wall=0; n_wall<rp->membership->nbits; n_wall++)
        {
          if (get_bit(rp->membership, n_wall))
          {
            struct wall *w = objp->wall_p[n_wall];
	    struct surface_grid *sg=w->grid;
	    if (sg!=NULL) {
              for (unsigned int n_tile=0; n_tile<sg->n_tiles; n_tile++)
              {
                if (sg->mol[n_tile]==NULL)
                {
                  tiles[n_slot] = &(sg->mol[n_tile]);
                  idx[n_slot] = n_tile;
                  walls[n_slot++] = w;
                }
	      }
	    }
	  }
        }

        /* distribute desired number of effector sites */
        /* for each effector type to add */
        for (struct eff_dat *effdp=rp->eff_dat_head;
             effdp !=NULL;
             effdp = effdp->next)
        {
          if (effdp->quantity_type == EFFNUM) {
            struct species *eff=effdp->eff;
            short orientation;
            unsigned int n_set = effdp->quantity;
            unsigned int n_clear = n_free_eff - n_set;

            /* Compute orientation */
            if (effdp->orientation > 0) orientation = 1;
            else if (effdp->orientation < 0) orientation = -1;
            else orientation = 0;

            /* Clamp n_set to number of available slots (w/ warning). */
            if (n_set > n_free_eff)
            {
              mcell_warn("Number of %s effectors to place (%d) exceeds number of free effector tiles (%d) in region %s[%s].\n"
                         "  Effectors %s placed on all available effector sites.",
                         eff->sym->name,
                         n_set,
                         n_free_eff,
                         rp->parent->sym->name,
                         rp->region_last_name,
                         eff->sym->name);
              n_set = n_free_eff;
              n_clear=0;
            }

            eff->population+=n_set;

            no_printf("distribute %d of effector %s\n",n_set,eff->sym->name);
            no_printf("n_set = %d  n_clear = %d  n_free_eff = %d\n",n_set,n_clear,n_free_eff);

            /* if filling more than half the free tiles
               init all with bread_crumbs
               choose which tiles to free again
               and then convert remaining bread_crumbs to actual molecules */
            if (n_set > n_free_eff/2) {
              no_printf("filling more than half the free tiles: init all with bread_crumb\n");
              for (unsigned int j=0;j<n_free_eff;j++) {
                *tiles[j]=bread_crumb;
              }

              no_printf("choose which tiles to free again\n");
              for (unsigned int j=0;j<n_clear;j++) {

                /* Loop until we find a vacant tile. */
                while (1) {
                  int slot_num = (int) (rng_dbl(world->rng)*n_free_eff);
                  if (*tiles[slot_num]==bread_crumb) {
                    *tiles[slot_num]=NULL;
                    break;
                  }
                }
              }

              no_printf("convert remaining bread_crumbs to actual molecules\n");
              for (unsigned int j=0;j<n_free_eff;j++) {
                if (*tiles[j]==bread_crumb) {
                  struct vector2 s_pos;
                  struct vector3 pos3d;
                  struct grid_molecule *mol;
                  if (world->randomize_gmol_pos) grid2uv_random(walls[j]->grid,idx[j],&s_pos);
                  else grid2uv(walls[j]->grid,idx[j],&s_pos);
                  uv2xyz(&s_pos, walls[j], &pos3d);
                  gsv = find_subvolume(&pos3d, gsv);

                  mol=(struct grid_molecule *) CHECKED_MEM_GET(gsv->local_storage->gmol, "grid molecule");
                  *tiles[j]=mol;
                  mol->t=0;
                  mol->t2=0;
                  mol->birthday=0;
                  mol->properties=eff;
                  mol->birthplace=walls[j]->birthplace->gmol;
                  mol->grid_index=idx[j];
                  mol->s_pos.u = s_pos.u;
                  mol->s_pos.v = s_pos.v;
                  if (orientation == 0)
                    mol->orient = (rng_uint(world->rng)&1) ? 1 : -1;
                  else
                    mol->orient = orientation;
                  mol->cmplx = NULL;
                  mol->grid=walls[j]->grid;
                  mol->flags=TYPE_GRID|ACT_NEWBIE|IN_SCHEDULE|IN_SURFACE;
                  if (mol->properties->space_step > 0) mol->flags |= ACT_DIFFUSE;
                  if (trigger_unimolecular(eff->hashval,(struct abstract_molecule *)mol)!=NULL
                      || (eff->flags&CAN_GRIDWALL)!=0 ) {
                    mol->flags|=ACT_REACT;
                  }
                  if ((mol->properties->flags&COUNT_ENCLOSED) != 0) mol->flags |= COUNT_ME;

                  if ((mol->properties->flags & (COUNT_CONTENTS|COUNT_ENCLOSED)) != 0)
                    count_region_from_scratch((struct abstract_molecule*)mol,NULL,1,NULL,NULL,mol->t);

                  if (schedule_add(gsv->local_storage->timer, mol))
                    mcell_allocfailed("Failed to add volume molecule '%s' to scheduler.", mol->properties->sym->name);
                }
              }
            }
            else {  /* just fill only the tiles we need */
              no_printf("fill only the tiles we need\n");
              for (unsigned int j=0;j<n_set;j++) {

                /* Loop until we find a vacant tile. */
                while (1) {
                  int slot_num = (int) (rng_dbl(world->rng)*n_free_eff);
                  if (*tiles[slot_num]==NULL) {
                    struct vector2 s_pos;
                    struct vector3 pos3d;
                    struct grid_molecule *mol;
                    if (world->randomize_gmol_pos) grid2uv_random(walls[slot_num]->grid,idx[slot_num],&s_pos);
                    else grid2uv(walls[slot_num]->grid,idx[slot_num],&s_pos);
                    uv2xyz(&s_pos, walls[slot_num], &pos3d);
                    gsv = find_subvolume(&pos3d, gsv);

                    mol=(struct grid_molecule *)CHECKED_MEM_GET(gsv->local_storage->gmol, "grid molecule");
                    *tiles[slot_num]=mol;
                    mol->t=0;
                    mol->t2=0;
                    mol->birthday=0;
                    mol->properties=eff;
                    mol->birthplace=walls[slot_num]->birthplace->gmol;
                    mol->grid_index=idx[slot_num];
                    mol->s_pos.u = s_pos.u;
                    mol->s_pos.v = s_pos.v;
                    mol->cmplx = NULL;
                    if (orientation == 0)
                      mol->orient = (rng_uint(world->rng) & 1) ? 1 : -1;
                    else
                      mol->orient = orientation;

                    mol->grid=walls[slot_num]->grid;
                    mol->flags=TYPE_GRID|ACT_NEWBIE|IN_SCHEDULE|IN_SURFACE;
                    if (mol->properties->space_step > 0) mol->flags |= ACT_DIFFUSE;
                    if (trigger_unimolecular(eff->hashval,(struct abstract_molecule *)mol)!=NULL
                        || (eff->flags&CAN_GRIDWALL)!=0) {
                      mol->flags|=ACT_REACT;
                    }

                    if ((mol->properties->flags & (COUNT_CONTENTS|COUNT_ENCLOSED)) != 0)
                      count_region_from_scratch((struct abstract_molecule*)mol,NULL,1,NULL,NULL,mol->t);

                    if (schedule_add(gsv->local_storage->timer, mol))
                      mcell_allocfailed("Failed to add volume molecule '%s' to scheduler.", mol->properties->sym->name);
                    break;
                  }
                }
              }
            }

            if(n_clear > 0)
            { 
              struct grid_molecule ***tiles_tmp;
              unsigned int *idx_tmp;
              struct wall **walls_tmp;

              /* allocate memory to hold array of pointers to remaining free tiles */
              tiles_tmp = CHECKED_MALLOC_ARRAY(struct grid_molecule **, n_clear, "effector placement tiles array");
              idx_tmp   = CHECKED_MALLOC_ARRAY(unsigned int,            n_clear, "effector placement indices array");
              walls_tmp = CHECKED_MALLOC_ARRAY(struct wall *,           n_clear, "effector placement walls array");

              n_slot = 0;
              for (unsigned int n_eff=0; n_eff<n_free_eff; n_eff++) {
                if (*tiles[n_eff] == NULL)
                {
                  tiles_tmp[n_slot] = tiles[n_eff];
                  idx_tmp[n_slot] = idx[n_eff];
                  walls_tmp[n_slot++] = walls[n_eff];
                }
              }
              /* free original array of pointers to all free tiles */
              free(tiles);
              free(idx);
              free(walls);
              tiles=tiles_tmp;
              idx=idx_tmp;
              walls=walls_tmp;
              n_free_eff=n_free_eff-n_set;
            }

            /* update n_occupied for each effector grid */
            for (int n_wall=0; n_wall<rp->membership->nbits; n_wall++)
            {
              if (get_bit(rp->membership, n_wall))
              {
                struct surface_grid *sg = objp->wall_p[n_wall]->grid;
                if (sg!=NULL)
                {
                  sg->n_occupied=0;
                  for (unsigned int n_tile=0; n_tile<sg->n_tiles; ++ n_tile)
                  {
                    if (sg->mol[n_tile] != NULL)
                      sg->n_occupied++;
                  }
                }
              }
            }
          }
        }
        /* free array of pointers to all free tiles */
        if (tiles!=NULL) {
          free(tiles);
        }
        if (idx!=NULL) {
          free(idx);
        }
        if (walls!=NULL) {
          free(walls);
        }
      } /* end if(world->chkpt_init) */
    }
    no_printf("Done initialize effectors by number.\n");
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
static struct void_list* rel_expr_grab_obj(struct release_evaluator *root, struct mem_helper *voidmem)
{
  struct void_list *vl = NULL;
  struct void_list *vr = NULL;
  
  if (root->left != NULL)
  {
    if (root->op&REXP_LEFT_REGION)
    {
      vl = CHECKED_MEM_GET(voidmem, "temporary list for region release");
      if (vl==NULL) return NULL;
      vl->data = ((struct region*)(root->left))->parent;
      vl->next=NULL;
    }
    else vl = rel_expr_grab_obj(root->left,voidmem);
  }
  if (root->right != NULL)
  {
    if (root->op&REXP_RIGHT_REGION)
    {
      vr = CHECKED_MEM_GET(voidmem, "temporary list for region release");
      if (vr==NULL) return NULL;
      vr->data = ((struct region*)(root->right))->parent;
      vr->next=NULL;
    }
    else vr = rel_expr_grab_obj(root->right,voidmem);
  }
  
  if (vl==NULL)
  {
    if (vr==NULL) return NULL;
    return vr;
  }
  else if (vr==NULL)
  {
    return vl;
  }
  else
  {
    struct void_list *vp;
    
    for (vp=vl;vp->next!=NULL;vp=vp->next) {}
    
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
static struct object** find_unique_rev_objects(struct release_evaluator *root, int *n)
{
  struct object **o_array;
  struct void_list *vp,*vq;
  struct mem_helper *voidmem;
  int n_unique;
  
  voidmem = create_mem(sizeof(struct void_list), 1024);
  if (voidmem == NULL)
    mcell_allocfailed("Failed to create temporary list memory pool.");
  
  vp = rel_expr_grab_obj(root,voidmem);
  if (vp==NULL) return NULL;
  
  vp = void_list_sort(vp);
  
  for (n_unique=1,vq=vp ; vq!=NULL && vq->next!=NULL ; vq=vq->next , n_unique++)
  {
    while (vq->data == vq->next->data)
    {
      vq->next = vq->next->next;
      if (vq->next==NULL) break;
    }
  }
  
  if (vq==NULL) n_unique--;
  *n = n_unique;
  
  o_array = CHECKED_MALLOC_ARRAY(struct object *, n_unique, "object array for region release");
  vq = vp;
  for (unsigned int n_obj=0; vq!=NULL; vq = vq->next, ++ n_obj)
    o_array[n_obj] = (struct object*)vq->data;

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
static int eval_rel_region_expr(struct release_evaluator *expr,
                                int n,
                                struct object **objs,
                                struct bit_array **result,
                                int *n_refinements)
{
  char bit_op;
  
  if (expr->left!=NULL)
  {
    if (expr->op&REXP_INCLUSION)
    {
      if (expr->right==NULL) return 1;  /* Should always have two arguments */
      if (eval_rel_region_expr(expr->left,n,objs,result,n_refinements))
        return 1;
      *n_refinements += 1;
      /* Just ignore right-hand argument; we'll mark that we should look at it later */
    }
    else
    {
      if (expr->op&REXP_LEFT_REGION)
      {
        int pos = void_array_search((void**)objs, (int)n,((struct region*)(expr->left))->parent);
        result[pos] = duplicate_bit_array( ((struct region*)(expr->left))->membership );
        if (result[pos]==NULL) return 1;
      }
      else
      {
        if (eval_rel_region_expr(expr->left,n,objs,result,n_refinements))
          return 1;
      }
      
      if (expr->right==NULL)
      {
	if (expr->op&REXP_NO_OP) return 0;
	else return 1;
      }
      
      if (expr->op&REXP_RIGHT_REGION)
      {
        int pos = void_array_search((void**)objs, (int)n,((struct region*)(expr->right))->parent);
        if (result[pos] == NULL)
        {
          result[pos] = duplicate_bit_array( ((struct region*)(expr->right))->membership );
          if (result[pos]==NULL) return 1;
        }
	else
	{
	  if (expr->op&REXP_UNION) bit_op = '|';
	  else if (expr->op&REXP_SUBTRACTION) bit_op = '-';
	  else if (expr->op&REXP_INTERSECTION) bit_op = '&';
	  else return 1;
  
          bit_operation(result[pos],((struct region*)(expr->right))->membership,bit_op);
	}
      }
      else
      {
	struct bit_array *res2[n];
        for (int i=0;i<n;i++) res2[i]=NULL;

        if (eval_rel_region_expr(expr->right,n,objs,res2,n_refinements))
          return 1;

        for (int i=0;i<n;i++)
	{
	  if (res2[i]==NULL) continue;
	  if (result[i]==NULL) result[i] = res2[i];
	  else
	  {
	    if (expr->op&REXP_UNION) bit_op = '|';
	    else if (expr->op&REXP_SUBTRACTION) bit_op = '-';
	    else if (expr->op&REXP_INTERSECTION) bit_op = '&';
	    else return 1;
	    
	    bit_operation(result[i],res2[i],bit_op);
	    free_bit_array(res2[i]);
	  }
	}
      }
    }
  }
  else return 1;  /* Left should always have something! */
  
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
                                   struct release_region_data *rrd)
{
  rrd->owners = find_unique_rev_objects(rrd->expression , &(rrd->n_objects));
  if (rrd->owners == NULL)
    mcell_error("No objects were found matching the 2-D region release request for release site '%s'.",
                rsop->name);
  
  rrd->in_release = CHECKED_MALLOC_ARRAY(struct bit_array *,
                                         rrd->n_objects,
                                         "region membership array for 2D region release");
  for (int n_object=0; n_object<rrd->n_objects; ++ n_object)
    rrd->in_release[n_object] = NULL;
  
  rrd->refinement = 0;
  if (eval_rel_region_expr(rrd->expression,rrd->n_objects,rrd->owners,rrd->in_release,&rrd->refinement))
    mcell_error("Could not evaluate region expression for release site '%s'.", rsop->name);

  for (int n_object=0; n_object<rrd->n_objects; n_object++)
  {
    if (rrd->owners[n_object] == NULL)
      mcell_internal_error("Object %d of %d in region expression for release site '%s' was not found!",
                           n_object+1,
                           rrd->n_objects,
                           rsop->name);
  }
  
  rrd->walls_per_obj = CHECKED_MALLOC_ARRAY(int,
                                            rrd->n_objects,
                                            "wall counts for 2D region release");
  
  rrd->n_walls_included=0;
  for (int n_object=0; n_object<rrd->n_objects; ++ n_object)
  {
    if (rrd->in_release[n_object] == NULL) rrd->walls_per_obj[n_object]=0;
    else rrd->walls_per_obj[n_object] = count_bits(rrd->in_release[n_object]);
    rrd->n_walls_included += rrd->walls_per_obj[n_object];
  }
  
  rrd->cum_area_list = CHECKED_MALLOC_ARRAY(double,
                                            rrd->n_walls_included,
                                            "cumulative area list for 2D region release");
  rrd->wall_index = CHECKED_MALLOC_ARRAY(int,
                                         rrd->n_walls_included,
                                         "wall indices for 2D region release");
  rrd->obj_index = CHECKED_MALLOC_ARRAY(int,
                                        rrd->n_walls_included,
                                        "object indices for 2D region release");

  unsigned int n_wall_overall = 0;
  for (int n_object=0; n_object<rrd->n_objects; ++ n_object)
  {
    if (rrd->walls_per_obj[n_object]==0) continue;
    int owner_type = rrd->owners[n_object]->object_type;
    if (owner_type != POLY_OBJ && owner_type != BOX_OBJ)
      mcell_internal_error("Found a region on an object which is neither a box nor a polygon (type=%d).", owner_type);

    struct polygon_object *po = (struct polygon_object*)
          (rrd->owners[n_object]->contents);
    const unsigned int n_walls = po->n_walls;
    for (unsigned int n_wall=0; n_wall<n_walls; ++ n_wall)
    {
      if (get_bit(rrd->in_release[n_object], n_wall))
      {
        rrd->cum_area_list[n_wall_overall] = rrd->owners[n_object]->wall_p[n_wall]->area;
        rrd->obj_index[n_wall_overall] = n_object;
        rrd->wall_index[n_wall_overall] = n_wall;
        ++ n_wall_overall;
      }
    }
  }

  for (int n_wall=1; n_wall<rrd->n_walls_included; n_wall++)
  {
    rrd->cum_area_list[n_wall] += rrd->cum_area_list[n_wall-1];
  }
  
  return 0;
}


/***************************************************************************
create_region_bbox:
  In: a region
  Out: pointer to a 2-element array contining the LLF and URB corners of
       a bounding box around the region, or NULL if out of memory.
***************************************************************************/
static struct vector3* create_region_bbox(struct region *r)
{
  struct vector3 *bbox =
        CHECKED_MALLOC_ARRAY(struct vector3, 2, "region bounding box");

  int found_first_wall = 0;
  for (int n_wall=0; n_wall<r->membership->nbits; ++ n_wall)
  {
    if (get_bit(r->membership, n_wall))
    {
      if (! found_first_wall)
      {
        bbox[0].x = bbox[1].x = r->parent->wall_p[n_wall]->vert[0]->x;
        bbox[0].y = bbox[1].y = r->parent->wall_p[n_wall]->vert[0]->y;
        bbox[0].z = bbox[1].z = r->parent->wall_p[n_wall]->vert[0]->z;
        found_first_wall = 1;
      }
      for (unsigned int n_vert=0; n_vert<3; ++ n_vert)
      {
        struct vector3 *v = r->parent->wall_p[n_wall]->vert[n_vert];
        if (bbox[0].x > v->x) bbox[0].x = v->x;
        else if (bbox[1].x < v->x) bbox[1].x = v->x;
        if (bbox[0].y > v->y) bbox[0].y = v->y;
        else if (bbox[1].y < v->y) bbox[1].y = v->y;
        if (bbox[0].z > v->z) bbox[0].z = v->z;
        else if (bbox[1].z < v->z) bbox[1].z = v->z;
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
                                struct vector3 *llf,
                                struct vector3 *urb)
{
  struct region *r;
  
  if (expr->left!=NULL)
  {
    if (expr->op&REXP_LEFT_REGION)
    {
      r = (struct region*)(expr->left);
      if (r->manifold_flag==MANIFOLD_UNCHECKED)
      {
        if (is_manifold(r)) r->manifold_flag = IS_MANIFOLD;
        else
          mcell_error("Cannot release a 3D molecule inside the unclosed region '%s'.", r->sym->name);
      }
      
      if (r->bbox==NULL)
        r->bbox = create_region_bbox(r);
      
      llf->x = r->bbox[0].x;
      llf->y = r->bbox[0].y;
      llf->z = r->bbox[0].z;
      urb->x = r->bbox[1].x;
      urb->y = r->bbox[1].y;
      urb->z = r->bbox[1].z;
    }
    else
    {
      if (eval_rel_region_bbox(expr->left,llf,urb))
        return 1;
    }
    
    if (expr->right==NULL)
    {
      if (expr->op&REXP_NO_OP) return 0;
      else
        mcell_internal_error("Right subtree of release expression is unexpectedly NULL.");
    }
    
    if (expr->op&REXP_SUBTRACTION) return 0;
    else
    {
      struct vector3 llf2;
      struct vector3 urb2;
      
      if (expr->op&REXP_RIGHT_REGION)
      {
        r = (struct region*)(expr->right);
        if (r->manifold_flag==MANIFOLD_UNCHECKED)
        {
          if (is_manifold(r)) r->manifold_flag = IS_MANIFOLD;
          else
            mcell_error("Cannot release a 3D molecule inside the unclosed region '%s'.", r->sym->name);
        }
        
        if (r->bbox==NULL)
          r->bbox = create_region_bbox(r);

        llf2.x = r->bbox[0].x;
        llf2.y = r->bbox[0].y;
        llf2.z = r->bbox[0].z;
        urb2.x = r->bbox[1].x;
        urb2.y = r->bbox[1].y;
        urb2.z = r->bbox[1].z;
      }
      else
      {
        if (eval_rel_region_bbox(expr->right,&llf2,&urb2))
          return 1;
      }
      
      if (expr->op&REXP_UNION)
      {
        if (llf->x > llf2.x) llf->x = llf2.x;
        if (llf->y > llf2.y) llf->y = llf2.y;
        if (llf->z > llf2.z) llf->z = llf2.z;
        if (urb->x < urb2.x) urb->x = urb2.x;
        if (urb->y < urb2.y) urb->y = urb2.y;
        if (urb->z < urb2.z) urb->z = urb2.z;
      }
      else if (expr->op&(REXP_INTERSECTION|REXP_INCLUSION))
      {
        if (llf->x < llf2.x) llf->x = llf2.x;
        if (llf->y < llf2.y) llf->y = llf2.y;
        if (llf->z < llf2.z) llf->z = llf2.z;
        if (urb->x > urb2.x) urb->x = urb2.x;
        if (urb->y > urb2.y) urb->y = urb2.y;
        if (urb->z > urb2.z) urb->z = urb2.z;
      }
      else
        mcell_internal_error("Release expression contains an unknown or unexpected operator: (%d).",
                             expr->op);
    }
  }
  else
    mcell_internal_error("Left subtree of release expression is unexpectedly NULL.");
  
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
static int init_rel_region_data_3d(struct release_region_data *rrd)
{
  rrd->n_walls_included = 0;
  
  if (eval_rel_region_bbox(rrd->expression,&(rrd->llf),&(rrd->urb)))
    return 1;

  if (rrd->llf.x >= rrd->urb.x ||
      rrd->llf.y >= rrd->urb.y ||
      rrd->llf.z >= rrd->urb.z)
  {
    return -1;  /* Special signal to print out "nothing in here" error msg */
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
static void output_relreg_eval_tree(FILE *f,
                                    char *prefix,
                                    char cA,
                                    char cB,
                                    struct release_evaluator *expr)
{
  size_t l = strlen(prefix);
  char my_op;
  
  if (expr->op&REXP_NO_OP)
  {
    fprintf(f,"%s >%s\n",prefix,((struct region*)(expr->left))->sym->name);
  }
  else
  {
    char prefixA[l+3];
    char prefixB[l+3];
    strncpy(prefixA,prefix,l);
    strncpy(prefixB,prefix,l);
    prefixA[l] = cA;
    prefixB[l] = cB;
    prefixA[l+1] = prefixB[l+1] = ' ';
    prefixA[l+2] = prefixB[l+2] = 0;
    
    if (expr->op&REXP_LEFT_REGION)
    {
      fprintf(f,"%s >%s\n",prefix,((struct region*)(expr->left))->sym->name);
    }
    else
    {
      output_relreg_eval_tree(f,prefixA,' ','|',expr->left);
    }
    
    my_op = '?';
    if (expr->op & REXP_UNION) my_op = '+';
    else if (expr->op & REXP_INTERSECTION) my_op = '*';
    else if (expr->op & REXP_SUBTRACTION) my_op = '-';
    else if (expr->op & REXP_INCLUSION) my_op = ':';
    
    fprintf(f,"%s%c\n",prefix,my_op);
    
    if (expr->op&REXP_RIGHT_REGION)
    {
      fprintf(f,"%s >%s\n",prefix,((struct region*)(expr->left))->sym->name);
    }
    else
    {
      output_relreg_eval_tree(f,prefixA,'|',' ',expr->right);
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

int init_releases(void)
{
  struct release_event_queue *req;
  struct abstract_element *ae;
  struct schedule_helper *sh;
  int i;
  
  for (sh=world->releaser ; sh!=NULL ; sh=sh->next_scale)
  {
    for (i=-1;i<sh->buf_len;i++)
    {
      for ( ae = (i==-1)?sh->current:sh->circ_buf_head[i] ; ae!=NULL ; ae=ae->next )
      {
        req = (struct release_event_queue*)ae;
        switch ((int) req->release_site->release_shape)
        {
          case SHAPE_REGION:
            if (req->release_site->mol_type == NULL)
              mcell_error("Molecule type was not specified for the region release site '%s'.",
                          req->release_site->name);
            if ((req->release_site->mol_type->flags & NOT_FREE) == 0)
            {
              switch (init_rel_region_data_3d(req->release_site->region_data))
              {
                case 0:
                  break;

                case -1:
                  mcell_warn("Region release site '%s' is empty!  Ignoring!  Evaluation tree:\n",
                             req->release_site->name);
                  output_relreg_eval_tree(mcell_get_error_file()," ",' ',' ',req->release_site->region_data->expression);
                  req->release_site->release_number_method=CONSTNUM;
                  req->release_site->release_number=0;
                  break;

                default:
                  mcell_error("Unexpected error while initializing 3-D region releases for release site '%s'.",
                              req->release_site->name);
                  break;
              }
            }
            else
            {
              if (init_rel_region_data_2d(req->release_site, req->release_site->region_data))
                mcell_error("Unexpected error while initializing 2-D region releases for release site '%s'.",
                            req->release_site->name);
            }
            break;

          case SHAPE_LIST:
            if (req->release_site->mol_list == NULL)
              mcell_error("Molecule positions for the LIST release site '%s' are not specified.", req->release_site->name);
            break;

          case SHAPE_SPHERICAL:
          case SHAPE_CUBIC:
          case SHAPE_ELLIPTIC:
          case SHAPE_RECTANGULAR:
          case SHAPE_SPHERICAL_SHELL:
            /* geometrical release sites */
            if (req->release_site->mol_type == NULL)
              mcell_error("Molecule type for the release site '%s' is not specified.", req->release_site->name);
            if (req->release_site->diameter == NULL)
              mcell_error("Diameter for the geometrical shape release site '%s' is not specified.", req->release_site->name);
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
  In: species
  Out: None.  If species is a surface (surface class) and special reactions 
       like TRANSPARENT, REFLECTIVE or ABSORPTIVE are defined for it,
       the reactions report is printed out.
***************************************************************************/
void publish_special_reactions_report(struct species *sp)
{
   struct name_orient *no;
   FILE *log_file;
   struct species *spec;
   int i;
   /* orientation of GENERIC_MOLECULE */
   int generic_mol_orient;
   /* flags */
   int refl_mols_generic_mol = 0;
   int transp_mols_generic_mol = 0;
   int absorb_mols_generic_mol = 0;

   /* Below I will employ following set of rules for printing out
      the relative orientations of surface classes and molecules;
      1. The orientation of surface class is always printed out as {1}.
      2. The orientation of molecule is printed out 
         as {1} when it is positive, {-1} when it is negative, 
         and {0} when it is zero, or absent.
   */

   struct name_list *nl_head = NULL, *nl, *nnext;
 
   /* create name list of 3D species */
   for(i = 0; i < world->n_species; i++)
   {
      spec = world->species_list[i];
      if ((spec == world->g_mol)  ||  (spec == world->g_surf)) continue;
      if(spec->flags & ON_GRID) continue;
      if(spec->flags & IS_SURFACE) continue;

      nl = CHECKED_MALLOC_STRUCT(struct name_list, "name_list");
      nl->name = my_strcat(spec->sym->name, NULL);
      nl->prev = NULL; /* we will use only FORWARD feature */

      if(nl_head == NULL)
      {
         nl->next = NULL;
         nl_head = nl;
      }else{
         nl->next = nl_head;
         nl_head = nl;
      }
   }


   log_file = mcell_get_log_file();

   if(sp->refl_mols != NULL)
   {
      fprintf(log_file, "Surfaces with surface class \"%s{1}\" are REFLECTIVE for molecules  ", sp->sym->name);
      /* search for GENERIC_MOLECULE */
      for(no = sp->refl_mols; no != NULL; no = no->next)
      {
        if(strcmp(no->name, "GENERIC_MOLECULE") == 0)
        {
           generic_mol_orient = no->orient;
           refl_mols_generic_mol = 1;
           break;
        }
      }

      if(refl_mols_generic_mol)
      {
         for(nl = nl_head; nl != NULL; nl = nl->next)
         {
            fprintf(log_file, "%s{%d}", nl->name, generic_mol_orient);
            if(nl->next != NULL) fprintf(log_file, ", ");
            else fprintf(log_file, ".");
         }
      }else{ 
         for(no = sp->refl_mols; no != NULL; no = no->next)
         {
            fprintf(log_file, "%s{%d}", no->name, no->orient);
            if(no->next != NULL) fprintf(log_file, ", ");
            else fprintf(log_file, ".");
         }
      }
      fprintf(log_file, "\n");
   }
  
   if(sp->transp_mols != NULL)
   {
      fprintf(log_file, "Surfaces with surface class \"%s{1}\" are TRANSPARENT for molecules  ", sp->sym->name);
      /* search for GENERIC_MOLECULE */
      for(no = sp->transp_mols; no != NULL; no = no->next)
      {
        if(strcmp(no->name, "GENERIC_MOLECULE") == 0)
        {
           generic_mol_orient = no->orient;
           transp_mols_generic_mol = 1;
           break;
        }
      }

      if(transp_mols_generic_mol)
      {
         for(nl = nl_head; nl != NULL; nl = nl->next)
         {
            fprintf(log_file, "%s{%d}", nl->name, generic_mol_orient);
            if(nl->next != NULL) fprintf(log_file, ", ");
            else fprintf(log_file, ".");
         }
      }else{ 
         for(no = sp->transp_mols; no != NULL; no = no->next)
         {
            fprintf(log_file, "%s{%d}", no->name, no->orient);
            if(no->next != NULL) fprintf(log_file, ", ");
            else fprintf(log_file, ".");
         }
      }
      fprintf(log_file, "\n");
   }

   if(sp->absorb_mols != NULL)
   {
      fprintf(log_file, "Surfaces with surface class \"%s{1}\" are ABSORPTIVE for molecules  ", sp->sym->name);
      /* search for GENERIC_MOLECULE */
      for(no = sp->absorb_mols; no != NULL; no = no->next)
      {
        if(strcmp(no->name, "GENERIC_MOLECULE") == 0)
        {
           generic_mol_orient = no->orient;
           absorb_mols_generic_mol = 1;
           break;
        }
      }

      if(absorb_mols_generic_mol)
      {
         for(nl = nl_head; nl != NULL; nl = nl->next)
         {
            fprintf(log_file, "%s{%d}", nl->name, generic_mol_orient);
            if(nl->next != NULL) fprintf(log_file, ", ");
            else fprintf(log_file, ".");
         }
      }else{ 
         for(no = sp->absorb_mols; no != NULL; no = no->next)
         {
            fprintf(log_file, "%s{%d}", no->name, no->orient);
            if(no->next != NULL) fprintf(log_file, ", ");
            else fprintf(log_file, ".");
         }
      }
      fprintf(log_file, "\n");
   }

   if((sp->refl_mols != NULL) || (sp->transp_mols != NULL) || (sp->absorb_mols != NULL))
   {
      fprintf(log_file, "\n");
   }

   /* remove name list */
   while(nl_head != NULL)
   {
     nnext = nl_head->next;
     free(nl_head);
     nl_head = nnext;
   }
}
