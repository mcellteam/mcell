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

#ifdef DEBUG
#define no_printf printf
#endif

extern struct volume *world;

/* Initialize the surface macromolecules on a given object */
static int init_complex_effectors(struct object *objp, struct region_list *head);

#define MICROSEC_PER_YEAR 365.25*86400.0*1e6

/* Sets default notification values */
int init_notifications()
{
  world->notify = (struct notifications*)malloc(sizeof(struct notifications));
  if (world->notify==NULL) {
     fprintf(world->err_file, "File '%s', Line  %ld: Out of memory.\n", __FILE__, (long)__LINE__);
     return 1;
  }
 
  /* Notifications */
  world->notify->progress_report = NOTIFY_FULL;
  world->notify->diffusion_constants = NOTIFY_BRIEF;
  world->notify->reaction_probabilities = NOTIFY_FULL;
  world->notify->time_varying_reactions = NOTIFY_FULL;
  world->notify->reaction_prob_notify = 0.0;
  world->notify->partition_location = NOTIFY_NONE;
  world->notify->box_triangulation = NOTIFY_NONE;
  world->notify->custom_iterations = NOTIFY_FULL;
  world->notify->custom_iteration_value = 0;  /* Ignored unless NOTIFY_CUSTOM set */
  world->notify->release_events = NOTIFY_FULL;
  world->notify->file_writes = NOTIFY_NONE;
  world->notify->final_summary = NOTIFY_FULL;
  world->notify->throughput_report = NOTIFY_FULL;
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
  
  if (world->log_freq!=-1) /* User set this */
  {
    world->notify->custom_iterations = NOTIFY_CUSTOM;
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
  if (wrld->volume_output_scheduler == NULL) {
    fprintf(wrld->err_file,"File '%s', Line %ld: Out of memory while creating volume_output_scheduler.\n", __FILE__, (long)__LINE__);
    exit(EXIT_FAILURE);
  }

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

        vo->t = vo->times[idx] * time_scale;
        vo->next_time = vo->times + idx;
      }

      /* Advance the next_time pointer */
      ++ vo->next_time;
    }

    if (schedule_add(world->volume_output_scheduler, vo)) {
      fprintf(world->err_file,"File %s, Line %ld: Out of memory while setting up volume output.\n", __FILE__, (long)__LINE__);
      exit(EXIT_FAILURE);
    }
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
  FILE *log_file, *file;
  struct sym_table *gp;
  struct output_block *obp,*obpn;
  struct output_set *set;
  int i;
  double f;
  int reactants_3D_present = 0; /* flag to check whether there are 3D reactants
                             (participants in the reactions
                              between 3D molecules) in the simulation */

  log_file=world->log_file;
#ifdef KELP
  if (world->procnum == 0) {
#endif
  fprintf(log_file,"MCell initializing simulation...\n");
  fflush(log_file);
#ifdef KELP
  }
#endif


  /* Initialize variables to reasonably safe values */

  world->curr_file=world->mdl_infile_name;

  /* by Erhan Gokcay 5/3/2002 ========================== */
  /* We can not initialize chkpt_infile anymore. It is set in mcell.c */
  /*  chkpt_infile=NULL; */
  /* =================================================== */

  world->chkpt_outfile=NULL;
  world->chkpt_iterations=0;
  world->chkpt_seq_num=0;

  /*world->chkpt_init=1; */  /* set in the main() */
  world->chkpt_flag=0;
  world->molecule_prefix_name=NULL;
  world->file_prefix_name=NULL;
  world->random_number_use=0;
  world->ray_voxel_tests=0;
  world->ray_polygon_tests=0;
  world->ray_polygon_colls=0;
  world->mol_mol_colls=0;
  world->mol_mol_mol_colls=0;
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
  
  world->viz_output_flag = 0; 
  world->use_expanded_list=1;
  world->randomize_gmol_pos=1;
  world->vacancy_search_dist2=0;
  world->surface_reversibility=0;
  world->volume_reversibility=0;
  world->n_reactions = 0;

  memset(&world->viz_state_info, 0, sizeof(world->viz_state_info));
  
  world->mcell_version = mcell_version;
  
  world->clamp_list = NULL;

  world->rng = malloc(sizeof(struct rng_state));
  if (world->rng==NULL)
  {
    fprintf(world->err_file,"File '%s', Line %ld: Out of memory, failed to allocate random number generator\n", __FILE__, (long)__LINE__);
    exit(EXIT_FAILURE);
  }
  if (world->seed_seq < 1 || world->seed_seq > INT_MAX) {
    fprintf(world->err_file,"File '%s', Line %ld: error, random sequence number not in range 1 to 2^31-1\n", __FILE__, (long)__LINE__);
    return(1);
  }
  rng_init(world->rng,world->seed_seq);
  fprintf(log_file,"MCell[%d]: random sequence %d\n",world->procnum,world->seed_seq);
  fflush(log_file);

  world->count_hashmask = COUNT_HASHMASK;
  world->count_hash = (struct counter**)malloc(sizeof(struct counter*)*(world->count_hashmask+1));
  if (world->count_hash == NULL)
  {
    fprintf(world->err_file,"File '%s', Line %ld: Out of memory while creating counter hash table\n", __FILE__, (long)__LINE__);
    exit(EXIT_FAILURE);
  }
  for (i=0;i<=world->count_hashmask;i++) world->count_hash[i] = NULL;
  
  world->oexpr_mem = create_mem_named(sizeof(struct output_expression),128,"output expression");
  if (world->oexpr_mem==NULL)
  {
    fprintf(world->err_file,"Out of memory while getting ready to store output expressions\n");
    return 1;
  }
  world->outp_request_mem = create_mem_named(sizeof(struct output_request),64,"output request");
  if (world->outp_request_mem==NULL)
  {
    fprintf(world->err_file,"Out of memory while getting ready to store lists of output commands\n");
    return 1;
  }
  world->counter_mem = create_mem_named(sizeof(struct counter),32,"counter");
  if (world->counter_mem==NULL)
  {
    fprintf(world->err_file,"Out of memory while getting ready to store reaction and molecule counts\n");
    return 1;
  }
  world->trig_request_mem = create_mem_named(sizeof(struct trigger_request),32,"trigger request");
  if (world->trig_request_mem==NULL)
  {
    fprintf(world->err_file,"Out of memory while getting ready to store output triggers\n");
    return 1;
  }
  world->magic_mem = create_mem_named(sizeof(struct magic_list),1024,"reaction-triggered release");
  if (world->magic_mem==NULL)
  {
    fprintf(world->err_file,"Out of memory while getting ready to store reaction-release list.\n");
    return 1;
  }

  if((world->main_sym_table=init_symtab(SYM_HASHSIZE)) == NULL){
    fprintf(world->err_file,"File '%s', Line %ld: initialization of symbol table failed\n", __FILE__, (long)__LINE__);
    return(1);
  }
	

  if ((gp=store_sym("WORLD_OBJ",OBJ,world->main_sym_table, NULL))==NULL) {
    fprintf(world->err_file,"File '%s', Line %ld: Out of memory while creating world root object\n", __FILE__, (long)__LINE__);
    return(1);
  }
  world->root_object=(struct object *)gp->value;
  world->root_object->object_type=META_OBJ;
  world->root_object->last_name="";

  if ((gp=store_sym("WORLD_INSTANCE",OBJ,world->main_sym_table, NULL))==NULL) {
    fprintf(world->err_file,"File '%s', Line %ld: Out of memory while creating world root instance.\n", __FILE__, (long)__LINE__);
    return(1);
  }
  world->root_instance=(struct object *)gp->value;
  world->root_instance->object_type=META_OBJ;
  world->root_instance->last_name="";

  if ((gp=store_sym("DEFAULT_RELEASE_PATTERN",RPAT,world->main_sym_table, NULL))
      ==NULL) {
    fprintf(world->err_file,"File '%s', Line %ld: Out of memory while creating default release pattern.\n", __FILE__, (long)__LINE__);
    return(1);
  }
  world->default_release_pattern=(struct release_pattern *)gp->value;
  world->default_release_pattern->delay=0;
  world->default_release_pattern->release_interval=FOREVER;
  world->default_release_pattern->train_interval=FOREVER;
  world->default_release_pattern->train_duration=FOREVER;
  world->default_release_pattern->number_of_trains=1;
   
  if ((gp=store_sym("GENERIC_MOLECULE",MOL,world->main_sym_table, NULL))
      ==NULL) {
    fprintf(world->err_file,"File '%s', Line %ld: Out of memory while creating generic molecule.\n", __FILE__, (long)__LINE__);
    return(1);
  }
  world->g_mol=(struct species *)gp->value;

  if ((gp=store_sym("GENERIC_SURFACE",MOL,world->main_sym_table, NULL))
      ==NULL) {
    fprintf(world->err_file,"File '%s', Line %ld: Out of memory while creating generic surface", __FILE__, (long)__LINE__);
    return(1);
  }
  world->g_surf=(struct species *)gp->value;
  world->g_surf->flags=IS_SURFACE;

  world->volume_output_head = NULL;

  world->output_block_head=NULL;
  world->output_request_head=NULL;
  world->viz_obj_head=NULL;
  world->viz_mode=-1;
  world->rk_mode_var=NULL;
  world->frame_data_head=NULL;

  world->releaser = create_scheduler(1.0,100.0,100,0.0);
  if(world->releaser == NULL){
	fprintf(world->err_file, "File '%s', Line %ld: Out of memory while creating releaser.\n", __FILE__, (long)__LINE__);
        exit(EXIT_FAILURE);
  }


  /* Parse the MDL file: */
  no_printf("Node %d parsing MDL file %s\n",world->procnum,world->mdl_infile_name);
  fflush(world->err_file);
  if (mdlparse_init(world)) {
    return(1);
  }
  no_printf("Done parsing MDL file: %s\n",world->mdl_infile_name);
  fflush(world->err_file);

  if(world->iterations == INT_MIN){
     fprintf( world->err_file, "Error: Total number of iterations is not specified either through the ITERATIONS keyword or through the command line option '-iterations'.\n");
     return 1;
  }

  /* Set up the array of species */
  if (init_species())
  {
    fprintf(world->err_file,"File '%s', Line %ld: error initializing species.\n", __FILE__, (long)__LINE__);
    return(1);
  }
  no_printf("Done setting up species.\n");


  /* Visualize all molecules if asked in "mdl" file */
  if((world->viz_mode == DREAMM_V3_MODE) || (world->viz_mode == DREAMM_V3_GROUPED_MODE))
  {
    if((world->viz_output_flag & VIZ_ALL_MOLECULES) != 0) {
       struct species *sp;
  
      for(i = 0; i < world->n_species; i++)
      {
         sp = world->species_list[i];
         if((sp->flags & IS_SURFACE) != 0) continue;
         if(strcmp(sp->sym->name, "GENERIC_MOLECULE") == 0) continue;  

         /* set viz_state to INCLUDE_OBJ for the molecule we want to visualize
             but will not assign state value */
         sp->viz_state = INCLUDE_OBJ;

      } 

    }
  }

 /* If there are no 3D molecules-reactants in the simulation
    set up the"use_expanded_list" flag to zero. */
  for(i = 0; i < world->n_species; i++)
  {
                        
        struct	species *sp = world->species_list[i];
        if(strcmp(sp->sym->name, "GENERIC_MOLECULE") == 0) continue;  
        if(strcmp(sp->sym->name, "GENERIC_SURFACE") == 0) continue;  
        
        if((sp->flags & (CAN_MOLMOL|CAN_MOLMOLMOL)) != 0){
		reactants_3D_present = 1;
                break;
        }
                      
  }

  if(reactants_3D_present == 0){
	world->use_expanded_list = 0;
  }

/* Instantiation Pass #1: Initialize the geometry */
  if (init_geom()) {
    fprintf(world->err_file,"File '%s', Line %ld: error initializing geometry.\n", __FILE__, (long)__LINE__);
    return(1);
  }
  
  no_printf("Done setting up geometry.\n");
  
/* Instantiation Pass #2: Partition geometry */
  if (init_partitions()) {
    fprintf(world->err_file,"File '%s', Line %ld: error initializing partitions.\n", __FILE__, (long)__LINE__);
    return(1);
  }
  
  if (distribute_world()) {
    fprintf(world->err_file,"File '%s', Line %ld: error moving geometry to partitions\n", __FILE__, (long)__LINE__);
    return(1);
  }
  
  if (sharpen_world()) {
    fprintf(world->err_file,"File '%s', Line %ld: error adding edges to geometry.\n", __FILE__, (long)__LINE__);
    return(1);
  }

/* Instantiation Pass #3: Initialize regions */
  if (prepare_counters())
  {
    fprintf(world->err_file,"File '%s' line %d: error while preparing count statements.\n",__FILE__,(int)__LINE__);
    return 1;
  }

  if (init_regions())
  {
    fprintf(world->err_file,"File '%s' line %d: error initializing object regions.\n",__FILE__,(int)__LINE__);
    return(1);
  }
  
  if (check_counter_geometry())
  {
    fprintf(world->err_file,"Count request did not have sensible geometry.\n");
    return 1;
  }
  
  if (world->place_waypoints_flag) {
    if (place_waypoints()) {
      fprintf(world->err_file,"File '%s', Line %ld: error storing waypoints.\n", __FILE__, (long)__LINE__);
      return(1);
    }
  }
  
  if (init_effectors())
  {
    fprintf(world->err_file,"File '%s', Line %ld: Error initializing effectors on regions.\n", __FILE__, (long)__LINE__);
    return 1;
  }
  
  if (init_releases())
  {
     if (world->place_waypoints_flag)
     {
         fprintf(world->err_file,"File '%s', Line %ld: Error initializing releases on regions\n", __FILE__, (long)__LINE__);
         return 1;
     }else{
         fprintf(world->err_file,"File '%s', Line %ld: Error initializing releases.\n", __FILE__, (long)__LINE__);
         return 1;
     }
  }

  if (world->chkpt_infile) {
    if ((world->chkpt_infs=fopen(world->chkpt_infile,"rb"))==NULL) {
      world->chkpt_seq_num=1;
    }
    else {
      fprintf(log_file,"MCell: reading from checkpoint file %s\n",world->chkpt_infile);
      if(read_chkpt(world->chkpt_infs)) {
	fprintf(world->err_file,"File '%s', Line %ld: error reading from checkpoint file %s.\n", __FILE__, (long)__LINE__, world->chkpt_infile);
	return(1);
      }
      fclose(world->chkpt_infs);
    }
  }
  else {
    world->chkpt_seq_num=1;
  }

  /* Initialize the frame data for the visualization and reaction output. */
  if (init_frame_data_list(&world->frame_data_head))
  {
    fprintf(world->err_file,"File '%s', Line %ld: Failed to initialize viz output.\n", __FILE__, (long)__LINE__);
    return 1;
  }

  /* Initialize the volume output */
  init_volume_data_output(world);

  world->count_scheduler = create_scheduler(1.0,100.0,100,world->start_time);
  if(world->count_scheduler == NULL){
	fprintf(world->err_file,"File '%s', Line %ld: Out of memory while creating count_scheduler.\n", __FILE__, (long)__LINE__);
        exit(EXIT_FAILURE);
  }

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
            file = fopen(set->outfile_name,"w");
            if (file==NULL)
            {
              fprintf(world->err_file,"Can't open output file %s\n",set->outfile_name);
              return 1;
            }
            fclose(file);
          }
          else if (obp->timer_type==OUTPUT_BY_ITERATION_LIST)
          {
            if(obp->time_now == NULL) continue;
            i = truncate_output_file(set->outfile_name,obp->t);
            if (i)
            {
              fprintf(world->err_file,"Failed to prepare output file %s to receive output\n",set->outfile_name);
              return 1;
            }
          }
          else if (obp->timer_type==OUTPUT_BY_TIME_LIST)
          {
            if(obp->time_now == NULL) continue;
            i = truncate_output_file(set->outfile_name,obp->t*world->time_unit);
            if (i)
            {
              fprintf(world->err_file,"Failed to prepare output file %s to receive output\n",set->outfile_name);
              return 1;
            }
          }
          else
          {
            i = truncate_output_file(set->outfile_name,obp->t*world->time_unit);
            if (i)
            {
              fprintf(world->err_file,"Failed to prepare output file %s to receive output\n",set->outfile_name);
              return 1;
            }
          }
        }
      }
    
    if (schedule_add(world->count_scheduler , obp))
    {
      fprintf(world->err_file,"File %s, Line %ld: Out of memory while setting up output.\n", __FILE__, (long)__LINE__);
      return 1;
    }
    obp = obpn;
  }

  no_printf("Done initializing simulation\n");
  fflush(log_file);
  return(0);
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
  
  for (i=0;i<SYM_HASHSIZE;i++)
  {
    for (gp = world->main_sym_table[i] ; gp != NULL ; gp = gp->next)
    {    
      if (gp->sym_type==MOL) count++;
    }
  }
  
  world->n_species = count;
  if((world->species_list = (struct species**)malloc(sizeof(struct species*)*world->n_species)) == NULL)
  {
	fprintf(world->err_file, "File '%s', Line %ld: Out of memory during species initialization.\n", __FILE__, (long)__LINE__);
        exit(EXIT_FAILURE);
  }
  count = 0;
  for (i=0;i<SYM_HASHSIZE;i++)
  {
    for (gp = world->main_sym_table[i] ; gp != NULL ; gp = gp->next)
    {    
      if (gp->sym_type==MOL)
      {
        s = (struct species*) gp->value;
        world->species_list[count] = s;
        world->species_list[count]->species_id = count;
        world->species_list[count]->chkpt_species_id = UINT_MAX;
/*        world->species_list[count]->hashval &= world->rx_hashsize-1; */
        world->species_list[count]->population = 0;
	world->species_list[count]->n_deceased = 0;
	world->species_list[count]->cum_lifetime = 0;
        if(world->species_list[count]->viz_state < 0){
        	world->species_list[count]->viz_state = EXCLUDE_OBJ;
        }
        if ( (s->flags & NOT_FREE) == 0 )
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
    Out: A freshly allocated storage with initialized memory pools, or NULL if
         memory allocation fails.  Program state remains valid upon failure of
         this function.
 *******************************************************************/
static struct storage *create_storage(int nsubvols)
{
  struct storage *shared_mem = NULL;
  if((shared_mem = (struct storage*)malloc(sizeof(struct storage))) == NULL)
  {
    fprintf(world->err_file,"File '%s', Line %ld: out of memory while initializing partitions.\n", __FILE__, (long)__LINE__);
    return NULL;
  }
  memset(shared_mem, 0, sizeof(struct storage));

  if (nsubvols < 8) nsubvols = 8;
  if (nsubvols > 4096) nsubvols = 4096;
  /* We should tune the algorithm for selecting allocation block sizes.  */
  /* XXX: Round up to power of 2?  Shouldn't matter, I think. */
  if ((shared_mem->list  = create_mem_named(sizeof(struct wall_list),nsubvols,"wall list")) == NULL) goto failure;
  if ((shared_mem->mol   = create_mem_named(sizeof(struct volume_molecule),nsubvols,"vol mol")) == NULL) goto failure;
  if ((shared_mem->gmol  = create_mem_named(sizeof(struct grid_molecule),nsubvols,"grid mol")) == NULL) goto failure;
  if ((shared_mem->face  = create_mem_named(sizeof(struct wall),nsubvols,"wall")) == NULL) goto failure;
  if ((shared_mem->join  = create_mem_named(sizeof(struct edge),nsubvols,"edge")) == NULL) goto failure;
  if ((shared_mem->tree  = create_mem_named(sizeof(struct vertex_tree),nsubvols,"vertex tree")) == NULL) goto failure;
  if ((shared_mem->grids = create_mem_named(sizeof(struct surface_grid),nsubvols,"surface grid")) == NULL) goto failure;
  if ((shared_mem->regl  = create_mem_named(sizeof(struct region_list),nsubvols,"region list")) == NULL) goto failure;
  if ((shared_mem->pslv  = create_mem_named(sizeof(struct per_species_list),32,"per species list")) == NULL) goto failure;
  shared_mem->coll = world->coll_mem;
  shared_mem->sp_coll = world->sp_coll_mem;
  shared_mem->tri_coll = world->tri_coll_mem;
  shared_mem->exdv = world->exdv_mem;

  if (world->chkpt_init)
  {
    if ((shared_mem->timer = create_scheduler(1.0,100.0,100,0.0)) == NULL)
      goto failure;
    shared_mem->current_time = 0.0;
  }

  if (world->time_step_max==0.0) shared_mem->max_timestep = MICROSEC_PER_YEAR;
  else
  {
    if (world->time_step_max < world->time_unit) shared_mem->max_timestep = 1.0;
    else shared_mem->max_timestep = world->time_step_max/world->time_unit;
  }

  return shared_mem;

failure:
  if (shared_mem)
  {
    if (shared_mem->timer) delete_scheduler(shared_mem->timer);
    if (shared_mem->list)  delete_mem(shared_mem->list);
    if (shared_mem->mol)   delete_mem(shared_mem->mol);
    if (shared_mem->gmol)  delete_mem(shared_mem->gmol);
    if (shared_mem->face)  delete_mem(shared_mem->face);
    if (shared_mem->join)  delete_mem(shared_mem->join);
    if (shared_mem->tree)  delete_mem(shared_mem->tree);
    if (shared_mem->grids) delete_mem(shared_mem->grids);
    if (shared_mem->regl)  delete_mem(shared_mem->regl);
    free(shared_mem);
  }

  fprintf(world->err_file,"File '%s', Line %ld: out of memory while initializing partitions.\n", __FILE__, (long)__LINE__);
  return NULL;
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
  if((world->waypoints = (struct waypoint*)malloc(sizeof(struct waypoint*)*world->n_waypoints)) == NULL)
  {
    fprintf(world->err_file,"File '%s', Line %ld: out of memory while initializing partitions.\n", __FILE__, (long)__LINE__);
    exit(EXIT_FAILURE);
  }

  /* Allocate the subvolumes */
  world->n_subvols = (world->nz_parts-1) * (world->ny_parts-1) * (world->nx_parts-1);
  printf("Creating %d subvolumes (%d,%d,%d per axis)\n",world->n_subvols,world->nx_parts-1,world->ny_parts-1,world->nz_parts-1);
  if ((world->subvol = (struct subvolume*)malloc(sizeof(struct subvolume)*world->n_subvols)) == NULL)
  {
        fprintf(world->err_file,"File '%s', Line %ld: Out of memory, trying to save intermediate results.\n", __FILE__, (long)__LINE__);
        int i = emergency_output();
        fprintf(world->err_file,"Fatal error: out of memory while partition initialization.\nAttempt to write intermediate results had %d errors.\n", i);
        exit(EXIT_FAILURE);
  }

  /* Decide how fine-grained to make the memory subdivisions */
  int subdivisions_per_storage = 0;
  char *env = getenv("SUBDIVISIONS_PER_STORAGE");
  if (env != NULL)
    subdivisions_per_storage = atoi(env);
  if (subdivisions_per_storage == 0)
    subdivisions_per_storage = 14;

  /* Allocate the data structures which are shared between storages */
  if ((world->coll_mem  = create_mem_named(sizeof(struct collision),128,"collision")) == NULL)
  {
    fprintf(world->err_file,"File '%s', Line %ld: Out of memory.\n", __FILE__, (long)__LINE__);
    exit(EXIT_FAILURE);
  }
  if ((world->sp_coll_mem  = create_mem_named(sizeof(struct sp_collision),128,"sp collision")) == NULL)
  {
    fprintf(world->err_file,"File '%s', Line %ld: Out of memory.\n", __FILE__, (long)__LINE__);
    exit(EXIT_FAILURE);
  }
  if ((world->tri_coll_mem  = create_mem_named(sizeof(struct tri_collision),128,"tri collision")) == NULL)
  {
    fprintf(world->err_file,"File '%s', Line %ld: Out of memory.\n", __FILE__, (long)__LINE__);
    exit(EXIT_FAILURE);
  }
  if ((world->exdv_mem  = create_mem_named(sizeof(struct exd_vertex),64,"exact disk vertex")) == NULL)
  {
    fprintf(world->err_file,"File '%s', Line %ld: Out of memory.\n", __FILE__, (long)__LINE__);
    exit(EXIT_FAILURE);
  }

  /* How many storage subdivisions along each axis? */
  int nx = (world->nx_parts + (subdivisions_per_storage) - 2) / (subdivisions_per_storage);
  int ny = (world->ny_parts + (subdivisions_per_storage) - 2) / (subdivisions_per_storage);
  int nz = (world->nz_parts + (subdivisions_per_storage) - 2) / (subdivisions_per_storage);

  /* Create memory pool for storages */
  if((world->storage_allocator = create_mem_named(sizeof(struct storage_list),nx*ny*nz,"storage allocator")) == NULL)
  {
	fprintf(world->err_file,"File '%s', Line %ld: Out of memory, trying to save intermediate results.\n", __FILE__, (long)__LINE__);
        int i = emergency_output();
	fprintf(world->err_file,"Fatal error: out of memory while partition initialization.\nAttempt to write intermediate results had %d errors.\n", i);
        exit(EXIT_FAILURE);
  }

  /* Allocate the storages */
  struct storage *shared_mem[nx*ny*nz];
  int cx = 0, cy = 0, cz = 0;
  for (i=0; i<nx*ny*nz; ++i)
  {
    /* Determine the number of subvolumes included in this subdivision */
    int xd = subdivisions_per_storage, yd = subdivisions_per_storage, zd = subdivisions_per_storage;
    if (cx == nx-1)
      xd = (world->nx_parts - 1) % subdivisions_per_storage;
    if (cy == ny-1)
      yd = (world->ny_parts - 1) % subdivisions_per_storage;
    if (cz == nz-1)
      zd = (world->nz_parts - 1) % subdivisions_per_storage;
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
      exit(EXIT_FAILURE);

    /* Add to the storage list */
    struct storage_list *l = (struct storage_list*)mem_get(world->storage_allocator);
    if (l == NULL)
    {
      fprintf(world->err_file,"File '%s', Line %ld: Out of memory.\n", __FILE__, (long)__LINE__);
      exit(EXIT_FAILURE);
    }
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
    int shidx = (i / (subdivisions_per_storage)) + nx * (j / (subdivisions_per_storage) + ny * (k / (subdivisions_per_storage)));
    sv->local_storage = shared_mem[shidx];
  }
  
  return 0;
}

/**
 * Initializes the geometry of the world.
 * Calls instance_obj() to instantiate all physical objects.
 * (Meta objects, box objects, polygon objects and release sites)
 * Populates viz_obj list vizp
 */
int init_geom(void)
{
  FILE *log_file;
  double tm[4][4];
  double vol_infinity;
  
  no_printf("Initializing physical objects\n");
  log_file=world->log_file;
  vol_infinity=sqrt(DBL_MAX)/4;
  world->bb_llf.x=vol_infinity;
  world->bb_llf.y=vol_infinity;
  world->bb_llf.z=vol_infinity;
  world->bb_urb.x=-vol_infinity;
  world->bb_urb.y=-vol_infinity;
  world->bb_urb.z=-vol_infinity;
  init_matrix(tm);
  
  if(compute_bb(world->root_instance,tm,NULL)){
     return 1;
  }
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
    fprintf(log_file,"MCell: world bounding box in microns =\n");
    fprintf(log_file,"         [ %.9g %.9g %.9g ] [ %.9g %.9g %.9g ]\n",
      world->bb_llf.x*world->length_unit,world->bb_llf.y*world->length_unit,
      world->bb_llf.z*world->length_unit,world->bb_urb.x*world->length_unit,
      world->bb_urb.y*world->length_unit,world->bb_urb.z*world->length_unit);
  }

  world->n_walls=world->root_instance->n_walls;
  world->n_verts=world->root_instance->n_verts;
  no_printf("World object contains %d walls and %d vertices\n",
    world->n_walls,world->n_verts);
  
  if (instance_obj(world->root_instance,tm,NULL,NULL)) {
    return(1);
  }

  return(0);
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
int instance_obj(struct object *objp, double (*im)[4], struct viz_obj *vizp, char *sub_name)
{
  FILE *log_file;
  struct object *child_objp;
  double tm[4][4];
  unsigned short l,m,n;
  char *tmp_name;
  int i;

  log_file=world->log_file;
  l=4;
  m=4;
  n=4;
  mult_matrix(objp->t_matrix,im,tm,l,m,n);
  if (vizp==NULL) {
    vizp=objp->viz_obj;
  }

  if (sub_name!=NULL) { 
    if (strcmp(sub_name,"")==0) {
      tmp_name=strdup("");
      if(tmp_name == NULL){
	fprintf(world->err_file, "File '%s', Line %ld: Out of memory, trying to save intermediate results.\n", __FILE__, (long)__LINE__);
        i = emergency_output();
	fprintf(world->err_file,"Fatal error: out of memory while object instantiation.\nAttempt to write intermediate results had %d errors.\n", i);
        exit(EXIT_FAILURE);
      }else{}
    }
    else {
      tmp_name=my_strcat(sub_name,"."); 
      if(tmp_name == NULL){
	fprintf(world->err_file, "File '%s', Line %ld: Out of memory, trying to save intermediate results.\n", __FILE__, (long)__LINE__);
        i = emergency_output();
	fprintf(world->err_file, "Fatal error: out of memory while object instantiation.\nAttempt to write intermediate results had %d errors.\n", i);
        exit(EXIT_FAILURE);
      }else{}
                  
    }
    sub_name=my_strcat(tmp_name,objp->last_name);    
    if(sub_name == NULL){
	fprintf(world->err_file, "File '%s', Line %ld: Out of memory, trying to save intermediate results.\n", __FILE__, (long)__LINE__);
        i = emergency_output();
	fprintf(world->err_file, "Fatal error: out of memory while object instantiation.\nAttempt to write intermediate results had %d errors.\n", i);
        exit(EXIT_FAILURE);
    }else{}
    free((void *)tmp_name);
  }
  else {
    sub_name=strdup(objp->last_name);    
    if(sub_name == NULL){
	fprintf(world->err_file, "File '%s', Line %ld: Out of memory, trying to save intermediate results.\n", __FILE__, (long)__LINE__);
        i = emergency_output();
	fprintf(world->err_file, "Fatal error: out of memory while object instantiation.\nAttempt to write intermediate results had %d errors.\n", i);
        exit(EXIT_FAILURE);
    }else{}
  }
  
  switch (objp->object_type) {
  case META_OBJ:
    no_printf("Meta object %s instanced\n",sub_name);
    fflush(log_file);
    child_objp=objp->first_child;
    while (child_objp!=NULL) {
      if (instance_obj(child_objp,tm,vizp,sub_name)) {
        return(1);
      }
      child_objp=child_objp->next;
    }
    break;
  case REL_SITE_OBJ:
    no_printf("Release site object %s instanced\n",sub_name);
    fflush(log_file);
    if (instance_release_site(objp,tm)) {
      return(1);
    }
    break;
  case BOX_OBJ:
    no_printf("Box object %s instanced\n",sub_name);
    fflush(log_file);
    if (instance_polygon_object(objp,tm,vizp,sub_name)) {
      return(1);
    }
    break;
  case POLY_OBJ:
    no_printf("Polygon list object %s instanced\n",sub_name);
    fflush(log_file);
    if (instance_polygon_object(objp,tm,vizp,sub_name)) {
      return(1);
    }
    break;
  }

  free((void *)sub_name);
  return(0);
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
  FILE *log_file;
  struct release_site_obj *rsop;
  struct release_event_queue *reqp;
  int i,j;

  log_file=world->log_file;
  rsop=(struct release_site_obj *)objp->contents;
  
  no_printf("Instancing release site object %s\n",objp->sym->name);
  fflush(log_file);
  if (rsop->release_prob==MAGIC_PATTERN_PROBABILITY)
  {
    struct magic_list *ml;
    struct rxn_pathname *rxpn;
    
    ml = (struct magic_list*)mem_get(world->magic_mem);
    if (ml==NULL)
    {
      fprintf(world->err_file,"Internal error in source file %s on line %d\n  Out of memory while preparing release site.\n",__FILE__,__LINE__);
      return 1;
    }
    
    ml->data = rsop;
    ml->type = magic_release;
    
    rxpn = (struct rxn_pathname*)rsop->pattern;
    ml->next=rxpn->magic;
    rxpn->magic=ml;
    
    /* Region releases need to be in release queue to get initialized */
    /* Release code itself is smart enough to ignore MAGIC_PATTERNs */
    if (rsop->release_shape==SHAPE_REGION) 
    {
      reqp = (struct release_event_queue*)malloc(sizeof(struct release_event_queue));
      if (reqp==NULL)
      {
        fprintf(world->err_file,"File '%s', Line %ld: Out of memory while instantiating release site.\n", __FILE__, (long)__LINE__);
        return 1;
      }
    
      reqp->release_site=rsop;
      reqp->event_time=0;
      reqp->train_counter=0;
      reqp->train_high_time=0;
      if (schedule_add(world->releaser,reqp))
      {
        fprintf(world->err_file,"File '%s', Line %ld: Out of memory while scheduling releases.\n", __FILE__, (long)__LINE__);
        return(1);
      }
    }
  }
  else
  {
    reqp = (struct release_event_queue*)malloc(sizeof(struct release_event_queue));
    if (reqp==NULL)
    {
      fprintf(world->err_file,"File '%s', Line %ld: Out of memory while instantiating release site.\n", __FILE__, (long)__LINE__);
      return 1;
    }
  
    reqp->release_site=rsop;
    reqp->event_time=rsop->pattern->delay;
    reqp->train_counter=0;
    reqp->train_high_time=rsop->pattern->delay;
    for (i=0;i<4;i++) for (j=0;j<4;j++) reqp->t_matrix[i][j]=im[i][j];
  
    /* Schedule the release event */
    if (schedule_add(world->releaser,reqp))
    {
      fprintf(world->err_file,"File '%s', Line %ld: Out of memory while scheduling releases.\n", __FILE__, (long)__LINE__);
      return(1);
    }
  
    if(rsop->pattern->train_duration > rsop->pattern->train_interval)
    {
      fprintf(world->err_file,"File '%s', Line %ld: Error - Release pattern train duration is greater than train interval\n", __FILE__, (long)__LINE__);
      return 1;
    }
  }
  
  no_printf("Done instancing release site object %s\n",objp->sym->name);

  fflush(log_file);
  return(0);
}



/**
 * Computes the bounding box for the entire simulation world.
 * Does things recursively in a manner similar to instance_obj().
 */
int compute_bb(struct object *objp, double (*im)[4], char *sub_name)
{
  FILE *log_file;
  struct object *child_objp;
  double tm[4][4];
  unsigned short l,m,n;
  char *tmp_name;
  int i;

  log_file=world->log_file;
  l=4;
  m=4;
  n=4;
  mult_matrix(objp->t_matrix,im,tm,l,m,n);

  if (sub_name!=NULL) { 
    if (strcmp(sub_name,"")==0) {
      tmp_name=strdup("");
      if(tmp_name == NULL){
		fprintf(world->err_file,"File '%s', Line %ld: Out of memory, trying to save intermediate results.\n", __FILE__, (long)__LINE__);
        	i = emergency_output();
		fprintf(world->err_file, "Fatal error: out of memory during bounding box computation.\nAttempt to write intermediate results had %d errors.\n", i);
        	exit(EXIT_FAILURE);
      }else{}
    }
    else {
      tmp_name=my_strcat(sub_name,".");              
      if(tmp_name == NULL){
		fprintf(world->err_file, "File '%s', Line %ld: Out of memory, trying to save intermediate results.\n", __FILE__, (long)__LINE__);
        	i = emergency_output();
		fprintf(world->err_file, "Fatal error: out of memory while bounding box computation.\nAttempt to write intermediate results had %d errors.\n", i);
        	exit(EXIT_FAILURE);
      }else{}
    }
    sub_name=my_strcat(tmp_name,objp->last_name);    
    if(sub_name == NULL){
		fprintf(world->err_file, "File '%s', Line %ld: Out of memory, trying to save intermediate results.\n", __FILE__, (long)__LINE__);
        	i = emergency_output();
		fprintf(world->err_file, "Fatal error: out of memory while bounding box computation.\nAttempt to write intermediate results had %d errors.\n", i);
        	exit(EXIT_FAILURE);
    }else{}
    free((void *)tmp_name);
  }
  else {
    sub_name=strdup(objp->last_name);    
    if(sub_name == NULL){
		fprintf(world->err_file, "File '%s', Line %ld: Out of memory, trying to save intermediate results.\n", __FILE__, (long)__LINE__);
        	i = emergency_output();
		fprintf(world->err_file, "Fatal error: out of memory while bounding box computation.\nAttempt to write intermediate results had %d errors.\n", i);
        	exit(EXIT_FAILURE);
    }else{}
  }

  switch (objp->object_type) {
  case META_OBJ:
    no_printf("Bounding box of Meta object %s is computed\n",sub_name);
    fflush(log_file);
    child_objp=objp->first_child;
    while (child_objp!=NULL) {
      if (compute_bb(child_objp,tm,sub_name)) {
        return(1);
      }
      child_objp=child_objp->next;
    }
    break;
  case REL_SITE_OBJ:
    no_printf("Bounding box of Release site object %s is computed\n",sub_name);
    fflush(log_file);
    if (compute_bb_release_site(objp,tm)) {
      return(1);
    }
    break;
  case BOX_OBJ:
    no_printf("Bounding box of Box object %s is computed\n",sub_name);
    fflush(log_file);
    if (compute_bb_polygon_object(objp,tm,sub_name)) {
      return(1);
    }
    break;
  case POLY_OBJ:
    no_printf("Bounding box of Polygon list object %s is computed\n",sub_name);
    fflush(log_file);
    if (compute_bb_polygon_object(objp,tm,sub_name)) {
      return(1);
    }
    break;
  }

  free((void *)sub_name);
  return(0);
}



/**
 * Updates the bounding box of the world based on the size
 * and location of a release site.
 * Used by compute_bb().
 */
int compute_bb_release_site(struct object *objp, double (*im)[4])
{
  struct release_site_obj *rsop;
  double location[1][4];
  unsigned short l,m,n;
  double diam_x, diam_y, diam_z; /* diameters of the release_site */ 
 
  rsop=(struct release_site_obj *)objp->contents;
  
  if (rsop->release_shape == SHAPE_REGION) return 0;

  l=1;
  m=4;
  n=4;

  if(rsop->location == NULL){
     fprintf(world->err_file, "ERROR: location is not specified for the geometrical shape release site.\n");
     return 1;
  }

  location[0][0]=rsop->location->x;
  location[0][1]=rsop->location->y;
  location[0][2]=rsop->location->z;
  location[0][3]=1.0;
  mult_matrix(location,im,location,l,m,n);
  
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

  return(0);
}



/**
 * Updates the bounding box of the world based on the size
 * and location of a polygon_object.
 * Used by compute_bb().
 */
int compute_bb_polygon_object(struct object *objp, double (*im)[4], char *full_name)
{
  struct polygon_object *pop;
  struct ordered_poly *opp;
  double p[1][4];
  int i,n_verts,n_walls;
  unsigned short l,m,n;

  pop=(struct polygon_object *)objp->contents;
  n_walls=pop->n_walls;
  l=1;
  m=4;
  n=4;

  opp=(struct ordered_poly *)pop->polygon_data;
  n_verts=opp->n_verts;

  for (i=0;i<n_verts;i++) {
    p[0][0]=opp->vertex[i].x;
    p[0][1]=opp->vertex[i].y;
    p[0][2]=opp->vertex[i].z;
    p[0][3]=1.0;
    mult_matrix(p,im,p,l,m,n);
    if (p[0][0]<world->bb_llf.x) world->bb_llf.x = p[0][0];
    if (p[0][1]<world->bb_llf.y) world->bb_llf.y = p[0][1];
    if (p[0][2]<world->bb_llf.z) world->bb_llf.z = p[0][2];
    if (p[0][0]>world->bb_urb.x) world->bb_urb.x = p[0][0];
    if (p[0][1]>world->bb_urb.y) world->bb_urb.y = p[0][1];
    if (p[0][2]>world->bb_urb.z) world->bb_urb.z = p[0][2];
  }

  return(0);
}


/**
 * Instantiates a polygon_object.
 * Creates walls from a template polygon_object or box object
 * as defined in the MDL file after applying the necessary geometric
 * transformations (scaling, rotation and translation).
 * <br>
 */
int instance_polygon_object(struct object *objp, double (*im)[4], struct viz_obj *vizp, char *full_name)
{
// #define INIT_VERTEX_NORMALS
// Uncomment to compute vertex normals
  FILE *log_file;
  struct polygon_object *pop;
  struct ordered_poly *opp;
  struct vector3 *v,**vp;
  struct wall *w,**wp;
  struct viz_child *vcp;
  double p[1][4];
#ifdef INIT_VERTEX_NORMALS
  struct vector3 *vertex_normal;
  double origin[1][4];
#endif
  double total_area;
  int i,n_verts,n_walls,index_0,index_1,index_2;
  unsigned int degenerate_count;
  unsigned short l,m,n;
  byte compute_vertex_normals;

  log_file=world->log_file;
  pop=(struct polygon_object *)objp->contents;
  n_walls=pop->n_walls;
  n_verts=pop->n_verts;
  l=1;
  m=4;
  n=4;
  total_area=0;

/* Allocate and initialize walls and vertices */
    w=(struct wall *)malloc(n_walls*sizeof(struct wall));
    wp=(struct wall **)malloc(n_walls*sizeof(struct wall *));
    v=(struct vector3 *)malloc(n_verts*sizeof(struct vector3));
    vp=(struct vector3 **)malloc(n_verts*sizeof(struct vector3 *)); 
    if (w==NULL || wp==NULL || v==NULL || vp==NULL)
    {
      fprintf(world->err_file,"File '%s', Line %ld: Out of memory while instantiating polygon object.  Quitting.\n", __FILE__, (long)__LINE__);
      exit(EXIT_FAILURE);
    }
    objp->walls=w;
    objp->wall_p=wp;
    objp->verts=v;
    objp->vert_p=vp;

    opp=(struct ordered_poly *)pop->polygon_data;

    compute_vertex_normals=0;

/* If we want vertex normals we'll have to add a place to store them
   in struct object.
*/
#ifdef INIT_VERTEX_NORMALS
    if (opp->normal!=NULL) {
      compute_vertex_normals=1;
    }
#endif
   if(vizp!=NULL)
   {
     if((world->viz_mode == DREAMM_V3_MODE) || (world->viz_mode == DREAMM_V3_GROUPED_MODE) || (objp->viz_state!=NULL))
     {

      if ((vcp=(struct viz_child *)malloc
           (sizeof(struct viz_child)))==NULL) {
		fprintf(world->err_file, "File '%s', Line %ld: Out of memory, trying to save intermediate results.\n", __FILE__, (long)__LINE__);
        	int i = emergency_output();
		fprintf(world->err_file, "Fatal error: out of memory while instantiation of polygon object.\nAttempt to write intermediate results had %d errors.\n", i);
        	exit(EXIT_FAILURE);
      }
      vcp->obj = objp;
      vcp->next = vizp->viz_child_head;
      vizp->viz_child_head = vcp;
    }
  }  

  for (i=0;i<n_verts;i++) {
    vp[i]=&v[i];
    p[0][0]=opp->vertex[i].x;
    p[0][1]=opp->vertex[i].y;
    p[0][2]=opp->vertex[i].z;
    p[0][3]=1.0;
    mult_matrix(p,im,p,l,m,n);
    v[i].x=p[0][0];
    v[i].y=p[0][1];
    v[i].z=p[0][2];

#ifdef INIT_VERTEX_NORMALS
    if (compute_vertex_normals) {
      p[0][0]=opp->normal[i].x;
      p[0][1]=opp->normal[i].y;
      p[0][2]=opp->normal[i].z;
      p[0][3]=1.0;
      origin[0][0]=0;
      origin[0][1]=0;
      origin[0][2]=0;
      origin[0][3]=1.0;
      mult_matrix(p,im,p,l,m,n);
      mult_matrix(origin,im,origin,l,m,n);
      vertex_normal[i].x=p[0][0]-origin[0][0];
      vertex_normal[i].y=p[0][1]-origin[0][1];
      vertex_normal[i].z=p[0][2]-origin[0][2];
      normalize(&vertex_normal[i]);
    }
#endif
  }
  
  degenerate_count=0;
  for (i=0;i<n_walls;i++) {
    if (!get_bit(pop->side_removed,i)) {
      wp[i]=&w[i];
      index_0=opp->element[i].vertex_index[0];
      index_1=opp->element[i].vertex_index[1];
      index_2=opp->element[i].vertex_index[2];

      init_tri_wall(objp,i,vp[index_0],vp[index_1],vp[index_2]);
      total_area+=wp[i]->area;

      if (wp[i]->area==0) {
        if (world->notify->degenerate_polys != WARN_COPE)
        {
          if (world->notify->degenerate_polys==WARN_ERROR)
          {
            log_file = world->err_file;
            fprintf(log_file,"\nError -- ");
          }
          else fprintf(log_file,"\nWarning -- ");
          
          fprintf(log_file,"Degenerate polygon found and automatically removed: %s %d\n\n",objp->sym->name,i);
          fprintf(log_file,"  Vertex 0: %.5e %.5e %.5e\n",vp[index_0]->x,vp[index_0]->y,vp[index_0]->z);
          fprintf(log_file,"  Vertex 1: %.5e %.5e %.5e\n",vp[index_1]->x,vp[index_1]->y,vp[index_1]->z);
          fprintf(log_file,"  Vertex 2: %.5e %.5e %.5e\n",vp[index_2]->x,vp[index_2]->y,vp[index_2]->z);
          
          if (world->notify->degenerate_polys==WARN_ERROR) return 1;
        }
        set_bit(pop->side_removed,i,1);
        objp->n_walls_actual--;
        degenerate_count++;
        wp[i]=NULL;
      }
    }
    else {
      wp[i]=NULL;
    }
  }
  if (degenerate_count) mdl_remove_gaps_from_regions(objp);
  
  objp->total_area=total_area;
  
#ifdef DEBUG    
  printf("n_walls = %d\n", n_walls);
  printf("n_walls_actual = %d\n", objp->n_walls_actual);
#endif

  return(0);
}


/********************************************************************
 init_regions:

    Traverse the world initializing regions on each object.

    In:  none
    Out: 0 on success, 1 on failure
 *******************************************************************/
int init_regions()
{
  if (world->clamp_list!=NULL) init_clamp_lists();

  if (instance_obj_regions(world->root_instance,NULL)) {
    return(1);
  }

  return(0);
}


/* First part of concentration clamp initialization. */
/* After this, list is grouped by surface class. */
/* Second part (list of objects) happens with regions. */
void init_clamp_lists()
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
 * them.  "sub_name" contains the name of the parent object of the one we're
 * instancing as we traverse down into the object hierarchy 
 */
int instance_obj_regions(struct object *objp,char *sub_name)
{
  FILE *log_file;
  struct object *child_objp;
  char *tmp_name;

  log_file=world->log_file;

  if (sub_name!=NULL) { 
    if (strcmp(sub_name,"")==0) {
      tmp_name=strdup("");
      if (tmp_name == NULL) {
		fprintf(world->err_file, "File %s, Line %ld: Out of memory, trying to save intermediate results.\n", __FILE__, (long)__LINE__);
        	int i = emergency_output();
		fprintf(world->err_file, "Fatal error: out of memory while instantiation of object regions.\nAttempt to write intermediate results had %d errors.\n", i);
        	exit(EXIT_FAILURE);
      }else{}
    }
    else {
      tmp_name=my_strcat(sub_name,".");              
      if (tmp_name == NULL) {
		fprintf(world->err_file, "File %s, Line %ld: Out of memory, trying to save intermediate results.\n", __FILE__, (long)__LINE__);
        	int i = emergency_output();
		fprintf(world->err_file, "Fatal error: out of memory while instantiation of object regions.\nAttempt to write intermediate results had %d errors.\n", i);
        	exit(EXIT_FAILURE);
      }else{}
    }
    sub_name=my_strcat(tmp_name,objp->last_name);    
    if (sub_name == NULL) {
		fprintf(world->err_file, "File %s, Line %ld: Out of memory, trying to save intermediate results.\n", __FILE__, (long)__LINE__);
        	int i = emergency_output();
		fprintf(world->err_file, "Fatal error: out of memory while instantiation of object regions.\nAttempt to write intermediate results had %d errors.\n", i);
        	exit(EXIT_FAILURE);
    }else{}
    free((void *)tmp_name);
  }
  else {
    sub_name=strdup(objp->last_name);    
    if (sub_name == NULL) {
		fprintf(world->err_file, "File %s, Line %ld: Out of memory, trying to save intermediate results.\n", __FILE__, (long)__LINE__);
        	int i = emergency_output();
		fprintf(world->err_file, "Fatal error: out of memory while instantiation of object regions.\nAttempt to write intermediate results had %d errors.\n", i);
        	exit(EXIT_FAILURE);
    }else{}
  }

  switch (objp->object_type) {
  case META_OBJ:
    no_printf("Initializing regions in meta object: %s\n",sub_name);
    fflush(log_file);
    child_objp=objp->first_child;
    while (child_objp!=NULL) {
      if (instance_obj_regions(child_objp,sub_name)) {
        return(1);
      }
      child_objp=child_objp->next;
    }
    break;
  case REL_SITE_OBJ:
    break;
  case BOX_OBJ:
    no_printf("Initializing regions in box object: %s\n",sub_name);
    fflush(log_file);
    if (init_wall_regions(objp,sub_name)) {
      return(1);
    }
    break;
  case POLY_OBJ:
    no_printf("Initializing regions in polygon list object: %s\n",sub_name);
    fflush(log_file);
    if (init_wall_regions(objp,sub_name)) {
      return(1);
    }
    break;
  }

  free((void *)sub_name);
  return(0);
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
int init_wall_regions(struct object *objp, char *full_name)
{
  FILE *log_file;
  struct polygon_object *pop;
  struct wall *w;
  struct region *rp;
  struct region_list *rlp,*wrlp;
  int i,n_walls;

  log_file=world->log_file;

  pop=(struct polygon_object *)objp->contents;
  n_walls=pop->n_walls;

  no_printf("Processing %d regions in polygon list object: %s\n",objp->num_regions,full_name);  

  /* prepend a copy of eff_dat for each element referenced in each region
     of this object to the eff_prop list for the referenced element */
  rlp=objp->regions;
  
  for (rlp=objp->regions ; rlp!=NULL ; rlp=rlp->next)
  {
    rp=rlp->reg;
    if (rp->membership==NULL)
    {
      fprintf(world->err_file,"File '%s', Line %ld: Internal error. incomplete region information for %s\n", __FILE__, (long)__LINE__, rp->sym->name);
      return 1;
    }

    for (i=0;i<rp->membership->nbits;i++)
    {
      if (get_bit(rp->membership,i))
      {
	/* prepend this region to wall region list of i_th wall only if the region is used in counting */
	w=objp->wall_p[i];
	rp->area += w->area;
	if (rp->surf_class!=NULL) w->surf_class = rp->surf_class;

	if ((rp->flags & COUNT_SOME_MASK) != 0)
	{  
	  wrlp = (struct region_list *)mem_get(w->birthplace->regl);
	  if (wrlp==NULL)
	  {
	    fprintf(world->err_file,"File '%s', Line %ld: Out of memory, can't place regions on geometry.\n", __FILE__, (long)__LINE__);
	    return 1;
	  }
	  wrlp->reg=rp;
	  wrlp->next=w->counting_regions;
	  w->counting_regions=wrlp;
	  w->flags|=rp->flags;
	}
      }
    }

  } /*end loop over all regions in object */
  for (i=0;i<n_walls;i++)
  {
    if (get_bit(pop->side_removed,i)) continue; 
     
    w=objp->wall_p[i];
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
    
    for (i=0;i<n_walls;i++)
    {
      if (get_bit(pop->side_removed,i)) continue;
      if (objp->wall_p[i]->surf_class != world->g_surf)
      {
        for (ccd=world->clamp_list ; ccd!=NULL ; ccd=ccd->next)
        {
          if (objp->wall_p[i]->surf_class == ccd->surf_class)
          {
            if (ccd->objp!=objp)
            {
              if (ccd->objp==NULL) ccd->objp=objp;
              else if (ccd->next_obj != NULL && ccd->next_obj->objp==objp) ccd=ccd->next_obj;
              else
              {
                temp = (struct ccn_clamp_data*)malloc(sizeof(struct ccn_clamp_data));
                if (temp==NULL)
                {
                  fprintf(world->err_file,"File '%s', Line %ld: Out of memory assembling concentration clamp data.\n", __FILE__, (long)__LINE__);
                  return 1;
                }
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
              {
                fprintf(world->err_file,"File '%s', Line %ld: Out of memory assembling concentration clamp data.\n", __FILE__, (long)__LINE__);
                return 1;
              }
              set_all_bits(ccd->sides,0);
            }
            set_bit(ccd->sides,i,1);
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
        
        ccd->side_idx = (int*)malloc(ccd->n_sides*sizeof(int));
        ccd->cum_area = (double*)malloc(ccd->n_sides*sizeof(double));
        if (ccd->side_idx==NULL || ccd->cum_area==NULL)
        {
          fprintf(world->err_file,"File '%s', Line %ld: Out of memory assembling concentration clamp data.\n", __FILE__, (long)__LINE__);
          return 1;
        }
        
        j=0;
        for (i=0;i<n_walls;i++)
        {
          if (get_bit(ccd->sides,i))
          {
            ccd->side_idx[j] = i;
            ccd->cum_area[j] = objp->wall_p[i]->area;
            j++;
          }
        }
        if (j!=ccd->n_sides)
        {
          fprintf(world->err_file,"File '%s', Line %ld: Miscounted the number of walls for concentration clamp\n  on object %s\n  surface class %s\n", __FILE__, (long)__LINE__, objp->sym->name,ccd->surf_class->sym->name);
          return 1;
        }
        
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

  return(0);
}



/********************************************************************
 init_effectors:

    Traverse the world placing grid molecules.

    In:  none
    Out: 0 on success, 1 on failure
 *******************************************************************/
int init_effectors()
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

  return(0);
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
  FILE *log_file;
  struct polygon_object *pop;
  struct wall *w;
  struct eff_dat *effdp,*dup_effdp,**eff_prop;
  struct region *rp;
  struct region_list *rlp,*rlp2,*reg_eff_num_head,*complex_head;
  int i,n_walls;
  byte reg_eff_num;
  byte complex_eff;

  log_file=world->log_file;

  pop=(struct polygon_object *)objp->contents;
  n_walls=pop->n_walls;

   
  /* allocate scratch storage to hold effector info for each wall */
  if ((eff_prop=(struct eff_dat **)malloc(n_walls*sizeof(struct eff_dat *)))==NULL)
  {
    fprintf(world->err_file,"File '%s', Line %ld:  Out of memory, can't create space for molecules on a region.\n", __FILE__, (long)__LINE__);
    return 1;
  }

  for (i=0;i<n_walls;i++) eff_prop[i]=NULL; 

  /* prepend a copy of eff_dat for each element referenced in each region
     of this object to the eff_prop list for the referenced element */
  reg_eff_num_head=NULL;

  /* List of regions which need macromol processing */
  complex_head=NULL;
  
  rlp=objp->regions;
  for (rlp=objp->regions ; rlp!=NULL ; rlp=rlp->next)
  {
    rp=rlp->reg;
    reg_eff_num=0;
    complex_eff=0;

    for (i=0;i<rp->membership->nbits;i++)
    {
      if (get_bit(rp->membership,i))
      {
	w=objp->wall_p[i];

	/* prepend region eff data for this region to eff_prop for i_th wall */
        for ( effdp=rp->eff_dat_head ; effdp!=NULL ; effdp=effdp->next )
	{
          if (effdp->eff->flags & IS_COMPLEX)
            complex_eff = 1;
          else if (effdp->quantity_type==EFFDENS)
	  {
	    if ((dup_effdp=(struct eff_dat *)malloc(sizeof(struct eff_dat)))==NULL)
	    {
	      fprintf(world->err_file,"File '%s', Line %ld: Out of memory, can't create space for molecules on a region.\n", __FILE__, (long)__LINE__);
	      return 1;
	    }
	    dup_effdp->eff=effdp->eff;
	    dup_effdp->quantity_type=effdp->quantity_type;
	    dup_effdp->quantity=effdp->quantity;
	    dup_effdp->orientation=effdp->orientation;
	    dup_effdp->next=eff_prop[i];
	    eff_prop[i]=dup_effdp;
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
	      if ((dup_effdp=(struct eff_dat *)malloc(sizeof(struct eff_dat)))==NULL)
	      {
		fprintf(world->err_file,"File '%s', Line %ld: Out of memory, can't create space for molecules on a region.\n", __FILE__, (long)__LINE__);
		return 1;
	      }
	      dup_effdp->eff=effdp->eff;
	      dup_effdp->quantity_type=effdp->quantity_type;
	      dup_effdp->quantity=effdp->quantity;
	      dup_effdp->orientation=effdp->orientation;
	      dup_effdp->next=eff_prop[i];
	      eff_prop[i]=dup_effdp;
	    }
	    else reg_eff_num=1;
	  }
	}

      }
    } /* done checking each wall */
    
    if (complex_eff)
    {
      if ((rlp2=(struct region_list *)malloc(sizeof(struct region_list)))==NULL)
      {
	fprintf(world->err_file,"File '%s', Line %ld: Out of memory, can't place regions on geometry.\n", __FILE__, (long)__LINE__);
	return 1;
      }
      rlp2->reg=rp;
      rlp2->next=complex_head;
      complex_head=rlp2;
    }
    else if (reg_eff_num)
    {
      if ((rlp2=(struct region_list *)malloc(sizeof(struct region_list)))==NULL)
      {
	fprintf(world->err_file,"File '%s', Line %ld: Out of memory, can't place regions on geometry.\n", __FILE__, (long)__LINE__);
	return 1;
      }
      rlp2->reg=rp;
      rlp2->next=reg_eff_num_head;
      reg_eff_num_head=rlp2;
    }
  } /*end for (... ; rlp != NULL ; ...) */

  /* Place macromolecular complexes, if any */
  if (complex_head!=NULL)
  {
    if (init_complex_effectors(objp, complex_head))
    {
      return 1;
    }
    
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
  for (i=0;i<n_walls;i++)
  {
    if (!get_bit(pop->side_removed,i))
    {
      if (eff_prop[i]!=NULL)
      {
	if (init_effectors_by_density(objp->wall_p[i],eff_prop[i])) return 1;
      }
    }
  }

  /* Place regular (non-macro) molecules by number */
  if (reg_eff_num_head!=NULL)
  {
    if (init_effectors_by_number(objp,reg_eff_num_head)) {
      return(1);
    }
    
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
  for (i=0;i<n_walls;i++)
  {
    if (eff_prop[i]!=NULL)
    {
      effdp=eff_prop[i];
      while(effdp!=NULL)
      {
	dup_effdp=effdp;
	effdp=effdp->next;
	free(dup_effdp);
      }
    }
  }
  free(eff_prop);

    
  return(0);
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
  int grid_idx;
  double p;
  short orient = effdp->orientation;

  /* Pick orientation */
  if (orient == 0)
  {
    orient = (rng_uint(world->rng) & 1) ? 1 : -1;
    if (world->notify->final_summary == NOTIFY_FULL)
    {
      world->random_number_use++;
    }
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
  /* How many times do we try each placement before we give up? */
  static const int MAX_TRIES = 25;

  if (nwalls <= 0)
    return 1;

  double max_weight = weights[nwalls - 1];
  int n_failures = 0;
  int n_total = n_to_place;
  while (n_to_place > 0)
  {
    int num_tries = MAX_TRIES;
    int chosen_wall = 0;
    double p = rng_dbl(world->rng) * max_weight;

    /* Pick a wall */
    if (world->notify->final_summary == NOTIFY_FULL)
      world->random_number_use++;
    chosen_wall = bisect_high(weights, nwalls, p);

    /* Try to find a spot for the release */
    while (-- num_tries >= 0)
    {
      if (! init_effectors_place_complex(walls[chosen_wall], rp, effdp))
        break;
    }

    if (num_tries >= 0)
      -- n_to_place;
    else
    {
      /* XXX: What criteria? */
      if (++ n_failures > 100)
      {
        fprintf(world->log_file,"\nMCell: Warning -- Unable to place some surface complexes of species '%s' (placed %d of %d).\n\n", effdp->eff->sym->name, n_total - n_to_place, n_total);
        fflush(world->log_file);
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
    int i;
    for (i=0; i<rp->membership->nbits; ++i)
    {
      if (get_bit(rp->membership, i))
      {
        w = objp->wall_p[i];
        if (w->grid == NULL  &&  create_grid(w, NULL))
          return 1;

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
        n_to_place = (int) (effdp->quantity * total_area + 0.5);
      else
      {
        fprintf(world->err_file, "File '%s', Line %ld: Unknown effector quantity type, trying to save intermediate results.\n", __FILE__, (long)__LINE__);
        int i = emergency_output();
        fprintf(world->err_file, "Fatal error: Unknown effector quantity type.\nAttempt to write intermediate results had %d errors.\n", i);
        exit(EXIT_FAILURE);
      }

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
  FILE *log_file;
  struct object *objp;
  struct species **eff;
  struct surface_grid *sg;
  struct eff_dat *effdp;
  struct grid_molecule *mol;
  short *orientation;
  unsigned int i,j,n,nr,n_occupied;
  int p_index;
  double rand,*prob,area,tot_prob,tot_density;
  struct subvolume *gsv = NULL;

  log_file=world->log_file;

  no_printf("Initializing effectors by density...\n");
  fflush(log_file);

  if (create_grid(w,NULL)) {
    return(1);
  }
  sg=w->grid;
  objp=w->parent_object;

  nr=0;
  effdp=effdp_head;
  while (effdp!=NULL) {
    nr++;
    effdp=effdp->next;
  }

  if ((eff=(struct species **)malloc(nr*sizeof(struct species *)))==NULL) {
			fprintf(world->err_file, "File '%s', Line %ld: Out of memory, trying to save intermediate results.\n", __FILE__, (long)__LINE__);
        		int i = emergency_output();
			fprintf(world->err_file, "Fatal error: out of memory while effectors by density intialization.\nAttempt to write intermediate results had %d errors.\n", i);
        		exit(EXIT_FAILURE);
  }
  if ((prob=(double *)malloc(nr*sizeof(double)))==NULL) {
			fprintf(world->err_file, "File '%s', Line %ld: Out of memory, trying to save intermediate results.\n", __FILE__, (long)__LINE__);
        		int i = emergency_output();
			fprintf(world->err_file, "Fatal error: out of memory while effectors by density intialization.\nAttempt to write intermediate results had %d errors.\n", i);
        		exit(EXIT_FAILURE);
  }
  if ((orientation=(short*)malloc(nr*sizeof(short)))==NULL) {
			fprintf(world->err_file, "File '%s', Line %ld: Out of memory, trying to save intermediate results.\n", __FILE__, (long)__LINE__);
        		int i = emergency_output();
			fprintf(world->err_file, "Fatal error: out of memory while effectors by density intialization.\nAttempt to write intermediate results had %d errors.\n", i);
        		exit(EXIT_FAILURE);
  }

  for (i=0;i<nr;i++) {
    eff[i]=NULL;
    prob[i]=0.0;
    orientation[i]=0;
  }

  n=sg->n_tiles;
  area=w->area;
  objp->n_tiles+=n;
  no_printf("Initializing %d effectors...\n",n);
  no_printf("  Area = %.9g\n",area);
  no_printf("  Grid_size = %d\n",sg->n);
  no_printf("  Number of effector types in wall = %d\n",nr);

  i=0;
  tot_prob=0;
  tot_density=0;
  effdp=effdp_head;
  while (effdp!=NULL) {
    no_printf("  Adding effector %s to wall at density %.9g\n",effdp->eff->sym->name,effdp->quantity);
    tot_prob+=(area*effdp->quantity)/(n*world->grid_density);
    prob[i]=tot_prob;
    if (effdp->orientation > 0) orientation[i] = 1;
    else if (effdp->orientation < 0) orientation[i] = -1;
    else{ 
       orientation[i] = (rng_uint(world->rng)&1)?1:-1;
       if(world->notify->final_summary == NOTIFY_FULL){
           world->random_number_use++;
       }
    }
    eff[i++]=effdp->eff;
    tot_density+=effdp->quantity;
    effdp=effdp->next;
  }

  if (tot_density>world->grid_density) {
    fprintf(log_file,"\nMCell: Warning -- Total effector density too high: %f.  Filling all available effector sites.\n\n",tot_density);
    fflush(log_file);
/*
    return(1);
*/
  }

  n_occupied=0;
  if (world->chkpt_init) {
    for (i=0;i<n;i++) {
      if (sg->mol[i] != NULL)
        continue;

      j=0;
      p_index=-1;
      rand = rng_dbl(world->rng);
      if(world->notify->final_summary == NOTIFY_FULL){
         world->random_number_use++;
      }
      while (j<nr && p_index==-1) {
        if (rand<=prob[j++]) {
          p_index=j-1;
        }
      }
      if (p_index!=-1) {
        struct vector2 s_pos;
        struct vector3 pos3d;
        n_occupied++;
        eff[p_index]->population++;

	if (world->randomize_gmol_pos) grid2uv_random(sg,i,&s_pos);
	else grid2uv(sg,i,&s_pos);
        uv2xyz(&s_pos, w, &pos3d);
        gsv = find_subvolume(&pos3d, gsv);
        mol=(struct grid_molecule *)mem_get(gsv->local_storage->gmol);
        if(mol == NULL){
			fprintf(world->err_file, "File '%s', Line %ld: Out of memory, trying to save intermediate results.\n", __FILE__, (long)__LINE__);
        		int i = emergency_output();
			fprintf(world->err_file, "Fatal error: out of memory while effectors by density intialization.\nAttempt to write intermediate results had %d errors.\n", i);
        		exit(EXIT_FAILURE);
        }
        sg->mol[i]=mol;
        mol->t=0;
        mol->t2=0;
	mol->birthday = 0;
        mol->cmplx = NULL;
        mol->s_pos.u = s_pos.u;
        mol->s_pos.v = s_pos.v;
        mol->properties=eff[p_index];
        mol->birthplace=w->birthplace->gmol;
        mol->grid_index=i;
        mol->orient=orientation[p_index];
        mol->grid=sg;

        mol->flags=TYPE_GRID|ACT_NEWBIE|IN_SCHEDULE|IN_SURFACE;
	if (mol->properties->space_step>0) mol->flags |= ACT_DIFFUSE;
        if ( trigger_unimolecular(eff[p_index]->hashval,(struct abstract_molecule *)mol)!=NULL
	     || (eff[p_index]->flags&CAN_GRIDWALL)!=0 ) {
          mol->flags|=ACT_REACT;
        }
        if ((mol->properties->flags&COUNT_ENCLOSED) != 0) mol->flags |= COUNT_ME;

        if ((mol->properties->flags & (COUNT_CONTENTS|COUNT_ENCLOSED)) != 0)
          count_region_from_scratch((struct abstract_molecule*)mol,NULL,1,NULL,NULL,mol->t);
      
        if (schedule_add(gsv->local_storage->timer,mol)){ 
		fprintf(world->err_file, "File '%s', Line %ld: Out of memory, trying to save intermediate results.\n", __FILE__, (long)__LINE__);
        	int i = emergency_output();
		fprintf(world->err_file, "Fatal error: out of memory while effectors by density intialization.\nAttempt to write intermediate results had %d errors.\n", i);
        	exit(EXIT_FAILURE);
         }
      }
    }
  }

  sg->n_occupied=n_occupied;
  objp->n_occupied_tiles+=n_occupied;

  for (i=0;i<nr;i++) {
    no_printf("Total number of effector %s = %d\n",eff[i]->sym->name,eff[i]->population);
  }

  free(eff);
  free(prob);
  free(orientation);

  no_printf("Done initializing %d effectors by density\n",n_occupied);
  fflush(log_file);

	
  return(0);
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
  FILE *log_file;
  struct polygon_object *pop;
  struct species *eff;
  struct grid_molecule ***tiles,***tiles_tmp;
  struct grid_molecule gmol,*bread_crumb,*mol;
  struct region_list *rlp; 
  struct region *rp;
/*  struct element_list *elp; */
  struct surface_grid *sg;
  struct eff_dat *effdp;
  struct wall **walls,**walls_tmp,*w;
  short orientation;
  unsigned int *index,*index_tmp;
  unsigned int n_free_eff,n_set,n_clear;
  unsigned int i,j,k;
  byte done;
  struct subvolume *gsv=NULL;
  struct vector3 pos3d;

    no_printf("Initializing effectors by number...\n");

    log_file=world->log_file;
    pop=(struct polygon_object *)objp->contents;

    tiles=NULL;
    tiles_tmp=NULL;
    index=NULL;
    index_tmp=NULL;
    walls=NULL;
    walls_tmp=NULL;
    bread_crumb=&gmol;

    /* traverse region list and add effector sites by number to whole regions
       as appropriate */
    rlp=reg_eff_num_head;
    while (rlp!=NULL) {
      rp=rlp->reg;
        /* initialize effector grids in region as needed and */
        /* count total number of free effector sites in region */
        n_free_eff=0;
	for (i=0;i<rp->membership->nbits;i++)
	{
	  if (get_bit(rp->membership,i))
	  {
	    w=objp->wall_p[i];
	    if (create_grid(w,NULL)) {
	      return(1);
	    }
	    sg=w->grid;
	    n_free_eff=n_free_eff+(sg->n_tiles-sg->n_occupied);
          }
        }
        no_printf("Number of free effector tiles in region %s = %d\n",rp->sym->name,n_free_eff);  
        fflush(stdout);

        if(n_free_eff == 0) {
              fprintf(log_file,"\nMCell: Warning -- Number of free effector tiles in region %s = %d\n", rp->sym->name, n_free_eff);
              rlp = rlp->next;   
              continue;
        }
 
      if (world->chkpt_init) {  /* only needed for denovo initiliazation */
        /* allocate memory to hold array of pointers to all free tiles */
        if ((tiles=(struct grid_molecule ***)malloc
           (n_free_eff*sizeof(struct grid_molecule **)))==NULL) {
			fprintf(world->err_file, "File '%s', Line %ld: Out of memory, trying to save intermediate results.\n", __FILE__, (long)__LINE__);
        		int i = emergency_output();
			fprintf(world->err_file, "Fatal error: out of memory while effectors by number intialization.\nAttempt to write intermediate results had %d errors.\n", i);
        		exit(EXIT_FAILURE);
        }
        if ((index=(unsigned int *)malloc
           (n_free_eff*sizeof(unsigned int)))==NULL) {
			fprintf(world->err_file, "File '%s', Line %ld: Out of memory, trying to save intermediate results.\n", __FILE__, (long)__LINE__);
        		int i = emergency_output();
			fprintf(world->err_file, "Fatal error: out of memory while effectors by number intialization.\nAttempt to write intermediate results had %d errors.\n", i);
        		exit(EXIT_FAILURE);
        }
        if ((walls=(struct wall **)malloc
           (n_free_eff*sizeof(struct wall *)))==NULL) {
			fprintf(world->err_file, "File '%s', Line %ld: Out of memory, trying to save intermediate results.\n", __FILE__, (long)__LINE__);
        		int i = emergency_output();
			fprintf(world->err_file, "Fatal error: out of memory while effectors by number intialization.\nAttempt to write intermediate results had %d errors.\n", i);
        		exit(EXIT_FAILURE);
        }
        /* initialize array of pointers to all free tiles */
        k=0;
	for (i=0;i<rp->membership->nbits;i++)
	{
	  if (get_bit(rp->membership,i))
	  {
	    w=objp->wall_p[i];
	    sg=w->grid;
	    if (sg!=NULL) {
	      for (j=0;j<sg->n_tiles;j++) {
		if (sg->mol[j]==NULL) {
		  tiles[k]=&(sg->mol[j]);
		  index[k]=j;
		  walls[k++]=w;
		}
	      }
	    }
	  }
        }
      } /* end while(world->chkpt_init) */


      /* distribute desired number of effector sites */
      /* for each effector type to add */
      effdp=rp->eff_dat_head;
      while (effdp!=NULL) {
        if (effdp->quantity_type==EFFNUM) {
          eff=effdp->eff;

          if (world->chkpt_init) {  /* only needed for denovo initiliazation */
	    if (effdp->orientation > 0) orientation = 1;
	    else if (effdp->orientation < 0) orientation = -1;
	    else orientation = 0;
  
            n_set=effdp->quantity;
            n_clear=n_free_eff-n_set;

            if (n_set > n_free_eff) {
              fprintf(log_file,"\nMCell: Warning -- Number of %s effectors to place (%d) exceeds number of free effector tiles (%d) in region %s[%s].",eff->sym->name,n_set,n_free_eff,rp->parent->sym->name,rp->region_last_name);
               fprintf(log_file, "  Effectors %s placed on all available effector sites.\n\n", eff->sym->name);
              n_set=n_free_eff;
              n_clear=0;
            }
            
            eff->population+=n_set;
            
            no_printf("distribute %d of effector %s\n",n_set,eff->sym->name);
            no_printf("n_set = %d  n_clear = %d  n_free_eff = %d\n",n_set,n_clear,n_free_eff);
            fflush(stdout);

            /* if filling more than half the free tiles
              init all with bread_crumbs
              choose which tiles to free again
              and then convert remaining bread_crumbs to actual molecules */
            if (n_set > n_free_eff/2) {
              no_printf("filling more than half the free tiles: init all with bread_crumb\n");
              fflush(stdout);
              for (j=0;j<n_free_eff;j++) {
                *tiles[j]=bread_crumb;
              }

              no_printf("choose which tiles to free again\n");
              fflush(stdout);
              for (j=0;j<n_clear;j++) {
                done=0;
                while (!done) {
		  k = (int) (rng_dbl(world->rng)*n_free_eff);
                  if(world->notify->final_summary == NOTIFY_FULL){
                      world->random_number_use++;
                  }
                  if (*tiles[k]==bread_crumb) {
                    *tiles[k]=NULL;
                    done=1;
                  }
                }
              }

              no_printf("convert remaining bread_crumbs to actual molecules\n");
              fflush(stdout);
              for (j=0;j<n_free_eff;j++) {
                if (*tiles[j]==bread_crumb) {
                  struct vector2 s_pos;
		  if (world->randomize_gmol_pos) grid2uv_random(walls[j]->grid,index[j],&s_pos);
		  else grid2uv(walls[j]->grid,index[j],&s_pos);
                  uv2xyz(&s_pos, walls[j], &pos3d);
                  gsv = find_subvolume(&pos3d, gsv);

                  mol=(struct grid_molecule *)
                    mem_get(gsv->local_storage->gmol);
                  if (mol == NULL){
			fprintf(world->err_file, "File '%s', Line %ld: Out of memory, trying to save intermediate results.\n", __FILE__, (long)__LINE__);
        		int i = emergency_output();
			fprintf(world->err_file, "Fatal error: out of memory while effectors by number initialization.\nAttempt to write intermediate results had %d errors.\n", i);
        		exit(EXIT_FAILURE);
                  }
                  *tiles[j]=mol;
                  mol->t=0;
                  mol->t2=0;
		  mol->birthday=0;
                  mol->properties=eff;
                  mol->birthplace=walls[j]->birthplace->gmol;
                  mol->grid_index=index[j];
                  mol->s_pos.u = s_pos.u;
                  mol->s_pos.v = s_pos.v;
                  /*mol->orient = (orientation==0) ? ((rng_uint(world->rng)&1)?1:-1) : orientation; */
                  if(orientation == 0){
                     if((rng_uint(world->rng)&1)){
                        mol->orient = 1;
                     }else{
                        mol->orient = -1;
                     }
                     if(world->notify->final_summary == NOTIFY_FULL){
                        world->random_number_use++;
                     }
                  }else{
                     mol->orient = orientation;
                  }
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
      
                  if ( schedule_add(gsv->local_storage->timer, mol) ){ 
			fprintf(world->err_file, "File '%s', Line %ld: Out of memory, trying to save intermediate results.\n", __FILE__, (long)__LINE__);
        		int i = emergency_output();
			fprintf(world->err_file, "Fatal error: out of memory while effectors by number initialization.\nAttempt to write intermediate results had %d errors.\n", i);
        		exit(EXIT_FAILURE);
		   }
                }
              }
            }
            else {  /* just fill only the tiles we need */
              no_printf("fill only the tiles we need\n");
              fflush(stdout);
              for (j=0;j<n_set;j++) {
                done=0;
                while (!done) {
		  k = (int) (rng_dbl(world->rng)*n_free_eff);
                  if(world->notify->final_summary == NOTIFY_FULL){
                      world->random_number_use++;
                  }
                  if (*tiles[k]==NULL) {
                    struct vector2 s_pos;
                    if (world->randomize_gmol_pos) grid2uv_random(walls[k]->grid,index[k],&s_pos);
                    else grid2uv(walls[k]->grid,index[k],&s_pos);
                    uv2xyz(&s_pos, walls[k], &pos3d);
                    gsv = find_subvolume(&pos3d, gsv);

                    mol=(struct grid_molecule *)mem_get
                      (gsv->local_storage->gmol);
                    if (mol == NULL){
			fprintf(world->err_file, "File '%s', Line %ld: Out of memory, trying to save intermediate results.\n", __FILE__, (long)__LINE__);
        		int i = emergency_output();
			fprintf(world->err_file, "Fatal error: out of memory while effectors by number initialization.\nAttempt to write intermediate results had %d errors.\n", i);
        		exit(EXIT_FAILURE);
                    }
                    *tiles[k]=mol;
                    mol->t=0;
                    mol->t2=0;
		    mol->birthday=0;                                                                                    
                    mol->properties=eff;
                    mol->birthplace=walls[k]->birthplace->gmol;
                    mol->grid_index=index[k];
                    mol->s_pos.u = s_pos.u;
                    mol->s_pos.v = s_pos.v;
                    mol->cmplx = NULL;
                    /* mol->orient = (orientation==0) ? ((rng_uint(world->rng)&1)?1:-1) : orientation; */
                    if(orientation == 0){
                       if((rng_uint(world->rng)&1)){
                          mol->orient = 1;
                       }else{
                          mol->orient = -1;
                       }
                       if(world->notify->final_summary == NOTIFY_FULL){
                          world->random_number_use++;
                       }
                    }else{
                       mol->orient = orientation;
                    }

                    mol->grid=walls[k]->grid;
                    mol->flags=TYPE_GRID|ACT_NEWBIE|IN_SCHEDULE|IN_SURFACE;
		    if (mol->properties->space_step > 0) mol->flags |= ACT_DIFFUSE;
                      if (trigger_unimolecular(eff->hashval,(struct abstract_molecule *)mol)!=NULL
		          || (eff->flags&CAN_GRIDWALL)!=0) {
                      mol->flags|=ACT_REACT;
                    }
                  
                    if ((mol->properties->flags & (COUNT_CONTENTS|COUNT_ENCLOSED)) != 0)
                      count_region_from_scratch((struct abstract_molecule*)mol,NULL,1,NULL,NULL,mol->t);
      
                    if ( schedule_add(gsv->local_storage->timer,mol) ){ 
			fprintf(world->err_file, "File '%s', Line %ld: Out of memory, trying to save intermediate results.\n", __FILE__, (long)__LINE__);
        		int i = emergency_output();
			fprintf(world->err_file, "Fatal error: out of memory while effectors by number initialization.\nAttempt to write intermediate results had %d errors.\n", i);
        		exit(EXIT_FAILURE);
                     }
                    done=1;
                  }
                }
              }
            }
         
      if(n_clear > 0)
      { 
        /* allocate memory to hold array of pointers to remaining free tiles */
            if ((tiles_tmp=(struct grid_molecule ***)malloc
                 (n_clear*sizeof(struct grid_molecule **)))==NULL) {
			fprintf(world->err_file, "File '%s', Line %ld: Out of memory, trying to save intermediate results.\n", __FILE__, (long)__LINE__);
        		int i = emergency_output();
			fprintf(world->err_file, "Fatal error: out of memory while effectors by number initialization.\nAttempt to write intermediate results had %d errors.\n", i);
        		exit(EXIT_FAILURE);
            }
            if ((index_tmp=(unsigned int *)malloc
               (n_clear*sizeof(unsigned int)))==NULL) {
			fprintf(world->err_file, "File '%s', Line %ld: Out of memory, trying to save intermediate results.\n", __FILE__, (long)__LINE__);
        		int i = emergency_output();
			fprintf(world->err_file, "Fatal error: out of memory while effectors by number initialization.\nAttempt to write intermediate results had %d errors.\n", i);
        		exit(EXIT_FAILURE);
            }
            if ((walls_tmp=(struct wall **)malloc
               (n_clear*sizeof(struct wall *)))==NULL) {
			fprintf(world->err_file, "File '%s', Line %ld: Out of memory, trying to save intermediate results.\n", __FILE__, (long)__LINE__);
        		int i = emergency_output();
			fprintf(world->err_file, "Fatal error: out of memory while effectors by number initialization.\nAttempt to write intermediate results had %d errors.\n", i);
        		exit(EXIT_FAILURE);
            }
            k=0;
            for (i=0;i<n_free_eff;i++) {
              if (*tiles[i]==NULL) {
                tiles_tmp[k]=tiles[i];
                index_tmp[k]=index[i];
                walls_tmp[k++]=walls[i];
              }
            }
            /* free original array of pointers to all free tiles */
            free(tiles);
            free(index);
            free(walls);
            tiles=tiles_tmp;
            index=index_tmp;
            walls=walls_tmp;
            n_free_eff=n_free_eff-n_set;
         }

            /* update n_occupied for each effector grid */
	    for (i=0;i<rp->membership->nbits;i++)
	    {
	      if (get_bit(rp->membership,i))
	      {
		sg=objp->wall_p[i]->grid;
		if (sg!=NULL)
		{
		  sg->n_occupied=0;
		  for (j=0;j<sg->n_tiles;j++) {
		    if (sg->mol[j]!=NULL) {
		      sg->n_occupied++;
		    }
		  }
                }
              }
            }
          } /* end while(world->chkpt_init) */
        }
        effdp=effdp->next;
      }
      /* free array of pointers to all free tiles */
      if (tiles!=NULL) {
        free(tiles);
      }
      if (index!=NULL) {
        free(index);
      }
      if (walls!=NULL) {
        free(walls);
      }
      rlp=rlp->next;
    }
    no_printf("Done initialize effectors by number.\n");
    return(0);
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
struct void_list* rel_expr_grab_obj(struct release_evaluator *root,struct mem_helper *voidmem)
{
  struct void_list *vl = NULL;
  struct void_list *vr = NULL;
  
  if (root->left != NULL)
  {
    if (root->op&REXP_LEFT_REGION)
    {
      vl = mem_get(voidmem);
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
      vr = mem_get(voidmem);
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

struct object** find_unique_rev_objects(struct release_evaluator *root,int *n)
{
  struct object **o_array;
  struct void_list *vp,*vq;
  struct mem_helper *voidmem;
  int i;
  
  voidmem = create_mem(sizeof(struct void_list), 1024);
  
  vp = rel_expr_grab_obj(root,voidmem);
  if (vp==NULL) return NULL;
  
  vp = void_list_sort(vp);
  
  for (i=1,vq=vp ; vq!=NULL && vq->next!=NULL ; vq=vq->next , i++)
  {
    while (vq->data == vq->next->data)
    {
      vq->next = vq->next->next;
      if (vq->next==NULL) break;
    }
  }
  
  if (vq==NULL) i--;
  *n = i;
  
  o_array = (struct object**)malloc(i*sizeof(struct object*));
  if (o_array==NULL) return NULL;
  
  for (i=0,vq=vp ; vq!=NULL ; vq=vq->next,i++)
  {
    o_array[i] = (struct object*)vq->data;
  }
  
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

int eval_rel_region_expr(struct release_evaluator *expr,int n,struct object **objs,struct bit_array **result,int *n_refinements)
{
  int i;
  char bit_op;
  
  if (expr->left!=NULL)
  {
    if (expr->op&REXP_INCLUSION)
    {
      if (expr->right==NULL) return 1;  /* Should always have two arguments */
      i=eval_rel_region_expr(expr->left,n,objs,result,n_refinements);
      if (i==1) return 1;
      *n_refinements += 1;
      /* Just ignore right-hand argument; we'll mark that we should look at it later */
    }
    else
    {
      if (expr->op&REXP_LEFT_REGION)
      {
	i = void_array_search((void**)objs,n,((struct region*)(expr->left))->parent);
	result[i] = duplicate_bit_array( ((struct region*)(expr->left))->membership );
	if (result[i]==NULL) return 1;
      }
      else
      {
	i = eval_rel_region_expr(expr->left,n,objs,result,n_refinements);
	if (i) return 1;
      }
      
      if (expr->right==NULL)
      {
	if (expr->op&REXP_NO_OP) return 0;
	else return 1;
      }
      
      if (expr->op&REXP_RIGHT_REGION)
      {
	i = void_array_search((void**)objs,n,((struct region*)(expr->right))->parent);
	if (result[i]==NULL)
	{
	  result[i] = duplicate_bit_array( ((struct region*)(expr->right))->membership );
	  if (result[i]==NULL) return 1;
	}
	else
	{
	  if (expr->op&REXP_UNION) bit_op = '|';
	  else if (expr->op&REXP_SUBTRACTION) bit_op = '-';
	  else if (expr->op&REXP_INTERSECTION) bit_op = '&';
	  else return 1;
  
	  bit_operation(result[i],((struct region*)(expr->right))->membership,bit_op);
	}
      }
      else
      {
	struct bit_array *res2[n];
	for (i=0;i<n;i++) res2[i]=NULL;
	
	i = eval_rel_region_expr(expr->right,n,objs,res2,n_refinements);
	if (i) return 1;
	
	for (i=0;i<n;i++)
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

int init_rel_region_data_2d(struct release_region_data *rrd)
{
  int i,j,k;
  struct polygon_object *po;

  rrd->owners = find_unique_rev_objects(rrd->expression , &(rrd->n_objects));
  if (rrd->owners==NULL)
  {
    fprintf(world->err_file,"File '%s', Line %ld: Error, cannot find any objects for region release\n", __FILE__, (long)__LINE__);
    return 1;
  }
  
  rrd->in_release = (struct bit_array**)malloc(rrd->n_objects*sizeof(struct bit_array*));
  if (rrd->in_release==NULL)
  {
    fprintf(world->err_file,"File '%s', Line %ld: Out of memory creating region lists for 2D region releases\n", __FILE__, (long)__LINE__);
    return 1;
  }
  for (i=0;i<rrd->n_objects;i++) rrd->in_release[i]=NULL;
  
  rrd->refinement = 0;
  i = eval_rel_region_expr(rrd->expression,rrd->n_objects,rrd->owners,rrd->in_release,&rrd->refinement);
  if (i)
  {
    fprintf(world->err_file,"File '%s', Line %ld: Error, could not evaluate region expression.\n", __FILE__, (long)__LINE__);
    return 1;
  }
  for (i=0;i<rrd->n_objects;i++)
  {
    if (rrd->owners[i]==NULL)
    {
      fprintf(world->err_file,"File '%s', Line %ld: Object %d of %d in region expression was not found!\n", __FILE__, (long)__LINE__, i+1,rrd->n_objects);
      return 1;
    }
  }
  
  rrd->walls_per_obj = (int*)malloc(rrd->n_objects*sizeof(int));
  if (rrd->walls_per_obj==NULL)
  {
    fprintf(world->err_file,"File '%s', Line %ld: Error, out of memory creating wall counts for 2D region releases\n", __FILE__, (long)__LINE__);
    return 1;
  }
  
  rrd->n_walls_included=0;
  for (i=0;i<rrd->n_objects;i++)
  {
    if (rrd->in_release[i]==NULL) rrd->walls_per_obj[i]=0;
    else rrd->walls_per_obj[i] = count_bits(rrd->in_release[i]);
    rrd->n_walls_included += rrd->walls_per_obj[i];
  }
  
  rrd->cum_area_list = (double*)malloc(rrd->n_walls_included*sizeof(double));
  rrd->wall_index = (int*)malloc(rrd->n_walls_included*sizeof(int));
  rrd->obj_index = (int*)malloc(rrd->n_walls_included*sizeof(int));
  if (rrd->cum_area_list==NULL || rrd->wall_index==NULL || rrd->obj_index==NULL)
  {
    fprintf(world->err_file,"File '%s', Line %ld: Out of memory creating area lists for 2D region releases\n", __FILE__, (long)__LINE__);
    return 1;
  }
  
  j = 0;
  for (i=0;i<rrd->n_objects;i++)
  {
    if (rrd->walls_per_obj[i]==0) continue;
    k = rrd->owners[i]->object_type;
    if (k != POLY_OBJ && k != BOX_OBJ)
    {
      fprintf(world->err_file,"File '%s', Line %ld: found a region on something that isn't a box or polygon object?\n", __FILE__, (long)__LINE__);
      return 1;
    }
    po = (struct polygon_object*)(rrd->owners[i]->contents);
    for (k=0;k<po->n_walls;k++)
    {
      if (get_bit(rrd->in_release[i],k))
      {
        rrd->cum_area_list[j] = rrd->owners[i]->wall_p[k]->area;
        rrd->obj_index[j] = i;
        rrd->wall_index[j] = k;
        j++;
      }
    }
  }
  
  for (i=1;i<rrd->n_walls_included;i++)
  {
    rrd->cum_area_list[i] += rrd->cum_area_list[i-1];
  }
  
  return 0;
}


/***************************************************************************
create_region_bbox:
  In: a region
  Out: pointer to a 2-element array contining the LLF and URB corners of
       a bounding box around the region, or NULL if out of memory.
***************************************************************************/

struct vector3* create_region_bbox(struct region *r)
{
  int i,j,k;
  struct vector3 *bbox;
  struct vector3 *v;
  
  bbox = (struct vector3*) malloc(2*sizeof(struct vector3));
  if (bbox==NULL) {
    fprintf(world->err_file, "File '%s', Line %ld: Out of memory while creating region bounding box.\n", __FILE__, (long)__LINE__);
    return NULL;
  }
  
  j=0;
  for (i=0;i<r->membership->nbits;i++)
  {
    if (get_bit(r->membership,i))
    {
      if (!j)
      {
        bbox[0].x = bbox[1].x = r->parent->wall_p[i]->vert[0]->x;
        bbox[0].y = bbox[1].y = r->parent->wall_p[i]->vert[0]->y;
        bbox[0].z = bbox[1].z = r->parent->wall_p[i]->vert[0]->z;
      }
      for (k=0;k<3;k++)
      {
        v = r->parent->wall_p[i]->vert[k];
        if (bbox[0].x > v->x) bbox[0].x = v->x;
        else if (bbox[1].x < v->x) bbox[1].x = v->x;
        if (bbox[0].y > v->y) bbox[0].y = v->y;
        else if (bbox[1].y < v->y) bbox[1].y = v->y;
        if (bbox[0].z > v->z) bbox[0].z = v->z;
        else if (bbox[1].z < v->z) bbox[1].z = v->z;
      }
      j++;
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

int eval_rel_region_bbox(struct release_evaluator *expr,struct vector3 *llf,struct vector3 *urb)
{
  int i;
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
        {
          fprintf(world->err_file,"File '%s', Line %ld: Error--cannot release a 3D molecule inside an unclosed region\n", __FILE__, (long)__LINE__);
          return 1;
        }
      }
      
      if (r->bbox==NULL)
      {
        r->bbox = create_region_bbox(r);
        if (r->bbox==NULL) return 1;
      }
      
      llf->x = r->bbox[0].x;
      llf->y = r->bbox[0].y;
      llf->z = r->bbox[0].z;
      urb->x = r->bbox[1].x;
      urb->y = r->bbox[1].y;
      urb->z = r->bbox[1].z;
    }
    else
    {
      i = eval_rel_region_bbox(expr->left,llf,urb);
      if (i) return 1;
    }
    
    if (expr->right==NULL)
    {
      if (expr->op&REXP_NO_OP) return 0;
      else return 1;
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
          {
            fprintf(world->err_file,"File '%s', Line %ld: Error--cannot release a 3D molecule inside an unclosed region.\n", __FILE__, (long)__LINE__);
            return 1;
          }
        }
        
        if (r->bbox==NULL)
        {
          r->bbox = create_region_bbox(r);
          if (r->bbox==NULL) return 1;          
        }

        llf2.x = r->bbox[0].x;
        llf2.y = r->bbox[0].y;
        llf2.z = r->bbox[0].z;
        urb2.x = r->bbox[1].x;
        urb2.y = r->bbox[1].y;
        urb2.z = r->bbox[1].z;
      }
      else
      {
        i = eval_rel_region_bbox(expr->right,&llf2,&urb2);
        if (i) return 1;
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
      else return 1;
    }
  }
  else return 1;  /* Left should always have something! */
  
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

int init_rel_region_data_3d(struct release_region_data *rrd)
{
  int i;
  
  rrd->n_walls_included = 0;
  
  i = eval_rel_region_bbox(rrd->expression,&(rrd->llf),&(rrd->urb));

  if (i) return 1;

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

void output_relreg_eval_tree(FILE *f,char *prefix,char cA,char cB,struct release_evaluator *expr)
{
  int l = strlen(prefix);
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

int init_releases()
{
  struct release_event_queue *req;
  struct abstract_element *ae;
  struct schedule_helper *sh;
  int i,j;
  
  for (sh=world->releaser ; sh!=NULL ; sh=sh->next_scale)
  {
    for (i=-1;i<sh->buf_len;i++)
    {
      for ( ae = (i==-1)?sh->current:sh->circ_buf_head[i] ; ae!=NULL ; ae=ae->next )
      {
        req = (struct release_event_queue*)ae;
        if (req->release_site->release_shape == SHAPE_REGION)
        {
          if (req->release_site->mol_type == NULL){
              fprintf(world->err_file,"ERROR: molecule name is not specified for the region release site.\n");

          }
          if ((req->release_site->mol_type->flags & NOT_FREE) == 0)
          {
            j = init_rel_region_data_3d(req->release_site->region_data);
            if (j==-1)
            {
              fprintf(world->err_file,"File '%s', Line %ld: Region release site is empty!  Ignoring!  Evaluation tree:\n", __FILE__, (long)__LINE__);
              output_relreg_eval_tree(world->err_file," ",' ',' ',req->release_site->region_data->expression);
              req->release_site->release_number_method=CONSTNUM;
              req->release_site->release_number=0;
            }
            else if (j)
	    {
	      fprintf(world->err_file,"File '%s', Line %ld: Error initializing region release for molecule %s\nEvaluation tree (3D release):\n", __FILE__,(long)__LINE__, req->release_site->mol_type->sym->name);
	      output_relreg_eval_tree(world->err_file," ",' ',' ',req->release_site->region_data->expression);
	      return 1;
	    }
          }
          else
          {
            j = init_rel_region_data_2d(req->release_site->region_data);
            if (j)
	    { 
	      fprintf(world->err_file,"File '%s', Line %ld: Error initializing region release for molecule %s\nEvaluation tree (2D release):\n",__FILE__, (long)__LINE__, req->release_site->mol_type->sym->name);
	      output_relreg_eval_tree(world->err_file," ",' ',' ',req->release_site->region_data->expression);
	      return 1;
	    }
          }
        }else if (req->release_site->release_shape != SHAPE_LIST){
            if(req->release_site->mol_type == NULL){
               fprintf(world->err_file, "ERROR: molecule name is not specified for the release site.\n");
               return 1;
            }
            if(req->release_site->diameter == NULL){
               fprintf(world->err_file, "ERROR: diameter for the geometrical shape release site is not specified.\n");
               return 1;
            }
        }else{
           /* this check should be for SHAPE_LIST release only */
            if(req->release_site->mol_list == NULL){
               fprintf(world->err_file, "ERROR: molecule positions for the SHAPE_LIST release site are not specified.\n");
               return 1;
            }
  

        }
      }
    }
  }
  
  return 0;
}

