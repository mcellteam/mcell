
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

#include "mcell_structs.h"
#include "strfunc.h"
#include "vector.h"
#include "rng.h"
/*
#include "geom_util.h"
*/
#include "sym_table.h"
#include "init.h"

#ifdef DEBUG
#define no_printf printf
#endif

extern struct volume *world;

/**
 * Prints out author and grant credits.
 * When invoked with the -info option, this function prints to log_file
 * a header like:

 <pre>
  MCell (tm) Version 2.50x  06/05/2001  Running on dalton.salk.edu

  Copyright (C) 1997,1998,1999 by The Salk Institute & Cornell University
  Co-authored by Thomas M. Bartol Jr. & Joel R. Stiles
  Acknowledgements:
    The authors thank Edwin E. Salpeter for input on theory and
    algorithm design, Terrence J. Sejnowski for development input
    and support (NSF Grant IBN-9603611), and Miriam M. Salpeter
    for fostering quantitative experimental applications.
    Additional support from NIH Grant K08NS01776 (Joel R. Stiles).
 </pre>

 */
void init_credits(void)
{
  FILE *log_file;
  unsigned int seed;
  double ran_vec[100];
  time_t the_time;
  char *institute[2],*author[2];
 
  log_file=world->log_file;
  time(&the_time);
  seed=(unsigned int)the_time;
  ran4_init(&seed);
  ran4(&seed,ran_vec,100,1.0);
  if (ran_vec[99]<0.5) {
    institute[0]=my_strdup("The Salk Institute");
    institute[1]=my_strdup("& Cornell University");
    author[0]=my_strdup("Thomas M. Bartol Jr.");
    author[1]=my_strdup("& Joel R. Stiles");
  }
  else {
    institute[0]=my_strdup("Cornell University");
    institute[1]=my_strdup("& The Salk Institute");
    author[0]=my_strdup("Joel R. Stiles");
    author[1]=my_strdup("& Thomas M. Bartol Jr.");
  }
  
    fprintf(log_file,"  Copyright (C) 1997,1998,1999 by %s %s\n",institute[0],institute[1]);
    fprintf(log_file,"  Co-authored by %s %s\n",author[0],author[1]);
    fprintf(log_file,"  Acknowledgements:\n");
    fprintf(log_file,"    The authors thank Edwin E. Salpeter for input on theory and\n");
    fprintf(log_file,"    algorithm design, Terrence J. Sejnowski for development input\n");
    fprintf(log_file,"    and support (NSF Grant IBN-9603611), and Miriam M. Salpeter\n");
    fprintf(log_file,"    for fostering quantitative experimental applications.\n");
    fprintf(log_file,"    Additional support from NIH Grant K08NS01776 (Joel R. Stiles).\n\n");

  return;
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
  FILE *log_file;
  struct sym_table *gp;
  double fact;
  int i;
  int *intp;
#include "seed_array.h"

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

  /* ============ by Erhan Gokcay 5/8/2002 =========== */
  /* This is set in the main program. */
  /*chkpt_init=1; */
  /* ================================================= */
  world->chkpt_flag=0;
  world->molecule_prefix_name=NULL;
  world->random_number_use=0;
  world->ray_voxel_tests=0;
  world->ray_polygon_tests=0;
  world->diffusion_steps=0;
  world->sim_elapsed_time=0;
  world->chkpt_elapsed_time=0;
  world->chkpt_elapsed_time_start=0;
  world->it_time=0;
  world->time_unit=0;
  world->start_time=0;
  world->current_time=0;
  world->current_start_time=0;
  world->effector_grid_density=1;
  world->length_unit=1.0/sqrt(world->effector_grid_density);
  world->max_diffusion_step=0;
  world->radial_directions=16384;
  world->radial_subdivisions=1024;
  world->fully_random=0;
  world->num_directions=world->radial_directions;
  world->r_num_directions=1.0/world->num_directions;
  world->r_step=NULL;
  world->d_step=NULL;
  world->n_release_events=0;

  if ((world->factorial_r=(double *)malloc(101*sizeof(double)))==NULL) {
    fprintf(log_file,"MCell: could not store factorial array\n");
    return(1);
  }
  world->factorial_r[0]=1;
  world->factorial_r[1]=1;
  fact=1;
  for (i=2;i<101;i++) {
    fact=fact*i;
    world->factorial_r[i]=1.0/fact;
  }
  if (world->seed_seq < 1 || world->seed_seq > 3000) {
    fprintf(log_file,"MCell: error random sequence number not in range 1 to 3000\n");
    return(1);
  }

  world->seed=seed_array[world->seed_seq-1];
  world->init_seed = world->seed;
  fprintf(log_file,"MCell[%d]: random sequence: %d  seed: %d\n", world->procnum,world->seed_seq,world->seed);
  fflush(log_file);
  ran4_init(&world->seed);

  world->main_sym_table=init_symtab(HASHSIZE);

  if ((gp=store_sym("WORLD_OBJ",OBJ,world->main_sym_table))==NULL) {
    fprintf(log_file,"MCell: could not store world root object\n");
    return(1);
  }
  world->root_object=(struct object *)gp->value;
  world->root_object->object_type=META_OBJ;
  world->root_object->last_name="";

  if ((gp=store_sym("WORLD_INSTANCE",OBJ,world->main_sym_table))==NULL) {
    fprintf(log_file,"MCell: could not store world root instance\n");
    return(1);
  }
  world->root_instance=(struct object *)gp->value;
  world->root_instance->object_type=META_OBJ;
  world->root_instance->last_name="";

  if ((gp=store_sym("DEFAULT_RELEASE_PATTERN",RPAT,world->main_sym_table))
      ==NULL) {
    fprintf(log_file,"MCell: cannot store default release pattern");
    return(1);
  }
  world->default_release_pattern=(struct release_pattern *)gp->value;
  world->default_release_pattern->delay=0;
  world->default_release_pattern->release_interval=2;
  world->default_release_pattern->train_interval=1;
  world->default_release_pattern->train_duration=1;
  world->default_release_pattern->number_of_trains=1;

  world->count_list=NULL;
  world->output_list=NULL;
  world->release_event_queue_head=NULL;
  world->tot_mols=0;
  world->viz_obj_head=NULL;
  world->viz_mode=0;
  world->frame_data_head=NULL;

  if ((world->count_zero=(struct count_list *)malloc
       (sizeof(struct count_list)))==NULL) {
    fprintf(log_file,"MCell: cannot store counter data\n");
    return(1); 
  }
  if (!(intp=(int *)malloc(sizeof(int)))) {
    fprintf(log_file,"MCell: cannot store counter data\n");
    return(1);
  }
  *intp=0;
  world->count_zero->freq=1;
  world->count_zero->frame_index=1;
  world->count_zero->n_output=0;
  world->count_zero->reset_flag=0;
  world->count_zero->update_flag=0;
  world->count_zero->index_type=TIME_STAMP_VAL;
  world->count_zero->n_data=1;
  world->count_zero->data_type=INT;
  world->count_zero->temp_data=(void *)intp;
  world->count_zero->final_data=(void *)intp;
  world->count_zero->operand1=NULL;
  world->count_zero->operand2=NULL;
  world->count_zero->oper='\0';
  world->count_zero->next=world->count_list;
  world->count_list=world->count_zero;
  
  world->releaser = create_scheduler(1.0,100.0,100,0.0);

  /* Parse the MDL file: */
  no_printf("Node %d parsing MDL file %s\n",world->procnum,world->mdl_infile_name);
  fflush(stderr);
  if (mdlparse_init(world)) {
    fprintf(log_file,"MCell: error parsing file: %s\n",world->curr_file);
    return(1);
  }
  no_printf("Done parsing MDL file: %s\n",world->mdl_infile_name);
  fflush(stderr);
  
  /* Set up the array of species */
  if (init_species())
  {
    fprintf(log_file,"MCell: error initializing species\n");
    return(1);
  }
  no_printf("Done setting up species.\n");

  /* Initialize the geometry */
  if (init_geom()) {
    fprintf(log_file,"MCell: error initializing geometry\n");
    return(1);
  }
  
  no_printf("Done setting up geometry.\n");
  
  if (init_partitions()) {
    fprintf(log_file,"MCell: error initializing partitions.\n");
    return(1);
  }
  


  /* Decompose the space into subvolumes */
/*
  if (decompose_volume(volume,wall_head)) {
    fprintf(log_file,"MCell: error decomposing volume\n");
    return(1);
  }
  resolve_contiguity(cmprt_head);
  if (init_effector_table(wall_head)) {
    fprintf(log_file,"MCell: error initializing effector table\n");
    return(1);
  }
  fflush(log_file);
  if (init_rx_table(rx_head)) {
    fprintf(log_file,"MCell: error initializing rx table\n");
    return(1);
  }
  if (init_release_event_table(release_event_queue_head)) {
    fprintf(log_file,"MCell: error initializing release event table\n");
    return(1);
  }
  

  if (chkpt_infile) {
    if ((chkpt_infs=fopen(chkpt_infile,"rb"))==NULL) {
      chkpt_seq_num=1;
    }
    else {
      fprintf(log_file,"MCell: reading from checkpoint file %s\n",chkpt_infile);
      if(read_chkpt(chkpt_infs)) {
	fprintf(log_file,"MCell: error reading from checkpoint file %s\n",chkpt_infile);
	return(1);
      }
      fclose(chkpt_infs);
    }
  }
  else {
    chkpt_seq_num=1;
  }
*/

  /**
   *Initialize the frame date list for the visualization 
   *and reaction output.
   **/
/*
  init_frame_data_list(frame_data_head);
  init_reaction_list(reaction_data_head);
*/

  no_printf("Done initializing simulation\n");
  fflush(log_file);
  return(0);
}


int init_species(void)
{
  int i;
  int count = 0;
  struct sym_table *gp;
  
/* FIX ME WHEN WE HAVE REACTIONS */
  world->n_reactions = 0;
  world->reaction_hash = (struct rxn**)malloc(sizeof(struct rxn*));
  world->hashsize = 1;
  world->reaction_hash[0] = (struct rxn*)malloc(sizeof(struct rxn));
  world->reaction_hash[0]->next = NULL;
  world->reaction_hash[0]->n_reactants = 1;
  world->reaction_hash[0]->n_pathways = 1;
  world->reaction_hash[0]->product_idx = (u_int*)malloc(sizeof(int));
  world->reaction_hash[0]->product_idx[0] = 1;
  world->reaction_hash[0]->cum_rates = (double*)malloc(sizeof(double));
  world->reaction_hash[0]->cum_rates[0] = 1e-100;
  world->reaction_hash[0]->cat_rates = world->reaction_hash[0]->cum_rates;
  world->reaction_hash[0]->players = (struct species**)malloc(2*sizeof(struct species*));
  world->reaction_hash[0]->players[0] = NULL;
  world->reaction_hash[0]->players[1] = NULL;
  world->reaction_hash[0]->geometries = (short*)malloc(2*sizeof(short));
  world->reaction_hash[0]->geometries[0] = 0;
  world->reaction_hash[0]->geometries[1] = 0;
  world->reaction_hash[0]->geometries = (short*)malloc(2*sizeof(short));
  world->reaction_hash[0]->fates = (byte*)malloc(sizeof(byte));
  world->reaction_hash[0]->fates[0] = 0;
/* END FIX ME*/
  
  for (i=0;i<HASHSIZE;i++)
  {
    for (gp = world->main_sym_table[i] ; gp != NULL ; gp = gp->next)
    {    
      if (gp->sym_type==MOL) count++;
    }
  }
  
  world->n_species = count;  printf("Found %d species!\n",world->n_species);
  world->species_list = (struct species**)malloc(sizeof(struct species*)*world->n_species);

  count = 0;
  for (i=0;i<HASHSIZE;i++)
  {
    for (gp = world->main_sym_table[i] ; gp != NULL ; gp = gp->next)
    {    
      if (gp->sym_type==MOL)
      {
        world->species_list[count] = (struct species*) gp->value;
        world->species_list[count]->hashval &= world->hashsize-1;
      }
    }
  }
  
  return 0;
}



/* This is just a placeholder for now--make one giant partition. */
int init_partitions(void)
{
  int i,j,k;
  struct subvolume *sv;
  
  world->n_axis_partitions = 2;
  world->x_partitions = (double*)malloc(sizeof(double)*world->n_axis_partitions);
  world->y_partitions = (double*)malloc(sizeof(double)*world->n_axis_partitions);
  world->z_partitions = (double*)malloc(sizeof(double)*world->n_axis_partitions);
  world->x_partitions[0] = - GIGANTIC;
  world->x_partitions[1] = GIGANTIC;
  world->y_partitions[0] = - GIGANTIC;
  world->y_partitions[1] = GIGANTIC;
  world->z_partitions[0] = - GIGANTIC;
  world->z_partitions[1] = GIGANTIC;
  
  
  world->n_fine_partitions = world->n_axis_partitions;
  world->x_fineparts = world->x_partitions;
  world->y_fineparts = world->y_partitions;
  world->z_fineparts = world->z_partitions;
  
  world->n_waypoints = 1;
  world->waypoints = (struct waypoint*)malloc(sizeof(struct waypoint*)*world->n_waypoints);
  
  world->n_subvols = world->n_waypoints;
  world->subvol = (struct subvolume*)malloc(sizeof(struct subvolume*)*world->n_subvols);
  for (i=0;i<world->n_axis_partitions-1;i++)
  for (j=0;j<world->n_axis_partitions-1;j++)
  for (k=0;k<world->n_axis_partitions-1;k++)
  {
    sv = & (world->subvol[k + world->n_axis_partitions*(j + k*world->n_axis_partitions)]);
    sv->wall_head = NULL;
    sv->wall_tail = NULL;
    sv->wall_count = 0;
    sv->mol_head = NULL;
    sv->mol_count = 0;
    
    sv->index = -1;
    
    sv->llf.x = i;
    sv->llf.y = j;
    sv->llf.z = k;
    sv->urb.x = i+1;
    sv->urb.y = j+1;
    sv->urb.z = k+1;


    
    sv->is_bsp = 0;
    
    sv->neighbor[0] = NULL;
    sv->neighbor[1] = NULL;
    sv->neighbor[2] = NULL;
    sv->neighbor[3] = NULL;
    sv->neighbor[4] = NULL;
    sv->neighbor[5] = NULL;
    
    sv->mem = (struct storage*)malloc(sizeof(struct storage));
    
    sv->mem->list = create_mem(sizeof(struct wall_list),50);
    sv->mem->mol  = create_mem(sizeof(struct molecule),50);
    sv->mem->smol  = create_mem(sizeof(struct surface_molecule),50);
    sv->mem->gmol  = create_mem(sizeof(struct grid_molecule),50);
    sv->mem->wall = create_mem(sizeof(struct wall),50);
    sv->mem->coll = create_mem(sizeof(struct collision),50);
    
    sv->mem->timer = create_scheduler(1.0,100.0,100,0.0);
    sv->mem->current_time = 0.0;
    sv->mem->max_timestep = 1000.0;
  }
  
  world->binning = 0;
  world->lookup = NULL;
  
  world->collide_hashmask = 0xFFFF;
  world->collide_hash = (struct counter**)malloc(sizeof(struct counter*)*(world->collide_hashmask+1));
  
  return 0;
}



/**
 * Initializes the geometry of the world.
 * Calls instance_obj() to instantiate all physical objects.
 * (Meta objects, box objects, polygon objects and release sites)
 * Populates viz_obj list vizp and lig_count_ref list lcrp.
 */
int init_geom(void)
{
  FILE *log_file;
  double tm[4][4];
  double vol_infinity;
  struct release_event_queue *req,*rqn;

  no_printf("Initializing physical objects\n");
  log_file=world->log_file;
  vol_infinity=sqrt(DBL_MAX)/4;
  world->bb_min.x=vol_infinity;
  world->bb_min.y=vol_infinity;
  world->bb_min.z=vol_infinity;
  world->bb_max.x=-vol_infinity;
  world->bb_max.y=-vol_infinity;
  world->bb_max.z=-vol_infinity;
  init_matrix(tm);
  
  compute_bb(world->root_instance,tm,NULL);
  if (world->procnum == 0) {
    fprintf(log_file,"MCell: world bounding box =\n");
    fprintf(log_file,"         [ %.9g %.9g %.9g ] [ %.9g %.9g %.9g ]\n",
      world->bb_min.x,world->bb_min.y,world->bb_min.z,world->bb_max.x,
      world->bb_max.y,world->bb_max.z);
  }
  
  if (instance_obj(world->root_instance,tm,NULL,NULL,NULL)) {
    return(1);
  }
  
/* Stick the queue of release events in a scheduler */
  req = world->release_event_queue_head;  
  while(req != NULL)
  {
    rqn = req->next;
    schedule_add(world->releaser , req);
    req = rqn;
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
int instance_obj(struct object *objp, double (*im)[4], struct viz_obj *vizp, struct lig_count_ref *lcrp, char *sub_name)
{
  FILE *log_file;
  struct object *child_objp;
  double tm[4][4];
  unsigned short l,m,n;
  char *tmp_name;

  log_file=world->log_file;
  l=4;
  m=4;
  n=4;
  mult_matrix(objp->t_matrix,im,tm,l,m,n);
  if (vizp==NULL) {
    vizp=objp->viz_obj;
  }
  if (lcrp==NULL) {
    lcrp=objp->lig_count_ref;
  }
  if (sub_name!=NULL) { 
    if (strcmp(sub_name,"")==0) {
      tmp_name=my_strdup("");
    }
    else {
      tmp_name=my_strcat(sub_name,".");              
    }
    sub_name=my_strcat(tmp_name,objp->last_name);    
    free((void *)tmp_name);
  }
  else {
    sub_name=my_strdup(objp->last_name);    
  }

  switch (objp->object_type) {
  case META_OBJ:
    no_printf("Meta object %s instanced\n",sub_name);
    fflush(log_file);
    child_objp=objp->first_child;
    while (child_objp!=NULL) {
      if (instance_obj(child_objp,tm,vizp,lcrp,sub_name)) {
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
/*
  case BOX_OBJ:
    no_printf("Box object %s instanced\n",sub_name);
    fflush(log_file);
    if (instance_polygon_object(objp,tm,vizp,lcrp,sub_name)) {
      return(1);
    }
    break;
  case POLY_OBJ:
    no_printf("Polygon list object %s instanced\n",sub_name);
    fflush(log_file);
    if (instance_polygon_object(objp,tm,vizp,lcrp,sub_name)) {
      return(1);
    }
    break;
*/
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
  double location[1][4];
  unsigned short l,m,n;

  log_file=world->log_file;
  rsop=(struct release_site_obj *)objp->contents;
  l=1;
  m=4;
  n=4;
  location[0][0]=rsop->location->x;
  location[0][1]=rsop->location->y;
  location[0][2]=rsop->location->z;
  location[0][3]=1;
  mult_matrix(location,im,location,l,m,n);

/*
  if ( ((location[0][0] >= volume->x_partitions[0])
			  && (location[0][0] <= volume->x_partitions[volume->n_x_subvol]))
    && ((location[0][1] >= volume->y_partitions[0])
			  && (location[0][1] <= volume->y_partitions[volume->n_y_subvol]))
    && ((location[0][2] >= volume->z_partitions[0])
			  && (location[0][2] <= volume->z_partitions[volume->n_z_subvol])) ){
*/

	  no_printf("Instancing release site object %s\n",objp->sym->name);
	  fflush(log_file);
	  if ((reqp=(struct release_event_queue *)malloc
				  (sizeof(struct release_event_queue)))==NULL) {
		printf ("MCell: cannot store release event queue\n");
		return(1);
	  }

	  reqp->location.x=location[0][0]/world->length_unit;
	  reqp->location.y=location[0][1]/world->length_unit;
	  reqp->location.z=location[0][2]/world->length_unit;

	  reqp->release_site=rsop;
	  reqp->event_type=TRAIN_HIGH_EVENT;
	  reqp->event_time=rsop->pattern->delay;
	  reqp->event_counter=0;
	  reqp->train_high_time=0;
	  reqp->index=world->n_release_events++;
	  reqp->next=world->release_event_queue_head;
	  world->release_event_queue_head=reqp;

	  no_printf("Done instancing release site object %s\n",objp->sym->name);
/*
  } else {
	no_printf("Ignoring release site object %s\n",objp->sym->name);
  }
*/

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

  log_file=world->log_file;
  l=4;
  m=4;
  n=4;
  mult_matrix(objp->t_matrix,im,tm,l,m,n);
  if (sub_name!=NULL) { 
    if (strcmp(sub_name,"")==0) {
      tmp_name=my_strdup("");
    }
    else {
      tmp_name=my_strcat(sub_name,".");              
    }
    sub_name=my_strcat(tmp_name,objp->last_name);    
    free((void *)tmp_name);
  }
  else {
    sub_name=my_strdup(objp->last_name);    
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
/*
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
*/
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

  rsop=(struct release_site_obj *)objp->contents;

  l=1;
  m=4;
  n=4;
  location[0][0]=rsop->location->x;
  location[0][1]=rsop->location->y;
  location[0][2]=rsop->location->z;
  location[0][3]=1.0;
  mult_matrix(location,im,location,l,m,n);
  if (location[0][0]<world->bb_min.x) {
    world->bb_min.x=location[0][0];
  }
  if (location[0][1]<world->bb_min.y) {
    world->bb_min.y=location[0][1];
  }
  if (location[0][2]<world->bb_min.z) {
    world->bb_min.z=location[0][2];
  }
  if (location[0][0]>world->bb_max.x) {
    world->bb_max.x=location[0][0];
  }
  if (location[0][1]>world->bb_max.y) {
    world->bb_max.y=location[0][1];
  }
  if (location[0][2]>world->bb_max.z) {
    world->bb_max.z=location[0][2];
  }

  return(0);
}






/* **************************************************************** */

# if 0






/**
 * Instantiates a polygon_object.
 * Creates a new polygon object from a template polygon_object or box object
 * as defined in the MDL file after applying the necessary geometric
 * transformations (scaling, rotation and translation).
 * <br>
 * <b>Note:</b> Box objects (box_poly) are also instantiated using this
 * function. Ultimately, even boxes are turned into polygon_object's.
 */
int instance_polygon_object(struct object *objp, double (*im)[4], struct viz_obj *vizp, struct lig_count_ref *obj_lcrp, char *full_name)
{
  struct polygon_object *pop;
  struct ordered_poly *opp;
  struct box_poly *bpp;
  struct vertex_list *vlp;
  struct vector3 *corner,**face,*rect_vert[4],**vert,*vp1,*vp2,v1,v2,vx,vy,vz;
  struct vector3 *vertex_normal,**face_vertex_normal;
  struct vector3 ab,bc,ac;
  struct cmprt_data *cdp;
  struct cmprt_data_list *cdlp;
  struct wall *wp;
  struct wall_list *wlp;
  struct effector *ep;
  struct rx *rx;
  struct eff_dat *effdp,*dup_effdp,**eff_prop;
  struct lig_count_ref *poly_lcrp;
  struct region *rp;
  struct region_list *rlp,*rlp2,*reg_eff_num_head,*rlp3,*reg_count_head;
  struct region_list *rlp4, *lig_hit_count;
  struct lig_hit_counter **lig_hit;
  struct element_list *elp;
  struct reg_counter_ref *rcrp;
  struct reg_counter_ref_list *rcrlp;
  struct counter_hash_table **countertab;
  double dx,dy,dz;
  double p[1][4],origin[1][4];
  double area,total_area;
  int i,j,k,cnt,n_verts,n_polys,n_element_verts,index_0,index_1,index_n;
  int done;
  unsigned short l,m,n;
  char *obj_name;
  byte compute_vertex_normals,reg_eff_num,reg_count_num;
  byte lig_hit_flag;
  byte throw_away = 0;
  byte wall_on_region_flag;

  pop=(struct polygon_object *)objp->contents;
  n_polys=pop->n_polys;
  l=1;
  m=4;
  n=4;
  total_area=0;
  obj_name=my_strdup(full_name);


   /* There's two cases here --
   * 	POLY_OBJ - an object comprised of triangular polygons.
   * 	BOX_OBJ - object that's a cuboid box; it's got 8 corners.
   */
  switch (objp->object_type) {

  case POLY_OBJ:

	/* Allocate and initialize memory */
    compute_vertex_normals=0;
    if (pop->list_type==ORDERED_POLY) {
      opp=(struct ordered_poly *)pop->polygon_data;
      if (opp->normal!=NULL) {
        compute_vertex_normals=1;
      }
      n_verts=opp->n_verts;
    }
    if ((cdp=(struct cmprt_data *)malloc
	 (sizeof(struct cmprt_data)))==NULL) {
      return(1);
    }
    if ((cdp->lig_count=(int *)malloc
	 ((1+n_ligand_types)*sizeof(int)))==NULL){
      return(1);
    }
    if ((cdp->conc=(double *)malloc
	 ((1+n_ligand_types)*sizeof(double)))==NULL){
      return(1);
    }
    for (i=0;i<1+n_ligand_types;i++) {
      cdp->lig_count[i]=0;
      cdp->conc[i]=0;
    }
    if ((cdp->corner=(struct vector3 *)malloc
	 (n_verts*sizeof(struct vector3)))==NULL){
      return(1);
    }
    cdp->vertex_normal=NULL;
    if (compute_vertex_normals) {
      if ((cdp->vertex_normal=(struct vector3 *)malloc
	   (n_verts*sizeof(struct vector3)))==NULL){
        return(1);
      }
    }
    if ((cdp->normal=(struct vector3 *)malloc
	 (n_polys*sizeof(struct vector3)))==NULL){
      return(1);
    }
    if ((cdp->wall=(struct wall **)malloc
	 (n_polys*sizeof(struct wall *)))==NULL){
      return(1);
    }
    if ((cdp->neighbor=(struct cmprt_data **)malloc
	 (n_polys*sizeof(struct cmprt_data *)))==NULL){
      return(1);
    }
    for (i=0;i<n_polys;i++) {
      cdp->wall[i]=NULL;
      cdp->neighbor[i]=NULL;
    }
    cdp->sym=objp->sym;
    cdp->full_name=obj_name;
    cdp->cmprt_type=pop->list_type;
    cdp->fully_closed=pop->fully_closed;
    cdp->vm=0;
    cdp->n_corners=n_verts;
    cdp->n_walls=n_polys;
    cdp->wall_list=NULL;
    cdp->next=cmprt_head;
    cmprt_head=cdp;
    if (vizp!=NULL) {
      if ((cdlp=(struct cmprt_data_list *)malloc
	   (sizeof(struct cmprt_data_list)))==NULL) {
        return(1);
      }
      cdlp->cmprt_data = cdp;
      cdlp->next = vizp->cmprt_data_list;
      vizp->cmprt_data_list = cdlp;
    }
    corner=cdp->corner;
    vertex_normal=cdp->vertex_normal;
 
    for (i=0;i<n_verts;i++) {
      if (pop->list_type==ORDERED_POLY) {
        p[0][0]=opp->vertex[i]->x;
        p[0][1]=opp->vertex[i]->y;
        p[0][2]=opp->vertex[i]->z;
        p[0][3]=1.0;
      }
      mult_matrix(p,im,p,l,m,n);
      corner[i].x=p[0][0]/length_unit;
      corner[i].y=p[0][1]/length_unit;
      corner[i].z=p[0][2]/length_unit;

      if (compute_vertex_normals) {
        if (pop->list_type==ORDERED_POLY) {
          p[0][0]=opp->normal[i]->x;
          p[0][1]=opp->normal[i]->y;
          p[0][2]=opp->normal[i]->z;
          p[0][3]=1.0;
        }
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
    }

    if ((eff_prop=(struct eff_dat **)malloc
       (n_polys*sizeof(struct eff_dat *)))==NULL) {
      return(1);
    }

    for (i=0;i<n_polys;i++) {
      if (pop->list_type==ORDERED_POLY) {
        n_element_verts=opp->element_data[i].n_verts;
        index_0=opp->element_data[i].vertex_index[0];
        index_1=opp->element_data[i].vertex_index[1];
        index_n=opp->element_data[i].vertex_index[n_element_verts-1];
      }
      vectorize(&corner[index_0],&corner[index_1],&v1);
      vectorize(&corner[index_0],&corner[index_n],&v2);
      cross_prod(&v1,&v2,&cdp->normal[i]);
      if (vect_length(&v1)==0
          || vect_length(&v2)==0
          || vect_length(&cdp->normal[i])==0) {
        fprintf(log_file,"\nMCell: Warning -- Degenerate polygon found and automatically removed: %s %d\n\n",objp->sym->name,i);
        pop->side_stat[i]=0;
      }
      else {
        normalize(&cdp->normal[i]);
      }
      eff_prop[i]=NULL;
    }

    poly_lcrp=pop->lig_count_ref;
    while (poly_lcrp!=NULL) {
      if (strcmp(full_name,poly_lcrp->full_name)==0) {
        poly_lcrp->count_list->temp_data=(void *)&cdp->lig_count[poly_lcrp->type];
      }
      poly_lcrp=poly_lcrp->next;
    }

    /* prepend a copy of eff_dat for each element referenced in each region
       of this object to the eff_prop list for the referenced element */
    reg_eff_num_head=NULL;
    reg_count_head=NULL;
    lig_hit_count=NULL;
    rlp=objp->region_list;
    
    while (rlp!=NULL) {
      rp=rlp->region;
      reg_eff_num=0;
      reg_count_num=0;
      lig_hit_flag=0;
      rcrlp=rp->reg_counter_ref_list;
      if (rcrlp!=NULL) {
	reg_count_num=1;
      }
      lig_hit=rp->lig_hit_counter;
      if (lig_hit!=NULL) {
	lig_hit_flag=1;
      }
      elp=rp->element_list;
      while (elp!=NULL) {
        for (i=elp->begin;i<=elp->end;i++) {
          if (pop->side_stat[i]) {
            effdp=rp->eff_dat;
            while (effdp!=NULL) {
              if (effdp->quantity_type==EFFDENS) {
                if ((dup_effdp=(struct eff_dat *)malloc
                     (sizeof(struct eff_dat)))==NULL){
                  return(1);
                }
                dup_effdp->rx=effdp->rx;
                dup_effdp->quantity_type=effdp->quantity_type;
                dup_effdp->quantity=effdp->quantity;
                dup_effdp->orient=effdp->orient;
                dup_effdp->next=eff_prop[i];
                eff_prop[i]=dup_effdp;
              }
              else {
                reg_eff_num=1;
              }
              effdp=effdp->next;
            }
          }
        }
        elp=elp->next;
      }
      if (reg_eff_num) {
        if ((rlp2=(struct region_list *)malloc
             (sizeof(struct region_list)))==NULL){
          return(1);
        }
        rlp2->region=rp;
        rlp2->next=reg_eff_num_head;
        reg_eff_num_head=rlp2;
      }
      if (reg_count_num) {
        if ((rlp3=(struct region_list *)malloc
             (sizeof(struct region_list)))==NULL){
          return(1);
        }
        rlp3->region=rp;
        rlp3->next=reg_count_head;
        reg_count_head=rlp3;
      }
      if (lig_hit_flag) {
        if ((rlp4=(struct region_list *)malloc
             (sizeof(struct region_list)))==NULL){
          return(1);
        }
        rlp4->region=rp;
        rlp4->next=lig_hit_count;
        lig_hit_count=rlp4;
      }

      rlp=rlp->next;
    }
    k=0;
    for (i=0;i<n_polys;i++) {
      if (pop->side_stat[i]) {
          face_vertex_normal=NULL;
          if (pop->list_type==ORDERED_POLY) {
            n_element_verts=opp->element_data[i].n_verts;
            if ((face=(struct vector3 **)malloc
	         (n_element_verts*sizeof(struct vector3 *)))==NULL) {
              return(1);
            }
            if (compute_vertex_normals) {
              if ((face_vertex_normal=(struct vector3 **)malloc
	           (n_element_verts*sizeof(struct vector3 *)))==NULL) {
                return(1);
              }
            }
			throw_away = 0;
            for (j=0;j<n_element_verts;j++) {
              index_0=opp->element_data[i].vertex_index[j];
              face[j]=&corner[index_0];
              if (compute_vertex_normals) {
                face_vertex_normal[j]=&vertex_normal[index_0];
              }
			  
			  /* Check if any of the vertices is in the region for
			   * for this node */
              if ( ((face[j][0].x <= volume->x_partitions[0])
                   || (face[j][0].x >= volume->x_partitions[volume->n_x_subvol]))
               || ((face[j][0].y <= volume->y_partitions[0])
                   || (face[j][0].y >= volume->y_partitions[volume->n_y_subvol]))
               || ((face[j][0].z <= volume->z_partitions[0])
                   || (face[j][0].z >= volume->z_partitions[volume->n_z_subvol])) ){
                throw_away++;
			  }
            }

			/* If the entire wall is outside the region for this node, skip it */
            if (throw_away == n_element_verts) {
			  cdp->n_walls--;
			  free(face);
			  free(face_vertex_normal);
			  continue;
			}

          }
	  if ((wp=init_wall(face,face_vertex_normal,n_element_verts))==NULL) {
	    return(1);
	  }
	  if ((wlp=(struct wall_list *)malloc
	       (sizeof(struct wall_list)))==NULL) {
	    return(1);
	  }
	  wlp->wall=wp;
	  wlp->next=cdp->wall_list;
	  cdp->wall_list=wlp;
          pop->cmprt_side_map[i]=k;
	  cdp->wall[k++]=wp;
	  wp->next_wall=wall_head;
	  wall_head=wp;
          if (wp->n_vert==3) {
            wp->wall_shape=TRI_POLY;
          }
          else {
            wp->wall_shape=GEN_POLY;
          }

          vectorize(wp->vert[0],wp->vert[1],&ab);
          vectorize(wp->vert[1],wp->vert[2],&bc);
          if (wp->wall_shape==TRI_POLY) {
            vectorize(wp->vert[0],wp->vert[2],&ac);
            cross_prod(&ab,&ac,&v1);
            area=0.5*vect_length(&v1);
            wp->area=area;
            total_area=total_area+wp->area;
          }
          else {
            /* todo: compute area of general polygon */
            wp->area=0;
            total_area=total_area+wp->area;
          }

	  wp->parent_object=objp;   /* Lin-Wei Wu 6/25/03 */
	  wp->side=i;
	  wp->parent=cdp;
	  wp->wall_type=pop->lig_prop[i];
	  wp->effectors=NULL;
          if (objp->viz_state!=NULL) {
	    wp->viz_state=objp->viz_state[i];
          }
          else {
	    wp->viz_state=EXCLUDE_OBJ;
          }

	  /* Try to find out which region contains this wall
	   * Only when there are regions in this object and 
	   * there are region counters for this object.
	   * Else we set it to NULL 6/26/03*/
	  rlp=objp->region_list;
	  if ((rlp!=NULL)&&((reg_count_head!=NULL)||(lig_hit_count!=NULL))) {
	    wp->region_list=init_region_list_for_wall(rlp,i); 
  	  }
	  else {
	    wp->region_list=NULL;
	  }
	  
	  /* end of new searching */
	  if ((effdp=eff_prop[i])!=NULL) {
	    if (init_effectors_by_density(wp,effdp)) {
	      return(1);
	    }
	  }
      }
    }
	  if (lig_hit_count!=NULL) {
	    rlp=lig_hit_count;
	    while(rlp!=NULL) {
	      rlp4=rlp;
	      rlp=rlp->next;
	      free(rlp4);
	    }
	  }
  
    if (reg_eff_num_head!=NULL) {
      if (init_effectors_by_number(pop,cdp,reg_eff_num_head)) {
        return(1);
      }
      /* free region list created to hold regions populated by number */
      rlp=reg_eff_num_head;
      while(rlp!=NULL) {
        rlp2=rlp;
        rlp=rlp->next;
        free(rlp2);
      }
    }

    /*Initialize reg_counter after all the effectors in the object 
     *had been placed with concern of the region overlapping
    */
    if (reg_count_head!=NULL) {
      if (!chkpt_infile) {
        if (init_region_counter(pop,cdp,reg_count_head)) {
          return(1);
        }
      }
      /* allocate memory for counter_hash_table, and 
       * save the reg_counter_ref_list to the table
       */       
      if (init_counter_hash_table(objp,reg_count_head)) {
	return(1);
      }
      /* free region_list created by the initialization of region counter*/
      rlp=reg_count_head;
      while(rlp!=NULL) {
        rlp3=rlp;
        rlp=rlp->next;
        free(rlp3);
      }
    }

   
    /* free eff_prop array and contents */
    for (i=0;i<n_polys;i++) {
      if (eff_prop[i]!=NULL) {
        effdp=eff_prop[i];
        while(effdp!=NULL) {
          dup_effdp=effdp;
          effdp=effdp->next;
          free(dup_effdp);
        }
      }
    }
    free(eff_prop);

    break;
  case BOX_OBJ:

	/* Allocate memory and initialize */
    if ((cdp=(struct cmprt_data *)malloc
	 (sizeof(struct cmprt_data)))==NULL) {
      return(1);
    }
    if ((cdp->lig_count=(int *)malloc
	 ((1+n_ligand_types)*sizeof(int)))==NULL){
      return(1);
    }
    if ((cdp->conc=(double *)malloc
	 ((1+n_ligand_types)*sizeof(double)))==NULL){
      return(1);
    }
    for (i=0;i<1+n_ligand_types;i++) {
      cdp->lig_count[i]=0;
      cdp->conc[i]=0;
    }
    if ((cdp->corner=(struct vector3 *)malloc
	 (8*sizeof(struct vector3)))==NULL){
      return(1);
    }
    if ((cdp->normal=(struct vector3 *)malloc
	 (n_polys*sizeof(struct vector3)))==NULL){
      return(1);
    }
    if ((cdp->wall=(struct wall **)malloc
	 (n_polys*sizeof(struct wall *)))==NULL){
      return(1);
    }
    if ((cdp->neighbor=(struct cmprt_data **)malloc
	 (n_polys*sizeof(struct cmprt_data *)))==NULL){
      return(1);
    }
    for (i=0;i<n_polys;i++) {
      cdp->wall[i]=NULL;
      cdp->neighbor[i]=NULL;
    }

	/* Setup the metadata */
    cdp->sym=objp->sym;
    cdp->full_name=obj_name;
    cdp->cmprt_type=pop->list_type;
    cdp->fully_closed=pop->fully_closed;
    cdp->vm=0;
    cdp->n_corners=8;
    cdp->n_walls=n_polys;
    cdp->wall_list=NULL;
    cdp->next=cmprt_head;
    cmprt_head=cdp;
    if (vizp!=NULL) {
      if ((cdlp=(struct cmprt_data_list *)malloc
	   (sizeof(struct cmprt_data_list)))==NULL) {
        return(1);
      }
      cdlp->cmprt_data = cdp;
      cdlp->next = vizp->cmprt_data_list;
      vizp->cmprt_data_list = cdlp;
    }
    corner=cdp->corner;

    bpp=(struct box_poly *)pop->polygon_data;
    vp1=bpp->llf;
    vp2=bpp->urb;
    cube_corners(vp1,vp2,corner);

    for (i=0;i<8;i++) {
      p[0][0]=corner[i].x;
      p[0][1]=corner[i].y;
      p[0][2]=corner[i].z;
      p[0][3]=1.0;
      mult_matrix(p,im,p,l,m,n);
      corner[i].x=p[0][0]/length_unit;
      corner[i].y=p[0][1]/length_unit;
      corner[i].z=p[0][2]/length_unit;

	  if ( ((corner[i].x <= volume->x_partitions[0])
             || (corner[i].x >= volume->x_partitions[volume->n_x_subvol]))
        || ((corner[i].y <= volume->y_partitions[0])
             || (corner[i].y >= volume->y_partitions[volume->n_y_subvol]))
        || ((corner[i].z <= volume->z_partitions[0])
             || (corner[i].z >= volume->z_partitions[volume->n_z_subvol])) ){
		throw_away++;
	  }
		
    }

	/* Completely outside this node... throw away */
	if (throw_away == 8) {
	  free(cdp->lig_count);
	  free(cdp->conc);
	  free(cdp->corner);
	  free(cdp->normal);
	  free(cdp->wall);
	  free(cdp->neighbor);
	  cmprt_head=cdp->next;
	  free(cdp);
	  return(0);
	}

    if ((eff_prop=(struct eff_dat **)malloc
       (n_polys*sizeof(struct eff_dat *)))==NULL) {
      return(1);
    }

    vectorize(&corner[0],&corner[4],&vx);
    vectorize(&corner[0],&corner[2],&vy);
    vectorize(&corner[0],&corner[1],&vz);
    cdp->volume=vect_length(&vx)*vect_length(&vy)*vect_length(&vz);
	
    for (i=0;i<6;i++) {
      cube_face(corner,rect_vert,i);
      vectorize(rect_vert[0],rect_vert[1],&v1);
      vectorize(rect_vert[0],rect_vert[3],&v2);
      cross_prod(&v1,&v2,&cdp->normal[i]);
      normalize(&cdp->normal[i]);
      eff_prop[i]=NULL;
    }

    poly_lcrp=pop->lig_count_ref;
    while (poly_lcrp!=NULL) {
      if (strcmp(full_name,poly_lcrp->full_name)==0) {
        poly_lcrp->count_list->temp_data=(void *)&cdp->lig_count[poly_lcrp->type];
      }
      poly_lcrp=poly_lcrp->next;
    }

    /* prepend a copy of eff_dat for each element referenced in each region
       of this object to the eff_prop list for the referenced element */
    reg_eff_num_head=NULL;
    reg_count_head=NULL;
    lig_hit_count=NULL;
    rlp=objp->region_list;
    no_printf("Traversing regions on object %s:\n",objp->sym->name);
    while (rlp!=NULL) {
      rp=rlp->region;
      rcrlp=rp->reg_counter_ref_list;
      no_printf("  checking region %s\n",rp->sym->name);
      reg_eff_num=0;
      reg_count_num=0;
      lig_hit_flag=0;
      if (rcrlp!=NULL) {
	reg_count_num=1;
      }
      lig_hit=rp->lig_hit_counter;
      if (lig_hit!=NULL) {
	lig_hit_flag=1;
      }	
      elp=rp->element_list;
      while (elp!=NULL) {
        for (i=elp->begin;i<=elp->end;i++) {
          if (pop->side_stat[i]) {
            effdp=rp->eff_dat;
            while (effdp!=NULL) {
              if (effdp->quantity_type==EFFDENS) {
                if ((dup_effdp=(struct eff_dat *)malloc
                     (sizeof(struct eff_dat)))==NULL){
                  return(1);
                }
                dup_effdp->rx=effdp->rx;
                dup_effdp->quantity_type=effdp->quantity_type;
                dup_effdp->quantity=effdp->quantity;
                dup_effdp->orient=effdp->orient;
                no_printf("    will add effector state %s to side %d at density = %.9g\n",dup_effdp->rx->sym->name,i,dup_effdp->quantity);
                no_printf("    this effdp = %p\n",(void *)dup_effdp);
                no_printf("    next effdp = %p\n",(void *)eff_prop[i]);
                dup_effdp->next=eff_prop[i];
                eff_prop[i]=dup_effdp;
              }
              else {
                reg_eff_num=1;
                no_printf("    add effector state %s to side %d by number = %.9g\n",effdp->rx->sym->name,i,effdp->quantity);
              }
              effdp=effdp->next;
            }
          }
        }
        elp=elp->next;
      }
      if (reg_eff_num) {
        if ((rlp2=(struct region_list *)malloc
             (sizeof(struct region_list)))==NULL){
          return(1);
        }
        rlp2->region=rp;
        rlp2->next=reg_eff_num_head;
        reg_eff_num_head=rlp2;
      }
      if (reg_count_num) {
        if ((rlp3=(struct region_list *)malloc
             (sizeof(struct region_list)))==NULL){
          return(1);
        }
        rlp3->region=rp;
        rlp3->next=reg_count_head;
        reg_count_head=rlp3;
      }
      if (lig_hit_flag) {
        if ((rlp4=(struct region_list *)malloc
             (sizeof(struct region_list)))==NULL){
          return(1);
        }
        rlp4->region=rp;
        rlp4->next=lig_hit_count;
        lig_hit_count=rlp4;
      }

      rlp=rlp->next;
    }

    for (i=0;i<6;i++) {
      if (pop->side_stat[i]) {
          if ((face=(struct vector3 **)malloc
	       (4*sizeof(struct vector3 *)))==NULL) {
            return(1);
          }
          cube_face(corner,face,i);
	  if ((wp=init_wall(face,NULL,4))==NULL) {
	    return(1);
	  }
	  if ((wlp=(struct wall_list *)malloc
	       (sizeof(struct wall_list)))==NULL) {
	    return(1);
	  }
	  wlp->wall=wp;
	  wlp->next=cdp->wall_list;
	  cdp->wall_list=wlp;
	  cdp->wall[i]=wp;
          pop->cmprt_side_map[i]=i;
	  wp->next_wall=wall_head;
	  wall_head=wp;
          wp->wall_shape=RECT_POLY;

          area=wp->length_first*wp->length_last;
          wp->area=area;
          total_area=total_area+wp->area;
	  
	  wp->parent_object=objp;   
	  wp->side=i;
	  wp->parent=cdp;
	  wp->wall_type=pop->lig_prop[i];
	  wp->effectors=NULL;
          if (objp->viz_state!=NULL) {
	    wp->viz_state=objp->viz_state[i];
          }
          else {
	    wp->viz_state=EXCLUDE_OBJ;
          }

	  /* Try to find out which region contains this wall
	   * Only when there are regions defined in this object 
	   * AND region counters exist for this object.
	   * ELSE we set it NULL. Lin-Wei Wu 6/26/03*/
	  rlp=objp->region_list;
	  if ((rlp!=NULL)&&((reg_count_head!=NULL)||(lig_hit_count!=NULL))) {
	    wp->region_list=init_region_list_for_wall(rlp,i);
	  }
	  else {
	    wp->region_list=NULL;
	  }
	  if ((effdp=eff_prop[i])!=NULL) {
	    if (init_effectors_by_density(wp,effdp)) {
	      return(1);
	    }
	  }
        }
    }
    	  if (lig_hit_count!=NULL) {
	    rlp=lig_hit_count;
	    while (rlp!=NULL) {
	      rlp4=rlp;
	      rlp=rlp->next;
	      free(rlp4);
	    }
	  }

    if (reg_eff_num_head!=NULL) {
      if (init_effectors_by_number(pop,cdp,reg_eff_num_head)) {
        return(1);
      }
      /* free region list created to hold regions populated by number */
      rlp=reg_eff_num_head;
      while(rlp!=NULL) {
        rlp2=rlp;
        rlp=rlp->next;
        free(rlp2);
      }
    }


    /*Initialize reg_counter after all the effectors in the object 
     *had been placed with concern of the region overlapping
    */
    
    if (reg_count_head!=NULL) {
      if (!chkpt_infile) {
        if (init_region_counter(pop,cdp,reg_count_head)) {
          return(1);
        }
      }
      /* allocate memory for counter_hash_table, and 
       * save the reg_counter_ref_list to the table
       */       
      if (init_counter_hash_table(objp,reg_count_head)) {
	return(1);
      }
      /* free region_list created for region counter purpose */
      rlp=reg_count_head;
      while(rlp!=NULL) {
        rlp3=rlp;
        rlp=rlp->next;
        free(rlp3);
      }
    }
    
    
    
    /* free eff_prop array and contents */
    for (i=0;i<n_polys;i++) {
      if (eff_prop[i]!=NULL) {
        effdp=eff_prop[i];
        while(effdp!=NULL) {
          dup_effdp=effdp;
          effdp=effdp->next;
          free(dup_effdp);
        }
      }
    }
    free(eff_prop);

    break;
  }

#ifdef KELP
  cdp->sym->ref_count--;
  if (!cdp->sym->ref_count) {	/* Done with the geometry information */
	destroy_sym_value(cdp->sym);	/* free up memory */
  }
#endif

  return(0);
}



/**
 ** Initialize region counters for the RX_STATE and save 
 ** reg_counter_ref_list to the counter_hash_table of the object
 ** Start from MCell 2.68
 */
int init_region_counter(struct polygon_object *pop, struct cmprt_data *cdp, struct region_list *reg_count_head)
{
  struct rx *rx;
  struct region_list *rlp; 
  struct region *rp;
  struct element_list *elp;
  struct effector *ep;
  struct reg_counter_ref_list *rcrlp;
  struct reg_counter_ref *rcrp;
  int i,j;
  
    /* traverse region list and add effector sites by number to whole regions
       as appropriate */
    rlp=reg_count_head;
    no_printf(" Initialzing region counter ...\n");
    fflush(log_file);
    while (rlp!=NULL) {
      rp=rlp->region;
      elp=rp->element_list;
      /* check effector numbers for each region of the object*/
      while (elp!=NULL) {
	for (i=elp->begin;i<=elp->end;i++) {
	  if (pop->side_stat[i]) {
	    ep=cdp->wall[pop->cmprt_side_map[i]]->effectors;
	    if (ep!=NULL) {
	      for (j=0;j<ep->n_tiles;j++) {
		rx=ep->tiles[j];
		rcrlp=rp->reg_counter_ref_list;
		while (rcrlp!=NULL) {
		  rcrp=rcrlp->reg_counter_ref;
		  if ((rx!=NULL)&&(rx==rcrp->state)&&(rcrp->count_type==RX_STATE)) {
		    rcrp->counter++;
		  }
		  rcrlp=rcrlp->next;
		}
	      }
	    }
	  }
	}
	elp=elp->next;
      }
      rlp=rlp->next;
    }
    no_printf(" Done initializing region counter\n");
    fflush(log_file);
  return(0);
}


/* Initialize the region_list that contains this wall 
 * Add from MCell 2.68, designed for the region counters.
 */
struct region_list  *init_region_list_for_wall(struct region_list *rlp, int index) {
  struct element_list *elp;
  struct region *rp;
  struct region_list *rlp1, *wall_region_head;
  byte wall_on_region_flag;
  
  wall_region_head=NULL; 
 
  while (rlp!=NULL) {
    wall_on_region_flag=0; 
    rp=rlp->region;
    elp=rp->element_list;
    while (elp!=NULL) {
      if ((index<=elp->end)&&(index>=elp->begin)) {
	wall_on_region_flag=1;
      }
      elp=elp->next;
    }

    if (wall_on_region_flag) {
      if ((rlp1=(struct region_list *)malloc(sizeof(struct region_list)))==NULL) {
	mdlerror("Can not save region_list for wall");
	return(NULL);
      }
      rlp1->region=rp;
      rlp1->next=wall_region_head;
      wall_region_head=rlp1;
    }
    rlp=rlp->next;
  }
  return(wall_region_head);

  /* free memory not needed any more*/  
  if (wall_region_head!=NULL) {
    rlp=wall_region_head;
    while (rlp!=NULL) {
      rlp1=rlp;
      rlp=rlp->next;
      free(rlp1);
    }
  }
}


/*
   ** Initialize the counter_hash_table of the object if 
   ** there are any reg_counters defined in this object, 
   ** and save the reg_counter_ref_list to the table.
   ** The table will be used for future counter searching.
   ** Start from MCell 2.68
 */  

int init_counter_hash_table (struct object *objp, struct region_list *reg_count_head) {

  struct counter_hash_table **countertab;
  struct counter_hash_table *chtp,*prev,*hash_table_head;
  struct region_list *rlp;
  struct region *rp;
  struct reg_counter_ref_list *rcrlp, *rcrlp2, *rcrlp_tmp;
  struct reg_counter_ref *rcrp, *rcrp2;
  char *rx_name;
  int col, done, i;
  unsigned short hashval;

  rlp=reg_count_head;
  rcrlp2=NULL;
  rcrp2=NULL;
  hash_table_head=NULL;
  no_printf("\t Initializing counter hash table for object %s\n",objp->sym->name);
  fflush(log_file);
  /* Initialize the object counter hash table if it is empty*/
  if (objp->counter_hash_table==NULL) {
    if ((countertab=(struct counter_hash_table **)malloc(COUNTER_HASH*sizeof(struct counter_hash_table *)))==NULL) {
      mdlerror("Can not save counter table to object");
      return(1);
    }
    countertab=init_countertab(COUNTER_HASH);
    objp->counter_hash_table=countertab;
  }
  
  /* save each reg_counter_ref_list to the table by the hash value of 
   * the rx state of the counter and sort by the address of the region
   */  
  while (rlp!=NULL) {
    rp=rlp->region;
    rcrlp=rp->reg_counter_ref_list;
    while (rcrlp!=NULL) {
      rcrp=rcrlp->reg_counter_ref;
      rx_name=rcrp->state->sym->name;
      hashval=hash(rx_name)&0x0000000f;
      /* Store counter in table, but check if it is already saved */
      if ((chtp=retrieve_counter(rx_name,rcrp, objp->counter_hash_table))==NULL) {
       if ((chtp=store_counter(rx_name,rcrlp, objp->counter_hash_table))==NULL) {
	fprintf(log_file,"Cannot store counter in table: %s\n",rx_name);
	return(NULL);
      }
      }
      else {
	warning(" Counter already defined");
	return(NULL);
      }
      objp->counter_hash_table[hashval]=chtp;      
      rcrlp=rcrlp->next;
    }
  rlp=rlp->next;
  }
  /*Print the counter_hash_table for debuging*/
  /*  
  for (i=0; i<COUNTER_HASH; i++) {
    hash_table_head=objp->counter_hash_table[i];
    while (hash_table_head!=NULL) {
      rcrlp2=(struct reg_counter_ref_list *)hash_table_head->value;
      if (rcrlp2!=NULL) {
	rcrp2=rcrlp2->reg_counter_ref;
	no_printf("\t\t col %d, counter %s, region %s = %d, type %d, counter %d;\n",i,rcrp2->state->sym->name,rcrp2->parent->sym->name,(int)rcrp2->parent, rcrp2->count_type,rcrp2->counter);
	fflush(log_file);
       }
      hash_table_head=hash_table_head->next;
     }
  }
  */  	
  no_printf("\t end of counter hash table initializing\n");
  fflush(log_file); 
  return(0);
}


/**
 * Updates the bounding box of the world based on the size
 * and location of a polygon_object (box_poly or polygon_object).
 * Used by compute_bb().
 */
int compute_bb_polygon_object(struct object *objp, double (*im)[4], char *full_name)
{
  struct polygon_object *pop;
  struct ordered_poly *opp;
  struct box_poly *bpp;
  struct vector3 corner[8],*vp1,*vp2;
  double p[1][4];
  int i,j,cnt,n_verts,n_polys,n_element_verts,index_0,index_1,index_n;
  unsigned short l,m,n;

  pop=(struct polygon_object *)objp->contents;
  n_polys=pop->n_polys;
  l=1;
  m=4;
  n=4;

  switch (objp->object_type) {
  case POLY_OBJ:
    if (pop->list_type==ORDERED_POLY) {
      opp=(struct ordered_poly *)pop->polygon_data;
      n_verts=opp->n_verts;
    }
 
    for (i=0;i<n_verts;i++) {
      if (pop->list_type==ORDERED_POLY) {
        p[0][0]=opp->vertex[i]->x;
        p[0][1]=opp->vertex[i]->y;
        p[0][2]=opp->vertex[i]->z;
        p[0][3]=1.0;
      }
      mult_matrix(p,im,p,l,m,n);
      if (p[0][0]<bb_min.x) {
        bb_min.x=p[0][0];
      }
      if (p[0][1]<bb_min.y) {
        bb_min.y=p[0][1];
      }
      if (p[0][2]<bb_min.z) {
        bb_min.z=p[0][2];
      }
      if (p[0][0]>bb_max.x) {
        bb_max.x=p[0][0];
      }
      if (p[0][1]>bb_max.y) {
        bb_max.y=p[0][1];
      }
      if (p[0][2]>bb_max.z) {
        bb_max.z=p[0][2];
      }
    }

/*
    fprintf(log_file,"MCell: Total area of physical object %s = %.9g microns^2\n",objp->sym->name,total_area*length_unit*length_unit);
*/
    break;
  case BOX_OBJ:
    bpp=(struct box_poly *)pop->polygon_data;
    vp1=bpp->llf;
    vp2=bpp->urb;
    cube_corners(vp1,vp2,corner);

    for (i=0;i<8;i++) {
      p[0][0]=corner[i].x;
      p[0][1]=corner[i].y;
      p[0][2]=corner[i].z;
      p[0][3]=1.0;
      mult_matrix(p,im,p,l,m,n);
      if (p[0][0]<bb_min.x) {
        bb_min.x=p[0][0];
      }
      if (p[0][1]<bb_min.y) {
        bb_min.y=p[0][1];
      }
      if (p[0][2]<bb_min.z) {
        bb_min.z=p[0][2];
      }
      if (p[0][0]>bb_max.x) {
        bb_max.x=p[0][0];
      }
      if (p[0][1]>bb_max.y) {
        bb_max.y=p[0][1];
      }
      if (p[0][2]>bb_max.z) {
        bb_max.z=p[0][2];
      }
    }

/*
    fprintf(log_file,"MCell: Total area of physical object %s = %.9g microns^2\n",objp->sym->name,total_area*length_unit*length_unit);
*/
    break;
  }

  return(0);
}

/**
 * Computes the corners of the cuboid whose bounding box is
 * specified by the two points p1 and p2.
 * @param p1 a vector3 specifying one corner of the cuboid.
 * @param p2 a vector3 specifying the other diagonal corner of the cuboid.
 * @param corner an array of 8 vector3s which consitute the 8 corners
 *               of the cuboid.
 */
void cube_corners(struct vector3 *p1, struct vector3 *p2, struct vector3 *corner)
{
  double dx,dy,dz;

  dx=p2->x-p1->x;
  dy=p2->y-p1->y;
  dz=p2->z-p1->z;
  if (dx<0) {
    swap_double(&p1->x,&p2->x);
  }
  if (dy<0) {
    swap_double(&p1->y,&p2->y);
  }
  if (dz<0) {
    swap_double(&p1->z,&p2->z);
  }
  corner[0].x=p1->x;
  corner[0].y=p1->y;
  corner[0].z=p1->z;
  corner[1].x=p1->x;
  corner[1].y=p1->y;
  corner[1].z=p2->z;
  corner[2].x=p1->x;
  corner[2].y=p2->y;
  corner[2].z=p1->z;
  corner[3].x=p1->x;
  corner[3].y=p2->y;
  corner[3].z=p2->z;
  corner[4].x=p2->x;
  corner[4].y=p1->y;
  corner[4].z=p1->z;
  corner[5].x=p2->x;
  corner[5].y=p1->y;
  corner[5].z=p2->z;
  corner[6].x=p2->x;
  corner[6].y=p2->y;
  corner[6].z=p1->z;
  corner[7].x=p2->x;
  corner[7].y=p2->y;
  corner[7].z=p2->z;
}


/**
 * Sets up the corners of one face of a cuboid.
 * Assigns vector3's corresponding to the corners of a cuboid
 * from 'corner' to the corresponding faces of the cuboid.
 * <br>
 * Builds the faces using the right hand rule.
 * @param corner an array of 8 vector3s which consitute the 8 corners
 *               of the cuboid.
 * @param face an array 4 vector3's that get assigned pointers to the
 *             corresponding 4 corner vector3s
 * @param i an integer that controls which face is computed.
 *          i == TP (top)
 *          i == BOT (bottom)
 *          i == FRNT (front)
 *          i == BCK (back)
 *          i == LFT (left)
 *          i == RT (right)
 */
void cube_face(struct vector3 *corner, struct vector3 **face, int i)
{
      /* Build each face using right-hand rule */
      if (i==TP) {
        /* top face */
        face[0]=&corner[1];
        face[1]=&corner[5];
        face[2]=&corner[7];
        face[3]=&corner[3];
      }
      else if (i==BOT) {
        /* bottom face */
        face[0]=&corner[0];
        face[1]=&corner[2];
        face[2]=&corner[6];
        face[3]=&corner[4];
      }
      else if (i==FRNT) {
        /* front face */
        face[0]=&corner[0];
        face[1]=&corner[4];
        face[2]=&corner[5];
        face[3]=&corner[1];
      }
      else if (i==BCK) {
        /* back face */
        face[0]=&corner[2];
        face[1]=&corner[3];
        face[2]=&corner[7];
        face[3]=&corner[6];
      }
      else if (i==LFT) {
        /* left face */
        face[0]=&corner[0];
        face[1]=&corner[1];
        face[2]=&corner[3];
        face[3]=&corner[2];
      }
      else if (i==RT) {
        /* right face */
        face[0]=&corner[4];
        face[1]=&corner[6];
        face[2]=&corner[7];
        face[3]=&corner[5];
      }
      return;
}


/**
 * Sets up the corners of all the faces of a cuboid.
 * Assigns vector3's corresponding to the corners of a cuboid
 * from 'corner' to the corresponding faces of the cuboid.
 * <br>
 * Builds the faces using the right hand rule.
 * <br>
 * Similar to cube_face(). [Note difference in spelling]
 */
void cube_faces(struct vector3 *corner, struct vector3 *(*face)[4])
{

  /* Build each face using right-hand rule */
  /* top face */
  face[0][0]=&corner[1];
  face[0][1]=&corner[5];
  face[0][2]=&corner[7];
  face[0][3]=&corner[3];
 
  /* bottom face */
  face[1][0]=&corner[0];
  face[1][1]=&corner[2];
  face[1][2]=&corner[6];
  face[1][3]=&corner[4];
 
  /* front face */
  face[2][0]=&corner[0];
  face[2][1]=&corner[4];
  face[2][2]=&corner[5];
  face[2][3]=&corner[1];
 
  /* back face */
  face[3][0]=&corner[2];
  face[3][1]=&corner[3];
  face[3][2]=&corner[7];
  face[3][3]=&corner[6];
 
  /* left face */
  face[4][0]=&corner[0];
  face[4][1]=&corner[1];
  face[4][2]=&corner[3];
  face[4][3]=&corner[2];
 
  /* right face */
  face[5][0]=&corner[4];
  face[5][1]=&corner[6];
  face[5][2]=&corner[7];
  face[5][3]=&corner[5];

  return;
}

int init_effector_grid(struct wall *wp)
{
  struct effector *ep;
  struct rx **tiles;
  struct vector3 v1,u_axis,v_axis,p_b,p_c,p_d,i_axis,j_axis;
  struct vector3 ab,bc,ac,step_u,step_v,diagonal,p0,p1,p2;
  unsigned short *dsp,*psp;
  signed char *orp;
  int *tsp;
  unsigned int i,j,k,l,m,n,p,nr,nl,ir,jr,n_ligs,n_occupied;
  int p_index,rx_index,subvol;
  int grid_size,uu,vv,vv_max;
  byte grid_shape;
  double area;
  double diag_x,diag_y,r_slope,width,u_width;
  double u_factor,u_factor_2,v_factor,v_val,binding_factor;

  no_printf("Creating new effector grid...\n");

  vectorize(wp->vert[0],wp->vert[1],&ab);
  vectorize(wp->vert[1],wp->vert[2],&bc);
  if (wp->wall_shape==RECT_POLY) {
    grid_shape=RECTANGULAR;
    area=wp->area;
    no_printf("Area = %g\n",area);
    grid_size=(int) ceil(sqrt(area/2));
    if (grid_size==0 || grid_size>46340) {
      fprintf(log_file,"MCell: too many tiles in effector grid: %d\n",grid_size);
      fflush(log_file);
      return(1);
    }
    n=2*grid_size*grid_size;
    binding_factor=n/area;
    no_printf("binding_factor = %g\n",binding_factor);
    r_slope=(wp->length_first/wp->length_last);
    diag_x=wp->length_first;
    diag_y=wp->length_last;
    u_factor=grid_size/wp->length_first;
    u_factor_2=wp->length_first/grid_size;
    v_factor=grid_size/wp->length_last;
    u_axis.x=1;
    u_axis.y=0;
    u_axis.z=0;
    v_axis.x=0;
    v_axis.y=1;
    v_axis.z=0;
    vectorize(wp->vert[0],wp->vert[1],&i_axis);
    normalize(&i_axis);
    vectorize(wp->vert[0],wp->vert[3],&j_axis);
    normalize(&j_axis);
  }
  else {
    grid_shape=TRIANGULAR;
    vectorize(wp->vert[0],wp->vert[2],&ac);
    area=wp->area;
    no_printf("Area = %g\n",area);
    grid_size=(int) ceil(sqrt(area));
    if (grid_size==0 || grid_size>65536) {
      fprintf(log_file,"MCell: too many tiles in effector grid: %d\n",grid_size);
      fflush(log_file);
      return(1);
    }
    n=grid_size*grid_size;
    binding_factor=n/area;
    no_printf("binding_factor = %g\n",binding_factor);
    p_b.x=wp->length_first;
    p_b.y=0;
    p_b.z=0;
    width=dot_prod(&ab,&ac)/wp->length_first;
    p_c.x=width;
    p_c.y=2*area/wp->length_first;
    p_c.z=0;
    r_slope=(p_c.x/p_c.y);
    diag_x=p_c.x;
    diag_y=p_c.y;
    p_d.x=wp->vert[0]->x+(width*ab.x/wp->length_first);
    p_d.y=wp->vert[0]->y+(width*ab.y/wp->length_first);
    p_d.z=wp->vert[0]->z+(width*ab.z/wp->length_first);
    vectorize(wp->vert[0],wp->vert[1],&i_axis);
    normalize(&i_axis);
    vectorize(&p_d,wp->vert[2],&j_axis);
    normalize(&j_axis);
    u_axis.x=p_c.y-p_b.y;
    u_axis.y=p_b.x-p_c.x;
    u_axis.z=0;
    normalize(&u_axis);
    u_width=dot_prod(&p_b,&u_axis);
    u_factor=grid_size/u_width;
    u_factor_2=wp->length_first/grid_size;
    v_factor=grid_size/p_c.y;
    vectorize(&p_b,&p_c,&v_axis);
    v_axis.x=v_axis.x/grid_size;
    v_axis.y=v_axis.y/grid_size;
    v_axis.z=0;

  }
  step_u.x=ab.x/grid_size;
  step_u.y=ab.y/grid_size;
  step_u.z=ab.z/grid_size;
  step_v.x=bc.x/grid_size;
  step_v.y=bc.y/grid_size;
  step_v.z=bc.z/grid_size;
  vectorize(&step_u,&step_v,&diagonal);

  no_printf("Initializing effector grid size %d ...\n",n);
  fflush(log_file);

  if ((tiles=(struct rx **)malloc(n*sizeof(struct rx *)))==NULL) {
    fprintf(log_file,"MCell: cannot store tiles for effector grid: %d\n",n);
    fflush(log_file);
    return(1);
  }
  if ((tsp=(int *)malloc(n*sizeof(int)))==NULL) {
    fprintf(log_file,"MCell: cannot store time stamp data for effector grid: %d\n",n);
    fflush(log_file);
    return(1);
  }
  if ((dsp=(unsigned short *)malloc(n*sizeof(unsigned short)))==NULL) {
    fprintf(log_file,"MCell: cannot store desired state data for effector grid: %d\n",n);
    fflush(log_file);
    return(1);
  }
  if ((orp=(signed char *)malloc(n*sizeof(signed char)))==NULL) {
    fprintf(log_file,"MCell: cannot store orientation data\n");
    fflush(log_file);
    return(1);
  }

  for (i=0;i<n;i++) {
    tiles[i]=NULL;
    tsp[i]=(-INT_MAX);
    dsp[i]=0;
    orp[i]=0;
  }

  no_printf("Done initializing effector grid size %d\n",n);
  fflush(log_file);

  if ((ep=(struct effector *)malloc(sizeof(struct effector)))==NULL) {
    return(1);
  }
  wp->effectors=ep;
  wp->effectors->grid_shape=grid_shape;
  wp->effectors->grid_size=grid_size;
  wp->effectors->u_axis.x=u_axis.x;
  wp->effectors->u_axis.y=u_axis.y;
  wp->effectors->v_axis.x=v_axis.x;
  wp->effectors->v_axis.y=v_axis.y;
  wp->effectors->i_axis.x=i_axis.x;
  wp->effectors->i_axis.y=i_axis.y;
  wp->effectors->i_axis.z=i_axis.z;
  wp->effectors->j_axis.x=j_axis.x;
  wp->effectors->j_axis.y=j_axis.y;
  wp->effectors->j_axis.z=j_axis.z;
  wp->effectors->r_slope=r_slope;
  wp->effectors->diag.x=diag_x;
  wp->effectors->diag.y=diag_y;
  wp->effectors->u_factor=u_factor;
  wp->effectors->u_factor_2=u_factor_2;
  wp->effectors->v_factor=v_factor;
  wp->effectors->r_u_factor=1/u_factor;
  wp->effectors->r_v_factor=1/v_factor;
  wp->effectors->binding_factor=binding_factor;
  wp->effectors->n_tiles=n;
  wp->effectors->n_occupied=0;
  
  wp->effectors->i=0; /* deprecated */
  wp->effectors->j=0; /* deprecated */

  wp->effectors->tiles=tiles;
  wp->effectors->time_stamp=tsp;
  wp->effectors->desired_state=dsp;
  wp->effectors->prev_state=NULL;
  wp->effectors->index=n_effector_grids++;
  wp->effectors->n_types=0;
  wp->effectors->orient=orp;
  wp->effectors->set=0;
  wp->effectors->wall=wp;
	
  no_printf("Done creating new effector grid.\n");

  return(0);
}


int init_effectors_by_number(struct polygon_object *pop, struct cmprt_data *cdp, struct region_list *reg_eff_num_head)
{
  struct rx *rx,***tiles,***tiles_tmp;
  struct region_list *rlp; 
  struct region *rp;
  struct element_list *elp;
  struct effector *ep;
  struct eff_dat *effdp;
  struct wall **walls,**walls_tmp,*wp;
  struct ligand *lp;
  struct vector3 ab,bc,step_u,step_v,diagonal,p0,p1,p2;
  unsigned short **dsp,**dsp_tmp,desired_state;
  unsigned short *psp;
  byte prev_state_flag,grid_shape;
  signed char **orient,**orient_tmp,orientation;
  double rand[1];
  double fuzz,tiny_x,tiny_y,tiny_z,r1,r2,r3,v_val;
  unsigned int *index,*index_tmp;
  unsigned int n_free_eff,n_set,n_clear;
  unsigned int i,j,k,m,p,nl,n_ligs;
  int grid_size;
  int subvol,uu,vv;
  byte done;

    no_printf("Initializing effectors by number...\n");

/* Place initially bound ligands on appropriate pole of effector */
/*
    fuzz=POLE*0.000001;
*/
    fuzz=0.0;

    nl=1+n_ligand_types;
    n_ligs=0;

    tiles=NULL;
    tiles_tmp=NULL;
    index=NULL;
    index_tmp=NULL;
    orient=NULL;
    orient_tmp=NULL;
    dsp=NULL;
    dsp_tmp=NULL;
    walls=NULL;
    walls_tmp=NULL;
    /* traverse region list and add effector sites by number to whole regions
       as appropriate */
    rlp=reg_eff_num_head;
    while (rlp!=NULL) {
      rp=rlp->region;
        /* initialize effector grids in region as needed and */
        /* count total number of free effector sites in region */
        n_free_eff=0;
        elp=rp->element_list;
        while (elp!=NULL) {
          for (i=elp->begin;i<=elp->end;i++) {
            if (pop->side_stat[i]) {
              wp=cdp->wall[pop->cmprt_side_map[i]];
              ep=wp->effectors;
              if (ep==NULL) {
                if (init_effector_grid(wp)) {
                  return(1);
                }
                ep=wp->effectors;
              }
              n_free_eff=n_free_eff+(ep->n_tiles-ep->n_occupied);
            }
          }
          elp=elp->next;
        }
        no_printf("Number of free effector tiles in region %s = %d\n",rp->sym->name,n_free_eff);
        fflush(stdout);
      if (chkpt_init) {  /* only needed for denovo initiliazation */
        /* allocate memory to hold array of pointers to all free tiles */
        if ((tiles=(struct rx ***)malloc
           (n_free_eff*sizeof(struct rx **)))==NULL) {
          return(1);
        }
        if ((index=(unsigned int *)malloc
           (n_free_eff*sizeof(unsigned int)))==NULL) {
          return(1);
        }
        if ((orient=(signed char **)malloc
           (n_free_eff*sizeof(signed char *)))==NULL) {
          return(1);
        }
        if ((dsp=(unsigned short **)malloc(n_free_eff*sizeof(unsigned short *)))==NULL) {
          return(1);
        }
        if ((walls=(struct wall **)malloc
           (n_free_eff*sizeof(struct wall *)))==NULL) {
          return(1);
        }
        /* initialize array of pointers to all free tiles */
        k=0;
        elp=rp->element_list;
        while (elp!=NULL) {
          for (i=elp->begin;i<=elp->end;i++) {
            if (pop->side_stat[i]) {
              ep=cdp->wall[pop->cmprt_side_map[i]]->effectors;
              if (ep!=NULL) {
                for (j=0;j<ep->n_tiles;j++) {
                  if (ep->tiles[j]==NULL) {
                    tiles[k]=&(ep->tiles[j]);
                    index[k]=j;
                    orient[k]=&(ep->orient[j]);
                    dsp[k]=&(ep->desired_state[j]);
                    walls[k++]=ep->wall;
                  }
                }
              }
            }
          }
          elp=elp->next;
        }
      } /* chkpt_init */

      /* distribute desired number of effector sites */
      /* for each effector type to add */
      effdp=rp->eff_dat;
      prev_state_flag=0;
      while (effdp!=NULL) {
        if (effdp->quantity_type==EFFNUM) {
          rx=effdp->rx;
          prev_state_flag=prev_state_flag||rx->parent_rx->prev_state_flag;

          if (chkpt_init) {  /* only needed for denovo initiliazation */
	    orientation=effdp->orient;
            desired_state=rx->state_index;
  
            n_set=effdp->quantity;
            n_clear=n_free_eff-n_set;
            rx->count+=n_set;

            if (n_set > n_free_eff) {
              fprintf(log_file,"\nMCell: Warning -- Number of %s effectors to place (%d) exceeds number of free effector tiles (%d) in region %s[%s]\n\n",rx->sym->name,n_set,n_free_eff,rp->parent->sym->name,rp->region_last_name);
              n_set=n_free_eff;
              n_clear=0;
            }
            no_printf("distribute %d of effector %s\n",n_set,rx->sym->name);
            no_printf("n_set = %d  n_clear = %d  n_free_eff = %d\n",n_set,n_clear,n_free_eff);
            fflush(stdout);

            /* if filling more than half the free tiles then init all to full */
            /* and choose which tiles to free again */
            if (n_set > n_free_eff/2) {
              no_printf("filling more than half the free tiles: init all to full\n");
              fflush(stdout);
              for (j=0;j<n_free_eff;j++) {
                *tiles[j]=rx;
                *orient[j]=orientation;
                *dsp[j]=desired_state;
              }
              no_printf("choose which tiles to free again\n");
              fflush(stdout);
              for (j=0;j<n_clear;j++) {
                done=0;
                while (!done) {
                  random_number_use++;
                  ran4(&seed,rand,1,n_free_eff);
                  k=rand[0];
                  if (*tiles[k]!=NULL) {
                    *tiles[k]=NULL;
                    *orient[k]=0;
                    *dsp[k]=0;
                    done=1;
                  }
                }
              }
              /* create bound molecules as needed on newly filled tiles */
              for (j=0;j<n_free_eff;j++) {
                if (*tiles[j]!=NULL) {
                  rx=*tiles[j];
                    for (m=1;m<nl;m++) {
	              for (p=0;p<rx->bound_ligands[m];p++) {
	                if ((lp=(struct ligand *)malloc
	                   (sizeof(struct ligand)))==NULL) {
	                  fprintf(log_file,"MCell: cannot store molecule\n");
	                  return(1);
	                }
	                lp->next_ligand=ligand_table[m]->top;
	                ligand_table[m]->top=lp;
	                tot_mols++;
	                ligand_table[m]->n_mols++;
                        n_ligs++;
	                lp->lig_num=ligand_table[m]->lig_index++;
	                
                        wp=walls[j];
                        tiny_x=fuzz*wp->normal.x;
                        tiny_y=fuzz*wp->normal.y;
                        tiny_z=fuzz*wp->normal.z;
                        r1=wp->vert[0]->x;
                        r2=wp->vert[0]->y;
                        r3=wp->vert[0]->z;
  
                        vectorize(wp->vert[0],wp->vert[1],&ab);
                        vectorize(wp->vert[1],wp->vert[2],&bc);
                        grid_size=wp->effectors->grid_size;
                        grid_shape=wp->effectors->grid_shape;
                        step_u.x=ab.x/grid_size;
                        step_u.y=ab.y/grid_size;
                        step_u.z=ab.z/grid_size;
                        step_v.x=bc.x/grid_size;
                        step_v.y=bc.y/grid_size;
                        step_v.z=bc.z/grid_size;
                        vectorize(&step_u,&step_v,&diagonal);
                        if (grid_shape==TRIANGULAR) {
                          uu=(int) sqrt((double) index[j]);
                          vv=index[j]-(uu*uu);
                        }
                        else {
                          uu=index[j]/(2*grid_size);
                          vv=index[j]-(uu*2*grid_size);
                        }
                        v_val=vv/2;
                        p0.x=(uu+1)*step_u.x;
                        p0.y=(uu+1)*step_u.y;
                        p0.z=(uu+1)*step_u.z;
                        p1.x=p0.x+(v_val*step_v.x);
                        p1.y=p0.y+(v_val*step_v.y);
                        p1.z=p0.z+(v_val*step_v.z);
                        if (vv%2==0) {
                          p2.x=p1.x+(0.33*diagonal.x);
                          p2.y=p1.y+(0.33*diagonal.y);
                          p2.z=p1.z+(0.33*diagonal.z);
                        }
                        else {
                          p2.x=p1.x+(0.66*diagonal.x);
                          p2.y=p1.y+(0.66*diagonal.y);
                          p2.z=p1.z+(0.66*diagonal.z);
                        }
  
	                /* put ligand in middle of effector */
	                lp->pos.x=r1+p2.x+tiny_x;
	                lp->pos.y=r2+p2.y+tiny_y;
	                lp->pos.z=r3+p2.z+tiny_z;
	                
	                lp->effector=wp->effectors;
	                lp->index=index[j];
                        subvol=find_subvol(volume,&lp->pos,-1);
                        if (subvol<0) {
	                  fprintf(log_file,"MCell: invalid subvolume referenced while initializing effectors\n");
                          return(1);
                        }
            
	                lp->subvol=subvol;
	              }
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
                  random_number_use++;
                  ran4(&seed,rand,1,n_free_eff);
                  k=rand[0];
                  if (*tiles[k]==NULL) {
                    *tiles[k]=rx;
                    *orient[k]=orientation;
                    *dsp[k]=desired_state;

                    /* create bound molecules as needed */
                    for (m=1;m<nl;m++) {
	              for (p=0;p<rx->bound_ligands[m];p++) {
	                if ((lp=(struct ligand *)malloc
	                   (sizeof(struct ligand)))==NULL) {
	                  fprintf(log_file,"MCell: cannot store molecule\n");
	                  return(1);
	                }
	                lp->next_ligand=ligand_table[m]->top;
	                ligand_table[m]->top=lp;
	                tot_mols++;
	                ligand_table[m]->n_mols++;
                        n_ligs++;
	                lp->lig_num=ligand_table[m]->lig_index++;
	                
                        wp=walls[k];
                        tiny_x=fuzz*wp->normal.x;
                        tiny_y=fuzz*wp->normal.y;
                        tiny_z=fuzz*wp->normal.z;
                        r1=wp->vert[0]->x;
                        r2=wp->vert[0]->y;
                        r3=wp->vert[0]->z;
  
                        vectorize(wp->vert[0],wp->vert[1],&ab);
                        vectorize(wp->vert[1],wp->vert[2],&bc);
                        grid_size=wp->effectors->grid_size;
                        grid_shape=wp->effectors->grid_shape;
                        step_u.x=ab.x/grid_size;
                        step_u.y=ab.y/grid_size;
                        step_u.z=ab.z/grid_size;
                        step_v.x=bc.x/grid_size;
                        step_v.y=bc.y/grid_size;
                        step_v.z=bc.z/grid_size;
                        vectorize(&step_u,&step_v,&diagonal);
                        if (grid_shape==TRIANGULAR) {
                          uu=(int) sqrt((double) index[k]);
                          vv=index[k]-(uu*uu);
                        }
                        else {
                          uu=index[k]/(2*grid_size);
                          vv=index[k]-(uu*2*grid_size);
                        }
                        v_val=vv/2;
                        p0.x=(uu+1)*step_u.x;
                        p0.y=(uu+1)*step_u.y;
                        p0.z=(uu+1)*step_u.z;
                        p1.x=p0.x+(v_val*step_v.x);
                        p1.y=p0.y+(v_val*step_v.y);
                        p1.z=p0.z+(v_val*step_v.z);
                        if (vv%2==0) {
                          p2.x=p1.x+(0.33*diagonal.x);
                          p2.y=p1.y+(0.33*diagonal.y);
                          p2.z=p1.z+(0.33*diagonal.z);
                        }
                        else {
                          p2.x=p1.x+(0.66*diagonal.x);
                          p2.y=p1.y+(0.66*diagonal.y);
                          p2.z=p1.z+(0.66*diagonal.z);
                        }
  
	                /* put ligand in middle of effector */
	                lp->pos.x=r1+p2.x+tiny_x;
	                lp->pos.y=r2+p2.y+tiny_y;
	                lp->pos.z=r3+p2.z+tiny_z;
	                
	                lp->effector=wp->effectors;
	                lp->index=index[k];
                        subvol=find_subvol(volume,&lp->pos,-1);
                        if (subvol<0) {
	                  fprintf(log_file,"MCell: invalid subvolume referenced while initializing effectors\n");
                          return(1);
                        }
            
	                lp->subvol=subvol;
	              }
                    }

                    done=1;
                  }
                }
              }
            }
        /* allocate memory to hold array of pointers to remaining free tiles */
            if ((tiles_tmp=(struct rx ***)malloc
                 (n_clear*sizeof(struct rx **)))==NULL) {
              return(1);
            }
            if ((index_tmp=(unsigned int *)malloc
               (n_clear*sizeof(unsigned int)))==NULL) {
              return(1);
            }
            if ((orient_tmp=(signed char **)malloc
                 (n_clear*sizeof(signed char *)))==NULL) {
              return(1);
            }
            if ((dsp_tmp=(unsigned short **)malloc
                 (n_clear*sizeof(unsigned short *)))==NULL) {
              return(1);
            }
            if ((walls_tmp=(struct wall **)malloc
               (n_clear*sizeof(struct wall *)))==NULL) {
              return(1);
            }
            k=0;
            for (i=0;i<n_free_eff;i++) {
              if (*tiles[i]==NULL) {
                tiles_tmp[k]=tiles[i];
                index_tmp[k]=index[i];
                orient_tmp[k]=orient[i];
                dsp_tmp[k]=dsp[i];
                walls_tmp[k++]=walls[i];
              }
            }
            /* free original array of pointers to all free tiles */
            free(tiles);
            free(index);
            free(orient);
            free(dsp);
            free(walls);
            tiles=tiles_tmp;
            index=index_tmp;
            orient=orient_tmp;
            dsp=dsp_tmp;
            walls=walls_tmp;
            n_free_eff=n_free_eff-n_set;

            /* update n_occupied for each effector grid */
            elp=rp->element_list;
            while (elp!=NULL) {
              for (i=elp->begin;i<=elp->end;i++) {
                if (pop->side_stat[i]) {
                  ep=cdp->wall[pop->cmprt_side_map[i]]->effectors;
                  if (ep!=NULL) {
		    ep->n_occupied=0;
                    for (j=0;j<ep->n_tiles;j++) {
                      if (ep->tiles[j]!=NULL) {
                        ep->n_occupied++;
                      }
                    }
                  }
                }
              }
              elp=elp->next;
            }
          } /* chkpt_init */
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
      if (orient!=NULL) {
        free(orient);
      }
      if (dsp!=NULL) {
        free(dsp);
      }
      if (walls!=NULL) {
        free(walls);
      }
      /* create arrays to hold prev_state as necessary */
      if (prev_state_flag) {
        elp=rp->element_list;
        while (elp!=NULL) {
          for (i=elp->begin;i<=elp->end;i++) {
            if (pop->side_stat[i]) {
              ep=cdp->wall[pop->cmprt_side_map[i]]->effectors;
              if (ep!=NULL) {
		psp=ep->prev_state;
                if (psp==NULL) {
                  if ((psp=(unsigned short *)malloc(ep->n_tiles*sizeof(unsigned short)))==NULL) {
                    return(1);
                  }
                  for (j=0;j<ep->n_tiles;j++) {
                    psp[j]=0;
                  }
                  ep->prev_state=psp;
                }
              }
            }
          }
          elp=elp->next;
        }
      }
      rlp=rlp->next;
    }
  no_printf("Done initialize effectors by number.\n");
  return(0);
}


int init_effectors_by_density(struct wall *wp, struct eff_dat *effdp_head)
{
  struct rx **tiles,*rx[NUM_ADD_EFFECTORS];
  struct effector *ep;
  struct eff_dat *effdp;
  struct vector3 v1,u_axis,v_axis,p_b,p_c,p_d,i_axis,j_axis;
  struct vector3 ab,bc,ac,step_u,step_v,diagonal,p0,p1,p2;
  unsigned short *dsp,*psp;
  signed char *orp,orientation[NUM_ADD_EFFECTORS];
  int *tsp;
  unsigned int i,j,k,l,m,n,p,nr,nl,ir,jr,n_ligs,n_occupied;
  unsigned int n_ran;
  int p_index,rx_index,subvol;
  int grid_size,uu,vv,vv_max;
  byte grid_shape;
  byte prev_state_flag;
  double rand[1],prob[NUM_ADD_EFFECTORS],area,tot_prob,tot_density;
  double fuzz,tiny_x,tiny_y,tiny_z,length0,length3,r1,r2,r3,i1,i2,i3,j1,j2,j3;
  double diag_x,diag_y,r_slope,width,u_width;
  double u_factor,u_factor_2,v_factor,v_val,binding_factor;
  struct ligand *lp;

  no_printf("Initializing effectors by density...\n");

  ep=wp->effectors;
  if (ep==NULL) {
    if (init_effector_grid(wp)) {
      return(1);
    }
  }
  ep=wp->effectors;

/* Place initially bound ligands on appropriate pole of effector */
/*
  fuzz=POLE*0.000001;
*/
  fuzz=0.0;
  tiny_x=fuzz*wp->normal.x;
  tiny_y=fuzz*wp->normal.y;
  tiny_z=fuzz*wp->normal.z;

  n=ep->n_tiles;
  no_printf("Initializing %d effectors...\n",n);
  fflush(log_file);
  area=wp->area;
  vectorize(wp->vert[0],wp->vert[1],&ab);
  vectorize(wp->vert[1],wp->vert[2],&bc);
  grid_size=ep->grid_size;
  grid_shape=ep->grid_shape;
  step_u.x=ab.x/grid_size;
  step_u.y=ab.y/grid_size;
  step_u.z=ab.z/grid_size;
  step_v.x=bc.x/grid_size;
  step_v.y=bc.y/grid_size;
  step_v.z=bc.z/grid_size;
  vectorize(&step_u,&step_v,&diagonal);

  no_printf("  Area = %.9g\n",area);
  no_printf("  grid_size = %d\n",grid_size);

  r1=wp->vert[0]->x;
  r2=wp->vert[0]->y;
  r3=wp->vert[0]->z;

  nl=1+n_ligand_types;
  n_ligs=0;

  tiles=ep->tiles;
  dsp=ep->desired_state;
  tsp=ep->time_stamp;
  orp=ep->orient;

  for (k=0;k<NUM_ADD_EFFECTORS;k++) {
    rx[k]=NULL;
    prob[k]=0.0;
    orientation[k]=0;
  }

  nr=0;
  tot_prob=0;
  tot_density=0;
  prev_state_flag=0;
  effdp=effdp_head;
  while (effdp!=NULL) {
    no_printf("  adding effector state %s to wall at density %.9g\n",effdp->rx->sym->name,effdp->quantity);
    tot_prob+=(area*effdp->quantity)/(n*effector_grid_density);
    prob[nr]=tot_prob;
    orientation[nr]=effdp->orient;
    rx[nr++]=effdp->rx;
    tot_density+=effdp->quantity;
    prev_state_flag=prev_state_flag||effdp->rx->parent_rx->prev_state_flag;
    effdp=effdp->next;
  }

  no_printf("Number of effector types in wall = %d\n",nr);
  fflush(log_file);
  if (tot_density>effector_grid_density) {
    fprintf(log_file,"\nMCell: Warning -- Total effector density too high: %f\n\n",tot_density);
    fflush(log_file);
/*
    return(1);
*/
  }

  if (prev_state_flag) {
    if ((psp=(unsigned short *)malloc(n*sizeof(unsigned short)))==NULL) {
      fprintf(log_file,"MCell: cannot store previous state data for effector grid: %d\n",n);
      fflush(log_file);
      return(1);
    }
  }
  else {
    psp=NULL;
  }

  k=0;
  n_occupied=0;
  for (uu=0;uu<grid_size;uu++) {
    p0.x=(uu+1)*step_u.x;
    p0.y=(uu+1)*step_u.y;
    p0.z=(uu+1)*step_u.z;
    if (grid_shape==RECTANGULAR) {
      vv_max=2*grid_size;
    }
    else {
      vv_max=(2*uu)+1;
    }
    v_val=0;
    for (vv=0;vv<vv_max;vv++) {
      p1.x=p0.x+(v_val*step_v.x);
      p1.y=p0.y+(v_val*step_v.y);
      p1.z=p0.z+(v_val*step_v.z);
      if (vv%2==0) {
        p2.x=p1.x+(0.33*diagonal.x);
        p2.y=p1.y+(0.33*diagonal.y);
        p2.z=p1.z+(0.33*diagonal.z);
      }
      else {
        p2.x=p1.x+(0.66*diagonal.x);
        p2.y=p1.y+(0.66*diagonal.y);
        p2.z=p1.z+(0.66*diagonal.z);
        v_val++;
      }
      if (prev_state_flag) {
        psp[k]=0;
      }
      
      if (chkpt_init) {
        l=0;
        p_index=-1;
        random_number_use++;
        ran4(&seed,rand,1,1.0);
        while (l<nr && p_index==-1) {
          if (rand[0]<=prob[l++]) {
	    p_index=l-1;
          }
        }
        if (p_index!=-1) {
          n_occupied++;
          tiles[k]=rx[p_index];
          dsp[k]=rx[p_index]->state_index;
          orp[k]=orientation[p_index];
          rx[p_index]->count++;
          for (m=1;m<nl;m++) {
	    for (p=0;p<rx[p_index]->bound_ligands[m];p++) {
	      if ((lp=(struct ligand *)malloc
	           (sizeof(struct ligand)))==NULL) {
	        fprintf(log_file,"MCell: cannot store molecule\n");
	        return(1);
	      }
	      lp->next_ligand=ligand_table[m]->top;
	      ligand_table[m]->top=lp;
	      tot_mols++;
	      ligand_table[m]->n_mols++;
              n_ligs++;
	      lp->lig_num=ligand_table[m]->lig_index++;
	      
	      /* put ligand in middle of effector */
	      lp->pos.x=r1+p2.x+tiny_x;
	      lp->pos.y=r2+p2.y+tiny_y;
	      lp->pos.z=r3+p2.z+tiny_z;
	      
	      lp->effector=ep;
	      lp->index=k;
              subvol=find_subvol(volume,&lp->pos,-1);
              if (subvol<0) {
	        fprintf(log_file,"MCell: invalid subvolume referenced while initializing effectors\n");
                return(1);
              }
  
	      lp->subvol=subvol;
	    }
          }
        }
      }
      k++;
    }
  }
  no_printf("  final index = %d\n",k);

  for (k=0;k<nr;k++) {
    no_printf("Number of effector type %d = %d\n",k+1,rx[k]->count);
  }
  no_printf("Number of molecules initialized on effector sites = %d\n",n_ligs);
  no_printf("Done initializing %d effectors by density\n",n_ran);
  fflush(log_file);

  ep->n_occupied=n_occupied;
  
  ep->prev_state=psp;
  ep->n_types=nr;
	
  return(0);
}


/**
 * Creates an array pointers to reference effector elements of the
 * wp linked list (walls) directly.
 * This effectively flattens the linked list.
 */
int init_effector_table(struct wall *wp)
{
  if (n_effector_grids>0) {
    if ((effector_table=(struct effector **)malloc
	 ((n_effector_grids)*sizeof(struct effector *)))==NULL) {
      return(1);
    }
    while (wp) {
      if (wp->effectors) {
	effector_table[wp->effectors->index]=wp->effectors;
       }
      wp=wp->next_wall;
    }
  }
  return(0);
}


/**
 * Creates an array pointers to reference elements of the
 * parent_rx linked list directly.
 * This effectively flattens the linked list.
 */
int init_rx_table(struct parent_rx *prxp)
{
  
  if (n_rx_types>0) {
    if ((rx_table=(struct parent_rx **)malloc
         (n_rx_types*sizeof(struct parent_rx *)))==NULL) {
      return(1);
    }
    while (prxp) {
      rx_table[prxp->rx_index]=prxp;
      prxp=prxp->next;
    }
  }
  return(0);
}

/**
 * Creates an array pointers to reference elements of the
 * release_event_queue linked list directly.
 * This effectively flattens the linked list.
 */
int init_release_event_table(struct release_event_queue *reqp)
{
  
  if (n_release_events>0) {
    if ((release_event_table=(struct release_event_queue **)malloc
         (n_release_events*sizeof(struct release_event_queue *)))==NULL) {
      return(1);
    }
    while (reqp) {
      release_event_table[reqp->index]=reqp;
      reqp=reqp->next;
    }
  }
  return(0);
}

/**
 * Frees up the memory used by the value of a symbol table entry
 * after we're done with it.
 * Currently deals with OBJ::POLY_OBJ::{BOX_POLY|ORDERED_POLY} types.
 */
void destroy_sym_value(struct sym_table *sym) {
  struct object *objp;
  struct polygon_object *pop;
  struct ordered_poly *opp;
  struct box_poly *bpp;
  struct vector3 *vect3;
  int i;

no_printf("Starting GARBAGE COLLECTING\n");
no_printf("sym: %s\n", sym->name);
  switch (sym->sym_type) {
	case OBJ:
      objp = (struct object *)sym->value;

      switch (objp->object_type) {
		case POLY_OBJ:
		  pop = (struct polygon_object *)objp->obj;
		  switch (pop->list_type) {
			case ORDERED_POLY:
              opp = (struct ordered_poly *)pop->polygon_data;
			  for (i=0;i<opp->n_verts;i++) {
				free(opp->vertex[i]);
			  }
			  free(opp->vertex);
			  if (opp->normal != NULL) {
				  for (i=0;i<opp->n_verts;i++) {
					free(opp->normal[i]);
				  }
				  free(opp->normal);
			  }
			  free(opp);
			break;

			case BOX_POLY:
			  bpp = (struct box_poly *)pop->polygon_data;
			  free(bpp->llf);
			  free(bpp->urb);
			  free(bpp);
			break;
		  }
		  free(pop);
		break;
	  }
	break;
  }

no_printf("Done GARBAGE COLLECTING\n");
fflush(stdout);
}




/* **************************************************************** */

#endif



