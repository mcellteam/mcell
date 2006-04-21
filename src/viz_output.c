
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "mcell_structs.h"
#include "grid_util.h"
#include "sched_util.h"
#include "viz_output.h"
#include <string.h>
#include "strfunc.h"
#include "util.h"

extern struct volume *world;
/* last iteration number in the simulation */
static long long last_iteration = 0; 
static int obj_to_show_number; /* number of viz_obj objects in the world */
static int eff_to_show_number; /* number of types of effectors */ 
static int mol_to_show_number; /* number of types of 3D mol's */

/* these arrays are used to hold values of main_index for objects */
u_int *surf_states = NULL;
u_int *surf_pos = NULL;
u_int *surf_con = NULL;
u_int *surf_field_indices = NULL;
u_int *eff_pos = NULL;
u_int *eff_orient = NULL;
u_int *eff_states = NULL;
u_int *eff_field_indices = NULL;
u_int *mol_pos = NULL;
u_int *mol_orient = NULL;
u_int *mol_states = NULL;
u_int *mol_field_indices = NULL;
/* array of viz_objects names */
char **obj_names = NULL;
/* array of volume_molecule names */
char **mol_names = NULL;
/* array of surface_molecule (effectors) names */
char **eff_names = NULL;
/* 2-dimensional array of regions names for each object. 
   First index refers to the object, second - to the regions. */
char  ***region_names = NULL; 
/* array that holds number of regions for each object */
int *obj_num_regions = NULL;
/* 2-dimensional array of regions indices for each object. 
   First index refers to the object, second - to the regions indices. */
int  **surf_region_values = NULL; 

/**************************************************************************
update_frame_data_list:
	In: struct frame_data_list * fdlp
        Out: Nothing. Calls output visualization functions if necessary.
	     Updates value of the current iteration step and pointer
             to the current iteration in the linked list
**************************************************************************/
void update_frame_data_list(struct frame_data_list *fdlp)
{
  FILE *log_file;

  log_file=world->log_file;
  
  while (fdlp!=NULL) {
    if(world->it_time==fdlp->viz_iterationll)
    {
      
      switch (world->viz_mode)
      {
	case DX_MODE:
          output_dx_objects(fdlp);
	  break;
	case DREAMM_V3_MODE:
          output_dreamm_objects(fdlp);
	  break;
	case RK_MODE:
          output_rk_custom(fdlp);
	  break;
	case ASCII_MODE:
	  output_ascii_molecules(fdlp);
	  break;
	case NO_VIZ_MODE:
	default:
	  /* Do nothing for vizualization */
	  break;
      }
      fdlp->curr_viz_iteration=fdlp->curr_viz_iteration->next;
      if (fdlp->curr_viz_iteration!=NULL) {
	switch (fdlp->list_type) {
	  case OUTPUT_BY_ITERATION_LIST:
	  	  fdlp->viz_iterationll=(long long)fdlp->curr_viz_iteration->value; 
	          break;
	  case OUTPUT_BY_TIME_LIST:
	          fdlp->viz_iterationll=(long long)(fdlp->curr_viz_iteration->value/world->time_unit + ROUND_UP);
	          break;
          default:
                  fprintf(log_file,"MCell: error - wrong frame_data_list list_type %d\n", fdlp->list_type);
                  break;
	}
      }
    }

    fdlp=fdlp->next;
  }
  return;
}

/********************************************************************* 
init_frame_data:
   In: struct frame_data_list* fdlp
   Out: nothing.  Initializes  frame_data_list structure. 
   	Sets the value of the current iteration step to the start value.
   	Sets the number of iterations. 
        Initializes parameters used in 
       'output_dreamm_objects_some_frame_data()'. 
***********************************************************************/
void init_frame_data_list(struct frame_data_list *fdlp)
{
  struct num_expr_list *nelp;
  int done, ii, jj;
  struct species *sp;
  int mol_orient_frame_present = 0;
  int mol_pos_frame_present = 0;
  int reg_data_frame_present = 0;
  int mesh_geometry_frame_present = 0;
  FILE *log_file;

  log_file=world->log_file;

  while (fdlp!=NULL) {
 
    fdlp->viz_iterationll=-1;
    fdlp->n_viz_iterations=0;
    nelp=fdlp->iteration_list;
    done=0;
    switch (fdlp->list_type) {
    case OUTPUT_BY_ITERATION_LIST:
      while (nelp!=NULL) {
	fdlp->n_viz_iterations++;
	if (!done) {
	  if (nelp->value>=world->start_time) {
             fdlp->viz_iterationll=(long long)nelp->value;
             fdlp->curr_viz_iteration=nelp;
             done=1;
	  }
	}
        
        if((long long)(nelp->value) > last_iteration){
             last_iteration = (long long)(nelp->value);
        }
         
	nelp=nelp->next;
      }
      break;
    case OUTPUT_BY_TIME_LIST:
      while (nelp!=NULL) {
	fdlp->n_viz_iterations++;
	if (!done) {
	  if (nelp->value>=world->current_start_real_time) {
	    fdlp->viz_iterationll=(long long)(nelp->value/world->time_unit+ROUND_UP);
	    fdlp->curr_viz_iteration=nelp;
	    done=1;
	  }
        }
        if((long long)(nelp->value/world->time_unit + ROUND_UP) > last_iteration){
             last_iteration = (long long)(nelp->value);
        }
	nelp=nelp->next;
      }
      break;
     default:
         fprintf(stderr,"MCell: error - wrong frame_data_list list_type %d\n", fdlp->list_type);
         break;
    }
        
    if(fdlp->type == MOL_ORIENT) {
        mol_orient_frame_present = 1;
    }
    if(fdlp->type == MOL_POS) {
        mol_pos_frame_present = 1;
    }
    if(fdlp->type == MESH_GEOMETRY) {
        mesh_geometry_frame_present = 1;
    }
    if(fdlp->type == REG_DATA) {
        reg_data_frame_present = 1;
    }

    fdlp=fdlp->next;
  }

  if((mol_orient_frame_present) & (!mol_pos_frame_present)){
     fprintf(log_file, "The input file contains ORIENTATIONS but not POSITIONS statement in the MOLECULES block. The molecules cannot be visualized.\n");
  }
  if((reg_data_frame_present) & (!mesh_geometry_frame_present)){
     fprintf(log_file, "The input file contains REGION_DATA but not GEOMETRY statement in the MESHES block. The meshes cannot be visualized.\n");
  }

  /* find out the number of objects to visualize. */
  struct viz_obj *vizp = world->viz_obj_head;
  while(vizp != NULL){
        struct viz_child *vcp = vizp->viz_child_head;
        while(vcp != NULL){
	   obj_to_show_number++;

           vcp = vcp->next;
        } /* end while (vcp) */
	vizp = vizp->next;
  } /* end while (vizp) */

  /* fill the obj_num_regions array */
  if (obj_to_show_number==0) obj_num_regions = (int*)malloc(sizeof(int));
  else obj_num_regions = (int *)malloc(obj_to_show_number * sizeof(int)); 
  if(obj_num_regions == NULL){
     fprintf(world->err_file, "File %s, Line %ld: memory allocation error.\n", __FILE__, (long)__LINE__);
     exit(1);
  }
  vizp = world->viz_obj_head;
  ii = 0;
  while(vizp != NULL){
        struct viz_child *vcp = vizp->viz_child_head;
        while(vcp != NULL){
           /* find out the number of regions in the object */ 
           obj_num_regions[ii] = vcp->obj->num_regions - 1;
           ii++;
           vcp = vcp->next;
        } /* end while (vcp) */
	vizp = vizp->next;
  } /* end while (vizp) */

  
  /* declare and initialize 2-dimensional array that will hold
     indices of the "region_data" objects which are components
     of the field object.
     Here first index points to the object and the second index
     - to the array of indices for this object regions */
  if (obj_to_show_number==0) surf_region_values = (int**)malloc(sizeof(int*));
  else surf_region_values = (int **)malloc(obj_to_show_number * sizeof(int*));
  if(surf_region_values == NULL){
     fprintf(world->err_file, "File %s, Line %ld: memory allocation error.\n", __FILE__, (long)__LINE__);
       exit(1);
  }
  /* initialize the array */
  for(ii = 0; ii < obj_to_show_number; ii++)
  {
     surf_region_values[ii] = NULL;

  }
  int n_regs; /* number of regions for the object */
  for(ii = 0; ii < obj_to_show_number; ii++)
  {
    n_regs = obj_num_regions[ii];
    if(n_regs > 0){
       surf_region_values[ii] = (int *)malloc(n_regs*sizeof(int));
       if(surf_region_values[ii] == NULL){
          fprintf(world->err_file, "File %s, Line %ld: memory allocation error.\n", __FILE__, (long)__LINE__);
          exit(1);
       }
       for(jj = 0; jj < n_regs; jj++)
       {
          surf_region_values[ii][jj] = 0;
       }
    }
      
  }


  /* find out number of effectors and 3D molecules to visualize */
  for(ii = 0; ii < world->n_species; ii++)
  {
     sp = world->species_list[ii];
     if((sp->flags & IS_SURFACE) != 0) continue;
     if(strcmp(sp->sym->name, "GENERIC_MOLECULE") == 0) continue;
     if(((sp->flags & ON_GRID) == ON_GRID) && (sp->viz_state > 0)){ 
	  eff_to_show_number++;
     }else if(((sp->flags & NOT_FREE) == 0) && (sp->viz_state > 0)) {
          mol_to_show_number++;

     }    
  }
 

  return;
}

/*************************************************************************
output_dx_objects:
	In: struct frame_data_list *fdlp
	Out: 0 on success, 1 on error; output visualization files (*.dx)
             are written.
**************************************************************************/

int output_dx_objects(struct frame_data_list *fdlp)
{
  FILE *log_file;
  FILE *wall_verts_header;
  FILE *wall_states_header;
  FILE *eff_pos_header;
  FILE *eff_states_header;
  FILE *mol_pos_header;
  FILE *mol_states_header;
  struct viz_obj *vizp;
  struct viz_child *vcp;
  struct wall *w,**wp;
  struct surface_grid *sg;
  struct grid_molecule *gmol;
  struct species **species_list;
  struct object *objp;
  struct polygon_object *pop;
  struct ordered_poly *opp;
  struct element_data *edp;
  struct vector3 p0;
  struct storage_list *slp;
  struct storage *sp;
  struct schedule_helper *shp;
  struct abstract_molecule *amp;
  struct molecule *molp,***viz_molp = NULL;
  float v1,v2,v3;
  u_int n_tiles,spec_id,*viz_mol_count = NULL;
  u_int mol_pos_index,mol_pos_field_index,mol_pos_group_index;
  u_int mol_states_index,mol_states_group_index;
  int ii,jj;
  int vi1,vi2,vi3;
  /*int vi4;*/
  int num;
  long long viz_iterationll;
  long long n_viz_iterations;
  /* int first_viz_iteration; */
  int pos_count,state_count,element_data_count;
  int effector_state_pos_count;
  int state;
  int viz_type;
  unsigned int index,n_eff;
  int word;
  byte *word_p;
  byte viz_eff,viz_mol,viz_surf;
  byte viz_surf_or_eff;
  byte viz_eff_pos,viz_eff_states;
  byte viz_mol_pos,viz_mol_states;
  byte viz_surf_pos,viz_surf_states;
  char file_name[1024];
  char my_byte_order[8];

  u_int n_species = world->n_species;
  
  log_file=world->log_file;
  no_printf("Viz output in DX mode...\n");

  word_p=(unsigned char *)&word;
  word=0x04030201;

  if (word_p[0]==1) {
    sprintf(my_byte_order,"lsb");
  }
  else {
    sprintf(my_byte_order,"msb");
  }

  viz_iterationll=fdlp->viz_iterationll;
  n_viz_iterations=fdlp->n_viz_iterations;
  /* first_viz_iteration=(viz_iteration==fdlp->iteration_list->value); */
  
 /* here is the check to prevent writing twice the same info 
    - at the end of one checkpoint and the beginning of the next checkpoint.
   */
  if((world->chkpt_flag) && (world->start_time > 0)){
     if (world->it_time == world->start_time){
         if(viz_iterationll % (world->chkpt_iterations) == 0){
             return 0;
         }
     }
  }

  viz_type=fdlp->type;
  viz_eff=((viz_type==ALL_FRAME_DATA) || (viz_type==EFF_POS)
               || (viz_type==EFF_STATES));
  viz_mol=((viz_type==ALL_FRAME_DATA) || (viz_type==MOL_POS)
               || (viz_type==MOL_STATES));
  viz_surf=((viz_type==ALL_FRAME_DATA) || (viz_type==SURF_POS)
               || (viz_type==SURF_STATES));
  viz_surf_or_eff=(viz_surf || viz_eff);


  viz_eff_pos=((viz_type==ALL_FRAME_DATA) || (viz_type==EFF_POS));
  viz_eff_states=((viz_type==ALL_FRAME_DATA) || (viz_type==EFF_STATES));
  viz_mol_pos=((viz_type==ALL_FRAME_DATA) || (viz_type==MOL_POS));
  viz_mol_states=((viz_type==ALL_FRAME_DATA) || (viz_type==MOL_STATES));
  viz_surf_pos=((viz_type==ALL_FRAME_DATA) || (viz_type==SURF_POS));
  viz_surf_states=((viz_type==ALL_FRAME_DATA) || (viz_type==SURF_STATES));


/* dump walls and effectors: */
  if (viz_surf_or_eff) {
  vizp = world->viz_obj_head;
  while(vizp!=NULL) {

    wall_verts_header=NULL;
    wall_states_header=NULL;
    eff_pos_header=NULL;
    eff_states_header=NULL;

    if (viz_surf_pos) {
      sprintf(file_name,"%s.mesh_elements.%lld.dx",
               vizp->name,viz_iterationll);
      if ((wall_verts_header=fopen(file_name,"wb"))==NULL) {
        fprintf(world->err_file,"File %s, Line %ld: error cannot open mesh elements file %s\n", __FILE__, (long)__LINE__, file_name);
        return(1);
      }
    }

    if (viz_surf_states) {
      sprintf(file_name,"%s.mesh_element_states.%lld.dx",
               vizp->name,viz_iterationll);
      if ((wall_states_header=fopen(file_name,"wb"))==NULL) {
        fprintf(world->err_file,"File %s, Line %ld: error cannot open mesh element states file %s\n", __FILE__, (long)__LINE__, file_name);
        return(1);
      }
    }

    if (viz_eff_pos) {
      sprintf(file_name,"%s.effector_site_positions.%lld.dx",vizp->name,viz_iterationll);
      if ((eff_pos_header=fopen(file_name,"wb"))==NULL) {
        fprintf(world->err_file,"File %s, Line %ld: error cannot open effector position file %s\n", __FILE__, (long)__LINE__, file_name);
        return(1);
      }
    }
  
    if (viz_eff_states) {
      sprintf(file_name,"%s.effector_site_states.%lld.dx",vizp->name,viz_iterationll);
      if ((eff_states_header=fopen(file_name,"wb"))==NULL) {
        fprintf(world->err_file,"File %s, Line %ld: error cannot open effector states file %s\n", __FILE__, (long)__LINE__, file_name);
        return(1);
      }
    }

    /* Traverse all visualized compartments */
    /* output mesh element positions and connections */
    /* output effector positions */
    /* output effector states */
    /* do not output effector normals header or normals data yet */
/*
    wall_verts_index = 1;
    wall_states_index = 1;
    wall_data_index = 1;
*/
/*
    poly_pos_count = 0;
*/
/*
    eff_pos_index = 1;
    eff_states_index = 1;
*/
    pos_count = 0;
    state_count = 0;
    effector_state_pos_count = 0;
    vcp = vizp->viz_child_head;
    while(vcp!=NULL) {
      objp = vcp->obj;
      pop=(struct polygon_object *)objp->contents;
      if (objp->object_type==BOX_OBJ) {
#if 0
        if (viz_surf_pos || viz_surf_states)
        {
          element_data_count=0.5*objp->n_walls_actual;
        
          if (viz_surf_pos && !viz_surf_states) {
            fprintf(wall_verts_header,
              "object \"%s.positions\" class array type float rank 1 shape 3 items %d %s binary data follows\n",
              objp->sym->name,objp->n_verts,my_byte_order);
            /* output box vertices */
            for (ii=0;ii<objp->n_verts;ii++) {
              v1 = world->length_unit*objp->verts[ii].x;
              v2 = world->length_unit*objp->verts[ii].y;
              v3 = world->length_unit*objp->verts[ii].z;
              fwrite(&v1,sizeof v1,1,wall_verts_header);
              fwrite(&v2,sizeof v2,1,wall_verts_header);
              fwrite(&v3,sizeof v3,1,wall_verts_header);
            }
            fprintf(wall_verts_header,
              "\nattribute \"dep\" string \"positions\"\n#\n");
/*
            wall_verts_index++;
*/
    
            /* output box wall connections */
            fprintf(wall_verts_header,
              "object \"%s.connections\" class array type int rank 1 shape 4 items %d %s binary data follows\n",
              objp->sym->name,element_data_count,my_byte_order);
            for (ii=0;ii<objp->n_walls;ii+=2) {
              if (pop->side_stat[ii]) {
                switch (ii) {
                  case TP:
                    vi1=3;
                    vi2=7;
                    vi3=1;
                    vi4=5;
                    fwrite(&vi1,sizeof vi1,1,wall_verts_header);
                    fwrite(&vi2,sizeof vi2,1,wall_verts_header);
                    fwrite(&vi3,sizeof vi3,1,wall_verts_header);
                    fwrite(&vi4,sizeof vi4,1,wall_verts_header);
                  break;
                  case BOT:
                    vi1=0;
                    vi2=4;
                    vi3=2;
                    vi4=6;
                    fwrite(&vi1,sizeof vi1,1,wall_verts_header);
                    fwrite(&vi2,sizeof vi2,1,wall_verts_header);
                    fwrite(&vi3,sizeof vi3,1,wall_verts_header);
                    fwrite(&vi4,sizeof vi4,1,wall_verts_header);
                  break;
                  case FRNT:
                    vi1=4;
                    vi2=0;
                    vi3=5;
                    vi4=1;
                    fwrite(&vi1,sizeof vi1,1,wall_verts_header);
                    fwrite(&vi2,sizeof vi2,1,wall_verts_header);
                    fwrite(&vi3,sizeof vi3,1,wall_verts_header);
                    fwrite(&vi4,sizeof vi4,1,wall_verts_header);
                  break;
                  case BCK:
                    vi1=2;
                    vi2=6;
                    vi3=3;
                    vi4=7;
                    fwrite(&vi1,sizeof vi1,1,wall_verts_header);
                    fwrite(&vi2,sizeof vi2,1,wall_verts_header);
                    fwrite(&vi3,sizeof vi3,1,wall_verts_header);
                    fwrite(&vi4,sizeof vi4,1,wall_verts_header);
                  break;
                  case LFT:
                    vi1=0;
                    vi2=2;
                    vi3=1;
                    vi4=3;
                    fwrite(&vi1,sizeof vi1,1,wall_verts_header);
                    fwrite(&vi2,sizeof vi2,1,wall_verts_header);
                    fwrite(&vi3,sizeof vi3,1,wall_verts_header);
                    fwrite(&vi4,sizeof vi4,1,wall_verts_header);
                  break;
                  case RT:
                    vi1=6;
                    vi2=4;
                    vi3=7;
                    vi4=5;
                    fwrite(&vi1,sizeof vi1,1,wall_verts_header);
                    fwrite(&vi2,sizeof vi2,1,wall_verts_header);
                    fwrite(&vi3,sizeof vi3,1,wall_verts_header);
                    fwrite(&vi4,sizeof vi4,1,wall_verts_header);
                  break;
                }
              }
            }
            fprintf(wall_verts_header,
              "\nattribute \"ref\" string \"positions\"\n");
            fprintf(wall_verts_header,
              "attribute \"element type\" string \"quads\"\n#\n");
/*
            wall_verts_index++;
*/
          }

          if (viz_surf_pos && viz_surf_states) {
            fprintf(wall_verts_header,
              "object \"%s.positions\" class array type float rank 1 shape 3 items %d %s binary data follows\n",
              objp->sym->name,objp->n_verts,my_byte_order);
            /* output box vertices */
            for (ii=0;ii<objp->n_verts;ii++) {
              v1 = world->length_unit*objp->verts[ii].x;
              v2 = world->length_unit*objp->verts[ii].y;
              v3 = world->length_unit*objp->verts[ii].z;
              fwrite(&v1,sizeof v1,1,wall_verts_header);
              fwrite(&v2,sizeof v2,1,wall_verts_header);
              fwrite(&v3,sizeof v3,1,wall_verts_header);
            }
            fprintf(wall_verts_header,
              "\nattribute \"dep\" string \"positions\"\n#\n");
/*
            wall_verts_index++;
*/
    
            /* output box wall connections */
            fprintf(wall_verts_header,
              "object \"%s.connections\" class array type int rank 1 shape 4 items %d %s binary data follows\n",
              objp->sym->name,element_data_count,my_byte_order);
            fprintf(wall_states_header,
              "object \"%s.states\" class array type int rank 0 items %d ascii data follows\n",
              objp->sym->name,element_data_count);
            for (ii=0;ii<objp->n_walls;ii+=2) {
              if (pop->side_stat[ii]) {
                switch (ii) {
                  case TP:
                    vi1=3;
                    vi2=7;
                    vi3=1;
                    vi4=5;
                    fwrite(&vi1,sizeof vi1,1,wall_verts_header);
                    fwrite(&vi2,sizeof vi2,1,wall_verts_header);
                    fwrite(&vi3,sizeof vi3,1,wall_verts_header);
                    fwrite(&vi4,sizeof vi4,1,wall_verts_header);
                  break;
                  case BOT:
                    vi1=0;
                    vi2=4;
                    vi3=2;
                    vi4=6;
                    fwrite(&vi1,sizeof vi1,1,wall_verts_header);
                    fwrite(&vi2,sizeof vi2,1,wall_verts_header);
                    fwrite(&vi3,sizeof vi3,1,wall_verts_header);
                    fwrite(&vi4,sizeof vi4,1,wall_verts_header);
                  break;
                  case FRNT:
                    vi1=4;
                    vi2=0;
                    vi3=5;
                    vi4=1;
                    fwrite(&vi1,sizeof vi1,1,wall_verts_header);
                    fwrite(&vi2,sizeof vi2,1,wall_verts_header);
                    fwrite(&vi3,sizeof vi3,1,wall_verts_header);
                    fwrite(&vi4,sizeof vi4,1,wall_verts_header);
                  break;
                  case BCK:
                    vi1=2;
                    vi2=6;
                    vi3=3;
                    vi4=7;
                    fwrite(&vi1,sizeof vi1,1,wall_verts_header);
                    fwrite(&vi2,sizeof vi2,1,wall_verts_header);
                    fwrite(&vi3,sizeof vi3,1,wall_verts_header);
                    fwrite(&vi4,sizeof vi4,1,wall_verts_header);
                  break;
                  case LFT:
                    vi1=0;
                    vi2=2;
                    vi3=1;
                    vi4=3;
                    fwrite(&vi1,sizeof vi1,1,wall_verts_header);
                    fwrite(&vi2,sizeof vi2,1,wall_verts_header);
                    fwrite(&vi3,sizeof vi3,1,wall_verts_header);
                    fwrite(&vi4,sizeof vi4,1,wall_verts_header);
                  break;
                  case RT:
                    vi1=6;
                    vi2=4;
                    vi3=7;
                    vi4=5;
                    fwrite(&vi1,sizeof vi1,1,wall_verts_header);
                    fwrite(&vi2,sizeof vi2,1,wall_verts_header);
                    fwrite(&vi3,sizeof vi3,1,wall_verts_header);
                    fwrite(&vi4,sizeof vi4,1,wall_verts_header);
                  break;
                }
                state=objp->viz_state[ii];
                fprintf(wall_states_header,"%d\n",state);
              }
            }
            fprintf(wall_verts_header,
              "\nattribute \"ref\" string \"positions\"\n");
            fprintf(wall_verts_header,
              "attribute \"element type\" string \"quads\"\n#\n");
/*
            wall_verts_index++;
*/
            fprintf(wall_states_header,"\nattribute \"dep\" string \"connections\"\n#\n");
/*
            wall_states_index++;
*/

          }

          if (viz_surf_states && !viz_surf_pos) {
            fprintf(wall_states_header,
              "object \"%s.states\" class array type int rank 0 items %d ascii data follows\n",
              objp->sym->name,element_data_count);
            for (ii=0;ii<objp->n_walls;ii+=2) {
              if (!get_bit(pop->side_removed,ii)) {
                state=objp->viz_state[ii];
                fprintf(wall_states_header,"%d\n",state);
              }
            }
          
          }
          fprintf(wall_states_header,"\nattribute \"dep\" string \"connections\"\n#\n");
/*
          wall_states_index++;
*/
        }
#endif
      }

      if (objp->object_type==POLY_OBJ || objp->object_type==BOX_OBJ) {
        if (viz_surf && (viz_surf_pos || viz_surf_states))
        {
          opp=(struct ordered_poly *)pop->polygon_data;
          edp=opp->element;
          element_data_count=objp->n_walls_actual;

          if (viz_surf_pos && !viz_surf_states) {

	    fprintf(wall_verts_header,
              "object \"%s.positions\" class array type float rank 1 shape 3 items %d %s binary data follows\n",
              objp->sym->name,objp->n_verts,my_byte_order);
            /* output polyhedron vertices */
            for (ii=0;ii<objp->n_verts;ii++) {
              v1 = world->length_unit*objp->verts[ii].x;
              v2 = world->length_unit*objp->verts[ii].y;
              v3 = world->length_unit*objp->verts[ii].z;
	      fwrite(&v1,sizeof v1,1,wall_verts_header);
	      fwrite(&v2,sizeof v2,1,wall_verts_header);
	      fwrite(&v3,sizeof v3,1,wall_verts_header);
            }
            fprintf(wall_verts_header,
              "\nattribute \"dep\" string \"positions\"\n#\n");
/*
	    wall_verts_index++;
*/

            /* output polygon element connections */
	    fprintf(wall_verts_header,
              "object \"%s.connections\" class array type int rank 1 shape 3 items %d %s binary data follows\n",
              objp->sym->name,element_data_count,my_byte_order);
            for (ii=0;ii<objp->n_walls;ii++) {
              if (!get_bit(pop->side_removed,ii)) {
                for (jj=0;jj<edp[ii].n_verts-2;jj++) {
                  vi1=edp[ii].vertex_index[0];
                  vi2=edp[ii].vertex_index[jj+1];
                  vi3=edp[ii].vertex_index[jj+2];
	          fwrite(&vi1,sizeof vi1,1,wall_verts_header);
	          fwrite(&vi2,sizeof vi2,1,wall_verts_header);
	          fwrite(&vi3,sizeof vi3,1,wall_verts_header);
                }
              }
            }
            fprintf(wall_verts_header,
              "\nattribute \"ref\" string \"positions\"\n");
            fprintf(wall_verts_header,
              "attribute \"element type\" string \"triangles\"\n#\n");
/*
	    wall_verts_index++;
*/


          }
    
          if (viz_surf_pos && viz_surf_states) {
	    fprintf(wall_verts_header,
              "object \"%s.positions\" class array type float rank 1 shape 3 items %d %s binary data follows\n",
              objp->sym->name,objp->n_verts,my_byte_order);
            /* output polyhedron vertices */
            for (ii=0;ii<objp->n_verts;ii++) {
              v1 = world->length_unit*objp->verts[ii].x;
              v2 = world->length_unit*objp->verts[ii].y;
              v3 = world->length_unit*objp->verts[ii].z;
	      fwrite(&v1,sizeof v1,1,wall_verts_header);
	      fwrite(&v2,sizeof v2,1,wall_verts_header);
	      fwrite(&v3,sizeof v3,1,wall_verts_header);
            }
            fprintf(wall_verts_header,
              "\nattribute \"dep\" string \"positions\"\n#\n");
/*
	    wall_verts_index++;
*/

            /* output polygon element connections */
	    fprintf(wall_verts_header,
              "object \"%s.connections\" class array type int rank 1 shape 3 items %d %s binary data follows\n",
              objp->sym->name,element_data_count,my_byte_order);
            fprintf(wall_states_header,
              "object \"%s.states\" class array type int rank 0 items %d ascii data follows\n",
              objp->sym->name,element_data_count);
            for (ii=0;ii<objp->n_walls;ii++) {
              if (!get_bit(pop->side_removed,ii)) {
                for (jj=0;jj<edp[ii].n_verts-2;jj++) {
                  vi1=edp[ii].vertex_index[0];
                  vi2=edp[ii].vertex_index[jj+1];
                  vi3=edp[ii].vertex_index[jj+2];
	          fwrite(&vi1,sizeof vi1,1,wall_verts_header);
	          fwrite(&vi2,sizeof vi2,1,wall_verts_header);
	          fwrite(&vi3,sizeof vi3,1,wall_verts_header);
                }
		state=objp->viz_state[ii];
		fprintf(wall_states_header,"%d\n",state);
              }
            }
            fprintf(wall_verts_header,
              "\nattribute \"ref\" string \"positions\"\n");
            fprintf(wall_verts_header,
              "attribute \"element type\" string \"triangles\"\n#\n");
/*
	    wall_verts_index++;
*/
            fprintf(wall_states_header,"\nattribute \"dep\" string \"connections\"\n#\n");
/*
            wall_states_index++;
*/
          }

          if (viz_surf_states && !viz_surf_pos) {
          fprintf(wall_states_header,
            "object \"%s.states\" class array type int rank 0 items %d ascii data follows\n",
            objp->sym->name,element_data_count);
          for (ii=0;ii<objp->n_walls;ii++) {
            if (!get_bit(pop->side_removed,ii)) {
              state=objp->viz_state[ii];
              fprintf(wall_states_header,"%d\n",state);
            }
          }
          fprintf(wall_states_header,"\nattribute \"dep\" string \"connections\"\n#\n");
/*
          wall_states_index++;
*/
          }
        }
      }

      wp = objp->wall_p;
      no_printf("Traversing walls in object %s\n",objp->sym->name);
      if (viz_eff) {
        n_eff=0;
        for (ii=0;ii<objp->n_walls;ii++) {
          w = wp[ii];
          if (w!=NULL) {
	    sg = w->effectors;
            if (sg!=NULL) {
              for (index=0;index<sg->n_tiles;index++) {
	        gmol=sg->mol[index];
	        if (gmol!=NULL) {
	          state=sg->mol[index]->properties->viz_state;
	        }
	        else {
	          state=EXCLUDE_OBJ;
	        }
                if (state!=EXCLUDE_OBJ) {
                  n_eff++;
                }
              }
            }
          }
        }
        no_printf("Dumping %d effectors...\n",n_eff);
        fflush(log_file);
        if (viz_eff_pos) {
          if (n_eff) {
            fprintf(eff_pos_header,
              "object \"%s.pos_and_norm\" class array type float rank 2 shape 2 3 items %d %s binary data follows\n",
              objp->sym->name,n_eff,my_byte_order);
/*
            eff_pos_index++;
*/
          }
          else {
            fprintf(eff_pos_header,
              "object \"%s.pos_and_norm\" array",objp->sym->name);
/*
            eff_pos_index++;
*/
          }
        }

        if (viz_eff_states) {
          if (n_eff) {
            fprintf(eff_states_header,
              "object \"%s.states\" class array type int rank 0 items %d ascii data follows\n",
              objp->sym->name,n_eff);
/*
            eff_states_index++;
*/
          }
          else {
            fprintf(eff_states_header,
              "object \"%s.states\" array\n",objp->sym->name);
/*
            eff_states_index++;
*/
          }
        }
      }
      wp = objp->wall_p;
      for (ii=0;ii<objp->n_walls;ii++) {
        w = wp[ii];
        if (w!=NULL) {
	  sg = w->effectors;

          /* dump the effectors */
          if (viz_eff) {
            if (sg!=NULL) {
              n_tiles=sg->n_tiles;
              for (index=0;index<n_tiles;index++) {
                grid2xyz(sg,index,&p0);
	        gmol=sg->mol[index];
	        if (gmol!=NULL) {
	          state=sg->mol[index]->properties->viz_state;
	        }
	        else {
	          state=EXCLUDE_OBJ;
	        }
      
                if (state!=EXCLUDE_OBJ) {
                  if (viz_eff_states) {
                   fprintf(eff_states_header,"%d\n",state);
                  }
                  if (viz_eff_pos) {
	           v1=world->length_unit*p0.x;
	           v2=world->length_unit*p0.y;
	           v3=world->length_unit*p0.z;
	           fwrite(&v1,sizeof v1,1,eff_pos_header);
	           fwrite(&v2,sizeof v2,1,eff_pos_header);
	           fwrite(&v3,sizeof v3,1,eff_pos_header);
                   v1=(float)((w->normal.x)*(gmol->orient));
                   v2=(float)((w->normal.y)*(gmol->orient));
                   v3=(float)((w->normal.z)*(gmol->orient));
                   fwrite(&v1,sizeof v1,1,eff_pos_header);
                   fwrite(&v2,sizeof v2,1,eff_pos_header);
                   fwrite(&v3,sizeof v3,1,eff_pos_header);
                  }
                }
              }
            }
          }
        }
      }
      if (viz_eff_pos) {
        fprintf(eff_pos_header,
          "\n#\n");
      }
      if (viz_eff_states) {
        fprintf(eff_states_header,
          "attribute \"dep\" string \"positions\"\n#\n");
      }
      vcp = vcp->next;
    }

    /* output effector states null object */
    if (viz_eff) {
      if (viz_eff_states) {
        fprintf(eff_states_header,"object \"null_object\" array\n");
        fprintf(eff_states_header,
          "  attribute \"dep\" string \"positions\"\n\n");
      }
    }

    /* output surface states null object */
    if (viz_surf_states) {
      if (wall_states_header!=NULL) {
        fprintf(wall_states_header,"object \"null_object\" array\n");
        fprintf(wall_states_header,
          "  attribute \"dep\" string \"connections\"\n\n");
      }
    }

    /* output effector and surface positions field objects */
    vcp = vizp->viz_child_head;
    while(vcp!=NULL) {
      objp = vcp->obj;

      /* effectors */
      if (viz_eff) {
        /* effector positions */
        if (viz_eff_pos) {
          fprintf(eff_pos_header,"object \"%s.field\" field\n",objp->sym->name);
/*
          fprintf(eff_pos_header,
            "  component \"positions\" \"%s.positions\"\n",objp->sym->name);
*/
          fprintf(eff_pos_header,
            "  component \"data\" \"%s.pos_and_norm\"\n\n",objp->sym->name);
        }
      }

      /* surfaces */
      if (viz_surf) {
        /* surface positions */
        if (viz_surf_pos) {
          if (wall_verts_header!=NULL) {
            fprintf(wall_verts_header,
              "object \"%s.field\" field\n",objp->sym->name);
            fprintf(wall_verts_header,
              "  component \"positions\" \"%s.positions\"\n",objp->sym->name);
            fprintf(wall_verts_header,
              "  component \"connections\" \"%s.connections\"\n\n",
              objp->sym->name);
          }
        }
      }
      vcp = vcp->next;
    }

    /* output effector positions null objects */
    if (viz_eff) {
      if (viz_eff_pos) {
        fprintf(eff_pos_header,"object \"null_pos_and_norm\" array\n\n");
/*
        fprintf(eff_pos_header,"object \"null_data\" array\n\n");
*/
        fprintf(eff_pos_header,"object \"null_object\" field\n");
/*
        fprintf(eff_pos_header,"  component \"positions\" \"null_positions\"\n");
*/
        fprintf(eff_pos_header,"  component \"data\" \"null_pos_and_norm\"\n\n");
      }
    }

    /* output surface positions null objects */
    if (viz_surf_pos) {
      if (wall_verts_header!=NULL) {
        fprintf(wall_verts_header,"object \"null_positions\" array\n\n");
        fprintf(wall_verts_header,"object \"null_connections\" array\n\n");
        fprintf(wall_verts_header,"object \"null_object\" field\n");
        fprintf(wall_verts_header,"  component \"positions\" \"null_positions\"\n");
        fprintf(wall_verts_header,"  component \"connections\" \"null_connections\"\n\n");
      }
    }

    /* output effector group objects and null members*/
    /* group object name is full MCell parent object name */
    if (viz_eff) {
      /* effector positions */
      if (viz_eff_pos) {
        fprintf(eff_pos_header,"object \"%s\" group\n",vizp->full_name);
        fprintf(eff_pos_header,
          "  member \"null_object (default)\" \"null_object\"\n");
      }
      /* effector states */
      if (viz_eff_states) {
        fprintf(eff_states_header,"object \"%s\" group\n",vizp->full_name);
        fprintf(eff_states_header,
          "  member \"null_object (default)\" \"null_object\"\n");
      }
    }

    /* output surface group objects and null members*/
    /* group object name is full MCell parent object name */
    if (viz_surf) {
      /* surface positions */
      if (viz_surf_pos) {
        if (wall_verts_header!=NULL) {
          fprintf(wall_verts_header,"object \"%s\" group\n",vizp->full_name);
          fprintf(wall_verts_header,
            "  member \"null_object (default)\" \"null_object\"\n");
        }
      }
      /* surface states */
      if (viz_surf_states) {
        if (wall_states_header!=NULL) {
          fprintf(wall_states_header,"object \"%s\" group\n",vizp->full_name);
          fprintf(wall_states_header,
            "  member \"null_object (default)\" \"null_object\"\n");
        }
      }
    }

    /* output group object members */
    /* member name is full MCell child object name */
    vcp = vizp->viz_child_head;
    while(vcp!=NULL) {
      objp = vcp->obj;

      /* effectors */
      if (viz_eff) {
        /* effector positions */
        if (viz_eff_pos) {
          fprintf(eff_pos_header,
            "  member \"%s\" \"%s.field\"\n",objp->sym->name,objp->sym->name);
        }
        /* effector states */
        if (viz_eff_states) {
          fprintf(eff_states_header,
            "  member \"%s.states\" \"%s.states\"\n",objp->sym->name,objp->sym->name);
        }
      }

      /* surfaces */
      if (viz_surf) {
        /* surface positions */
        if (viz_surf_pos) {
          if (wall_verts_header!=NULL) {
            fprintf(wall_verts_header,
              "  member \"%s\" \"%s.field\"\n",objp->sym->name,objp->sym->name);
          }
        }
        /* surface states */
        if (viz_surf_states) {
          if (wall_states_header!=NULL) {
            fprintf(wall_states_header,
              "  member \"%s.states\" \"%s.states\"\n",objp->sym->name,objp->sym->name);
          }
        }
      }
      vcp = vcp->next;
    }

    if (wall_verts_header!=NULL) {
      fclose(wall_verts_header);
    }
    if (wall_states_header!=NULL) {
      fclose(wall_states_header);
    }
    if (eff_pos_header!=NULL) {
      fclose(eff_pos_header);
    }
    if (eff_states_header!=NULL) {
      fclose(eff_states_header);
    }
    vizp = vizp->next;
  }

  } /* end viz_surf_or_eff */

/* dump diffusible molecules: */
  if (viz_mol) {
    mol_pos_header = NULL;
    mol_states_header = NULL;
    mol_pos_index=1;
    mol_states_index=1;
    pos_count=0;
    state_count=0;
    if (viz_mol_pos) {
      sprintf(file_name,"%s.molecule_positions.%lld.dx",world->molecule_prefix_name,viz_iterationll);
      if ((mol_pos_header=fopen(file_name,"wb"))==NULL) {
        fprintf(world->err_file,"File %s, Line %ld: error cannot open molecule positions header file %s\n", __FILE__, (long)__LINE__, file_name);
        return(1);
      }
    }
    if (viz_mol_states) {
      sprintf(file_name,"%s.molecule_states.%lld.dx",world->molecule_prefix_name,viz_iterationll);
      if ((mol_states_header=fopen(file_name,"wb"))==NULL) {
        fprintf(world->err_file,"File %s, Line %ld: error cannot open molecule states header file %s\n", __FILE__, (long)__LINE__, file_name);
        return(1);
      }
    }

    species_list=world->species_list;
    n_species=world->n_species;

    if ((viz_molp=(struct molecule ***)malloc(n_species*sizeof(struct molecule **)))==NULL) {
        fprintf(world->err_file,"File %s, Line %ld: memory allocation error.\n", __FILE__, (long)__LINE__);
      return(1);
    }
    if ((viz_mol_count=(u_int *)malloc(n_species*sizeof(u_int)))==NULL) {
        fprintf(world->err_file,"File %s, Line %ld: memory allocation error.\n", __FILE__, (long)__LINE__);
      return(1);
    }

    for (ii=0;ii<n_species;ii++) {

      spec_id=species_list[ii]->species_id;
      viz_molp[spec_id]=NULL;
      viz_mol_count[spec_id]=0;

      if (species_list[ii]->viz_state!=EXCLUDE_OBJ) {
        num=species_list[ii]->population;
        if (num>0) {
          if ((viz_molp[spec_id]=(struct molecule **)malloc
            (num*sizeof(struct molecule *)))==NULL) {
                fprintf(world->err_file,"File %s, Line %ld: memory allocation error.\n", __FILE__, (long)__LINE__);
                return(1);
          }
        }
      }
    }

    slp=world->storage_head;
    while (slp!=NULL) {
      sp=slp->store;
      shp=sp->timer;
      while (shp!=NULL) {

        for (ii=0;ii<shp->buf_len;ii++) {
          amp=(struct abstract_molecule *)shp->circ_buf_head[ii];
          while (amp!=NULL) {
            if ((amp->properties!=NULL) && (amp->flags&TYPE_3D)==TYPE_3D) {
              molp=(struct molecule *)amp;
              if (molp->properties->viz_state!=EXCLUDE_OBJ) {
                spec_id=molp->properties->species_id;
                if (viz_mol_count[spec_id]<molp->properties->population) {
                  viz_molp[spec_id][viz_mol_count[spec_id]++]=molp;
                }
                else {
                  fprintf(log_file,"MCell: molecule count disagreement!!\n");
                  fprintf(log_file,"  Species %s  population = %d  count = %d\n",molp->properties->sym->name,molp->properties->population,viz_mol_count[spec_id]);
                }
              }
            }
            amp=amp->next;
          }
        }

        amp=(struct abstract_molecule *)shp->current;
        while (amp!=NULL) {
          if ((amp->properties!=NULL) && (amp->flags&TYPE_3D)==TYPE_3D) {
            molp=(struct molecule *)amp;
            if (molp->properties->viz_state!=EXCLUDE_OBJ) {
              spec_id=molp->properties->species_id;
              if (viz_mol_count[spec_id]<molp->properties->population) {
                viz_molp[spec_id][viz_mol_count[spec_id]++]=molp;
              }
              else {
                fprintf(log_file,"MCell: molecule count disagreement!!\n");
                fprintf(log_file,"  Species %s  population = %d  count = %d\n",molp->properties->sym->name,molp->properties->population,viz_mol_count[spec_id]);
              }
            }
          }
          amp=amp->next;
        }
        
        shp=shp->next_scale;
      }

      slp=slp->next;
    }

    for (ii=0;ii<n_species;ii++) {
      spec_id=species_list[ii]->species_id;
      num=viz_mol_count[spec_id];
      state=species_list[ii]->viz_state;
      if (state!=EXCLUDE_OBJ
          && num!=species_list[ii]->population
          && ((species_list[ii]->flags & NOT_FREE)==0)) {
        fprintf(log_file,"MCell: molecule count disagreement!!\n");
        fprintf(log_file,"  Species %s  population = %d  count = %d\n",species_list[ii]->sym->name,species_list[ii]->population,num);
      }

      if (num>0 && ((species_list[ii]->flags & NOT_FREE)==0)) {
        if (viz_mol_pos) {
          fprintf(mol_pos_header,"object \"%d\" class array type float rank 1 shape 3 items %d %s binary data follows\n",mol_pos_index,num,my_byte_order);
          for (jj=0;jj<num;jj++) {
            molp=viz_molp[spec_id][jj];
	    v1=world->length_unit*molp->pos.x;
	    v2=world->length_unit*molp->pos.y;
	    v3=world->length_unit*molp->pos.z;
	    fwrite(&v1,sizeof v1,1,mol_pos_header);
	    fwrite(&v2,sizeof v2,1,mol_pos_header);
	    fwrite(&v3,sizeof v3,1,mol_pos_header);
          }
          fprintf(mol_pos_header,"\n#\n");
          mol_pos_index++;
        }
        if (viz_mol_states) {
          fprintf(mol_states_header,"object \"%d\"\n  constantarray type int items %d\n",mol_states_index,num);
          fprintf(mol_states_header,"  data follows %d\n",state);
          fprintf(mol_states_header,"  attribute \"dep\" string \"positions\"\n\n");
          mol_states_index++;
        }
      }
/* output empty arrays for zero molecule counts here */
      else if (num==0 && ((species_list[ii]->flags & NOT_FREE)==0)) {
        if (viz_mol_pos) {
          fprintf(mol_pos_header,"object \"%d\" array\n\n",mol_pos_index);
          mol_pos_index++;
        }
        if (viz_mol_states) {
          fprintf(mol_states_header,"object \"%d\" array\n",mol_states_index);
          fprintf(mol_states_header,"  attribute \"dep\" string \"positions\"\n\n");
          mol_states_index++;
        }
      }
    }

/* build fields and groups here */
    if (viz_mol_pos) {
      mol_pos_field_index=0;
      for (ii=0;ii<n_species;ii++) {
        if ((species_list[ii]->flags & NOT_FREE)==0) {
          fprintf(mol_pos_header,
            "object \"%d\" field\n",mol_pos_index+mol_pos_field_index);
          fprintf(mol_pos_header,
            "  component \"positions\" \"%d\"\n\n",1+mol_pos_field_index);
          mol_pos_field_index++;
        }
      }
      fprintf(mol_pos_header,"object \"null_positions\" array\n\n");
      fprintf(mol_pos_header,"object \"null_object\" field\n");
      fprintf(mol_pos_header,"  component \"positions\" \"null_positions\"\n\n");
      fprintf(mol_pos_header,"object \"%d\" group\n",
        2*mol_pos_index-1);
      fprintf(mol_pos_header,"  member \"null_object (default)\" \"null_object\"\n");
      mol_pos_group_index=0;
      for (ii=0;ii<n_species;ii++) {
        if ((species_list[ii]->flags & NOT_FREE)==0) {
          fprintf(mol_pos_header,"  member \"%s\" \"%d\"\n",species_list[ii]->sym->name,mol_pos_index+mol_pos_group_index);
          mol_pos_group_index++;
        }
      }
      fclose(mol_pos_header);
    }
    if (viz_mol_states) {
      fprintf(mol_states_header,"object \"null_object\" array\n");
      fprintf(mol_states_header,"  attribute \"dep\" string \"positions\"\n\n");
      fprintf(mol_states_header,"object \"%d\" group\n",mol_states_index);
      fprintf(mol_states_header,"  member \"null_object (default)\" \"null_object\"\n");
      mol_states_group_index=0;
      for (ii=0;ii<n_species;ii++) {
        if ((species_list[ii]->flags & NOT_FREE)==0) {
          fprintf(mol_states_header,"  member \"%s\" \"%d\"\n",species_list[ii]->sym->name,1+mol_states_group_index);
          mol_states_group_index++;
        }
      }
      fclose(mol_states_header);
    }

  }

  if (viz_molp != NULL) {
    for (ii=0;ii<n_species;ii++) {
      if (viz_molp[ii]!=NULL) {
        free(viz_molp[ii]);
      }
    }
    
    free(viz_molp);
  }
  if (viz_mol_count != NULL){
    free (viz_mol_count);
  }
  return(0);
}

/*************************************************************************
output_dreamm_objects:
	In: struct frame_data_list *fdlp
	Out: 0 on success, 1 on error; output visualization files (*.dx)
             in dreamm group format are written.
        NB! NOT TESTED FOR SURFACE_MOLECULES!!!
**************************************************************************/

int output_dreamm_objects(struct frame_data_list *fdlp)
{
  FILE *log_file;
  FILE *master_header = NULL;
  FILE *mesh_pos_data = NULL;  /* data file for wall vertices */
  FILE *mesh_states_data = NULL; /* data file for wall states */
  FILE *region_data = NULL; /* data file for region's data */
  FILE *mol_pos_data = NULL;	/* data file for molecule positions */
  FILE *mol_states_data = NULL; /* data file for molecule states */
  FILE *mol_orient_data = NULL; /* data file for molecule orientations */
  FILE *frame_numbers_data = NULL; /* data file for frame numbers */
  FILE *time_values_data = NULL; /* data file for time_values */
  struct viz_obj *vizp = NULL;
  struct viz_child *vcp = NULL;
  struct surface_grid *sg;
  struct wall *w,**wp;
  struct species **species_list;
  struct species *specp;
  struct object *objp;
  struct polygon_object *pop;
  struct ordered_poly *opp;
  struct element_data *edp;
  struct vector3 p0;
  struct storage_list *slp;
  struct storage *sp;
  struct schedule_helper *shp;
  struct abstract_molecule *amp;
  struct grid_molecule *gmol;
  struct molecule *molp,***viz_molp = NULL;       /* for 3D molecules */
  struct region *rp;
  struct region_list *rlp;
  float v1,v2,v3;
  u_int spec_id = 0, *viz_mol_count = NULL, *viz_grid_mol_count = NULL;
  static u_int main_index = 1;
  static u_int series_index = 0;
  int ii,jj;
  int vi1,vi2,vi3;
  u_int num;
  u_int viz_iteration;
  u_int n_viz_iterations;
  int element_data_count;
  int state;
  int viz_type;
  unsigned int index; 
  /* indices used in arrays "u_int *surf_pos", etc. */
  int surf_pos_index = 0; 
  int surf_con_index = 0;
  int surf_states_index = 0;
  int surf_region_values_index = 0;
  int surf_obj_region_values_index = 0;
  int eff_states_index = 0;
  int eff_pos_index = 0;
  int eff_orient_index = 0;
  int mol_states_index = 0;
  int mol_pos_index = 0;
  int mol_orient_index = 0;
  int word;
  byte *word_p;
  byte viz_mol_pos_flag = 0, viz_mol_states_flag = 0;	/* flags */
  byte viz_mol_orient_flag = 0, viz_region_data_flag = 0;   /* flags */
  byte viz_mol_all_data_flag = 0, viz_surf_all_data_flag = 0;	/* flags */
  byte viz_surf_pos_flag = 0, viz_surf_states_flag = 0;	/* flags */
  char file_name[1024];
  char *ch_ptr = NULL; /* pointer used to extract data file name */
  char *grid_mol_name = NULL; /* points to the name of the grid molecule */
  char mesh_pos_name[1024]; /* meshes vertices data file name */
  char mesh_states_name[1024]; /* meshes states data file name */
  char region_viz_data_name[1024]; /* region_viz_data file name */
  char mol_pos_name[1024]; /* molecule positions data file name */
  char mol_states_name[1024]; /* molecule states data file name */
  char mol_orient_name[1024]; /* molecule orientations data file name */
  char frame_numbers_name[1024]; /* frame numbers data file name */ 
  char time_values_name[1024]; /* time values data file name */
  char buffer[100]; /* used to write 'frame_data' object information */
  /* used to write combined group information */
  static u_int member_meshes_iteration = UINT_MAX;
  static u_int member_molecules_iteration = UINT_MAX;
  static u_int member_effectors_iteration = UINT_MAX;
  /* linked list that stores data for the 'frame_data' object */
  static struct infinite_string_array frame_data_series_list;
  static u_int frame_data_series_count = 0; /* count elements in frame_data_series_list array.*/
  /* linked lists that stores data for the 'frame_numbers' object */
  static struct infinite_uint_array frame_numbers_meshes;
  static struct infinite_uint_array frame_numbers_vol_mols;
  static struct infinite_uint_array frame_numbers_surf_mols;
  static struct infinite_uint_array time_values;

  static u_int frame_numbers_meshes_count = 0; /* count elements in 
                                           frame_numbers_meshes array  */
  static u_int frame_numbers_vol_mols_count = 0; /* count elements in 
                                           frame_numbers_vol_mols array  */
  static u_int frame_numbers_surf_mols_count = 0; /* count elements in 
                                           frame_numbers_surf_mols array  */
  static u_int time_values_count = 0; /* count elements in 
                                           time_values array  */
  char my_byte_order[8];  /* shows binary ordering ('lsb' or 'msb') */
  static int mesh_pos_byte_offset = 0;  /* defines position of the object data
                                  in the mesh positions binary data file */
  static int mesh_states_byte_offset = 0; /* defines position of the object data
                                    in the mesh states binary file */
  static int region_data_byte_offset = 0; /* defines position of the 
                            object data in the region_data binary file */
  static int region_data_byte_offset_prev = 0; /* defines position of the 
                            object data in the region_data binary file */
  static int mol_pos_byte_offset = 0; /*defines position of the object data 
                           in the molecule positions binary file */
  static int mol_orient_byte_offset = 0; /*defines position of the object data 
                           in the molecule orientations binary file */
  static int mol_states_byte_offset = 0; /* defines position of the object data
                              in the molecule states binary file. */
  static int frame_numbers_byte_offset = 0; /* defines position of the frame 
                          numbers data in the frame_numbers binary file */
  static int time_values_byte_offset = 0; /* defines position of the time 
                          values data in the time_values binary file */
  int mol_pos_byte_offset_prev = 0; /* used when defining position  
                           in the molecule positions binary file */
  int mol_orient_byte_offset_prev = 0; /* used when defining position                               in the molecule orientations binary file */
  int mol_states_byte_offset_prev = 0; /* used when defining position                                  in the molecule states binary file. */
  /* placeholders for the members of the combined group */
  static long long mesh_group_index = 0;
  static long long mol_group_index = 0;
  static long long eff_group_index = 0;
  int combined_group_index = 0;

  /* counts number of times master_header file was opened.*/
  static int count_master_header = 0;
  /* counts number of times mol_pos_data file was opened.*/
  static int count_mol_pos_data = 0;
  /* counts number of times mol_orient_data file was opened. */
  static int count_mol_orient_data = 0;
  /* counts number of times mol_states_data file was opened. */
  static int count_mol_states_data = 0;
  /* counts number of times mesh_pos_data file was opened.*/
  static int count_mesh_pos_data = 0;
  /* counts number of times mesh_states_data file was opened. */
  static int count_mesh_states_data = 0;
  /* counts number of times region_values_data file was opened. */
  static int count_region_data = 0;

  /* points to the values of the current iteration steps
     for certain frame types. */
  static long long curr_surf_pos_iteration_step = -1;
  static long long curr_region_data_iteration_step = -1;
  static long long curr_mol_pos_iteration_step = -1;
  static long long curr_mol_orient_iteration_step = -1;

  /* points to the special case in iteration steps
     when both values for GEOMETRY and REGION_VIZ_VALUES
     are equal.  E.g. for the cases when GEOMETRY = [0,200], 
     REGION_VIZ_VALUES = [0,100, 200,300] 
     special_surf_iteration_step = [0,200].Same for (MOL_POS,MOL_ORIENT). 
  */
   
  static long long special_surf_iteration_step = -1;
   
  static long long special_mol_iteration_step = -1;
  
  /* counts number of this function executions
     for special_surf_iteration_step/special_mol_iteration_step cases. */
  
  static int special_surf_frames_counter = 0;
  
  static int special_mol_frames_counter = 0;

  struct frame_data_list * fdl_ptr; 
  u_int n_species = world->n_species;
  species_list=world->species_list;

  log_file=world->log_file;
  no_printf("Viz output in DREAMM_V3 mode...\n");
  
  if(world->file_prefix_name == NULL) {
   	fprintf(world->err_file, "File %s, Line %ld: Inside VIZ_DATA_OUTPUT block the required keyword FILENAME_PREFIX is missing.\n", __FILE__, (long)__LINE__);
   	exit(1);
  }

  word_p=(unsigned char *)&word;
  word=0x04030201;

  if (word_p[0]==1) {
    sprintf(my_byte_order,"lsb");
  }
  else {
    sprintf(my_byte_order,"msb");
  }
  
  /*initialize infinite arrays. */
  ia_init(&frame_numbers_meshes);
  ia_init(&frame_numbers_vol_mols);
  ia_init(&frame_numbers_surf_mols);
  ia_init(&time_values);
  ia_init(&frame_data_series_list);

  viz_iteration = (u_int)(fdlp->viz_iterationll);
  n_viz_iterations = (u_int)(fdlp->n_viz_iterations);

  viz_type=fdlp->type;

  /* here is the check to prevent writing twice the same info 
    - at the end of one checkpoint and the beginning of the next checkpoint.
   */
  if((world->chkpt_flag) && (world->start_time > 0)){
     if (world->it_time == world->start_time){
         if(viz_iteration % (world->chkpt_iterations) == 0){
             return 0;
         }
     }
  }

  /* initialize flags */
  viz_mol_pos_flag = (viz_type==MOL_POS);
  viz_mol_orient_flag = (viz_type==MOL_ORIENT);
  viz_mol_all_data_flag = (viz_type==ALL_MOL_DATA);
  if(viz_mol_all_data_flag){
       viz_mol_pos_flag = viz_mol_orient_flag = 1;
  }
  if(viz_mol_pos_flag){
     if((world->viz_output_flag & VIZ_MOLECULES_STATES) != 0){
     	viz_mol_states_flag = 1;
     }
  } 
  
  viz_surf_pos_flag = (viz_type==MESH_GEOMETRY);
  viz_region_data_flag = (viz_type==REG_DATA);
  viz_surf_all_data_flag = (viz_type==ALL_MESH_DATA);
  if(viz_surf_all_data_flag){
     viz_surf_pos_flag = viz_region_data_flag = 1;
  }
  if(viz_surf_pos_flag){
     if((world->viz_output_flag & VIZ_SURFACE_STATES) != 0){
     	viz_surf_states_flag = 1;
     }
  }


  /* these arrays are used to hold values of main_index for objects */
  if((viz_surf_pos_flag) || (viz_region_data_flag))
  {
     if(obj_to_show_number > 0)
     {
  	if((surf_states == NULL) && (viz_surf_states_flag)){ 
     		if ((surf_states=(u_int *)malloc(obj_to_show_number*sizeof(u_int)))==NULL)      {
                        fprintf(world->err_file, "File %s, Line %ld: memory allocation error.\n", __FILE__, (long)__LINE__);
         		return (1);
     		}
        	for(ii = 0; ii < obj_to_show_number; ii++){
			surf_states[ii] = 0;
     		}
   	}
   	if((surf_pos == NULL) && (viz_surf_pos_flag)){
     		if ((surf_pos=(u_int *)malloc(obj_to_show_number*sizeof(u_int)))==NULL) {
                        fprintf(world->err_file, "File %s, Line %ld: memory allocation error.\n", __FILE__, (long)__LINE__);
         		return (1);
     		}
     		for(ii = 0; ii < obj_to_show_number; ii++){
			surf_pos[ii] = 0;
     		}
   	}
   	if((surf_con == NULL) && (viz_surf_pos_flag)){
      		if ((surf_con=(u_int *)malloc(obj_to_show_number*sizeof(u_int)))==NULL) {
                        fprintf(world->err_file, "File %s, Line %ld: memory allocation error.\n", __FILE__, (long)__LINE__);
         		return (1);
      		}
      		for(ii = 0; ii < obj_to_show_number; ii++){
			surf_con[ii] = 0;
      		}
        }
   	if(surf_field_indices == NULL){
      		if ((surf_field_indices=(u_int *)malloc(obj_to_show_number*sizeof(u_int)))==NULL) {
                        fprintf(world->err_file, "File %s, Line %ld: memory allocation error.\n", __FILE__, (long)__LINE__);
         		return (1);
      		}
      		for(ii = 0; ii < obj_to_show_number; ii++){
			surf_field_indices[ii] = 0;
      		}
        }
    

   	/* initialize array of viz_objects names */
   	if(obj_names == NULL){
      		if ((obj_names = (char **)malloc(obj_to_show_number*sizeof(char *)))==NULL) {
                        fprintf(world->err_file, "File %s, Line %ld: memory allocation error.\n", __FILE__, (long)__LINE__);
         		return (1);
      		}
      		for(ii = 0; ii < obj_to_show_number; ii++){
			obj_names[ii] = NULL;
      		}
                /* create an array of region names */
   	        if(region_names == NULL){
      	            if ((region_names = (char ***)malloc(obj_to_show_number*sizeof(char **)))==NULL) { 
                        fprintf(world->err_file, "File %s, Line %ld: memory allocation error.\n", __FILE__, (long)__LINE__);
         	        return (1);
      	            }
                    /* initialize it */
                    for(ii = 0; ii < obj_to_show_number; ii++)
                    {
                       region_names[ii] = NULL;
                    }
                 }          
      
     	        ii = 0;
     	        vizp = world->viz_obj_head;
     	        while(vizp != NULL){
		   vcp = vizp->viz_child_head;
        	   while(vcp != NULL){
         	       obj_names[ii] = my_strdup(vcp->obj->sym->name);
                       if(obj_names[ii] == NULL){
                            fprintf(world->err_file, "File %s, Line %ld: memory allocation error.\n", __FILE__, (long)__LINE__);
         		    return (1);
         	        }
                
                        int n_regs;   /* number of regions in the object */
                        /* subtract the default region ALL */
                        n_regs = vcp->obj->num_regions - 1;
                        if(n_regs > 0){
		            region_names[ii] = (char **)malloc(n_regs*sizeof(char *));
                            if(region_names[ii] == NULL){ 
                               fprintf(world->err_file, "File %s, Line %ld: memory allocation error.\n", __FILE__, (long)__LINE__);
         		       return (1);
      		            }
                        }

                        jj = 0;
 
                        for(rlp = vcp->obj->regions; rlp != NULL; rlp = rlp->next){
                            rp = rlp->reg;
                            if(strcmp(rp->region_last_name, "ALL") == 0) continue;
                            
                            region_names[ii][jj] = my_strdup(rp->region_last_name);
                            if(region_names[ii][jj] == NULL)
                            { 
                               fprintf(world->err_file, "File %s, Line %ld: memory allocation error.\n", __FILE__, (long)__LINE__);
         		       return (1);
      		            }
                            jj++;
                             
                        }
                        
                        ii++;
         		vcp = vcp->next;
       		    } /* end while (vcp) */
       		    vizp = vizp->next;
    	         } /* end while (vizp) */

   	}
       
       } /* end if (obj_to_show_number > 0) */
    }  /* end if (viz_surf_pos_flag  || viz_region_data_flag) */



     /* these arrays are used to hold values of main_index for effectors 
      and 3D molecules */
   if(viz_mol_pos_flag || viz_mol_orient_flag)
   {
     if(eff_to_show_number > 0)
     {
     	if((eff_states == NULL) && (viz_mol_states_flag)){ 
     		if((eff_states=(u_int *)malloc(eff_to_show_number*sizeof(u_int)))==NULL)        {
                        fprintf(world->err_file, "File %s, Line %ld: memory allocation error.\n", __FILE__, (long)__LINE__);
         		return (1);
     		}
     		for(ii = 0; ii < eff_to_show_number; ii++){
			eff_states[ii] = 0;
     		}
     	}
     	if((eff_pos == NULL) && (viz_mol_pos_flag)){
     		if ((eff_pos=(u_int *)malloc(eff_to_show_number*sizeof(u_int)))==NULL) 	{
                        fprintf(world->err_file, "File %s, Line %ld: memory allocation error.\n", __FILE__, (long)__LINE__);
         		return (1);
     		}
     		for(ii = 0; ii < eff_to_show_number; ii++){
			eff_pos[ii] = 0;
     		}
     	}
     	if((eff_orient == NULL) && (viz_mol_orient_flag)){
     		if ((eff_orient = (u_int *)malloc(eff_to_show_number*sizeof(u_int)))==NULL) 	{
                        fprintf(world->err_file, "File %s, Line %ld: memory allocation error.\n", __FILE__, (long)__LINE__);
         		return (1);
     		}
     		for(ii = 0; ii < eff_to_show_number; ii++){
			eff_orient[ii] = 0;
     		}
     	}
     	if(eff_field_indices == NULL){
     		if ((eff_field_indices = (u_int *)malloc(eff_to_show_number*sizeof(u_int)))==NULL) 	{
                        fprintf(world->err_file, "File %s, Line %ld: memory allocation error.\n", __FILE__, (long)__LINE__);
         		return (1);
     		}
     		for(ii = 0; ii < eff_to_show_number; ii++){
			eff_field_indices[ii] = 0;
     		}
     	}
   	/* initialize array of grid_mol's names */
   	if(eff_names == NULL){

      		if ((eff_names = (char **)malloc(eff_to_show_number*sizeof(char *)))==NULL) {
                        fprintf(world->err_file, "File %s, Line %ld: memory allocation error.\n", __FILE__, (long)__LINE__);
         		return (1);
      		}
      		for(ii = 0; ii < eff_to_show_number; ii++){
			eff_names[ii] = NULL;
      		}
   	
                index = 0;
        
                for(ii = 0; ii < world->n_species; ii++)
                {
     	           specp = world->species_list[ii];
     	           if((specp->flags & IS_SURFACE) != 0) continue;
     	           if(strcmp(specp->sym->name, "GENERIC_MOLECULE") == 0) continue;
                   if(((specp->flags & ON_GRID) == ON_GRID) && (specp->viz_state > 0)){ 
                      eff_names[index] = my_strdup(specp->sym->name);
                      if(eff_names[index] == NULL){
                         fprintf(world->err_file, "File %s, Line %ld: memory allocation error.\n", __FILE__, (long)__LINE__);
                         return (1);
                      }
                      index++;
                    } 
                 }
        }
        index = 0;
     } /* end if (eff_to_show_number > 0) */
   }  /* end if(viz_eff_pos_flag  */
   
   if(viz_mol_pos_flag || viz_mol_orient_flag)
   {
     if(mol_to_show_number > 0)
     {
     	if((mol_states == NULL) && (viz_mol_states_flag)){ 
     		if((mol_states=(u_int *)malloc(mol_to_show_number*sizeof(u_int)))==NULL)        {
                         fprintf(world->err_file, "File %s, Line %ld: memory allocation error.\n", __FILE__, (long)__LINE__);
         		return (1);
     		}
     		for(ii = 0; ii < mol_to_show_number; ii++){
			mol_states[ii] = 0;
     		}
     	}
     	if((mol_pos == NULL) && (viz_mol_pos_flag)) {
     		if ((mol_pos=(u_int *)malloc(mol_to_show_number*sizeof(u_int)))==NULL) 	{
                         fprintf(world->err_file, "File %s, Line %ld: memory allocation error.\n", __FILE__, (long)__LINE__);
         		return (1);
     		}
     		for(ii = 0; ii < mol_to_show_number; ii++){
			mol_pos[ii] = 0;
     		}
     	}
     	if((mol_orient == NULL) && (viz_mol_orient_flag)){
     		if ((mol_orient = (u_int *)malloc(mol_to_show_number*sizeof(u_int)))==NULL) 	{
                         fprintf(world->err_file, "File %s, Line %ld: memory allocation error.\n", __FILE__, (long)__LINE__);
         		return (1);
     		}
     		for(ii = 0; ii < mol_to_show_number; ii++){
			mol_orient[ii] = 0;
     		}
     	}
     	if(mol_field_indices == NULL){
     		if ((mol_field_indices = (u_int *)malloc(mol_to_show_number*sizeof(u_int)))==NULL) 	{
                         fprintf(world->err_file, "File %s, Line %ld: memory allocation error.\n", __FILE__, (long)__LINE__);
         		return (1);
     		}
     		for(ii = 0; ii < mol_to_show_number; ii++){
			mol_field_indices[ii] = 0;
     		}
     	}
        
   	/* initialize array of mol's names */
   	if(mol_names == NULL){
      		if ((mol_names = (char **)malloc(mol_to_show_number*sizeof(char *)))==NULL) {
                         fprintf(world->err_file, "File %s, Line %ld: memory allocation error.\n", __FILE__, (long)__LINE__);
         		return (1);
      		}
      		for(ii = 0; ii < mol_to_show_number; ii++){
			mol_names[ii] = NULL;
      		}
   	
        	index = 0;
        	for(ii = 0; ii < world->n_species; ii++)
        	{
     	   	   specp = world->species_list[ii];
     	   	   if((specp->flags & IS_SURFACE) != 0) continue;
     	   	   if(strcmp(specp->sym->name, "GENERIC_MOLECULE") == 0) continue;
           	   if(((specp->flags & NOT_FREE) == 0) && (specp->viz_state > 0)){ 
                       mol_names[index] = my_strdup(specp->sym->name);
                       if(mol_names[index] == NULL){
                         fprintf(world->err_file, "File %s, Line %ld: memory allocation error.\n", __FILE__, (long)__LINE__);
                         return (1);
                       }
                       index++;
                   } 
                }
        }
        index = 0;

     } /* end if (mol_to_show_number > 0) */
  }  /* end if(viz_mol_pos_flag || viz_mol_orient_flag) */

  /* Open master header file. */
  if(world->chkpt_flag){
    sprintf(file_name,"%s.%u.dx",world->file_prefix_name, world->chkpt_seq_num);
  }else{
    sprintf(file_name,"%s.dx",world->file_prefix_name);
  }

  if(count_master_header == 0){
      if ((master_header=fopen(file_name,"w"))==NULL) {
           fprintf(world->err_file, "File %s, Line %ld: cannot open master header file %s.\n", __FILE__, (long)__LINE__,file_name);
           return(1);
      }
      count_master_header++;

  }else{
      if ((master_header=fopen(file_name,"a"))==NULL) {
           fprintf(world->err_file, "File %s, Line %ld: cannot open master header file %s.\n", __FILE__, (long)__LINE__,file_name);
           return(1);
      }
      count_master_header++;

  }

  if (viz_mol_pos_flag) {
     if(world->chkpt_flag){
       sprintf(file_name,"%s.molecule_positions.%u.bin",world->file_prefix_name,               world->chkpt_seq_num);
     }else{
       sprintf(file_name,"%s.molecule_positions.bin",world->file_prefix_name);
     }
     /* remove the folder name from the molecule_positions data file name */
     ch_ptr = strrchr(file_name, '/');
     ++ch_ptr;
     strcpy(mol_pos_name, ch_ptr);
     
     if (count_mol_pos_data == 0){
        if ((mol_pos_data=fopen(file_name,"wb"))==NULL) {
           fprintf(world->err_file, "File %s, Line %ld: cannot open file %s.\n", __FILE__, (long)__LINE__,file_name);
           return(1);
        }else{}
        count_mol_pos_data++;
     }else{
        if ((mol_pos_data=fopen(file_name,"ab"))==NULL) {
              fprintf(world->err_file, "File %s, Line %ld: cannot open file %s.\n", __FILE__, (long)__LINE__,file_name);
              return(1);
        }
        count_mol_pos_data++;
    }
  }

  if(viz_mol_orient_flag)
  {
    if(world->chkpt_flag){
      sprintf(file_name,"%s.molecule_orientations.%u.bin",world->file_prefix_name, world->chkpt_seq_num);
    }else{
      sprintf(file_name,"%s.molecule_orientations.bin",world->file_prefix_name);
    }
     /* remove the folder name from the molecule_positions data file name */
     ch_ptr = strrchr(file_name, '/');
     ++ch_ptr;
     strcpy(mol_orient_name, ch_ptr);

     if (count_mol_orient_data == 0){
        if ((mol_orient_data=fopen(file_name,"wb"))==NULL) {
           fprintf(world->err_file, "File %s, Line %ld: cannot open file %s.\n", __FILE__, (long)__LINE__,file_name);
           return(1);
        }
        count_mol_orient_data++;
     }else{
        if ((mol_orient_data=fopen(file_name,"ab"))==NULL) {
              fprintf(world->err_file, "File %s, Line %ld: cannot open file %s.\n", __FILE__, (long)__LINE__,file_name);
              return(1);
        }
        count_mol_orient_data++;

     }
   
  }


    if (viz_mol_states_flag) {
      if(world->chkpt_flag){
         sprintf(file_name,"%s.molecule_states.%u.bin",world->file_prefix_name, world->chkpt_seq_num);
      }else{
         sprintf(file_name,"%s.molecule_states.bin",world->file_prefix_name);
      }
       /* remove the folder name from the molecule_states data file name */
       ch_ptr = strrchr(file_name, '/');
       ++ch_ptr;
       strcpy(mol_states_name, ch_ptr);
     
       if (count_mol_states_data == 0){
            if ((mol_states_data = fopen(file_name,"wb"))==NULL) {
                   fprintf(world->err_file, "File %s, Line %ld: cannot open file %s.\n", __FILE__, (long)__LINE__,file_name);
           	   return(1);
            }
            count_mol_states_data++;
       }else{
           if ((mol_states_data = fopen(file_name,"ab"))==NULL) {
                   fprintf(world->err_file, "File %s, Line %ld: cannot open file %s.\n", __FILE__, (long)__LINE__,file_name);
           	   return(1);
            }
            count_mol_states_data++;
       }
    }
    
      if (viz_surf_pos_flag) {
        if(world->chkpt_flag){
      	   sprintf(file_name,"%s.mesh_positions.%u.bin",world->file_prefix_name, world->chkpt_seq_num);
        }else{
      	   sprintf(file_name,"%s.mesh_positions.bin",world->file_prefix_name);
        }
     
     	/* remove the folder name from the mesh_positions data file name */
     	ch_ptr = strrchr(file_name, '/');
     	++ch_ptr;
     	strcpy(mesh_pos_name, ch_ptr);

       if (count_mesh_pos_data == 0){
      	 if ((mesh_pos_data=fopen(file_name,"wb"))==NULL) {
                fprintf(world->err_file, "File %s, Line %ld: cannot open file %s.\n", __FILE__, (long)__LINE__,file_name);
                return(1);
          }
          count_mesh_pos_data++;
       }else{
      	   if ((mesh_pos_data=fopen(file_name,"ab"))==NULL) {
               fprintf(world->err_file, "File %s, Line %ld: cannot open file %s.\n", __FILE__, (long)__LINE__,file_name);
               return(1);
           }
           count_mesh_pos_data++;
       }

      }


      if (viz_surf_states_flag) {
         if(world->chkpt_flag){
           sprintf(file_name,"%s.mesh_states.%u.bin", world->file_prefix_name, world->chkpt_seq_num);
         }else{
           sprintf(file_name,"%s.mesh_states.bin", world->file_prefix_name);
         }
     
        /* remove the folder name from the mesh_states data file name */
        ch_ptr = strrchr(file_name, '/');
        ++ch_ptr;
        strcpy(mesh_states_name, ch_ptr);

       if (count_mesh_states_data == 0){
           if ((mesh_states_data=fopen(file_name,"wb"))==NULL) {
                  fprintf(world->err_file, "File %s, Line %ld: cannot open file %s.\n", __FILE__, (long)__LINE__,file_name);
                  return(1);
            }
            count_mesh_states_data++;
       }else{
          if ((mesh_states_data=fopen(file_name,"ab"))==NULL) {
              fprintf(world->err_file, "File %s, Line %ld: cannot open file %s.\n", __FILE__, (long)__LINE__,file_name);
              return(1);
          }
          count_mesh_states_data++;
       }
      }

      if (viz_region_data_flag) {
         if(world->chkpt_flag){
            sprintf(file_name,"%s.region_indices.%u.bin", world->file_prefix_name, world->chkpt_seq_num);
         }else{
            sprintf(file_name,"%s.region_indices.bin", world->file_prefix_name);
         }
     
        /* remove the folder name from the region values data file name */
        ch_ptr = strrchr(file_name, '/');
        ++ch_ptr;
        strcpy(region_viz_data_name, ch_ptr);

       if (count_region_data == 0){
          if ((region_data=fopen(file_name,"wb"))==NULL) {
              fprintf(world->err_file, "File %s, Line %ld: cannot open file %s.\n", __FILE__, (long)__LINE__,file_name);
              return(1);
           }
           count_region_data++;
       }else{
          if ((region_data=fopen(file_name,"ab"))==NULL) {
                fprintf(world->err_file, "File %s, Line %ld: cannot open file %s.\n", __FILE__, (long)__LINE__,file_name);
                return(1);
          }
          count_region_data++;
       }
      }


    /* find out the values of the current iteration steps for
       (GEOMETRY, REG_DATA),  
       (MOL_POS, MOL_ORIENT) frames combinations */
   
    fdl_ptr = world->frame_data_head;
    while(fdl_ptr != NULL){
        if(fdl_ptr->type == MESH_GEOMETRY){
             curr_surf_pos_iteration_step = fdl_ptr->viz_iterationll;
        }if(fdl_ptr->type == REG_DATA){
             curr_region_data_iteration_step = fdl_ptr->viz_iterationll;
        }if(fdl_ptr->type == ALL_MESH_DATA){
             curr_surf_pos_iteration_step = fdl_ptr->viz_iterationll;
             curr_region_data_iteration_step = fdl_ptr->viz_iterationll;
        }else if(fdl_ptr->type == MOL_POS){
             curr_mol_pos_iteration_step = fdl_ptr->viz_iterationll;
	}else if(fdl_ptr->type == MOL_ORIENT) 
        {
             curr_mol_orient_iteration_step = fdl_ptr->viz_iterationll;
        }else if(fdl_ptr->type == ALL_MOL_DATA){
             curr_mol_pos_iteration_step = fdl_ptr->viz_iterationll;
             curr_mol_orient_iteration_step = fdl_ptr->viz_iterationll;
        }	
        fdl_ptr = fdl_ptr->next;
    }
  
    /* If the values of the current iteration steps for REG_VIZ_VALUES and
       MESH_GEOMETRY are equal set this value to the 
       special_surf_iteration_step.
       Do the same for the (MOL_POS,MOL_ORIENT).
    */
     
    if(curr_region_data_iteration_step == curr_surf_pos_iteration_step){
	special_surf_iteration_step = curr_surf_pos_iteration_step;
    }
    
    if(curr_mol_orient_iteration_step == curr_mol_pos_iteration_step){
	special_mol_iteration_step = curr_mol_orient_iteration_step;
    }
    
    /* check for the special_iteration_step  */  
    if(viz_surf_pos_flag || viz_region_data_flag){	    
	    if(fdlp->viz_iterationll == special_surf_iteration_step){
		special_surf_frames_counter++;
            }
    }


/* dump walls */
  if((viz_surf_pos_flag) || (viz_region_data_flag)) {
     vizp = world->viz_obj_head;
     
     while(vizp!=NULL) {

    /* Traverse all visualized compartments 
       output mesh element positions and connections */
    vcp = vizp->viz_child_head;

    while(vcp!=NULL) {
      objp = vcp->obj;
      if(objp->viz_state == NULL) continue;
      pop=(struct polygon_object *)objp->contents;

      if (objp->object_type==POLY_OBJ || objp->object_type==BOX_OBJ) {
          opp=(struct ordered_poly *)pop->polygon_data;
          edp=opp->element;
          element_data_count=objp->n_walls_actual;

          if (viz_surf_pos_flag) {
            fprintf(master_header,
              "object %d class array type float rank 1 shape 3 items %d %s binary data file %s,%d # %s.positions #\n",
              main_index,objp->n_verts,my_byte_order,mesh_pos_name, mesh_pos_byte_offset, objp->sym->name);
            fprintf(master_header,
              "\tattribute \"dep\" string \"positions\"\n\n");
            surf_pos[surf_pos_index] = main_index;
            surf_pos_index++;
            main_index++;

            /* output polyhedron vertices */
            for (ii=0;ii<objp->n_verts;ii++) {
              v1 = (float)(world->length_unit*objp->verts[ii].x);
              v2 = (float)(world->length_unit*objp->verts[ii].y);
              v3 = (float)(world->length_unit*objp->verts[ii].z);
              fwrite(&v1,sizeof v1,1,mesh_pos_data);
              fwrite(&v2,sizeof v2,1,mesh_pos_data);
              fwrite(&v3,sizeof v3,1,mesh_pos_data);
              mesh_pos_byte_offset += (sizeof(v1) + sizeof(v2) + sizeof(v3));
            }

            /* output polygon element connections */
            fprintf(master_header,
              "object %d class array type int rank 1 shape 3 items %d %s binary data file %s,%d # %s.connections #\n",
              main_index,element_data_count,my_byte_order, mesh_pos_name, mesh_pos_byte_offset, objp->sym->name);
            fprintf(master_header,
              "\tattribute \"ref\" string \"positions\"\n");
            fprintf(master_header,
              "\tattribute \"element type\" string \"triangles\"\n\n");

            for (ii=0;ii<objp->n_walls;ii++) {
              if (!get_bit(pop->side_removed,ii)) {
                for (jj=0;jj<edp[ii].n_verts-2;jj++) {
                  vi1=edp[ii].vertex_index[0];
                  vi2=edp[ii].vertex_index[jj+1];
                  vi3=edp[ii].vertex_index[jj+2];
	          fwrite(&vi1,sizeof vi1,1,mesh_pos_data);
	          fwrite(&vi2,sizeof vi2,1,mesh_pos_data);
	          fwrite(&vi3,sizeof vi3,1,mesh_pos_data);
                  mesh_pos_byte_offset += (sizeof(vi1) + sizeof(vi2) + sizeof(vi3));
                }
              }
            }
            surf_con[surf_con_index] = main_index;
            main_index++;
            surf_con_index++;

          } /* end viz_surf_pos_flag for POLY_OBJ */
 
          if (viz_surf_states_flag) {
            fprintf(master_header,
              "object %d class array type int rank 0 items %d %s binary data file %s,%d # %s.states #\n", main_index, element_data_count, my_byte_order, mesh_states_name, mesh_states_byte_offset, objp->sym->name);
            fprintf(master_header,"\tattribute \"dep\" string \"connections\"\n\n");
           surf_states[surf_states_index] = main_index;
           surf_states_index++;
           main_index++;

            for (ii=0;ii<objp->n_walls;ii++) {
               if (!get_bit(pop->side_removed,ii)) {
                 state=objp->viz_state[ii];
                 fwrite(&state,sizeof (state),1,mesh_states_data);
                 mesh_states_byte_offset += sizeof(state);
               }
            }

          } 	/* end viz_surf_states_flag for POLY_OBJ */

          if(viz_region_data_flag && (objp->num_regions > 1))
          {
              surf_region_values_index = 0;

              for(rlp = objp->regions; rlp != NULL; rlp = rlp->next)
              {
                  rp = rlp->reg;
                  if(strcmp(rp->region_last_name, "ALL") == 0) continue; 
                  
                  /* number of walls in the region */
                  int region_walls_number = 0; 
                  
                  for(jj = 0; jj < objp->n_walls; jj++)
                  {
                    int n = objp->wall_p[jj]->side;
                    if(get_bit(rp->membership,n))
                    {
                        fwrite(&n,sizeof (n),1,region_data);
                        region_data_byte_offset += sizeof(n);
                        region_walls_number++;
                    }

                  }

                  fprintf(master_header,
                        "object %d class array type int rank 0 items %d %s binary data file %s,%d # %s.region_data #\n", main_index, region_walls_number, my_byte_order, region_viz_data_name, region_data_byte_offset_prev, objp->sym->name);
                  fprintf(master_header,"\tattribute \"ref\" string \"connections\"\n");  
                  fprintf(master_header,"\tattribute \"identity\" string \"region_indices\"\n");
                  fprintf(master_header,"\tattribute \"name\" string \"%s\"\n", rp->region_last_name);
                  if(rp->region_viz_value > 0){
                      fprintf(master_header,"\tattribute \"viz_value\" value %d\n", rp->region_viz_value);
                  }

                  region_data_byte_offset_prev = region_data_byte_offset;
                  surf_region_values[surf_obj_region_values_index][surf_region_values_index] = main_index;
                  surf_region_values_index++;
                  main_index++;
                  fprintf(master_header, "\n\n");

              } /* end for */
              

           } /* end if(region_data_flag) */


       }	/* end POLY_OBJ */

       surf_obj_region_values_index++; 

      vcp = vcp->next;
      }
      
      vizp = vizp->next;
      }
    } /* end (viz_surf_pos_flag || viz_surf_states_flag || viz_region_data_flag) for vizp */
    
    /* build fields here */
      int show_meshes = 0;
   if(obj_to_show_number > 0)
   {
      if(fdlp->type == ALL_MESH_DATA) {
           show_meshes = 1;
      }
      else if((viz_surf_pos_flag || viz_region_data_flag) && ((special_surf_frames_counter == 0) || (special_surf_frames_counter == 2))){
		show_meshes = 1;
      }
    }


    if(show_meshes)
    {
         for(ii = 0; ii < obj_to_show_number; ii++)
         {
             if(obj_names[ii] != NULL){
                fprintf(master_header,
                    "object %d field   # %s #\n",main_index, obj_names[ii]);
             }
             if(surf_pos[ii] > 0){
             	fprintf(master_header,
                 "\tcomponent \"positions\" value %d\n",surf_pos[ii]);
             }
             if(surf_con[ii] > 0){
             	fprintf(master_header,
                 "\tcomponent \"connections\" value %d\n",surf_con[ii]);
             }
             if(surf_states != NULL){
                if(surf_states[ii] > 0){
             	   fprintf(master_header,
                	"\tcomponent \"state_values\" value %d\n",surf_states[ii]);
	        }
             }
             if(surf_region_values[ii] != NULL)
             {
                for(jj = 0; jj < obj_num_regions[ii]; jj++)
                {
                   if(surf_region_values[ii][jj] > 0){
             	      fprintf(master_header,
                	"\tcomponent \"%s\" value %d\n",region_names[ii][jj], surf_region_values[ii][jj]);
	           }
                }
             }
             fprintf(master_header, "\n");
             surf_field_indices[ii] = main_index;

             main_index++;
           }     
     }
    
    /* create a group object for all meshes. */
    if(show_meshes)
    {
        vizp = world->viz_obj_head;
        if(vizp != NULL){

        	fprintf(master_header,"object \"%s%u\" group # %s #\n","meshes_", viz_iteration, "meshes");
        	mesh_group_index = main_index;
        	main_index++;
                member_meshes_iteration = viz_iteration;

                ii = 0;
        	while(vizp != NULL) {
         		vcp = vizp->viz_child_head;
         		while(vcp != NULL){
                          if(surf_field_indices[ii] > 0){
             			fprintf(master_header,"\tmember \"%s\" value %d\n",vcp->obj->sym->name,surf_field_indices[ii]);
                          }
                          ii++;
		       	  vcp = vcp->next;
         		}
         		vizp = vizp->next;
        	}       
        	fprintf(master_header, "\n"); 
      	}  /* end (if vizp) */

        /* store iteration_number for meshes */
        ia_uint_store(&frame_numbers_meshes, frame_numbers_meshes_count, viz_iteration);
        frame_numbers_meshes_count++;
  
   } 

      
    /* Visualize molecules. */



/* dump grid molecules. */
    /* check for the special_iteration_step  */ 
    if(viz_mol_pos_flag || viz_mol_orient_flag){	    
	    if(fdlp->viz_iterationll == special_mol_iteration_step){
		special_mol_frames_counter++;
            }
    }


 if(viz_mol_pos_flag || viz_mol_orient_flag){	    

       /* create references to the numbers of grid molecules of each name. */
       if ((viz_grid_mol_count=(u_int *)malloc(n_species*sizeof(u_int)))==NULL) {
         return(1);
       }

       /* perform initialization */
      for (ii = 0; ii < n_species; ii++){
         spec_id = species_list[ii]->species_id;
         viz_grid_mol_count[spec_id]=0;

      } 


   for (ii = 0; ii < n_species; ii++)
   {
     specp = species_list[ii];
     if((specp->flags & ON_GRID) == 0) continue; 
     if(specp->viz_state == EXCLUDE_OBJ) continue; 
 
        grid_mol_name = specp->sym->name;
        spec_id = specp->species_id; 

        if(viz_mol_pos_flag){
      	   mol_pos_byte_offset_prev = mol_pos_byte_offset;
        }
        if(viz_mol_orient_flag){
      	   mol_orient_byte_offset_prev = mol_orient_byte_offset;
        }
        if(viz_mol_states_flag){
      	   mol_states_byte_offset_prev = mol_states_byte_offset;
        }
        /* Write binary data files. */
        vizp = world->viz_obj_head;
        while(vizp != NULL) {
         vcp = vizp->viz_child_head;
         while(vcp != NULL){
	     objp = vcp->obj;
             wp = objp->wall_p;
             no_printf("Traversing walls in object %s\n",objp->sym->name);
             for (jj=0;jj<objp->n_walls;jj++) {
                w = wp[jj];
                if (w!=NULL) {
	           sg = w->effectors;
                   if (sg!=NULL) {
                     for (index=0;index<sg->n_tiles;index++) {
                        grid2xyz(sg,index,&p0);
	                gmol=sg->mol[index];
	                if ((gmol!=NULL) && (viz_mol_states_flag)){
	                   state=sg->mol[index]->properties->viz_state;
                        }
                        if (gmol != NULL) {
                            if(spec_id == gmol->properties->species_id){
                                  viz_grid_mol_count[spec_id]++;

                               if(viz_mol_pos_flag){
                                   /* write positions information */
	           		   v1=(float)(world->length_unit*p0.x);
	           		   v2=(float)(world->length_unit*p0.y);
	           		   v3=(float)(world->length_unit*p0.z);
	           		   fwrite(&v1,sizeof v1,1,mol_pos_data);
	               		   fwrite(&v2,sizeof v2,1,mol_pos_data);
	           		   fwrite(&v3,sizeof v3,1,mol_pos_data);
                                   mol_pos_byte_offset += (sizeof(v1) + sizeof(v2) + sizeof(v3));
                               }
                               if(viz_mol_orient_flag){
                                   /* write orientations information */
                   		   v1=(float)((w->normal.x)*(gmol->orient));
                   		   v2=(float)((w->normal.y)*(gmol->orient));
                   		   v3=(float)((w->normal.z)*(gmol->orient));
                   		   fwrite(&v1,sizeof v1,1,mol_orient_data);
                   		   fwrite(&v2,sizeof v2,1,mol_orient_data);
                   		   fwrite(&v3,sizeof v3,1,mol_orient_data);
                                   mol_orient_byte_offset += (sizeof(v1) + sizeof(v2) + sizeof(v3));
                  		}
                          } /* end if strcmp */
                           
                        } /* end if (gmol)*/
                     } /* end for */
                   } /* end if (sg) */
                 } /* end if (w)*/
              } /* end for */
              vcp = vcp->next;
            } /* end while (vcp) */
                                              
            vizp = vizp->next;
         } /* end while (vzp) */

        num = viz_grid_mol_count[spec_id];
        
        if(viz_mol_pos_flag)
        {
           if(num > 0)
           {   
        	fprintf(master_header,"object %d class array type float rank 1 shape 3 items %d %s binary data file %s,%d # %s positions #\n",main_index,num,my_byte_order, mol_pos_name, mol_pos_byte_offset_prev, grid_mol_name);
        	fprintf(master_header,"\tattribute \"dep\" string \"positions\"\n\n");
            }else{
                /* output empty arrays for zero molecule counts here */
                fprintf(master_header,"object %d array   # %s positions #\n",main_index, grid_mol_name);
            }
            eff_pos[eff_pos_index] = main_index;
            eff_pos_index++;
            main_index++;
         }
         
         if(viz_mol_orient_flag)
         {
           if(num > 0)
           {
        	fprintf(master_header,"object %d class array type float rank 1 shape 3 items %d %s binary data file %s,%d   # %s orientations #\n",main_index,num,my_byte_order, mol_orient_name, mol_orient_byte_offset_prev, grid_mol_name);
        	fprintf(master_header,"\tattribute \"dep\" string \"positions\"\n\n");
            }else{
                /* output empty arrays for zero molecule counts here */
                fprintf(master_header,"object %d array   # %s orientations #\n",main_index, grid_mol_name);
            }
            eff_orient[eff_orient_index] = main_index;
            eff_orient_index++;
            main_index++;
        }

        if (viz_mol_states_flag) {
          if(num > 0)
          {
            /* write states information. */
            fwrite(&state,sizeof state,1,mol_states_data);
            mol_states_byte_offset += (sizeof state);
        
	    fprintf(master_header,"object %d class constantarray type int items %d %s binary data file %s,%d  # %s states #\n",main_index,num, my_byte_order,mol_states_name,mol_states_byte_offset_prev, grid_mol_name);

            fprintf(master_header,"\tattribute \"dep\" string \"positions\"\n\n");
	  }else{
             /* output empty arrays for zero molecule counts here */
             fprintf(master_header,"object %d array   # %s states #\n",main_index, grid_mol_name);
          }
          eff_states[eff_states_index] = main_index;
          eff_states_index++;
          main_index++;
       }
       if(num == 0){
         fprintf(master_header, "\n");
       }
   } /* end for loop */

 } /* end if(viz_mol_pos_flag || viz_mol_orient_flag) */

/* build fields for grid molecules here */
  int show_effectors = 0;

  if(eff_to_show_number > 0){
      if(fdlp->type == ALL_MOL_DATA) {
           show_effectors = 1;
      }
      else if((viz_mol_pos_flag || viz_mol_orient_flag) && ((special_mol_frames_counter ==0) || (special_mol_frames_counter == 2))){
		show_effectors = 1;
      }
  }
 
  if(show_effectors){


       for(ii = 0; ii < eff_to_show_number; ii++)
       {
             if(eff_names[ii] != NULL){
                fprintf(master_header,
                    "object %d field   # %s #\n",main_index, eff_names[ii]);
             }
             if(eff_pos[ii] > 0){
             	fprintf(master_header,
                 "\tcomponent \"positions\" value %d\n",eff_pos[ii]);
             }
             if(eff_orient[ii] > 0){     
             	fprintf(master_header,
                 	"\tcomponent \"data\" value %d # orientations #\n",eff_orient[ii]);
             }
             if(viz_mol_states_flag)
             {
                if(eff_states[ii] > 0){
             	   fprintf(master_header,
                	"\tcomponent \"state_values\" value %d\n",eff_states[ii]);
                }
             }
             fprintf(master_header, "\n");
             eff_field_indices[ii] = main_index;
             main_index++;
      }
  }




/* dump 3D molecules: */

 if(viz_mol_pos_flag || viz_mol_orient_flag){	 
   
        /* create references to the molecules. */
        if ((viz_molp=(struct molecule ***)malloc(n_species*sizeof(struct molecule **)))==NULL) {
                fprintf(world->err_file, "File %s, Line %ld: memory allocation error.\n", __FILE__, (long)__LINE__);
      		return(1);
    	}
    	if ((viz_mol_count=(u_int *)malloc(n_species*sizeof(u_int)))==NULL) {
                fprintf(world->err_file, "File %s, Line %ld: memory allocation error.\n", __FILE__, (long)__LINE__);
      		return(1);
    	}
  
    for (ii=0;ii<n_species;ii++) {
      /* perform initialization */
      spec_id=species_list[ii]->species_id;
      viz_molp[spec_id]=NULL;
      viz_mol_count[spec_id]=0;

      if (species_list[ii]->viz_state == EXCLUDE_OBJ)  continue;

      num=species_list[ii]->population;
      if ((num>0) && (species_list[ii]->flags & NOT_FREE) == 0) {
          /* create references for 3D molecules */
          if ((viz_molp[spec_id]=(struct molecule **)malloc
            (num*sizeof(struct molecule *)))==NULL) {
                fprintf(world->err_file, "File %s, Line %ld: memory allocation error.\n", __FILE__, (long)__LINE__);
                return(1);
          }
      }
    } /* end for */

    slp=world->storage_head;
    while (slp!=NULL) {
      sp=slp->store;
      shp=sp->timer;
      while (shp!=NULL) {
        
        for (ii=0;ii<shp->buf_len;ii++) {
          amp=(struct abstract_molecule *)shp->circ_buf_head[ii];
          while (amp!=NULL) {
            if ((amp->properties!=NULL) && (amp->flags & TYPE_3D) == TYPE_3D){
              molp=(struct molecule *)amp;
              if (molp->properties->viz_state!=EXCLUDE_OBJ) {
                spec_id=molp->properties->species_id;
                if (viz_mol_count[spec_id]<molp->properties->population) {
                  viz_molp[spec_id][viz_mol_count[spec_id]++]=molp;
                }
                else {
                  fprintf(log_file,"MCell: molecule count disagreement!!\n");
                  fprintf(log_file,"  Species %s  population = %d  count = %d\n",molp->properties->sym->name,molp->properties->population,viz_mol_count[spec_id]);
                } /* end if/else */
              }  /* end if(molp) */
            } /* end if(amp) */
            amp=amp->next;
          } /* end while  (amp) */
        } /* end for */
        
        amp=(struct abstract_molecule *)shp->current;
        while (amp!=NULL) {
          if ((amp->properties!=NULL) && (amp->flags & TYPE_3D) == TYPE_3D) {
            molp=(struct molecule *)amp;
            if (molp->properties->viz_state!=EXCLUDE_OBJ) {
              spec_id=molp->properties->species_id;
              if (viz_mol_count[spec_id]<molp->properties->population) {
                viz_molp[spec_id][viz_mol_count[spec_id]++]=molp;
              }
              else {
                fprintf(log_file,"MCell: molecule count disagreement!!\n");
                fprintf(log_file,"  Species %s  population = %d  count = %d\n",molp->properties->sym->name,molp->properties->population,viz_mol_count[spec_id]);
              }
            }

           }
          amp=amp->next;
        }
        
        shp=shp->next_scale;
      } /* end (while (shp) */

      slp=slp->next;
    } /* end (while slp) */


    for (ii=0;ii<n_species;ii++) {
      
      if(species_list[ii]->viz_state == EXCLUDE_OBJ) continue; 
      spec_id=species_list[ii]->species_id;
      
      if(viz_mol_count[spec_id] > 0)
      {
         num=viz_mol_count[spec_id];
         if(num!=species_list[ii]->population
             && ((species_list[ii]->flags & NOT_FREE)==0)) {
             fprintf(log_file,"MCell: molecule count disagreement!!\n");
             fprintf(log_file,"  Species %s  population = %d  count = %d\n",species_list[ii]->sym->name,species_list[ii]->population,num);
         }
      }

      if (viz_mol_count[spec_id]>0 && ((species_list[ii]->flags & NOT_FREE) == 0)) {
        /* here are 3D diffusing molecules */
        num = viz_mol_count[spec_id];
       if(viz_mol_pos_flag)
       { 	  
	  mol_pos_byte_offset_prev = mol_pos_byte_offset;
          for (jj=0;jj<num;jj++) {
            molp=viz_molp[spec_id][jj];
	    v1=(float)(world->length_unit*molp->pos.x);
	    v2=(float)(world->length_unit*molp->pos.y);
	    v3=(float)(world->length_unit*molp->pos.z);
	    fwrite(&v1,sizeof v1,1,mol_pos_data);
	    fwrite(&v2,sizeof v2,1,mol_pos_data);
	    fwrite(&v3,sizeof v3,1,mol_pos_data);
            mol_pos_byte_offset += (sizeof(v1) + sizeof(v2) + sizeof(v3));
          }
          fprintf(master_header,"object %d class array type float rank 1 shape 3 items %d %s binary data file %s,%d  # %s positions #\n",main_index,num,my_byte_order, mol_pos_name, mol_pos_byte_offset_prev, species_list[ii]->sym->name);
          fprintf(master_header,"\tattribute \"dep\" string \"positions\"\n\n");
          mol_pos[mol_pos_index] = main_index;
          mol_pos_index++;
          main_index++;
        }

        if(viz_mol_orient_flag)
        {  
          /* write molecule orientations information. 
             for 3D molecules we use default orientation [0,0,1] */
      	  mol_orient_byte_offset_prev = mol_orient_byte_offset;
          v1 = 0;
          v2 = 0;
          v3 = 1;
          fwrite(&v1,sizeof v1,1,mol_orient_data);
          fwrite(&v2,sizeof v2,1,mol_orient_data);
          fwrite(&v3,sizeof v3,1,mol_orient_data);
          mol_orient_byte_offset += (sizeof(v1) + sizeof(v2) + sizeof(v3));
          fprintf(master_header,"object %d class constantarray type float rank 1 shape 3 items %d %s binary data file %s,%d # %s orientations #\n",main_index,num, my_byte_order, mol_orient_name, mol_orient_byte_offset_prev, species_list[ii]->sym->name);
          fprintf(master_header,"\tattribute \"dep\" string \"positions\"\n\n");
          mol_orient[mol_orient_index] = main_index;
          mol_orient_index++;
          main_index++;
        }

        if(viz_mol_states_flag)
        {
          /* write molecule states information. */ 
      	  mol_states_byte_offset_prev = mol_states_byte_offset;
          fwrite(&state,sizeof state,1,mol_states_data);
          mol_states_byte_offset += (sizeof state);
          fprintf(master_header,"object %d class constantarray type int items %d %s binary data file %s,%d  # %s states #\n",main_index,num, my_byte_order,mol_states_name,mol_states_byte_offset_prev, species_list[ii]->sym->name);
          fprintf(master_header,"\tattribute \"dep\" string \"positions\"\n\n");
          mol_states[mol_states_index] = main_index;
          mol_states_index++;
          main_index++;
        }
      }
      /* output empty arrays for zero molecule counts here */
      else if ((viz_mol_count[spec_id]==0) && ((species_list[ii]->flags & NOT_FREE) == 0)) {
        if(viz_mol_pos_flag)
        {
          fprintf(master_header,"object %d array   # %s positions #\n",main_index, species_list[ii]->sym->name);
          mol_pos[mol_pos_index] = main_index;
          mol_pos_index++;
          main_index++;
        }
        if(viz_mol_orient_flag)
        {
          fprintf(master_header,"object %d array   # %s orientations #\n",main_index, species_list[ii]->sym->name);
          mol_orient[mol_orient_index] = main_index;
          mol_orient_index++;
          main_index++;
        }
       
       if(viz_mol_states_flag)
       {
          fprintf(master_header,"object %d array   # %s states #\n",main_index, species_list[ii]->sym->name);
          mol_states[mol_states_index] = main_index;
          mol_states_index++;
          main_index++;
          
       }
	 fprintf(master_header,"\n");
      } /* end else if */
    }
  } /* end if((viz_mol_pos_flag) || (viz_mol_orient_flag)) */
 

/* build fields here */
      int show_molecules = 0;
   if(mol_to_show_number > 0)
   {
      if(fdlp->type == ALL_MOL_DATA) {
           show_molecules = 1;
      }
      else if((viz_mol_pos_flag || viz_mol_orient_flag) && ((special_mol_frames_counter ==0) || (special_mol_frames_counter == 2))){
		show_molecules = 1;
      }
   }

    if(show_molecules)
    {
      for (ii=0; ii<mol_to_show_number; ii++) {
               if(mol_names[ii] != NULL){
                  fprintf(master_header,
                      "object %d field   # %s #\n",main_index, mol_names[ii]);
                }
                if(mol_pos[ii] > 0) {
             	   fprintf(master_header,
                 	"\tcomponent \"positions\" value %d\n",mol_pos[ii]);
                }
                if(mol_orient[ii] > 0){
             	   fprintf(master_header,
                 	"\tcomponent \"data\" value %d # orientations #\n",mol_orient[ii]);
                }
             if(viz_mol_states_flag)
             {
                if(mol_states[ii] > 0){ 
             	   fprintf(master_header,
                	"\tcomponent \"state_values\" value %d\n",mol_states[ii]);
                }
             }
             fprintf(master_header, "\n");
             mol_field_indices[ii] = main_index;
             main_index++;
      }
    }

      /* create group objects for molecules/effectors */

        if(show_molecules)
        {
                fprintf(master_header,"object \"%s%u\" group # %s #\n", "volume_molecules_", viz_iteration, "volume molecules"); 
        	mol_group_index = main_index;
                member_molecules_iteration = viz_iteration;
      		for (ii=0;ii<mol_to_show_number;ii++) {
                   if(mol_field_indices[ii] > 0){
          		fprintf(master_header,"\tmember \"%s\" value %d\n",mol_names[ii], mol_field_indices[ii]);
                    }
          	}
                main_index++;
            /* store iteration_number for volume molecules */
            ia_uint_store(&frame_numbers_vol_mols,frame_numbers_vol_mols_count, 
                                viz_iteration);
            frame_numbers_vol_mols_count++;
      	}

        fprintf(master_header, "\n"); 

 
        if(show_effectors)
        {
                fprintf(master_header,"object \"%s%u\" group # %s #\n","surface_molecules_", viz_iteration, "surface molecules");
               eff_group_index = main_index;
               member_effectors_iteration = viz_iteration;
      	       for (ii=0;ii<eff_to_show_number;ii++) {
                   if(eff_field_indices[ii] > 0){
                       fprintf(master_header,"\tmember \"%s\" value %d\n",eff_names[ii], eff_field_indices[ii]);
                   }
               }
               main_index++;
            /* store iteration_number for surface molecules */
            ia_uint_store(&frame_numbers_surf_mols,
                      frame_numbers_surf_mols_count, viz_iteration);
            frame_numbers_surf_mols_count++;
      	}
        fprintf(master_header, "\n"); 

      /* create combined group object for meshes and molecules */
      int show_combined_group = 1;
      if(((fdlp->viz_iterationll == curr_surf_pos_iteration_step) ||
	  (fdlp->viz_iterationll == curr_region_data_iteration_step)) && 
        (mesh_group_index == 0)) 
      {
            	show_combined_group = 0;
      }
      if(((fdlp->viz_iterationll == curr_mol_pos_iteration_step) ||
	 (fdlp->viz_iterationll == curr_mol_orient_iteration_step))
             && ((mol_group_index == 0) && (eff_group_index == 0))){ 
            	show_combined_group = 0;
      }

      if(show_combined_group)
      {
        fprintf(master_header,"object %d group\n",main_index);
        combined_group_index = main_index;
        if(member_meshes_iteration != UINT_MAX){  
          	fprintf(master_header,"\tmember \"meshes\" value \"%s%u\"\n", "meshes_",member_meshes_iteration);
	     	mesh_group_index = 0; 

      	}
        if(member_molecules_iteration != UINT_MAX) {  
          	fprintf(master_header,"\tmember \"volume_molecules\" value \"%s%u\"\n",  "volume_molecules_", member_molecules_iteration); 
		mol_group_index = 0; 
        }
       if(member_effectors_iteration != UINT_MAX) {  
          	fprintf(master_header,"\tmember \"surface_molecules\" value \"%s%u\"\n",  "surface_molecules_", member_effectors_iteration);
		eff_group_index = 0; 
        }
      	
        fprintf(master_header,"\n");
      }


     /* create entry into "frame_data" object. */
      if(show_combined_group){

        /* create an entry into a 'frame_data' object. */
      	sprintf(buffer, "\tmember %d value %d position %u\n", series_index, combined_group_index, viz_iteration);
	ia_string_store(&frame_data_series_list, frame_data_series_count, buffer);
        frame_data_series_count++;
      	series_index++;
      	main_index++;

      	fprintf(master_header, "\n\n");

        if(special_surf_frames_counter == 2){
		special_surf_frames_counter = 0;
        }
        if(special_mol_frames_counter == 2){
		special_mol_frames_counter = 0;
        }
          

	/* put value of viz_iteration into the time_values array */ 
      ia_uint_store(&time_values, time_values_count,viz_iteration);  
        time_values_count++;

     }

     /* flag that signals the need to write "footers" for 
        the master file */
     int time_to_write_footers = 0;
     /* flag pointing to another frame with certain condition */
     int found = 0;
     struct frame_data_list *fdlp_temp;
     long long next_iteration_step = 0;  /* next iteration for this frame */
     long long next_iteration_step_temp = INT_MAX; /*next iteration for other frame */
     

     if(fdlp->curr_viz_iteration->next != NULL){
        next_iteration_step = (long long)(fdlp->curr_viz_iteration->next->value);
     }

     if(world->chkpt_flag){
        /* check whether it is the last frame */
        

        if((world->it_time == (world->start_time + world->chkpt_iterations)) ||
                 (world->it_time == last_iteration)){


             /* look ahead to find out whether there are 
                other frames to be output */
             for(fdlp_temp = fdlp->next; fdlp_temp != NULL; fdlp_temp = fdlp_temp->next){
                if((fdlp_temp->viz_iterationll == (world->start_time + world->chkpt_iterations)) || (fdlp_temp->viz_iterationll == last_iteration)){
                   found = 1;
                   break;
                }
             }
             
             if(!found){
                /* this is the last frame */
                    time_to_write_footers = 1; 
             }
        }
        else if(next_iteration_step > (world->start_time + world->chkpt_iterations)) {
             /* look ahead to find out whether next_iteration_step
                for other frames is less than 'world->chkpt_iterations'
             */
             for(fdlp_temp = fdlp->next; fdlp_temp != NULL; fdlp_temp = fdlp_temp->next){
                if((fdlp_temp->curr_viz_iteration != NULL) && 
                    (fdlp_temp->curr_viz_iteration->next != NULL)){
                    next_iteration_step_temp = (long long)(fdlp_temp->curr_viz_iteration->next->value);
                   if(next_iteration_step_temp  < world->start_time + world->chkpt_iterations){
                      found = 1;
                      break;
                   }
                }
             }
             if(!found){
                /* this is the last frame */
                    time_to_write_footers = 1; 
             }
        }else{}

    }else{
          if(world->it_time == last_iteration){

             /* look ahead to find out whether there are 
                other frames to be output */
             for(fdlp_temp = fdlp->next; fdlp_temp != NULL; fdlp_temp = fdlp_temp->next){
                if(fdlp_temp->viz_iterationll == last_iteration){
                   found = 1;
                   break;
                }
             }
             if(!found){
                /* this is the last frame */
                    time_to_write_footers = 1; 
             }
           }

    }   


     if(time_to_write_footers)
     {
       
        u_int elem1;
        double t_value;
        int extra_elems;
        int frame_numbers_count;

	/* write 'frame_numbers' object. */
        if(world->chkpt_flag){
      	   sprintf(file_name,"%s.frame_numbers.%u.bin",                                    world->file_prefix_name, world->chkpt_seq_num);
        }else{
      	   sprintf(file_name,"%s.frame_numbers.bin",world->file_prefix_name);
        }

     /* remove the folder name from the frame_numbers data file name */
     	ch_ptr = strrchr(file_name, '/');
     	++ch_ptr;
     	strcpy(frame_numbers_name, ch_ptr);

      	if ((frame_numbers_data=fopen(file_name,"wb"))==NULL) {
            fprintf(world->err_file, "File %s, Line %ld: cannot open file %s.\n", __FILE__, (long)__LINE__, file_name);
           return(1);
        }

        /* write "time_values" object. */
        if(world->chkpt_flag){
     	    sprintf(file_name,"%s.time_values.%u.bin",world->file_prefix_name, world->chkpt_seq_num);
        }else{
     	    sprintf(file_name,"%s.time_values.bin",world->file_prefix_name);
        }

     	/* remove the folder name from the time_values data file name */
     	ch_ptr = strrchr(file_name, '/');
     	++ch_ptr;
     	strcpy(time_values_name, ch_ptr);

      	if ((time_values_data=fopen(file_name,"wb"))==NULL) {
            fprintf(world->err_file, "File %s, Line %ld: cannot open file %s.\n", __FILE__, (long)__LINE__, file_name);
            return(1);
        }

        if(frame_numbers_meshes_count > 0)
        {
           extra_elems = frame_numbers_vol_mols_count - frame_numbers_meshes_count;
           if(extra_elems > 0){
              /* pad the frame_numbers_meshes array with the last 
                 element so that it will have the same number of elements 
                 as frame_numbers_vol_mols array */
              elem1 = ia_uint_get(&frame_numbers_meshes, frame_numbers_meshes_count - 1);
              
		for(ii = 0; ii < extra_elems; ii++){
                   ia_uint_store(&frame_numbers_meshes, frame_numbers_meshes_count + ii, elem1);
                }
           }
        }

        if(frame_numbers_vol_mols_count > 0)
        {
           extra_elems = frame_numbers_meshes_count - frame_numbers_vol_mols_count;
           if(extra_elems > 0){
              /* pad the frame_numbers_vol_mols array with the last 
                 element so that it will have the same number of elements 
                 as frame_numbers_meshes array */
              elem1 = ia_uint_get(&frame_numbers_vol_mols, frame_numbers_vol_mols_count - 1);
              
		for(ii = 0; ii < extra_elems; ii++){
                   ia_uint_store(&frame_numbers_vol_mols, frame_numbers_vol_mols_count + ii, elem1);
                }
           }
        }

        if(frame_numbers_surf_mols_count > 0)
        {
           extra_elems = frame_numbers_meshes_count - frame_numbers_surf_mols_count;
           if(extra_elems > 0){
              /* pad the frame_numbers_surf_mols array with the last 
                 element so that it will have the same number of elements 
                 as frame_numbers_meshes array */
              elem1 = ia_uint_get(&frame_numbers_surf_mols, frame_numbers_surf_mols_count - 1);
              
		for(ii = 0; ii < extra_elems; ii++){
                   ia_uint_store(&frame_numbers_surf_mols, frame_numbers_surf_mols_count + ii, elem1);
                }
           }
        }

      if(frame_numbers_vol_mols_count > frame_numbers_meshes_count){
         frame_numbers_count = frame_numbers_vol_mols_count;
      }else{
         frame_numbers_count = frame_numbers_meshes_count;
      }
       
      fprintf(master_header,"object \"frame_numbers\" class array  type unsigned int rank 1 shape 3 items %u %s binary data file %s,%d\n",frame_numbers_count, my_byte_order,frame_numbers_name, frame_numbers_byte_offset);
      for(ii = 0; ii < frame_numbers_count; ii++){
                elem1 = ia_uint_get(&frame_numbers_meshes, ii);
                if(elem1 == UINT_MAX) 
                {
                   fprintf(world->err_file, "File %s, Line %ld: ia_uint_get() tries to access uninitialized data.\n", __FILE__, (long)__LINE__);
                   return 1;
                }
                fwrite(&elem1, sizeof(elem1),1,frame_numbers_data);

                elem1 = ia_uint_get(&frame_numbers_vol_mols, ii);
                if(elem1 == UINT_MAX) 
                {
                   fprintf(world->err_file, "File %s, Line %ld: ia_uint_get() tries to access uninitialized data.\n", __FILE__, (long)__LINE__);
                   return 1;
                }
                fwrite(&elem1, sizeof(elem1),1,frame_numbers_data);

                elem1 = ia_uint_get(&frame_numbers_surf_mols, ii);
                if(elem1 == UINT_MAX) 
                {
                   fprintf(world->err_file, "File %s, Line %ld: ia_uint_get() tries to access uninitialized data.\n", __FILE__, (long)__LINE__);
                   return 1;
                }
                fwrite(&elem1, sizeof(elem1),1,frame_numbers_data);

        }
     	fprintf(master_header, "\n\n");



        if(time_values_count > 0)
        {
        	fprintf(master_header,"object \"time_values\" class array  type double rank 0 items %u %s binary data file %s,%d\n",time_values_count, my_byte_order,time_values_name, time_values_byte_offset);
												for(ii = 0; ii < time_values_count; ii++){
                	elem1 = ia_uint_get(&time_values, ii);
                        if(elem1 == UINT_MAX) 
                        {
                          fprintf(world->err_file, "File %s, Line %ld: ia_uint_get() tries to access uninitialized data.\n", __FILE__, (long)__LINE__);
                          return 1;
                        }
                        t_value = elem1*world->time_unit;
                        fwrite(&(t_value), sizeof(t_value),1,time_values_data);
                 }
		fprintf(master_header, "\n\n");
	}
                /* write 'frame_data' object. */
    		fprintf(master_header,"object \"%s\" class series\n", "frame_data");
        	if(frame_data_series_count > 0)
        	{
        	   char *elem;
		   for(ii = 0; ii < frame_data_series_count; ii++){
                       elem = ia_string_get(&frame_data_series_list, ii);
                       if(elem == NULL){
                          fprintf(world->err_file, "File %s, Line %ld: ia_string_get() tries to access uninitialized data.\n", __FILE__, (long)__LINE__);
                        }
			fprintf(master_header, "\t%s", elem);
        	    }
                 }
	         fprintf(master_header, "\n\n");
    } /* end if(time_to_write_footers) */

 
   /* free allocated memory */
   if (viz_molp != NULL) {
    for (ii=0;ii<n_species;ii++) {
      if (viz_molp[ii]!=NULL) {
        free(viz_molp[ii]);
      }
    }
    
    free(viz_molp);
  }
   
  if (viz_mol_count != NULL) {
    free (viz_mol_count);
  }


  if (viz_grid_mol_count != NULL) {
    free (viz_grid_mol_count);
  }
  
    if(master_header != NULL){
    	fclose(master_header);
    }
    if(mol_pos_data != NULL){
    	fclose(mol_pos_data);
    }
    if(mol_orient_data != NULL){
    	fclose(mol_orient_data);
    }
    if(mol_states_data != NULL){
    	fclose(mol_states_data);
    }
    if(mesh_pos_data != NULL){
    	fclose(mesh_pos_data);
    }
    if(mesh_states_data != NULL){
    	fclose(mesh_states_data);
    }
    if(region_data != NULL){
    	fclose(region_data);
    }
    if(frame_numbers_data != NULL){
        fclose(frame_numbers_data);
    }
    if(time_values_data != NULL){
        fclose(time_values_data);
    }



  return(0);

}


/************************************************************************ 
output_rk_custom:
Rex Kerr's personal visualization mode output function 
*************************************************************************/

int output_rk_custom(struct frame_data_list *fdlp)
{
  FILE *log_file;
  FILE *custom_file;
  char cf_name[1024];
  struct storage_list *slp;
  struct schedule_helper *shp;
  struct abstract_element *aep;
  struct abstract_molecule *amp;
  struct molecule *mp;
  struct grid_molecule *gmp;
  short orient = 0;
  struct species *target;
  
  int i,j,k;
  double d;

  int id;
  struct vector3 where;
  
  no_printf("Output in CUSTOM_RK mode...\n");
  log_file = world->log_file;
  
  if (world->rk_mode_var==NULL) return output_ascii_molecules(fdlp);
  
  if ((fdlp->type==ALL_FRAME_DATA) || (fdlp->type==MOL_POS) || (fdlp->type==MOL_STATES))
  {
    sprintf(cf_name,"%s.rk.dat",world->molecule_prefix_name);
    custom_file = fopen(cf_name,(world->rk_mode_var->n_written)?"a+":"w");
    if (!custom_file)
    {
      fprintf(log_file,"Couldn't open file %s for viz output.\n",cf_name);
      return 1;
    }
    else no_printf("Writing to file %s\n",cf_name);
    
    world->rk_mode_var->n_written++;
    fprintf(custom_file,"%lld",fdlp->viz_iterationll);
    for (j=0;j<world->n_species;j++)
    {
      target = world->species_list[j];
      if (target==NULL) continue;
      if (target->viz_state==EXCLUDE_OBJ) continue;
      for (k=0;k<world->rk_mode_var->n_bins;k++) world->rk_mode_var->bins[k] = 0;
      id = target->viz_state;
      
      fprintf(custom_file,"\t%d",id);
      
      for (slp = world->storage_head ; slp != NULL ; slp = slp->next)
      {
	for (shp = slp->store->timer ; shp != NULL ; shp = shp->next_scale)
	{
	  for (i=-1;i<shp->buf_len;i++)
	  {
	    if (i<0) aep = shp->current;
	    else aep = shp->circ_buf_head[i];
	    
	    for (aep=(i<0)?shp->current:shp->circ_buf_head[i] ; aep!=NULL ; aep=aep->next)
	    {
	      amp = (struct abstract_molecule*)aep;
	      if (amp->properties != target) continue;
  
	      if ((amp->properties->flags & NOT_FREE)==0)
	      {
		mp = (struct molecule*)amp;
		where.x = mp->pos.x;
		where.y = mp->pos.y;
		where.z = mp->pos.z;
	      }
	      else if ((amp->properties->flags & ON_GRID)!=0)
	      {
		gmp = (struct grid_molecule*)amp;
		uv2xyz(&(gmp->s_pos),gmp->grid->surface,&where);
		orient = gmp->orient;
	      }
	      else continue;
	      
	      d = dot_prod(&where , world->rk_mode_var->direction);
	      
	      k = bin(world->rk_mode_var->parts,world->rk_mode_var->n_bins-1,d);
	      world->rk_mode_var->bins[k]++;
	    }
	  }
	}
      }
      for (i=k=0 ; k<world->rk_mode_var->n_bins ; k++) i+=world->rk_mode_var->bins[k];
      if (i!=target->population) printf("Wanted to bin %d but found %d instead\n",target->population,i);
      for (k=0 ; k<world->rk_mode_var->n_bins ; k++) fprintf(custom_file," %d",world->rk_mode_var->bins[k]);
    }
    fprintf(custom_file,"\n");
    fclose(custom_file);
  }
  
  return 0;
}


/************************************************************************ 
output_ascii_molecules:
In: a frame data list (internal viz output data structure)
Out: 0 on success, 1 on failure.  The positions of molecules are output
     in exponential floating point notation (with 8 decimal places)
*************************************************************************/

int output_ascii_molecules(struct frame_data_list *fdlp)
{
  FILE *log_file;
  FILE *custom_file;
  char cf_name[1024];
  char cf_format[256];
  struct storage_list *slp;
  struct schedule_helper *shp;
  struct abstract_element *aep;
  struct abstract_molecule *amp;
  struct molecule *mp;
  struct grid_molecule *gmp;
  short orient = 0;
  
  int ndigits,i;
  long long lli;

  int id;
  struct vector3 where;
  
  no_printf("Output in ASCII mode (molecules only)...\n");
  log_file = world->log_file;
  
  if ((fdlp->type==ALL_FRAME_DATA) || (fdlp->type==MOL_POS) || (fdlp->type==MOL_STATES))
  {
    lli = 10;
    for (ndigits = 1 ; lli <= world->iterations && ndigits<20 ; lli*=10 , ndigits++) {}
    sprintf(cf_format,"%%s.ascii.%%0%dlld.dat",ndigits);
    sprintf(cf_name,cf_format,world->molecule_prefix_name,fdlp->viz_iterationll);
    custom_file = fopen(cf_name,"w");
    if (!custom_file)
    {
      fprintf(log_file,"Couldn't open file %s for viz output.\n",cf_name);
      return 1;
    }
    else no_printf("Writing to file %s\n",cf_name);
    
    for (slp = world->storage_head ; slp != NULL ; slp = slp->next)
    {
      for (shp = slp->store->timer ; shp != NULL ; shp = shp->next_scale)
      {
        for (i=-1;i<shp->buf_len;i++)
        {
          if (i<0) aep = shp->current;
          else aep = shp->circ_buf_head[i];
          
          for (aep=(i<0)?shp->current:shp->circ_buf_head[i] ; aep!=NULL ; aep=aep->next)
          {
            amp = (struct abstract_molecule*)aep;
            if (amp->properties == NULL) continue;
            if (amp->properties->viz_state == EXCLUDE_OBJ) continue;

            id = amp->properties->viz_state;
            
            if ((amp->properties->flags & NOT_FREE)==0)
            {
              mp = (struct molecule*)amp;
              where.x = mp->pos.x;
              where.y = mp->pos.y;
              where.z = mp->pos.z;
            }
            else if ((amp->properties->flags & ON_GRID)!=0)
            {
              gmp = (struct grid_molecule*)amp;
              uv2xyz(&(gmp->s_pos),gmp->grid->surface,&where);
              orient = gmp->orient;
            }
            else continue;
            
            where.x *= world->length_unit;
            where.y *= world->length_unit;
            where.z *= world->length_unit;
            fprintf(custom_file,"%d %15.8e %15.8e %15.8e %2d\n",id,where.x,where.y,where.z,orient);
          }
        }
      }
    }
    fclose(custom_file);
  }
  
  return 0;
}


/* **************************************************************** */

#if 0





int output_rayshade_objects(struct frame_data_list *fdlp)
{
FILE *wall_header, *wall_verts;
FILE *eff_header, *eff_states, *eff_pos;
FILE *lig_header, *lig_states, *lig_pos;
struct viz_obj *vizp;
struct cmprt_data_list *cdlp;
struct cmprt_data *cdp;
struct wall_list *wlp;
struct wall *wp;
struct effector *rp;
struct rx *rxp;
struct ligand_info *ligip;
struct ligand *ligp;
struct vector3 glyph_axis,rotation_axis;
struct vector3 ab,bc,step_u,step_v,diagonal;
struct vector3 p0,p1,p2;
double i,j,j_max,i_max,i1,i2,i3,j1,j2,j3;
double rad,dot,length,rotation_angle,final_angle,v_val;
float v1,v2,v3,n1,n2,n3;
int ii;
int uu,vv,vv_max,grid_size;
int num;
int viz_iteration,n_viz_iterations;
int first_viz_iteration;
int wall_index,eff_index,lig_index,r_index,orient;
int pos_count,state_count,lig_count;
int state;
int viz_type;
unsigned int index,n_eff;
byte viz_eff,viz_lig,viz_surf;
byte viz_surf_or_eff;
byte viz_eff_pos,viz_eff_states;
byte viz_lig_pos,viz_lig_states;
byte viz_surf_pos,viz_surf_states;
byte grid_shape;
char file_name[1024];



viz_iteration=(int)(fdlp->viz_iteration/time_unit + ROUND_UP);
n_viz_iterations=fdlp->n_viz_iterations;
first_viz_iteration=(viz_iteration==fdlp->iteration_list->value);

viz_type=fdlp->type;
viz_eff=((viz_type==ALL_FRAME_DATA) || (viz_type==EFF_POS)
             || (viz_type==EFF_STATES));
viz_lig=((viz_type==ALL_FRAME_DATA) || (viz_type==MOL_POS)
             || (viz_type==MOL_STATES));
viz_surf=((viz_type==ALL_FRAME_DATA) || (viz_type==SURF_POS)
             || (viz_type==SURF_STATES));
viz_surf_or_eff=(viz_surf || viz_eff);


viz_eff_pos=((viz_type==ALL_FRAME_DATA) || (viz_type==EFF_POS));
viz_eff_states=((viz_type==ALL_FRAME_DATA) || (viz_type==EFF_STATES));
viz_lig_pos=((viz_type==ALL_FRAME_DATA) || (viz_type==MOL_POS));
viz_lig_states=((viz_type==ALL_FRAME_DATA) || (viz_type==MOL_STATES));
viz_surf_pos=((viz_type==ALL_FRAME_DATA) || (viz_type==SURF_POS));
viz_surf_states=((viz_type==ALL_FRAME_DATA) || (viz_type==SURF_STATES));


rad=360/(2*MY_PI);

glyph_axis.x = 0;
glyph_axis.y = 0;
glyph_axis.z = 1;

/* dump walls and effectors: */
  if (viz_surf_or_eff) {
    vizp = viz_obj_head;
    while(vizp!=NULL) {
      if (viz_surf) {
        sprintf(file_name,"%s.surfaces.%d.ray",vizp->name,viz_iteration);
        if ((wall_header=fopen(file_name,"w"))==NULL) {
          fprintf(log_file,"MCell: error cannot open rayshade surface file %s\n",file_name);
          return(1);
        }
        fprintf(wall_header,"/* surfaces */\n");
      }

      if (viz_eff) {
        sprintf(file_name,"%s.effectors.%d.ray",vizp->name,viz_iteration);
        if ((eff_header=fopen(file_name,"w"))==NULL) {
          fprintf(log_file,"MCell: error cannot open effector header file %s\n",file_name);
          return(1);
        }
        fprintf(eff_header,"/* effectors */\n");
      }

      wall_index = 3;
      eff_index = 2;
      pos_count = 0;
      state_count = 0;
      cdlp = vizp->cmprt_data_list;
      while(cdlp!=NULL) {
        cdp = cdlp->cmprt_data;
        no_printf("Traversing surfaces in compartment %s\n",cdp->sym->name);
        fflush(log_file);
        wlp = cdp->wall_list;
        while (wlp!=NULL) {
          wp = wlp->wall;
 
          /* dump the effectors */
          if (viz_eff) {
            rp = wp->effectors;
            if (rp!=NULL) {
              n_eff=0;
              for (index=0;index<rp->n_tiles;index++) {
	        rxp=rp->tiles[index];
	        if (rxp!=NULL) {
	          state=rxp->viz_state;
	        }
	        else {
	          state=EXCLUDE_OBJ;
	        }
                if (state!=EXCLUDE_OBJ) {
                  n_eff++;
                }
              }
              no_printf("Dumping %d effectors...\n",n_eff);
              fflush(log_file);
              grid_shape=rp->grid_shape;
              grid_size=rp->grid_size;
              vectorize(wp->vert[0],wp->vert[1],&ab);
              vectorize(wp->vert[1],wp->vert[2],&bc);
              step_u.x=ab.x/grid_size;
              step_u.y=ab.y/grid_size;
              step_u.z=ab.z/grid_size;
              step_v.x=bc.x/grid_size;
              step_v.y=bc.y/grid_size;
              step_v.z=bc.z/grid_size;
              vectorize(&step_u,&step_v,&diagonal);
  
	      cross_prod(&glyph_axis,&(wp->normal),&rotation_axis);
	      dot=dot_prod(&glyph_axis,&(wp->normal))*0.999999999999;
	      rotation_angle = acos(dot)*rad;
/*
	      length = vect_length(&rotation_axis)*0.999999999999;
	      rotation_angle = asin(length)*rad;
*/
  
              index=0;
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
	          rxp=rp->tiles[index];
	          if (rxp!=NULL) {
	            state=rxp->viz_state;
	            orient=(rp->orient[rxp->rx_index]-1)/-2;
	          }
	          else {
	            state=EXCLUDE_OBJ;
	          }
      
	          if (state!=EXCLUDE_OBJ) {
	            final_angle=rotation_angle-(orient*180);
	            v1=length_unit*(wp->vert[0]->x+p2.x);
	            v2=length_unit*(wp->vert[0]->y+p2.y);
	            v3=length_unit*(wp->vert[0]->z+p2.z);
	            fprintf(eff_header,"object effector_mat_%d effector_glyph\n",state);
		    if (final_angle!=0) {
	              fprintf(eff_header,"rotate %.6g %.6g %.6g %.6g\n",rotation_axis.x,
		        rotation_axis.y,rotation_axis.z,final_angle);
	            }
	            fprintf(eff_header,"translate %.6g %.6g %.6g\n",v1,v2,v3);
	          }
                  index++;
                }
              }
            }
          }

          /* dump the wall */
          if (viz_surf) {
	    if (wp->viz_state!=EXCLUDE_OBJ) {
              if (wp->n_vert==3) {
                fprintf(wall_header,"triangle surface_mat_%d\n",wp->viz_state);
                for (ii=0;ii<wp->n_vert;ii++) {
                  v1 = length_unit*wp->vert[ii]->x;
                  v2 = length_unit*wp->vert[ii]->y;
                  v3 = length_unit*wp->vert[ii]->z;
                  fprintf(wall_header,"%.6g %.6g %.6g\n",v1,v2,v3);
                  if (wp->vert_normal!=NULL) {
                    n1 = wp->vert_normal[ii]->x;
                    n2 = wp->vert_normal[ii]->y;
                    n3 = wp->vert_normal[ii]->z;
                    fprintf(wall_header,"%.6g %.6g %.6g\n",n1,n2,n3);
                  }
                }
              }
              else {
                fprintf(wall_header,"poly surface_mat_%d\n",wp->viz_state);
                for (ii=0;ii<wp->n_vert;ii++) {
                  v1 = length_unit*wp->vert[ii]->x;
                  v2 = length_unit*wp->vert[ii]->y;
                  v3 = length_unit*wp->vert[ii]->z;
                  fprintf(wall_header,"%.6g %.6g %.6g\n",v1,v2,v3);
                }
              }
	    }
          }
          wlp = wlp->next;
        }

        cdlp = cdlp->next;
      }

      if (viz_surf) {
        fclose(wall_header);
      }
      if (viz_eff) {
        fclose(eff_header);
      }
      vizp = vizp->next;
    }
  }

/* dump ligands: */
  if (viz_lig) {
    if (n_ligand_types>0) {
      sprintf(file_name,"ligands.%d.ray",viz_iteration);
      if ((lig_header=fopen(file_name,"w"))==NULL) {
        fprintf(log_file,"MCell: error cannot open molecule header file %s\n",file_name);
        return(1);
      }
  
      fprintf(lig_header,"/* ligands */\n");
    }

    lig_index=2;
    pos_count=0;
    state_count=0;
    for (ii=1;ii<1+n_ligand_types;ii++) {
      ligip=ligand_table[ii];
      ligp=ligip->top;
      num=ligip->lig_index;
      lig_count=0;
      if (num>0) {
        while (ligp!=NULL) {
  
          rp=ligp->effector;
          if (rp!=NULL) {
            r_index=ligp->index;
            state=rp->tiles[r_index]->viz_state;
          }
	  else {
	    state=ligip->viz_state;
	  }

	  if (state!=EXCLUDE_OBJ) {
	    v1=length_unit*ligp->pos.x;
	    v2=length_unit*ligp->pos.y;
	    v3=length_unit*ligp->pos.z;
	    fprintf(lig_header,"object ligand_mat_%d ligand_glyph\n",state);
	    fprintf(lig_header,"translate %.6g %.6g %.6g\n",v1,v2,v3);
	  }

  	  lig_count++;
	  ligp=ligp->next_ligand;
        }
      }
    }

    if (n_ligand_types>0) {
      fclose(lig_header);
    }
  }

  no_printf("Done traversing.\n");
  fflush(log_file);
  return(0);
}

int output_radiance_objects(struct frame_data_list *fdlp)
{
FILE *wall_header, *wall_verts;
FILE *eff_header, *eff_states, *eff_pos;
FILE *lig_header, *lig_states, *lig_pos;
struct viz_obj *vizp;
struct cmprt_data_list *cdlp;
struct cmprt_data *cdp;
struct wall_list *wlp;
struct wall *wp;
struct effector *rp;
struct rx *rxp;
struct ligand_info *ligip;
struct ligand *ligp;
struct vector3 glyph_axis,rotation_axis;
struct vector3 ab,bc,step_u,step_v,diagonal;
struct vector3 p0,p1,p2;
double i,j,j_max,i_max,i1,i2,i3,j1,j2,j3;
double rad,dot,length,rotation_angle,final_angle,v_val;
float v1,v2,v3,n1,n2,n3;
int ii;
int uu,vv,vv_max,grid_size;
int num;
int viz_iteration,n_viz_iterations;
int first_viz_iteration;
int wall_index,eff_index,lig_index,r_index,orient;
int pos_count,state_count,lig_count,eff_count;
int state;
int viz_type;
unsigned int index,n_eff;
byte viz_eff,viz_lig,viz_surf;
byte viz_surf_or_eff;
byte viz_eff_pos,viz_eff_states;
byte viz_lig_pos,viz_lig_states;
byte viz_surf_pos,viz_surf_states;
byte grid_shape;
char file_name[1024];

viz_iteration=fdlp->viz_iteration;
n_viz_iterations=fdlp->n_viz_iterations;
first_viz_iteration=(viz_iteration==fdlp->iteration_list->value);

viz_type=fdlp->type;
viz_eff=((viz_type==ALL_FRAME_DATA) || (viz_type==EFF_POS)
             || (viz_type==EFF_STATES));
viz_lig=((viz_type==ALL_FRAME_DATA) || (viz_type==MOL_POS)
             || (viz_type==MOL_STATES));
viz_surf=((viz_type==ALL_FRAME_DATA) || (viz_type==SURF_POS)
             || (viz_type==SURF_STATES));
viz_surf_or_eff=(viz_surf || viz_eff);


viz_eff_pos=((viz_type==ALL_FRAME_DATA) || (viz_type==EFF_POS));
viz_eff_states=((viz_type==ALL_FRAME_DATA) || (viz_type==EFF_STATES));
viz_lig_pos=((viz_type==ALL_FRAME_DATA) || (viz_type==MOL_POS));
viz_lig_states=((viz_type==ALL_FRAME_DATA) || (viz_type==MOL_STATES));
viz_surf_pos=((viz_type==ALL_FRAME_DATA) || (viz_type==SURF_POS));
viz_surf_states=((viz_type==ALL_FRAME_DATA) || (viz_type==SURF_STATES));


rad=360/(2*MY_PI);

glyph_axis.x = 0;
glyph_axis.y = 0;
glyph_axis.z = 1;

eff_count=0;

/* dump walls and effectors: */
  if (viz_surf_or_eff) {
    vizp = viz_obj_head;
    while(vizp!=NULL) {
      if (viz_surf) {
        sprintf(file_name,"%s.surfaces.%d.tmesh",vizp->name,viz_iteration);
        if ((wall_header=fopen(file_name,"w"))==NULL) {
          fprintf(log_file,"MCell: error cannot open radiance surface file %s\n",file_name);
          return(1);
        }
        fprintf(wall_header,"# surfaces\n");
      }

      if (viz_eff) {
        sprintf(file_name,"%s.effectors.%d.rad",vizp->name,viz_iteration);
        if ((eff_header=fopen(file_name,"w"))==NULL) {
          fprintf(log_file,"MCell: error cannot open effector header file %s\n",file_name);
          return(1);
        }
        fprintf(eff_header,"# effectors\n");
      }

      wall_index = 3;
      eff_index = 2;
      pos_count = 0;
      state_count = 0;
      cdlp = vizp->cmprt_data_list;
      while(cdlp!=NULL) {
        cdp = cdlp->cmprt_data;
        no_printf("Traversing surfaces in compartment %s\n",cdp->sym->name);
        fflush(log_file);
        wlp = cdp->wall_list;
        while (wlp!=NULL) {
          wp = wlp->wall;
 
          /* dump the effectors */
          if (viz_eff) {
            rp = wp->effectors;
            if (rp!=NULL) {
              n_eff=0;
              for (index=0;index<rp->n_tiles;index++) {
	        rxp=rp->tiles[index];
	        if (rxp!=NULL) {
	          state=rxp->viz_state;
	        }
	        else {
	          state=EXCLUDE_OBJ;
	        }
                if (state!=EXCLUDE_OBJ) {
                  n_eff++;
                }
              }
              no_printf("Dumping %d effectors...\n",n_eff);
              fflush(log_file);
              grid_shape=rp->grid_shape;
              grid_size=rp->grid_size;
              vectorize(wp->vert[0],wp->vert[1],&ab);
              vectorize(wp->vert[1],wp->vert[2],&bc);
              step_u.x=ab.x/grid_size;
              step_u.y=ab.y/grid_size;
              step_u.z=ab.z/grid_size;
              step_v.x=bc.x/grid_size;
              step_v.y=bc.y/grid_size;
              step_v.z=bc.z/grid_size;
              vectorize(&step_u,&step_v,&diagonal);
  
	      cross_prod(&glyph_axis,&(wp->normal),&rotation_axis);
	      dot=dot_prod(&glyph_axis,&(wp->normal))*0.999999999999;
	      rotation_angle = acos(dot)*rad;
/*
	      length = vect_length(&rotation_axis)*0.999999999999;
	      rotation_angle = asin(length)*rad;
*/
  
              index=0;
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
	          rxp=rp->tiles[index];
	          if (rxp!=NULL) {
	            state=rxp->viz_state;
	            orient=(rp->orient[rxp->rx_index]-1)/-2;
	          }
	          else {
	            state=EXCLUDE_OBJ;
	          }
      
	          if (state!=EXCLUDE_OBJ) {
	            final_angle=rotation_angle-(orient*180);
	            v1=length_unit*(wp->vert[0]->x+p2.x);
	            v2=length_unit*(wp->vert[0]->y+p2.y);
	            v3=length_unit*(wp->vert[0]->z+p2.z);
                    eff_count++;
	            fprintf(eff_header,"effector_mat_%d instance effector_glyph.%d\n",state,eff_count);
                    final_angle=0;
		    if (final_angle!=0) {
	              fprintf(eff_header,"rotate %.6g %.6g %.6g %.6g\n",
                        rotation_axis.x,rotation_axis.y,rotation_axis.z,
                        final_angle);
	            }
                    else {
                      fprintf(eff_header,"5 effector_glyph.oct ");
                    }
	            fprintf(eff_header,"-t %.6g %.6g %.6g\n0\n0\n",v1,v2,v3);
	          }
                  index++;
                }
              }
            }
          }

          /* dump the wall */
          if (viz_surf) {
	    if (wp->viz_state!=EXCLUDE_OBJ) {
              if (wp->n_vert==3) {
                fprintf(wall_header,"m surface_mat_%d\n",wp->viz_state);
                for (ii=0;ii<wp->n_vert;ii++) {
                  v1 = length_unit*wp->vert[ii]->x;
                  v2 = length_unit*wp->vert[ii]->y;
                  v3 = length_unit*wp->vert[ii]->z;
                  fprintf(wall_header,"v %d %.6g %.6g %.6g",ii,v1,v2,v3);
                  if (wp->vert_normal!=NULL) {
                    n1 = wp->vert_normal[ii]->x;
                    n2 = wp->vert_normal[ii]->y;
                    n3 = wp->vert_normal[ii]->z;
                    fprintf(wall_header," n %.6g %.6g %.6g",n1,n2,n3);
                  }
                  fprintf(wall_header,"\n");
                  fprintf(wall_header,"t 0 1 2\n");
                }
              }
              else {
/*
                fprintf(wall_header,"poly surface_mat_%d\n",wp->viz_state);
                for (ii=0;ii<wp->n_vert;ii++) {
                  v1 = length_unit*wp->vert[ii]->x;
                  v2 = length_unit*wp->vert[ii]->y;
                  v3 = length_unit*wp->vert[ii]->z;
                  fprintf(wall_header,"%.6g %.6g %.6g\n",v1,v2,v3);
                }
*/
              }
	    }
          }
          wlp = wlp->next;
        }

        cdlp = cdlp->next;
      }

      if (viz_surf) {
        fclose(wall_header);
      }
      if (viz_eff) {
        fclose(eff_header);
      }
      vizp = vizp->next;
    }
  }

/* dump ligands: */
  if (viz_lig) {
    if (n_ligand_types>0) {
      sprintf(file_name,"ligands.%d.rad",viz_iteration);
      if ((lig_header=fopen(file_name,"w"))==NULL) {
        fprintf(log_file,"MCell: error cannot open ligand header file %s\n",file_name);
        return(1);
      }
  
      fprintf(lig_header,"# ligands\n");
    }

    lig_index=2;
    pos_count=0;
    state_count=0;
    for (ii=1;ii<1+n_ligand_types;ii++) {
      ligip=ligand_table[ii];
      ligp=ligip->top;
      num=ligip->lig_index;
      lig_count=0;
      if (num>0) {
        while (ligp!=NULL) {
  
          rp=ligp->effector;
          if (rp!=NULL) {
            r_index=ligp->index;
            state=rp->tiles[r_index]->viz_state;
          }
	  else {
	    state=ligip->viz_state;
	  }

	  if (state!=EXCLUDE_OBJ) {
	    v1=length_unit*ligp->pos.x;
	    v2=length_unit*ligp->pos.y;
	    v3=length_unit*ligp->pos.z;
	    fprintf(lig_header,"ligand_mat_%d instance ligand_glyph\n",state);
	    fprintf(lig_header,"5 ligand_glyph.oct -t %.6g %.6g %.6g\n",
              v1,v2,v3);
	  }

  	  lig_count++;
	  ligp=ligp->next_ligand;
        }
      }
    }

    if (n_ligand_types>0) {
      fclose(lig_header);
    }
  }

  no_printf("Done traversing.\n");
  fflush(log_file);
  return(0);
}


int output_povray_objects(struct frame_data_list *fdlp)
{
FILE *wall_header, *wall_verts;
FILE *eff_header, *eff_states, *eff_pos;
FILE *lig_header, *lig_states, *lig_pos;
struct viz_obj *vizp;
struct cmprt_data_list *cdlp;
struct cmprt_data *cdp;
struct wall_list *wlp;
struct wall *wp;
struct effector *rp;
struct rx *rxp;
struct ligand_info *ligip;
struct ligand *ligp;
struct vector3 glyph_axis,rotation_axis,scale,translate;
struct vector3 ab,bc,step_u,step_v,diagonal;
struct vector3 p0,p1,p2;
struct state_list *surf_state_p,*eff_state_p,*lig_state_p;
double i,j,j_max,i_max,i1,i2,i3,j1,j2,j3;
double rad,dot,axis_length,length,rotation_angle,final_angle,v_val;
double t_matrix[4][4];
float v1,v2,v3,n1,n2,n3;
int ii,jj;
int uu,vv,vv_max,grid_size;
int num;
int viz_iteration,n_viz_iterations;
int first_viz_iteration;
int wall_index,eff_index,lig_index,r_index,orient;
int pos_count,state_count,lig_count;
int state,surf_state,eff_state,lig_state;
int viz_type;
int vert_order[2][3];
unsigned int index,n_eff;
byte viz_eff,viz_lig,viz_surf;
byte viz_surf_or_eff;
byte viz_eff_pos,viz_eff_states;
byte viz_lig_pos,viz_lig_states;
byte viz_surf_pos,viz_surf_states;
byte grid_shape;
char file_name[1024];

viz_iteration=fdlp->viz_iteration;
n_viz_iterations=fdlp->n_viz_iterations;
first_viz_iteration=(viz_iteration==fdlp->iteration_list->value);

viz_type=fdlp->type;
viz_eff=((viz_type==ALL_FRAME_DATA) || (viz_type==EFF_POS)
             || (viz_type==EFF_STATES));
viz_lig=((viz_type==ALL_FRAME_DATA) || (viz_type==MOL_POS)
             || (viz_type==MOL_STATES));
viz_surf=((viz_type==ALL_FRAME_DATA) || (viz_type==SURF_POS)
             || (viz_type==SURF_STATES));
viz_surf_or_eff=(viz_surf || viz_eff);


viz_eff_pos=((viz_type==ALL_FRAME_DATA) || (viz_type==EFF_POS));
viz_eff_states=((viz_type==ALL_FRAME_DATA) || (viz_type==EFF_STATES));
viz_lig_pos=((viz_type==ALL_FRAME_DATA) || (viz_type==MOL_POS));
viz_lig_states=((viz_type==ALL_FRAME_DATA) || (viz_type==MOL_STATES));
viz_surf_pos=((viz_type==ALL_FRAME_DATA) || (viz_type==SURF_POS));
viz_surf_states=((viz_type==ALL_FRAME_DATA) || (viz_type==SURF_STATES));


rad=360/(2*MY_PI);

glyph_axis.x = 0;
glyph_axis.y = 0;
glyph_axis.z = 1;

vert_order[0][0] = 0;
vert_order[0][1] = 1;
vert_order[0][2] = 2;
    
vert_order[1][0] = 0;
vert_order[1][1] = 2;
vert_order[1][2] = 3;

scale.x=1;
scale.y=1;
scale.z=1;

/* dump walls and effectors: */
  if (viz_surf_or_eff) {
    vizp = viz_obj_head;
    while(vizp!=NULL) {
      if (viz_surf) {
        sprintf(file_name,"%s.surfaces.%d.pov",vizp->name,viz_iteration);
        if ((wall_header=fopen(file_name,"w"))==NULL) {
          fprintf(log_file,"MCell: error cannot open povray surface file %s\n",file_name);
          return(1);
        }
        fprintf(wall_header,"/* surfaces in povray-3.1 format */\n");
        fprintf(wall_header,"#include \"surface_materials.pov\"\n");
      }

      if (viz_eff) {
        sprintf(file_name,"%s.effectors.%d.pov",vizp->name,viz_iteration);
        if ((eff_header=fopen(file_name,"w"))==NULL) {
          fprintf(log_file,"MCell: error cannot open effector header file %s\n",file_name);
          return(1);
        }
        fprintf(eff_header,"/* effectors in povray-3.1 format */\n");
        fprintf(eff_header,"#include \"effector_materials.pov\"\n");
        fprintf(eff_header,"#include \"effector_glyph.pov\"\n");
      }

      wall_index = 3;
      eff_index = 2;
      pos_count = 0;
      state_count = 0;
      surf_state_p=surf_state_head;
      eff_state_p=eff_state_head;
      while (surf_state_p!=NULL || eff_state_p!=NULL) { 
        if (surf_state_p!=NULL) {
          surf_state=surf_state_p->state;
          if (viz_surf) {
/*
	    fprintf(wall_header,
              "ReadArchive \"surface_mat_%d.pov\"\n",surf_state);
*/
          }
        }
        else {
          surf_state=EXCLUDE_OBJ;
        }
        if (eff_state_p!=NULL) { 
          eff_state=eff_state_p->state;
          if (viz_eff) {
/*
	    fprintf(eff_header,
              "ReadArchive \"effector_mat_%d.pov\"\n",eff_state);
*/
          }
        }
        else {
          eff_state=EXCLUDE_OBJ;
        }
      cdlp = vizp->cmprt_data_list;
      while(cdlp!=NULL) {
        cdp = cdlp->cmprt_data;
        no_printf("Traversing surfaces in compartment %s\n",cdp->sym->name);
        fflush(log_file);
        wlp = cdp->wall_list;
        while (wlp!=NULL) {
          wp = wlp->wall;
 
          /* dump the effectors */
          if (viz_eff) {
            rp = wp->effectors;
            if (rp!=NULL) {
              n_eff=0;
              for (index=0;index<rp->n_tiles;index++) {
	        rxp=rp->tiles[index];
	        if (rxp!=NULL) {
	          state=rxp->viz_state;
	        }
	        else {
	          state=EXCLUDE_OBJ;
	        }
                if (state!=EXCLUDE_OBJ) {
                  n_eff++;
                }
              }
              no_printf("Dumping %d effectors...\n",n_eff);
              fflush(log_file);
              grid_shape=rp->grid_shape;
              grid_size=rp->grid_size;
              vectorize(wp->vert[0],wp->vert[1],&ab);
              vectorize(wp->vert[1],wp->vert[2],&bc);
              step_u.x=ab.x/grid_size;
              step_u.y=ab.y/grid_size;
              step_u.z=ab.z/grid_size;
              step_v.x=bc.x/grid_size;
              step_v.y=bc.y/grid_size;
              step_v.z=bc.z/grid_size;
              vectorize(&step_u,&step_v,&diagonal);

	      cross_prod(&glyph_axis,&(wp->normal),&rotation_axis);
	      dot=dot_prod(&glyph_axis,&(wp->normal))*0.999999999999;
	      rotation_angle = acos(dot)*rad;
	      axis_length = vect_length(&rotation_axis);
/*
	      length = vect_length(&rotation_axis)*0.999999999999;
	      rotation_angle = asin(length)*rad;
*/
  
              index=0;
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
	          rxp=rp->tiles[index];
	          if (rxp!=NULL) {
	            state=rxp->viz_state;
	            orient=(rp->orient[rxp->rx_index]-1)/-2;
	          }
	          else {
	            state=EXCLUDE_OBJ;
	          }
      
	          if (state==eff_state && state!=EXCLUDE_OBJ) {
	            final_angle=rotation_angle-(orient*180);
	            v1=length_unit*(wp->vert[0]->x+p2.x);
	            v2=length_unit*(wp->vert[0]->y+p2.y);
	            v3=length_unit*(wp->vert[0]->z+p2.z);
                    translate.x=v1;
                    translate.y=v2;
                    translate.z=v3;
                    fprintf(eff_header,"object{effector_glyph ");
                    fprintf(eff_header,"texture{effector_mat_%d}\n",state);
		    if (final_angle!=0 && axis_length>0) {
                      tform_matrix(&scale,&translate,&rotation_axis,
                        final_angle,t_matrix);
                      fprintf(eff_header,
                        "matrix<%.9g,%.9g,%.9g,\n%.9g,%.9g,%.9g,\n%.9g,%.9g,%.9g,\n%.9g,%.9g,%.9g>no_shadow}\n",
                        t_matrix[0][0],t_matrix[0][1],t_matrix[0][2],
                        t_matrix[1][0],t_matrix[1][1],t_matrix[1][2],
                        t_matrix[2][0],t_matrix[2][1],t_matrix[2][2],
                        t_matrix[3][0],t_matrix[3][1],t_matrix[3][2]);
	            }
                    else {
                      fprintf(eff_header,
                        "translate<%.9g,%.9g,%.9g>no_shadow}\n",v1,v2,v3);
                    }
	          }
                  index++;
                }
              }
            }
          }

          /* dump the wall */
          if (viz_surf) {
            state=wp->viz_state;
	    if (state==surf_state && state!=EXCLUDE_OBJ) {
              if (wp->wall_shape==RECT_POLY) {
                for (ii=0;ii<2;ii++) {
                  fprintf(wall_header,"triangle {\n");
                  for (jj=0;jj<3;jj++) {
                    if (jj>0) {
                      fprintf(wall_header,",\n");
                    }
                    v1 = length_unit*wp->vert[vert_order[ii][jj]]->x;
                    v2 = length_unit*wp->vert[vert_order[ii][jj]]->y;
                    v3 = length_unit*wp->vert[vert_order[ii][jj]]->z;
                    fprintf(wall_header,"<%.9g,%.9g,%.9g>",v1,v2,v3);
                  }
                  fprintf(wall_header,"\ntexture{surface_mat_%d}no_shadow}\n",
                          wp->viz_state);
                }
              }
              else {
                if (wp->vert_normal==NULL) {
                  fprintf(wall_header,"triangle {\n");
                  for (jj=0;jj<3;jj++) {
                    if (jj>0) {
                      fprintf(wall_header,",\n");
                    }
                    v1 = length_unit*wp->vert[vert_order[0][jj]]->x;
                    v2 = length_unit*wp->vert[vert_order[0][jj]]->y;
                    v3 = length_unit*wp->vert[vert_order[0][jj]]->z;
                    fprintf(wall_header,"<%.9g,%.9g,%.9g>",v1,v2,v3);
                  }
                }
                else {
                  fprintf(wall_header,"smooth_triangle {\n");
                  for (jj=0;jj<3;jj++) {
                    if (jj>0) {
                      fprintf(wall_header,",\n");
                    }
                    v1 = length_unit*wp->vert[vert_order[0][jj]]->x;
                    v2 = length_unit*wp->vert[vert_order[0][jj]]->y;
                    v3 = length_unit*wp->vert[vert_order[0][jj]]->z;
                    n1 = wp->vert_normal[vert_order[0][jj]]->x;
                    n2 = wp->vert_normal[vert_order[0][jj]]->y;
                    n3 = wp->vert_normal[vert_order[0][jj]]->z;
                    fprintf(wall_header,"<%.9g,%.9g,%.9g>,<%.9g,%.9g,%.9g>",
                      v1,v2,v3,n1,n2,n3);
                  }
                }
                fprintf(wall_header,"\ntexture{surface_mat_%d}no_shadow}\n",
                        wp->viz_state);
              }
	    }
          }

          wlp = wlp->next;
        }

        cdlp = cdlp->next;
      }
        if (surf_state_p!=NULL) {
          surf_state_p=surf_state_p->next;
        }
        if (eff_state_p!=NULL) { 
          eff_state_p=eff_state_p->next;
        }
      }

      if (viz_surf) {
        fclose(wall_header);
      }
      if (viz_eff) {
        fclose(eff_header);
      }
      vizp = vizp->next;
    }
  }

/* dump ligands: */
  if (viz_lig) {
    if (n_ligand_types>0) {
      sprintf(file_name,"ligands.%d.pov",viz_iteration);
      if ((lig_header=fopen(file_name,"w"))==NULL) {
        fprintf(log_file,"MCell: error cannot open ligand header file %s\n",file_name);
        return(1);
      }
      fprintf(lig_header,"/* ligands in povray-3.1 format */\n");
      fprintf(lig_header,"#include \"ligand_materials.pov\"\n");
      fprintf(lig_header,"#include \"ligand_glyph.pov\"\n");
    }

    lig_index=2;
    pos_count=0;
    state_count=0;
    lig_state_p=lig_state_head;
    while (lig_state_p!=NULL) { 
      lig_state=lig_state_p->state;
      for (ii=1;ii<1+n_ligand_types;ii++) {
        ligip=ligand_table[ii];
        ligp=ligip->top;
        num=ligip->lig_index;
        lig_count=0;
        if (num>0) {
          while (ligp!=NULL) {
    
            rp=ligp->effector;
            if (rp!=NULL) {
              r_index=ligp->index;
              state=rp->tiles[r_index]->viz_state;
            }
	    else {
	      state=ligip->viz_state;
	    }
  
	    if (state==lig_state && state!=EXCLUDE_OBJ) {
	      v1=length_unit*ligp->pos.x;
	      v2=length_unit*ligp->pos.y;
	      v3=length_unit*ligp->pos.z;
              fprintf(lig_header,"object{ligand_glyph ");
              fprintf(lig_header,"texture{ligand_mat_%d}\n",state);
              fprintf(lig_header,"translate<%.9g,%.9g,%.9g>no_shadow}\n",v1,v2,v3);
	    }
  
  	    lig_count++;
	    ligp=ligp->next_ligand;
          }
        }
      }
      lig_state_p=lig_state_p->next;
    }

    if (n_ligand_types>0) {
      fclose(lig_header);
    }
  }

  no_printf("Done traversing.\n");
  fflush(log_file);
  return(0);
}


int output_renderman_objects(struct frame_data_list *fdlp)
{
FILE *wall_header, *wall_verts;
FILE *eff_header, *eff_states, *eff_pos;
FILE *lig_header, *lig_states, *lig_pos;
struct viz_obj *vizp;
struct cmprt_data_list *cdlp;
struct cmprt_data *cdp;
struct wall_list *wlp;
struct wall *wp;
struct effector *rp;
struct rx *rxp;
struct ligand_info *ligip;
struct ligand *ligp;
struct vector3 glyph_axis,rotation_axis;
struct vector3 ab,bc,step_u,step_v,diagonal;
struct vector3 p0,p1,p2;
struct state_list *surf_state_p,*eff_state_p,*lig_state_p;
double i,j,j_max,i_max,i1,i2,i3,j1,j2,j3;
double rad,dot,axis_length,rotation_angle,final_angle,v_val,eff_bbox;
float v1,v2,v3,n1,n2,n3;
int ii;
int uu,vv,vv_max,grid_size;
int num;
int viz_iteration,n_viz_iterations;
int first_viz_iteration;
int wall_index,eff_index,lig_index,r_index,orient;
int pos_count,state_count,lig_count;
int state,surf_state,eff_state,lig_state;
int viz_type;
unsigned int index,n_eff;
byte viz_eff,viz_lig,viz_surf;
byte viz_surf_or_eff;
byte viz_eff_pos,viz_eff_states;
byte viz_lig_pos,viz_lig_states;
byte viz_surf_pos,viz_surf_states;
byte grid_shape;
char file_name[1024];

viz_iteration=fdlp->viz_iteration;
n_viz_iterations=fdlp->n_viz_iterations;
first_viz_iteration=(viz_iteration==fdlp->iteration_list->value);

viz_type=fdlp->type;
viz_eff=((viz_type==ALL_FRAME_DATA) || (viz_type==EFF_POS)
             || (viz_type==EFF_STATES));
viz_lig=((viz_type==ALL_FRAME_DATA) || (viz_type==MOL_POS)
             || (viz_type==MOL_STATES));
viz_surf=((viz_type==ALL_FRAME_DATA) || (viz_type==SURF_POS)
             || (viz_type==SURF_STATES));
viz_surf_or_eff=(viz_surf || viz_eff);


viz_eff_pos=((viz_type==ALL_FRAME_DATA) || (viz_type==EFF_POS));
viz_eff_states=((viz_type==ALL_FRAME_DATA) || (viz_type==EFF_STATES));
viz_lig_pos=((viz_type==ALL_FRAME_DATA) || (viz_type==MOL_POS));
viz_lig_states=((viz_type==ALL_FRAME_DATA) || (viz_type==MOL_STATES));
viz_surf_pos=((viz_type==ALL_FRAME_DATA) || (viz_type==SURF_POS));
viz_surf_states=((viz_type==ALL_FRAME_DATA) || (viz_type==SURF_STATES));


rad=360/(2*MY_PI);
eff_bbox=0.5*length_unit;

glyph_axis.x = 0;
glyph_axis.y = 0;
glyph_axis.z = 1;

/* dump walls and effectors: */
  if (viz_surf_or_eff) {
    vizp = viz_obj_head;
    while(vizp!=NULL) {
      if (viz_surf) {
        sprintf(file_name,"%s.surfaces.%d.rib",vizp->name,viz_iteration);
        if ((wall_header=fopen(file_name,"w"))==NULL) {
          fprintf(log_file,"MCell: error cannot open renderman surface file %s\n",file_name);
          return(1);
        }
        fprintf(wall_header,"# surfaces\n");
      }

      if (viz_eff) {
        sprintf(file_name,"%s.effectors.%d.rib",vizp->name,viz_iteration);
        if ((eff_header=fopen(file_name,"w"))==NULL) {
          fprintf(log_file,"MCell: error cannot open effector header file %s\n",file_name);
          return(1);
        }
        fprintf(eff_header,"# effectors\n");
      }

      wall_index = 3;
      eff_index = 2;
      pos_count = 0;
      state_count = 0;
      surf_state_p=surf_state_head;
      eff_state_p=eff_state_head;
      while (surf_state_p!=NULL || eff_state_p!=NULL) { 
        if (surf_state_p!=NULL) {
          surf_state=surf_state_p->state;
          if (viz_surf) {
	    fprintf(wall_header,
              "ReadArchive \"surface_mat_%d.rib\"\n",surf_state);
          }
        }
        else {
          surf_state=EXCLUDE_OBJ;
        }
        if (eff_state_p!=NULL) { 
          eff_state=eff_state_p->state;
          if (viz_eff) {
	    fprintf(eff_header,
              "ReadArchive \"effector_mat_%d.rib\"\n",eff_state);
          }
        }
        else {
          eff_state=EXCLUDE_OBJ;
        }
      cdlp = vizp->cmprt_data_list;
      while(cdlp!=NULL) {
        cdp = cdlp->cmprt_data;
        no_printf("Traversing surfaces in compartment %s\n",cdp->sym->name);
        fflush(log_file);
        wlp = cdp->wall_list;
        while (wlp!=NULL) {
          wp = wlp->wall;
 
          /* dump the effectors */
          if (viz_eff) {
            rp = wp->effectors;
            if (rp!=NULL) {
              n_eff=0;
              for (index=0;index<rp->n_tiles;index++) {
	        rxp=rp->tiles[index];
	        if (rxp!=NULL) {
	          state=rxp->viz_state;
	        }
	        else {
	          state=EXCLUDE_OBJ;
	        }
                if (state!=EXCLUDE_OBJ) {
                  n_eff++;
                }
              }
              no_printf("Dumping %d effectors...\n",n_eff);
              fflush(log_file);
              grid_shape=rp->grid_shape;
              grid_size=rp->grid_size;
              vectorize(wp->vert[0],wp->vert[1],&ab);
              vectorize(wp->vert[1],wp->vert[2],&bc);
              step_u.x=ab.x/grid_size;
              step_u.y=ab.y/grid_size;
              step_u.z=ab.z/grid_size;
              step_v.x=bc.x/grid_size;
              step_v.y=bc.y/grid_size;
              step_v.z=bc.z/grid_size;
              vectorize(&step_u,&step_v,&diagonal);
  
	      cross_prod(&glyph_axis,&(wp->normal),&rotation_axis);
	      dot=dot_prod(&glyph_axis,&(wp->normal))*0.999999999999;
	      rotation_angle = acos(dot)*rad;
	      axis_length = vect_length(&rotation_axis);
/*
	      axis_length = vect_length(&rotation_axis)*0.999999999999;
	      rotation_angle = asin(axis_length)*rad;
*/
  
              index=0;
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
	          rxp=rp->tiles[index];
	          if (rxp!=NULL) {
	            state=rxp->viz_state;
	            orient=(rp->orient[rxp->rx_index]-1)/-2;
	          }
	          else {
	            state=EXCLUDE_OBJ;
	          }
      
	          if (state==eff_state && state!=EXCLUDE_OBJ) {
	            final_angle=rotation_angle-(orient*180);
	            v1=length_unit*(wp->vert[0]->x+p2.x);
	            v2=length_unit*(wp->vert[0]->y+p2.y);
	            v3=length_unit*(wp->vert[0]->z+p2.z);
	            fprintf(eff_header,"AttributeBegin\n");
	            fprintf(eff_header,"Translate %.9g %.9g %.9g\n",v1,v2,v3);
		    if (final_angle!=0 && axis_length>0) {
	              fprintf(eff_header,"Rotate %.9g %.9g %.9g %.9g\n",final_angle,
	                rotation_axis.x,rotation_axis.y,rotation_axis.z);
	            }
/*
	            fprintf(eff_header,"ObjectInstance 1\n");
*/
                    fprintf(eff_header,"Procedural \"DelayedReadArchive\" [ \"effector_glyph_%d.rib\" ] [ %.4g %.4g %.4g %.4g %.4g %.4g ]\n",
                      state,-eff_bbox,eff_bbox,-eff_bbox,eff_bbox,
                      -eff_bbox,eff_bbox);
	            fprintf(eff_header,"AttributeEnd\n");
	          }
                  index++;
                }
              }
            }
          }

          /* dump the wall */
          if (viz_surf) {
            state=wp->viz_state;
	    if (state==surf_state && state!=EXCLUDE_OBJ) {
	      fprintf(wall_header,"AttributeBegin\n");
	      fprintf(wall_header,"Polygon\n");
	      fprintf(wall_header,"  \"P\" [");
              for (ii=0;ii<wp->n_vert;ii++) {
                v1 = length_unit*wp->vert[ii]->x;
                v2 = length_unit*wp->vert[ii]->y;
                v3 = length_unit*wp->vert[ii]->z;
                fprintf(wall_header," %.9g %.9g %.9g ",v1,v2,v3);
              }
	      fprintf(wall_header,"  ]\n");
              if (wp->vert_normal!=NULL) {
	        fprintf(wall_header,"  \"N\" [");
                for (ii=0;ii<wp->n_vert;ii++) {
                    n1 = wp->vert_normal[ii]->x;
                    n2 = wp->vert_normal[ii]->y;
                    n3 = wp->vert_normal[ii]->z;
                    fprintf(wall_header," %.9g %.9g %.9g ",n1,n2,n3);
                }
	        fprintf(wall_header,"  ]\n");
              }
	      fprintf(wall_header,"AttributeEnd\n");
	    }
          }

          wlp = wlp->next;
        }

        cdlp = cdlp->next;
      }
        if (surf_state_p!=NULL) {
          surf_state_p=surf_state_p->next;
        }
        if (eff_state_p!=NULL) { 
          eff_state_p=eff_state_p->next;
        }
      }

      if (viz_surf) {
        fclose(wall_header);
      }
      if (viz_eff) {
        fclose(eff_header);
      }
      vizp = vizp->next;
    }
  }

/* dump ligands: */
  if (viz_lig) {
    if (n_ligand_types>0) {
      sprintf(file_name,"ligands.%d.rib",viz_iteration);
      if ((lig_header=fopen(file_name,"w"))==NULL) {
        fprintf(log_file,"MCell: error cannot open ligand header file %s\n",file_name);
        return(1);
      }
  
      fprintf(lig_header,"# ligands\n");
    }

    lig_index=2;
    pos_count=0;
    state_count=0;
    lig_state_p=lig_state_head;
    while (lig_state_p!=NULL) { 
      lig_state=lig_state_p->state;
      fprintf(lig_header,"ReadArchive \"ligand_mat_%d.rib\"\n",lig_state);
      for (ii=1;ii<1+n_ligand_types;ii++) {
        ligip=ligand_table[ii];
        ligp=ligip->top;
        num=ligip->lig_index;
        lig_count=0;
        if (num>0) {
          while (ligp!=NULL) {
    
            rp=ligp->effector;
            if (rp!=NULL) {
              r_index=ligp->index;
              state=rp->tiles[r_index]->viz_state;
            }
	    else {
	      state=ligip->viz_state;
	    }
  
	    if (state==lig_state && state!=EXCLUDE_OBJ) {
	      v1=length_unit*ligp->pos.x;
	      v2=length_unit*ligp->pos.y;
	      v3=length_unit*ligp->pos.z;
	      fprintf(lig_header,"AttributeBegin\n");
	      fprintf(lig_header,"Translate %.9g %.9g %.9g\n",v1,v2,v3);
	      fprintf(lig_header,"ObjectInstance 1\n");
	      fprintf(lig_header,"AttributeEnd\n");
	    }
  
  	    lig_count++;
	    ligp=ligp->next_ligand;
          }
        }
      }
      lig_state_p=lig_state_p->next;
    }

    if (n_ligand_types>0) {
      fclose(lig_header);
    }
  }

  no_printf("Done traversing.\n");
  fflush(log_file);
  return(0);
}



int output_irit_objects(struct frame_data_list *fdlp)
{
FILE *wall_header, *wall_verts;
FILE *eff_header, *eff_states, *eff_pos;
FILE *lig_header, *lig_states, *lig_pos;
struct viz_obj *vizp;
struct cmprt_data_list *cdlp;
struct cmprt_data *cdp;
struct wall_list *wlp;
struct wall *wp;
struct effector *rp;
struct rx *rxp;
struct ligand_info *ligip;
struct ligand *ligp;
struct vector3 glyph_axis,rotation_axis;
struct vector3 ab,bc,step_u,step_v,diagonal;
struct vector3 p0,p1,p2;
double rad,length,rotation_angle,final_angle,v_val;
float v1,v2,v3;
int uu,vv,vv_max,grid_size;
int ii;
int num;
int viz_iteration;
int wall_index,eff_index,lig_index,r_index,orient;
int pos_count,state_count,lig_count;
int state;
unsigned int index,n_eff;
byte grid_shape;
char file_name[1024];

viz_iteration=fdlp->viz_iteration;

rad=360/(2*MY_PI);

glyph_axis.x = 0;
glyph_axis.y = 0;
glyph_axis.z = 1;

/* dump walls and effectors: */
  vizp = viz_obj_head;
  while(vizp!=NULL) {
    sprintf(file_name,"%s.surfaces.%d.dat",vizp->name,viz_iteration);
    if ((wall_header=fopen(file_name,"w"))==NULL) {
      fprintf(log_file,"MCell: error cannot open IRIT surface file %s\n",file_name);
      return(1);
    }

    sprintf(file_name,"%s.effectors.%d.dat",vizp->name,viz_iteration);
    if ((eff_header=fopen(file_name,"w"))==NULL) {
      fprintf(log_file,"MCell: error cannot open IRIT effector file %s\n",file_name);
      return(1);
    }

    fprintf(wall_header,"The Surfaces\n");
    fprintf(wall_header,"[OBJECT THE_SURFACES\n");

    fprintf(eff_header,"The Receptors\n");
    fprintf(eff_header,"[OBJECT THE_EFFECTORS\n");

    wall_index = 3;
    eff_index = 2;
    pos_count = 0;
    state_count = 0;
    cdlp = vizp->cmprt_data_list;
    while(cdlp!=NULL) {
      cdp = cdlp->cmprt_data;
      wlp = cdp->wall_list;
      no_printf("Traversing surfaces in compartment %s\n",cdp->sym->name);
      while (wlp!=NULL) {
        wp = wlp->wall;
	rp = wp->effectors;

        /* dump the effectors */
        if (rp!=NULL) {
          n_eff=0;
          for (index=0;index<rp->n_tiles;index++) {
	    rxp=rp->tiles[index];
	    if (rxp!=NULL) {
	      state=rxp->viz_state;
	    }
	    else {
	      state=EXCLUDE_OBJ;
	    }
            if (state!=EXCLUDE_OBJ) {
              n_eff++;
            }
          }
          no_printf("Dumping %d effectors...\n",n_eff);
          fflush(log_file);
          grid_shape=rp->grid_shape;
          grid_size=rp->grid_size;
          if (n_eff>0) {
            fprintf(eff_header,"\t[POINTLIST %d\n",n_eff);
          }
          vectorize(wp->vert[0],wp->vert[1],&ab);
          vectorize(wp->vert[1],wp->vert[2],&bc);
          step_u.x=ab.x/grid_size;
          step_u.y=ab.y/grid_size;
          step_u.z=ab.z/grid_size;
          step_v.x=bc.x/grid_size;
          step_v.y=bc.y/grid_size;
          step_v.z=bc.z/grid_size;
          vectorize(&step_u,&step_v,&diagonal);

	  cross_prod(&glyph_axis,&(wp->normal),&rotation_axis);
	  length = vect_length(&rotation_axis)*0.999999999999;
	  rotation_angle = asin(length)*rad;

          index=0;
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
	      rxp=rp->tiles[index];
	      if (rxp!=NULL) {
	        state=rxp->viz_state;
	        orient=(rp->orient[rxp->rx_index]-1)/-2;
	      }
	      else {
	        state=EXCLUDE_OBJ;
	      }
    
	      if (state!=EXCLUDE_OBJ) {
	        final_angle=rotation_angle-(orient*180);
	        v1=length_unit*(wp->vert[0]->x+p2.x);
	        v2=length_unit*(wp->vert[0]->y+p2.y);
	        v3=length_unit*(wp->vert[0]->z+p2.z);
                fprintf(eff_header,"\t\t[ %g %g %g ]\n",v1,v2,v3);
	      }
              index++;
            }
          }
          if (n_eff>0) {
            fprintf(eff_header,"\t]\n");
          }
        }

        /* dump the wall */
	if (wp->viz_state!=EXCLUDE_OBJ) {
	  fprintf(wall_header,"\t[POLYGON %d\n",wp->n_vert);
          for (ii=0;ii<wp->n_vert;ii++) {
            v1 = length_unit*wp->vert[ii]->x;
            v2 = length_unit*wp->vert[ii]->y;
            v3 = length_unit*wp->vert[ii]->z;
            fprintf(wall_header,"\t\t[ %g %g %g ]\n",v1,v2,v3);
          }
	  fprintf(wall_header,"\t]\n");
	}
        wlp = wlp->next;
      }

      cdlp = cdlp->next;
    }

    fprintf(wall_header,"]\n");
    fprintf(eff_header,"]\n");
    fclose(wall_header);
    fclose(eff_header);
    vizp = vizp->next;
  }

/* dump ligands: */
  if (n_ligand_types>0) {
    sprintf(file_name,"ligands.%d.dat",viz_iteration);
    if ((lig_header=fopen(file_name,"w"))==NULL) {
      fprintf(log_file,"MCell: error cannot open IRIT ligand file %s\n",file_name);
      return(1);
    }

    fprintf(lig_header,"The Ligands\n");
    fprintf(lig_header,"[OBJECT THE_LIGANDS\n");
  }

  lig_index=2;
  pos_count=0;
  state_count=0;
  for (ii=1;ii<1+n_ligand_types;ii++) {
    ligip=ligand_table[ii];
    ligp=ligip->top;
    num=ligip->lig_index;
    lig_count=0;
    if (num>0) {
      while (ligp!=NULL) {

        rp=ligp->effector;
        if (rp!=NULL) {
          r_index=ligp->index;
          state=rp->tiles[r_index]->viz_state;
        }
	else {
	  state=ligip->viz_state;
	}

	if (state!=EXCLUDE_OBJ) {
	  lig_count++;
	}

	ligp=ligp->next_ligand;
      }
    }

    ligp=ligip->top;
    if (lig_count>0) {
      fprintf(lig_header,"\t[POINTLIST %d\n",lig_count);
      while (ligp!=NULL) {

        rp=ligp->effector;
        if (rp!=NULL) {
          r_index=ligp->index;
          state=rp->tiles[r_index]->viz_state;
        }
	else {
	  state=ligip->viz_state;
	}

	if (state!=EXCLUDE_OBJ) {
	  v1=length_unit*ligp->pos.x;
	  v2=length_unit*ligp->pos.y;
	  v3=length_unit*ligp->pos.z;
	  fprintf(lig_header,"\t\t[ %g %g %g ]\n",v1,v2,v3);
	}

	ligp=ligp->next_ligand;
      }
      fprintf(lig_header,"\t]\n");
    }
  }

  if (n_ligand_types>0) {
    fprintf(lig_header,"]\n");
    fclose(lig_header);
  }

  no_printf("IRIT output done.\n");
  return(0);
}


int output_mcell_objects(struct frame_data_list *fdlp)
{
FILE *wall_header, *wall_verts;
FILE *eff_header, *eff_states, *eff_pos;
FILE *lig_header, *lig_states, *lig_pos;
struct viz_obj *vizp;
struct cmprt_data_list *cdlp;
struct cmprt_data *cdp;
struct wall_list *wlp;
struct wall *wp;
struct effector *rp;
struct rx *rxp;
struct ligand_info *ligip;
struct ligand *ligp;
struct vector3 glyph_axis,rotation_axis;
struct vector3 ab,bc,step_u,step_v,diagonal;
struct vector3 p0,p1,p2,et_or,et0,et1,et2;
struct state_list *surf_state_p,*eff_state_p,*lig_state_p;
double i,j,j_max,i_max,i1,i2,i3,j1,j2,j3;
double rad,dot,axis_length,rotation_angle,final_angle,v_val;
float v1,v2,v3,n1,n2,n3;
int ii;
int uu,vv,vv_max,grid_size;
int num;
int viz_iteration,n_viz_iterations;
int first_viz_iteration;
int wall_index,eff_index,lig_index,r_index,orient;
int pos_count,state_count,lig_count;
int state,surf_state,eff_state,lig_state;
int viz_type;
unsigned int index,n_eff;
byte viz_eff,viz_lig,viz_surf;
byte viz_surf_or_eff;
byte viz_eff_pos,viz_eff_states;
byte viz_lig_pos,viz_lig_states;
byte viz_surf_pos,viz_surf_states;
byte grid_shape;
char file_name[1024];

viz_iteration=fdlp->viz_iteration;
n_viz_iterations=fdlp->n_viz_iterations;
first_viz_iteration=(viz_iteration==fdlp->iteration_list->value);

viz_type=fdlp->type;
viz_eff=((viz_type==ALL_FRAME_DATA) || (viz_type==EFF_POS)
             || (viz_type==EFF_STATES));
viz_lig=((viz_type==ALL_FRAME_DATA) || (viz_type==MOL_POS)
             || (viz_type==MOL_STATES));
viz_surf=((viz_type==ALL_FRAME_DATA) || (viz_type==SURF_POS)
             || (viz_type==SURF_STATES));
viz_surf_or_eff=(viz_surf || viz_eff);


viz_eff_pos=((viz_type==ALL_FRAME_DATA) || (viz_type==EFF_POS));
viz_eff_states=((viz_type==ALL_FRAME_DATA) || (viz_type==EFF_STATES));
viz_lig_pos=((viz_type==ALL_FRAME_DATA) || (viz_type==MOL_POS));
viz_lig_states=((viz_type==ALL_FRAME_DATA) || (viz_type==MOL_STATES));
viz_surf_pos=((viz_type==ALL_FRAME_DATA) || (viz_type==SURF_POS));
viz_surf_states=((viz_type==ALL_FRAME_DATA) || (viz_type==SURF_STATES));


rad=360/(2*MY_PI);

glyph_axis.x = 0;
glyph_axis.y = 0;
glyph_axis.z = 1;

/* dump walls and effectors: */
  if (viz_surf_or_eff) {
    vizp = viz_obj_head;
    while(vizp!=NULL) {
      if (viz_surf) {
        sprintf(file_name,"%s.surfaces.%d.viz",vizp->name,viz_iteration);
        if ((wall_header=fopen(file_name,"w"))==NULL) {
          fprintf(log_file,"MCell: error cannot open renderman surface file %s\n",file_name);
          return(1);
        }
      }

      if (viz_eff) {
        sprintf(file_name,"%s.effectors.%d.viz",vizp->name,viz_iteration);
        if ((eff_header=fopen(file_name,"w"))==NULL) {
          fprintf(log_file,"MCell: error cannot open effector header file %s\n",file_name);
          return(1);
        }
      }

      wall_index = 3;
      eff_index = 2;
      pos_count = 0;
      state_count = 0;
      surf_state_p=surf_state_head;
      eff_state_p=eff_state_head;
      while (surf_state_p!=NULL || eff_state_p!=NULL) { 
        if (surf_state_p!=NULL) {
          surf_state=surf_state_p->state;
          if (viz_surf) {
	    fprintf(wall_header,
              "state %d\n",surf_state);
          }
        }
        else {
          surf_state=EXCLUDE_OBJ;
        }
        if (eff_state_p!=NULL) { 
          eff_state=eff_state_p->state;
          if (viz_eff) {
	    fprintf(eff_header,
              "state %d\n",eff_state);
          }
        }
        else {
          eff_state=EXCLUDE_OBJ;
        }
      cdlp = vizp->cmprt_data_list;
      while(cdlp!=NULL) {
        cdp = cdlp->cmprt_data;
        no_printf("Traversing surfaces in compartment %s\n",cdp->sym->name);
        fflush(log_file);
        wlp = cdp->wall_list;
        while (wlp!=NULL) {
          wp = wlp->wall;
 
          /* dump the effectors */
          if (viz_eff) {
            rp = wp->effectors;
            if (rp!=NULL) {
              n_eff=0;
              for (index=0;index<rp->n_tiles;index++) {
	        rxp=rp->tiles[index];
	        if (rxp!=NULL) {
	          state=rxp->viz_state;
	        }
	        else {
	          state=EXCLUDE_OBJ;
	        }
                if (state!=EXCLUDE_OBJ) {
                  n_eff++;
                }
              }
              no_printf("Dumping %d effectors...\n",n_eff);
              fflush(log_file);
              grid_shape=rp->grid_shape;
              grid_size=rp->grid_size;
              vectorize(wp->vert[0],wp->vert[1],&ab);
              vectorize(wp->vert[1],wp->vert[2],&bc);
              step_u.x=ab.x/grid_size;
              step_u.y=ab.y/grid_size;
              step_u.z=ab.z/grid_size;
              step_v.x=bc.x/grid_size;
              step_v.y=bc.y/grid_size;
              step_v.z=bc.z/grid_size;
              vectorize(&step_u,&step_v,&diagonal);
  
	      cross_prod(&glyph_axis,&(wp->normal),&rotation_axis);
	      dot=dot_prod(&glyph_axis,&(wp->normal))*0.999999999999;
	      rotation_angle = acos(dot)*rad;
	      axis_length = vect_length(&rotation_axis);
/*
	      axis_length = vect_length(&rotation_axis)*0.999999999999;
	      rotation_angle = asin(axis_length)*rad;
*/
  
              index=0;
              for (uu=0;uu<grid_size;uu++) {
                p0.x=(uu+1)*step_u.x;
                p0.y=(uu+1)*step_u.y;
                p0.z=(uu+1)*step_u.z;
                et_or.x=uu*step_u.x;
                et_or.y=uu*step_u.y;
                et_or.z=uu*step_u.z;
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
                  et0.x=et_or.x+(v_val*step_v.x);
                  et0.y=et_or.y+(v_val*step_v.y);
                  et0.z=et_or.z+(v_val*step_v.z);
                  if (vv%2==0) {
                    p2.x=p1.x+(0.33*diagonal.x);
                    p2.y=p1.y+(0.33*diagonal.y);
                    p2.z=p1.z+(0.33*diagonal.z);
                    et1.x=et0.x+step_u.x;
                    et1.y=et0.y+step_u.y;
                    et1.z=et0.z+step_u.z;
                    et2.x=et1.x+step_v.x;
                    et2.y=et1.y+step_v.y;
                    et2.z=et1.z+step_v.z;
                  }
                  else {
                    p2.x=p1.x+(0.66*diagonal.x);
                    p2.y=p1.y+(0.66*diagonal.y);
                    p2.z=p1.z+(0.66*diagonal.z);
                    et1.x=et0.x+step_u.x+step_v.x;
                    et1.y=et0.y+step_u.y+step_v.y;
                    et1.z=et0.z+step_u.z+step_v.z;
                    et2.x=et0.x+step_v.x;
                    et2.y=et0.y+step_v.y;
                    et2.z=et0.z+step_v.z;
                    v_val++;
                  }
	          rxp=rp->tiles[index];
	          if (rxp!=NULL) {
	            state=rxp->viz_state;
	            orient=(rp->orient[rxp->rx_index]-1)/-2;
	          }
	          else {
	            state=EXCLUDE_OBJ;
	          }
      
	          if (state==eff_state && state!=EXCLUDE_OBJ) {
	            final_angle=rotation_angle-(orient*180);
	            v1=length_unit*(wp->vert[0]->x+p2.x);
	            v2=length_unit*(wp->vert[0]->y+p2.y);
	            v3=length_unit*(wp->vert[0]->z+p2.z);
	            fprintf(eff_header,"location %.9g %.9g %.9g\n",v1,v2,v3);
	            v1=length_unit*(wp->vert[0]->x+et0.x);
	            v2=length_unit*(wp->vert[0]->y+et0.y);
	            v3=length_unit*(wp->vert[0]->z+et0.z);
	            fprintf(eff_header,"polygon %.9g %.9g %.9g\n",v1,v2,v3);
	            v1=length_unit*(wp->vert[0]->x+et1.x);
	            v2=length_unit*(wp->vert[0]->y+et1.y);
	            v3=length_unit*(wp->vert[0]->z+et1.z);
	            fprintf(eff_header," %.9g %.9g %.9g\n",v1,v2,v3);
	            v1=length_unit*(wp->vert[0]->x+et2.x);
	            v2=length_unit*(wp->vert[0]->y+et2.y);
	            v3=length_unit*(wp->vert[0]->z+et2.z);
	            fprintf(eff_header," %.9g %.9g %.9g\n",v1,v2,v3);
	          }
                  index++;
                }
              }
            }
          }

          /* dump the wall */
          if (viz_surf) {
            state=wp->viz_state;
	    if (state==surf_state && state!=EXCLUDE_OBJ) {
	      fprintf(wall_header,"polygon");
              for (ii=0;ii<wp->n_vert;ii++) {
                v1 = length_unit*wp->vert[ii]->x;
                v2 = length_unit*wp->vert[ii]->y;
                v3 = length_unit*wp->vert[ii]->z;
                fprintf(wall_header," %.9g %.9g %.9g\n",v1,v2,v3);
              }
	    }
          }

          wlp = wlp->next;
        }

        cdlp = cdlp->next;
      }
        if (surf_state_p!=NULL) {
          surf_state_p=surf_state_p->next;
        }
        if (eff_state_p!=NULL) { 
          eff_state_p=eff_state_p->next;
        }
      }

      if (viz_surf) {
        fclose(wall_header);
      }
      if (viz_eff) {
        fclose(eff_header);
      }
      vizp = vizp->next;
    }
  }

/* dump ligands: */
  if (viz_lig) {
    if (n_ligand_types>0) {
      sprintf(file_name,"ligands.%d.viz",viz_iteration);
      if ((lig_header=fopen(file_name,"w"))==NULL) {
        fprintf(log_file,"MCell: error cannot open ligand header file %s\n",file_name);
        return(1);
      }
    }

    lig_index=2;
    pos_count=0;
    state_count=0;
    lig_state_p=lig_state_head;
    while (lig_state_p!=NULL) { 
      lig_state=lig_state_p->state;
      fprintf(lig_header,"state %d\n",lig_state);
      for (ii=1;ii<1+n_ligand_types;ii++) {
        ligip=ligand_table[ii];
        ligp=ligip->top;
        num=ligip->lig_index;
        lig_count=0;
        if (num>0) {
          while (ligp!=NULL) {
    
            rp=ligp->effector;
            if (rp!=NULL) {
              r_index=ligp->index;
              state=rp->tiles[r_index]->viz_state;
            }
	    else {
	      state=ligip->viz_state;
	    }
  
	    if (state==lig_state && state!=EXCLUDE_OBJ) {
	      v1=length_unit*ligp->pos.x;
	      v2=length_unit*ligp->pos.y;
	      v3=length_unit*ligp->pos.z;
	      fprintf(lig_header,"location %.9g %.9g %.9g\n",v1,v2,v3);
	    }
  
  	    lig_count++;
	    ligp=ligp->next_ligand;
          }
        }
      }
      lig_state_p=lig_state_p->next;
    }

    if (n_ligand_types>0) {
      fclose(lig_header);
    }
  }

  no_printf("Done traversing.\n");
  fflush(log_file);
  return(0);
}

int output_voxel_image(struct frame_data_list *fdlp)
{
FILE *wall_header, *wall_verts;
FILE *eff_header, *eff_states, *eff_pos;
FILE *lig_header, *lig_states, *lig_pos;
struct viz_obj *vizp;
struct cmprt_data_list *cdlp;
struct cmprt_data *cdp;
struct wall_list *wlp;
struct wall *wp;
struct effector *rp;
struct rx *rxp;
struct ligand_info *ligip;
struct ligand *ligp;
struct vector3 glyph_axis,rotation_axis;
struct vector3 ab,bc,step_u,step_v,diagonal;
struct vector3 p0,p1,p2;
struct state_list *surf_state_p,*eff_state_p,*lig_state_p;
double rad,length,rotation_angle,final_angle,v_val;
double umin,umax,vmin,vmax,w_val,wmin,wmax,voxel_uv,voxel_w;
double v1,v2,v3;
int iv1,iv2,iv3;
int uu,vv,vv_max,grid_size;
int i,j;
int num;
int viz_iteration;
int viz_type;
int wall_index,eff_index,r_index,orient;
int pos_count,state_count,lig_count,eff_count;
int state,surf_state,eff_state,lig_state;
unsigned int index,n_eff;
unsigned int image_x,image_y,num_pixels,*image_buffer;
byte viz_eff,viz_lig,viz_surf;
byte viz_surf_or_eff;
byte viz_eff_pos,viz_eff_states;
byte viz_lig_pos,viz_lig_states;
byte viz_surf_pos,viz_surf_states;
byte grid_shape;
byte i_val;
char *rx_name;
char file_name[1024];

  no_printf("voxel image output...\n");

  viz_iteration=fdlp->viz_iteration;

  viz_type=fdlp->type;
  viz_eff=((viz_type==ALL_FRAME_DATA) || (viz_type==EFF_POS)
             || (viz_type==EFF_STATES));
  viz_lig=((viz_type==ALL_FRAME_DATA) || (viz_type==MOL_POS)
             || (viz_type==MOL_STATES));
  viz_surf=((viz_type==ALL_FRAME_DATA) || (viz_type==SURF_POS)
             || (viz_type==SURF_STATES));
  viz_surf_or_eff=(viz_surf || viz_eff);


  viz_eff_pos=((viz_type==ALL_FRAME_DATA) || (viz_type==EFF_POS));
  viz_eff_states=((viz_type==ALL_FRAME_DATA) || (viz_type==EFF_STATES));
  viz_lig_pos=((viz_type==ALL_FRAME_DATA) || (viz_type==MOL_POS));
  viz_lig_states=((viz_type==ALL_FRAME_DATA) || (viz_type==MOL_STATES));
  viz_surf_pos=((viz_type==ALL_FRAME_DATA) || (viz_type==SURF_POS));
  viz_surf_states=((viz_type==ALL_FRAME_DATA) || (viz_type==SURF_STATES));

  voxel_w=0.50;
  voxel_uv=0.134;
  umin=-1.072;
  umax=1.072;
  vmin=-1.072;
  vmax=1.072;
/*
  voxel_uv=0.025;
  umin=-2.1;
  umax=2.1;
  vmin=-2.5;
  vmax=0.5;
*/
  w_val=0;
  wmin=w_val-(voxel_w/2);
  wmax=w_val+(voxel_w/2);

  image_x=(unsigned int)((umax-umin)/voxel_uv);
  image_y=(unsigned int)((vmax-vmin)/voxel_uv);
  num_pixels=image_x*image_y;

  if ((image_buffer=(unsigned int *)malloc
                    (num_pixels*sizeof(unsigned int)))==NULL) {
    printf ("MCell: cannot store voxel image\n");
    return(1);
  }

/* voxelize effectors: */
eff_state_p=eff_state_head;
while (eff_state_p!=NULL) { 
  eff_state=eff_state_p->state;
  eff_count=0;

  if (viz_eff) {
    no_printf("Voxelizing effector state %d\n",eff_state);
    rx_name=eff_state_p->name;
    for (j=0;j<num_pixels;j++) {
      image_buffer[j]=0;
    }
    vizp = viz_obj_head;
    while(vizp!=NULL) {


      cdlp = vizp->cmprt_data_list;
      while(cdlp!=NULL) {
        cdp = cdlp->cmprt_data;
        no_printf("Traversing surfaces in compartment %s\n",cdp->sym->name);
        fflush(log_file);
        wlp = cdp->wall_list;
        while (wlp!=NULL) {
          wp = wlp->wall;
 
          /* dump the effectors */
          rp = wp->effectors;
          if (rp!=NULL) {
            n_eff=0;
            for (index=0;index<rp->n_tiles;index++) {
	      rxp=rp->tiles[index];
	      if (rxp!=NULL) {
	        state=rxp->viz_state;
	      }
	      else {
	        state=EXCLUDE_OBJ;
	      }
              if (state==eff_state) {
                rx_name=rxp->sym->name;
                n_eff++;
              }
            }

            no_printf("Voxelizing %d effectors...\n",n_eff);
            if (n_eff>0) {
              fflush(log_file);
              grid_shape=rp->grid_shape;
              grid_size=rp->grid_size;
              vectorize(wp->vert[0],wp->vert[1],&ab);
              vectorize(wp->vert[1],wp->vert[2],&bc);
              step_u.x=ab.x/grid_size;
              step_u.y=ab.y/grid_size;
              step_u.z=ab.z/grid_size;
              step_v.x=bc.x/grid_size;
              step_v.y=bc.y/grid_size;
              step_v.z=bc.z/grid_size;
              vectorize(&step_u,&step_v,&diagonal);
  
              index=0;
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
	          rxp=rp->tiles[index];
	          if (rxp!=NULL) {
	            state=rxp->viz_state;
	          }
	          else {
	            state=EXCLUDE_OBJ;
	          }
      
                  if (state==eff_state && state!=EXCLUDE_OBJ) {
	            v1=length_unit*(wp->vert[0]->x+p2.x);
	            v2=length_unit*(wp->vert[0]->y+p2.y);
	            v3=length_unit*(wp->vert[0]->z+p2.z);
                    if (v1>=umin && v1<=umax) {
                      if (v3>=vmin && v3<=vmax) {
                        if (v2>=wmin && v2<=wmax) {
	                  iv1=(int)((v1-umin)/voxel_uv);
	                  iv3=(image_y-1)-(int)((v3-vmin)/voxel_uv);
                          eff_count++;
                          image_buffer[iv3*image_x+iv1]++;
                        } 
                      } 
                    } 
	          }
                  index++;
                }
              }
            }
          }
          wlp = wlp->next;
        }

        cdlp = cdlp->next;
      }

      vizp = vizp->next;
    }
  }
  no_printf("MCell: writing file voxel_image.effector_%s.state_%d.%d.pgm\n",
    rx_name,eff_state,viz_iteration);
  sprintf(file_name,"voxel_image.effector_%s.state_%d.%d.pgm",
    rx_name,eff_state,viz_iteration);
  if ((eff_header=fopen(file_name,"w"))==NULL) {
    fprintf(log_file,"MCell: error cannot open voxel image file %s\n",
    file_name);
    return(1);
  }
  
  fprintf(eff_header,"P2\n%d %d\n255\n",image_x,image_y);
  for (j=0;j<num_pixels;j++) {
    i_val=image_buffer[j];
    fprintf(eff_header,"%d\n",i_val);
  }
  fclose(eff_header);

  eff_state_p = eff_state_p->next;
}


/* voxelize ligands: */

  lig_state_p=lig_state_head;
  while (lig_state_p!=NULL) { 
    lig_state=lig_state_p->state;
    rx_name=lig_state_p->name;
    for (i=1;i<1+n_ligand_types;i++) {
      for (j=0;j<num_pixels;j++) {
        image_buffer[j]=0;
      }
      ligip=ligand_table[i];
      ligp=ligip->top;
      num=ligip->lig_index;
      if (num>0) {
        no_printf("ligand = %s state = %d\n",ligip->sym->name,state);
  
        lig_count=0;
        while (ligp!=NULL) {
  
          rp=ligp->effector;
          if (rp!=NULL) {
            r_index=ligp->index;
            state=rp->tiles[r_index]->viz_state;
          }
	  else {
	    state=ligip->viz_state;
	  }
          no_printf("ligand = %s state = %d\n",ligip->sym->name,state);
  
	  if (state==lig_state && state!=EXCLUDE_OBJ) {
	    v1=length_unit*ligp->pos.x;
	    v2=length_unit*ligp->pos.y;
	    v3=length_unit*ligp->pos.z;
            if (v1>=umin && v1<=umax) {
              if (v3>=vmin && v3<=vmax) {
                if (v2>=wmin && v2<=wmax) {
	          iv1=(int)((v1-umin)/voxel_uv);
	          iv3=(image_y-1)-(int)((v3-vmin)/voxel_uv);
                  lig_count++;
                  image_buffer[iv3*image_x+iv1]++;
                } 
              } 
            } 
	  }
  
	  ligp=ligp->next_ligand;
        }
      }
      state=ligip->viz_state;
      if (state==lig_state && state!=EXCLUDE_OBJ) {
      
        no_printf("ligand = %s state = %d\n",ligip->sym->name,state);
        no_printf("MCell: writing file voxel_image.ligand_%s.state_%d.%d.pgm\n",
          rx_name,lig_state,viz_iteration);
        sprintf(file_name,"voxel_image.ligand_%s.state_%d.%d.pgm",
          rx_name,lig_state,viz_iteration);
        if ((lig_header=fopen(file_name,"w"))==NULL) {
          fprintf(log_file,"MCell: error cannot open voxel image file %s\n",
            file_name);
          return(1);
        }
  
        fprintf(lig_header,"P2\n%d %d\n255\n",image_x,image_y);
        for (j=0;j<num_pixels;j++) {
          i_val=image_buffer[j];
          fprintf(lig_header,"%d\n",i_val);
        }
    
        fclose(lig_header);
      }
    }
    lig_state_p=lig_state_p->next;
  }

  free(image_buffer);

  no_printf("voxel image output done.\n");
  return(0);
}


int output_voxel_volume(struct frame_data_list *fdlp)
{
FILE *wall_header, *wall_verts;
FILE *eff_header, *eff_states, *eff_pos;
FILE *lig_header, *lig_states, *lig_pos;
struct viz_obj *vizp;
struct cmprt_data_list *cdlp;
struct cmprt_data *cdp;
struct wall_list *wlp;
struct wall *wp;
struct effector *rp;
struct rx *rxp;
struct ligand_info *ligip;
struct ligand *ligp;
struct vector3 glyph_axis,rotation_axis;
struct vector3 ab,bc,step_u,step_v,diagonal;
struct vector3 p0,p1,p2;
struct state_list *surf_state_p,*eff_state_p,*lig_state_p;
double rad,length,rotation_angle,final_angle,v_val;
double umin,umax,vmin,vmax,wmin,wmax,voxel_uvw;
double v1,v2,v3;
double u,v,w;
int iv1,iv2,iv3;
int uu,vv,vv_max,grid_size;
unsigned int i,j,k,ii,jj,kk;
int num;
int viz_iteration;
int viz_type;
int wall_index,eff_index,r_index,orient;
int pos_count,state_count,lig_count,eff_count;
int state,surf_state,eff_state,lig_state;
unsigned int index,n_eff;
unsigned int vol_x,vol_y,vol_z,vol_xy,num_voxels,*voxel_buffer;
byte viz_eff,viz_lig,viz_surf;
byte viz_surf_or_eff;
byte viz_eff_pos,viz_eff_states;
byte viz_lig_pos,viz_lig_states;
byte viz_surf_pos,viz_surf_states;
byte grid_shape;
byte i_val;
char *rx_name;
char file_name[1024];

  no_printf("voxel volume output...\n");

  viz_iteration=fdlp->viz_iteration;

  viz_type=fdlp->type;
  viz_eff=((viz_type==ALL_FRAME_DATA) || (viz_type==EFF_POS)
             || (viz_type==EFF_STATES));
  viz_lig=((viz_type==ALL_FRAME_DATA) || (viz_type==MOL_POS)
             || (viz_type==MOL_STATES));
  viz_surf=((viz_type==ALL_FRAME_DATA) || (viz_type==SURF_POS)
             || (viz_type==SURF_STATES));
  viz_surf_or_eff=(viz_surf || viz_eff);


  viz_eff_pos=((viz_type==ALL_FRAME_DATA) || (viz_type==EFF_POS));
  viz_eff_states=((viz_type==ALL_FRAME_DATA) || (viz_type==EFF_STATES));
  viz_lig_pos=((viz_type==ALL_FRAME_DATA) || (viz_type==MOL_POS));
  viz_lig_states=((viz_type==ALL_FRAME_DATA) || (viz_type==MOL_STATES));
  viz_surf_pos=((viz_type==ALL_FRAME_DATA) || (viz_type==SURF_POS));
  viz_surf_states=((viz_type==ALL_FRAME_DATA) || (viz_type==SURF_STATES));

  voxel_uvw=0.025;
  umin=-0.25;
  umax=0.25;
  vmin=-0.25;
  vmax=0.25;
  wmin=-0.5;
  wmax=0.0;

  vol_x=(unsigned int)((umax-umin)/voxel_uvw);
  vol_y=(unsigned int)((vmax-vmin)/voxel_uvw);
  vol_z=(unsigned int)((wmax-wmin)/voxel_uvw);

  vol_xy=vol_x*vol_y;
  num_voxels=vol_x*vol_y*vol_z;

  if ((voxel_buffer=(unsigned int *)malloc
                    (num_voxels*sizeof(unsigned int)))==NULL) {
    printf ("MCell: cannot store voxel image\n");
    return(1);
  }

/* voxelize effectors: */
eff_state_p=eff_state_head;
while (eff_state_p!=NULL) { 
  eff_state=eff_state_p->state;
  eff_count=0;

  if (viz_eff) {
    no_printf("Voxelizing effector state %d\n",eff_state);
    rx_name=eff_state_p->name;
    for (j=0;j<num_voxels;j++) {
      voxel_buffer[j]=0;
    }
    vizp = viz_obj_head;
    while(vizp!=NULL) {


      cdlp = vizp->cmprt_data_list;
      while(cdlp!=NULL) {
        cdp = cdlp->cmprt_data;
        no_printf("Traversing surfaces in compartment %s\n",cdp->sym->name);
        fflush(log_file);
        wlp = cdp->wall_list;
        while (wlp!=NULL) {
          wp = wlp->wall;
 
          /* dump the effectors */
          rp = wp->effectors;
          if (rp!=NULL) {
            n_eff=0;
            for (index=0;index<rp->n_tiles;index++) {
	      rxp=rp->tiles[index];
	      if (rxp!=NULL) {
	        state=rxp->viz_state;
	      }
	      else {
	        state=EXCLUDE_OBJ;
	      }
              if (state==eff_state) {
                rx_name=rxp->sym->name;
                n_eff++;
              }
            }

            no_printf("Voxelizing %d effectors...\n",n_eff);
            if (n_eff>0) {
              fflush(log_file);
              grid_shape=rp->grid_shape;
              grid_size=rp->grid_size;
              vectorize(wp->vert[0],wp->vert[1],&ab);
              vectorize(wp->vert[1],wp->vert[2],&bc);
              step_u.x=ab.x/grid_size;
              step_u.y=ab.y/grid_size;
              step_u.z=ab.z/grid_size;
              step_v.x=bc.x/grid_size;
              step_v.y=bc.y/grid_size;
              step_v.z=bc.z/grid_size;
              vectorize(&step_u,&step_v,&diagonal);
  
              index=0;
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
	          rxp=rp->tiles[index];
	          if (rxp!=NULL) {
	            state=rxp->viz_state;
	          }
	          else {
	            state=EXCLUDE_OBJ;
	          }
      
                  if (state==eff_state && state!=EXCLUDE_OBJ) {
	            v1=length_unit*(wp->vert[0]->x+p2.x);
	            v2=length_unit*(wp->vert[0]->y+p2.y);
	            v3=length_unit*(wp->vert[0]->z+p2.z);
                    if (v1>=umin && v1<=umax) {
                      if (v2>=vmin && v2<=vmax) {
                        if (v3>=wmin && v3<=wmax) {
	                  iv1=(int)((v1-umin)/voxel_uvw);
	                  iv2=(int)((v2-vmin)/voxel_uvw);
	                  iv3=(int)((v3-wmin)/voxel_uvw);
                          eff_count++;
                          voxel_buffer[iv3*vol_xy+iv2*vol_x+iv1]++;
                        } 
                      } 
                    } 
	          }
                  index++;
                }
              }
            }
          }
          wlp = wlp->next;
        }

        cdlp = cdlp->next;
      }

      vizp = vizp->next;
    }
  }
  no_printf("MCell: writing file voxel_volume.effector_%s.state_%d.%d.vox\n",
    rx_name,eff_state,viz_iteration);
  sprintf(file_name,"voxel_volume.effector_%s.state_%d.%d.vox",
    rx_name,eff_state,viz_iteration);
  if ((eff_header=fopen(file_name,"w"))==NULL) {
    fprintf(log_file,"MCell: error cannot open voxel volume file %s\n",
    file_name);
    return(1);
  }

  for (i=0;i<vol_x;i++) {
    u=i*voxel_uvw+umin;
    for (j=0;j<vol_y;j++) {
      v=j*voxel_uvw+vmin;
      jj=j*vol_x;
      for (k=0;k<vol_z;k++) {
        w=k*voxel_uvw+wmin;
        kk=k*vol_xy;
        ii=kk+jj+i;
        fprintf(eff_header,"%.9g %.9g %.9g %u\n",u,v,w,voxel_buffer[ii]);
      }
    }
  }
  fclose(eff_header);
  
  eff_state_p = eff_state_p->next;
}


/* voxelize ligands: */

  lig_state_p=lig_state_head;
  while (lig_state_p!=NULL) { 
    lig_state=lig_state_p->state;
    rx_name=lig_state_p->name;
    for (i=1;i<1+n_ligand_types;i++) {
      for (j=0;j<num_voxels;j++) {
        voxel_buffer[j]=0;
      }
      ligip=ligand_table[i];
      ligp=ligip->top;
      num=ligip->lig_index;
      if (num>0) {
        no_printf("ligand = %s state = %d\n",ligip->sym->name,state);
  
        lig_count=0;
        while (ligp!=NULL) {
  
          rp=ligp->effector;
          if (rp!=NULL) {
            r_index=ligp->index;
            state=rp->tiles[r_index]->viz_state;
          }
	  else {
	    state=ligip->viz_state;
	  }
          no_printf("ligand = %s state = %d\n",ligip->sym->name,state);
  
	  if (state==lig_state && state!=EXCLUDE_OBJ) {
	    v1=length_unit*ligp->pos.x;
	    v2=length_unit*ligp->pos.y;
	    v3=length_unit*ligp->pos.z;
            if (v1>=umin && v1<=umax) {
              if (v2>=vmin && v2<=vmax) {
                if (v3>=wmin && v3<=wmax) {
	          iv1=(int)((v1-umin)/voxel_uvw);
	          iv2=(int)((v2-vmin)/voxel_uvw);
	          iv3=(int)((v3-wmin)/voxel_uvw);
                  lig_count++;
                  voxel_buffer[iv3*vol_xy+iv2*vol_x+iv1]++;
                } 
              } 
            } 
	  }
  
	  ligp=ligp->next_ligand;
        }
      }
      state=ligip->viz_state;
      if (state==lig_state && state!=EXCLUDE_OBJ) {
      
        no_printf("ligand = %s state = %d\n",ligip->sym->name,state);
        no_printf("MCell: writing file voxel_volume.ligand_%s.state_%d.%d.vox\n",
          rx_name,lig_state,viz_iteration);
        sprintf(file_name,"voxel_volume.ligand_%s.state_%d.%d.vox",
          rx_name,lig_state,viz_iteration);
        if ((lig_header=fopen(file_name,"w"))==NULL) {
          fprintf(log_file,"MCell: error cannot open voxel volume file %s\n",
            file_name);
          return(1);
        }
  
        for (i=0;i<vol_x;i++) {
          u=i*voxel_uvw+umin;
          for (j=0;j<vol_y;j++) {
            v=j*voxel_uvw+vmin;
            jj=j*vol_x;
            for (k=0;k<vol_z;k++) {
              w=k*voxel_uvw+wmin;
              kk=k*vol_xy;
              ii=kk+jj+i;
              fprintf(lig_header,"%.9g %.9g %.9g %u\n",u,v,w,voxel_buffer[ii]);
            }
          }
        }

        fclose(lig_header);
      }
    }
    lig_state_p=lig_state_p->next;
  }

  free(voxel_buffer);

  no_printf("voxel volume output done.\n");
  return(0);
}





#endif /* if 0 */

/* **************************************************************** */


