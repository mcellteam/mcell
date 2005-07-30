
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
/* serial # of the frame_iteration in the linked list.
 each frame may have several iterations */
static int frames_iterations_count = 0;
/* total # of frames_iterations in the linked list*/	
static int total_frames_iterations = 0; 
static int obj_to_show_number; /* number of viz_obj objects in the world */
static int eff_to_show_number; /* number of types of effectors */ 
static int mol_to_show_number; /* number of types of 3D mol's */
/* arrays used in the function 'output_dreamm_objects_some_frame_data()' */
u_int *surf_states = NULL;
u_int *surf_pos = NULL;
u_int *surf_con = NULL;
char **obj_names = NULL;
u_int *eff_pos = NULL;
u_int *eff_orient = NULL;
u_int *eff_states = NULL;
char **eff_names = NULL;
u_int *mol_pos = NULL;
u_int *mol_orient = NULL;
u_int *mol_states = NULL;
char **mol_names = NULL;

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
          if(fdlp->type == ALL_FRAME_DATA) output_dreamm_objects_all_frame_data(fdlp);
          else
	  {
            frames_iterations_count++;
            output_dreamm_objects_some_frame_data(fdlp);
	  }
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
  int done, ii;
  struct species **species_list;

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
	nelp=nelp->next;
      }
      break;
    case OUTPUT_BY_TIME_LIST:
      while (nelp!=NULL) {
	fdlp->n_viz_iterations++;
	if (!done) {
	  if (nelp->value>=world->current_start_time) {
	    fdlp->viz_iterationll=(long long)(nelp->value/world->time_unit+ROUND_UP);
	    fdlp->curr_viz_iteration=nelp;
	    done=1;
	  }
      }
	nelp=nelp->next;
      }
      break;
     default:
         fprintf(stderr,"MCell: error - wrong frame_data_list list_type %d\n", fdlp->list_type);
         break;
    }
    fdlp=fdlp->next;
  }

  struct frame_data_list *fdl_ptr = world->frame_data_head;
  /* Count total number of frames miltiplied by the number
     of iterations per frame in the linked list of frames. */  
  while(fdl_ptr != NULL){
	total_frames_iterations += (fdl_ptr->n_viz_iterations);
        fdl_ptr = fdl_ptr->next;
  }

  /* find out the number of objects to visualize. */
  struct viz_obj *vizp = world->viz_obj_head;
  while(vizp != NULL){
        struct viz_child *vcp = vizp->viz_child_head;
        while(vcp != NULL){
		obj_to_show_number++;
                vcp = vcp->next;
        }
	vizp = vizp->next;
  }
  /* find out number of effectors and 3D molecules to visualize */
  u_int n_species = world->n_species;
  species_list=world->species_list;

  for (ii=0;ii<n_species;ii++) {

          if(species_list[ii]->viz_state == EXCLUDE_OBJ) continue;
          if((species_list[ii]->flags & ON_GRID) == ON_GRID) 
          { 
                /* calculate number of effectors */
	        eff_to_show_number++;
          }else if((species_list[ii]->flags & NOT_FREE) == 0) {
                /* calculate number of 3D molecules */
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
  int n_viz_iterations;
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
        fprintf(log_file,"MCell: error cannot open mesh elements file %s\n",
               file_name);
        return(1);
      }
    }

    if (viz_surf_states) {
      sprintf(file_name,"%s.mesh_element_states.%lld.dx",
               vizp->name,viz_iterationll);
      if ((wall_states_header=fopen(file_name,"wb"))==NULL) {
        fprintf(log_file,"MCell: error cannot open mesh element states file %s\n",
               file_name);
        return(1);
      }
    }

    if (viz_eff_pos) {
      sprintf(file_name,"%s.effector_site_positions.%lld.dx",vizp->name,viz_iterationll);
      if ((eff_pos_header=fopen(file_name,"wb"))==NULL) {
        fprintf(log_file,"MCell: error cannot open effector position file %s\n",file_name);
        return(1);
      }
    }
  
    if (viz_eff_states) {
      sprintf(file_name,"%s.effector_site_states.%lld.dx",vizp->name,viz_iterationll);
      if ((eff_states_header=fopen(file_name,"wb"))==NULL) {
        fprintf(log_file,"MCell: error cannot open effector states file %s\n",file_name);
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
                   v1=w->normal.x;
                   v2=w->normal.y;
                   v3=w->normal.z;
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
        fprintf(log_file,"MCell: error cannot open molecule positions header file %s\n",file_name);
        return(1);
      }
    }
    if (viz_mol_states) {
      sprintf(file_name,"%s.molecule_states.%lld.dx",world->molecule_prefix_name,viz_iterationll);
      if ((mol_states_header=fopen(file_name,"wb"))==NULL) {
        fprintf(log_file,"MCell: error cannot open molecule states header file %s\n",file_name);
        return(1);
      }
    }

    species_list=world->species_list;
    n_species=world->n_species;

    if ((viz_molp=(struct molecule ***)malloc(n_species*sizeof(struct molecule **)))==NULL) {
      return(1);
    }
    if ((viz_mol_count=(u_int *)malloc(n_species*sizeof(u_int)))==NULL) {
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
output_dreamm_objects_some_frame_data:
	In: struct frame_data_list *fdlp
	Out: 0 on success, 1 on error; output visualization files (*.dx)
             in dreamm group format are written.
        NB! NOT TESTED FOR SURFACE_MOLECULES!!!
**************************************************************************/

int output_dreamm_objects_some_frame_data(struct frame_data_list *fdlp)
{
  FILE *log_file;
  FILE *master_header = NULL;
  FILE *mesh_pos_data = NULL;  /* data file for wall vertices */
  FILE *mesh_states_data = NULL; /* data file for wall states */
  FILE *mol_pos_data = NULL;	/* data file for molecule positions */
  FILE *mol_states_data = NULL; /* data file for molecule states */
  FILE *mol_orient_data = NULL; /* data file for molecule orientations */
  struct viz_obj *vizp = NULL;
  struct viz_child *vcp = NULL;
  struct surface_grid *sg;
  struct wall *w,**wp;
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
  struct grid_molecule *gmol;
  struct molecule *molp,***viz_molp = NULL;       /* for 3D molecules */
  float v1,v2,v3;
  u_int spec_id, *viz_mol_count = NULL, *viz_grid_mol_count = NULL;
  static u_int main_index = 1;
  u_int group_index = 0;
  static u_int series_index = 0;
  int ii,jj;
  int vi1,vi2,vi3;
  /* int vi4;*/
  int num;
  long long viz_iterationll;
  int n_viz_iterations;
  int element_data_count;
  int state;
  int viz_type;
  unsigned int index; 
  /* indices used in arrays "u_int *surf_pos", etc. */
  int surf_pos_index = 0; 
  int surf_con_index = 0;
  int surf_states_index = 0;
  int eff_states_index = 0;
  int eff_pos_index = 0;
  int eff_orient_index = 0;
  int mol_states_index = 0;
  int mol_pos_index = 0;
  int mol_orient_index = 0;
  int word;
  byte *word_p;
  byte viz_eff_pos_flag,viz_eff_states_flag;	/* flags */
  byte viz_mol_pos_flag,viz_mol_states_flag;	/* flags */
  byte viz_mol_pos_states_flag;	/* flag */
  byte viz_surf_pos_flag,viz_surf_states_flag;	/* flags */
  char file_name[1024];
  char *ch_ptr = NULL; /* pointer used to extract data file name */
  char *grid_mol_name = NULL; /* points to the name of the grid molecule */
  char mesh_pos_name[1024]; /* wall vertices data file name */
  char mesh_states_name[1024]; /* wall states data file name */
  char mol_pos_name[1024]; /* molecule_positions data file name */
  char mol_states_name[1024]; /* molecule_positions data file name */
  char mol_orient_name[1024]; /* molecule orientations data file name */
  char buffer[100]; /* used to write 'frame_data' object information */
  /* linked list that stores data for the 'frame_data' object */
  static struct infinite_string_array frame_data_list;
  static int frame_data_count = 0; /* count elements in frame_data_list array.*/
  /* linked lists that stores data for the 'frame_numbers' object */
  static struct infinite_int_array frame_numbers_pos;
  static struct infinite_int_array frame_numbers_mesh;
  static struct infinite_int_array frame_numbers_mol;
  static struct infinite_int_array frame_numbers_eff;
  static int frame_numbers_count = 0; /* count elements in frame_numbers_pos array. */
  char my_byte_order[8];  /* shows binary ordering ('lsb' or 'msb') */
  static int mesh_pos_byte_offset = 0;  /* defines position of the object data
                                  in the mesh positions binary data file */
  static int mesh_states_byte_offset = 0; /* defines position of the object data
                                    in the mesh states binary file */
  static int mol_pos_byte_offset = 0; /*defines position of the object data 
                           in the molecule positions binary file */
  static int mol_orient_byte_offset = 0; /*defines position of the object data 
                           in the molecule orientations binary file */
  static int mol_states_byte_offset = 0; /* defines position of the object data
                              in the molecule states binary file. */
  int mol_pos_byte_offset_prev = 0; /*defines position of the object data 
                           in the effector positions binary file */
  int mol_orient_byte_offset_prev = 0; /*defines position of the object da                            ta in the effector orientations binary file */
  int mol_states_byte_offset_prev = 0; /* defines position of the object                                           data in the effector states binary file. */
  
  /* linked list that stores members of the group for surfaces */
  struct num_expr_list *surf_head = NULL;  
  struct num_expr_list *surf_end = NULL;
  /* linked list that stores members of the group for 3D molecules */
  struct num_expr_list *mol_head = NULL;  
  struct num_expr_list *mol_end = NULL;
  /* linked list that stores members of the group for effectors */
  struct num_expr_list *eff_head = NULL;  
  struct num_expr_list *eff_end = NULL;

  struct num_expr_list *curr_ptr = NULL;
  struct num_expr_list *newNode = NULL;
  struct num_expr_list *temp = NULL;
  /* placeholders for the members of the combined group */
  static int mesh_group_index = 0;
  static int mol_group_index = 0;
  static int eff_group_index = 0;
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

  /* points to the values of the current iteration steps
     for certain frame types. */
  static long long curr_surf_pos_iteration_step = -1;
  static long long curr_surf_states_iteration_step = -1;
  static long long curr_eff_pos_iteration_step = -1;
  static long long curr_eff_states_iteration_step = -1;
  static long long curr_mol_pos_states_iteration_step = -1;

  /* points to the special case in iteration steps
     when both values for SURF_POS and SURF_STATES
     are equal.  E.g. SURF_POS = [0,200], SURF_STATES = [0,100, 200,300]
     special_surf_iteration_step = [0,200]
     Same for (EFF_POS,EFF_STATES) 
  */
  static long long special_surf_iteration_step = -1;
  static long long special_eff_iteration_step = -1;
  /* counts number of this function executions
     for special_surf_iteration_step cases. */
  static int special_surf_frames_counter = 0;
  static int special_eff_frames_counter = 0;
  

  struct frame_data_list * fdl_ptr; 
  u_int n_species = world->n_species;
  species_list=world->species_list;

  log_file=world->log_file;
  no_printf("Viz output in DREAMM_V3 mode...\n");

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

  viz_type=fdlp->type;

  if(viz_type == ALL_FRAME_DATA) {
   	fprintf(log_file, "Wrong viz_type. Received %d type.  Expected other type.\n", viz_type);
   	return 1;
  }
  
  /* initialize flags */
  viz_eff_pos_flag = (viz_type==EFF_POS);
  viz_eff_states_flag = (viz_type==EFF_STATES);
  viz_mol_pos_flag = (viz_type==MOL_POS);
  viz_mol_states_flag = (viz_type==MOL_STATES);
  viz_mol_pos_states_flag = (viz_type==MOL_POS_STATES);
  viz_surf_pos_flag = (viz_type==SURF_POS);
  viz_surf_states_flag = (viz_type==SURF_STATES);

  if(viz_mol_pos_flag || viz_mol_states_flag) {
   	fprintf(log_file, "Wrong vizualization option for molecules. Expected type such as MOLECULE_POSITIONS_STATES.\n");
   	return 1;
  }

  /* these arrays are used to hold values of main_index for objects */
  if(viz_surf_pos_flag || viz_surf_states_flag)
  {
     if(obj_to_show_number > 0)
     {
  	if(surf_states == NULL){ 
     		if ((surf_states=(u_int *)malloc(obj_to_show_number*sizeof(u_int)))==NULL)      {
         		fprintf(stderr,"MCell: memory allocation error\n");
         		return (1);
     		}
        	for(ii = 0; ii < obj_to_show_number; ii++){
			surf_states[ii] = 0;
     		}
   	}
   	if(surf_pos == NULL){
     		if ((surf_pos=(u_int *)malloc(obj_to_show_number*sizeof(u_int)))==NULL) {
         		fprintf(stderr,"MCell: memory allocation error\n");
         		return (1);
     		}
     		for(ii = 0; ii < obj_to_show_number; ii++){
			surf_pos[ii] = 0;
     		}
   	}
   	if(surf_con == NULL){
      		if ((surf_con=(u_int *)malloc(obj_to_show_number*sizeof(u_int)))==NULL) {
         		fprintf(stderr,"MCell: memory allocation error\n");
         		return (1);
      		}
      		for(ii = 0; ii < obj_to_show_number; ii++){
			surf_con[ii] = 0;
      		}
   	}
   	/* initialize array of viz_objects names */
   	if(obj_names == NULL){
      		if ((obj_names = (char **)malloc(obj_to_show_number*sizeof(char *)))==NULL) {
         		fprintf(stderr,"MCell: memory allocation error\n");
         		return (1);
      		}
      		for(ii = 0; ii < obj_to_show_number; ii++){
			obj_names[ii] = NULL;
      		}
   	}
  
     	ii = 0;
     	vizp = world->viz_obj_head;
     	while(vizp != NULL){
		vcp = vizp->viz_child_head;
        	while(vcp != NULL){
         		if(obj_names != NULL){
         		   obj_names[ii] = my_strdup(vcp->obj->sym->name);
                           if(obj_names[ii] == NULL){
         		       fprintf(stderr,"MCell: memory allocation error\n");
         		       return (1);
         		   }
                           ii++;
         		}
         		vcp = vcp->next;
       		}
       		vizp = vizp->next;
    	}
       } /* end if (obj_to_show_number > 0) */
    }  /* end if (viz_surf_pos_flag || viz_surf_states_flag) */

     /* these arrays are used to hold values of main_index for effectors 
      and 3D molecules */
   if(viz_eff_pos_flag || viz_eff_states_flag)
   {
     if(eff_to_show_number > 0)
     {
     	if(eff_states == NULL){ 
     		if((eff_states=(u_int *)malloc(eff_to_show_number*sizeof(u_int)))==NULL)        {
         		fprintf(stderr,"MCell: memory allocation error\n");
         		return (1);
     		}
     		for(ii = 0; ii < eff_to_show_number; ii++){
			eff_states[ii] = 0;
     		}
     	}
     	if(eff_pos == NULL){
     		if ((eff_pos=(u_int *)malloc(eff_to_show_number*sizeof(u_int)))==NULL) 	{
         		fprintf(stderr,"MCell: memory allocation error\n");
         		return (1);
     		}
     		for(ii = 0; ii < eff_to_show_number; ii++){
			eff_pos[ii] = 0;
     		}
     	}
     	if(eff_orient == NULL){
     		if ((eff_orient = (u_int *)malloc(eff_to_show_number*sizeof(u_int)))==NULL) 	{
         		fprintf(stderr,"MCell: memory allocation error\n");
         		return (1);
     		}
     		for(ii = 0; ii < eff_to_show_number; ii++){
			eff_orient[ii] = 0;
     		}
     	}
   	/* initialize array of grid_mol's names */
   	if(eff_names == NULL){
      		if ((eff_names = (char **)malloc(eff_to_show_number*sizeof(char *)))==NULL) {
         		fprintf(stderr,"MCell: memory allocation error\n");
         		return (1);
      		}
      		for(ii = 0; ii < eff_to_show_number; ii++){
			eff_names[ii] = NULL;
      		}
   	}
        index = 0;
        for (ii=0;ii<n_species;ii++) {

          if(species_list[ii]->viz_state == EXCLUDE_OBJ) continue;
          if((species_list[ii]->flags & ON_GRID) == ON_GRID) 
          { 
             eff_names[index] = my_strdup(species_list[ii]->sym->name);
             if(eff_names[index] == NULL){
                 fprintf(stderr,"MCell: memory allocation error\n");
                 return (1);
             }
             index++;
          }
        }
        index = 0;
     } /* end if (eff_to_show_number > 0) */
   }  /* end if(viz_eff_pos_flag || viz_eff_states_flag) */
   
   if(viz_mol_pos_states_flag )
   {
     if(mol_to_show_number > 0)
     {
     	if(mol_states == NULL){ 
     		if((mol_states=(u_int *)malloc(mol_to_show_number*sizeof(u_int)))==NULL)        {
         		fprintf(stderr,"MCell: memory allocation error\n");
         		return (1);
     		}
     		for(ii = 0; ii < mol_to_show_number; ii++){
			mol_states[ii] = 0;
     		}
     	}
     	if(mol_pos == NULL){
     		if ((mol_pos=(u_int *)malloc(mol_to_show_number*sizeof(u_int)))==NULL) 	{
         		fprintf(stderr,"MCell: memory allocation error\n");
         		return (1);
     		}
     		for(ii = 0; ii < mol_to_show_number; ii++){
			mol_pos[ii] = 0;
     		}
     	}
     	if(mol_orient == NULL){
     		if ((mol_orient = (u_int *)malloc(mol_to_show_number*sizeof(u_int)))==NULL) 	{
         		fprintf(stderr,"MCell: memory allocation error\n");
         		return (1);
     		}
     		for(ii = 0; ii < mol_to_show_number; ii++){
			mol_orient[ii] = 0;
     		}
     	}
   	/* initialize array of mol's names */
   	if(mol_names == NULL){
      		if ((mol_names = (char **)malloc(mol_to_show_number*sizeof(char *)))==NULL) {
         		fprintf(stderr,"MCell: memory allocation error\n");
         		return (1);
      		}
      		for(ii = 0; ii < mol_to_show_number; ii++){
			mol_names[ii] = NULL;
      		}
   	}
        index = 0;
        for (ii=0;ii<n_species;ii++) {

          if(species_list[ii]->viz_state == EXCLUDE_OBJ) continue;
          
          if((species_list[ii]->flags & NOT_FREE) == 0) 
          { 
             mol_names[index] = my_strdup(species_list[ii]->sym->name);
             if(mol_names[index] == NULL){
                 fprintf(stderr,"MCell: memory allocation error\n");
                 return (1);
             }
             index++;
          }
        }
        index = 0;

     } /* end if (mol_to_show_number > 0) */
  }  /* end if(viz_mol_pos_flag || viz_mol_states_flag) */

  /* Open master header file. */
  sprintf(file_name,"%s.master_header.dx",world->molecule_prefix_name);

  if(count_master_header == 0){
      if ((master_header=fopen(file_name,"w"))==NULL) {
           fprintf(log_file,"MCell: error cannot open master header file %s\n",file_name);
           return(1);
      }else{}
      count_master_header++;

  }else{
      if ((master_header=fopen(file_name,"a"))==NULL) {
        fprintf(log_file,"MCell: error cannot open master header file %s\n",file_name);
        return(1);
      }else{}
      count_master_header++;
  }


  if (viz_mol_pos_states_flag || viz_eff_pos_flag) {
      sprintf(file_name,"%s.molecule_positions.bin",world->molecule_prefix_name);
     /* remove the folder name from the molecule_positions data file name */
     ch_ptr = strrchr(file_name, '/');
     ++ch_ptr;
     strcpy(mol_pos_name, ch_ptr);
     
     if (count_mol_pos_data == 0){
        if ((mol_pos_data=fopen(file_name,"wb"))==NULL) {
           fprintf(log_file,"MCell: error cannot open molecule positions data file %s\n",file_name);
           return(1);
        }else{}
        count_mol_pos_data++;
     }else{
        if ((mol_pos_data=fopen(file_name,"ab"))==NULL) {
           fprintf(log_file,"MCell: error cannot open molecule positions data file %s\n",file_name);
           return(1);
        }else{}
        count_mol_pos_data++;

    }

      sprintf(file_name,"%s.molecule_orientations.bin",world->molecule_prefix_name);
     /* remove the folder name from the molecule_positions data file name */
     ch_ptr = strrchr(file_name, '/');
     ++ch_ptr;
     strcpy(mol_orient_name, ch_ptr);

     if (count_mol_orient_data == 0){
        if ((mol_orient_data=fopen(file_name,"wb"))==NULL) {
           fprintf(log_file,"MCell: error cannot open molecule orientations data file %s\n",file_name);
           return(1);
        }else{}
        count_mol_orient_data++;
     }else{
        if ((mol_orient_data=fopen(file_name,"ab"))==NULL) {
           fprintf(log_file,"MCell: error cannot open molecule orientations data file %s\n",file_name);
           return(1);
        }else{}
        count_mol_orient_data++;

     }

  }

    if (viz_mol_pos_states_flag || viz_eff_states_flag) {
       sprintf(file_name,"%s.molecule_states.bin",world->molecule_prefix_name);
       /* remove the folder name from the molecule_states data file name */
       ch_ptr = strrchr(file_name, '/');
       ++ch_ptr;
       strcpy(mol_states_name, ch_ptr);
     
       if (count_mol_states_data == 0){
        	if ((mol_states_data = fopen(file_name,"wb"))==NULL) {
           		fprintf(log_file,"MCell: error cannot open molecule states data file %s\n",file_name);
           		return(1);
         	}else{}
         	count_mol_states_data++;
       }else{
        	if ((mol_states_data = fopen(file_name,"ab"))==NULL) {
           		fprintf(log_file,"MCell: error cannot open molecule states data file %s\n",file_name);
           		return(1);
         	}else{}
         	count_mol_states_data++;
       }
    }
    
    /* find out the values of the current iteration steps for
       the (SURF_POS,SURF_STATES), (EFF_POS, EFF_STATES), 
       (MOL_POS_STATES) frames combinations */
    fdl_ptr = world->frame_data_head;
    while(fdl_ptr != NULL){
	if(fdl_ptr->type == SURF_STATES) 
        {
             curr_surf_states_iteration_step = fdl_ptr->viz_iterationll;
        }else if(fdl_ptr->type == SURF_POS){
             curr_surf_pos_iteration_step = fdl_ptr->viz_iterationll;
	}else if(fdl_ptr->type == EFF_STATES) 
        {
             curr_eff_states_iteration_step = fdl_ptr->viz_iterationll;
        }else if(fdl_ptr->type == EFF_POS){
             curr_eff_pos_iteration_step = fdl_ptr->viz_iterationll;
	}else if(fdl_ptr->type == MOL_POS_STATES) 
        {
             curr_mol_pos_states_iteration_step = fdl_ptr->viz_iterationll;
        }	
        fdl_ptr = fdl_ptr->next;
    }
    /* If the values of the current iteration steps for SURF_POS and
       SURF_STATES are equal set this value to the special_surf_iteration_step.
       Do the same for the (EFF_POS,EFF_STATES).
    */
    if(curr_surf_states_iteration_step == curr_surf_pos_iteration_step){
	special_surf_iteration_step = curr_surf_states_iteration_step;
    }
    if(curr_eff_states_iteration_step == curr_eff_pos_iteration_step){
	special_eff_iteration_step = curr_eff_states_iteration_step;
    }



/* dump walls */
  if (viz_surf_pos_flag || viz_surf_states_flag) {
     vizp = world->viz_obj_head;
     
     while(vizp!=NULL) {
   
      if (viz_surf_pos_flag) {
      	sprintf(file_name,"%s.mesh_positions.bin",vizp->name);
     
     	/* remove the folder name from the mesh_positions data file name */
     	ch_ptr = strrchr(file_name, '/');
     	++ch_ptr;
     	strcpy(mesh_pos_name, ch_ptr);

       if (count_mesh_pos_data == 0){
      	 if ((mesh_pos_data=fopen(file_name,"wb"))==NULL) {
           fprintf(log_file,"MCell: error cannot open mesh positions file %s\n",
               file_name);
           return(1);
         }else{}
         count_mesh_pos_data++;
       }else{
      	 if ((mesh_pos_data=fopen(file_name,"ab"))==NULL) {
           fprintf(log_file,"MCell: error cannot open mesh positions file %s\n",
               file_name);
           return(1);
         }else{}
         count_mesh_pos_data++;
       }

      }


      if (viz_surf_states_flag) {
         sprintf(file_name,"%s.mesh_states.bin", vizp->name);
     
        /* remove the folder name from the mesh_states data file name */
        ch_ptr = strrchr(file_name, '/');
        ++ch_ptr;
        strcpy(mesh_states_name, ch_ptr);

       if (count_mesh_states_data == 0){
        if ((mesh_states_data=fopen(file_name,"wb"))==NULL) {
           fprintf(log_file,"MCell: error cannot open mesh states file %s\n",
               file_name);
           return(1);
        }else{}
         count_mesh_states_data++;
       }else{
          if ((mesh_states_data=fopen(file_name,"ab"))==NULL) {
             fprintf(log_file,"MCell: error cannot open mesh states file %s\n",
               file_name);
             return(1);
          }else{}
          count_mesh_states_data++;
       }
      }

    /* Traverse all visualized compartments 
       output mesh element positions and connections */
    

    vcp = vizp->viz_child_head;

    /* check for the special_iteration_step  */  
    if(viz_surf_pos_flag || viz_surf_states_flag){	    
	    if(fdlp->viz_iterationll == special_surf_iteration_step){
		special_surf_frames_counter++;
            }
    }
    

    while(vcp!=NULL) {
      objp = vcp->obj;
      pop=(struct polygon_object *)objp->contents;
      if (objp->object_type==BOX_OBJ) {
#if 0	
        if (viz_surf_pos_flag || viz_surf_states_flag)
        {
          element_data_count=0.5*objp->n_walls_actual;
        
          if (viz_surf_pos_flag) {
             
	      fprintf(master_header, "object %d class array type float rank 1 shape 3 items %d %s binary data file %s,%d # %s.positions #\n",
              main_index,objp->n_verts,my_byte_order,mesh_pos_name, mesh_pos_byte_offset, objp->sym->name);
            fprintf(master_header,
              "\tattribute \"dep\" string \"positions\"\n\n");
            // output box vertices 
            for (ii=0;ii<objp->n_verts;ii++) {
              v1 = world->length_unit*objp->verts[ii].x;
              v2 = world->length_unit*objp->verts[ii].y;
              v3 = world->length_unit*objp->verts[ii].z;
              fwrite(&v1,sizeof v1,1,mesh_pos_data);
              fwrite(&v2,sizeof v2,1,mesh_pos_data);
              fwrite(&v3,sizeof v3,1,mesh_pos_data);
              mesh_pos_byte_offset += (sizeof(v1) + sizeof(v2) + sizeof(v3));
            }
            surf_pos[surf_pos_index] = main_index;
            surf_pos_index++;
            main_index++;
          
    
            /* output box wall connections */
            fprintf(master_header,
              "object %d class array type int rank 1 shape 4 items %d %s binary data file %s,%d # %s.connections #\n",
              main_index,element_data_count,my_byte_order, mesh_pos_name, mesh_pos_byte_offset, objp->sym->name);
            fprintf(master_header,
              "\tattribute \"ref\" string \"positions\"\n");
            fprintf(master_header,
              "\tattribute \"element type\" string \"quads\"\n\n");
             

            for (ii=0;ii<objp->n_walls;ii+=2) {
              if (!get_bit(pop->side_stat,ii)) {
                switch (ii) {
                  case TP:
                    vi1=3;
                    vi2=7;
                    vi3=1;
                    vi4=5;
                    fwrite(&vi1,sizeof vi1,1,mesh_pos_data);
                    fwrite(&vi2,sizeof vi2,1,mesh_pos_data);
                    fwrite(&vi3,sizeof vi3,1,mesh_pos_data);
                    fwrite(&vi4,sizeof vi4,1,mesh_pos_data);
                    mesh_pos_byte_offset += (sizeof(vi1) + sizeof(vi2) + sizeof(vi3) + sizeof(vi4));
                  break;
                  case BOT:
                    vi1=0;
                    vi2=4;
                    vi3=2;
                    vi4=6;
                    fwrite(&vi1,sizeof vi1,1,mesh_pos_data);
                    fwrite(&vi2,sizeof vi2,1,mesh_pos_data);
                    fwrite(&vi3,sizeof vi3,1,mesh_pos_data);
                    fwrite(&vi4,sizeof vi4,1,mesh_pos_data);
                    mesh_pos_byte_offset += (sizeof(vi1) + sizeof(vi2) + sizeof(vi3) + sizeof(vi4));
                  break;
                  case FRNT:
                    vi1=4;
                    vi2=0;
                    vi3=5;
                    vi4=1;
                    fwrite(&vi1,sizeof vi1,1,mesh_pos_data);
                    fwrite(&vi2,sizeof vi2,1,mesh_pos_data);
                    fwrite(&vi3,sizeof vi3,1,mesh_pos_data);
                    fwrite(&vi4,sizeof vi4,1,mesh_pos_data);
                    mesh_pos_byte_offset += (sizeof(vi1) + sizeof(vi2) + sizeof(vi3) + sizeof(vi4));
                  break;
                  case BCK:
                    vi1=2;
                    vi2=6;
                    vi3=3;
                    vi4=7;
                    fwrite(&vi1,sizeof vi1,1,mesh_pos_data);
                    fwrite(&vi2,sizeof vi2,1,mesh_pos_data);
                    fwrite(&vi3,sizeof vi3,1,mesh_pos_data);
                    fwrite(&vi4,sizeof vi4,1,mesh_pos_data);
                    mesh_pos_byte_offset += (sizeof(vi1) + sizeof(vi2) + sizeof(vi3) + sizeof(vi4));
                  break;
                  case LFT:
                    vi1=0;
                    vi2=2;
                    vi3=1;
                    vi4=3;
                    fwrite(&vi1,sizeof vi1,1,mesh_pos_data);
                    fwrite(&vi2,sizeof vi2,1,mesh_pos_data);
                    fwrite(&vi3,sizeof vi3,1,mesh_pos_data);
                    fwrite(&vi4,sizeof vi4,1,mesh_pos_data);
                    mesh_pos_byte_offset += (sizeof(vi1) + sizeof(vi2) + sizeof(vi3) + sizeof(vi4));
                  break;
                  case RT:
                    vi1=6;
                    vi2=4;
                    vi3=7;
                    vi4=5;
                    fwrite(&vi1,sizeof vi1,1,mesh_pos_data);
                    fwrite(&vi2,sizeof vi2,1,mesh_pos_data);
                    fwrite(&vi3,sizeof vi3,1,mesh_pos_data);
                    fwrite(&vi4,sizeof vi4,1,mesh_pos_data);
                    mesh_pos_byte_offset += (sizeof(vi1) + sizeof(vi2) + sizeof(vi3) + sizeof(vi4));
                  break;
                }
              }
            }
            surf_con[surf_con_index] = main_index;
            main_index++;
            surf_con_index++;


          } /* end viz_surf_pos_flag */

          if (viz_surf_states_flag) {

            fprintf(master_header,
              "object %d class array type int rank 0 items %d %s binary data file %s,%d # *%s.states #\n", main_index, element_data_count, my_byte_order, mesh_states_name, mesh_states_byte_offset, objp->sym->name);
            fprintf(master_header,"\tattribute \"dep\" string \"connections\"\n\n");
           surf_states[surf_states_index] = main_index;
           surf_states_index++;
           main_index++;
 
	    for (ii=0;ii<objp->n_walls;ii+=2) {
              if (!get_bit(pop->side_removed,ii)) {
                state=objp->viz_state[ii];
                fwrite(&state,sizeof (state),1,mesh_states_data);
                mesh_states_byte_offset += sizeof(state);
              }
            }

 
          } /* end (viz_surf_states_flag) */

        } /* end (viz_surf_pos_flag || viz_surf_states_flag) for BOX_OBJ */
#endif
      }  /* end BOX_OBJ */



      if (objp->object_type==POLY_OBJ || objp->object_type==BOX_OBJ) {
        if (viz_surf_pos_flag || viz_surf_states_flag)
        {
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
              v1 = world->length_unit*objp->verts[ii].x;
              v2 = world->length_unit*objp->verts[ii].y;
              v3 = world->length_unit*objp->verts[ii].z;
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
              "attribute \"element type\" string \"triangles\"\n\n");

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


        } /* end (viz_surf_pos_flag || viz_surf_states_flag) for POLY_OBJ */
      }	/* end POLY_OBJ */
      vcp = vcp->next;
      }
      
      vizp = vizp->next;
      }
    } /* end (viz_surf_pos_flag || viz_surf_states_flag) for vizp */

    
    /* build fields here */
      int show_meshes = 0;
      if((viz_surf_pos_flag || viz_surf_states_flag) && ((special_surf_frames_counter == 0) || (special_surf_frames_counter == 2))){
		show_meshes = 1;
      }
 
    if(show_meshes)
    {
         for(ii = 0; ii < obj_to_show_number; ii++)
         {
             fprintf(master_header,
                "object %d field   # %s #\n",main_index, obj_names[ii]);
             if(surf_pos[ii] > 0){
             	fprintf(master_header,
                 "\tcomponent \"positions\" value %d\n",surf_pos[ii]);
             }
             if(surf_con[ii] > 0){
             	fprintf(master_header,
                 "\tcomponent \"connections\" value %d\n",surf_con[ii]);
             }
             if(surf_states[ii] > 0){
             	fprintf(master_header,
                	"\tcomponent \"state_values\" value %d\n",surf_states[ii]);
	     }
             fprintf(master_header, "\n");
	     /* put number of this field into the linked list 
                of the field numbers for meshes */
             newNode = (struct num_expr_list *)malloc(sizeof (struct num_expr_list));
	     if(newNode == NULL) return 1;
             newNode->value = (double)main_index;
             newNode->next = NULL;    

             if((surf_head == NULL) && (surf_end == NULL)){
		surf_head = surf_end = newNode;
	     }else{
		surf_end->next = newNode;
		surf_end = newNode;
             }
             main_index++;
           }     
     }
    
    /* create a group object for all meshes. */
    if(show_meshes)
    {
        vizp = world->viz_obj_head;
        if(vizp != NULL){

        	fprintf(master_header,"object %d group # %s #\n",main_index, "meshes");
        	mesh_group_index = main_index;
        	main_index++;

        	curr_ptr = surf_head;
        	while(vizp != NULL) {
         		vcp = vizp->viz_child_head;
              
         		while(vcp != NULL){
                		if(curr_ptr != NULL){
                   	    		group_index = (int)curr_ptr->value;
                		}else{
		   	    		break;
                		}
             			fprintf(master_header,"\tmember \"%s\" value %d\n",vcp->obj->sym->name,group_index);
             	       		curr_ptr = curr_ptr->next;		
		       		vcp = vcp->next;
         		}
         		vizp = vizp->next;
        	}       
        	fprintf(master_header, "\n\n"); 
      	}  /* end (if vizp) */
   } 

      
    /* Visualize molecules. */

/* dump grid molecules. */
    /* check for the special_iteration_step  */  
    if(viz_eff_pos_flag || viz_eff_states_flag){	    
	    if(fdlp->viz_iterationll == special_eff_iteration_step){
		special_eff_frames_counter++;
            }
    }
    
    /* create references to the numbers of grid molecules of each name. */
    if ((viz_grid_mol_count=(u_int *)malloc(n_species*sizeof(u_int)))==NULL) {
      return(1);
    }

   for (ii = 0; ii < n_species; ii++)
   {
      /* perform initialization */
      spec_id=species_list[ii]->species_id;
      viz_grid_mol_count[spec_id]=0;

     if((species_list[ii]->flags & ON_GRID) == 0) continue; 
     if(species_list[ii]->viz_state == EXCLUDE_OBJ) continue; 
        
     grid_mol_name = species_list[ii]->sym->name;
    
     if(viz_eff_pos_flag || viz_eff_states_flag){

      	mol_pos_byte_offset_prev = mol_pos_byte_offset;
      	mol_orient_byte_offset_prev = mol_orient_byte_offset;
      	mol_states_byte_offset_prev = mol_states_byte_offset;
      
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
	                if (gmol!=NULL) {
	                   state=sg->mol[index]->properties->viz_state;
	                }else{
			   state = EXCLUDE_OBJ;
                        }
                        if (state!=EXCLUDE_OBJ) {
                          spec_id = gmol->properties->species_id;
                          viz_grid_mol_count[spec_id]++;
                          if(strcmp(gmol->properties->sym->name,
                            grid_mol_name) == 0){
                  		if (viz_eff_pos_flag) {
                                   /* write positions information */
	           		   v1=world->length_unit*p0.x;
	           		   v2=world->length_unit*p0.y;
	           		   v3=world->length_unit*p0.z;
	           		   fwrite(&v1,sizeof v1,1,mol_pos_data);
	               		   fwrite(&v2,sizeof v2,1,mol_pos_data);
	           		   fwrite(&v3,sizeof v3,1,mol_pos_data);
                                   mol_pos_byte_offset += (sizeof(v1) + sizeof(v2) + sizeof(v3));
                                   /* write orientations information */
                   		   v1=w->normal.x;
                   		   v2=w->normal.y;
                   		   v3=w->normal.z;
                   		   fwrite(&v1,sizeof v1,1,mol_orient_data);
                   		   fwrite(&v2,sizeof v2,1,mol_orient_data);
                   		   fwrite(&v3,sizeof v3,1,mol_orient_data);
                                   mol_orient_byte_offset += (sizeof(v1) + sizeof(v2) + sizeof(v3));

                  		}
                          }
                           
                        } /* end if */
                     } /* end for */
                   } /* end if */
                 } /* end if */
              } /* end for */
              vcp = vcp->next;
            } /* end while (vcp) */
                                              
            vizp = vizp->next;
         } /* end while (vzp) */

        num = viz_grid_mol_count[ii];
        
        if(viz_eff_pos_flag)
        {   
        	fprintf(master_header,"object %d class array type float rank 1 shape 3 items %d %s binary data file %s,%d # %s positions #\n",main_index,num,my_byte_order, mol_pos_name, mol_pos_byte_offset_prev, grid_mol_name);
        	fprintf(master_header,"\tattribute \"dep\" string \"positions\"\n\n");
                eff_pos[eff_pos_index] = main_index;
                eff_pos_index++;
        	main_index++;
   
        	fprintf(master_header,"object %d class array type float rank 1 shape 3 items %d %s binary data file %s,%d   # %s orientations #\n",main_index,num,my_byte_order, mol_orient_name, mol_orient_byte_offset_prev, species_list[ii]->sym->name);
        	fprintf(master_header,"\tattribute \"dep\" string \"positions\"\n\n");
                eff_orient[eff_orient_index] = main_index;
                eff_orient_index++;
        	main_index++;
         }

        if (viz_eff_states_flag) {
            /* write states information. */
            fwrite(&state,sizeof state,1,mol_states_data);
            mol_states_byte_offset += (sizeof state);
        
	    fprintf(master_header,"object %d class constantarray type int items %d %s binary data file %s,%d  # %s states #\n",main_index,num, my_byte_order,mol_states_name,mol_states_byte_offset_prev, species_list[ii]->sym->name);

            fprintf(master_header,"\tattribute \"dep\" string \"positions\"\n\n");
	     eff_states[eff_states_index] = main_index;
             eff_states_index++;
             main_index++;
        }
  
      } /* end if(viz_eff_pos_flag || viz_eff_states_flag) */

   } /* end for */


/* build fields for grid molecules here */
      int show_effectors = 0;

      if((viz_eff_pos_flag || viz_eff_states_flag) && ((special_eff_frames_counter ==0) || (special_eff_frames_counter == 2))){
		show_effectors = 1;
      }

  if(show_effectors){

       for(ii = 0; ii < eff_to_show_number; ii++)
       {

             fprintf(master_header,
                "object %d field   # %s #\n",main_index, eff_names[ii]);
             if(eff_pos[ii] > 0){
             	fprintf(master_header,
                 "\tcomponent \"positions\" value %d\n",eff_pos[ii]);
             }
             if(eff_orient[ii] > 0){     
             	fprintf(master_header,
                 	"\tcomponent \"data\" value %d # orientations #\n",eff_orient[ii]);
             }
             if(eff_states[ii] > 0){
             	fprintf(master_header,
                	"\tcomponent \"state_values\" value %d\n",eff_states[ii]);
             }
             fprintf(master_header, "\n");
             /* put number of this field into the linked list 
                of the field numbers for grid molecules */
             newNode = (struct num_expr_list *)malloc(sizeof (struct num_expr_list));
	     if(newNode == NULL) {
		fprintf(log_file, "Memory allocation error\n");
		return 1;
             }
             newNode->value = (double)main_index;
             newNode->next = NULL;    

             if((eff_head == NULL) && (eff_end == NULL)){
		eff_head = eff_end = newNode;
	     }else{
		eff_end->next = newNode;
		eff_end = newNode;
             }
             main_index++;
      }
  }



/* dump 3D and surface molecules: */

  if(viz_mol_pos_states_flag){	 
   
        /* create references to the molecules. */
        if ((viz_molp=(struct molecule ***)malloc(n_species*sizeof(struct molecule **)))==NULL) {
		fprintf(log_file, "Memory allocation error\n");
      		return(1);
    	}
    	if ((viz_mol_count=(u_int *)malloc(n_species*sizeof(u_int)))==NULL) {
		fprintf(log_file, "Memory allocation error\n");
      		return(1);
    	}
    }

  if(viz_mol_pos_states_flag){	 
  
    for (ii=0;ii<n_species;ii++) {
      /* perform initialization */
      spec_id=species_list[ii]->species_id;
      viz_molp[spec_id]=NULL;
      viz_mol_count[spec_id]=0;

      if (species_list[ii]->viz_state!=EXCLUDE_OBJ) {

        num=species_list[ii]->population;
        if ((num>0) && (species_list[ii]->flags & NOT_FREE) == 0) {
          /* create references for 3D molecules */
          if ((viz_molp[spec_id]=(struct molecule **)malloc
            (num*sizeof(struct molecule *)))==NULL) {
		fprintf(log_file, "Memory allocation error\n");
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
                }
              }
            }
            amp=amp->next;
          } /* end while */
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
      }

      slp=slp->next;
    }

    for (ii=0;ii<n_species;ii++) {
      
      if(species_list[ii]->viz_state == EXCLUDE_OBJ) continue; 
      spec_id=species_list[ii]->species_id;
      
      if(viz_mol_count[spec_id] > 0)
      {
         num=viz_mol_count[spec_id];
         if (state!=EXCLUDE_OBJ
             && num!=species_list[ii]->population
             && ((species_list[ii]->flags & NOT_FREE)==0) 
             && ((species_list[ii]->flags & ON_SURFACE)==0)) {
             fprintf(log_file,"MCell: molecule count disagreement!!\n");
             fprintf(log_file,"  Species %s  population = %d  count = %d\n",species_list[ii]->sym->name,species_list[ii]->population,num);
         }
      }
      if (viz_mol_count[spec_id]>0 && ((species_list[ii]->flags & NOT_FREE) == 0) && ((species_list[ii]->flags & ON_SURFACE) == 0)) {
        /* here are 3D diffusing molecules */
        num = viz_mol_count[spec_id];
      	  
	  mol_pos_byte_offset_prev = mol_pos_byte_offset;
          for (jj=0;jj<num;jj++) {
            molp=viz_molp[spec_id][jj];
	    v1=world->length_unit*molp->pos.x;
	    v2=world->length_unit*molp->pos.y;
	    v3=world->length_unit*molp->pos.z;
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
          
          /* write molecule orientations information. 
             for 3D molecules we put the number of items in the array as 0. */
      	  mol_orient_byte_offset_prev = mol_orient_byte_offset;
          int num_orient_mol = 0;
          fprintf(master_header,"object %d class array type float rank 1 shape 3 items %d  # %s orientations #\n\n",main_index,num_orient_mol, species_list[ii]->sym->name);
          mol_orient[mol_orient_index] = main_index;
          mol_orient_index++;
          main_index++;

      	  mol_states_byte_offset_prev = mol_states_byte_offset;
          fwrite(&state,sizeof state,1,mol_states_data);
          mol_states_byte_offset += (sizeof state);
          fprintf(master_header,"object %d class constantarray type int items %d %s binary data file %s,%d  # %s states #\n",main_index,num, my_byte_order,mol_states_name,mol_states_byte_offset_prev, species_list[ii]->sym->name);
          fprintf(master_header,"\tattribute \"dep\" string \"positions\"\n\n");
          mol_states[mol_states_index] = main_index;
          mol_states_index++;
          main_index++;
      }
      /* output empty arrays for zero molecule counts here */
      else if ((viz_mol_count[spec_id]==0) && ((species_list[ii]->flags & NOT_FREE) == 0)) {

          fprintf(master_header,"object %d array   # %s positions #\n",main_index, species_list[ii]->sym->name);
          mol_pos[mol_pos_index] = main_index;
          mol_pos_index++;
          main_index++;
          fprintf(master_header,"object %d array   # %s orientations #\n",main_index, species_list[ii]->sym->name);
          mol_orient[mol_orient_index] = main_index;
          mol_orient_index++;
          main_index++;


          fprintf(master_header,"object %d array   # %s states #\n",main_index, species_list[ii]->sym->name);
          mol_states[mol_states_index] = main_index;
          mol_states_index++;
          main_index++;
          
	 fprintf(master_header,"\n");
      }
    }
  } /* end if(viz_mol_pos_states_flag) */
 

/* build fields here */
      int show_molecules = 0;
      if(viz_mol_pos_states_flag ){
		show_molecules = 1;
      } 

    if(show_molecules)
    {
      for (ii=0; ii<mol_to_show_number; ii++) {

             fprintf(master_header,
                "object %d field   # %s #\n",main_index, mol_names[ii]);
             if(mol_pos[ii] > 0) {
             	fprintf(master_header,
                 	"\tcomponent \"positions\" value %d\n",mol_pos[ii]);
             }
             if(mol_orient[ii] > 0){
             	fprintf(master_header,
                 	"\tcomponent \"data\" value %d # orientations #\n",mol_orient[ii]);
             }
             if(mol_states[ii] > 0){ 
             	fprintf(master_header,
                	"\tcomponent \"state_values\" value %d\n",mol_states[ii]);
             }
             fprintf(master_header, "\n");
             /* put number of this field into the linked list 
                of the field numbers for molecules */
             newNode = (struct num_expr_list *)malloc(sizeof (struct num_expr_list));
	     if(newNode == NULL) return 1;
             newNode->value = (double)main_index;
             newNode->next = NULL;    

             if((mol_head == NULL) && (mol_end == NULL)){
		mol_head = mol_end = newNode;
	     }else{
		mol_end->next = newNode;
		mol_end = newNode;
             }
             main_index++;
      }
    }

      /* create group object for molecules/effectors */
      if(show_molecules || show_effectors) {

        if(show_molecules){
                fprintf(master_header,"object %d group # %s #\n",main_index, "3D molecules");
        	mol_group_index = main_index;
        }else if(show_effectors){
                fprintf(master_header,"object %d group # %s #\n",main_index, "effectors");
               eff_group_index = main_index;
        }
        main_index++;

        if(show_effectors)
        {
                curr_ptr = eff_head;
      		for (ii=0;ii<eff_to_show_number;ii++) {
                	if(curr_ptr != NULL){
                		group_index = (int)curr_ptr->value;
                	}else{
				break;
                        }
          		fprintf(master_header,"\tmember \"%s\" value %d\n",eff_names[ii], group_index);
             		curr_ptr = curr_ptr->next;		
          	}
      	}
        if(show_molecules)
        {
                curr_ptr = mol_head;
      		for (ii=0;ii<mol_to_show_number;ii++) {
                	if(curr_ptr != NULL){
                		group_index = (int)curr_ptr->value;
                	}else{
				break;
                        }
          		fprintf(master_header,"\tmember \"%s\" value %d\n",mol_names[ii], group_index);
             		curr_ptr = curr_ptr->next;		
          	}
      	}
        fprintf(master_header, "\n\n"); 
      }	
 

      /* create combined group object for meshes and molecules */
      int show_combined_group = 1;
      if(((fdlp->viz_iterationll == curr_surf_pos_iteration_step) ||
	  (fdlp->viz_iterationll == curr_surf_states_iteration_step)) && 
         (mesh_group_index == 0))
      {
            	show_combined_group = 0;
      }
      if(((fdlp->viz_iterationll == curr_eff_pos_iteration_step) ||
	 (fdlp->viz_iterationll == curr_eff_states_iteration_step))
             && (eff_group_index == 0)){
            	show_combined_group = 0;
      }
      if((fdlp->viz_iterationll == curr_mol_pos_states_iteration_step) 
             && (mol_group_index == 0)){
            	show_combined_group = 0;
      }
      

      if(show_combined_group)
      {
        fprintf(master_header,"object %d group\n",main_index);
        combined_group_index = main_index;
      	if(mesh_group_index > 0){
          	fprintf(master_header,"\tmember \"meshes\" value %d\n",mesh_group_index);
		ia_int_store(&frame_numbers_mesh, frame_numbers_count,mesh_group_index);
		mesh_group_index = 0;

      	}else{
		ia_int_store(&frame_numbers_mesh, frame_numbers_count,-1);
        }
      	if(mol_group_index > 0){
          	fprintf(master_header,"\tmember \"molecules\" value %d\n",mol_group_index);
		ia_int_store(&frame_numbers_mol, frame_numbers_count,mol_group_index);
		mol_group_index = 0;
        }else{
		ia_int_store(&frame_numbers_mol, frame_numbers_count,-1);
        }
        if(eff_group_index > 0){
          	fprintf(master_header,"\tmember \"effectors\" value %d\n",eff_group_index);
		ia_int_store(&frame_numbers_eff, frame_numbers_count,eff_group_index);
		eff_group_index = 0;
        }else{
		ia_int_store(&frame_numbers_eff, frame_numbers_count,-1);
        }
      	fprintf(master_header,"\n");
      }
     /* create entry into "frame_data" object. */
      if(show_combined_group){

        /* create an entry into an 'frame_data' object. */
      	sprintf(buffer, "\tmember %d value %d position %lld\n", series_index, combined_group_index, viz_iterationll);
	ia_string_store(&frame_data_list, frame_data_count, buffer);
        frame_data_count++;
      	series_index++;
      	main_index++;

      	fprintf(master_header, "\n\n");
        
        if(special_surf_frames_counter == 2){
		special_surf_frames_counter = 0;
        }
        if(special_eff_frames_counter == 2){
		special_eff_frames_counter = 0;
        }
 
	/* put value of viz_iteration into the frame_numbers_pos */ 
	ia_int_store(&frame_numbers_pos, frame_numbers_count,(int)viz_iterationll);  /*FIXME--viz_iteration is a long long!*/
        frame_numbers_count++;

     }

     if(frames_iterations_count == total_frames_iterations)
     {
	/* write 'frame_data' object. */
    	fprintf(master_header,"object \"frame_data\" class series\n");
        if(frame_data_count > 0)
        {
        	char *elem;
		for(ii = 0; ii < frame_data_count; ii++){
                	elem = ia_string_get(&frame_data_list, ii);
			fprintf(master_header, "\t%s", elem);
        	}
        }
	fprintf(master_header, "\n\n");

	/* write 'frame_numbers' object. */
        if(frame_numbers_count > 0)
        {
        	int elem1, elem2, elem3, elem4;
        	fprintf(master_header,"object \"frame_numbers\" class array  type int rank 1 shape 4 items %d data follows\n",frame_numbers_count);
		for(ii = 0; ii < frame_numbers_count; ii++){
                	elem1 = ia_int_get(&frame_numbers_pos, ii);
                	elem2 = ia_int_get(&frame_numbers_mesh, ii);
                	elem3 = ia_int_get(&frame_numbers_mol, ii);
                	elem4 = ia_int_get(&frame_numbers_eff, ii);
			fprintf(master_header, "\t%d\t%d\t%d\t%d\n", elem1, elem2, elem3, elem4);
        	}
		fprintf(master_header, "\n\n");
        }

        /* free allocated memory */
        if(surf_states != NULL) free(surf_states);
        if(surf_pos != NULL) free(surf_pos);
        if(surf_con != NULL) free(surf_con);
        if(obj_names != NULL) {
            for (ii=0;ii<obj_to_show_number;ii++) {
               if (obj_names[ii]!=NULL) {
                   free(obj_names[ii]);
               }
            }
            free(obj_names);
        }
        if(eff_states != NULL) free(eff_states);
        if(eff_pos != NULL) free(eff_pos);
        if(eff_names != NULL) {
            for (ii=0;ii<eff_to_show_number;ii++) {
               if (eff_names[ii]!=NULL) {
                   free(eff_names[ii]);
               }
            }
            free(eff_names);
        }
        if(mol_states != NULL) free(mol_states);
        if(mol_pos != NULL) free(mol_pos);
        if(mol_names != NULL) {
            for (ii=0;ii<mol_to_show_number;ii++) {
               if (mol_names[ii]!=NULL) {
                   free(mol_names[ii]);
               }
            }
            free(mol_names);
        }

        /* free linked list for frame_numbers. */
  	struct infinite_int_array *arr_int_ptr, *temp_int_ptr;
  	arr_int_ptr = (&frame_numbers_pos)->next;
  	while(arr_int_ptr != NULL){
		temp_int_ptr = arr_int_ptr;
		arr_int_ptr = arr_int_ptr->next;
		free(temp_int_ptr);
  	}     
  	arr_int_ptr = (&frame_numbers_mesh)->next;
  	while(arr_int_ptr != NULL){
		temp_int_ptr = arr_int_ptr;
		arr_int_ptr = arr_int_ptr->next;
		free(temp_int_ptr);
  	}     
  	arr_int_ptr = (&frame_numbers_mol)->next;
  	while(arr_int_ptr != NULL){
		temp_int_ptr = arr_int_ptr;
		arr_int_ptr = arr_int_ptr->next;
		free(temp_int_ptr);
  	}     
  	arr_int_ptr = (&frame_numbers_eff)->next;
  	while(arr_int_ptr != NULL){
		temp_int_ptr = arr_int_ptr;
		arr_int_ptr = arr_int_ptr->next;
		free(temp_int_ptr);
  	}     

  	/* free linked list for frame_data. */
  	struct infinite_string_array *arr_string_ptr, *temp_string_ptr;
  	arr_string_ptr = (&frame_data_list);
        for(ii = 0; ii < BLOCK_SIZE; ii++){
	    if(arr_string_ptr->data[ii] != NULL){
		free(arr_string_ptr->data[ii]);
            }
        }
  	arr_string_ptr = (&frame_data_list)->next;
  	while(arr_string_ptr != NULL){
           for(ii = 0; ii < BLOCK_SIZE; ii++){
	       if(arr_string_ptr->data[ii] != NULL){
		   free(arr_string_ptr->data[ii]);
               }
           }
	   temp_string_ptr = arr_string_ptr;
	   arr_string_ptr = arr_string_ptr->next;
	   free(temp_string_ptr);
  	}
    } 
     
 
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
  
  /* free linked list for surfaces. */
  if(surf_head != NULL){
    curr_ptr = surf_head;
    while(curr_ptr != NULL){
        temp = curr_ptr;
        curr_ptr = curr_ptr->next;
        free(temp);
    }
    surf_head = surf_end = NULL;
  }      

  /* free linked list for effectors. */
  if(eff_head != NULL){
    curr_ptr = eff_head;
    while(curr_ptr != NULL){
        temp = curr_ptr;
        curr_ptr = curr_ptr->next;
        free(temp);
    }
    eff_head = eff_end = NULL;
  }      
  
  /* free linked list for molecules. */
  if(mol_head != NULL){
    curr_ptr = mol_head;
    while(curr_ptr != NULL){
        temp = curr_ptr;
        curr_ptr = curr_ptr->next;
        free(temp);
    }
    mol_head = mol_end = NULL;
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


  return(0);

}

/*************************************************************************
output_dreamm_objects_all_frame_data:
	In: struct frame_data_list *fdlp
	Out: 0 on success, 1 on error; output visualization files (*.dx)
             in dreamm group format are written.
        NB! NOT TESTED FOR SURFACE_MOLECULES!!!
**************************************************************************/

int output_dreamm_objects_all_frame_data(struct frame_data_list *fdlp)
{
  FILE *log_file;
  FILE *master_header = NULL;
  FILE *mesh_pos_data = NULL;  /* data file for wall vertices */
  FILE *mesh_states_data = NULL; /* data file for wall states */
  FILE *mol_pos_data = NULL;	/* data file for molecule positions */
  FILE *mol_states_data = NULL; /* data file for molecule states */
  FILE *mol_orient_data = NULL; /* data file for molecule orientations */
  struct viz_obj *vizp = NULL;
  struct viz_child *vcp = NULL;
  struct surface_grid *sg;
  struct wall *w,**wp;
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
  struct grid_molecule *gmol;
  struct molecule *molp,***viz_molp = NULL;       /* for 3D molecules */
  float v1,v2,v3;
  u_int spec_id, *viz_mol_count, *viz_grid_mol_count;
  static u_int main_index = 1;
  u_int group_index = 0;
  static u_int series_index = 0;
  int ii,jj;
  int vi1,vi2,vi3;
  /* int vi4; */
  int num;
  long long viz_iterationll;
  int n_viz_iterations;
  int element_data_count;
  int state;
  int viz_type;
  unsigned int index;
  /* indices used in arrays "u_int *surf_pos", etc. */
  int surf_pos_index = 0; 
  int surf_con_index = 0;
  int surf_states_index = 0;
  int eff_states_index = 0;
  int eff_pos_index = 0;
  int eff_orient_index = 0;
  int mol_states_index = 0;
  int mol_pos_index = 0;
  int mol_orient_index = 0;
  int word;
  byte *word_p;
  char file_name[1024];
  char *ch_ptr = NULL; /* pointer used to extract data file name */
  char *grid_mol_name = NULL; /* points to the name of the grid molecule */
  char mesh_pos_name[1024]; /* wall vertices data file name */
  char mesh_states_name[1024]; /* wall states data file name */
  char mol_pos_name[1024]; /* molecule_positions data file name */
  char mol_states_name[1024]; /* molecule_positions data file name */
  char mol_orient_name[1024]; /* molecule orientations data file name */
  char buffer[100]; /* used to write 'frame_data' object information */
  /* linked list that stores data for the 'frame_data' object */
  static struct infinite_string_array frame_data_list;
  static int frame_data_count = 0; /* count elements in frame_data_list array. */
  /* linked lists that stores data for the 'frame_numbers' object */
  static struct infinite_int_array frame_numbers_pos;
  static struct infinite_int_array frame_numbers_mesh;
  static struct infinite_int_array frame_numbers_mol;
  static struct infinite_int_array frame_numbers_eff;
  static int frame_numbers_count = 0; /* count elements in frame_numbers_pos array. */
  char my_byte_order[8];  /* shows binary ordering ('lsb' or 'msb') */
  static int mesh_pos_byte_offset = 0;  /* defines position of the object data
                                  in the mesh positions binary data file */
  static int mesh_states_byte_offset = 0; /* defines position of the object data
                                    in the mesh states binary file */
  static int mol_pos_byte_offset = 0; /*defines position of the object data 
                           in the molecule positions binary file */
  static int mol_orient_byte_offset = 0; /*defines position of the object data 
                           in the molecule orientations binary file */
  static int mol_states_byte_offset = 0; /* defines position of the object data
                              in the molecule states binary file. */
  int mol_pos_byte_offset_prev = 0; /*defines position of the object data 
                           in the effector positions binary file */
  int mol_orient_byte_offset_prev = 0; /*defines position of the object data                             in the effector orientations binary file */
  int mol_states_byte_offset_prev = 0; /* defines position of the object                          d ata in the effector states binary file. */

  /* linked list that stores members of the group for surfaces */
  struct num_expr_list *surf_head = NULL;  
  struct num_expr_list *surf_end = NULL;
  /* linked list that stores members of the group for 3D molecules */
  struct num_expr_list *mol_head = NULL;  
  struct num_expr_list *mol_end = NULL;
  /* linked list that stores members of the group for effectors */
  struct num_expr_list *eff_head = NULL;  
  struct num_expr_list *eff_end = NULL;
  
  struct num_expr_list *curr_ptr = NULL;
  struct num_expr_list *newNode = NULL;
  struct num_expr_list *temp = NULL;
  /* placeholders for the members of the combined group */
  int mesh_group_index = 0;
  int mol_group_index = 0;
  int eff_group_index = 0;
  int combined_group_index = 0;
  int first_iteration = 0;      /* flag signaling the first iteration
                                   in the iteration list */


  u_int n_species = world->n_species;
  species_list=world->species_list;
  
  /* these arrays are used to hold values of main_index for objects */
  if(obj_to_show_number > 0)
  {
  	if(surf_states == NULL){ 
		if ((surf_states=(u_int *)malloc(obj_to_show_number*sizeof(u_int)))==NULL)      {
         		fprintf(stderr,"MCell: memory allocation error\n");
         		return (1);
     		}
        	for(ii = 0; ii < obj_to_show_number; ii++){
			surf_states[ii] = 0;
     		}
   	}
   	if(surf_pos == NULL){
     		if ((surf_pos=(u_int *)malloc(obj_to_show_number*sizeof(u_int)))==NULL) {
         		fprintf(stderr,"MCell: memory allocation error\n");
         		return (1);
     		}
     		for(ii = 0; ii < obj_to_show_number; ii++){
			surf_pos[ii] = 0;
     		}
   	}
   	if(surf_con == NULL){
      		if ((surf_con=(u_int *)malloc(obj_to_show_number*sizeof(u_int)))==NULL) {
         		fprintf(stderr,"MCell: memory allocation error\n");
         		return (1);
      		}
      		for(ii = 0; ii < obj_to_show_number; ii++){
			surf_con[ii] = 0;
      		}
   	}
   	/* initialize array of viz_objects names */
   	if(obj_names == NULL){
      		if ((obj_names = (char **)malloc(obj_to_show_number*sizeof(char *)))==NULL) {
         		fprintf(stderr,"MCell: memory allocation error\n");
         		return (1);
      		}
      		for(ii = 0; ii < obj_to_show_number; ii++){
			obj_names[ii] = NULL;
      		}
   	}
  
     	ii = 0;
     	vizp = world->viz_obj_head;
     	while(vizp != NULL){
		vcp = vizp->viz_child_head;
        	while(vcp != NULL){
         		if(obj_names != NULL){
         		   obj_names[ii] = my_strdup(vcp->obj->sym->name);
                           if(obj_names[ii] == NULL){
         		       fprintf(stderr,"MCell: memory allocation error\n");
         		       return (1);
         		   }
                           ii++;
         		}
         		vcp = vcp->next;
       		}
       		vizp = vizp->next;
    	}
  } /* end if (obj_to_show_number > 0) */

  /* these arrays are used to hold values of main_index for effectors 
     and 3D molecules */
  if(eff_to_show_number > 0)
  {
     	if(eff_states == NULL){ 
     		if((eff_states=(u_int *)malloc(eff_to_show_number*sizeof(u_int)))==NULL)        {
         		fprintf(stderr,"MCell: memory allocation error\n");
         		return (1);
     		}
     		for(ii = 0; ii < eff_to_show_number; ii++){
			eff_states[ii] = 0;
     		}
     	}
     	if(eff_pos == NULL){
     		if ((eff_pos=(u_int *)malloc(eff_to_show_number*sizeof(u_int)))==NULL) 	{
         		fprintf(stderr,"MCell: memory allocation error\n");
         		return (1);
     		}
     		for(ii = 0; ii < eff_to_show_number; ii++){
			eff_pos[ii] = 0;
     		}
     	}
     	if(eff_orient == NULL){
     		if ((eff_orient = (u_int *)malloc(eff_to_show_number*sizeof(u_int)))==NULL) 	{
         		fprintf(stderr,"MCell: memory allocation error\n");
         		return (1);
     		}
     		for(ii = 0; ii < eff_to_show_number; ii++){
			eff_orient[ii] = 0;
     		}
     	}
   	/* initialize array of grid_mol's names */
   	if(eff_names == NULL){
      		if ((eff_names = (char **)malloc(eff_to_show_number*sizeof(char *)))==NULL) {
         		fprintf(stderr,"MCell: memory allocation error\n");
         		return (1);
      		}
      		for(ii = 0; ii < eff_to_show_number; ii++){
			eff_names[ii] = NULL;
      		}
   	}
        index = 0;
        for (ii=0;ii<n_species;ii++) {

          if(species_list[ii]->viz_state == EXCLUDE_OBJ) continue;
          if((species_list[ii]->flags & ON_GRID) == ON_GRID) 
          { 
             eff_names[index] = my_strdup(species_list[ii]->sym->name);
             if(eff_names[index] == NULL){
                 fprintf(stderr,"MCell: memory allocation error\n");
                 return (1);
             }
             index++;
          }
        }
        index = 0;
   } /* end if (eff_to_show_number > 0) */
   
  if(mol_to_show_number > 0)
  {
     	if(mol_states == NULL){ 
     		if((mol_states=(u_int *)malloc(mol_to_show_number*sizeof(u_int)))==NULL)        {
         		fprintf(stderr,"MCell: memory allocation error\n");
         		return (1);
     		}
     		for(ii = 0; ii < mol_to_show_number; ii++){
			mol_states[ii] = 0;
     		}
     	}
     	if(mol_pos == NULL){
     		if ((mol_pos=(u_int *)malloc(mol_to_show_number*sizeof(u_int)))==NULL) 	{
         		fprintf(stderr,"MCell: memory allocation error\n");
         		return (1);
     		}
     		for(ii = 0; ii < mol_to_show_number; ii++){
			mol_pos[ii] = 0;
     		}
     	}
     	if(mol_orient == NULL){
     		if ((mol_orient = (u_int *)malloc(mol_to_show_number*sizeof(u_int)))==NULL) 	{
         		fprintf(stderr,"MCell: memory allocation error\n");
         		return (1);
     		}
     		for(ii = 0; ii < mol_to_show_number; ii++){
			mol_orient[ii] = 0;
     		}
     	}
   	/* initialize array of mol's names */
   	if(mol_names == NULL){
      		if ((mol_names = (char **)malloc(mol_to_show_number*sizeof(char *)))==NULL) {
         		fprintf(stderr,"MCell: memory allocation error\n");
         		return (1);
      		}
      		for(ii = 0; ii < mol_to_show_number; ii++){
			mol_names[ii] = NULL;
      		}
   	}
        index = 0;
        for (ii=0;ii<n_species;ii++) {

          if(species_list[ii]->viz_state == EXCLUDE_OBJ) continue;
          
          if((species_list[ii]->flags & NOT_FREE) == 0) 
          { 
             mol_names[index] = my_strdup(species_list[ii]->sym->name);
             if(mol_names[index] == NULL){
                 fprintf(stderr,"MCell: memory allocation error\n");
                 return (1);
             }
             index++;
          }
        }
        index = 0;

  } /* end if (mol_to_show_number > 0) */

  /* check whether it is a first iteration for a given 'fdlp' */
  if(first_iteration == 0){
  	if(world->it_time <= fdlp->iteration_list->value){
  		first_iteration = 1;
  	}
  }

  log_file=world->log_file;
  no_printf("Viz output in DREAMM_V3 mode all frame data...\n");

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

  viz_type=fdlp->type;
  if(viz_type != ALL_FRAME_DATA) {
   	fprintf(log_file, "Wrong viz_type. Expected ALL_FRAME_DATA.  Received %d\n", viz_type);
   	return 1;
  }


  /* Open master header file. */
  sprintf(file_name,"%s.master_header.dx",world->molecule_prefix_name);

  if(first_iteration == 1){
      if ((master_header=fopen(file_name,"w"))==NULL) {
           fprintf(log_file,"MCell: error cannot open master header file %s\n",file_name);
           return(1);
      }

  }else{
      if ((master_header=fopen(file_name,"a"))==NULL) {
        fprintf(log_file,"MCell: error cannot open master header file %s\n",file_name);
        return(1);
      }
  }
  

  sprintf(file_name,"%s.molecule_positions.bin",world->molecule_prefix_name);
  /* remove the folder name from the molecule_positions data file name */
  ch_ptr = strrchr(file_name, '/');
  ++ch_ptr;
  strcpy(mol_pos_name, ch_ptr);
     

  if (first_iteration == 1){
        if ((mol_pos_data=fopen(file_name,"wb"))==NULL) {
           fprintf(log_file,"MCell: error cannot open molecule positions data file %s\n",file_name);
           return(1);
        }else{}
   }else{
        if ((mol_pos_data=fopen(file_name,"ab"))==NULL) {
           fprintf(log_file,"MCell: error cannot open molecule positions data file %s\n",file_name);
           return(1);
        }else{}

   }

     sprintf(file_name,"%s.molecule_orientations.bin",world->molecule_prefix_name);
     /* remove the folder name from the molecule_positions data file name */
     ch_ptr = strrchr(file_name, '/');
     ++ch_ptr;
     strcpy(mol_orient_name, ch_ptr);
     

     if (first_iteration == 1){
        if ((mol_orient_data=fopen(file_name,"wb"))==NULL) {
           fprintf(log_file,"MCell: error cannot open molecule orientations data file %s\n",file_name);
           return(1);
        }else{}
     }else{
        if ((mol_orient_data=fopen(file_name,"ab"))==NULL) {
           fprintf(log_file,"MCell: error cannot open molecule orientations data file %s\n",file_name);
           return(1);
        }else{}

     }


       sprintf(file_name,"%s.molecule_states.bin",world->molecule_prefix_name);
       /* remove the folder name from the molecule_states data file name */
       ch_ptr = strrchr(file_name, '/');
       ++ch_ptr;
       strcpy(mol_states_name, ch_ptr);

     
       if (first_iteration == 1){
        if ((mol_states_data = fopen(file_name,"wb"))==NULL) {
           fprintf(log_file,"MCell: error cannot open molecule states data file %s\n",file_name);
           return(1);
         }else{}
       }else{
        if ((mol_states_data = fopen(file_name,"ab"))==NULL) {
           fprintf(log_file,"MCell: error cannot open molecule states data file %s\n",file_name);
           return(1);
         }else{}

       }
    

/* dump walls */
     vizp = world->viz_obj_head;
     
     while(vizp!=NULL) {
   
      	sprintf(file_name,"%s.mesh_positions.bin",vizp->name);
     
     	/* remove the folder name from the mesh_positions data file name */
     	ch_ptr = strrchr(file_name, '/');
     	++ch_ptr;
     	strcpy(mesh_pos_name, ch_ptr);
     

       if (first_iteration == 1){
      	 if ((mesh_pos_data=fopen(file_name,"wb"))==NULL) {
           fprintf(log_file,"MCell: error cannot open mesh positions file %s\n",
               file_name);
           return(1);
         }else{}
       }else{
      	 if ((mesh_pos_data=fopen(file_name,"ab"))==NULL) {
           fprintf(log_file,"MCell: error cannot open mesh positions file %s\n",
               file_name);
           return(1);
         }else{}
       }


         sprintf(file_name,"%s.mesh_states.bin", vizp->name);
     
        /* remove the folder name from the mesh_states data file name */
        ch_ptr = strrchr(file_name, '/');
        ++ch_ptr;
        strcpy(mesh_states_name, ch_ptr);


       if (first_iteration == 1){
        if ((mesh_states_data=fopen(file_name,"wb"))==NULL) {
           fprintf(log_file,"MCell: error cannot open mesh states file %s\n",
               file_name);
           return(1);
        }else{}
       }else{
          if ((mesh_states_data=fopen(file_name,"ab"))==NULL) {
             fprintf(log_file,"MCell: error cannot open mesh states file %s\n",
               file_name);
             return(1);
          }else{}
       }

    /* Traverse all visualized compartments 
       output mesh element positions and connections */
    

    vcp = vizp->viz_child_head;

    while(vcp!=NULL) {
      objp = vcp->obj;
      pop=(struct polygon_object *)objp->contents;
      if (objp->object_type==BOX_OBJ) {
#if 0
          
            element_data_count=0.5*objp->n_walls_actual;
        
               
            fprintf(master_header,
              "object %d class array type float rank 1 shape 3 items %d %s binary data file %s,%d # %s.positions #\n",
              main_index,objp->n_verts,my_byte_order,mesh_pos_name, mesh_pos_byte_offset, objp->sym->name);
            fprintf(master_header,
              "\tattribute \"dep\" string \"positions\"\n\n");

            /* output box vertices */
            for (ii=0;ii<objp->n_verts;ii++) {
              v1 = world->length_unit*objp->verts[ii].x;
              v2 = world->length_unit*objp->verts[ii].y;
              v3 = world->length_unit*objp->verts[ii].z;
              fwrite(&v1,sizeof v1,1,mesh_pos_data);
              fwrite(&v2,sizeof v2,1,mesh_pos_data);
              fwrite(&v3,sizeof v3,1,mesh_pos_data);
              mesh_pos_byte_offset += (sizeof(v1) + sizeof(v2) + sizeof(v3));
            }
            surf_pos[surf_pos_index] = main_index;
            surf_pos_index++;
	    main_index++;
    
            /* output box wall connections */
            fprintf(master_header,
              "object %d class array type int rank 1 shape 4 items %d %s binary data file %s,%d # %s.connections #\n",
              main_index,element_data_count,my_byte_order, mesh_pos_name, mesh_pos_byte_offset, objp->sym->name);
            fprintf(master_header,
              "\tattribute \"ref\" string \"positions\"\n");
            fprintf(master_header,
              "\tattribute \"element type\" string \"quads\"\n\n");
            surf_con[surf_con_index] = main_index;
            surf_con_index++;
            main_index++;


            fprintf(master_header,
              "object %d class array type int rank 0 items %d %s binary data file %s,%d # %s.states #\n", main_index, element_data_count, my_byte_order, mesh_states_name, mesh_states_byte_offset, objp->sym->name);
            fprintf(master_header,"\tattribute \"dep\" string \"connections\"\n\n");
            surf_states[surf_states_index] = main_index;
            surf_states_index++;
            main_index++;
            
            for (ii=0;ii<objp->n_walls;ii+=2) {
              if (!get_bit(pop->side_removed,ii)) {
                switch (ii) {
                  case TP:
                    vi1=3;
                    vi2=7;
                    vi3=1;
                    vi4=5;
                    fwrite(&vi1,sizeof vi1,1,mesh_pos_data);
                    fwrite(&vi2,sizeof vi2,1,mesh_pos_data);
                    fwrite(&vi3,sizeof vi3,1,mesh_pos_data);
                    fwrite(&vi4,sizeof vi4,1,mesh_pos_data);
                    mesh_pos_byte_offset += (sizeof(vi1) + sizeof(vi2) + sizeof(vi3) + sizeof(vi4));
                  break;
                  case BOT:
                    vi1=0;
                    vi2=4;
                    vi3=2;
                    vi4=6;
                    fwrite(&vi1,sizeof vi1,1,mesh_pos_data);
                    fwrite(&vi2,sizeof vi2,1,mesh_pos_data);
                    fwrite(&vi3,sizeof vi3,1,mesh_pos_data);
                    fwrite(&vi4,sizeof vi4,1,mesh_pos_data);
                    mesh_pos_byte_offset += (sizeof(vi1) + sizeof(vi2) + sizeof(vi3) + sizeof(vi4));
                  break;
                  case FRNT:
                    vi1=4;
                    vi2=0;
                    vi3=5;
                    vi4=1;
                    fwrite(&vi1,sizeof vi1,1,mesh_pos_data);
                    fwrite(&vi2,sizeof vi2,1,mesh_pos_data);
                    fwrite(&vi3,sizeof vi3,1,mesh_pos_data);
                    fwrite(&vi4,sizeof vi4,1,mesh_pos_data);
                    mesh_pos_byte_offset += (sizeof(vi1) + sizeof(vi2) + sizeof(vi3) + sizeof(vi4));
                  break;
                  case BCK:
                    vi1=2;
                    vi2=6;
                    vi3=3;
                    vi4=7;
                    fwrite(&vi1,sizeof vi1,1,mesh_pos_data);
                    fwrite(&vi2,sizeof vi2,1,mesh_pos_data);
                    fwrite(&vi3,sizeof vi3,1,mesh_pos_data);
                    fwrite(&vi4,sizeof vi4,1,mesh_pos_data);
                    mesh_pos_byte_offset += (sizeof(vi1) + sizeof(vi2) + sizeof(vi3) + sizeof(vi4));
                  break;
                  case LFT:
                    vi1=0;
                    vi2=2;
                    vi3=1;
                    vi4=3;
                    fwrite(&vi1,sizeof vi1,1,mesh_pos_data);
                    fwrite(&vi2,sizeof vi2,1,mesh_pos_data);
                    fwrite(&vi3,sizeof vi3,1,mesh_pos_data);
                    fwrite(&vi4,sizeof vi4,1,mesh_pos_data);
                    mesh_pos_byte_offset += (sizeof(vi1) + sizeof(vi2) + sizeof(vi3) + sizeof(vi4));
                  break;
                  case RT:
                    vi1=6;
                    vi2=4;
                    vi3=7;
                    vi4=5;
                    fwrite(&vi1,sizeof vi1,1,mesh_pos_data);
                    fwrite(&vi2,sizeof vi2,1,mesh_pos_data);
                    fwrite(&vi3,sizeof vi3,1,mesh_pos_data);
                    fwrite(&vi4,sizeof vi4,1,mesh_pos_data);
                    mesh_pos_byte_offset += (sizeof(vi1) + sizeof(vi2) + sizeof(vi3) + sizeof(vi4));
                  break;
                }
                state=objp->viz_state[ii];
                fwrite(&state,sizeof (state),1,mesh_states_data);
                mesh_states_byte_offset += sizeof(state);
              }
            }
#endif
      }


      if (objp->object_type==POLY_OBJ || objp->object_type==BOX_OBJ) {
          opp=(struct ordered_poly *)pop->polygon_data;
          edp=opp->element;
          element_data_count=objp->n_walls_actual;


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
              v1 = world->length_unit*objp->verts[ii].x;
              v2 = world->length_unit*objp->verts[ii].y;
              v3 = world->length_unit*objp->verts[ii].z;
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
            surf_con[surf_con_index] = main_index;
            surf_con_index++;
            main_index++;


            fprintf(master_header,
              "object %d class array type int rank 0 items %d %s binary data file %s,%d # %s.states #\n", main_index, element_data_count, my_byte_order, mesh_states_name, mesh_states_byte_offset, objp->sym->name);
            fprintf(master_header,"\tattribute \"dep\" string \"connections\"\n\n");
           surf_states[surf_states_index] = main_index;
           surf_states_index++;
           main_index++;


            for (ii=0;ii<objp->n_walls;ii++) {
              if (!get_bit(pop->side_removed,ii)) {
                for (jj=0;jj<edp[ii].n_verts-2;jj++) {
                  vi1=edp[ii].vertex_index[0];
                  vi2=edp[ii].vertex_index[jj+1];
                  vi3=edp[ii].vertex_index[jj+2];
	          fwrite(&vi1,sizeof vi1,1,mesh_pos_data);
	          fwrite(&vi2,sizeof vi2,1,mesh_pos_data);
	          fwrite(&vi3,sizeof vi3,1,mesh_pos_data);
                  mesh_pos_byte_offset += (sizeof(v1) + sizeof(v2) + sizeof(v3));
                }
              state=objp->viz_state[ii];
              fwrite(&state,sizeof (state),1,mesh_states_data);
              mesh_states_byte_offset += sizeof(state);
              }
            }
      }
      vcp = vcp->next;
      }
      
      vizp = vizp->next;
      }
  
    /* create field objects */
    for(ii = 0; ii < obj_to_show_number; ii++)
    {
             fprintf(master_header,
                "object %d field   # %s #\n",main_index, obj_names[ii]);
             fprintf(master_header,
                 "\tcomponent \"positions\" value %d\n",surf_pos[ii]);
             fprintf(master_header,
                 "\tcomponent \"connections\" value %d\n",surf_con[ii]);
             fprintf(master_header,
                "\tcomponent \"state_values\" value %d\n\n",surf_states[ii]);
	     
	     /* put number of this field into the linked list 
                of the field numbers for meshes */
             newNode = (struct num_expr_list *)malloc(sizeof (struct num_expr_list));
	     if(newNode == NULL) return 1;
             newNode->value = (double)main_index;
             newNode->next = NULL;    

             if((surf_head == NULL) && (surf_end == NULL)){
		surf_head = surf_end = newNode;
	     }else{
		surf_end->next = newNode;
		surf_end = newNode;
             }
             main_index++;
    }     

      
    /* create a group object for all meshes. */
    vizp = world->viz_obj_head;
    if(vizp != NULL){

        	fprintf(master_header,"object %d group # %s #\n",main_index, "meshes");
        	mesh_group_index = main_index;
        	main_index++;

        	curr_ptr = surf_head;
        	while(vizp != NULL) {
         		vcp = vizp->viz_child_head;
              
         		while(vcp != NULL){
                		if(curr_ptr != NULL){
                   	    		group_index = (int)curr_ptr->value;
                		}else{
		   	    		break;
                		}
             			fprintf(master_header,"\tmember \"%s\" value %d\n",vcp->obj->sym->name,group_index);
             	       		curr_ptr = curr_ptr->next;		
		       		vcp = vcp->next;
       		}
       		vizp = vizp->next;
       	}       
       	fprintf(master_header, "\n\n"); 
   }  /* end (if vizp) */
    

   /* Visualize molecules. */

/* dump grid molecules. */
    
    /* create references to the numbers of grid molecules of each name. */
    if ((viz_grid_mol_count=(u_int *)malloc(n_species*sizeof(u_int)))==NULL) {
      return(1);
    }

   for (ii = 0; ii < n_species; ii++)
   {
      /* perform initialization */
      spec_id=species_list[ii]->species_id;
      viz_grid_mol_count[spec_id]=0;

     if((species_list[ii]->flags & ON_GRID) == 0) continue; 
     if(species_list[ii]->viz_state == EXCLUDE_OBJ) continue; 
        
     grid_mol_name = species_list[ii]->sym->name;
    
     mol_pos_byte_offset_prev = mol_pos_byte_offset;
     mol_orient_byte_offset_prev = mol_orient_byte_offset;
     mol_states_byte_offset_prev = mol_states_byte_offset;
      
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
	                if (gmol!=NULL) {
	                   state=sg->mol[index]->properties->viz_state;
	                }else{
			   state = EXCLUDE_OBJ;
                        }
                        if (state!=EXCLUDE_OBJ) {
                          spec_id = gmol->properties->species_id;
                          viz_grid_mol_count[spec_id]++;
                          if(strcmp(gmol->properties->sym->name,
                            grid_mol_name) == 0){
                                   /* write positions information */
	           		   v1=world->length_unit*p0.x;
	           		   v2=world->length_unit*p0.y;
	           		   v3=world->length_unit*p0.z;
	           		   fwrite(&v1,sizeof v1,1,mol_pos_data);
	               		   fwrite(&v2,sizeof v2,1,mol_pos_data);
	           		   fwrite(&v3,sizeof v3,1,mol_pos_data);
                                   mol_pos_byte_offset += (sizeof(v1) + sizeof(v2) + sizeof(v3));
                                   /* write orientations information */
                   		   v1=w->normal.x;
                   		   v2=w->normal.y;
                   		   v3=w->normal.z;
                   		   fwrite(&v1,sizeof v1,1,mol_orient_data);
                   		   fwrite(&v2,sizeof v2,1,mol_orient_data);
                   		   fwrite(&v3,sizeof v3,1,mol_orient_data);
                                   mol_orient_byte_offset += (sizeof(v1) + sizeof(v2) + sizeof(v3));

                          }
                           
                        } /* end if */
                     } /* end for */
                   } /* end if */
                 } /* end if */
              } /* end for */
              vcp = vcp->next;
            } /* end while (vcp) */
                                              
            vizp = vizp->next;
         } /* end while (vzp) */

        num = viz_grid_mol_count[ii];
           
        fprintf(master_header,"object %d class array type float rank 1 shape 3 items %d %s binary data file %s,%d # %s positions #\n",main_index,num,my_byte_order, mol_pos_name, mol_pos_byte_offset_prev, grid_mol_name);
        fprintf(master_header,"\tattribute \"dep\" string \"positions\"\n\n");
        eff_pos[eff_pos_index] = main_index;
        eff_pos_index++;
        main_index++;
   
        fprintf(master_header,"object %d class array type float rank 1 shape 3 items %d %s binary data file %s,%d   # %s orientations #\n",main_index,num,my_byte_order, mol_orient_name, mol_orient_byte_offset_prev, species_list[ii]->sym->name);
        fprintf(master_header,"\tattribute \"dep\" string \"positions\"\n\n");
        eff_orient[eff_orient_index] = main_index;
        eff_orient_index++;
        main_index++;

        /* write states information. */
        fwrite(&state,sizeof state,1,mol_states_data);
        mol_states_byte_offset += (sizeof state);

        fprintf(master_header,"object %d class constantarray type int items %d %s binary data file %s,%d  # %s states #\n",main_index,num, my_byte_order,mol_states_name,mol_states_byte_offset_prev, species_list[ii]->sym->name);

        fprintf(master_header,"\tattribute \"dep\" string \"positions\"\n\n");
        eff_states[eff_states_index] = main_index;
        eff_states_index++;
        main_index++;
  

   } /* end for */

/* build fields for grid molecules here */

    for(ii = 0; ii < eff_to_show_number; ii++)
    {

             fprintf(master_header,
                "object %d field   # %s #\n",main_index, eff_names[ii]);
             fprintf(master_header,
                 "\tcomponent \"positions\" value %d\n",eff_pos[ii]);
             fprintf(master_header,
                 "\tcomponent \"data\" value %d # orientations #\n",eff_orient[ii]);
             fprintf(master_header,
                "\tcomponent \"state_values\" value %d\n\n",eff_states[ii]);
             /* put number of this field into the linked list 
                of the field numbers for grid molecules */
             newNode = (struct num_expr_list *)malloc(sizeof (struct num_expr_list));
	     if(newNode == NULL) return 1;
             newNode->value = (double)main_index;
             newNode->next = NULL;    

             if((eff_head == NULL) && (eff_end == NULL)){
		eff_head = eff_end = newNode;
	     }else{
		eff_end->next = newNode;
		eff_end = newNode;
             }
             main_index++;
   }

      /* create group object for effectors */

        if(eff_to_show_number > 0){
                fprintf(master_header,"object %d group # %s #\n",main_index, "effectors");
               eff_group_index = main_index;
               main_index++;

                curr_ptr = eff_head;
      		for (ii=0;ii<eff_to_show_number;ii++) {
                	if(curr_ptr != NULL){
                		group_index = (int)curr_ptr->value;
                	}else{
				break;
                        }
          		fprintf(master_header,"\tmember \"%s\" value %d\n",eff_names[ii], group_index);
             		curr_ptr = curr_ptr->next;		
          	}
                fprintf(master_header, "\n\n"); 
      	}




/* dump 3D and surface molecules: */

    /* create references to the molecules. */
    if ((viz_molp=(struct molecule ***)malloc(n_species*sizeof(struct molecule **)))==NULL) {
      return(1);
    }
    if ((viz_mol_count=(u_int *)malloc(n_species*sizeof(u_int)))==NULL) {
      return(1);
    }

    for (ii=0;ii<n_species;ii++) {
      /* perform initialization */
      spec_id=species_list[ii]->species_id;
      viz_molp[spec_id]=NULL;
      viz_mol_count[spec_id]=0;

      if (species_list[ii]->viz_state!=EXCLUDE_OBJ) {

        num=species_list[ii]->population;
        if ((num>0) && (species_list[ii]->flags & NOT_FREE) == 0) {
          /* create references for 3D molecules */
          if ((viz_molp[spec_id]=(struct molecule **)malloc
            (num*sizeof(struct molecule *)))==NULL) {
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
                }
              }
            }
            amp=amp->next;
          } /* end while */
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
      }

      slp=slp->next;
    }

    for (ii=0;ii<n_species;ii++) {
      
      if(species_list[ii]->viz_state == EXCLUDE_OBJ) continue; 
      spec_id=species_list[ii]->species_id;
      
      if(viz_mol_count[spec_id] > 0)
      {
         num=viz_mol_count[spec_id];
         if (state!=EXCLUDE_OBJ
             && num!=species_list[ii]->population
             && ((species_list[ii]->flags & NOT_FREE)==0) 
             && ((species_list[ii]->flags & ON_SURFACE)==0)) {
             fprintf(log_file,"MCell: molecule count disagreement!!\n");
             fprintf(log_file,"  Species %s  population = %d  count = %d\n",species_list[ii]->sym->name,species_list[ii]->population,num);
         }
      }
      if (viz_mol_count[spec_id]>0 && ((species_list[ii]->flags & NOT_FREE) == 0) && ((species_list[ii]->flags & ON_SURFACE) == 0)) {

        /* here are 3D diffusing molecules */
        num = viz_mol_count[spec_id];

        mol_pos_byte_offset_prev = mol_pos_byte_offset;
        for (jj=0;jj<num;jj++) {
            molp=viz_molp[spec_id][jj];
	    v1=world->length_unit*molp->pos.x;
	    v2=world->length_unit*molp->pos.y;
	    v3=world->length_unit*molp->pos.z;
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
          
          /* write molecule orientations information. 
             for 3D molecules we put the number of items in the array as 0.*/
      	  mol_orient_byte_offset_prev = mol_orient_byte_offset;
          int num_orient_mol = 0;
          fprintf(master_header,"object %d class array type float rank 1 shape 3 items %d  # %s orientations #\n\n",main_index,num_orient_mol, species_list[ii]->sym->name);
          mol_orient[mol_orient_index] = main_index;
          mol_orient_index++;
          main_index++;

      	  mol_states_byte_offset_prev = mol_states_byte_offset;
          fwrite(&state,sizeof state,1,mol_states_data);
          mol_states_byte_offset += (sizeof state);
          fprintf(master_header,"object %d class constantarray type int items %d %s binary data file %s,%d  # %s states #\n",main_index,num, my_byte_order,mol_states_name,mol_states_byte_offset_prev, species_list[ii]->sym->name);
          fprintf(master_header,"\tattribute \"dep\" string \"positions\"\n\n");
          mol_states[mol_states_index] = main_index;
          mol_states_index++;
          main_index++;
      }
      /* output empty arrays for zero molecule counts here */
      else if ((viz_mol_count[spec_id]==0) && ((species_list[ii]->flags & NOT_FREE) == 0)) {

          fprintf(master_header,"object %d array   # %s positions #\n",main_index, species_list[ii]->sym->name);
          mol_pos[mol_pos_index] = main_index;
          mol_pos_index++;
          main_index++;
          fprintf(master_header,"object %d array   # %s orientations #\n",main_index, species_list[ii]->sym->name);
          mol_orient[mol_orient_index] = main_index;
          mol_orient_index++;
          main_index++;

          fprintf(master_header,"object %d array   # %s states #\n",main_index, species_list[ii]->sym->name);
          mol_states[mol_states_index] = main_index;
          mol_states_index++;
          main_index++;
      }
    }

/* build fields and groups here */
      for (ii=0; ii<mol_to_show_number; ii++) {

             fprintf(master_header,
                "object %d field   # %s #\n",main_index, mol_names[ii]);
             fprintf(master_header,
                 "\tcomponent \"positions\" value %d\n",mol_pos[ii]);
             fprintf(master_header,
                 "\tcomponent \"data\" value %d # orientations #\n",mol_orient[ii]);
             fprintf(master_header,
                "\tcomponent \"state_values\" value %d\n\n",mol_states[ii]);
             /* put number of this field into the linked list 
                of the field numbers for molecules */
             newNode = (struct num_expr_list *)malloc(sizeof (struct num_expr_list));
	     if(newNode == NULL) return 1;
             newNode->value = (double)main_index;
             newNode->next = NULL;    

             if((mol_head == NULL) && (mol_end == NULL)){
		mol_head = mol_end = newNode;
	     }else{
		mol_end->next = newNode;
		mol_end = newNode;
             }
             main_index++;
      }
        
      /* create group object for molecules */
        if(mol_to_show_number > 0){
                fprintf(master_header,"object %d group # %s #\n",main_index, "3D molecules");
        	mol_group_index = main_index;
                main_index++;

                curr_ptr = mol_head;
      		for (ii=0;ii<mol_to_show_number;ii++) {
                	if(curr_ptr != NULL){
                		group_index = (int)curr_ptr->value;
                	}else{
				break;
                        }
          		fprintf(master_header,"\tmember \"%s\" value %d\n",mol_names[ii], group_index);
             		curr_ptr = curr_ptr->next;		
          	}
                fprintf(master_header, "\n\n"); 
      	}

 

      /* create combined group object for meshes and molecules */
        fprintf(master_header,"object %d group\n",main_index);
        combined_group_index = main_index;
      	if(mesh_group_index > 0){
          	fprintf(master_header,"\tmember \"meshes\" value %d\n",mesh_group_index);
		ia_int_store(&frame_numbers_mesh, frame_numbers_count, mesh_group_index);
		mesh_group_index = 0;

      	}else{
		ia_int_store(&frame_numbers_mesh, frame_numbers_count, -1);
        }
        if(eff_group_index > 0){
          	fprintf(master_header,"\tmember \"effectors\" value %d\n",eff_group_index);
		ia_int_store(&frame_numbers_eff, frame_numbers_count, eff_group_index);
		eff_group_index = 0;
        }else{
		ia_int_store(&frame_numbers_eff, frame_numbers_count, -1);
        } 
      	if(mol_group_index > 0){
          	fprintf(master_header,"\tmember \"molecules\" value %d\n",mol_group_index);
		ia_int_store(&frame_numbers_mol, frame_numbers_count, mol_group_index);
		mol_group_index = 0;
        }else{
		ia_int_store(&frame_numbers_mol, frame_numbers_count, -1);
        }
      	fprintf(master_header,"\n");



     /* create entry into "frame_data" object. */
      	sprintf(buffer, "\tmember %d value %d position %lld\n", series_index, combined_group_index, viz_iterationll);
	ia_string_store(&frame_data_list, frame_data_count, buffer);
        frame_data_count++;

        series_index++;
      	main_index++;
      	fprintf(master_header, "\n\n");
     
	/* put value of viz_iteration into the frame_numbers_pos */ 
	ia_int_store(&frame_numbers_pos, frame_numbers_count,(int)viz_iterationll);  /*FIXME--viz_iteration is a long long! */
        frame_numbers_count++;

     if(fdlp->curr_viz_iteration->next == NULL)
     {
	/* write 'frame_data' object. */
    	fprintf(master_header,"object \"frame_data\" class series\n");
        if(frame_data_count > 0)
        {
        	char *elem;
		for(ii = 0; ii < frame_data_count; ii++){
                	elem = ia_string_get(&frame_data_list, ii);
			fprintf(master_header, "\t%s", elem);
        	}
        }
	fprintf(master_header, "\n\n");

	/* write 'frame_numbers' object. */
        if(frame_numbers_count > 0)
        {
        	int elem1, elem2, elem3, elem4;
        	fprintf(master_header,"object \"frame_numbers\" class array  type int rank 0 items %d data follows\n",frame_numbers_count);
		for(ii = 0; ii < frame_numbers_count; ii++){
                	elem1 = ia_int_get(&frame_numbers_pos, ii);
                	elem2 = ia_int_get(&frame_numbers_mesh, ii);
                	elem3 = ia_int_get(&frame_numbers_mol, ii);
                	elem4 = ia_int_get(&frame_numbers_eff, ii);
			fprintf(master_header, "\t%d\t%d\t%d\t%d\n", elem1, elem2, elem3, elem4);
        	}
		fprintf(master_header, "\n\n");
        }

        if(surf_states != NULL) free(surf_states);
        if(surf_pos != NULL) free(surf_pos);
        if(surf_con != NULL) free(surf_con);
        if(obj_names != NULL) {
            for (ii=0;ii<obj_to_show_number;ii++) {
               if (obj_names[ii]!=NULL) {
                   free(obj_names[ii]);
               }
            }
            free(obj_names);
        }
        if(eff_states != NULL) free(eff_states);
        if(eff_pos != NULL) free(eff_pos);
        if(eff_names != NULL) {
            for (ii=0;ii<eff_to_show_number;ii++) {
               if (eff_names[ii]!=NULL) {
                   free(eff_names[ii]);
               }
            }
            free(eff_names);
        }
        if(mol_states != NULL) free(mol_states);
        if(mol_pos != NULL) free(mol_pos);
        if(mol_names != NULL) {
            for (ii=0;ii<mol_to_show_number;ii++) {
               if (mol_names[ii]!=NULL) {
                   free(mol_names[ii]);
               }
            }
            free(mol_names);
        }
        
	/* free linked list for frame_numbers. */
  	struct infinite_int_array *arr_int_ptr, *temp_int_ptr;
  	arr_int_ptr = (&frame_numbers_pos)->next;
  	while(arr_int_ptr != NULL){
		temp_int_ptr = arr_int_ptr;
		arr_int_ptr = arr_int_ptr->next;
		free(temp_int_ptr);
  	}     
  	arr_int_ptr = (&frame_numbers_mesh)->next;
  	while(arr_int_ptr != NULL){
		temp_int_ptr = arr_int_ptr;
		arr_int_ptr = arr_int_ptr->next;
		free(temp_int_ptr);
  	}     
  	arr_int_ptr = (&frame_numbers_mol)->next;
  	while(arr_int_ptr != NULL){
		temp_int_ptr = arr_int_ptr;
		arr_int_ptr = arr_int_ptr->next;
		free(temp_int_ptr);
  	}     
  	arr_int_ptr = (&frame_numbers_eff)->next;
  	while(arr_int_ptr != NULL){
		temp_int_ptr = arr_int_ptr;
		arr_int_ptr = arr_int_ptr->next;
		free(temp_int_ptr);
  	}     

  	/* free linked list for frame_data. */
  	struct infinite_string_array *arr_string_ptr, *temp_string_ptr;
  	arr_string_ptr = (&frame_data_list);
        for(ii = 0; ii < BLOCK_SIZE; ii++){
	    if(arr_string_ptr->data[ii] != NULL){
		free(arr_string_ptr->data[ii]);
            }
        }
  	arr_string_ptr = (&frame_data_list)->next;
  	while(arr_string_ptr != NULL){
           for(ii = 0; ii < BLOCK_SIZE; ii++){
	       if(arr_string_ptr->data[ii] != NULL){
		   free(arr_string_ptr->data[ii]);
               }
           }
	   temp_string_ptr = arr_string_ptr;
	   arr_string_ptr = arr_string_ptr->next;
	   free(temp_string_ptr);
  	}
     }


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
  
  /* free linked list for surfaces. */
  if(surf_head != NULL){
    curr_ptr = surf_head;
    while(curr_ptr != NULL){
        temp = curr_ptr;
        curr_ptr = curr_ptr->next;
        free(temp);
    }
    surf_head = surf_end = NULL;
  }      

  /* free linked list for effectors. */
  if(eff_head != NULL){
    curr_ptr = eff_head;
    while(curr_ptr != NULL){
        temp = curr_ptr;
        curr_ptr = curr_ptr->next;
        free(temp);
    }
    eff_head = eff_end = NULL;
  }      
  
  /* free linked list for molecules */
  if(mol_head != NULL){
    curr_ptr = mol_head;
    while(curr_ptr != NULL){
        temp = curr_ptr;
        curr_ptr = curr_ptr->next;
        free(temp);
    }
    mol_head = mol_end = NULL;
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
  struct surface_molecule *smp;
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
	      else if ((amp->properties->flags & ON_SURFACE)!=0)
	      {
		smp = (struct surface_molecule*)amp;
		where.x = smp->pos.x;
		where.y = smp->pos.y;
		where.z = smp->pos.z;
		orient = smp->orient;
	      }
	      else if ((amp->properties->flags & ON_GRID)!=0)
	      {
		gmp = (struct grid_molecule*)amp;
		grid2xyz(gmp->grid,gmp->grid_index,&where);
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
  struct surface_molecule *smp;
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
    sprintf(cf_format,"%%s.rk.%%0%dlld.dat",ndigits);
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
            else if ((amp->properties->flags & ON_SURFACE)!=0)
            {
              smp = (struct surface_molecule*)amp;
              where.x = smp->pos.x;
              where.y = smp->pos.y;
              where.z = smp->pos.z;
              orient = smp->orient;
            }
            else if ((amp->properties->flags & ON_GRID)!=0)
            {
              gmp = (struct grid_molecule*)amp;
              grid2xyz(gmp->grid,gmp->grid_index,&where);
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


