
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "mcell_structs.h"
#include "grid_util.h"
#include "sched_util.h"
#include "viz_output.h"


extern struct volume *world;


void update_frame_data_list(struct frame_data_list *fdlp)
{
  FILE *log_file;

  log_file=world->log_file;
  
  while (fdlp!=NULL) {
    if(world->it_time==fdlp->viz_iteration) {
      if(world->viz_mode==DX_MODE) {
        output_dx_objects(fdlp);
      }
/*
      if(world->viz_mode==RADIANCE_MODE) {
        output_radiance_objects(fdlp);
      }
      if(world->viz_mode==RAYSHADE_MODE) {
        output_rayshade_objects(fdlp);
      }
      if(world->viz_mode==POVRAY_MODE) {
        output_povray_objects(fdlp);
      }
      if(world->viz_mode==RENDERMAN_MODE) {
        output_renderman_objects(fdlp);
      }
      if(world->viz_mode==IRIT_MODE) {
        output_irit_objects(fdlp);
      }
      if(world->viz_mode==MCELL_MODE) {
        output_mcell_objects(fdlp);
      }
      if(world->voxel_image_mode==1) {
        output_voxel_image(fdlp);
      }
      if(world->voxel_volume_mode==1) {
        output_voxel_volume(fdlp);
      }
*/
      fdlp->curr_viz_iteration=fdlp->curr_viz_iteration->next;
      if (fdlp->curr_viz_iteration!=NULL) {
	switch (fdlp->list_type) {
	case FRAME_NUMBER:
	  fdlp->viz_iteration=(int)fdlp->curr_viz_iteration->value; 
	  break;
	case REAL_TIME:
	  fdlp->viz_iteration=(int)(fdlp->curr_viz_iteration->value/world->time_unit + ROUND_UP);
	  break;
	}
      }
    }
    fdlp=fdlp->next;
  }
  return;
}


void init_frame_data_list(struct frame_data_list *fdlp)
{
  struct num_expr_list *nelp;
  int done;

  while (fdlp!=NULL) {
    fdlp->viz_iteration=-1;
    fdlp->n_viz_iterations=0;
    nelp=fdlp->iteration_list;
    done=0;
    switch (fdlp->list_type) {
    case FRAME_NUMBER:
      while (nelp!=NULL) {
	fdlp->n_viz_iterations++;
	if (!done) {
	  if (nelp->value>=world->start_time) {
          fdlp->viz_iteration=nelp->value;
          fdlp->curr_viz_iteration=nelp;
          done=1;
	  }
	}
	nelp=nelp->next;
      }
      break;
    case REAL_TIME:
      while (nelp!=NULL) {
	fdlp->n_viz_iterations++;
	if (!done) {
	  if (nelp->value>=world->current_start_time) {
	    fdlp->viz_iteration=(int)(nelp->value/world->time_unit+ROUND_UP);
	    fdlp->curr_viz_iteration=nelp;
	    done=1;
	  }
      }
	nelp=nelp->next;
      }
      break;
    }
    fdlp=fdlp->next;
  }
  return;
}


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
struct molecule *molp,***viz_molp;
float v1,v2,v3;
u_int n_tiles,n_species,spec_id,*viz_mol_count;
u_int mol_pos_index,mol_pos_field_index,mol_pos_group_index;
u_int mol_states_index,mol_states_group_index;
int ii,jj;
int vi1,vi2,vi3,vi4;
int num;
int viz_iteration,n_viz_iterations;
int first_viz_iteration;
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

viz_iteration=fdlp->viz_iteration;
n_viz_iterations=fdlp->n_viz_iterations;
first_viz_iteration=(viz_iteration==fdlp->iteration_list->value);

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
      sprintf(file_name,"%s.mesh_elements.%d.dx",
               vizp->name,viz_iteration);
      if ((wall_verts_header=fopen(file_name,"wb"))==NULL) {
        fprintf(log_file,"MCell: error cannot open mesh elements file %s\n",
               file_name);
        return(1);
      }
    }

    if (viz_surf_states) {
      sprintf(file_name,"%s.mesh_element_states.%d.dx",
               vizp->name,viz_iteration);
      if ((wall_states_header=fopen(file_name,"wb"))==NULL) {
        fprintf(log_file,"MCell: error cannot open mesh element states file %s\n",
               file_name);
        return(1);
      }
    }

    if (viz_eff_pos) {
      sprintf(file_name,"%s.effector_site_positions.%d.dx",vizp->name,viz_iteration);
      if ((eff_pos_header=fopen(file_name,"wb"))==NULL) {
        fprintf(log_file,"MCell: error cannot open effector position file %s\n",file_name);
        return(1);
      }
    }
  
    if (viz_eff_states) {
      sprintf(file_name,"%s.effector_site_states.%d.dx",vizp->name,viz_iteration);
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
        if (viz_surf_pos || viz_surf_states) {
          element_data_count=0.5*objp->n_walls;
        }
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
            if (pop->side_stat[ii]) {
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

      if (objp->object_type==POLY_OBJ) {
        if (viz_surf) {
          if (viz_surf_pos || viz_surf_states) {
            opp=(struct ordered_poly *)pop->polygon_data;
            edp=opp->element_data;
            element_data_count=objp->n_walls;
          }
          if (viz_surf_pos && !viz_surf_states) {

	    fprintf(wall_verts_header,
              "object \"%s.positions\" class array type float rank 1 shape 3 items %d %s binary data follows\n",
              objp->sym->name,objp->n_verts,my_byte_order);
            /* output polyhedron vertices */
            for (ii=0;ii<objp->n_walls;ii++) {
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
              if (pop->side_stat[ii]) {
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
              if (pop->side_stat[ii]) {
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
            if (pop->side_stat[ii]) {
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
      sprintf(file_name,"%s.molecule_positions.%d.dx",world->molecule_prefix_name,viz_iteration);
      if ((mol_pos_header=fopen(file_name,"wb"))==NULL) {
        fprintf(log_file,"MCell: error cannot open molecule positions header file %s\n",file_name);
        return(1);
      }
    }
    if (viz_mol_states) {
      sprintf(file_name,"%s.molecule_states.%d.dx",world->molecule_prefix_name,viz_iteration);
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

  for (ii=0;ii<n_species;ii++) {
    if (viz_molp[ii]!=NULL) {
      free(viz_molp[ii]);
    }
  }

  if (viz_molp!=NULL) {
    free(viz_molp);
  }

  return(0);
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


