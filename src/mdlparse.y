%{
  #include <stdio.h> 
  #include <stdlib.h>
  #include <stdarg.h>
  #include <string.h> 
  #include <time.h> 
  #include <math.h>
  #include <float.h>
  #include <limits.h>
  #include <sys/errno.h>
  #include "vector.h"
  #include "strfunc.h"
  #include "mcell_structs.h"
  #include "mem_util.h"
  #include "sym_table.h"
  #include "diffuse_util.h"
  #include "mdlparse_util.h"
  #include "mdlparse.h"

/*
  #include "chkpt.h"
  #include "init.h"
  #include "geom_util.h"
*/

  #define YYPARSE_PARAM mdlparse_params
  #define YYLEX_PARAM mdlparse_params

  #define mdlpvp ((struct mdlparse_vars *)YYPARSE_PARAM)
  #define volp mdlpvp->vol

  #ifdef DEBUG
  #define no_printf printf
  #endif

%}


%union {
int tok;
double dbl;
char *str;
struct sym_table *sym;
struct evnt *evnt;
struct vector3 *vec3;
struct num_expr_list *nel;
struct object *obj;
struct object_list *objl;
struct count_list *cnt;
} 


%{
  #include "mdllex.flex.c"
%}


%pure_parser

%name-prefix="mdl"
%output="mdlparse.bison.c"


%token <tok> ABS
%token <tok> ABSORPTIVE
%token <tok> ACOS
%token <tok> ADD_EFFECTOR
%token <tok> ALL
%token <tok> ALL_ELEMENTS
%token <tok> ALL_EVENTS
%token <tok> ALL_HITS
%token <tok> ASIN
%token <tok> ATAN
%token <tok> BACK
%token <tok> BINDING_POLE
%token <tok> BOTH_POLES
%token <tok> BOTTOM
%token <tok> BOX
%token <tok> CEIL
%token <tok> CHARGE
%token <tok> CHECKPOINT_INFILE
%token <tok> CHECKPOINT_OUTFILE
%token <tok> CHECKPOINT_ITERATIONS
%token <tok> CONCENTRATION
%token <tok> COLOR
%token <tok> COLOR_SIDE
%token <tok> COLOR_EFFECTOR
%token <tok> COMPARTMENT
%token <tok> CORNERS
%token <tok> COS
%token <tok> COUNT
%token <tok> CUMULATE_FOR_EACH_TIME_STEP
%token <tok> DBL_ARROW
%token <tok> DIFFUSION_CONSTANT
%token <tok> DEFINE
%token <tok> DEFINE_EFFECTOR_SITE_POSITIONS
%token <tok> DEFINE_MOLECULE
%token <tok> DEFINE_REACTIONS
%token <tok> DEFINE_RELEASE_PATTERN
%token <tok> DEFINE_SURFACE_REGIONS
%token <tok> DELAY
%token <tok> DENSITY
%token <tok> DX
%token <tok> EFFECTOR
%token <tok> EFFECTOR_GRID_DENSITY
%token <tok> EFFECTOR_POSITIONS
%token <tok> EFFECTOR_STATE 
%token <tok> EFFECTOR_STATES
%token <tok> EITHER_POLE
%token <tok> ELEMENT
%token <tok> ELEMENT_CONNECTIONS
%token <tok> ELEMENT_LIST
%token <tok> EOF_TOK
%token <tok> EXP
%token <tok> EXPRESSION
%token <tok> FALSE
%token <tok> FCLOSE
%token <tok> FLOOR
%token <tok> FOPEN
%token <tok> FOR_EACH_EFFECTOR
%token <tok> FOR_EACH_MOLECULE
%token <tok> FOR_EACH_TIME_STEP
%token <tok> FPRINT_TIME
%token <tok> FPRINTF
%token <tok> FRONT
%token <tok> FULLY_CLOSED
%token <tok> GAUSSIAN_RELEASE_NUMBER
%token <tok> INCLUDE_FILE
%token <tok> INITIAL_EVENTS
%token <tok> INPUT_FILE
%token <tok> INSTANTIATE
%token <tok> INTEGER
%token <tok> INTERIM_EVENTS
%token <tok> IRIT
%token <tok> ITERATION_FRAME_DATA
%token <tok> ITERATION_LIST
%token <tok> ITERATIONS
%token <tok> FULLY_RANDOM
%token <tok> LFT_ARROW
%token <tok> LEFT
%token <tok> MOLECULE
%token <tok> MOLECULE_POSITIONS
%token <tok> MOLECULE_STATES
%token <tok> LOCATION
%token <tok> LOG
%token <tok> LOG10
%token <tok> MAX_TOK
%token <tok> MCELL_GENERIC
%token <tok> MEAN_DIAMETER
%token <tok> MEAN_NUMBER
%token <tok> MIN_TOK
%token <tok> MOD
%token <tok> MODE
%token <tok> MOLECULE_FILE_PREFIX
%token <tok> NAME
%token <tok> NEGATIVE_POLE
%token <tok> NO
%token <tok> NONE
%token <tok> NORMAL
%token <tok> NO_SPECIES
%token <tok> NUMBER
%token <tok> NUMBER_BOUND
%token <tok> NUMBER_OF_TRAINS
%token <tok> NUMBER_TO_RELEASE
%token <tok> OBJECT
%token <tok> OBJECT_FILE_DESIGNATORS
%token <tok> ORIENTATION
%token <tok> PARALLEL_PARTITION
%token <tok> PART
%token <tok> PARTITION_X
%token <tok> PARTITION_Y
%token <tok> PARTITION_Z
%token <tok> PARTS
%token <tok> PI_TOK
%token <tok> POLE_ORIENTATION
%token <tok> POLYGON_LIST
%token <tok> POSITIVE_BACK
%token <tok> POSITIVE_FRONT
%token <tok> POSITIVE_POLE
%token <tok> POVRAY
%token <tok> PRINT_TIME
%token <tok> PRINTF
%token <tok> RADIAL_DIRECTIONS
%token <tok> RADIAL_SUBDIVISIONS
%token <tok> RADIANCE
%token <tok> RAYSHADE
%token <tok> REACTION_DATA_OUTPUT
%token <tok> REACTION_GROUP
%token <tok> REAL
%token <tok> REFERENCE_STATE
%token <tok> REFLECTIVE
%token <tok> REFERENCE_DIFFUSION_CONSTANT
%token <tok> REGION
%token <tok> RELEASE_INTERVAL
%token <tok> RELEASE_PATTERN
%token <tok> RELEASE_PROBABILITY
%token <tok> REMOVE_ELEMENT
%token <tok> RENDERMAN
%token <tok> RIGHT
%token <tok> ROTATE
%token <tok> ROUND_OFF
%token <tok> RT_ARROW
%token <tok> SCALE
%token <tok> SEED
%token <tok> SIN
%token <tok> SITE_DIAMETER
%token <tok> SPECIFIED_EFFECTORS
%token <tok> SPECIFIED_MOLECULES
%token <tok> SPHERICAL_RELEASE_SITE
%token <tok> SPRINTF
%token <tok> SQRT
%token <tok> STANDARD_DEVIATION
%token <tok> STATE
%token <tok> STATE_VALUES
%token <tok> STEP
%token <tok> STR_VALUE
%token <tok> STRING_TO_NUM
%token <tok> SUM_OVER_ALL_EFFECTORS
%token <tok> SUM_OVER_ALL_MOLECULES
%token <tok> SUM_OVER_ALL_TIME_STEPS
%token <tok> SURFACE_POSITIONS
%token <tok> SURFACE_STATES
%token <tok> TAN
%token <tok> TIME_LIST
%token <tok> TIME_STEP
%token <tod> TIME_STEP_MAX
%token <tok> TO
%token <tok> TOP
%token <tok> TRAIN_DURATION
%token <tok> TRAIN_INTERVAL
%token <tok> TRANSFORM
%token <tok> TRANSLATE
%token <tok> TRANSPARENT
%token <tok> TRUE
%token <tok> UNLIMITED
%token <tok> VAR
%token <tok> VERTEX_LIST
%token <tok> VIZ_DATA_OUTPUT
%token <tok> VOLUME_DEPENDENT_RELEASE_NUMBER
%token <tok> VOXEL_IMAGE_MODE
%token <tok> VOXEL_VOLUME_MODE
%token <tok> WORLD
%token <tok> X_TOK
%token <tok> XY_TOK
%token <tok> XZ_TOK
%token <tok> XYZ_TOK
%token <tok> Y_TOK
%token <tok> YES
%token <tok> YZ_TOK
%token <tok> Z_TOK


%type <tok> mdl_format
%type <tok> mdl_stmt_list
%type <tok> mdl_stmt
%type <tok> time_def
%type <tok> time_max_def
%type <tok> iteration_def
%type <tok> grid_density_def
%type <tok> radial_directions_def
%type <tok> radial_subdivisions_def
%type <tok> assignment_stmt
%type <tok> molecule_def
%type <tok> chkpt_stmt
%type <tok> release_pattern_def
%type <tok> physical_object_def
%type <tok> instance_def
%type <tok> include_stmt
%type <tok> end_of_mdl_file

%type <tok> rx_net_def
%type <tok> rx_stmt
%type <tok> list_rx_stmts
%type <tok> rxn
%type <tok> list_rxns
%type <tok> rx_group_def
%type <tok> unimolecular_rxn
%type <tok> bimolecular_rxn
%type <tok> product
%type <tok> list_products

%type <sym> assign_var
%type <sym> existing_var_only
%type <sym> existing_str_var
%type <sym> existing_num_var
%type <sym> existing_array
%type <sym> new_molecule

%type <dbl> num_value
%type <dbl> num_expr
%type <dbl> num_expr_only
%type <dbl> intOrReal 
%type <dbl> diffusion_def
%type <dbl> reference_diffusion_def
%type <dbl> charge_def

%type <str> str_value
%type <str> str_expr
%type <str> str_expr_only
%type <str> file_name

%type <sym> existing_molecule
%type <sym> new_release_pattern
%type <sym> existing_release_pattern
%type <sym> object_def
%type <sym> new_object
%type <sym> existing_object
%type <sym> object_ref
%type <sym> existing_object_ref
%type <sym> meta_object_def
%type <sym> release_site_def

%type <sym> reactant

%type <nel> array_value

%type <vec3> point


/**********************

%type <tok> boolean
%type <tok> orientation_cmd
%type <tok> wall_prop
%type <tok> side
%type <tok> side_name
%type <tok> partition_dimension
%type <tok> orientation
%type <tok> r_spec
%type <tok> lig_spec
%type <tok> t_spec
%type <tok> event_spec
%type <tok> io_stmt
%type <tok> viz_mode_def
%type <tok> viz_output_def
%type <tok> list_viz_output_cmds
%type <tok> viz_output_cmd
%type <tok> viz_iteration_def 
%type <tok> viz_output_list_def
%type <tok> viz_time_def 
%type <tok> viz_object_prefixes_def
%type <tok> list_viz_object_prefixes
%type <tok> viz_object_prefix
%type <tok> viz_molecule_prefix_def
%type <tok> viz_state_values_def
%type <tok> list_viz_state_values
%type <tok> viz_state_value
%type <tok> partition_def
%type <tok> voxel_image_mode_def
%type <tok> voxel_volume_mode_def
%type <tok> surface_region_def
%type <tok> effector_site_def
%type <tok> output_def
%type <tok> add_effector
%type <tok> remove_side
%type <tok> wall_prop_cmd
%type <tok> fopen_stmt
%type <tok> fclose_stmt
%type <tok> printf_stmt
%type <tok> fprintf_stmt
%type <tok> print_time_stmt
%type <tok> fprint_time_stmt
%type <tok> sprintf_stmt
%type <tok> polygon_object_cmd
%type <tok> list_polygon_object_cmds
%type <tok> polarity
%type <tok> polarity_spec
%type <tok> viz_frame_data_def
%type <tok> list_frame_data_specs
%type <tok> frame_data_spec
%type <tok> frame_data_item
%type <tok> parallel_partition_def
%type <tok> partition_plane

%type <dbl> rate
%type <dbl> effector_quantity_cmd

%type <sym> transition
%type <sym> existing_molecule_or_reaction_state
%type <sym> existing_logicalOrPhysical
%type <sym> new_region
%type <sym> existing_region
%type <sym> box_def
%type <sym> polygon_def
%type <sym> reaction_state
%type <sym> existing_reaction_state
%type <sym> new_file_stream
%type <sym> existing_file_stream

%type <str> file_mode
%type <str> format_string

%type <cnt> count_expr
%type <cnt> count_value
%type <cnt> count_value_init

%type <evnt> event 


***************************************/



%right '='
%left '+' '-'
%left '*' '/'
%left '^'
%left UNARYMINUS

%%


mdl_format:
	mdl_stmt_list
;


mdl_stmt_list: mdl_stmt
	| mdl_stmt_list mdl_stmt
;
 

mdl_stmt: time_def
        | time_max_def
	| iteration_def
	| grid_density_def
	| radial_directions_def
	| radial_subdivisions_def
	| assignment_stmt
	| molecule_def
	| rx_net_def
	| chkpt_stmt
	| release_pattern_def
	| physical_object_def
	| instance_def
/*
	| partition_def
	| parallel_partition_def
	| surface_region_def
	| effector_site_def
	| viz_output_def
	| output_def
	| io_stmt
*/
	| include_stmt
	| end_of_mdl_file
;


end_of_mdl_file: EOF_TOK
{
	if (mdlpvp->include_stack_ptr==0) {
          no_printf("terminal end of file.  curr_file = %s\n",volp->curr_file); 
          no_printf("include_flag = %d\n",mdlpvp->include_flag); 
          fflush(stderr);
          if (!mdlpvp->include_flag) {
/*
            partition_volume(volume);
            build_ligand_table();
            if (mdlpvp->vol->procnum == 0) print_rx();
*/
            no_printf("\n***** End of input file %s reached.\n\n",volp->curr_file);
          }
	  return(0);
	}
	else {
          no_printf("intermediate end of file.  curr_file = %s\n",
            volp->curr_file); 
          no_printf("include_flag = %d\n",mdlpvp->include_flag); 
          fflush(stderr);
	  mdl_switch_to_buffer(mdlpvp->include_stack[--mdlpvp->include_stack_ptr]);
	  volp->curr_file=mdlpvp->include_filename[mdlpvp->include_stack_ptr];
          mdlpvp->cval=mdlpvp->cval_stack[mdlpvp->include_stack_ptr]=mdlpvp->cval;
          mdlpvp->cval_2=mdlpvp->cval_2_stack[mdlpvp->include_stack_ptr]=mdlpvp->cval_2;
          no_printf("Switching back to MDL file: %s\n",volp->curr_file);
	}
};


include_stmt: INCLUDE_FILE 
{
	mdlpvp->include_flag = 1;
}
	'=' file_name 
{
  mdlpvp->a_str=$<str>4;
  if (mdlpvp->include_stack_ptr>=MAX_INCLUDE_DEPTH) {
    sprintf(mdlpvp->mdl_err_msg,"%s %s","Includes nested too deeply:",mdlpvp->a_str);
    mdlerror(mdlpvp->mdl_err_msg,mdlpvp);
    return(1);
  }
  if ((mdlin=fopen(mdlpvp->a_str,"r"))==NULL) {
    sprintf(mdlpvp->mdl_err_msg,"%s %s","Cannot open include file:",mdlpvp->a_str);
    mdlerror(mdlpvp->mdl_err_msg,mdlpvp);
    return(1);
  }
  mdlpvp->include_stack[mdlpvp->include_stack_ptr]=YY_CURRENT_BUFFER;
  mdlpvp->cval_stack[mdlpvp->include_stack_ptr]=mdlpvp->cval;
  mdlpvp->cval_2_stack[mdlpvp->include_stack_ptr]=mdlpvp->cval_2;
  mdlpvp->cval=NULL;
  mdlpvp->cval_2=NULL;
  mdlpvp->include_filename[mdlpvp->include_stack_ptr++]=volp->curr_file;
  volp->curr_file=my_strdup(mdlpvp->a_str);
  yyless(0);
  yyclearin;
  mdl_switch_to_buffer(mdl_create_buffer(mdlin,YY_BUF_SIZE));
  no_printf("Now parsing MDL file: %s\n",volp->curr_file);
  fflush(stderr);
  mdlpvp->include_flag = 0;
};


assignment_stmt: assign_var '=' num_expr_only
{
  mdlpvp->gp=$<sym>1;
  if ((mdlpvp->gp->value=(void *)malloc(sizeof(double)))==NULL) {
    sprintf(mdlpvp->mdl_err_msg,"%s %s","Cannot store value in table:",mdlpvp->gp->name);
    mdlerror(mdlpvp->mdl_err_msg,mdlpvp);
    return(1);
  }
  mdlpvp->gp->sym_type=DBL;
  mdlpvp->dblp=(double *)mdlpvp->gp->value;
  *mdlpvp->dblp=$<dbl>3;
  no_printf("\n%s is equal to: %f\n",mdlpvp->gp->name,*(double *)mdlpvp->gp->value);
  fflush(stderr);
}
	| assign_var '=' str_expr_only
{
  mdlpvp->gp=$<sym>1;
  mdlpvp->gp->sym_type=STR;
  mdlpvp->gp->value=(void *)my_strdup($<str>3);
  no_printf("\n%s is equal to: %s\n",mdlpvp->gp->name,(char *)mdlpvp->gp->value);
  fflush(stderr);
}
	| assign_var '=' existing_var_only
{
  mdlpvp->gp=$<sym>1;
  mdlpvp->tp=$<sym>3;
  switch (mdlpvp->tp->sym_type) {
  case DBL:
    if ((mdlpvp->gp->value=(void *)malloc(sizeof(double)))==NULL) {
      sprintf(mdlpvp->mdl_err_msg,"%s %s","Cannot store value in table:",mdlpvp->gp->name);
      mdlerror(mdlpvp->mdl_err_msg,mdlpvp);
      return(1);
    }
    mdlpvp->gp->sym_type=DBL;
    mdlpvp->dblp=(double *)mdlpvp->gp->value;
    *mdlpvp->dblp=*(double *)mdlpvp->tp->value;
    no_printf("\n%s is equal to: %f\n",mdlpvp->gp->name,*(double *)mdlpvp->gp->value);
    fflush(stderr);
    break;
  case STR:
    mdlpvp->gp->sym_type=STR;
    mdlpvp->gp->value=(void *)my_strdup((char *)mdlpvp->tp->value);
    no_printf("\n%s is equal to: %s\n",mdlpvp->gp->name,(char *)mdlpvp->gp->value);
    fflush(stderr);
    break;
  case ARRAY:
    mdlpvp->gp->sym_type=ARRAY;
    mdlpvp->gp->value=mdlpvp->tp->value;
#ifdef DEBUG
    mdlpvp->elp=(struct num_expr_list *)mdlpvp->gp->value;
    no_printf("\n%s is equal to: [",mdlpvp->gp->name);
    fflush(stderr);
    while (mdlpvp->elp!=NULL) {
      no_printf("%f",mdlpvp->elp->value);
      fflush(stderr);
      mdlpvp->elp=mdlpvp->elp->next;
      if (mdlpvp->elp!=NULL) {
        no_printf(",");
        fflush(stderr);
      }
    }
    no_printf("]\n");
    fflush(stderr);
#endif
    break;
  }
}
	| assign_var '=' '[' 
{
  mdlpvp->num_pos=0;
  mdlpvp->el_head=NULL;
  mdlpvp->el_tail=NULL;
}
array_expr ']'
{
  mdlpvp->gp=$<sym>1;
  mdlpvp->gp->sym_type=ARRAY;
  mdlpvp->gp->value=(void *)mdlpvp->el_head;
#ifdef DEBUG
  mdlpvp->elp=(struct num_expr_list *)mdlpvp->gp->value;
  no_printf("\n%s is equal to: [",mdlpvp->gp->name);
  fflush(stderr);
  while (mdlpvp->elp!=NULL) {
    no_printf("%f",mdlpvp->elp->value);
    fflush(stderr);
    mdlpvp->elp=mdlpvp->elp->next;
    if (mdlpvp->elp!=NULL) {
      no_printf(",");
      fflush(stderr);
    }
  }
  no_printf("]\n");
  fflush(stderr);
#endif
};


array_value: '['
{
  mdlpvp->num_pos=0;
  mdlpvp->el_head=NULL;
  mdlpvp->el_tail=NULL;
}
array_expr ']'
{
#ifdef DEBUG
  mdlpvp->elp=mdlpvp->el_head;
  no_printf("\nArray expression: [");
  while (mdlpvp->elp!=NULL) {
    no_printf("%f",mdlpvp->elp->value);
    mdlpvp->elp=mdlpvp->elp->next;
    if (mdlpvp->elp!=NULL) {
      no_printf(",");
    }
  }
  no_printf("]\n");
#endif
  $$=mdlpvp->el_head;
}
	| existing_array
{
  mdlpvp->gp=$<sym>1;
  $$=(struct num_expr_list *)mdlpvp->gp->value;
};


existing_array: VAR
{
  if (mdlpvp->cval_2!=NULL) {
    mdlpvp->sym_name=mdlpvp->cval_2;
  }
  else {
    mdlpvp->sym_name=mdlpvp->cval;
  }
  if ((mdlpvp->gp=retrieve_sym(mdlpvp->sym_name,ARRAY,volp->main_sym_table))==NULL) {
    sprintf(mdlpvp->mdl_err_msg,"%s %s","Undefined array variable:",mdlpvp->sym_name);
    mdlerror(mdlpvp->mdl_err_msg,mdlpvp);
    if (mdlpvp->sym_name==mdlpvp->cval) {
      mdlpvp->cval=NULL;
    }
    else {
      mdlpvp->cval_2=NULL;
    }
    free((void *)mdlpvp->sym_name);
    return(1);
  }
  if (mdlpvp->sym_name==mdlpvp->cval) {
    mdlpvp->cval=NULL;
  }
  else {
    mdlpvp->cval_2=NULL;
  }
  free((void *)mdlpvp->sym_name);
#ifdef KELP
  mdlpvp->gp->ref_count++;
  no_printf("ref_count: %d\n",mdlpvp->gp->ref_count);
#endif
  $$=mdlpvp->gp;
};


assign_var: VAR
{
  if (mdlpvp->cval_2!=NULL) {
    mdlpvp->sym_name=mdlpvp->cval_2;
  }
  else {
    mdlpvp->sym_name=mdlpvp->cval;
  }
  if ((mdlpvp->gp=retrieve_sym(mdlpvp->sym_name,DBL,volp->main_sym_table))==NULL) {
    if ((mdlpvp->gp=retrieve_sym(mdlpvp->sym_name,STR,volp->main_sym_table))==NULL) {
      if ((mdlpvp->gp=retrieve_sym(mdlpvp->sym_name,ARRAY,volp->main_sym_table))==NULL) {
        if ((mdlpvp->gp=store_sym(mdlpvp->sym_name,TMP,volp->main_sym_table))==NULL) {
          sprintf(mdlpvp->mdl_err_msg,"%s %s","Cannot store variable in table:",mdlpvp->sym_name);
          mdlerror(mdlpvp->mdl_err_msg,mdlpvp);
	  if (mdlpvp->sym_name==mdlpvp->cval) {
	    mdlpvp->cval=NULL;
          }
	  else {
	    mdlpvp->cval_2=NULL;
	  }
          free((void *)mdlpvp->sym_name);
          return(1);
        }
      }
      else {
        if (mdlpvp->sym_name==mdlpvp->cval) {
          mdlpvp->cval=NULL;
        }
        else {
          mdlpvp->cval_2=NULL;
        }
        free((void *)mdlpvp->sym_name);
      }
    }
    else {
      if (mdlpvp->sym_name==mdlpvp->cval) {
        mdlpvp->cval=NULL;
      }
      else {
        mdlpvp->cval_2=NULL;
      }
      free((void *)mdlpvp->sym_name);
    }
  }
  else {
    if (mdlpvp->sym_name==mdlpvp->cval) {
      mdlpvp->cval=NULL;
    }
    else {
      mdlpvp->cval_2=NULL;
    }
    free((void *)mdlpvp->sym_name);
  }
  if (mdlpvp->sym_name==mdlpvp->cval) {
    mdlpvp->cval=NULL;
  }
  else {
    mdlpvp->cval_2=NULL;
  }
  $$=mdlpvp->gp;
};


existing_var_only: VAR
{
  if (mdlpvp->cval_2!=NULL) {
    mdlpvp->sym_name=mdlpvp->cval_2;
  }
  else {
    mdlpvp->sym_name=mdlpvp->cval;
  }
  if ((mdlpvp->gp=retrieve_sym(mdlpvp->sym_name,DBL,volp->main_sym_table))==NULL) {
    if ((mdlpvp->gp=retrieve_sym(mdlpvp->sym_name,STR,volp->main_sym_table))==NULL) {
      if ((mdlpvp->gp=retrieve_sym(mdlpvp->sym_name,ARRAY,volp->main_sym_table))==NULL) {
        sprintf(mdlpvp->mdl_err_msg,"%s %s","Undefined variable:",mdlpvp->sym_name);
        mdlerror(mdlpvp->mdl_err_msg,mdlpvp);
        if (mdlpvp->sym_name==mdlpvp->cval) {
          mdlpvp->cval=NULL;
        }
        else {
          mdlpvp->cval_2=NULL;
        }
        free((void *)mdlpvp->sym_name);
        return(1);
      }
    }
  }
  if (mdlpvp->sym_name==mdlpvp->cval) {
    mdlpvp->cval=NULL;
  }
  else {
    mdlpvp->cval_2=NULL;
  }
  free((void *)mdlpvp->sym_name);
#ifdef KELP
  mdlpvp->gp->ref_count++;
  no_printf("ref_count: %d\n",mdlpvp->gp->ref_count);
#endif
  $$=mdlpvp->gp;
};


num_expr: num_value
{
  $$=$<dbl>1;
}
	| '(' num_expr ')'
{
  $$=$<dbl>2;
}
	| EXP '(' num_expr ')'
{
  $$=exp($<dbl>3);
}
	| LOG '(' num_expr ')'
{
  $$=log($<dbl>3);
}
	| LOG10 '(' num_expr ')'
{
  $$=log10($<dbl>3);
}
	| MAX_TOK '(' num_expr ',' num_expr ')'
{
  mdlpvp->val_1=$<dbl>3;
  mdlpvp->val_2=$<dbl>5;
  if (mdlpvp->val_1 > mdlpvp->val_2) {
    $$=mdlpvp->val_1;
  }
  else {
    $$=mdlpvp->val_2;
  }
}
	| MIN_TOK '(' num_expr ',' num_expr ')'
{
  mdlpvp->val_1=$<dbl>3;
  mdlpvp->val_2=$<dbl>5;
  if (mdlpvp->val_1 < mdlpvp->val_2) {
    $$=mdlpvp->val_1;
  }
  else {
    $$=mdlpvp->val_2;
  }
}
	| ROUND_OFF '(' num_expr ',' num_expr ')'
{
  mdlpvp->val_1=$<dbl>3;
  mdlpvp->val_2=$<dbl>5;
  sprintf(mdlpvp->format_str,"%%.%dg",(int)mdlpvp->val_1);
  sprintf(mdlpvp->str_buf,mdlpvp->format_str,mdlpvp->val_2);
  $$=strtod(mdlpvp->str_buf,(char **)NULL);
}
	|  FLOOR '(' num_expr ')'
{
  $$=floor($<dbl>3);
}
	|  CEIL '(' num_expr ')'
{
  $$=ceil($<dbl>3);
}
	| SIN '(' num_expr ')'
{
  $$=sin($<dbl>3);
}
	| COS '(' num_expr ')'
{
  $$=cos($<dbl>3);
}
	| TAN '(' num_expr ')'
{
  $$=tan($<dbl>3);
}
	| ASIN '(' num_expr ')'
{
  $$=asin($<dbl>3);
}
	| ACOS '(' num_expr ')'
{
  $$=acos($<dbl>3);
}
	| ATAN '(' num_expr ')'
{
  $$=atan($<dbl>3);
}
	| SQRT '(' num_expr ')'
{
  $$=sqrt($<dbl>3);
}
	| ABS '(' num_expr ')'
{
  $$=fabs($<dbl>3);
}
	| MOD '(' num_expr ',' num_expr ')'
{
  $$=fmod($<dbl>3,$<dbl>5);
}
	| PI_TOK
{
  $$=MY_PI;
}
	| SEED
{
  $$=volp->seed_seq;
}
	| STRING_TO_NUM '(' str_expr ')'
{
  $$=strtod($<str>3,(char **)NULL);
  if (errno==ERANGE) {
    sprintf(mdlpvp->mdl_err_msg,"%s %s","Error converting string to number:",$<str>3);
    mdlerror(mdlpvp->mdl_err_msg,mdlpvp);
  }
}
	| num_expr '+' num_expr
{
  $$=$<dbl>1+$<dbl>3;
}
	| num_expr '-' num_expr
{
  $$=$<dbl>1-$<dbl>3;
}
	| num_expr '*' num_expr
{
  $$=$<dbl>1*$<dbl>3;
}
	| num_expr '/' num_expr
{
  $$=$<dbl>1/$<dbl>3;
}
	| num_expr '^' num_expr
{
  $$=pow($<dbl>1,$<dbl>3);
}
	| '-' num_expr %prec UNARYMINUS
{
  $$=-$<dbl>2;
};


num_expr_only: intOrReal
{
  $$=$<dbl>1;
}
	| '(' num_expr ')'
{
  $$=$<dbl>2;
}
	| EXP '(' num_expr ')'
{
  $$=exp($<dbl>3);
}
	| LOG '(' num_expr ')'
{
  $$=log($<dbl>3);
}
	| LOG10 '(' num_expr ')'
{
  $$=log10($<dbl>3);
}
	| MAX_TOK '(' num_expr ',' num_expr ')'
{
  mdlpvp->val_1=$<dbl>3;
  mdlpvp->val_2=$<dbl>5;
  if (mdlpvp->val_1 > mdlpvp->val_2) {
    $$=mdlpvp->val_1;
  }
  else {
    $$=mdlpvp->val_2;
  }
}
	| MIN_TOK '(' num_expr ',' num_expr ')'
{
  mdlpvp->val_1=$<dbl>3;
  mdlpvp->val_2=$<dbl>5;
  if (mdlpvp->val_1 < mdlpvp->val_2) {
    $$=mdlpvp->val_1;
  }
  else {
    $$=mdlpvp->val_2;
  }
}
	| ROUND_OFF '(' num_expr ',' num_expr ')'
{
  mdlpvp->val_1=$<dbl>3;
  mdlpvp->val_2=$<dbl>5;
  sprintf(mdlpvp->format_str,"%%.%dg",(int)mdlpvp->val_1);
  sprintf(mdlpvp->str_buf,mdlpvp->format_str,mdlpvp->val_2);
  $$=strtod(mdlpvp->str_buf,(char **)NULL);
}
	|  FLOOR '(' num_expr ')'
{
  $$=floor($<dbl>3);
}
	|  CEIL '(' num_expr ')'
{
  $$=ceil($<dbl>3);
}
	| SIN '(' num_expr ')'
{
  $$=sin($<dbl>3);
}
	| COS '(' num_expr ')'
{
  $$=cos($<dbl>3);
}
	| TAN '(' num_expr ')'
{
  $$=tan($<dbl>3);
}
	| ASIN '(' num_expr ')'
{
  $$=asin($<dbl>3);
}
	| ACOS '(' num_expr ')'
{
  $$=acos($<dbl>3);
}
	| ATAN '(' num_expr ')'
{
  $$=atan($<dbl>3);
}
	| SQRT '(' num_expr ')'
{
  $$=sqrt($<dbl>3);
}
	| ABS '(' num_expr ')'
{
  $$=fabs($<dbl>3);
}
	| MOD '(' num_expr ',' num_expr ')'
{
  $$=fmod($<dbl>3,$<dbl>5);
}
	| PI_TOK
{
  $$=MY_PI;
}
	| SEED
{
  $$=volp->seed_seq;
}
	| STRING_TO_NUM '(' str_expr ')'
{
  $$=strtod($<str>3,(char **)NULL);
  if (errno==ERANGE) {
    sprintf(mdlpvp->mdl_err_msg,"%s %s","Error converting string to number:",$<str>3);
    mdlerror(mdlpvp->mdl_err_msg,mdlpvp);
  }
}
	| num_expr '+' num_expr
{
  $$=$<dbl>1+$<dbl>3;
}
	| num_expr '-' num_expr
{
  $$=$<dbl>1-$<dbl>3;
}
	| num_expr '*' num_expr
{
  $$=$<dbl>1*$<dbl>3;
}
	| num_expr '/' num_expr
{
  $$=$<dbl>1/$<dbl>3;
}
	| num_expr '^' num_expr
{
  $$=pow($<dbl>1,$<dbl>3);
}
	| '-' num_expr %prec UNARYMINUS
{
  $$=-$<dbl>2;
};


 num_value:  intOrReal {$$=$<dbl>1;}  
	| existing_num_var
{
  mdlpvp->gp=$<sym>1;
  $$=*(double *)mdlpvp->gp->value;
};


existing_num_var: VAR
{
  if (mdlpvp->cval_2!=NULL) {
    mdlpvp->sym_name=mdlpvp->cval_2;
  }
  else {
    mdlpvp->sym_name=mdlpvp->cval;
  }
  if ((mdlpvp->gp=retrieve_sym(mdlpvp->sym_name,DBL,volp->main_sym_table))==NULL) {
    sprintf(mdlpvp->mdl_err_msg,"%s %s","Undefined numeric variable:",mdlpvp->sym_name);
    mdlerror(mdlpvp->mdl_err_msg,mdlpvp);
    if (mdlpvp->sym_name==mdlpvp->cval) {
      mdlpvp->cval=NULL;
    }
    else {
      mdlpvp->cval_2=NULL;
    }
    free((void *)mdlpvp->sym_name);
    return(1);
  }
  if (mdlpvp->sym_name==mdlpvp->cval) {
    mdlpvp->cval=NULL;
  }
  else {
    mdlpvp->cval_2=NULL;
  }
  free((void *)mdlpvp->sym_name);
#ifdef KELP
  mdlpvp->gp->ref_count++;
  no_printf("ref_count: %d\n",mdlpvp->gp->ref_count);
#endif
  $$=mdlpvp->gp;
};


intOrReal: INTEGER {$$=(double)mdlpvp->ival;}
	| REAL {$$=mdlpvp->rval;};

file_name: str_expr
{
  strcpy(mdlpvp->str_buf,$<str>1);
  $$=$<str>1;
};


str_expr: str_value
{
  $$=$<str>1;
}
	| str_expr '&' str_expr
{
  $$=my_strcat($<str>1,$<str>3);
};


str_value: STR_VALUE
{
  $$=strip_quotes(mdlpvp->strval);
  free(mdlpvp->strval);
}
	| INPUT_FILE
{
  $$=volp->mdl_infile_name;
}
	| existing_str_var
{
  mdlpvp->gp=$<sym>1;
  $$=(char *)mdlpvp->gp->value;
};


existing_str_var: VAR
{
  if (mdlpvp->cval_2!=NULL) {
    mdlpvp->sym_name=mdlpvp->cval_2;
  }
  else {
    mdlpvp->sym_name=mdlpvp->cval;
  }
  if ((mdlpvp->gp=retrieve_sym(mdlpvp->sym_name,STR,volp->main_sym_table))==NULL) {
    sprintf(mdlpvp->mdl_err_msg,"%s %s","Undefined text variable:",mdlpvp->sym_name);
    mdlerror(mdlpvp->mdl_err_msg,mdlpvp);
    if (mdlpvp->sym_name==mdlpvp->cval) {
      mdlpvp->cval=NULL;
    }
    else {
      mdlpvp->cval_2=NULL;
    }
    free((void *)mdlpvp->sym_name);
    return(1);
  }
  if (mdlpvp->sym_name==mdlpvp->cval) {
    mdlpvp->cval=NULL;
  }
  else {
    mdlpvp->cval_2=NULL;
  }
  free((void *)mdlpvp->sym_name);
#ifdef KELP
  mdlpvp->gp->ref_count++;
  no_printf("ref_count: %d\n",mdlpvp->gp->ref_count);
#endif
  $$=mdlpvp->gp;
};


str_expr_only: STR_VALUE
{
  $$=strip_quotes(mdlpvp->strval);
  free(mdlpvp->strval);
}
	| INPUT_FILE
{
  $$=volp->mdl_infile_name;
}
	| str_expr '&' str_expr
{
  $$=my_strcat($<str>1,$<str>3);
};


array_expr:  num_expr
{
  if ((mdlpvp->elp=(struct num_expr_list *)malloc
      (sizeof(struct num_expr_list)))==NULL) {
    mdlerror("Cannot store numerical array data");
    return(1);
  }
  mdlpvp->elp->value=$<dbl>1;
  if (mdlpvp->el_tail==NULL) {
    mdlpvp->el_tail=mdlpvp->elp;
  }
  mdlpvp->el_tail->next=mdlpvp->elp;
  mdlpvp->elp->next=NULL;
  mdlpvp->el_tail=mdlpvp->elp;
  if (mdlpvp->el_head==NULL) {
    mdlpvp->el_head=mdlpvp->elp;
  }
  mdlpvp->num_pos++;
}
	| array_expr ',' num_expr
{
  if ((mdlpvp->elp=(struct num_expr_list *)malloc
      (sizeof(struct num_expr_list)))==NULL) {
    mdlerror("Cannot store numerical array data");
    return(1);
  }
  mdlpvp->elp->value=$<dbl>3;
  if (mdlpvp->el_tail==NULL) {
    mdlpvp->el_tail=mdlpvp->elp;
  }
  mdlpvp->el_tail->next=mdlpvp->elp;
  mdlpvp->elp->next=NULL;
  mdlpvp->el_tail=mdlpvp->elp;
  if (mdlpvp->el_head==NULL) {
    mdlpvp->el_head=mdlpvp->elp;
  }
  mdlpvp->num_pos++;
};


time_def: TIME_STEP '=' num_expr
{
volp->time_unit=$<dbl>3;
if (volp->time_unit<0) {
  sprintf(mdlpvp->mdl_err_msg,"Time unit = %g\n\tSetting to %g\n",volp->time_unit,-volp->time_unit);
  mdl_warning(mdlpvp);
  volp->time_unit=-volp->time_unit;
}
no_printf("Time unit = %g\n",volp->time_unit);
fflush(stderr);
};

time_max_def: TIME_STEP_MAX '=' num_expr
{
volp->time_step_max = $<dbl>3;
if (volp->time_step_max<0) {
  sprintf(mdlpvp->mdl_err_msg,"Maximum time step = %g\n\tSetting to %g\n",volp->time_unit,-volp->time_step_max);
  mdl_warning(mdlpvp);
  volp->time_step_max=-volp->time_step_max;
}
no_printf("Maximum time step = %g\n",volp->time_step_max);
fflush(stderr);
};

iteration_def: ITERATIONS '=' num_expr
{
if (volp->iterations==0) {
  volp->iterations=(int) $<dbl>3;
}
no_printf("Iterations = %d\n",volp->iterations);
fflush(stderr);
};


radial_directions_def: RADIAL_DIRECTIONS '=' num_expr
{
  volp->radial_directions = (int) $<dbl>3;
  volp->num_directions=0;
  if (volp->d_step!=NULL) {
    free(volp->d_step);
  }
  if ((volp->d_step=init_d_step(volp->radial_directions,&volp->num_directions))==NULL) {
    mdlerror("Cannot store d_step data for molecule");
    return(1);
  }
  volp->r_num_directions=1.0/volp->num_directions;
  no_printf("desired radial directions = %d\n",volp->radial_directions);
  no_printf("actual radial directions = %d\n",volp->num_directions);
}
	| RADIAL_DIRECTIONS '=' FULLY_RANDOM
{
	volp->fully_random=1;
};


radial_subdivisions_def: RADIAL_SUBDIVISIONS '=' num_expr
{
  volp->radial_subdivisions = (int) $<dbl>3;
  if (volp->r_step!=NULL) {
    free(volp->r_step);
  }
  if ((volp->r_step=init_r_step(volp->radial_subdivisions))==NULL) {
    mdlerror("Cannot store r_step data for molecule");
    return(1);
  }
  no_printf("radial subdivisions = %d\n",volp->radial_subdivisions);
};


grid_density_def: EFFECTOR_GRID_DENSITY '=' num_expr
{
  volp->effector_grid_density=$<dbl>3;
  no_printf("Max density = %f\n",volp->effector_grid_density);
  volp->length_unit=1.0/sqrt(volp->effector_grid_density);
  no_printf("Length unit = %f\n",volp->length_unit);
  mdlpvp->mc_factor=1.0e11*volp->effector_grid_density*sqrt(MY_PI*volp->time_unit)/N_AV;
  mdlpvp->transport_mc_factor=6.2415e18*mdlpvp->mc_factor;
  fflush(stderr);
};


molecule_def: DEFINE_MOLECULE new_molecule '{'
	reference_diffusion_def
	diffusion_def
        charge_def
'}'
{
  mdlpvp->gp=$<sym>2;
  mdlpvp->specp=(struct species *)mdlpvp->gp->value;
  mdlpvp->specp->sym=mdlpvp->gp;
  /* points to the head of region counter list of the transition count for each ligand */
  mdlpvp->specp->D_ref=$<dbl>4;
  mdlpvp->specp->D=$<dbl>5;
  mdlpvp->specp->charge=(int) $<dbl>6;
  if (volp->time_unit==0) {
    sprintf(mdlpvp->mdl_err_msg,"%s %s","TIME_STEP not yet specified.  Cannot define molecule:",mdlpvp->specp->sym->name);
    mdlerror(mdlpvp->mdl_err_msg,mdlpvp);
    return(1);
  }
  if (mdlpvp->specp->D_ref==0) {
    mdlpvp->specp->D_ref=mdlpvp->specp->D;
  }
  mdlpvp->mc_factor=1.0e11*volp->effector_grid_density*sqrt(MY_PI*volp->time_unit)/N_AV;
  mdlpvp->transport_mc_factor=6.2415e18*mdlpvp->mc_factor;
  mdlpvp->specp->space_step=sqrt(4.0*1.0e8*mdlpvp->specp->D*volp->time_unit)/volp->length_unit;
  if (volp->r_step==NULL) {
    if ((volp->r_step=init_r_step(volp->radial_subdivisions))==NULL) {
      mdlerror("Cannot store r_step data for molecule");
      return(1);
    }
  }
  if (volp->d_step==NULL) {
    if ((volp->d_step=init_d_step(volp->radial_directions,&volp->num_directions))==NULL) {
      mdlerror("Cannot store d_step data for molecule");
      return(1);
    }
  }
/*
  if (mdlpvp->specp->space_step*r_step[volp->radial_subdivisions-1]>max_diffusion_step) {
    max_diffusion_step=mdlpvp->specp->space_step*r_step[volp->radial_subdivisions-1];
  }
*/
  mdlpvp->l_perp_bar=sqrt(4*1.0e8*mdlpvp->specp->D*volp->time_unit/MY_PI);
  mdlpvp->l_perp_rms=sqrt(2*1.0e8*mdlpvp->specp->D*volp->time_unit);
  mdlpvp->l_r_bar=2*mdlpvp->l_perp_bar;
  mdlpvp->l_r_rms=sqrt(6*1.0e8*mdlpvp->specp->D*volp->time_unit);
  if (volp->procnum == 0) {
    fprintf(volp->log_file,"\nMCell: Theoretical average diffusion distances for molecule %s:\n\n",mdlpvp->specp->sym->name);
    fprintf(volp->log_file,"\tl_r_bar = %.9g microns\n",mdlpvp->l_r_bar);
    fprintf(volp->log_file,"\tl_r_rms = %.9g microns\n",mdlpvp->l_r_rms);
    fprintf(volp->log_file,"\tl_perp_bar = %.9g microns\n",mdlpvp->l_perp_bar);
    fprintf(volp->log_file,"\tl_perp_rms = %.9g microns\n\n",mdlpvp->l_perp_rms);
  }
  no_printf("Molecule %s defined with D = %g\n",mdlpvp->specp->sym->name,mdlpvp->specp->D);
};


new_molecule: VAR
{
  if (mdlpvp->cval_2!=NULL) {
    mdlpvp->sym_name=mdlpvp->cval_2;
  }
  else {
    mdlpvp->sym_name=mdlpvp->cval;
  }
  if ((mdlpvp->gp=retrieve_sym(mdlpvp->sym_name,MOL,volp->main_sym_table))==NULL) {
    if ((mdlpvp->gp=store_sym(mdlpvp->sym_name,MOL,volp->main_sym_table))==NULL) {
      sprintf(mdlpvp->mdl_err_msg,"%s %s","Cannot store molecule in table:",mdlpvp->sym_name);
      mdlerror(mdlpvp->mdl_err_msg,mdlpvp);
      if (mdlpvp->sym_name==mdlpvp->cval) {
        mdlpvp->cval=NULL;
      }
      else {
        mdlpvp->cval_2=NULL;
      }
      free((void *)mdlpvp->sym_name);
      return(1);
    }
    else {
      mdlpvp->specp=(struct species *)mdlpvp->gp->value;
    }
  }
  else {
    sprintf(mdlpvp->mdl_err_msg,"%s %s","Molecule already defined:",mdlpvp->sym_name);
    mdlerror(mdlpvp->mdl_err_msg,mdlpvp);
    if (mdlpvp->sym_name==mdlpvp->cval) {
      mdlpvp->cval=NULL;
    }
    else {
      mdlpvp->cval_2=NULL;
    }
    free((void *)mdlpvp->sym_name);
    return(1);
  }
  if (mdlpvp->sym_name==mdlpvp->cval) {
    mdlpvp->cval=NULL;
  }
  else {
    mdlpvp->cval_2=NULL;
  }
  $$=mdlpvp->gp;
};


diffusion_def: DIFFUSION_CONSTANT '=' num_expr
{
  $$=$<dbl>3;
};


reference_diffusion_def: /* empty */
{
  $$=0;
}
	| REFERENCE_DIFFUSION_CONSTANT '=' num_expr
{
  $$=$<dbl>3;
};


charge_def: /* empty */
{
  $$=1;
}
	| CHARGE '=' num_expr
{
  $$=$<dbl>3;
};


chkpt_stmt: CHECKPOINT_INFILE '=' file_name
{
  if ( volp->chkpt_infile == NULL ) {
    volp->chkpt_infile=$<str>3;
    if ((mdlpvp->file=fopen(volp->chkpt_infile,"r"))==NULL) {
      volp->chkpt_init=1;
    }
    else {
      volp->chkpt_init=0;
      fclose(mdlpvp->file);
    }
    volp->chkpt_flag = 1;
  }
}
	| CHECKPOINT_OUTFILE '=' file_name
{
  volp->chkpt_outfile=$<str>3;
  volp->chkpt_flag = 1;
}
	| CHECKPOINT_ITERATIONS '=' num_expr
{
  volp->chkpt_iterations=(int) $<dbl>3;
  volp->chkpt_flag = 1;
};


release_pattern_def: DEFINE_RELEASE_PATTERN new_release_pattern '{'
{
	mdlpvp->gp=$<sym>2;
	mdlpvp->rpatp=(struct release_pattern *)mdlpvp->gp->value;
}
	list_req_release_pattern_cmds
	'}'
{
  no_printf("Release pattern %s defined:\n",mdlpvp->gp->name);
  no_printf("\tdelay = %f\n",mdlpvp->rpatp->delay);
  no_printf("\trelease_interval = %f\n",mdlpvp->rpatp->release_interval);
  no_printf("\ttrain_interval = %f\n",mdlpvp->rpatp->train_interval);
  no_printf("\ttrain_duration = %f\n",mdlpvp->rpatp->train_duration);
  no_printf("\tnumber_of_trains = %d\n",mdlpvp->rpatp->number_of_trains);
};

new_release_pattern: VAR
{
  if (mdlpvp->cval_2!=NULL) {
    mdlpvp->sym_name=mdlpvp->cval_2;
  }
  else {
    mdlpvp->sym_name=mdlpvp->cval;
  }
  if ((mdlpvp->gp=retrieve_sym(mdlpvp->sym_name,RPAT,volp->main_sym_table))==NULL) {
    if ((mdlpvp->gp=store_sym(mdlpvp->sym_name,RPAT,volp->main_sym_table))==NULL) {
      sprintf(mdlpvp->mdl_err_msg,"%s %s","Cannot store release pattern in table:",mdlpvp->sym_name);
      mdlerror(mdlpvp->mdl_err_msg,mdlpvp);
      if (mdlpvp->sym_name==mdlpvp->cval) {
        mdlpvp->cval=NULL;
      }
      else {
        mdlpvp->cval_2=NULL;
      }
      free((void *)mdlpvp->sym_name);
      return(1);
    }
  }
  else {
    sprintf(mdlpvp->mdl_err_msg,"%s %s","Release pattern already defined:",mdlpvp->sym_name);
    mdlerror(mdlpvp->mdl_err_msg,mdlpvp);
    if (mdlpvp->sym_name==mdlpvp->cval) {
      mdlpvp->cval=NULL;
    }
    else {
      mdlpvp->cval_2=NULL;
    }
    free((void *)mdlpvp->sym_name);
    return(1);
  }
  if (mdlpvp->sym_name==mdlpvp->cval) {
    mdlpvp->cval=NULL;
  }
  else {
    mdlpvp->cval_2=NULL;
  }
  $$=mdlpvp->gp;
};

existing_release_pattern: VAR
{
  if (mdlpvp->cval_2!=NULL) {
    mdlpvp->sym_name=mdlpvp->cval_2;
  }
  else {
    mdlpvp->sym_name=mdlpvp->cval;
  }
  if ((mdlpvp->gp=retrieve_sym(mdlpvp->sym_name,RPAT,volp->main_sym_table))==NULL) {
    sprintf(mdlpvp->mdl_err_msg,"%s %s","Undefined release pattern:",mdlpvp->sym_name);
    mdlerror(mdlpvp->mdl_err_msg,mdlpvp);
    if (mdlpvp->sym_name==mdlpvp->cval) {
      mdlpvp->cval=NULL;
    }
    else {
      mdlpvp->cval_2=NULL;
    }
    free((void *)mdlpvp->sym_name);
    return(1);
  }
  if (mdlpvp->sym_name==mdlpvp->cval) {
    mdlpvp->cval=NULL;
  }
  else {
    mdlpvp->cval_2=NULL;
  }
  free((void *)mdlpvp->sym_name);
#ifdef KELP
  mdlpvp->gp->ref_count++;
  no_printf("ref_count: %d\n",mdlpvp->gp->ref_count);
#endif
  $$=mdlpvp->gp;
};

list_req_release_pattern_cmds: req_release_pattern_cmd
	| list_req_release_pattern_cmds req_release_pattern_cmd
;

req_release_pattern_cmd:
	DELAY '=' num_expr
{
  mdlpvp->rpatp->delay=$<dbl>3/volp->time_unit;
}
	| RELEASE_INTERVAL '=' num_expr
{
  mdlpvp->rpatp->release_interval=$<dbl>3/volp->time_unit;
}
	| TRAIN_INTERVAL '=' num_expr
{
  mdlpvp->rpatp->train_interval=$<dbl>3/volp->time_unit;
}
	| TRAIN_DURATION '=' num_expr
{
  mdlpvp->rpatp->train_duration=$<dbl>3/volp->time_unit;
}
	| NUMBER_OF_TRAINS '=' num_expr
{
  mdlpvp->rpatp->number_of_trains=(int) $<dbl>3;
}
	| NUMBER_OF_TRAINS '=' UNLIMITED
{
  mdlpvp->rpatp->number_of_trains=INT_MAX;
};


physical_object_def: object_def
{
  mdlpvp->gp=$<sym>1;
  mdlpvp->objp=(struct object *)mdlpvp->gp->value;
  mdlpvp->curr_obj=volp->root_object;
  if (mdlpvp->curr_obj->first_child==NULL) {
    mdlpvp->curr_obj->first_child=mdlpvp->objp;
  }
  if (mdlpvp->curr_obj->last_child!=NULL) {
    mdlpvp->curr_obj->last_child->next=mdlpvp->objp;
  }
  mdlpvp->curr_obj->last_child=mdlpvp->objp;
  mdlpvp->objp->parent=mdlpvp->curr_obj;
  mdlpvp->objp->next=NULL;
};

object_def:
        meta_object_def
{
        $$=$<sym>1;
}
        | release_site_def
{
        $$=$<sym>1;
}
/*
        | box_def
{
        $$=$<sym>1;
}
        | polygon_def
{
        $$=$<sym>1;
}
*/
;

meta_object_def:
        new_object OBJECT '{'
{
        no_printf("Here 1\n");
        fflush(stderr);
	mdlpvp->gp=$<sym>1;
        mdlpvp->objp=(struct object *)mdlpvp->gp->value;
        mdlpvp->objp->object_type=META_OBJ;
        mdlpvp->objp->parent=mdlpvp->curr_obj;
        mdlpvp->curr_obj=mdlpvp->objp;
        no_printf("Here 2\n");
        fflush(stderr);
}
        list_objects
        list_opt_object_cmds
        '}'
{
        mdlpvp->curr_obj=mdlpvp->curr_obj->parent;
        if (mdlpvp->object_name_list_end->prev!=NULL) {
          mdlpvp->object_name_list_end=mdlpvp->object_name_list_end->prev;
        }
        else {
          mdlpvp->object_name_list_end->name=NULL;
        }
        $$=$<sym>1;
};

new_object: VAR
{
  if (mdlpvp->cval_2!=NULL) {
    mdlpvp->sym_name=mdlpvp->cval_2;
  }
  else {
    mdlpvp->sym_name=mdlpvp->cval;
  }
  if (mdlpvp->object_name_list==NULL) {
    if ((mdlpvp->object_name_list=(struct name_list *)malloc
        (sizeof(struct name_list)))==NULL) {
      sprintf(mdlpvp->mdl_err_msg,"%s %s","Cannot store object name:",mdlpvp->sym_name);
      mdlerror(mdlpvp->mdl_err_msg,mdlpvp);
      return(1);
    }
    mdlpvp->object_name_list->name=NULL;
    mdlpvp->object_name_list->prev=NULL;
    mdlpvp->object_name_list->next=NULL;
    mdlpvp->object_name_list_end=mdlpvp->object_name_list;
  }
  if ((mdlpvp->object_name_list_end=concat_obj_name(mdlpvp->object_name_list_end,mdlpvp->sym_name))==NULL) {
    sprintf(mdlpvp->mdl_err_msg,"%s %s","Cannot store object name:",mdlpvp->sym_name);
    return(1);
  }
  mdlpvp->obj_name=mdlpvp->object_name_list_end->name;
  if ((mdlpvp->gp=retrieve_sym(mdlpvp->obj_name,OBJ,volp->main_sym_table))==NULL) {
    if ((mdlpvp->gp=store_sym(mdlpvp->obj_name,OBJ,volp->main_sym_table))==NULL) {
      sprintf(mdlpvp->mdl_err_msg,"%s %s","Cannot store object in table:",mdlpvp->obj_name);
      mdlerror(mdlpvp->mdl_err_msg,mdlpvp);
      if (mdlpvp->sym_name==mdlpvp->cval) {
        mdlpvp->cval=NULL;
      }
      else {
        mdlpvp->cval_2=NULL;
      }
      free((void *)mdlpvp->sym_name);
      return(1);
    }
  }
  else {
    sprintf(mdlpvp->mdl_err_msg,"%s %s","Object already defined:",mdlpvp->obj_name);
    mdlerror(mdlpvp->mdl_err_msg,mdlpvp);
    if (mdlpvp->sym_name==mdlpvp->cval) {
      mdlpvp->cval=NULL;
    }
    else {
      mdlpvp->cval_2=NULL;
    }
    free((void *)mdlpvp->sym_name);
    return(1);
  }
  if (mdlpvp->sym_name==mdlpvp->cval) {
    mdlpvp->cval=NULL;
  }
  else {
    mdlpvp->cval_2=NULL;
  }
  mdlpvp->objp=(struct object *)mdlpvp->gp->value;
  mdlpvp->objp->last_name=mdlpvp->sym_name;
  no_printf("Creating new object: %s\n",mdlpvp->obj_name);
  fflush(stderr);
  $$=mdlpvp->gp;
};

list_opt_object_cmds: opt_object_cmd | list_opt_object_cmds opt_object_cmd;

opt_object_cmd:
        transformation
        | /* empty */
{
};

transformation:
        TRANSLATE '=' point
{
        mdlpvp->pntp1=$<vec3>3;
        init_matrix(mdlpvp->tm);
        translate_matrix(mdlpvp->tm,mdlpvp->tm,mdlpvp->pntp1);
        mult_matrix(mdlpvp->curr_obj->t_matrix,mdlpvp->tm,mdlpvp->curr_obj->t_matrix,4,4,4);
}
        | SCALE '=' point
{
        mdlpvp->pntp1=$<vec3>3;
        init_matrix(mdlpvp->tm);
        scale_matrix(mdlpvp->tm,mdlpvp->tm,mdlpvp->pntp1);
        mult_matrix(mdlpvp->curr_obj->t_matrix,mdlpvp->tm,mdlpvp->curr_obj->t_matrix,4,4,4);
}
        | ROTATE '=' point ',' num_expr
{
        mdlpvp->pntp1=$<vec3>3;
        init_matrix(mdlpvp->tm);
        rotate_matrix(mdlpvp->tm,mdlpvp->tm,mdlpvp->pntp1,$<dbl>5);
        mult_matrix(mdlpvp->curr_obj->t_matrix,mdlpvp->tm,mdlpvp->curr_obj->t_matrix,4,4,4);
};

list_objects:
        object_ref
{
  mdlpvp->gp=$<sym>1;
  mdlpvp->objp=(struct object *)mdlpvp->gp->value;
  if (mdlpvp->curr_obj->first_child==NULL) {
    mdlpvp->curr_obj->first_child=mdlpvp->objp;
  }
  if (mdlpvp->curr_obj->last_child!=NULL) {
    mdlpvp->curr_obj->last_child->next=mdlpvp->objp;
  }
  mdlpvp->curr_obj->last_child=mdlpvp->objp;
  mdlpvp->objp->parent=mdlpvp->curr_obj;
  mdlpvp->objp->next=NULL;
}
        | list_objects object_ref
{
  mdlpvp->gp=$<sym>2;
  mdlpvp->objp=(struct object *)mdlpvp->gp->value;
  if (mdlpvp->curr_obj->first_child==NULL) {
    mdlpvp->curr_obj->first_child=mdlpvp->objp;
  }
  if (mdlpvp->curr_obj->last_child!=NULL) {
    mdlpvp->curr_obj->last_child->next=mdlpvp->objp;
  }
  mdlpvp->curr_obj->last_child=mdlpvp->objp;
  mdlpvp->objp->parent=mdlpvp->curr_obj;
  mdlpvp->objp->next=NULL;
};

object_ref:
        existing_object_ref
{
  $$=$<sym>1;
}
        | object_def
{
  $$=$<sym>1;
};

existing_object_ref:
        new_object OBJECT existing_object '{'
{
  mdlpvp->gp=$<sym>1;
  mdlpvp->tp=$<sym>3;
  mdlpvp->objp=(struct object *)mdlpvp->gp->value;
  mdlpvp->objp2=(struct object *)mdlpvp->tp->value;
  mdlpvp->objp->object_type=mdlpvp->objp2->object_type;
  /* replicate all of object tree rooted at mdlpvp->objp2 */
  if (copy_object(volp,mdlpvp->curr_obj,mdlpvp->objp,mdlpvp->objp2,mdlpvp->mdl_err_msg)) {
    mdlerror(mdlpvp->mdl_err_msg,mdlpvp);
    return(1);
  }
  mdlpvp->objp->parent=mdlpvp->curr_obj;
  mdlpvp->curr_obj=mdlpvp->objp;
}
          list_opt_object_cmds
        '}'
{
  mdlpvp->curr_obj=mdlpvp->curr_obj->parent;
  if (mdlpvp->object_name_list_end->prev!=NULL) {
    mdlpvp->object_name_list_end=mdlpvp->object_name_list_end->prev;
  }
  else {
    mdlpvp->object_name_list_end->name=NULL;
  }
  $$=$<sym>1;
};

release_site_def: new_object SPHERICAL_RELEASE_SITE '{'
{
  mdlpvp->gp=$<sym>1;
  mdlpvp->objp=(struct object *)mdlpvp->gp->value;
  if ((mdlpvp->rsop=(struct release_site_obj *)malloc
              (sizeof(struct release_site_obj)))==NULL) {
    mdlerror("Cannot store release site data");
    return(1);
  }
  mdlpvp->rsop->location=NULL;
  mdlpvp->rsop->mol_type=NULL;
  mdlpvp->rsop->release_number_method=CONSTNUM;
  mdlpvp->rsop->release_number=0;
  mdlpvp->rsop->mean_number=0;
  mdlpvp->rsop->mean_diameter=0;
  mdlpvp->rsop->concentration=0;
  mdlpvp->rsop->standard_deviation=0;
  mdlpvp->rsop->diameter=0;
  mdlpvp->rsop->release_prob=1.0;
  mdlpvp->rsop->pattern=volp->default_release_pattern;
  mdlpvp->objp->object_type=REL_SITE_OBJ;
  mdlpvp->objp->contents=mdlpvp->rsop;
  mdlpvp->objp->parent=mdlpvp->curr_obj;
  mdlpvp->curr_obj=mdlpvp->objp;
}
	list_release_site_cmds
	list_opt_object_cmds
	'}'
{
  no_printf("Release site %s defined:\n",mdlpvp->curr_obj->sym->name);
  no_printf("\tLocation = [%f,%f,%f]\n",mdlpvp->rsop->location->x,mdlpvp->rsop->location->y,mdlpvp->rsop->location->z);
  mdlpvp->curr_obj=mdlpvp->curr_obj->parent;
  if (mdlpvp->object_name_list_end->prev!=NULL) {
    mdlpvp->object_name_list_end=mdlpvp->object_name_list_end->prev;
  }
  else {
    mdlpvp->object_name_list_end->name=NULL;
  }
  $$=$<sym>1;
};

list_release_site_cmds: release_site_cmd
	| list_release_site_cmds release_site_cmd
;

release_site_cmd:
	LOCATION '=' point
{
  mdlpvp->rsop->location=$<vec3>3;
}
	| MOLECULE '=' existing_molecule
{
  mdlpvp->gp=$<sym>3;
  mdlpvp->rsop->mol_type=(struct species *)mdlpvp->gp->value;
}
	| release_number_cmd
	| SITE_DIAMETER '=' num_expr
{
  mdlpvp->rsop->diameter=$<dbl>3/volp->length_unit;
}
	| RELEASE_PROBABILITY '=' num_expr
{
  mdlpvp->rsop->release_prob=$<dbl>3;
}
	| RELEASE_PATTERN '=' existing_release_pattern
{
  mdlpvp->gp=$<sym>3;
  mdlpvp->rsop->pattern=(struct release_pattern *)mdlpvp->gp->value;
};

release_number_cmd: constant_release_number_cmd
	| gaussian_release_number_cmd
	| volume_dependent_number_cmd
;

constant_release_number_cmd: NUMBER_TO_RELEASE '=' num_expr
{
  mdlpvp->rsop->release_number_method=CONSTNUM;
  mdlpvp->rsop->release_number=(int) $<dbl>3;
}
	| GAUSSIAN_RELEASE_NUMBER '{'
	MEAN_NUMBER '=' num_expr
	'}'
{
  mdlpvp->rsop->release_number_method=CONSTNUM;
  mdlpvp->rsop->release_number=(int) $<dbl>5;
};

gaussian_release_number_cmd: GAUSSIAN_RELEASE_NUMBER '{'
	MEAN_NUMBER '=' num_expr
        STANDARD_DEVIATION '=' num_expr
	'}'
{
  mdlpvp->rsop->release_number_method=GAUSSNUM;
  mdlpvp->rsop->mean_number=(int) $<dbl>5;
  mdlpvp->rsop->standard_deviation=$<dbl>8;
};

volume_dependent_number_cmd: VOLUME_DEPENDENT_RELEASE_NUMBER '{'
	MEAN_DIAMETER '=' num_expr
	STANDARD_DEVIATION '=' num_expr
	CONCENTRATION '=' num_expr
	'}'
{
  mdlpvp->rsop->release_number_method=VOLNUM;
  mdlpvp->rsop->mean_diameter=$<dbl>5;
  mdlpvp->rsop->standard_deviation=$<dbl>8;
  mdlpvp->rsop->concentration=$<dbl>11;
};


point: array_value
{
  mdlpvp->el_head=$<nel>1;
  if ((mdlpvp->pntp1=(struct vector3 *)malloc(sizeof(struct vector3)))==NULL) {
    mdlerror("Cannot store point data");
    return(1);
  }
  mdlpvp->pntp1->x=0;
  mdlpvp->pntp1->y=0;
  mdlpvp->pntp1->z=0;
  mdlpvp->elp=mdlpvp->el_head;
  if (mdlpvp->elp!=NULL) {
    mdlpvp->pntp1->x=mdlpvp->elp->value;
    mdlpvp->elp=mdlpvp->elp->next;
  }
  else {
    mdlerror("Three dimensional value required");
    return(1);
  }
  if (mdlpvp->elp!=NULL) {
    mdlpvp->pntp1->y=mdlpvp->elp->value;
    mdlpvp->elp=mdlpvp->elp->next;
  }
  else {
    mdlerror("Three dimensional value required");
    return(1);
  }
  if (mdlpvp->elp!=NULL) {
    mdlpvp->pntp1->z=mdlpvp->elp->value;
    mdlpvp->elp=mdlpvp->elp->next;
  }
  else {
    mdlerror("Three dimensional value required");
    return(1);
  }
  $$=mdlpvp->pntp1;
};


existing_molecule: VAR
{
  if (mdlpvp->cval_2!=NULL) {  
    mdlpvp->sym_name=mdlpvp->cval_2;
  }   
  else {
    mdlpvp->sym_name=mdlpvp->cval;
  } 

  if ((mdlpvp->gp=retrieve_sym(mdlpvp->sym_name,MOL,volp->main_sym_table))==NULL) {
    sprintf(mdlpvp->mdl_err_msg,"%s %s","Undefined molecule:",mdlpvp->sym_name);
    mdlerror(mdlpvp->mdl_err_msg,mdlpvp);
    if (mdlpvp->sym_name==mdlpvp->cval) {
      mdlpvp->cval=NULL;
    }
    else {
      mdlpvp->cval_2=NULL;
    }
    free((void *)mdlpvp->sym_name);
    return(1);
  }
  if (mdlpvp->sym_name==mdlpvp->cval) {
    mdlpvp->cval=NULL;
  }
  else {
    mdlpvp->cval_2=NULL;
  }
  free((void *)mdlpvp->sym_name);
#ifdef KELP
  mdlpvp->gp->ref_count++;
  no_printf("ref_count: %d\n",mdlpvp->gp->ref_count);
#endif
  $$=mdlpvp->gp;
};


instance_def: INSTANTIATE new_object OBJECT '{'
{
  mdlpvp->gp=$<sym>2;
  mdlpvp->objp=(struct object *)mdlpvp->gp->value;

  mdlpvp->curr_obj=volp->root_instance;
  if (mdlpvp->curr_obj->first_child==NULL) { 
    mdlpvp->curr_obj->first_child=mdlpvp->objp;
  }
  if (mdlpvp->curr_obj->last_child!=NULL) {
    mdlpvp->curr_obj->last_child->next=mdlpvp->objp;
  } 
  mdlpvp->curr_obj->last_child=mdlpvp->objp;
  mdlpvp->objp->parent=mdlpvp->curr_obj;
  mdlpvp->objp->next=NULL;

  mdlpvp->objp->object_type=META_OBJ;
  mdlpvp->objp->parent=mdlpvp->curr_obj;
  mdlpvp->curr_obj=mdlpvp->objp;
}
        list_objects
        list_opt_object_cmds
        '}'
{
  mdlpvp->curr_obj=volp->root_object;
  if (mdlpvp->object_name_list_end->prev!=NULL) {
    mdlpvp->object_name_list_end=mdlpvp->object_name_list_end->prev;
  }
  else {
    mdlpvp->object_name_list_end->name=NULL;
  }
};


existing_object: VAR
{
  if (mdlpvp->prefix_name!=NULL) {
    free((void *)mdlpvp->prefix_name);
    mdlpvp->prefix_name=NULL;
  }
  if (mdlpvp->cval_2!=NULL) {
    mdlpvp->sym_name=mdlpvp->cval_2;
  }
  else {
    mdlpvp->sym_name=mdlpvp->cval;
  }
  if ((mdlpvp->gp=retrieve_sym(get_first_name(mdlpvp->sym_name),OBJ,volp->main_sym_table))==NULL) {
    sprintf(mdlpvp->mdl_err_msg,"%s %s","Undefined object:",mdlpvp->sym_name);
    mdlerror(mdlpvp->mdl_err_msg,mdlpvp);
    if (mdlpvp->sym_name==mdlpvp->cval) {
      mdlpvp->cval=NULL;
    }
    else {
      mdlpvp->cval_2=NULL;
    }
    free((void *)mdlpvp->sym_name);
    return(1);
  }
  mdlpvp->top_objp=(struct object *)mdlpvp->gp->value;
  if ((mdlpvp->objp=find_full_name(mdlpvp->top_objp,mdlpvp->sym_name,NULL))==NULL) {
    sprintf(mdlpvp->mdl_err_msg,"%s %s","Undefined object:",mdlpvp->sym_name);
    mdlerror(mdlpvp->mdl_err_msg,mdlpvp);
    if (mdlpvp->sym_name==mdlpvp->cval) {
      mdlpvp->cval=NULL;
    }
    else {
      mdlpvp->cval_2=NULL;
    }
    free((void *)mdlpvp->sym_name);
    return(1);
  }
  if (mdlpvp->sym_name==mdlpvp->cval) {
    mdlpvp->cval=NULL;
  }
  else {
    mdlpvp->cval_2=NULL;
  }
  sprintf(mdlpvp->full_name,"%s",mdlpvp->sym_name);
  mdlpvp->prefix_name=get_prefix_name(mdlpvp->sym_name);
#ifdef KELP
  mdlpvp->objp->sym->ref_count++;
  no_printf("ref_count: %d\n",mdlpvp->objp->sym->ref_count);
#endif
  $$=mdlpvp->objp->sym;
};


rx_net_def: DEFINE_REACTIONS '{'
	list_rx_stmts
	'}'
;


list_rx_stmts: rx_stmt
	| list_rx_stmts rx_stmt
;


rx_stmt: rx_group_def
	| rxn
;


rx_group_def: REACTION_GROUP reaction_group_name '{'
	list_rxns
        '}'
;


reaction_group_name: VAR
{
};


list_rxns: rxn
	| list_rxns rxn
;


rxn: unimolecular_rxn
	| bimolecular_rxn
;


unimolecular_rxn: reactant RT_ARROW 
{
  mdlpvp->gp=$<sym>1;
  if ((mdlpvp->tp=retrieve_sym(mdlpvp->gp->name,RX,volp->main_sym_table))
      !=NULL) {
  }
  else if ((mdlpvp->tp=store_sym(mdlpvp->gp->name,RX,volp->main_sym_table))
      ==NULL) {
    sprintf(mdlpvp->mdl_err_msg,"%s %s %s","Cannot store reaction:",
      mdlpvp->gp->name," -> ... ");
    mdlerror(mdlpvp->mdl_err_msg,mdlpvp);
    return(1);
  }
  /* fill in reactant parts of struct rxn here */
  if (mdlpvp->path_mem==NULL) {
    if ((mdlpvp->path_mem=create_mem(sizeof(struct pathway),16384))==NULL) {
      sprintf(mdlpvp->mdl_err_msg,"%s %s %s","Cannot store reaction:",
        mdlpvp->gp->name," -> ... ");
      mdlerror(mdlpvp->mdl_err_msg,mdlpvp);
      return(1);
    } 
  }
  if (mdlpvp->prod_mem==NULL) {
    if ((mdlpvp->prod_mem=create_mem(sizeof(struct product),16384))==NULL) {
      sprintf(mdlpvp->mdl_err_msg,"%s %s %s","Cannot store reaction:",
        mdlpvp->gp->name," -> ... ");
      mdlerror(mdlpvp->mdl_err_msg,mdlpvp);
      return(1);
    } 
  }
  if ((mdlpvp->pathp=(struct pathway *)mem_get(mdlpvp->path_mem))==NULL) {
    sprintf(mdlpvp->mdl_err_msg,"%s %s %s","Cannot store reaction:",
      mdlpvp->gp->name," -> ... ");
    mdlerror(mdlpvp->mdl_err_msg,mdlpvp);
    return(1);
  }
  mdlpvp->rxnp=(struct rxn *)mdlpvp->tp->value;
  mdlpvp->rxnp->n_reactants=1;
  mdlpvp->rxnp->n_pathways++;
  mdlpvp->pathp->reactant1=(struct species *)mdlpvp->gp->value;
  mdlpvp->pathp->reactant2=NULL;
  mdlpvp->pathp->reactant3=NULL;
  mdlpvp->pathp->km=0;
  mdlpvp->pathp->kcat=0;
  mdlpvp->pathp->orientation1=mdlpvp->orient_class;
  mdlpvp->pathp->orientation2=0;
  mdlpvp->pathp->orientation3=0;
  mdlpvp->pathp->product_head=NULL;

  mdlpvp->pathp->next=mdlpvp->rxnp->pathway_head;
  mdlpvp->rxnp->pathway_head=mdlpvp->pathp;
  
}
	list_products fwd_rx_rate1
{
  mdlpvp->pathp->km=mdlpvp->fwd_km;
#ifdef DEBUG
  no_printf("Unimolecular reaction defined:\n");
  no_printf("  %s[%d] ->",mdlpvp->rxnp->pathway_head->reactant1->sym->name,
    mdlpvp->rxnp->pathway_head->orientation1);
  for (mdlpvp->prodp=mdlpvp->rxnp->pathway_head->product_head;
      mdlpvp->prodp!=NULL;mdlpvp->prodp=mdlpvp->prodp->next) {
    if (mdlpvp->prodp!=mdlpvp->rxnp->pathway_head->product_head) {
      no_printf(" +");
    }
    no_printf(" %s[%d]",mdlpvp->prodp->prod->sym->name,mdlpvp->prodp->orientation);
  }
  no_printf(" [%.9g]\n",mdlpvp->rxnp->pathway_head->km);
#endif
};


bimolecular_rxn: reactant '+'
{
  mdlpvp->orient_class1=mdlpvp->orient_class;
}
	reactant RT_ARROW 
{
  mdlpvp->stp1=$<sym>1;
  mdlpvp->stp2=$<sym>4;
  mdlpvp->orient_class2=mdlpvp->orient_class;
  mdlpvp->sym_name=concat_rx_name(mdlpvp->stp1->name,mdlpvp->stp2->name);
  if ((mdlpvp->stp3=retrieve_sym(mdlpvp->sym_name,RX,volp->main_sym_table))
      !=NULL) {
  }
  else if ((mdlpvp->stp3=store_sym(mdlpvp->sym_name,RX,volp->main_sym_table))
      ==NULL) {
    sprintf(mdlpvp->mdl_err_msg,"%s %s + %s -> ...","Cannot store reaction:",
      mdlpvp->stp1->name,mdlpvp->stp2->name);
    mdlerror(mdlpvp->mdl_err_msg,mdlpvp);
    return(1);
  }
  /* fill in reactant parts of struct rxn here */
  if (mdlpvp->path_mem==NULL) {
    if ((mdlpvp->path_mem=create_mem(sizeof(struct pathway),16384))==NULL) {
      sprintf(mdlpvp->mdl_err_msg,"%s %s + %s -> ...","Cannot store reaction:",
        mdlpvp->stp1->name,mdlpvp->stp2->name);
      mdlerror(mdlpvp->mdl_err_msg,mdlpvp);
      return(1);
    } 
  }
  if (mdlpvp->prod_mem==NULL) {
    if ((mdlpvp->prod_mem=create_mem(sizeof(struct product),16384))==NULL) {
      sprintf(mdlpvp->mdl_err_msg,"%s %s + %s -> ...","Cannot store reaction:",
        mdlpvp->stp1->name,mdlpvp->stp2->name);
      mdlerror(mdlpvp->mdl_err_msg,mdlpvp);
      return(1);
    } 
  }
  if ((mdlpvp->pathp=(struct pathway *)mem_get(mdlpvp->path_mem))==NULL) {
    sprintf(mdlpvp->mdl_err_msg,"%s %s + %s -> ...","Cannot store reaction:",
      mdlpvp->stp1->name,mdlpvp->stp2->name);
    mdlerror(mdlpvp->mdl_err_msg,mdlpvp);
    return(1);
  }
  mdlpvp->rxnp=(struct rxn *)mdlpvp->stp3->value;
  mdlpvp->rxnp->n_reactants=2;
  mdlpvp->rxnp->n_pathways++;
  mdlpvp->pathp->reactant1=(struct species *)mdlpvp->stp1->value;
  mdlpvp->pathp->reactant2=(struct species *)mdlpvp->stp2->value;
  mdlpvp->pathp->reactant3=NULL;
  mdlpvp->pathp->km=0;
  mdlpvp->pathp->kcat=0;
  mdlpvp->pathp->orientation1=mdlpvp->orient_class1;
  mdlpvp->pathp->orientation2=mdlpvp->orient_class2;
  mdlpvp->pathp->orientation3=0;
  mdlpvp->pathp->product_head=NULL;

  mdlpvp->pathp->next=mdlpvp->rxnp->pathway_head;
  mdlpvp->rxnp->pathway_head=mdlpvp->pathp;
  
}
	list_products fwd_rx_rate1
{
  mdlpvp->pathp->km=mdlpvp->fwd_km;
#ifdef DEBUG
  no_printf("Bimolecular reaction defined:\n");
  no_printf("  %s[%d] + %s[%d] ->",
    mdlpvp->rxnp->pathway_head->reactant1->sym->name,
    mdlpvp->rxnp->pathway_head->orientation1,
    mdlpvp->rxnp->pathway_head->reactant2->sym->name,
    mdlpvp->rxnp->pathway_head->orientation2);
  for (mdlpvp->prodp=mdlpvp->rxnp->pathway_head->product_head;
      mdlpvp->prodp!=NULL;mdlpvp->prodp=mdlpvp->prodp->next) {
    if (mdlpvp->prodp!=mdlpvp->rxnp->pathway_head->product_head) {
      no_printf(" +");
    }
    no_printf(" %s[%d]",mdlpvp->prodp->prod->sym->name,mdlpvp->prodp->orientation);
  }
  no_printf(" [%.9g]\n",mdlpvp->rxnp->pathway_head->km);
#endif
};


reactant: existing_molecule
{
  mdlpvp->orient_class=1;
}
	orientation_class
{
  $$=$<sym>1;
};


list_products: product
	| list_products '+' product
;


product: existing_molecule
{
  mdlpvp->orient_class=1;
}
	orientation_class
{
  mdlpvp->gp=$<sym>1;
  if ((mdlpvp->prodp=(struct product *)mem_get(mdlpvp->prod_mem))==NULL) {
    sprintf(mdlpvp->mdl_err_msg,"%s %s %s","Cannot store reaction:",
      mdlpvp->rxnp->sym->name," -> ... ");
    mdlerror(mdlpvp->mdl_err_msg,mdlpvp);
    return(1);
  }
  mdlpvp->prodp->prod=(struct species *)mdlpvp->gp->value;
  mdlpvp->prodp->orientation=mdlpvp->orient_class;

  mdlpvp->prodp->next=mdlpvp->pathp->product_head;
  mdlpvp->pathp->product_head=mdlpvp->prodp;
}
	| NO_SPECIES
{
};


orientation_class: /* empty */
	| list_orient_marks
;

list_orient_marks: list_head_marks
	| list_tail_marks
{
  mdlpvp->orient_class*=-1;
};


list_head_marks: head_mark
	| list_head_marks head_mark
;


head_mark: '\''
{
  mdlpvp->orient_class+=1;
};


list_tail_marks: tail_mark
	| list_tail_marks tail_mark
;


tail_mark: ','
{
  mdlpvp->orient_class+=1;
};


fwd_rx_rate1: '[' num_expr ']'
{
  mdlpvp->fwd_km=$<dbl>2;
}
	| '[' '>' num_expr ']'
{
  mdlpvp->fwd_km=$<dbl>3;
};


fwd_rx_rate2: '[' num_expr ',' num_expr ']'
{
  mdlpvp->fwd_km=$<dbl>2;
  mdlpvp->fwd_kcat=$<dbl>4;
}
	| '[' '>' num_expr ',' num_expr ']'
{
  mdlpvp->fwd_km=$<dbl>3;
  mdlpvp->fwd_kcat=$<dbl>5;
};





%%





/* **************************************************************** */

#if 0



viz_output_def: VIZ_DATA_OUTPUT '{'
	list_viz_output_cmds
	'}'
;


list_viz_output_cmds: viz_output_cmd
	| list_viz_output_cmds viz_output_cmd
;


viz_output_cmd:
	viz_mode_def
	| voxel_image_mode_def
	| voxel_volume_mode_def
	| viz_output_list_def
	| viz_frame_data_def
	| viz_molecule_prefix_def
	| viz_object_prefixes_def
	| viz_state_values_def
;

viz_mode_def: MODE '=' NONE
{
	volp->viz_mode = NO_VIZ_MODE; 
}
	| MODE '=' DX
{
	volp->viz_mode = DX_MODE; 
}
	| MODE '=' RADIANCE
{
	volp->viz_mode = RADIANCE_MODE;
}
	| MODE '=' RAYSHADE
{
	volp->viz_mode = RAYSHADE_MODE;
}
	| MODE '=' POVRAY
{
	volp->viz_mode = POVRAY_MODE;
}
	| MODE '=' RENDERMAN
{
	volp->viz_mode = RENDERMAN_MODE;
}
	| MODE '=' IRIT
{
	volp->viz_mode = IRIT_MODE;
}
	| MODE '=' MCELL_GENERIC
{
	volp->viz_mode = MCELL_MODE;
};

voxel_image_mode_def: VOXEL_IMAGE_MODE '=' boolean
{
  volp->voxel_image_mode = $<tok>3;
};

voxel_volume_mode_def: VOXEL_VOLUME_MODE '=' boolean
{
  volp->voxel_volume_mode = $<tok>3;
};

viz_output_list_def: viz_iteration_def | viz_time_def ;

viz_iteration_def: ITERATION_LIST '='
{
  mdlpvp->num_pos=0;
  mdlpvp->el_head=NULL;
  mdlpvp->el_tail=NULL;
}
	'[' list_range_specs ']'
{
  sort_num_expr_list(mdlpvp->el_head);
  if ((mdlpvp->fdlp=(struct frame_data_list *)malloc
        (sizeof(struct frame_data_list)))==NULL) {
    mdlerror("Cannot store iteration list data");
    return(1);
  }
  mdlpvp->fdlp->list_type=FRAME_NUMBER;
  mdlpvp->fdlp->type=ALL_FRAME_DATA;
  mdlpvp->fdlp->viz_iteration=-1;
  mdlpvp->fdlp->n_viz_iterations=0;
  mdlpvp->fdlp->iteration_list=mdlpvp->el_head;
  mdlpvp->fdlp->curr_viz_iteration=mdlpvp->el_head;
  mdlpvp->fdlp->next=volp->frame_data_head;
  volp->frame_data_head = mdlpvp->fdlp;  
};

list_range_specs: range_spec
	| list_range_specs ',' range_spec
;

range_spec: num_expr
{
  if ((mdlpvp->elp=(struct num_expr_list *)malloc
      (sizeof(struct num_expr_list)))==NULL) {
    mdlerror("Cannot store iteration list data");
    return(1);
  }
 
  mdlpvp->elp->value=$<dbl>1;
  if (mdlpvp->el_tail==NULL) {
    mdlpvp->el_tail=mdlpvp->elp;
  }
  mdlpvp->el_tail->next=mdlpvp->elp;
  mdlpvp->elp->next=NULL;
  mdlpvp->el_tail=mdlpvp->elp;
  if (mdlpvp->el_head==NULL) {
    mdlpvp->el_head=mdlpvp->elp;
  }
  mdlpvp->num_pos++;
}
	| '[' num_expr TO num_expr STEP num_expr ']'
{
  for (mdlpvp->tmp_dbl=$<dbl>2;mdlpvp->tmp_dbl<$<dbl>4;mdlpvp->tmp_dbl+=$<dbl>6) {
    if ((mdlpvp->elp=(struct num_expr_list *)malloc
        (sizeof(struct num_expr_list)))==NULL) {
      mdlerror("Cannot store iteration list data");
      return(1);
    }
  
    mdlpvp->elp->value=mdlpvp->tmp_dbl;
    if (mdlpvp->el_tail==NULL) {
      mdlpvp->el_tail=mdlpvp->elp;
    }
    mdlpvp->el_tail->next=mdlpvp->elp;
    mdlpvp->elp->next=NULL;
    mdlpvp->el_tail=mdlpvp->elp;
    if (mdlpvp->el_head==NULL) {
      mdlpvp->el_head=mdlpvp->elp;
    }
    mdlpvp->num_pos++;
  }
  /*
   ** Introduce EPSILON=1e-14 to reduce the computer fluctuation to 
   ** the double numbers, which will miss the last data of the range 
   ** expression.
   */
  if (fabs($<dbl>4-mdlpvp->tmp_dbl)<=EPSILON) {
    if ((mdlpvp->elp=(struct num_expr_list *)malloc
        (sizeof(struct num_expr_list)))==NULL) {
      mdlerror("Cannot store iteration list data");
      return(1);
    }
  
    mdlpvp->elp->value=mdlpvp->tmp_dbl;
    if (mdlpvp->el_tail==NULL) {
      mdlpvp->el_tail=mdlpvp->elp;
    }
    mdlpvp->el_tail->next=mdlpvp->elp;
    mdlpvp->elp->next=NULL;
    mdlpvp->el_tail=mdlpvp->elp;
    if (mdlpvp->el_head==NULL) {
      mdlpvp->el_head=mdlpvp->elp;
    }
    mdlpvp->num_pos++;
  }

};


viz_time_def: TIME_LIST '='
{
  mdlpvp->num_pos=0;
  mdlpvp->el_head=NULL;
  mdlpvp->el_tail=NULL;
}
	'[' list_range_specs ']'
{
  sort_num_expr_list(mdlpvp->el_head);
  if ((mdlpvp->fdlp=(struct frame_data_list *)malloc
        (sizeof(struct frame_data_list)))==NULL) {
    mdlerror("Cannot store iteration list data");
    return(1);
  }
  mdlpvp->fdlp->list_type=REAL_TIME;
  mdlpvp->fdlp->type=ALL_FRAME_DATA;
  mdlpvp->fdlp->viz_iteration=-1;
  mdlpvp->fdlp->n_viz_iterations=0;
  mdlpvp->fdlp->iteration_list=mdlpvp->el_head;
  mdlpvp->fdlp->curr_viz_iteration=mdlpvp->el_head;
  mdlpvp->fdlp->next=volp->frame_data_head;
  volp->frame_data_head = mdlpvp->fdlp;  
};


viz_frame_data_def: ITERATION_FRAME_DATA '{'
	list_frame_data_specs
	'}'
;

list_frame_data_specs: frame_data_spec
	| list_frame_data_specs frame_data_spec
;

frame_data_spec: frame_data_item '='
{
  mdlpvp->num_pos=0;
  mdlpvp->el_head=NULL;
  mdlpvp->el_tail=NULL;
}
	'[' list_range_specs ']'
{
  sort_num_expr_list(mdlpvp->el_head);
  if ((mdlpvp->fdlp=(struct frame_data_list *)malloc
        (sizeof(struct frame_data_list)))==NULL) {
    mdlerror("Cannot store iteration frame data");
    return(1);
  }
  mdlpvp->fdlp->list_type=FRAME_NUMBER;
  mdlpvp->fdlp->type=$<tok>1;
  mdlpvp->fdlp->viz_iteration=-1;
  mdlpvp->fdlp->n_viz_iterations=0;
  mdlpvp->fdlp->iteration_list=mdlpvp->el_head;
  mdlpvp->fdlp->curr_viz_iteration=mdlpvp->el_head;
  mdlpvp->fdlp->next=volp->frame_data_head;
  volp->frame_data_head = mdlpvp->fdlp;
};


frame_data_item:  ALL
{
  $$=ALL_FRAME_DATA;
}
	| EFFECTOR_POSITIONS
{
  $$=EFF_POS;
}
	| EFFECTOR_STATES
{
  $$=EFF_STATES;
}
	| MOLECULE_POSITIONS
{
  $$=MOL_POS;
}
	| MOLECULE_STATES
{
  $$=MOL_STATES;
}
	| SURFACE_POSITIONS
{
  $$=SURF_POS;
}
	| SURFACE_STATES
{
  $$=SURF_STATES;
};


viz_molecule_prefix_def: MOLECULE_FILE_PREFIX '=' str_expr
{
  volp->molecule_prefix_name = $<str>3;
};

viz_object_prefixes_def: OBJECT_FILE_DESIGNATORS '{'
	list_viz_object_prefixes
	'}'
;

list_viz_object_prefixes: viz_object_prefix
	| list_viz_object_prefixes viz_object_prefix
;

viz_object_prefix: existing_object '=' str_expr
{
  mdlpvp->gp=$<sym>1;
  mdlpvp->objp=(struct object *)mdlpvp->gp->value;
  if ((mdlpvp->vizp=(struct viz_obj *)malloc(sizeof(struct viz_obj)))==NULL) {
    mdlerror("Cannot store viz object");
    return(1);
  }
  mdlpvp->objp->viz_obj = mdlpvp->vizp;
  mdlpvp->vizp->name = $<str>3;
  mdlpvp->vizp->full_name = my_strdup(mdlpvp->full_name);
  mdlpvp->vizp->obj = mdlpvp->objp;
  mdlpvp->vizp->cmprt_data_list = NULL;
  mdlpvp->vizp->next = volp->viz_obj_head;
  volp->viz_obj_head = vizp;
};

viz_state_values_def: STATE_VALUES '{'
	list_viz_state_values
	'}'
;

list_viz_state_values: viz_state_value
	| list_viz_state_values viz_state_value
;

viz_state_value: existing_logicalOrPhysical '=' num_expr
{
  mdlpvp->gp=$<sym>1;
  mdlpvp->viz_state=(int) $<dbl>3;
  switch (mdlpvp->gp->sym_type) {
    case OBJ:
      mdlpvp->existing_state=0;  
      slp=surf_state_head;
      while (slp!=NULL && !existing_state) {
        existing_state=existing_state||(viz_state==slp->state);
        slp=slp->next;
      }
      if (!existing_state) {
        if ((slp=(struct state_list *)malloc 
            (sizeof(struct state_list)))==NULL) {
          mdlerror("MCell: error cannot store state list");
          return(1);
        }
        slp->state=viz_state;
        slp->name=mdlpvp->gp->name;
        slp->next=surf_state_head;
        surf_state_head=slp;
      }
      mdlpvp->objp=(struct object *)mdlpvp->gp->value;
      switch (mdlpvp->objp->object_type) {
        case META_OBJ:
          if (set_viz_state_value(mdlpvp->objp,viz_state)) {
            mdlerror("Cannot store viz state value for object");
            return(1);
          }
          break;
        case BOX_OBJ:
          pop=(struct polygon_object *)mdlpvp->objp->obj;
          if (mdlpvp->objp->viz_state==NULL) {
            if ((mdlpvp->objp->viz_state=(int *)malloc
                 (pop->n_polys*sizeof(int)))==NULL) {
              mdlerror("Cannot store viz state value for object");
              return(1);
            }
          }
          for (i=0;i<pop->n_polys;i++) {
            mdlpvp->objp->viz_state[i]=viz_state;
          }
          break;
        case POLY_OBJ:
          pop=(struct polygon_object *)mdlpvp->objp->obj;
          if (mdlpvp->objp->viz_state==NULL) {
            if ((mdlpvp->objp->viz_state=(int *)malloc
                 (pop->n_polys*sizeof(int)))==NULL) {
              mdlerror("Cannot store viz state value for object");
              return(1);
            }
          }
          for (i=0;i<pop->n_polys;i++) {
            mdlpvp->objp->viz_state[i]=viz_state;
          }
          break;
        default:
          mdlerror("Cannot set viz state value of this type of object");
          return(1);
          break;
      }
      break;
    case RX:
      existing_state=0;  
      slp=eff_state_head;
      while (slp!=NULL && !existing_state) {
        existing_state=existing_state||(viz_state==slp->state);
        slp=slp->next;
      }
      if (!existing_state) {
        if ((slp=(struct state_list *)malloc 
            (sizeof(struct state_list)))==NULL) {
          mdlerror("MCell: error cannot store state list");
          return(1);
        }
        slp->state=viz_state;
        slp->name=mdlpvp->gp->name;
        slp->next=eff_state_head;
        eff_state_head=slp;
      }
      rxp=(struct rx *)mdlpvp->gp->value;
      rxp->viz_state=viz_state;
      break;
    case MOL:
      existing_state=0;  
      slp=lig_state_head;
      while (slp!=NULL && !existing_state) {
        existing_state=existing_state||(viz_state==slp->state);
        slp=slp->next;
      }
      if (!existing_state) {
        if ((slp=(struct state_list *)malloc 
            (sizeof(struct state_list)))==NULL) {
          mdlerror("MCell: error cannot store state list");
          return(1);
        }
        slp->state=viz_state;
        slp->name=mdlpvp->gp->name;
        slp->next=lig_state_head;
        lig_state_head=slp;
      }
      ligip=(struct ligand_info *)mdlpvp->gp->value;
      ligip->viz_state=viz_state;
      break;
  }
}
	| existing_logicalOrPhysical
{
  element_head=NULL;
}
	'[' list_element_specs ']' '=' num_expr
{
  mdlpvp->gp=$<sym>1;
  viz_state=(int) $<dbl>7;
  switch (mdlpvp->gp->sym_type) {
    case OBJ:
      existing_state=0;  
      slp=surf_state_head;
      while (slp!=NULL && !existing_state) {
        existing_state=existing_state||(viz_state==slp->state);
        slp=slp->next;
      }
      if (!existing_state) {
        if ((slp=(struct state_list *)malloc 
            (sizeof(struct state_list)))==NULL) {
          mdlerror("MCell: error cannot store state list");
          return(1);
        }
        slp->state=viz_state;
        slp->name=mdlpvp->gp->name;
        slp->next=surf_state_head;
        surf_state_head=slp;
      }
      mdlpvp->objp=(struct object *)mdlpvp->gp->value;
      switch (mdlpvp->objp->object_type) {
        case META_OBJ:
          mdlerror("Cannot set viz state value of elements of this type of object");
          return(1);
          break;
        case BOX_OBJ:
          pop=(struct polygon_object *)mdlpvp->objp->obj;
          if (mdlpvp->objp->viz_state==NULL) {
            if ((mdlpvp->objp->viz_state=(int *)malloc
                 (pop->n_polys*sizeof(int)))==NULL) {
              mdlerror("Cannot store viz state value for object");
              return(1);
            }
            for (i=0;i<pop->n_polys;i++) {
              mdlpvp->objp->viz_state[i]=EXCLUDE_OBJ;
            }
          }
          elmlp=element_head;
          while (elmlp!=NULL) {
            if (elmlp->begin==ALL_SIDES) {
              elmlp->begin=0;
              elmlp->end=pop->n_polys-1;
            }
            if (elmlp->begin < 0 || elmlp->end > pop->n_polys-1) {
              mdlerror("Cannot set viz state value -- element out of range");
              return(1);
            }
            for (i=elmlp->begin;i<=elmlp->end;i++) {
              mdlpvp->objp->viz_state[i]=viz_state;
            }
	    elmlp=elmlp->next;
          }
          break;
        case POLY_OBJ:
          pop=(struct polygon_object *)mdlpvp->objp->obj;
          if (mdlpvp->objp->viz_state==NULL) {
            if ((mdlpvp->objp->viz_state=(int *)malloc
                 (pop->n_polys*sizeof(int)))==NULL) {
              mdlerror("Cannot store viz state value for object");
              return(1);
            }
            for (i=0;i<pop->n_polys;i++) {
              mdlpvp->objp->viz_state[i]=EXCLUDE_OBJ;
            }
          }
          elmlp=element_head;
          while (elmlp!=NULL) {
            if (elmlp->begin==ALL_SIDES) {
              elmlp->begin=0;
              elmlp->end=pop->n_polys-1;
            }
            if (elmlp->begin < 0 || elmlp->end > pop->n_polys-1) {
              mdlerror("Cannot set viz state value -- element out of range");
              return(1);
            }
            for (i=elmlp->begin;i<=elmlp->end;i++) {
              mdlpvp->objp->viz_state[i]=viz_state;
            }
	    elmlp=elmlp->next;
          }
          break;
        default:
          mdlerror("Cannot set viz state value of this type of object");
          return(1);
          break;
      }
      break;
    default:
      mdlerror("Cannot set viz state value of elements of this type of object");
      return(1);
      break;
  }
};

existing_logicalOrPhysical: VAR
{
  if (mdlpvp->cval_2!=NULL) {
    mdlpvp->sym_name=mdlpvp->cval_2;
  }
  else {
    mdlpvp->sym_name=mdlpvp->cval;
  }
  if ((mdlpvp->gp=retrieve_sym(mdlpvp->sym_name,RX,volp->main_sym_table))==NULL) {
    if ((mdlpvp->gp=retrieve_sym(mdlpvp->sym_name,MOL,volp->main_sym_table))==NULL) {
      if ((mdlpvp->gp=retrieve_sym(get_first_name(mdlpvp->sym_name),OBJ,volp->main_sym_table))==NULL) {
        sprintf(mdlpvp->mdl_err_msg,"%s %s","Undefined object:",mdlpvp->sym_name);
        mdlerror(mdlpvp->mdl_err_msg,mdlpvp);
        if (mdlpvp->sym_name==mdlpvp->cval) {
          mdlpvp->cval=NULL;
        }
        else {
          mdlpvp->cval_2=NULL;
        }
        free((void *)mdlpvp->sym_name);
        return(1);
      }
      mdlpvp->top_objp=(struct object *)mdlpvp->gp->value;
      if ((mdlpvp->objp=find_full_name(mdlpvp->top_objp,mdlpvp->sym_name,NULL))==NULL) {
        sprintf(mdlpvp->mdl_err_msg,"%s %s","Undefined object:",mdlpvp->sym_name);
        mdlerror(mdlpvp->mdl_err_msg,mdlpvp);
        if (mdlpvp->sym_name==mdlpvp->cval) {
          mdlpvp->cval=NULL;
        }
        else {
          mdlpvp->cval_2=NULL;
        }
        free((void *)mdlpvp->sym_name);
        return(1);
      }
      mdlpvp->gp=mdlpvp->objp->sym;
    }
  }
  if (mdlpvp->sym_name==mdlpvp->cval) {
    mdlpvp->cval=NULL;
  }
  else {
    mdlpvp->cval_2=NULL;
  }
  free((void *)mdlpvp->sym_name);
#ifdef KELP
  mdlpvp->gp->ref_count++;
  no_printf("ref_count: %d\n",mdlpvp->gp->ref_count);
#endif
  $$=mdlpvp->gp;
};

list_element_specs: element_spec
	| list_element_specs ',' element_spec
;

element_spec: num_expr
{
  if ((elmlp=(struct element_list *)malloc
             (sizeof(struct element_list)))==NULL) {
  }
  elmlp->begin=(unsigned int) ($<dbl>1);
  elmlp->end=elmlp->begin;
  elmlp->next=element_head;
  element_head=elmlp;
}
	| num_expr TO num_expr
{
  if ((elmlp=(struct element_list *)malloc
             (sizeof(struct element_list)))==NULL) {
  }
  elmlp->begin=(unsigned int) ($<dbl>1);
  elmlp->end=(unsigned int) ($<dbl>3);
  elmlp->next=element_head;
  element_head=elmlp;
}
	| side_name
{
  if ((elmlp=(struct element_list *)malloc
             (sizeof(struct element_list)))==NULL) {
  }
  elmlp->begin=$<tok>1;
  elmlp->end=elmlp->begin;
  elmlp->next=element_head;
  element_head=elmlp;
};

parallel_partition_def: PARALLEL_PARTITION '=' partition_plane
;

partition_plane: X_TOK
{
  parallel_x = 1;
}
	| XY_TOK
{
  parallel_x = 1;
  parallel_y = 1;
}
	| XZ_TOK
{
  parallel_x = 1;
  parallel_z = 1;
}
	| XYZ_TOK
{
  parallel_x = 1;
  parallel_y = 1;
  parallel_z = 1;
}
	| Y_TOK
{
  parallel_y = 1;
}
	| YZ_TOK
{
  parallel_y = 1;
  parallel_z = 1;
}
	| Z_TOK
{
  parallel_z = 1;
};

partition_def: partition_dimension '='
{
  mdlpvp->num_pos=0;
  mdlpvp->el_head=NULL;
  mdlpvp->el_tail=NULL;
}
        '[' list_range_specs ']' 
{
  if ((mdlpvp->dblp=(double *)malloc((mdlpvp->num_pos+2)*sizeof(double)))==NULL) {
    mdlerror("Cannot store volume partition data");
    return(1);
  }
  sort_num_expr_list(mdlpvp->el_head);
  i=1;
  mdlpvp->elp=mdlpvp->el_head;
  while(mdlpvp->elp!=NULL) {
    mdlpvp->dblp[i++]=mdlpvp->elp->value/volp->length_unit;
    mdlpvp->elp=mdlpvp->elp->next;
  }
  mdlpvp->dblp[0]=-vol_infinity;
  mdlpvp->dblp[mdlpvp->num_pos+1]=vol_infinity;
  switch ($<tok>1) {
  case X_PART:
    if (volume->x_partitions!=NULL) {
      free(volume->x_partitions);
    }
    volume->n_x_partitions=mdlpvp->num_pos;
    volume->x_partitions=mdlpvp->dblp;
    break;
  case Y_PART:
    if (volume->y_partitions!=NULL) {
      free(volume->y_partitions);
    }
    volume->n_y_partitions=mdlpvp->num_pos;
    volume->y_partitions=mdlpvp->dblp;
    break;
  case Z_PART:
    if (volume->z_partitions!=NULL) {
      free(volume->z_partitions);
    }
    volume->n_z_partitions=mdlpvp->num_pos;
    volume->z_partitions=mdlpvp->dblp;
    break;
  }
};

partition_dimension: PARTITION_X {$$=X_PART;} 
	| PARTITION_Y {$$=Y_PART;}
	| PARTITION_Z {$$=Z_PART;}
;





surface_region_def: DEFINE_SURFACE_REGIONS '{'
	list_object_region_def
	'}'
{
};

list_object_region_def: object_region_def
	| list_object_region_def object_region_def
;

object_region_def: OBJECT existing_object '{'
{
  mdlpvp->gp=$<sym>2;
  mdlpvp->objp=(struct object *)mdlpvp->gp->value;
  obj_name=mdlpvp->objp->sym->name;
  if (mdlpvp->objp->object_type!=BOX_OBJ && mdlpvp->objp->object_type!=POLY_OBJ) {
    sprintf(mdlpvp->mdl_err_msg,"%s %s","Cannot define region on non-surface object:",mdlpvp->objp->sym->name);
    mdlerror(mdlpvp->mdl_err_msg,mdlpvp);
    return(1);
  }
  region_list_head=mdlpvp->objp->region_list;
}
	list_region_def
	'}'
{
  mdlpvp->objp->region_list=region_list_head;
};

list_region_def: region_def
	| list_region_def region_def
;

region_def: REGION new_region '{'
{
  element_head=NULL;
}
	ELEMENT_LIST '=' '[' list_element_specs ']'
	'}'
{
  mdlpvp->gp=$<sym>2;
  rp=(struct region *)mdlpvp->gp->value;
  pop=(struct polygon_object *)mdlpvp->objp->obj;
  elmlp=element_head;
  while (elmlp!=NULL) {
    if (elmlp->begin==ALL_SIDES) {
      elmlp->begin=0;
      elmlp->end=pop->n_polys-1;
    }
    if (elmlp->begin < 0 || elmlp->end > pop->n_polys-1) {
    mdlerror("Cannot create region -- element out of range");
    return(1);
    }
    elmlp=elmlp->next;
  }
  mdlpvp->objp->num_regions++;
  rp->element_list=element_head;
};

new_region: VAR
{
  if (mdlpvp->cval_2!=NULL) {
    mdlpvp->sym_name=mdlpvp->cval_2;
  }
  else {
    mdlpvp->sym_name=mdlpvp->cval;
  }
  if ((rp=make_new_region(obj_name,mdlpvp->sym_name,mdlpvp->mdl_err_msg))==NULL) {
    mdlerror(mdlpvp->mdl_err_msg,mdlpvp);
    return(1);
  }
  if ((rlp=(struct region_list *)malloc(sizeof(struct region_list)))==NULL) {
    sprintf(mdlpvp->mdl_err_msg,"%s %s","Cannot store region name:",rp->sym->name);
    mdlerror(mdlpvp->mdl_err_msg,mdlpvp);
    return(1);
  }
  if (mdlpvp->sym_name==mdlpvp->cval) {
    mdlpvp->cval=NULL;
  }
  else {
    mdlpvp->cval_2=NULL;
  }
  rp->region_last_name=mdlpvp->sym_name;
  rp->parent=mdlpvp->curr_obj;
  rp->reg_counter_ref_list=NULL;
  rp->lig_hit_counter=NULL;
  rlp->region=rp;
  rlp->next=region_list_head;
  region_list_head=rlp;
  no_printf("Creating new region: %s\n",rp->sym->name);
  fflush(stderr);
  $$=rp->sym;
};


effector_site_def: DEFINE_EFFECTOR_SITE_POSITIONS '{'
	list_region_ref
	'}'
{
};

list_region_ref: region_ref
	| list_region_ref region_ref
;

region_ref: REGION existing_region '{'
{
  mdlpvp->gp=$<sym>2;
  rp=(struct region *)mdlpvp->gp->value;
}
	list_add_effector_state
	'}'
;

list_add_effector_state: add_effector_state
	| list_add_effector_state add_effector_state
;

add_effector_state: EFFECTOR_STATE existing_reaction_state '{'
	effector_quantity_cmd
	orientation_cmd
	'}'
{
  mdlpvp->tp=$<sym>2;
  rxp=(struct rx *)mdlpvp->tp->value;
  if ((effdp=(struct eff_dat *)malloc(sizeof(struct eff_dat)))==NULL) {
    sprintf(mdlpvp->mdl_err_msg,"%s %s","Cannot store object name:",mdlpvp->sym_name);
    return(1);
  }
  effdp->next=rp->eff_dat;
  rp->eff_dat=effdp;
  effdp->rx=rxp;
  effdp->quantity_type=effector_quantity_type;
  effdp->quantity=$<dbl>4;
  effdp->orient=(signed char) $<tok>5;
};

effector_quantity_cmd: DENSITY '=' num_expr
{
  $$=$<dbl>3;
  effector_quantity_type=EFFDENS;
}
	| NUMBER '=' num_expr
{
  $$=$<dbl>3;
  effector_quantity_type=EFFNUM;
};

existing_region: existing_object '[' VAR ']'
{
  if (mdlpvp->cval_2!=NULL) {
    mdlpvp->sym_name=mdlpvp->cval_2;
  }
  else {
    mdlpvp->sym_name=mdlpvp->cval;
  }
  mdlpvp->tp=$<sym>1;
  mdlpvp->objp=(struct object *)mdlpvp->tp->value;
  obj_name=mdlpvp->objp->sym->name;
  strncpy(temp,"",256);
  strncpy(temp,obj_name,254);
  strcat(temp,",");   
  region_name=my_strcat(temp,mdlpvp->sym_name);
  if ((mdlpvp->gp=retrieve_sym(region_name,REG,volp->main_sym_table))==NULL) {
    sprintf(mdlpvp->mdl_err_msg,"%s %s","Undefined region:",region_name);
    mdlerror(mdlpvp->mdl_err_msg,mdlpvp);
    if (mdlpvp->sym_name==mdlpvp->cval) {
      mdlpvp->cval=NULL;
    }
    else {
      mdlpvp->cval_2=NULL;
    }
    free((void *)mdlpvp->sym_name);
    return(1);
  }
  free((void *)region_name);
  if (mdlpvp->sym_name==mdlpvp->cval) {
    mdlpvp->cval=NULL;
  }
  else {
    mdlpvp->cval_2=NULL;
  }
  free((void *)mdlpvp->sym_name);
#ifdef KELP
  mdlpvp->gp->ref_count++;
  no_printf("ref_count: %d\n",mdlpvp->gp->ref_count);
#endif
  $$=mdlpvp->gp;
};








polygon_def: new_object POLYGON_LIST '{'
{
  mdlpvp->gp=$<sym>1;
  mdlpvp->objp=(struct object *)mdlpvp->gp->value;
  if ((pop=(struct polygon_object *)malloc
              (sizeof(struct polygon_object)))==NULL) {
    mdlerror("Cannot store box object data");
    return(1);
  }
  pop->lig_count_ref=NULL;
  pop->viz_state_ref=NULL;
  pop->list_type=ORDERED_POLY;
  pop->polygon_data=NULL;
  pop->n_polys=0;
  pop->fully_closed=0;
  pop->side_stat=NULL;
  pop->cmprt_side_map=NULL;
  pop->lig_prop=NULL;

  if ((opp=(struct ordered_poly *)malloc
              (sizeof(struct ordered_poly)))==NULL) {
    mdlerror("Cannot store polygon object data");
    return(1);
  }
  opp->vertex=NULL;
  opp->normal=NULL;
  opp->element_data=NULL;
  opp->n_verts=0;
  pop->polygon_data=(void *)opp;

  mdlpvp->objp->object_type=POLY_OBJ;
  mdlpvp->objp->contents=pop;
  mdlpvp->objp->parent=mdlpvp->curr_obj;
  mdlpvp->curr_obj=mdlpvp->objp;
}
	vertex_list_cmd
        element_connection_cmd
{
  pop->n_polys=n_polys;
  if ((pop->side_stat=(unsigned short *)malloc
              (pop->n_polys*sizeof(unsigned short)))==NULL) {
    mdlerror("Cannot store polygon object data");
    return(1);
  }
  if ((pop->cmprt_side_map=(unsigned int *)malloc
              (pop->n_polys*sizeof(unsigned int)))==NULL) {
    mdlerror("Cannot store polygon object data");
    return(1);
  }
  if ((pop->lig_prop=(byte **)malloc
              (pop->n_polys*sizeof(byte *)))==NULL) {
    mdlerror("Cannot store polygon object data");
    return(1);
  }
  for (i=0;i<pop->n_polys;i++) {
    if ((pop->lig_prop[i]=(byte *)malloc
         ((1+n_ligand_types)*sizeof(byte)))==NULL) {
      mdlerror("Cannot store polygon object data");
      return(1);
    }
    for (j=0;j<1+n_ligand_types;j++) {
      *(pop->lig_prop[i]+j)=RFLCT;
    }
    pop->side_stat[i]=1;
    pop->cmprt_side_map[i]=0;
  }

  opp->n_verts=n_verts;
  if ((opp->vertex=(struct vector3 **)malloc
              (opp->n_verts*sizeof(struct vector3 *)))==NULL) {
    mdlerror("Cannot store polygon object data");
    return(1);
  }

  vlp=vertex_head;
  if (vlp->normal!=NULL) {
    if ((opp->normal=(struct vector3 **)malloc
                (opp->n_verts*sizeof(struct vector3 *)))==NULL) {
      mdlerror("Cannot store polygon object data");
      return(1);
    }
  }
  for (i=0;i<n_verts;i++) {
    opp->vertex[i]=vlp->vertex;
    if (opp->normal!=NULL) {
      opp->normal[i]=vlp->normal;
    }
    vlp=vlp->next;
  }
  if ((edp=(struct element_data *)malloc
              (pop->n_polys*sizeof(struct element_data)))==NULL) {
    mdlerror("Cannot store polygon object data");
    return(1);
  }
  opp->element_data=edp;
  eclp=connection_head;
  for (i=0;i<pop->n_polys;i++) {
    edp[i].n_verts=eclp->n_verts;
    if ((intp=(int *)malloc(eclp->n_verts*sizeof(int)))==NULL) {
      mdlerror("Cannot store polygon object data");
      return(1);
    }
    edp[i].vertex_index=intp;
    mdlpvp->elp=eclp->connection_list;
    for (j=0;j<eclp->n_verts;j++) {
      edp[i].vertex_index[j]=(int)mdlpvp->elp->value;
      mdlpvp->elp=mdlpvp->elp->next;
    }
    eclp=eclp->next;
  }
}
        FULLY_CLOSED '=' boolean
	list_polygon_object_cmds
	list_opt_object_cmds
	'}'
{
  no_printf("Polygon list %s defined:\n",mdlpvp->curr_obj->sym->name);
  pop->fully_closed=0;
  mdlpvp->curr_obj=mdlpvp->curr_obj->parent;
  if (object_name_list_end->prev!=NULL) {
    object_name_list_end=object_name_list_end->prev;
  }
  else {
    object_name_list_end->name=NULL;
  }
  $$=$<sym>1;
};

vertex_list_cmd: VERTEX_LIST '{' 
{
  n_verts=0;
  vertex_head=NULL;
  vertex_tail=NULL;
}
	list_points 
        '}'
;

list_points: point
{
  pntp1=$<vec3>1;
  if ((vlp=(struct vertex_list *)malloc(sizeof(struct vertex_list)))==NULL) {
    mdlerror("Cannot store vertex data");
    return(1);
  }
  n_verts++;
  vlp->vertex=pntp1;
  vlp->normal=NULL;
  if (vertex_tail==NULL) {
    vertex_tail=vlp;
  }
  vertex_tail->next=vlp;
  vlp->next=NULL;
  vertex_tail=vlp;
  if (vertex_head==NULL) {
    vertex_head=vlp;
  }
}
        | point NORMAL point
{
  pntp1=$<vec3>1;
  pntp2=$<vec3>3;
  if ((vlp=(struct vertex_list *)malloc(sizeof(struct vertex_list)))==NULL) {
    mdlerror("Cannot store vertex data");
    return(1);
  }
  n_verts++;
  vlp->vertex=pntp1;
  vlp->normal=pntp2;
  if (vertex_tail==NULL) {
    vertex_tail=vlp;
  }
  vertex_tail->next=vlp;
  vlp->next=NULL;
  vertex_tail=vlp;
  if (vertex_head==NULL) {
    vertex_head=vlp;
  }
}
	| list_points point
{
  pntp1=$<vec3>2;
  if ((vlp=(struct vertex_list *)malloc(sizeof(struct vertex_list)))==NULL) {
    mdlerror("Cannot store vertex data");
    return(1);
  }
  n_verts++;
  vlp->vertex=pntp1;
  vlp->normal=NULL;
  if (vertex_tail==NULL) {
    vertex_tail=vlp;
  }
  vertex_tail->next=vlp;
  vlp->next=NULL;
  vertex_tail=vlp;
  if (vertex_head==NULL) {
    vertex_head=vlp;
  }
}
	| list_points point NORMAL point
{
  pntp1=$<vec3>2;
  pntp2=$<vec3>4;
  if ((vlp=(struct vertex_list *)malloc(sizeof(struct vertex_list)))==NULL) {
    mdlerror("Cannot store vertex data");
    return(1);
  }
  n_verts++;
  vlp->vertex=pntp1;
  vlp->normal=pntp2;
  if (vertex_tail==NULL) {
    vertex_tail=vlp;
  }
  vertex_tail->next=vlp;
  vlp->next=NULL;
  vertex_tail=vlp;
  if (vertex_head==NULL) {
    vertex_head=vlp;
  }
};

element_connection_cmd: ELEMENT_CONNECTIONS '{'
{
  n_polys=0;
  connection_head=NULL;
  connection_tail=NULL;
}
        list_arrays
        '}'
;

list_arrays: array_value
{
  mdlpvp->elp=$<nel>1;
  if ((eclp=(struct element_connection_list *)malloc
              (sizeof(struct element_connection_list)))==NULL) {
    mdlerror("Cannot store element connection data");
    return(1);
  }
  n_polys++;
  eclp->connection_list=mdlpvp->elp;
  eclp->n_verts=0;
  while (mdlpvp->elp!=NULL) {
    eclp->n_verts++;
    mdlpvp->elp=mdlpvp->elp->next;
  }
  if (connection_tail==NULL) {
    connection_tail=eclp;
  }
  connection_tail->next=eclp;
  eclp->next=NULL;
  connection_tail=eclp;
  if (connection_head==NULL) {
    connection_head=eclp;
  }
}
        | list_arrays array_value
{
  mdlpvp->elp=$<nel>2;
  if ((eclp=(struct element_connection_list *)malloc
              (sizeof(struct element_connection_list)))==NULL) {
    mdlerror("Cannot store element connection data");
    return(1);
  }
  n_polys++;
  eclp->connection_list=mdlpvp->elp;
  eclp->n_verts=0;
  while (mdlpvp->elp!=NULL) {
    eclp->n_verts++;
    mdlpvp->elp=mdlpvp->elp->next;
  }
  if (connection_tail==NULL) {
    connection_tail=eclp;
  }
  connection_tail->next=eclp;
  eclp->next=NULL;
  connection_tail=eclp;
  if (connection_head==NULL) {
    connection_head=eclp;
  }
};

boolean: TRUE { $$=1; }
	| FALSE { $$=0; }
	| YES { $$=1; }
	| NO { $$=0; }
;

box_def: new_object BOX '{'
{
  mdlpvp->gp=$<sym>1;
  mdlpvp->objp=(struct object *)mdlpvp->gp->value;
  if ((pop=(struct polygon_object *)malloc
              (sizeof(struct polygon_object)))==NULL) {
    mdlerror("Cannot store box object data");
    return(1);
  }
  pop->lig_count_ref=NULL;
  pop->viz_state_ref=NULL;
  pop->list_type=BOX_POLY;
  pop->polygon_data=NULL;
  pop->n_polys=0;
  pop->fully_closed=0;
  pop->side_stat=NULL;
  pop->cmprt_side_map=NULL;
  pop->lig_prop=NULL;

  if ((bpp=(struct box_poly *)malloc
              (sizeof(struct box_poly)))==NULL) {
    mdlerror("Cannot store box object data");
    return(1);
  }
  bpp->llf=NULL;
  bpp->urb=NULL;
  pop->polygon_data=(void *)bpp;

  pop->n_polys=6;
  if ((pop->side_stat=(unsigned short *)malloc
              (pop->n_polys*sizeof(unsigned short)))==NULL) {
    mdlerror("Cannot store box object data");
    return(1);
  }
  if ((pop->cmprt_side_map=(unsigned int *)malloc
              (pop->n_polys*sizeof(unsigned int)))==NULL) {
    mdlerror("Cannot store box object data");
    return(1);
  }
  if ((pop->lig_prop=(byte **)malloc
              (pop->n_polys*sizeof(byte *)))==NULL) {
    mdlerror("Cannot store box object data");
    return(1);
  }
  for (i=0;i<pop->n_polys;i++) {
    if ((pop->lig_prop[i]=(byte *)malloc
         ((1+n_ligand_types)*sizeof(byte)))==NULL) {
      mdlerror("Cannot store box object data");
      return(1);
    }
    for (j=0;j<1+n_ligand_types;j++) {
      *(pop->lig_prop[i]+j)=RFLCT;
    }
    pop->side_stat[i]=1;
    pop->cmprt_side_map[i]=0;
  }

  mdlpvp->objp->object_type=BOX_OBJ;
  mdlpvp->objp->contents=pop;
  mdlpvp->objp->parent=mdlpvp->curr_obj;
  mdlpvp->curr_obj=mdlpvp->objp;
}
	CORNERS '=' point ',' point
        FULLY_CLOSED '=' boolean
	list_polygon_object_cmds
	list_opt_object_cmds
	'}'
{
  bpp->llf=$<vec3>7;
  bpp->urb=$<vec3>9;

  pop->fully_closed=$<tok>12;

#ifdef DEBUG
  no_printf("Box %s defined:\n",mdlpvp->curr_obj->sym->name);
  pntp1=bpp->llf;
  no_printf("    LLF Corner = [ %f, %f, %f ]\n",pntp1->x,pntp1->y,pntp1->z);
  pntp1=bpp->urb;
  no_printf("    URB Corner = [ %f, %f, %f ]\n",pntp1->x,pntp1->y,pntp1->z);
#endif
  mdlpvp->curr_obj=mdlpvp->curr_obj->parent;
  if (object_name_list_end->prev!=NULL) {
    object_name_list_end=object_name_list_end->prev;
  }
  else {
    object_name_list_end->name=NULL;
  }
  $$=$<sym>1;
};

list_polygon_object_cmds: polygon_object_cmd
	| list_polygon_object_cmds polygon_object_cmd
	| /* empty */
{
};

polygon_object_cmd:
	remove_side
	| wall_prop_cmd
	| add_effector
;

remove_side: REMOVE_ELEMENT '=' side_name
{
  if ($<tok>3==ALL_SIDES) {
    for (i=0;i<pop->n_polys;i++) {
      pop->side_stat[i]=0;
    }
    no_printf("Element removed: ALL\n");
  }
  else {
    if ($<tok>3 > (pop->n_polys)-1) {
      mdlerror("Reference to non-existant element");
      return(1);
    }
    else {
      pop->side_stat[$<tok>3]=0;
      no_printf("Element removed: %d\n",$<tok>3);
    }
  }
  fflush(stderr);
}
	| REMOVE_ELEMENT '=' num_expr
{
  if ($<dbl>3 > (pop->n_polys)-1) {
    mdlerror("Reference to non-existant element");
    return(1);
  }
  else {
    pop->side_stat[(int)$<dbl>3]=0;
    no_printf("Element removed: %d\n",$<dbl>3);
  }
};

side: ELEMENT '=' side_name
{
	$$=$<tok>3;
}
	| ELEMENT '=' num_expr
{
  $$=(int) $<dbl>3;
};

side_name: TOP {$$=TP;}
	| BOTTOM {$$=BOT;}
	| FRONT {$$=FRNT;}
	| BACK {$$=BCK;}
	| LEFT {$$=LFT;}
	| RIGHT {$$=RT;}
	| ALL_ELEMENTS {$$=ALL_SIDES;}
;

wall_prop_cmd: wall_prop '{'
	MOLECULE '=' existing_molecule
	side
{
  mdlpvp->gp=$<sym>5;
  ligip=(struct ligand_info *)mdlpvp->gp->value;
  if ($<tok>6==ALL_SIDES) {
    for (i=0;i<pop->n_polys;i++) {
      pop->lig_prop[i][ligip->type]=$<tok>1;
    }
  }
  else {
    if ($<tok>6 > (pop->n_polys)-1) {
      mdlerror("Reference to non-existant element");
      return(1);
    }
    else {
      pop->lig_prop[$<tok>6][ligip->type]=$<tok>1;
    }
  }
  no_printf("Property of element %d for molecule %s = %d\n",
          $<tok>6,mdlpvp->gp->name,$<tok>1);
  fflush(stderr);
}
	'}'
; 

wall_prop: REFLECTIVE {$$=RFLCT;}
	| TRANSPARENT {$$=TRANSP;}
	| ABSORPTIVE {$$=SINK;}
;

add_effector: ADD_EFFECTOR '{'
	req_add_effector_cmds
	'}'
;

req_add_effector_cmds:  STATE '=' existing_reaction_state
        DENSITY '=' num_expr
        side
	orientation_cmd
{
        mdlpvp->tp=$<sym>3;
        rxp=(struct rx *)mdlpvp->tp->value;
	effdp=NULL;

        /* create automatic region to hold effector data*/
        sprintf(mdlpvp->str_buf,"%d",mdlpvp->curr_obj->num_regions++);
        mdlpvp->sym_name=my_strdup(mdlpvp->str_buf);
        if ((rp=make_new_region(obj_name,mdlpvp->sym_name,mdlpvp->mdl_err_msg))==NULL) {
          mdlerror(mdlpvp->mdl_err_msg,mdlpvp);
          return(1);
        }
        if ((rlp=(struct region_list *)malloc(sizeof(struct region_list)))==NULL) {
          sprintf(mdlpvp->mdl_err_msg,"%s %s","Cannot store region name:",rp->sym->name);
          mdlerror(mdlpvp->mdl_err_msg,mdlpvp);
          return(1);
        }
        rp->region_last_name=mdlpvp->sym_name;
        rp->parent=mdlpvp->curr_obj;
	rp->reg_counter_ref_list=NULL;
	rp->lig_hit_counter=NULL;
        rlp->region=rp;
        rlp->next=mdlpvp->curr_obj->region_list;
        mdlpvp->curr_obj->region_list=rlp;
        no_printf("Creating new region: %s\n",rp->sym->name);
        fflush(stderr);

        if ((elmlp=(struct element_list *)malloc
                   (sizeof(struct element_list)))==NULL) {
        }
        elmlp->begin=0;
        elmlp->end=0;
        elmlp->next=NULL;
        rp->element_list=elmlp;

        if ($<tok>7==ALL_SIDES) {
          if ((effdp=(struct eff_dat *)malloc
	     (sizeof(struct eff_dat)))==NULL) {
            mdlerror("Cannot store effector data for object");
            return(1);
          }
          effdp->rx=rxp;
          effdp->quantity_type=EFFDENS;
          effdp->quantity=$<dbl>6;
	  effdp->orient = (signed char) $<tok>8;
          effdp->next=NULL;
          rp->eff_dat=effdp;
          elmlp->begin=0;
          elmlp->end=pop->n_polys-1;
        }
        else {
          if ($<tok>7 > (pop->n_polys)-1) {
            mdlerror("Reference to non-existant element");
            return(1);
          }
          else {
	    if (pop->side_stat[$<tok>7]) {
              if ((effdp=(struct eff_dat *)malloc
	         (sizeof(struct eff_dat)))==NULL) {
                mdlerror("Cannot store effector data for object");
                return(1);
              }
              effdp->rx=rxp;
              effdp->quantity_type=EFFDENS;
              effdp->quantity=$<dbl>6;
	      effdp->orient = (signed char) $<tok>8;
              effdp->next=NULL;
              rp->eff_dat=effdp;
              elmlp->begin=$<tok>7;
              elmlp->end=elmlp->begin;
            }
            else {
              mdlerror("Cannot add effectors to removed element");
              return(1);
            }
	  }
        }
        no_printf("Add effector %s, density %f, side %d\n",
                mdlpvp->tp->name,$<dbl>6,$<tok>7);
        fflush(stderr);
};

orientation_cmd: /* empty */
{
	$$ = OUTWRD;
}
	| POLE_ORIENTATION '=' orientation
{
	$$ = $<tok>3;
};

orientation: POSITIVE_FRONT
{
	$$ = OUTWRD;
}
	| POSITIVE_BACK
{
	$$ = INWRD;
};

existing_reaction_state: VAR
{
  if (mdlpvp->cval_2!=NULL) {
    mdlpvp->sym_name=mdlpvp->cval_2;
  }       
  else {  
    mdlpvp->sym_name=mdlpvp->cval;  
  }       
  if ((mdlpvp->gp=retrieve_sym(mdlpvp->sym_name,RX,volp->main_sym_table))==NULL) {
    sprintf(mdlpvp->mdl_err_msg,"%s %s","Undefined reaction state:",mdlpvp->sym_name);
    mdlerror(mdlpvp->mdl_err_msg,mdlpvp);
    if (mdlpvp->sym_name==mdlpvp->cval) {
      mdlpvp->cval=NULL;
    }
    else {
      mdlpvp->cval_2=NULL;
    }
    free((void *)mdlpvp->sym_name);
    return(1);
  }
  if (mdlpvp->sym_name==mdlpvp->cval) {
    mdlpvp->cval=NULL;
  }
  else {
    mdlpvp->cval_2=NULL;
  }
  free((void *)mdlpvp->sym_name);
#ifdef KELP
  mdlpvp->gp->ref_count++;
  no_printf("ref_count: %d\n",mdlpvp->gp->ref_count);
#endif
  $$=mdlpvp->gp;
};


req_rx_cmds: list_sub_rxs
	reference_state_cmd
{
  rxp=ref_state;
  no_printf("Counting bound molecules for rx %s with %d states\n",rxp->sym->name,n_state);
  fflush(stderr);
  if ((prxp->rx_states=(struct rx **)malloc
       (n_state*sizeof(struct rx *)))==NULL) {
    mdlerror("Cannot store rx state data");
    return(1);
  }
  prxp->num_states=n_state;
  prxp->num_transitions=count_rx(ref_state,0);
  if (procnum == 0) fprintf(log_file,"\nMCell: Transition probabilities for reaction %s with %d states and %d transitions:\n\n",prxp->reaction_name,prxp->num_states,prxp->num_transitions);
  set_rx_rates(ref_state);
  if (procnum == 0) fprintf(log_file,"\n");
  no_printf("Total number of rx transitions = %d\n",prxp->num_transitions);
  no_printf("Done counting bound molecules\n");
  prxp->prev_state_flag=0;
};

reaction_state: VAR
{
  if (mdlpvp->cval_2!=NULL) {
    mdlpvp->sym_name=my_strcat(reaction_name,mdlpvp->cval_2);
    free(mdlpvp->cval_2);
    mdlpvp->cval_2=NULL;
  }
  else {
    mdlpvp->sym_name=my_strcat(reaction_name,mdlpvp->cval);
    free(mdlpvp->cval);
    mdlpvp->cval=NULL;
  }
  if ((mdlpvp->gp=store_sym(mdlpvp->sym_name,RX,volp->main_sym_table))==NULL) {
    sprintf(mdlpvp->mdl_err_msg,"%s %s","Cannot store reaction in table:",mdlpvp->sym_name);
    mdlerror(mdlpvp->mdl_err_msg,mdlpvp);
    free((void *)mdlpvp->sym_name);
    return(1);
  }
  rxp=(struct rx *)mdlpvp->gp->value;
  if (rxp->parent_rx==NULL) {
    rxp->parent_rx=prxp;
    rxp->hashval=hash(rxp->sym->name)&0x0000000f;
    rxp->state_index=n_state++;
  }
  $$=mdlpvp->gp;
};

reference_state_cmd: REFERENCE_STATE reaction_state '{'
{
  mdlpvp->gp=$<sym>2;
  rxp=(struct rx *)mdlpvp->gp->value;
  ref_state=rxp;
  no_printf("Setting bound molecules for reference state %s\n",rxp->sym->name);
}
        list_number_bound_cmds
	'}'
;

list_number_bound_cmds: number_bound_cmd
	| list_number_bound_cmds number_bound_cmd;

number_bound_cmd: existing_molecule NUMBER_BOUND '=' num_expr
{
  mdlpvp->gp=$<sym>1;
  ligip=(struct ligand_info *)mdlpvp->gp->value;
  ref_state->bound_ligands[ligip->type]=(unsigned char) $<dbl>4;
  no_printf("Number of molecules of type %s bound to reference state %s set to: %d\n",ligip->sym->name,ref_state->sym->name,ref_state->bound_ligands[ligip->type]);
};

list_sub_rxs: sub_rx | list_sub_rxs sub_rx;

sub_rx: reaction_state list_rx
{
	mdlpvp->gp=$<sym>1;
	rxp=(struct rx *)mdlpvp->gp->value;
	rxp->sym=mdlpvp->gp;
	rxp->hashval=hash(rxp->sym->name)&0x0000000f;
	rxp->parent_rx=prxp;
        no_printf("parser processing sub-rx: %s\n",rxp->sym->name);
        fflush(stderr);

	temp1=n_rxn;
        total_rate=0;
	for (i=0;i<temp1;i++) {
	  n=rxn[i].n_pathways;
	  for (j=0;j<n;j++) {
	    if ((rxn[i].event_type!=BIND)
                 && (rxn[i].event_type!=TRANSPORT)
                 && (rxn[i].event_type!=PRODUCE_POISSON)) {
	      total_rate+=rxn[i].pathway[j].rate;
            }
          }
        }
        pe=1.0-exp(-volp->time_unit*total_rate);
	for (i=0;i<temp1;i++) {
	  n=rxn[i].n_pathways;
	  for (j=0;j<n;j++) {
	    if ((rxn[i].event_type!=BIND)
                 && (rxn[i].event_type!=TRANSPORT)
                 && (rxn[i].event_type!=PRODUCE_POISSON)) {
              if (total_rate>0.0) {
	        pt=pe*(rxn[i].pathway[j].rate/total_rate);
              }
              else {
                pt=0;
              }
            }
          }
        }

	for (i=0;i<temp1;i++) {
	  n=rxn[i].n_pathways;
	  if ((rx_dat=(struct rx_info *)malloc
               (sizeof(struct rx_info)))==NULL) {
	    mdlerror("Cannot store rx_info");
	    return(1);
	  }
	  if ((rx_dat->rate=(double *)malloc
               (n*sizeof(double)))==NULL) {
	    mdlerror("Cannot store rate data");
	    return(1);
	  }
	  if ((rx_dat->rate_t=(struct t_func **)malloc
               (n*sizeof(struct t_func *)))==NULL) {
	    mdlerror("Cannot store rate data");
	    return(1);
	  }
	  if ((rx_dat->n_rate_t=(int *)malloc
	       (n*sizeof(int)))==NULL) {
	    mdlerror("Cannot store rate data");
	    return(1);
	  }
	  if ((rx_dat->rate_t_index=(int *)malloc
	       (n*sizeof(int)))==NULL) {
	    mdlerror("Cannot store rate data");
	    return(1);
	  }
	  if ((rx_dat->polarity=(signed char *)malloc
	       (n*sizeof(signed char)))==NULL) {
	    mdlerror("Cannot store polarity data");
	    return(1);
	  }
	  if ((rx_dat->transition_index=(int *)malloc
	       (n*sizeof(int)))==NULL) {
	    mdlerror("Cannot store transition index data");
	    return(1);
	  }
	  if ((rx_dat->next_state_index=(unsigned short *)malloc
	       (n*sizeof(unsigned short)))==NULL) {
	    mdlerror("Cannot store next state index data");
	    return(1);
	  }
	  if ((rx_dat->init_transition_count_dt=(int *)malloc
	       (n*sizeof(int)))==NULL) {
	    mdlerror("Cannot store transition count data");
	    return(1);
	  }
	  if ((rx_dat->init_transition_count_cum=(int *)malloc
	       (n*sizeof(int)))==NULL) {
	    mdlerror("Cannot store transition count data");
	    return(1);
	  }
	  if ((rx_dat->transition_count_dt=(int *)malloc
               (n*sizeof(int)))==NULL) {
	    mdlerror("Cannot store transition count data");
	    return(1);
	  }
	  if ((rx_dat->transition_count_cum=(int *)malloc
	       (n*sizeof(int)))==NULL) {
	    mdlerror("Cannot store transition count data");
	    return(1);
	  }
	  if ((rx_dat->next_rx=(struct rx **)malloc
               (n*sizeof(struct rx *)))==NULL) {
	    mdlerror("Cannot store next_rx data");
	    return(1);
	  }
	  switch (rxn[i].event_type) {
	    case BIND:
	      rx_type=&rxp->bind_rx[rxn[i].l_type];
	      rx_proc=bind_proc;
	      break;
	    case TRANSPORT:
	      rx_type=&rxp->transport_rx[rxn[i].l_type];
	      rx_proc=transport_proc;
	      break;
	    case DISSOC:
	      rx_type=&rxp->dissoc_rx[rxn[i].l_type];
	      rx_proc=dissoc_proc;
	      break;
	    case PRODUCE:
	      rx_type=&rxp->product_rx[rxn[i].l_type];
	      rx_proc=product_proc;
	      break;
	    case PRODUCE_POISSON:
	      rx_type=&rxp->product_poisson_rx[rxn[i].l_type];
	      rx_proc=product_poisson_proc;
	      break;
	    case DEGRADE:
	      rx_type=&rxp->degrade_rx[rxn[i].l_type];
	      rx_proc=degrade_proc;
	      break;
	    case ISOM:
	      rx_type=&(rxp->isom_rx);
	      rx_proc=isom_proc;
	      break;
	  }
	  *rx_type=rx_dat;
	  rx_dat->n_rates=n;
	  rx_dat->rx_proc=rx_proc;
	  for (j=0;j<n;j++) {
	    mdlpvp->gp=rxn[i].pathway[j].sym;
	    rx_dat->next_rx[j]=(struct rx *)mdlpvp->gp->value;
	    rx_dat->rate[j]=rxn[i].pathway[j].rate;
	    rx_dat->rate_t[j]=NULL;
	    rx_dat->n_rate_t[j]=0;
	    rx_dat->rate_t_index[j]=0;
            rate_file=rxn[i].pathway[j].rate_file;
            if (rate_file!=NULL) {
              if ((file=fopen(rate_file,"r"))==NULL) {
                sprintf(mdlpvp->mdl_err_msg,"error opening rate file: %s\n",rate_file);
                mdlerror(mdlpvp->mdl_err_msg,mdlpvp);
                return(1);
              }   
              count=0;
              while (fscanf(file,"%lf %lf",&mdlpvp->val_1,&mdlpvp->val_2)!=EOF) {
                count++;
              }
	      rx_dat->n_rate_t[j]=count;
              if ((rx_dat->rate_t[j]=(struct t_func *)malloc
                   (count*sizeof(struct t_func)))==NULL) {
                mdlerror("Cannot store subvolume partition data");
                return(1);
              }
              rewind(file);
              count=0;
              while (fscanf(file,"%lf %lf",&mdlpvp->val_1,&mdlpvp->val_2)!=EOF) {
	        (rx_dat->rate_t[j])[count].time=mdlpvp->val_1;
	        if (rxn[i].event_type==BIND || rxn[i].event_type==TRANSPORT) {
	          (rx_dat->rate_t[j])[count].value=mdlpvp->val_2*rxn[i].pathway[j].rate_file_scale_factor;
                }
                else {
	          (rx_dat->rate_t[j])[count].value=mdlpvp->val_2;
                }
                count++;
              }
              fclose(file);
            }
	    rx_dat->polarity[j]=rxn[i].pathway[j].polarity;
	    rx_dat->transition_index[j]=0;
	    rx_dat->next_state_index[j]=0;
	    rx_dat->init_transition_count_dt[j]=0;
	    rx_dat->init_transition_count_cum[j]=0;
	    rx_dat->transition_count_dt[j]=0;
	    rx_dat->transition_count_cum[j]=0;
	  }
	}
        no_printf("parser done processing sub-rx: %s\n",rxp->sym->name);
        fflush(stderr);
	n_rxn=0;
};

list_rx: rx | list_rx rx;

rx: '[' transition '{' rate ':' event polarity '}' ']'
{

	mdlpvp->gp=$<sym>2;
        rx_match=0;
        for (i=0;i<n_rxn;i++) {
          if (rxn[i].l_type==ep.l_type && rxn[i].event_type==ep.event_type) {
            rx_match=1;
            rx_num=i;
          }
        }
        if (!rx_match) {
          rx_num=n_rxn;
          rxn[rx_num].n_pathways=0;
	  n_rxn++;
        }
	rxn[rx_num].l_type=ep.l_type;
	rxn[rx_num].event_type=ep.event_type;
        path_index=rxn[rx_num].n_pathways;
	rxn[rx_num].pathway[path_index].sym=mdlpvp->gp;
	rxn[rx_num].pathway[path_index].polarity=$<tok>7;
	if (ep.event_type==BIND) {
          if (rxn[rx_num].pathway[path_index].polarity==BOTH_P) {
	    rxn[rx_num].pathway[path_index].rate=0.5*mdlpvp->mc_factor*$<dbl>4/sqrt(ep.d_const);
	    rxn[rx_num].pathway[path_index].rate_file_scale_factor=0.5*mdlpvp->mc_factor/sqrt(ep.d_const);
	    rxn[rx_num].pathway[path_index].rate_file=rate_file;
          }
          else {
	    rxn[rx_num].pathway[path_index].rate=mdlpvp->mc_factor*$<dbl>4/sqrt(ep.d_const);
	    rxn[rx_num].pathway[path_index].rate_file_scale_factor=mdlpvp->mc_factor/sqrt(ep.d_const);
	    rxn[rx_num].pathway[path_index].rate_file=rate_file;
          }
	}
	else if (ep.event_type==TRANSPORT) {
	  rxn[rx_num].pathway[path_index].rate=mdlpvp->transport_mc_factor*$<dbl>4/(ep.charge*sqrt(ep.d_const));
	  rxn[rx_num].pathway[path_index].rate_file_scale_factor=mdlpvp->transport_mc_factor/(ep.charge*sqrt(ep.d_const));
	  rxn[rx_num].pathway[path_index].rate_file=rate_file;
          printf("Transport probability = %.16g\n",rxn[rx_num].pathway[path_index].rate);
	}
	else {
	  rxn[rx_num].pathway[path_index].rate=$<dbl>4;
	  rxn[rx_num].pathway[path_index].rate_file=rate_file;
	}
	rxn[rx_num].n_pathways++;
	no_printf("event: %d %d",ep.l_type,ep.event_type);
	no_printf(" transition: %s",mdlpvp->gp->name);
	no_printf(" rate: %f = %22.19f\n",
	            $<dbl>4,rxn[rx_num].pathway[path_index].rate);
        no_printf(" polarity: %d = %d\n",
                    $<tok>7,rxn[rx_num].pathway[path_index].polarity);
	fflush(stderr);

}
	| '[' transition '{' rate '}' ']'
{

	mdlpvp->gp=$<sym>2;
	ep.l_type=0;
	ep.event_type=ISOM;
        rx_match=0;
        for (i=0;i<n_rxn;i++) {
          if (rxn[i].l_type==ep.l_type && rxn[i].event_type==ep.event_type) {
            rx_match=1;
            rx_num=i;
          }
        }
        if (!rx_match) {
          rx_num=n_rxn;
          rxn[rx_num].n_pathways=0;
	  n_rxn++;
        }
	rxn[rx_num].l_type=ep.l_type;
	rxn[rx_num].event_type=ep.event_type;
        path_index=rxn[rx_num].n_pathways;
	rxn[rx_num].pathway[path_index].sym=mdlpvp->gp;
	rxn[rx_num].pathway[path_index].rate=$<dbl>4;
	rxn[rx_num].pathway[path_index].rate_file=rate_file;
	rxn[rx_num].pathway[path_index].polarity=DEFAULT_P;
	rxn[rx_num].n_pathways++;
	no_printf("event: %d %d",ep.l_type,ep.event_type);
	no_printf(" transition: %s",mdlpvp->gp->name);
	no_printf(" rate: %f = %22.19f\n",
	            $<dbl>4,rxn[rx_num].pathway[path_index].rate);
        no_printf(" polarity: %d\n",rxn[rx_num].pathway[path_index].polarity);
	fflush(stderr);
};

polarity: ',' polarity_spec { $$=$<tok>2; }
	| /* empty */ { $$=DEFAULT_P; }
;

polarity_spec: POSITIVE_POLE { $$=POS_P; }
	| NEGATIVE_POLE { $$=NEG_P; }
	| BOTH_POLES { $$=BOTH_P; }
	| BINDING_POLE { $$=BINDING_P; }
	| EITHER_POLE { $$=EITHER_P; }
;

event: '+' existing_molecule
{
	mdlpvp->gp=$<sym>2;
	ep.l_type=((struct ligand_info *)mdlpvp->gp->value)->type;
	ep.d_const=((struct ligand_info *)mdlpvp->gp->value)->ref_d_const;
	ep.event_type=BIND;
}
	| '~' existing_molecule
{
	mdlpvp->gp=$<sym>2;
	ep.l_type=((struct ligand_info *)mdlpvp->gp->value)->type;
	ep.d_const=((struct ligand_info *)mdlpvp->gp->value)->ref_d_const;
	ep.charge=((struct ligand_info *)mdlpvp->gp->value)->charge;
	ep.event_type=TRANSPORT;
}
	| '-' existing_molecule
{
	mdlpvp->gp=$<sym>2;
	ep.l_type=((struct ligand_info *)mdlpvp->gp->value)->type;
	ep.event_type=DISSOC;
}
	| '*' existing_molecule
{
	mdlpvp->gp=$<sym>2;
	ep.l_type=((struct ligand_info *)mdlpvp->gp->value)->type;
	ep.event_type=PRODUCE;
}
	| '#' existing_molecule
{
	mdlpvp->gp=$<sym>2;
	ep.l_type=((struct ligand_info *)mdlpvp->gp->value)->type;
	ep.event_type=DEGRADE;
}
	| '@' existing_molecule
{
	mdlpvp->gp=$<sym>2;
	ep.l_type=((struct ligand_info *)mdlpvp->gp->value)->type;
	ep.event_type=PRODUCE_POISSON;
};

transition: '>' reaction_state 
{
	$$=$<sym>2;
};

rate: num_expr_only 
{
  $$=$<dbl>1;
  rate_file=NULL;
}
	| str_expr_only
{
  $$=0;
  rate_file=$<str>1;
}
	| existing_var_only
{
  mdlpvp->gp=$<sym>1;
  rate_file=NULL;
  switch (mdlpvp->gp->sym_type) {
  case DBL:
    $$=*(double *)mdlpvp->gp->value;
    break;
  case STR:
    $$=0;
    rate_file=my_strdup((char *)mdlpvp->gp->value);
    break;
  default:
    mdlerror("Invalid variable type referenced");
    return(1);
    break;
  }
};  


output_def: REACTION_DATA_OUTPUT '{'
{
  build_ligand_table();
}
       output_step_def
       list_count_cmds
{  if ((olp=(struct output_list *)
       malloc(sizeof(struct output_list)))==NULL) {
    mdlerror("Cannot store output list data");
    return(1);
  }
  n_reac_frames++;
  olp->freq=(int) (output_freq+ROUND_UP);
  olp->reaction_list=rolp;
  olp->counter_info=cilp;
  olp->n_output=buffersize;
  olp->id=n_reac_frames-1;
  olp->counter=0;
  olp->index=0;
  olp->out_type=reac_out_type;
  olp->next=output_list;
  output_list=olp;
}
       
'}'
       ;


output_step_def:  step_def 
                | reac_iteration_def 
                | reac_time_def 
                ;

step_def: STEP '=' num_expr
{
  output_freq=$<dbl>3;
  output_freq=output_freq/volp->time_unit;
  if (chkpt_iterations) {
    n_output=(int)(chkpt_iterations/output_freq+1);
    buffersize=min_count_buffer(chkpt_iterations-start_time,n_output,COUNTBUFFERSIZE);
  }
  else {
    n_output=(int)(volp->iterations/output_freq+1);
    buffersize=min_count_buffer(volp->iterations-start_time,n_output,COUNTBUFFERSIZE);
  }
  
  reac_out_type=FREQ_DATA;
  if (output_freq>volp->iterations) {
    sprintf(mdlpvp->mdl_err_msg,"Output frequency too high\n\tSetting output frequency to %f microseconds\n",volp->iterations*volp->time_unit/1.0e-6);
    mdl_warning(mdlpvp);
    output_freq=volp->iterations;
  }
  if (output_freq<1) {
    sprintf(mdlpvp->mdl_err_msg,"Output frequency too low\n\tSetting output frequency to %f microseconds\n",volp->time_unit/1.0e-6);
    mdl_warning(mdlpvp);
    output_freq=1;
  }
  no_printf("Output definition:\n");
  no_printf("Output frequency = %f\n",output_freq);
}
        | /* empty */
{
  output_freq=1;
  if (output_freq>volp->iterations) {
    sprintf(mdlpvp->mdl_err_msg,"Output frequency too high\n\tSetting output frequency to %f microseconds\n",volp->iterations*volp->time_unit/1.0e-6);
    mdl_warning(mdlpvp);
    output_freq=volp->iterations;
  }
  if (output_freq<1) {
    sprintf(mdlpvp->mdl_err_msg,"Output frequency too low\n\tSetting output frequency to %f microseconds\n",volp->time_unit/1.0e-6);
    mdl_warning(mdlpvp);
    output_freq=1;
  }
  no_printf("Output definition:\n");
  no_printf("Output frequency = %f\n",output_freq);
};


reac_iteration_def: ITERATION_LIST '='
{
  mdlpvp->num_pos=0;
  mdlpvp->el_head=NULL;
  mdlpvp->el_tail=NULL;
}
	'[' list_range_specs ']'
{
  n_output=mdlpvp->num_pos;
  reac_out_type=FRAME_DATA;

  /**
   * Compute the output buffersize.
   **/
  if (chkpt_iterations) {
    buffersize=min_count_buffer(chkpt_iterations-start_time+1,n_output,COUNTBUFFERSIZE);
  }
  else {
    buffersize=min_count_buffer(volp->iterations-start_time+1,n_output,COUNTBUFFERSIZE);
  }

  sort_num_expr_list(mdlpvp->el_head);
  
  /**
   * Construct the reaction list. 
   **/
  if ((rolp=(struct reaction_list *)
       malloc(sizeof(struct reaction_list)))==NULL) {
    mdlerror("Cannot store output list data");
    return(1);
  }
  
  if ((intp=(int *)malloc(buffersize*sizeof(int)))==NULL) {
    mdlerror("Cannot store integer data");
    return(1);
  }

  i=0;
  /*
  mdlpvp->elp=mdlpvp->el_head;
  while(mdlpvp->elp!=NULL) {
     if (mdlpvp->elp->value>volp->iterations) {
      mdlerror("Reaction output iteration number is too big!");
      return(1);
    }
    intp[i++]=mdlpvp->elp->value;
    mdlpvp->elp=mdlpvp->elp->next;
  }
  */
  for (i=0;i<=buffersize;i++) {
    intp[i]=0;
  }

  
  rolp->list_type=FRAME_NUMBER;
  rolp->n_reac_iterations=0;
  rolp->reac_iteration=-1;
  rolp->array=intp;
  rolp->iteration_list=mdlpvp->el_head;
  rolp->curr_reac_iteration=mdlpvp->el_head;
  rolp->next=reaction_data_head;
  reaction_data_head=rolp;
}
;

reac_time_def: TIME_LIST '='
{
  mdlpvp->num_pos=0;
  mdlpvp->el_head=NULL;
  mdlpvp->el_tail=NULL;
}
	'[' list_range_specs ']'
{
  n_output=mdlpvp->num_pos;
  reac_out_type=FRAME_DATA;

  /**
   * Compute the output buffersize.
   **/
  if (chkpt_iterations) {
    buffersize=min_count_buffer(chkpt_iterations-start_time+1,n_output,COUNTBUFFERSIZE);
  }
  else {
    buffersize=min_count_buffer(volp->iterations-start_time+1,n_output,COUNTBUFFERSIZE);
  }

  sort_num_expr_list(mdlpvp->el_head);
  
  /**
   * Construct the reaction list. 
   **/
  if ((rolp=(struct reaction_list *)
       malloc(sizeof(struct reaction_list)))==NULL) {
    mdlerror("Cannot store output list data");
    return(1);
  }
  
  if ((intp=(int *)malloc(buffersize*sizeof(int)))==NULL) {
    mdlerror("Cannot store integer data");
    return(1);
  }

  i=0;

  for (i=0;i<=buffersize;i++) {
    intp[i]=0;
  }
  rolp->list_type=REAL_TIME;
  rolp->n_reac_iterations=0;
  rolp->reac_iteration=-1;
  rolp->array=intp;
  rolp->iteration_list=mdlpvp->el_head;
  rolp->curr_reac_iteration=mdlpvp->el_head;
  rolp->next=reaction_data_head;
  reaction_data_head=rolp;
};



list_count_cmds: count_cmd | list_count_cmds count_cmd;

count_cmd: '{' count_expr '}' '=' '>' outfile_syntax
{
  clp=$<cnt>2;
  cilp->count_list=clp;
  cilp->next=counter_info;
  counter_info=cilp;
  }
;

count_expr: count_value
{
	$$=$<cnt>1;
}
| '(' count_expr ')' 
{
  $$=$<cnt>2;
}
| count_expr '+' count_expr
{
  if ((clp=(struct count_list *)malloc
       (sizeof(struct count_list)))==NULL) {
    mdlerror("Cannot store counter data");
    return(1);
  }
  clp->freq=(int) (output_freq+ROUND_UP);
  clp->frame_index=n_reac_frames;
  clp->n_output=buffersize;
  clp->reset_flag=0;
  clp->update_flag=0;
  clp->data_type=EXPR;
  clp->index_type=UNKNOWN;
  clp->n_data=0;
  clp->temp_data=NULL;
  clp->final_data=NULL;
  clp->operand1=$<cnt>1;
  clp->operand2=$<cnt>3;
  clp->oper='+';
  clp->next=count_list;
  count_list=clp;
  $$=clp;
}
        | count_expr '-' count_expr
{
  if ((clp=(struct count_list *)malloc
       (sizeof(struct count_list)))==NULL) {
    mdlerror("Cannot store counter data");
    return(1);
  }
  clp->freq=(int) (output_freq+ROUND_UP);
  clp->frame_index=n_reac_frames;
  clp->n_output=buffersize;
  clp->reset_flag=0;
  clp->data_type=EXPR;
  clp->update_flag=0;
  clp->index_type=UNKNOWN;
  clp->n_data=0;
  clp->temp_data=NULL;
  clp->final_data=NULL;
  clp->operand1=$<cnt>1;
  clp->operand2=$<cnt>3;
  clp->oper='-';
  clp->next=count_list;
  count_list=clp;
  $$=clp;
}
        | count_expr '*' count_expr
{
  if ((clp=(struct count_list *)malloc
       (sizeof(struct count_list)))==NULL) {
    mdlerror("Cannot store counter data");
    return(1);
  }
  clp->freq=(int) (output_freq+ROUND_UP);
  clp->frame_index=n_reac_frames;
  clp->n_output=buffersize;
  clp->reset_flag=0;
  clp->data_type=EXPR;
  clp->update_flag=0;
  clp->index_type=UNKNOWN;
  clp->n_data=0;
  clp->temp_data=NULL;
  clp->final_data=NULL;
  clp->operand1=$<cnt>1;
  clp->operand2=$<cnt>3;
  clp->oper='*';
  clp->next=count_list;
  count_list=clp;
  $$=clp;
}
        | count_expr '/' count_expr
{
  if ((clp=(struct count_list *)malloc
       (sizeof(struct count_list)))==NULL) {
    mdlerror("Cannot store counter data");
    return(1);
  }
  clp->freq=(int) (output_freq+ROUND_UP);
  clp->frame_index=n_reac_frames;
  clp->n_output=buffersize;
  clp->reset_flag=0;
  clp->data_type=EXPR;
  clp->update_flag=0;
  clp->index_type=UNKNOWN;
  clp->n_data=0;
  clp->temp_data=NULL;
  clp->final_data=NULL;
  clp->operand1=$<cnt>1;
  clp->operand2=$<cnt>3;
  clp->oper='/';
  clp->next=count_list;
  count_list=clp;
  $$=clp;
}
        | '-' count_expr %prec UNARYMINUS
{
  if ((clp=(struct count_list *)malloc
       (sizeof(struct count_list)))==NULL) {
    mdlerror("Cannot store counter data");
    return(1);
  }
  if ((clp2=(struct count_list *)malloc
       (sizeof(struct count_list)))==NULL) {
    mdlerror("Cannot store counter data");
    return(1);
  }
  if (!(mdlpvp->dblp=(double *)malloc
	(sizeof(double)))) {
    mdlerror("Cannot store counter data");
    return(1);
  }
  *mdlpvp->dblp=-1;
  clp2->freq=(int) (output_freq+ROUND_UP);
  clp2->frame_index=n_reac_frames;
  clp2->n_output=buffersize;
  clp2->reset_flag=0;
  clp2->update_flag=0;
  clp2->index_type=UNKNOWN;
  clp2->n_data=1;
  clp2->data_type=DBL;
  clp2->temp_data=(void *)mdlpvp->dblp;
  clp2->final_data=(void *)mdlpvp->dblp;
  clp2->operand1=NULL;
  clp2->operand2=NULL;
  clp2->oper='\0';
  clp2->next=count_list;
  count_list=clp2;
  
  clp->freq=(int) (output_freq+ROUND_UP);
  clp->frame_index=n_reac_frames;
  clp->n_output=buffersize;
  clp->reset_flag=0;
  clp->update_flag=0;
  clp->index_type=UNKNOWN;
  clp->data_type=EXPR;
  clp->n_data=0;
  clp->temp_data=NULL;
  clp->final_data=NULL;
  clp->operand1=clp2;
  clp->operand2=$<cnt>2;
  clp->oper='*';
  clp->next=count_list;
  count_list=clp;
  $$=clp;
}
;

count_value: COUNT '[' count_value_init count_syntax ']'
{
        $$=$<cnt>3;
}
        | EXPRESSION '[' num_expr ']'
{
  if ((clp=(struct count_list *)malloc
       (sizeof(struct count_list)))==NULL) {
    mdlerror("Cannot store counter data");
    return(1);
  }
  if (!(mdlpvp->dblp=(double *)malloc
	(sizeof(double)))) {
    mdlerror("Cannot store counter data");
    return(1);
  }
  *mdlpvp->dblp=$<dbl>3;
  clp->freq=(int) (output_freq+ROUND_UP);
  clp->frame_index=n_reac_frames;
  clp->n_output=buffersize;
  clp->reset_flag=0;
  clp->update_flag=0;
  clp->index_type=TIME_STAMP_VAL;
  clp->n_data=1;
  clp->data_type=DBL;
  clp->temp_data=(void *)mdlpvp->dblp;
  clp->final_data=(void *)mdlpvp->dblp;
  clp->operand1=NULL;
  clp->operand2=NULL;
  clp->oper='\0';
  clp->next=count_list;
  count_list=clp;
  $$=clp;
}
; 

count_value_init: /* empty */
{
  if ((clp=(struct count_list *)malloc
       (sizeof(struct count_list)))==NULL) {
    mdlerror("Cannot store counter data");
    return(1);
  }
  clp->freq=(int) (output_freq+ROUND_UP);
  clp->frame_index=n_reac_frames;
  clp->n_output=buffersize;
  clp->reset_flag=0;
  clp->update_flag=1;
  clp->index_type=TIME_STAMP_VAL;
  clp->next=count_list;
  count_list=clp;
  $$=clp;
}
	;

outfile_syntax: file_name
{
  if ((cilp=(struct counter_info *)malloc
       (sizeof(struct counter_info)))==NULL) {
    mdlerror("Cannot store output list data");
    return(1);
  }

  sprintf(cilp->outfile_name,"%s",$<str>1);
  no_printf("Output file set to %s\n",cilp->outfile_name); 
};

count_syntax: r_transition_syntax
        | r_state_or_ligand_diffusion_syntax
	| lig_hit_syntax
	| lig_transition_syntax
;

r_transition_syntax: existing_reaction_state
        '[' '>' existing_reaction_state ']'
	',' WORLD ',' r_spec ',' t_spec ',' event_spec
{
	no_printf("\nReaction transition syntax:\n");
	p1=$<sym>1;
	p2=$<sym>4;
	rxp1=(struct rx *)p1->value;
	rxp2=(struct rx *)p2->value;
	prxp=rxp1->parent_rx;
	transition_index=-1;
	found_rxip=NULL;
        for (i=0;i<1+n_ligand_types;i++) {
          mrxip[0]=rxp1->bind_rx[i];
          mrxip[1]=rxp1->transport_rx[i];
          mrxip[2]=rxp1->product_poisson_rx[i];
          mrxip[3]=rxp1->dissoc_rx[i];
          mrxip[4]=rxp1->product_rx[i];
          mrxip[5]=rxp1->degrade_rx[i];
          for (j=0;j<6;j++) {
            if (mrxip[j]!=NULL) {
              nrx=mrxip[j]->n_rates;
              for (k=0;k<nrx;k++) {
                if(rxp2==mrxip[j]->next_rx[k]) {
	          transition_index=mrxip[j]->transition_index[k];
		  found_transition=k;
		  found_rxip=mrxip[j];
                }
              }
            }
          }
        }
        mrxip[0]=rxp1->isom_rx;
        if (mrxip[0]!=NULL) {
          nrx=mrxip[0]->n_rates;
          for (k=0;k<nrx;k++) {
            if(rxp2==mrxip[0]->next_rx[k]) {
	      transition_index=mrxip[0]->transition_index[k];
	      found_transition=k;
	      found_rxip=mrxip[0];
            }
          }
        }
	if (transition_index==-1) {
	  mdlerror("Invalid state transition specified");
	  return(1);
	}

	switch ($<tok>9) {
	  case OVER_E:
	    switch ($<tok>11) {
	      case SUM:
	        sprintf(mdlpvp->mdl_err_msg,"Count transition over all time not yet implemented\n");
	        mdl_warning(mdlpvp);
	        break;
	      case DT:
		switch ($<tok>13) {
		case A_EVENTS:
                  if (volp->iterations<0) {
                    sprintf(mdlpvp->mdl_err_msg,"Iterations = %d\n\tSetting iterations to 0\n",volp->iterations);
	            mdl_warning(mdlpvp);
                    volp->iterations=0;
                  }
		  i1=buffersize;
		  if ((intp=(int *)malloc(i1*sizeof(int)))==NULL) {
		    mdlerror("Cannot store count data");
		    return(1);
		  }
	          for (i=0;i<i1;i++) {
	            intp[i]=0;
	          }
		  clp->reset_flag=1;
		  clp->final_data=(void *)intp;
		  clp->data_type=INT;
		  clp->n_data=i1;
		  clp->temp_data
		    =(void *)&(found_rxip->transition_count_dt[found_transition]);
		  break;
		case INIT_EVENTS:
	          prxp->prev_state_flag=1;
                  if (volp->iterations<0) {
                    sprintf(mdlpvp->mdl_err_msg,"Iterations = %d\n\tSetting iterations to 0\n",volp->iterations);
	            mdl_warning(mdlpvp);
                    volp->iterations=0;
                  }
		  i1=buffersize;
		  
		  if ((intp=(int *)malloc(i1*sizeof(int)))==NULL) {
		    mdlerror("Cannot store count data");
		    return(1);
		  }
	          for (i=0;i<i1;i++) {
	            intp[i]=0;
	          }
		  clp->reset_flag=1;
		  clp->final_data=(void *)intp;
		  clp->data_type=INT;
		  clp->n_data=i1;
		  clp->temp_data
		    =(void *)&(found_rxip->init_transition_count_dt[found_transition]);
		  break;
		}
                break;
	      case CUM:
	        sprintf(mdlpvp->mdl_err_msg,"Count transition cumulatively not yet implemented\n");
	        mdl_warning(mdlpvp);
	        break;
	    } 
	    break;
	}
}
	| existing_reaction_state '[' '>' existing_reaction_state ']'
	',' existing_region ',' r_spec ',' t_spec ',' event_spec
{
	no_printf("\nReaction transition on region syntax:\n");
	fflush(stderr);
	p1=$<sym>1;
	p2=$<sym>4;
	mdlpvp->gp=$<sym>7;
	rxp1=(struct rx *)p1->value;
	rxp2=(struct rx *)p2->value;
	rp=(struct region *)mdlpvp->gp->value;

	switch ($<tok>9) {
	  case OVER_E:
	    switch ($<tok>11) {
	      case SUM:
	        sprintf(mdlpvp->mdl_err_msg,"Count transition over all time on region not yet implemented\n");
	        mdl_warning(mdlpvp);
	        break;
	      case DT:
		switch ($<tok>13) {
		case A_EVENTS:
                  if (volp->iterations<0) {
                    sprintf(mdlpvp->mdl_err_msg,"Iterations = %d\n\tSetting iterations to 0\n",volp->iterations);
	            mdl_warning(mdlpvp);
                    volp->iterations=0;
                  }
		  i1=buffersize;
		  defined_counter_flag=0;
		  if ((intp=(int *)malloc(i1*sizeof(int)))==NULL) {
		    mdlerror("Cannot store count data");
		    return(1);
		  }
	          for (i=0;i<i1;i++) {
	            intp[i]=0;
	          }
		  
		  /* Check if this counter already defined,
		   * if it is true, point the count_list to the
		   * defined rcrp2->counter
		   */
		  rcrlp2=rp->reg_counter_ref_list;
		  if (rcrlp2!=NULL) {
		    while (rcrlp2!=NULL) {
		      rcrp2=rcrlp2->reg_counter_ref;
		      if (rcrp2->state==rxp1&&rcrp2->next_state==rxp2&&rcrp2->count_type==TRANSITIONS&&rcrp->count_method==DT) {
			clp->final_data=(void *)intp;
			clp->data_type=INT;
			clp->n_data=i1;
			clp->temp_data=(void *)&(rcrp2->counter);
			defined_counter_flag=1;
		      }
		      rcrlp2=rcrlp2->next;
		    }
		  }
		  /* If this region counter is not defined before, 
		   * allocate memory for this region counter and
		   * assign it to the region.
		   */
		  rcrlp2=rp->reg_counter_ref_list;
		  if (rcrlp2==NULL || defined_counter_flag==0) {
		    if ((rcrp=(struct reg_counter_ref *)malloc
			(sizeof(struct reg_counter_ref)))==NULL) {
		      mdlerror("Can not save data for region counter");
		      return(1);
		    }
		    if ((rcrlp=(struct reg_counter_ref_list *)malloc
			(sizeof(struct reg_counter_ref_list)))==NULL) {
		      mdlerror("Can not save data for region counter");
		      return(1);
		    }
		    reg_counter_ref_head=NULL;
		    rcrp->counter=0;
		    rcrp->state=rxp1;
		    rcrp->parent=rp;
		    rcrp->next_state=rxp2;
		    rcrp->count_type=TRANSITIONS;
		    rcrp->count_method=DT;
		    rcrp->transition_count_each=NULL;
		    rcrp->next=reg_counter_ref_head;
		    reg_counter_ref_head=rcrp;
  
		    rcrlp->reg_counter_ref=rcrp;
		    rcrlp->next=rp->reg_counter_ref_list;
		    rp->reg_counter_ref_list=rcrlp;
  
		    clp->reset_flag=1;
		    clp->final_data=(void *)intp;
		    clp->data_type=INT;
		    clp->n_data=i1;
		    clp->temp_data =(void *)&(rcrp->counter);
		  }
		  break;
		case INIT_EVENTS:
	          prxp->prev_state_flag=1;
                  if (volp->iterations<0) {
                    sprintf(mdlpvp->mdl_err_msg,"Iterations = %d\n\tSetting iterations to 0\n",volp->iterations);
	            mdl_warning(mdlpvp);
                    volp->iterations=0;
                  }
		  i1=buffersize;
		  defined_counter_flag=0;
		  
		  if ((intp=(int *)malloc(i1*sizeof(int)))==NULL) {
		    mdlerror("Cannot store count data");
		    return(1);
		  }
	          for (i=0;i<i1;i++) {
	            intp[i]=0;
	          }
		  
		  /* Check if this counter already defined,
		   * if it is true, point the count_list to the
		   * defined rcrp2->counter
		   */
		  rcrlp2=rp->reg_counter_ref_list;
		  if (rcrlp2!=NULL) {
		    while (rcrlp2!=NULL) {
		      rcrp2=rcrlp2->reg_counter_ref;
		      if (rcrp2->state==rxp1&&rcrp2->next_state==rxp2&&rcrp2->count_type==INIT_TRANS&&rcrp->count_method==DT) {
			clp->final_data=(void *)intp;
			clp->data_type=INT;
			clp->n_data=i1;
			clp->temp_data=(void *)&(rcrp2->counter);
			defined_counter_flag=1;
		      }
		      rcrlp2=rcrlp2->next;
		    }
		  }
		  /* If this region counter is not defined before, 
		   * allocate memory for this region counter and
		   * assign it to the region.
		   */
		  rcrlp2=rp->reg_counter_ref_list;
		  if (rcrlp2==NULL || defined_counter_flag==0) {
		    if ((rcrp=(struct reg_counter_ref *)malloc
			(sizeof(struct reg_counter_ref)))==NULL) {
		      mdlerror("Can not save data for region counter");
		      return(1);
		    }
		    if ((rcrlp=(struct reg_counter_ref_list *)malloc
			(sizeof(struct reg_counter_ref_list)))==NULL) {
		      mdlerror("Can not save data for region counter");
		      return(1);
		    }
		    reg_counter_ref_head=NULL;
		    rcrp->counter=0;
		    rcrp->state=rxp1;
		    rcrp->parent=rp;
		    rcrp->next_state=rxp2;
		    rcrp->count_type=INIT_TRANS;
		    rcrp->count_method=DT;
		    rcrp->transition_count_each=NULL;
		    rcrp->next=reg_counter_ref_head;
		    reg_counter_ref_head=rcrp;

	  
		    rcrlp->reg_counter_ref=rcrp;
		    rcrlp->next=rp->reg_counter_ref_list;
		    rp->reg_counter_ref_list=rcrlp;
  
		    clp->reset_flag=1;
		    clp->final_data=(void *)intp;
		    clp->data_type=INT;
		    clp->n_data=i1;
		    clp->temp_data=(void *)&(rcrp->counter);
		  }

		  break;
		}
                break;
	      case CUM:

		switch ($<tok>13) {
		case A_EVENTS:
                  if (volp->iterations<0) {
                    sprintf(mdlpvp->mdl_err_msg,"Iterations = %d\n\tSetting iterations to 0\n",volp->iterations);
	            mdl_warning(mdlpvp);
                    volp->iterations=0;
                  }
		  i1=buffersize;
		  defined_counter_flag=0;
		  if ((intp=(int *)malloc(i1*sizeof(int)))==NULL) {
		    mdlerror("Cannot store count data");
		    return(1);
		  }
	          for (i=0;i<i1;i++) {
	            intp[i]=0;
	          }
		  
		  /* Check if this counter already defined,
		   * if it is true, point the count_list to the
		   * defined rcrp2->counter
		   */
		  rcrlp2=rp->reg_counter_ref_list;
		  if (rcrlp2!=NULL) {
		    while (rcrlp2!=NULL) {
		      rcrp2=rcrlp2->reg_counter_ref;
		      if (rcrp2->state==rxp1&&rcrp2->next_state==rxp2&&rcrp2->count_type==TRANSITIONS&&rcrp->count_method==CUM) {
			clp->final_data=(void *)intp;
			clp->data_type=INT;
			clp->n_data=i1;
			clp->temp_data=(void *)&(rcrp2->counter);
			defined_counter_flag=1;
		      }
		      rcrlp2=rcrlp2->next;
		    }
		  }
		  /* If this region counter is not defined before, 
		   * allocate memory for this region counter and
		   * assign it to the region.
		   */
		  rcrlp2=rp->reg_counter_ref_list;
		  if (rcrlp2==NULL || defined_counter_flag==0) {
		    if ((rcrp=(struct reg_counter_ref *)malloc
			(sizeof(struct reg_counter_ref)))==NULL) {
		      mdlerror("Can not save data for region counter");
		      return(1);
		    }
		    if ((rcrlp=(struct reg_counter_ref_list *)malloc
			(sizeof(struct reg_counter_ref_list)))==NULL) {
		      mdlerror("Can not save data for region counter");
		      return(1);
		    }
		    reg_counter_ref_head=NULL;
		    rcrp->counter=0;
		    rcrp->state=rxp1;
		    rcrp->parent=rp;
		    rcrp->next_state=rxp2;
		    rcrp->count_type=TRANSITIONS;
		    rcrp->count_method=CUM;
		    rcrp->transition_count_each=NULL;
		    rcrp->next=reg_counter_ref_head;
		    reg_counter_ref_head=rcrp;

	  
		    rcrlp->reg_counter_ref=rcrp;
		    rcrlp->next=rp->reg_counter_ref_list;
		    rp->reg_counter_ref_list=rcrlp;

		    clp->index_type=TIME_STAMP_VAL;
		    clp->final_data=(void *)intp;
		    clp->data_type=INT;
		    clp->n_data=i1;
		    clp->temp_data=(void *)&(rcrp->counter);
		  }
		  break;
		case INIT_EVENTS:
	          prxp->prev_state_flag=1;
                  if (volp->iterations<0) {
                    sprintf(mdlpvp->mdl_err_msg,"Iterations = %d\n\tSetting iterations to 0\n",volp->iterations);
	            mdl_warning(mdlpvp);
                    volp->iterations=0;
                  }
		  i1=buffersize;
		  defined_counter_flag=0;
		  if ((intp=(int *)malloc(i1*sizeof(int)))==NULL) {
		    mdlerror("Cannot store count data");
		    return(1);
		  }
	          for (i=0;i<i1;i++) {
	            intp[i]=0;
	          }

		  /* Check if this counter already defined,
		   * if it is true, point the count_list to the
		   * defined rcrp2->counter
		   */
		  rcrlp2=rp->reg_counter_ref_list;
		  if (rcrlp2!=NULL) {
		    while (rcrlp2!=NULL) {
		      rcrp2=rcrlp2->reg_counter_ref;
		      if (rcrp2->state==rxp1&&rcrp2->next_state==rxp2&&rcrp2->count_type==INIT_TRANS&&rcrp->count_method==CUM) {
			clp->final_data=(void *)intp;
			clp->data_type=INT;
			clp->n_data=i1;
			clp->temp_data=(void *)&(rcrp2->counter);
			defined_counter_flag=1;
		      }
		      rcrlp2=rcrlp2->next;
		    }
		  }
		  /* If this region counter is not defined before, 
		   * allocate memory for this region counter and
		   * assign it to the region.
		   */
		  rcrlp2=rp->reg_counter_ref_list;
		  if (rcrlp2==NULL || defined_counter_flag==0) {
		    if ((rcrp=(struct reg_counter_ref *)malloc
			(sizeof(struct reg_counter_ref)))==NULL) {
		      mdlerror("Can not save data for region counter");
		      return(1);
		    }
		    if ((rcrlp=(struct reg_counter_ref_list *)malloc
			(sizeof(struct reg_counter_ref_list)))==NULL) {
		      mdlerror("Can not save data for region counter");
		      return(1);
		    }
		    reg_counter_ref_head=NULL;
		    rcrp->counter=0;
		    rcrp->state=rxp1;
		    rcrp->parent=rp;
		    rcrp->next_state=rxp2;
		    rcrp->count_type=INIT_TRANS;
		    rcrp->count_method=CUM;
		    rcrp->transition_count_each=NULL;
		    rcrp->next=reg_counter_ref_head;
		    reg_counter_ref_head=rcrp;

	  
		    rcrlp->reg_counter_ref=rcrp;
		    rcrlp->next=rp->reg_counter_ref_list;
		    rp->reg_counter_ref_list=rcrlp;

		    clp->index_type=TIME_STAMP_VAL;
		    clp->final_data=(void *)intp;
		    clp->data_type=INT;
		    clp->n_data=i1;
		    clp->temp_data=(void *)&(rcrp->counter);
		  }
		  break;
	        } 
	        break;
	      }
	    break;
	}
}
	;

r_state_or_ligand_diffusion_syntax: existing_molecule_or_reaction_state
	',' WORLD ',' FOR_EACH_TIME_STEP
{
	p1=$<sym>1;
        switch (p1->sym_type) {
        case RX:
	  no_printf("\nReaction state syntax:\n");
	  fflush(stderr);
	  rxp1=(struct rx *)p1->value;
          if (volp->iterations<0) {
            sprintf(mdlpvp->mdl_err_msg,"Iterations = %d\n\tSetting iterations to 0\n",volp->iterations);
	    mdl_warning(mdlpvp);
            volp->iterations=0;
          }
	  i1=buffersize;
	  if ((intp=(int *)malloc(i1*sizeof(int)))==NULL) {
	    mdlerror("Cannot store count data");
	    return(1);
          }
	  for (i=0;i<i1;i++) {
	    intp[i]=0;
	  }
          clp->final_data=(void *)intp;
          clp->data_type=INT;
          clp->n_data=i1;
	  clp->temp_data=(void *)&(rxp1->count);
	  no_printf("Counting reaction state %s\n",p1->name);
	  fflush(stderr);
          break;
        case MOL:
          no_printf("\nWorld molecule diffusion syntax:\n");
          fflush(stderr);
          ligip=(struct ligand_info *)p1->value;
  
          if (world_lig_count==NULL) {
            if ((world_lig_count=(int *)malloc
                 ((1+n_ligand_types)*sizeof(int)))==NULL) {
              mdlerror("Cannot store world molecule count data");
              return(1); 
            }
	    for (i=0;i<(1+n_ligand_types);i++) {
              world_lig_count[i]=0;
            }
          }
  
          if (volp->iterations<0) {
            sprintf(mdlpvp->mdl_err_msg,"Iterations = %d\n\tSetting iterations to 0\n",volp->iterations);
	    mdl_warning(mdlpvp);
            volp->iterations=0;
          }

	  i1=buffersize;
	  if ((intp=(int *)malloc(i1*sizeof(int)))==NULL) {
	    mdlerror("Cannot store count data");
	    return(1);
          }
	  for (i=0;i<i1;i++) {
	    intp[i]=0;
	  }
          clp->final_data=(void *)intp;
          clp->data_type=INT;
          clp->n_data=i1;
	  clp->temp_data=(void *)&world_lig_count[ligip->type];
  
	  no_printf("Counting molecule %s in WORLD\n",mdlpvp->gp->name);
	  fflush(stderr);
          break;
        }
}
	| existing_molecule_or_reaction_state
	  ',' existing_object ',' FOR_EACH_TIME_STEP
{
	p1=$<sym>1;
        switch (p1->sym_type) {
        case RX:
  
	  mdlerror("Counting of effector states within objects not yet implemented");
	  return(1);
	  break;
        case MOL:
	  no_printf("\nLigand diffusion syntax:\n");
	  fflush(stderr);
	  ligip=(struct ligand_info *)p1->value;
	  mdlpvp->tp=$<sym>3;
	  mdlpvp->objp=(struct object *)mdlpvp->tp->value;
          mdlpvp->objp2=mdlpvp->top_objp;
  
          if (volp->iterations<0) {
            sprintf(mdlpvp->mdl_err_msg,"Iterations = %d\n\tSetting iterations to 0\n",volp->iterations);
	    mdl_warning(mdlpvp);
            volp->iterations=0;
          }

          freq=(int) (output_freq+ROUND_UP);
	  clp->freq=freq;
	  clp->n_output=buffersize;
          clp->reset_flag=0;
	  clp->update_flag=0;
	  clp->data_type=EXPR;
	  clp->index_type=UNKNOWN;
	  clp->n_data=0;
	  clp->temp_data=NULL;
	  clp->final_data=NULL;
	  clp->operand1=NULL;
	  clp->operand2=NULL;
	  clp->oper='+';
	  if (build_lig_count_tree(mdlpvp->objp,mdlpvp->objp2,clp,ligip->type,buffersize,prefix_name)) {
	    mdlerror("Cannot store molecule count_list data");
	    return(1);
	  }

	  no_printf("Counting molecule %s in compartment %s\n",mdlpvp->gp->name,mdlpvp->tp->name);
	  fflush(stderr);
	  break;
	}
}
	| existing_molecule_or_reaction_state
	  ',' existing_region ',' FOR_EACH_TIME_STEP
{	  
	p1=$<sym>1;
	mdlpvp->gp=$<sym>3;
	switch(p1->sym_type) {
	case RX:
	  no_printf("Counting rx state on region syntax: \n");
	  fflush(stderr);
	  rxp1=(struct rx *)p1->value;
	  rp=(struct region*)mdlpvp->gp->value;
	  
          if (volp->iterations<0) {
            sprintf(mdlpvp->mdl_err_msg,"Iterations = %d\n\tSetting iterations to 0\n",volp->iterations);
	    mdl_warning(mdlpvp);
            volp->iterations=0;
          }
	  defined_counter_flag=0;
	  i1=buffersize;
	  if ((intp=(int *)malloc(buffersize*sizeof(int)))==NULL) {
	    mdlerror("Cannot store count data");
	    return(1);
          }
	  for (i=0;i<i1;i++) {
	    intp[i]=0;
	  }

	  /* Check if this counter already defined,
	   * if it is true, point the count_list to the
	   * defined rcrp2->counter
	   */
	  rcrlp2=rp->reg_counter_ref_list;
	  if (rcrlp2!=NULL) {
	    while (rcrlp2!=NULL) {
	      rcrp2=rcrlp2->reg_counter_ref;
	      if (rcrp2->state==rxp1&&rcrp2->count_type==RX_STATE) {
		clp->final_data=(void *)intp;
		clp->data_type=INT;
		clp->n_data=i1;
		clp->temp_data=(void *)&(rcrp2->counter);
		defined_counter_flag=1;
	      }
	      rcrlp2=rcrlp2->next;
	    }
	  }
	  /* If this region counter is not defined before, 
	   * allocate memory for this region counter and
	   * assign it to the region.
	   */
	  rcrlp2=rp->reg_counter_ref_list;
	  if (rcrlp2==NULL || defined_counter_flag==0){
	    if ((rcrp=(struct reg_counter_ref *)malloc
	  	(sizeof(struct reg_counter_ref)))==NULL) {
	      mdlerror("Can not save data for region counter");
	      return(1);
	    }
	    if ((rcrlp=(struct reg_counter_ref_list *)malloc
	  	(sizeof(struct reg_counter_ref_list)))==NULL) {
	      mdlerror("Can not save data for region counter");
	      return(1);
	    }
	    reg_counter_ref_head=NULL;
	    rcrp->counter=0;
	    rcrp->state=rxp1;
	    rcrp->parent=rp;
	    rcrp->next_state=NULL;
	    rcrp->count_type=RX_STATE;
	    rcrp->count_method=DT;
	    rcrp->transition_count_each=NULL;
	    rcrp->next=reg_counter_ref_head;
	    reg_counter_ref_head=rcrp;
	    
	    rcrlp->reg_counter_ref=rcrp;
	    rcrlp->next=rp->reg_counter_ref_list;
	    rp->reg_counter_ref_list=rcrlp;
	    
	    clp->final_data=(void *)intp;
	    clp->data_type=INT;
	    clp->n_data=i1;
	    clp->temp_data=(void *)&(rcrp->counter);
	  }

	  no_printf("Counting state %s in region %s\n",p1->name,mdlpvp->gp->name);
	  fflush(stderr);
	  
	  break;
	case MOL:
	  sprintf(mdlpvp->mdl_err_msg,"Count ligand on object surface region is not possible\n");
	  mdl_warning(mdlpvp);

	  return(1);
	  break;
	}
};

lig_hit_syntax:existing_molecule_or_reaction_state ',' existing_region ',' 
	  FOR_EACH_TIME_STEP ',' ALL_HITS
{
	p1=$<sym>1;
	mdlpvp->gp=$<sym>3;

	switch(p1->sym_type) {
	case RX:
	  sprintf(mdlpvp->mdl_err_msg,"This case doesn't exist!\n");
	  mdl_warning(mdlpvp);
	  return(1);
 	  break;
	case MOL:
	  ligip=(struct ligand_info *)p1->value;
	  rp=(struct region*)mdlpvp->gp->value;
	  
          if (volp->iterations<0) {
            sprintf(mdlpvp->mdl_err_msg,"Iterations = %d\n\tSetting iterations to 0\n",volp->iterations);
	    mdl_warning(mdlpvp);
            volp->iterations=0;
          }
	  i1=buffersize;
	  if ((intp=(int *)malloc(buffersize*sizeof(int)))==NULL) {
	    mdlerror("Cannot store count data");
	    return(1);
          }
	  for (i=0;i<i1;i++) {
	    intp[i]=0;
	  }

	  if (rp->lig_hit_counter==NULL) {
	    if ((rp->lig_hit_counter=(struct lig_hit_counter **)malloc((1+n_ligand_types)*sizeof(struct lig_hit_counter *)))==NULL) {
	        mdlerror("Cannot store ligand hit count data");
	        return(1);
	    }
	    for(i1=0;i1<1+n_ligand_types;i1++) {
	      rp->lig_hit_counter[i1]=NULL;
	    }
	  }
	  lig_hit=rp->lig_hit_counter[ligip->type];
	  if (rp->lig_hit_counter[ligip->type]==NULL) {
	    if ((rp->lig_hit_counter[ligip->type]=
	                  (struct lig_hit_counter *)malloc
	                  (sizeof(struct lig_hit_counter)))==NULL) {
	      mdlerror("Cannot store transition count data");
	      return(1);
	    }
	    rp->lig_hit_counter[ligip->type]->ligand=ligip;
	    rp->lig_hit_counter[ligip->type]->counter=0;

	    clp->data_type=INT;
	    clp->update_flag=1;
	    clp->reset_flag=1;
	    clp->index_type=TIME_STAMP_VAL;

	    clp->final_data=(void *)intp;
	    clp->data_type=INT;
	    clp->n_data=i1;
	    clp->temp_data=(void *)&(rp->lig_hit_counter[ligip->type]->counter);
	  }
	  else  {
	    sprintf(mdlpvp->mdl_err_msg,"Duplicated request for ligand hit counting!\n");
	    mdl_warning(mdlpvp);
	  }

	  no_printf("Counting ligand %s hits with  region %s\n",p1->name,mdlpvp->gp->name);
	  fflush(stderr);
	  break;
	}
}
;

existing_molecule_or_reaction_state: VAR
{
  if (mdlpvp->cval_2!=NULL) {
    mdlpvp->sym_name=mdlpvp->cval_2;
  }       
  else {  
    mdlpvp->sym_name=mdlpvp->cval;  
  }       
  if ((mdlpvp->gp=retrieve_sym(mdlpvp->sym_name,RX,volp->main_sym_table))==NULL) {
    if ((mdlpvp->gp=retrieve_sym(mdlpvp->sym_name,MOL,volp->main_sym_table))==NULL) {
      sprintf(mdlpvp->mdl_err_msg,"%s %s","Undefined molecule or reaction state:",mdlpvp->sym_name);
      mdlerror(mdlpvp->mdl_err_msg,mdlpvp);
      if (mdlpvp->sym_name==mdlpvp->cval) {
        mdlpvp->cval=NULL;
      }
      else {
        mdlpvp->cval_2=NULL;
      }
      free((void *)mdlpvp->sym_name);
      return(1);
    }
  }
  if (mdlpvp->sym_name==mdlpvp->cval) {
    mdlpvp->cval=NULL;
  }
  else {
    mdlpvp->cval_2=NULL;
  }
  free((void *)mdlpvp->sym_name);
#ifdef KELP
  mdlpvp->gp->ref_count++;
  no_printf("ref_count: %d\n",mdlpvp->gp->ref_count);
#endif
  $$=mdlpvp->gp;
};

lig_transition_syntax: existing_reaction_state
	'[' '>' existing_reaction_state ']'
        ',' WORLD ',' lig_spec ',' t_spec
{
	no_printf("\nLigand transition syntax:\n");
	fflush(stderr);
	p1=$<sym>1;
	p2=$<sym>4;
	rxp1=(struct rx *)p1->value;
	rxp2=(struct rx *)p2->value;
	prxp=rxp1->parent_rx;
	transition_index=-1;
	found_rxip=NULL;
        for (i=0;i<1+n_ligand_types;i++) {
          mrxip[0]=rxp1->bind_rx[i];
          mrxip[1]=rxp1->transport_rx[i];
          mrxip[2]=rxp1->product_poisson_rx[i];
          mrxip[3]=rxp1->dissoc_rx[i];
          mrxip[4]=rxp1->product_rx[i];
          mrxip[5]=rxp1->degrade_rx[i];
          for (j=0;j<6;j++) {
            if (mrxip[j]!=NULL) {
              nrx=mrxip[j]->n_rates;
              for (k=0;k<nrx;k++) {
                if(rxp2==mrxip[j]->next_rx[k]) {
	          transition_index=mrxip[j]->transition_index[k];
		  found_transition=k;
		  found_rxip=mrxip[j];
                  ligip=ligand_table[i];
                }
              }
            }
          }
        }
	if (transition_index==-1) {
	  mdlerror("Invalid state transition specified");
	  return(1);
	}
	no_printf("transition index = %d\n",transition_index);
	fflush(stderr);
	switch ($<tok>9) {
	case EACH_L:
	  switch ($<tok>11) {
	  case SUM:
	    if (ligip->transition_count_each==NULL) {
	      i1=n_rx_types;
	      if ((ligip->transition_count_each=
	                  (struct lig_transition_count **)malloc
	                  (i1*sizeof(struct lig_transition_count *)))==NULL) {
	        mdlerror("Cannot store transition count data");
	        return(1);
	      }
	      for(i1=0;i1<n_rx_types;i1++) {
	        ligip->transition_count_each[i1]=NULL;
	      }
	    }
	    if (ligip->transition_count_each[prxp->rx_index]==NULL) {
	      if ((ligip->transition_count_each[prxp->rx_index]=
	                  (struct lig_transition_count *)malloc
	                  (sizeof(struct lig_transition_count)))==NULL) {
	        mdlerror("Cannot store transition count data");
	        return(1);
	      }
	      ligip->transition_count_each[prxp->rx_index]->sum=NULL;
	      ligip->transition_count_each[prxp->rx_index]->dt=NULL;
	      ligip->transition_count_each[prxp->rx_index]->cum=NULL;
	    }
	    if (ligip->transition_count_each[prxp->rx_index]->sum==NULL) {
	      i1=prxp->num_transitions;
	      if ((ligip->transition_count_each[prxp->rx_index]->sum=
	                  (struct count_list **)malloc
	                  (i1*sizeof(struct count_list *)))==NULL) {
	        mdlerror("Cannot store transition count data");
	        return(1);
	      }
	      for(i1=0;i1<prxp->num_transitions;i1++) {
	        ligip->transition_count_each[prxp->rx_index]->sum[i1]=NULL;
	      }
	    }
	    if (ligip->transition_count_each[prxp->rx_index]->sum[transition_index]
	        !=NULL) {
	      sprintf(mdlpvp->mdl_err_msg,"Duplicated request to count a type of molecule transition\n");
	      mdl_warning(mdlpvp);
	    }
	    ligip->transition_count_each[prxp->rx_index]->sum[transition_index]=clp;
	    clp->data_type=INT;
	    clp->update_flag=0;
	    clp->index_type=INDEX_VAL;
	    clp->temp_data=NULL; /*temp_data not required in this case*/ 
	    clp->final_data=NULL; /*final_data is allocated at 
	                            time of vesicle release only*/ 
	    clp->n_data=0; /*n_data is known at time of vesicle release only*/
	    break;
	  case DT:
	    sprintf(mdlpvp->mdl_err_msg,"Count molecule FOR_EACH_MOLECULE,FOR_EACH_TIME_STEP not yet implemented\n");
	    mdl_warning(mdlpvp);
	    break;
	  case CUM:
	    sprintf(mdlpvp->mdl_err_msg,"Count molecule FOR_EACH_MOLECULE,CUMULATE_FOR_EACH_TIME_STEP not yet implemented\n");
	    mdl_warning(mdlpvp);
	    break;
	  } 
	  break;
	case OVER_L:
	  switch ($<tok>11) {
	  case SUM:
	    sprintf(mdlpvp->mdl_err_msg,"Count ligand transitions SUM_OVER_ALL_MOLECULES,SUM_OVER_ALL_TIME_STEPS not yet implemented\n");
	    mdl_warning(mdlpvp);
	    break;
	  case DT:
	    if (ligip->transition_count_all==NULL) {
	      i1=n_rx_types;
	      if ((ligip->transition_count_all=
	                  (struct lig_transition_count **)malloc
	                  (i1*sizeof(struct lig_transition_count *)))==NULL) {
	        mdlerror("Cannot store transition count data");
	        return(1);
	      }
	      for(i1=0;i1<n_rx_types;i1++) {
	        ligip->transition_count_all[i1]=NULL;
	      }
	    }
	    if (ligip->transition_count_all[prxp->rx_index]==NULL) {
	      if ((ligip->transition_count_all[prxp->rx_index]=
	                  (struct lig_transition_count *)malloc
	                  (sizeof(struct lig_transition_count)))==NULL) {
	        mdlerror("Cannot store transition count data");
	        return(1);
	      }
	      ligip->transition_count_all[prxp->rx_index]->sum=NULL;
	      ligip->transition_count_all[prxp->rx_index]->dt=NULL;
	      ligip->transition_count_all[prxp->rx_index]->cum=NULL;
	    }
	    if (ligip->transition_count_all[prxp->rx_index]->dt==NULL) {
	      i1=prxp->num_transitions;
	      if ((ligip->transition_count_all[prxp->rx_index]->dt=
	                  (struct count_list **)malloc
	                  (i1*sizeof(struct count_list *)))==NULL) {
	        mdlerror("Cannot store transition count data");
	        return(1);
	      }
	      for(i1=0;i1<prxp->num_transitions;i1++) {
	        ligip->transition_count_all[prxp->rx_index]->dt[i1]=NULL;
	      }
	    }
	    if (ligip->transition_count_all[prxp->rx_index]->dt[transition_index]
	        !=NULL) {
	      sprintf(mdlpvp->mdl_err_msg,"Duplicated request to count this type of molecule transition\n");
	      mdl_warning(mdlpvp);
	    }
	    ligip->transition_count_all[prxp->rx_index]->dt[transition_index]=clp;
	    clp->data_type=INT;
	    clp->update_flag=1;
	    clp->reset_flag=1;
	    clp->index_type=TIME_STAMP_VAL;

            if (volp->iterations<0) {
              sprintf(mdlpvp->mdl_err_msg,"Iterations = %d\n\tSetting iterations to 0\n",volp->iterations);
	      mdl_warning(mdlpvp);
              volp->iterations=0;
            }

	    i1=buffersize;
	    if ((intp=(int *)malloc(i1*sizeof(int)))==NULL) {
	      mdlerror("Cannot store count data");
	      return(1);
            }
	    for (i=0;i<i1;i++) {
	      intp[i]=0;
	    }
            clp->final_data=(void *)intp;
            clp->n_data=i1;
    
	    if ((intp=(int *)malloc(sizeof(int)))==NULL) {
	      mdlerror("Cannot store count data");
	      return(1);
            }
	    clp->temp_data=(void *)intp;
	    *(int *)clp->temp_data=0;
	    break;
	  case CUM:
	    if (ligip->transition_count_all==NULL) {
	      i1=n_rx_types;
	      if ((ligip->transition_count_all=
	                  (struct lig_transition_count **)malloc
	                  (i1*sizeof(struct lig_transition_count *)))==NULL) {
	        mdlerror("Cannot store transition count data");
	        return(1);
	      }
	      for(i1=0;i1<n_rx_types;i1++) {
	        ligip->transition_count_all[i1]=NULL;
	      }
	    }
	    if (ligip->transition_count_all[prxp->rx_index]==NULL) {
	      if ((ligip->transition_count_all[prxp->rx_index]=
	                  (struct lig_transition_count *)malloc
	                  (sizeof(struct lig_transition_count)))==NULL) {
	        mdlerror("Cannot store transition count data");
	        return(1);
	      }
	      ligip->transition_count_all[prxp->rx_index]->sum=NULL;
	      ligip->transition_count_all[prxp->rx_index]->dt=NULL;
	      ligip->transition_count_all[prxp->rx_index]->cum=NULL;
	    }
	    if (ligip->transition_count_all[prxp->rx_index]->cum==NULL) {
	      i1=prxp->num_transitions;
	      if ((ligip->transition_count_all[prxp->rx_index]->cum=
	                  (struct count_list **)malloc
	                  (i1*sizeof(struct count_list *)))==NULL) {
	        mdlerror("Cannot store transition count data");
	        return(1);
	      }
	      for(i1=0;i1<prxp->num_transitions;i1++) {
	        ligip->transition_count_all[prxp->rx_index]->cum[i1]=NULL;
	      }
	    }
	    if (ligip->transition_count_all[prxp->rx_index]->cum[transition_index]
	        !=NULL) {
	      sprintf(mdlpvp->mdl_err_msg,"Duplicated request to count this type of molecule transition\n");
	      mdl_warning(mdlpvp);
	    }
	    ligip->transition_count_all[prxp->rx_index]->cum[transition_index]=clp;
	    clp->data_type=INT;
	    clp->update_flag=1;
	    clp->reset_flag=0;
	    clp->index_type=TIME_STAMP_VAL;

            if (volp->iterations<0) {
              sprintf(mdlpvp->mdl_err_msg,"Iterations = %d\n\tSetting iterations to 0\n",volp->iterations);
	      mdl_warning(mdlpvp);
              volp->iterations=0;
            }

	    i1=buffersize;
	    if ((intp=(int *)malloc(i1*sizeof(int)))==NULL) {
	      mdlerror("Cannot store count data");
	      return(1);
            }
	    for (i=0;i<i1;i++) {
	      intp[i]=0;
	    }
            clp->final_data=(void *)intp;
            clp->n_data=i1;
    
	    if ((intp=(int *)malloc(sizeof(int)))==NULL) {
	      mdlerror("Cannot store count data");
	      return(1);
            }
	    clp->temp_data=(void *)intp;
	    *(int *)clp->temp_data=0;
	    break;
	  }
	  break;
	}
}
	| existing_reaction_state
	'[' '>' existing_reaction_state ']'
        ',' existing_region ',' lig_spec ',' t_spec
{
	no_printf("\nLigand transition over surface region syntax:\n");
	fflush(stderr);
	p1=$<sym>1;
	p2=$<sym>4;
	mdlpvp->gp=$<sym>7;
	rxp1=(struct rx *)p1->value;
	rxp2=(struct rx *)p2->value;
	rp=(struct region *)mdlpvp->gp->value;
	prxp=rxp1->parent_rx;
	transition_index=-1;
	found_rxip=NULL;
        for (i=0;i<1+n_ligand_types;i++) {
          mrxip[0]=rxp1->bind_rx[i];
          mrxip[1]=rxp1->transport_rx[i];
          mrxip[2]=rxp1->product_poisson_rx[i];
          mrxip[3]=rxp1->dissoc_rx[i];
          mrxip[4]=rxp1->product_rx[i];
          mrxip[5]=rxp1->degrade_rx[i];
          for (j=0;j<6;j++) {
            if (mrxip[j]!=NULL) {
              nrx=mrxip[j]->n_rates;
              for (k=0;k<nrx;k++) {
                if(rxp2==mrxip[j]->next_rx[k]) {
	          transition_index=mrxip[j]->transition_index[k];
		  found_transition=k;
		  found_rxip=mrxip[j];
                  ligip=ligand_table[i];
                }
              }
            }
          }
        }
	if (transition_index==-1) {
	  mdlerror("Invalid state transition specified");
	  return(1);
	}
	no_printf("transition index = %d\n",transition_index);
	fflush(stderr);
	switch ($<tok>9) {
	case EACH_L:
	  switch ($<tok>11) {
	  case SUM:
	  /* Check if this counter already defined,
	   * if it is true, point the count_list to the
	   * defined rcrp2->counter
	   */
	  defined_counter_flag=0;  
	  rcrlp2=rp->reg_counter_ref_list;
	  if (rcrlp2!=NULL) {
	    while (rcrlp2!=NULL) {
	      rcrp2=rcrlp2->reg_counter_ref;
	      if (rcrp2->state==rxp1&&rcrp2->next_state==rxp2&&rcrp2->count_type==MOL_TRANS_EACH&&rcrp2->count_method==SUM) {
		rcrp2->transition_count_each[prxp->rx_index]->sum[transition_index]=clp;

		clp->reset_flag=1;
		clp->update_flag=0;
		clp->index_type=INDEX_VAL;
		clp->data_type=INT;
		clp->n_data=0;
		clp->temp_data=NULL;
		clp->final_data=NULL;		
		defined_counter_flag=1;
	      }
	      rcrlp2=rcrlp2->next;
	    }
	  }
	  /* If this region counter is not defined before, 
	   * allocate memory for this region counter and
	   * assign it to the region.
	   */
	  rcrlp2=rp->reg_counter_ref_list;
	  if (rcrlp2==NULL || defined_counter_flag==0){
	    if ((rcrp=(struct reg_counter_ref *)malloc
		(sizeof(struct reg_counter_ref)))==NULL) {
	      mdlerror("Can not save data for region counter");
	      return(1);
	    }
	    if ((rcrlp=(struct reg_counter_ref_list *)malloc
		(sizeof(struct reg_counter_ref_list)))==NULL) {
	      mdlerror("Can not save data for region counter");
	      return(1);
	    }

	    reg_counter_ref_head=NULL;
	    rcrp->counter=0;
	    rcrp->state=rxp1;
	    rcrp->parent=rp;
	    rcrp->next_state=rxp2;
	    rcrp->count_type=MOL_TRANS_EACH;
	    rcrp->count_method=SUM;
	    rcrp->transition_count_each=NULL;
	    rcrp->next=reg_counter_ref_head;
	    reg_counter_ref_head=rcrp;


  	    
	    rcrlp->reg_counter_ref=rcrp;
	    rcrlp->next=rp->reg_counter_ref_list;
	    rp->reg_counter_ref_list=rcrlp;
  
	    /* Save reg_counter_ref_list for ligand transition 
	     * FOR_EACH_MOLECULE, SUM_OVER_ALL_TIME_STEPS to 
	     * ligip->region_transition_count_each
	     */  
	    if ((ligrcrlp=(struct reg_counter_ref_list *)malloc
		(sizeof(struct reg_counter_ref_list)))==NULL) {
		mdlerror("Can not save data for region counter");
		return(1);
	    }	      
	    ligrcrlp->reg_counter_ref=rcrp;
	    ligrcrlp->next=ligip->region_transition_count_each;
	    ligip->region_transition_count_each=ligrcrlp;
	    
	    if (rcrp->transition_count_each==NULL) {
	      i1=n_rx_types;
	      if ((rcrp->transition_count_each=
	                  (struct lig_transition_count **)malloc
	                  (i1*sizeof(struct lig_transition_count)))==NULL) {
	        mdlerror("Cannot store transition count data");
	        return(1);

	       }

	       for (i1=0;i1<n_rx_types;i1++) {
		 rcrp->transition_count_each[i1]=NULL;
	       }
	    }
	    if (rcrp->transition_count_each[prxp->rx_index]==NULL) {
	      if ((rcrp->transition_count_each[prxp->rx_index]=
	                  (struct lig_transition_count *)malloc
	                  (sizeof(struct lig_transition_count)))==NULL) {
	        mdlerror("Cannot store transition count data");
	        return(1);
	      }
	      rcrp->transition_count_each[prxp->rx_index]->sum=NULL;
	      rcrp->transition_count_each[prxp->rx_index]->dt=NULL;
	      rcrp->transition_count_each[prxp->rx_index]->cum=NULL;
	    }
	    if (rcrp->transition_count_each[prxp->rx_index]->sum==NULL) {
	      i1=prxp->num_transitions;
	      if ((rcrp->transition_count_each[prxp->rx_index]->sum=
	                  (struct count_list **)malloc
	                  (i1*sizeof(struct count_list *)))==NULL) {
	        mdlerror("Cannot store transition count data");
	        return(1);
	      }
	      for(i1=0;i1<prxp->num_transitions;i1++) {
	        rcrp->transition_count_each[prxp->rx_index]->sum[i1]=NULL;
	      }
	    }
	    if (rcrp->transition_count_each[prxp->rx_index]->sum[transition_index]
	        !=NULL) {
	      sprintf(mdlpvp->mdl_err_msg,"Duplicated request to count a type of molecule transition\n");
	      mdl_warning(mdlpvp);
	    }
	    rcrp->transition_count_each[prxp->rx_index]->sum[transition_index]=clp;

	    clp->reset_flag=1;
	    clp->update_flag=0;
	    clp->index_type=INDEX_VAL;
	    clp->data_type=INT;
	    clp->n_data=0;
	    clp->temp_data=NULL;
	    clp->final_data=NULL;
	  }
	    break;
	  case DT:
	    sprintf(mdlpvp->mdl_err_msg,"Count transition FOR_EACH_MOLECULE,FOR_EACH_TIME_STEP on region not yet implemented\n");
	    mdl_warning(mdlpvp);
	    break;
	  case CUM:
	    sprintf(mdlpvp->mdl_err_msg,"Count transition FOR_EACH_MOLECULE,CUMULATE_FOR_EACH_TIME_STEP on region not yet implemented\n");
	    mdl_warning(mdlpvp);
	    break;
	  } 
	  break;
	case OVER_L:
	  switch ($<tok>11) {
	  case SUM:
	    sprintf(mdlpvp->mdl_err_msg,"Count transitions SUM_OVER_ALL_MOLECULES,SUM_OVER_ALL_TIME_STEPS not yet implemented\n");
	    mdl_warning(mdlpvp);
	    break;
	  case DT:
	    if (volp->iterations<0) {
              sprintf(mdlpvp->mdl_err_msg,"Iterations = %d\n\tSetting iterations to 0\n",volp->iterations);
	      mdl_warning(mdlpvp);
              volp->iterations=0;
            }

	    i1=buffersize;
	    defined_counter_flag=0;
	    if ((intp=(int *)malloc(i1*sizeof(int)))==NULL) {
	      mdlerror("Cannot store count data");
	      return(1);
            }
	    for (i=0;i<i1;i++) {
	      intp[i]=0;
	    }

	  /* Check if this counter already defined,
	   * if it is true, point the count_list to the
	   * defined rcrp2->counter
	   */
	  rcrlp2=rp->reg_counter_ref_list;
	  if (rcrlp2!=NULL) {
	    while (rcrlp2!=NULL) {
	      rcrp2=rcrlp2->reg_counter_ref;
	      if (rcrp2->state==rxp1&&rcrp2->next_state==rxp2&&rcrp2->count_type==MOL_TRANS_ALL&&rcrp2->count_method==DT) {
		clp->data_type=INT;
		clp->update_flag=1;
		clp->reset_flag=1;
		clp->index_type=TIME_STAMP_VAL;	
		clp->final_data=(void *)intp;
		clp->n_data=i1;
		clp->temp_data=(void *)&(rcrp2->counter);
		defined_counter_flag=1;
	      }
	      rcrlp2=rcrlp2->next;
	    }
	  }
	  /* If this region counter is not defined before, 
	   * allocate memory for this region counter and
	   * assign it to the region.
	   */
	  rcrlp2=rp->reg_counter_ref_list;
	  if (rcrlp2==NULL || defined_counter_flag==0){
	    if ((rcrp=(struct reg_counter_ref *)malloc
		(sizeof(struct reg_counter_ref)))==NULL) {
	      mdlerror("Can not save data for region counter");
	      return(1);
	    }
	    if ((rcrlp=(struct reg_counter_ref_list *)malloc
		(sizeof(struct reg_counter_ref_list)))==NULL) {
	      mdlerror("Can not save data for region counter");
	      return(1);
	    }
	    reg_counter_ref_head=NULL;
	    rcrp->counter=0;
	    rcrp->state=rxp1;
	    rcrp->parent=rp;
	    rcrp->next_state=rxp2;
	    rcrp->count_type=MOL_TRANS_ALL;
	    rcrp->count_method=DT;
	    rcrp->transition_count_each=NULL;
	    rcrp->next=reg_counter_ref_head;
	    reg_counter_ref_head=rcrp;

	  
	    rcrlp->reg_counter_ref=rcrp;
	    rcrlp->next=rp->reg_counter_ref_list;
	    rp->reg_counter_ref_list=rcrlp;
  
	    
	    clp->data_type=INT;
	    clp->update_flag=1;
	    clp->reset_flag=1;
	    clp->index_type=TIME_STAMP_VAL;
            clp->final_data=(void *)intp;
            clp->n_data=i1;
	    clp->temp_data=(void *)&(rcrp->counter);
	  }
	    break;
	  case CUM:
            if (volp->iterations<0) {
              sprintf(mdlpvp->mdl_err_msg,"Iterations = %d\n\tSetting iterations to 0\n",volp->iterations);
	      mdl_warning(mdlpvp);
              volp->iterations=0;
            }

	    i1=buffersize;
	    defined_counter_flag=0;
	    if ((intp=(int *)malloc(i1*sizeof(int)))==NULL) {
	      mdlerror("Cannot store count data");
	      return(1);
            }
	    for (i=0;i<i1;i++) {
	      intp[i]=0;
	    }
	    
	    /* Check if this counter already defined,
	     * if it is true, point the count_list to the
	     * defined rcrp2->counter
	     */
	    rcrlp2=rp->reg_counter_ref_list;
	    if (rcrlp2!=NULL) {
	      while (rcrlp2!=NULL) {
		rcrp2=rcrlp2->reg_counter_ref;
		if (rcrp2->state==rxp1&&rcrp2->next_state==rxp2&&rcrp2->count_type==MOL_TRANS_ALL&&rcrp2->count_method==CUM) {
		  clp->final_data=(void *)intp;
		  clp->n_data=i1;
		  clp->temp_data=(void *)&(rcrp2->counter);
	    
		  clp->data_type=INT;
		  clp->update_flag=1;
		  clp->reset_flag=0;
		  clp->index_type=TIME_STAMP_VAL;
		  defined_counter_flag=1;
		}
		rcrlp2=rcrlp2->next;
	      }
	    }
	    /* If this region counter is not defined before, 
	     * allocate memory for this region counter and
	     * assign it to the region.
	     */
	    rcrlp2=rp->reg_counter_ref_list;
	    if (rcrlp2==NULL || defined_counter_flag==0){
	      if ((rcrp=(struct reg_counter_ref *)malloc
		(sizeof(struct reg_counter_ref)))==NULL) {
		mdlerror("Can not save data for region counter");
		return(1);
	      }
	      if ((rcrlp=(struct reg_counter_ref_list *)malloc
		(sizeof(struct reg_counter_ref_list)))==NULL) {
		mdlerror("Can not save data for region counter");
		return(1);
	      }
	      reg_counter_ref_head=NULL;
	      rcrp->counter=0;
	      rcrp->state=rxp1;
	      rcrp->parent=rp;
	      rcrp->next_state=rxp2;
	      rcrp->count_type=MOL_TRANS_ALL;
	      rcrp->count_method=CUM;
	      rcrp->transition_count_each=NULL;
	      rcrp->next=reg_counter_ref_head;
	      reg_counter_ref_head=rcrp;

	  
	      rcrlp->reg_counter_ref=rcrp;
	      rcrlp->next=rp->reg_counter_ref_list;
	      rp->reg_counter_ref_list=rcrlp;

	      clp->final_data=(void *)intp;
	      clp->n_data=i1;
	      clp->temp_data=(void *)&(rcrp->counter);
	    
	      clp->data_type=INT;
	      clp->update_flag=1;
	      clp->reset_flag=0;
	      clp->index_type=TIME_STAMP_VAL;
	    }
	    break;
	  }
	  break;
	}
}
	  
        ;

lig_spec: SUM_OVER_ALL_MOLECULES {$$=OVER_L;}
	| FOR_EACH_MOLECULE {$$=EACH_L;}
	| SPECIFIED_MOLECULES {$$=SPEC_L;}
	;

r_spec: SUM_OVER_ALL_EFFECTORS {$$=OVER_E;}
	| FOR_EACH_EFFECTOR {$$=EACH_E;}
	| SPECIFIED_EFFECTORS {$$=SPEC_E;}
	;

t_spec: SUM_OVER_ALL_TIME_STEPS {$$=SUM;}
	| FOR_EACH_TIME_STEP {$$=DT;}
	| CUMULATE_FOR_EACH_TIME_STEP {$$=CUM;}
	;

event_spec: ALL_EVENTS {$$=A_EVENTS;}
	| INITIAL_EVENTS {$$=INIT_EVENTS;}
	| INTERIM_EVENTS {$$=INTER_EVENTS;}
	;

io_stmt: fopen_stmt
	| fclose_stmt
	| printf_stmt
	| fprintf_stmt
	| sprintf_stmt
	| print_time_stmt
	| fprint_time_stmt
	;

fopen_stmt: new_file_stream FOPEN '(' file_name ',' file_mode ')'
{
mdlpvp->gp=$<sym>1;
filep=(struct file_stream *)mdlpvp->gp->value;
filep->name=$<str>4;
a_str=$<str>6;
if ((filep->stream=fopen(filep->name,a_str))==NULL) {
  sprintf(mdlpvp->mdl_err_msg,"%s %s","Cannot open file:",filep->name);
  mdlerror(mdlpvp->mdl_err_msg,mdlpvp);
  return(1);
}
}
	;

new_file_stream: VAR
{
  if (mdlpvp->cval_2!=NULL) {
    mdlpvp->sym_name=mdlpvp->cval_2;
  }
  else {
    mdlpvp->sym_name=mdlpvp->cval;
  }

  if ((mdlpvp->gp=store_sym(mdlpvp->sym_name,FSTRM,volp->main_sym_table))==NULL) {
    sprintf(mdlpvp->mdl_err_msg,"%s %s","Cannot store file stream in table:",mdlpvp->sym_name);
    mdlerror(mdlpvp->mdl_err_msg,mdlpvp);
    if (mdlpvp->sym_name==mdlpvp->cval) {
      mdlpvp->cval=NULL;
    }
    else {
      mdlpvp->cval_2=NULL;
    }
    free((void *)mdlpvp->sym_name);
    return(1);
  }

  if (mdlpvp->sym_name==mdlpvp->cval) {
    mdlpvp->cval=NULL;
  }
  else {
    mdlpvp->cval_2=NULL;
  }

  $$=mdlpvp->gp;
}
	;

file_mode: str_expr
{
a_str=$<str>1;
c=a_str[0];
if (c!='r' 
    && c!='w'
    && c!='a') {
  sprintf(mdlpvp->mdl_err_msg,"%s %s","Invalid file mode:",a_str);
  mdlerror(mdlpvp->mdl_err_msg,mdlpvp);
  free((void *)mdlpvp->cval);
  mdlpvp->cval=NULL;
  return(1);
}
$$=a_str;
}
	;

fclose_stmt: FCLOSE '(' existing_file_stream ')'
{
mdlpvp->gp=$<sym>3;
filep=(struct file_stream *)mdlpvp->gp->value;
if (fclose(filep->stream)!=0) {
  sprintf(mdlpvp->mdl_err_msg,"%s %s","Error closing file:",filep->name);
  mdlerror(mdlpvp->mdl_err_msg,mdlpvp);
  return(1);
}
}
	;

existing_file_stream: VAR
{
  if (mdlpvp->cval_2!=NULL) {
    mdlpvp->sym_name=mdlpvp->cval_2;
  }
  else {
    mdlpvp->sym_name=mdlpvp->cval;
  }

  if ((mdlpvp->gp=retrieve_sym(mdlpvp->sym_name,FSTRM,volp->main_sym_table))==NULL) {
    sprintf(mdlpvp->mdl_err_msg,"%s %s","Undefined file stream:",mdlpvp->sym_name);
    mdlerror(mdlpvp->mdl_err_msg,mdlpvp);
    if (mdlpvp->sym_name==mdlpvp->cval) {
      mdlpvp->cval=NULL;
    }
    else {
      mdlpvp->cval_2=NULL;
    }
    free((void *)mdlpvp->sym_name);
    return(1);
  }

  if (mdlpvp->sym_name==mdlpvp->cval) {
    mdlpvp->cval=NULL;
  }
  else {
    mdlpvp->cval_2=NULL;
  }
  free((void *)mdlpvp->sym_name);

#ifdef KELP
  mdlpvp->gp->ref_count++;
  no_printf("ref_count: %d\n",mdlpvp->gp->ref_count);
#endif

  $$=mdlpvp->gp;
}
	;

printf_stmt: PRINTF arg_list_init '(' format_string ',' list_args ')'
{
a_str=$<str>4;
if (my_fprintf(stderr,a_str,arg_list)) {
  sprintf(mdlpvp->mdl_err_msg,"%s %s","Could not print to stderr:",a_str);
  mdlerror(mdlpvp->mdl_err_msg,mdlpvp);
  return(1);
}
}
	| PRINTF arg_list_init '(' format_string ')'
{
a_str=$<str>4;
if (my_fprintf(stderr,a_str,NULL)) {
  sprintf(mdlpvp->mdl_err_msg,"%s %s","Could not print to stderr:",a_str);
  mdlerror(mdlpvp->mdl_err_msg,mdlpvp);
  return(1);
}
};

arg_list_init: /* empty */
{
	num_args=0;
};

format_string: str_expr
{
a_str=$<str>1;
rem_str=a_str;
strcpy(fmt_str,"");
while(rem_str!=NULL) {
  pos1=strcspn(rem_str,"\\");
  if(pos1==strlen(rem_str)) {  /* no \ found */
    strcat(fmt_str,rem_str);
    rem_str=NULL;
  }
  else {  /* found a \ */
    strncat(fmt_str,rem_str,pos1);
    c=rem_str[pos1+1];
    switch(c) {
    case 'n':
      strcat(fmt_str,"\n");
      break;
    case 't':
      strcat(fmt_str,"\t");
      break;
    case '\\':
      strcat(fmt_str,"\\");
      break;
    case '\"':
      strcat(fmt_str,"\"");
      break;
    }
    rem_str=rem_str+pos1+2;
  }
}
$$=fmt_str;
}
	;

list_args: num_expr_only
{
        if ((arg_list[num_args].arg_value=(void *)double_dup($<dbl>1))
             ==NULL) {
          sprintf(mdlpvp->mdl_err_msg,"%s","Could not store argument");
          mdlerror(mdlpvp->mdl_err_msg,mdlpvp);
          return(1);
        }
        arg_list[num_args++].arg_type=DBL;
}
	| list_args ',' num_expr_only
{
        if ((arg_list[num_args].arg_value=(void *)double_dup($<dbl>3))
             ==NULL) {
          sprintf(mdlpvp->mdl_err_msg,"%s","Could not store argument");
          mdlerror(mdlpvp->mdl_err_msg,mdlpvp);
          return(1);
        }
        arg_list[num_args++].arg_type=DBL;
}
	| str_expr_only
{
        arg_list[num_args].arg_value=(void *)$<str>1;
        arg_list[num_args++].arg_type=STR;
}
	| list_args ',' str_expr_only
{
        arg_list[num_args].arg_value=(void *)$<str>3;
        arg_list[num_args++].arg_type=STR;
}
	| existing_var_only
{
  mdlpvp->gp=$<sym>1;
  arg_list[num_args].arg_value=mdlpvp->gp->value;
  switch (mdlpvp->gp->sym_type) {
  case DBL:
    arg_list[num_args++].arg_type=DBL;
    break;
  case STR:
    arg_list[num_args++].arg_type=STR;
    break;
  default:
    mdlerror("Invalid variable type referenced");
    return(1);
    break;
  }
}
	| list_args ',' existing_var_only
{
  mdlpvp->gp=$<sym>3;
  arg_list[num_args].arg_value=mdlpvp->gp->value;
  switch (mdlpvp->gp->sym_type) {
  case DBL:
    arg_list[num_args++].arg_type=DBL;
    break;
  case STR:
    arg_list[num_args++].arg_type=STR;
    break;
  default:
    mdlerror("Invalid variable type referenced");
    return(1);
    break;
  }
}
	;

fprintf_stmt: FPRINTF arg_list_init '(' existing_file_stream ',' format_string ',' list_args ')' 
{
  mdlpvp->gp=$<sym>4;
  filep=(struct file_stream *)mdlpvp->gp->value;
  a_str=$<str>6;
  if (my_fprintf(filep->stream,a_str,arg_list)) {
    sprintf(mdlpvp->mdl_err_msg,"%s %s","Could not print to file:",filep->name);
    mdlerror(mdlpvp->mdl_err_msg,mdlpvp);
    return(1);
  }
}
	| FPRINTF arg_list_init '(' existing_file_stream ',' format_string ')'
{
  mdlpvp->gp=$<sym>4;
  filep=(struct file_stream *)mdlpvp->gp->value;
  a_str=$<str>6;
  if (my_fprintf(filep->stream,a_str,NULL)) {
    sprintf(mdlpvp->mdl_err_msg,"%s %s","Could not print to file:",filep->name);
    mdlerror(mdlpvp->mdl_err_msg,mdlpvp);
    return(1);
  }
}
	;

sprintf_stmt: SPRINTF arg_list_init '(' assign_var ',' format_string ',' list_args ')'
{ 
  mdlpvp->gp=$<sym>4;
  a_str=$<str>6;
  if (my_sprintf(mdlpvp->str_buf2,a_str,arg_list)) {
    sprintf(mdlpvp->mdl_err_msg,"%s %s","Could not sprintf to variable:",mdlpvp->gp->name);
    mdlerror(mdlpvp->mdl_err_msg,mdlpvp);
    return(1);
  }
  mdlpvp->gp->sym_type=STR;
  mdlpvp->gp->value=(void *)my_strdup(mdlpvp->str_buf2);
}
        | SPRINTF arg_list_init '(' assign_var ',' format_string ')'
{
  mdlpvp->gp=$<sym>4;
  a_str=$<str>6;
  if (my_sprintf(mdlpvp->str_buf,a_str,NULL)) {
    sprintf(mdlpvp->mdl_err_msg,"%s %s","Could not sprintf to variable:",mdlpvp->gp->name);
    mdlerror(mdlpvp->mdl_err_msg,mdlpvp);
    return(1);
  }
  mdlpvp->gp->sym_type=STR;
  mdlpvp->gp->value=(void *)my_strdup(mdlpvp->str_buf);
};

print_time_stmt: PRINT_TIME '(' format_string ')'
{
  a_str=$<str>3;
  the_time=time(NULL);
  strftime(time_str,128,a_str,localtime(&the_time));
  if (procnum == 0) fprintf(stderr,"%s",time_str);
};


fprint_time_stmt: FPRINT_TIME '(' existing_file_stream ',' format_string ')'
{
  mdlpvp->gp=$<sym>3;
  filep=(struct file_stream *)mdlpvp->gp->value;
  a_str=$<str>5;
  the_time=time(NULL);
  strftime(time_str,128,a_str,localtime(&the_time));
  if (fprintf(filep->stream,"%s",time_str)==EOF) {
    sprintf(mdlpvp->mdl_err_msg,"%s %s","Could not print to file:",filep->name);
    mdlerror(mdlpvp->mdl_err_msg,mdlpvp);
    return(1);
  }
};





#endif

/* ***************************************************************** */





/*  %%  */



/* Begin Bison Epilogue: */


void mdlerror(char *s,...)
{
  va_list ap;
  FILE *log_file;
  struct mdlparse_vars *mpvp;

  va_start(ap,s);
  mpvp=va_arg(ap,struct mdlparse_vars *);  
  va_end(ap);

  log_file=stderr;
  if (mpvp->vol->log_file!=NULL) {
    log_file=mpvp->vol->log_file;
  }

  if (mpvp->vol->procnum == 0) {
	fprintf(log_file,"MCell: error on line: %d of file: %s  %s\n",\
	        mpvp->line_num[mpvp->include_stack_ptr],mpvp->vol->curr_file,s);
	fflush(log_file);
  }
	return;
}


int mdlparse_init(struct volume *vol)
{
  struct mdlparse_vars *mpvp;
  FILE *mdl_infile;
  u_int i;

  if ((mpvp=(struct mdlparse_vars *)malloc
      (sizeof(struct mdlparse_vars)))==NULL) {
    fprintf(vol->log_file,"MCell: out of memory storing mdlparse vars\n");
    return(1);
  }

  mpvp->vol=vol;
  
  mpvp->cval=NULL;
  mpvp->cval_2=NULL;
  mpvp->strval=NULL;
  mpvp->ival=0;
  mpvp->rval=0;
  mpvp->gp=NULL;
  mpvp->tp=NULL;
  mpvp->sym_name=NULL;
  mpvp->val_1=0;
  mpvp->val_2=0;
  mpvp->tmp_dbl=0;
  mpvp->dblp=NULL;
  mpvp->num_pos=0;
  mpvp->elp=NULL;
  mpvp->el_head=NULL;
  mpvp->el_tail=NULL;
  mpvp->include_stack_ptr=0;
  mpvp->include_flag=0;
  mpvp->curr_obj=vol->root_object;
  mpvp->path_mem=NULL;
  mpvp->prod_mem=NULL;

  for(i=0;i<MAX_INCLUDE_DEPTH;i++) {
    mpvp->line_num[i]=1;
  }

  no_printf("Opening file %s\n",vol->mdl_infile_name);
  fflush(stderr);
  if ((mdl_infile=fopen(vol->mdl_infile_name,"r"))==NULL) {
    fprintf(vol->log_file,"MCell: error opening file: %s\n",
      vol->mdl_infile_name);
    return(1);
  }

  mdlrestart(mdl_infile);
  if (mdlparse((void *)mpvp)) {
     return(1);
  }
  fclose(mdl_infile);

  free(mpvp);

  return(0);
}


