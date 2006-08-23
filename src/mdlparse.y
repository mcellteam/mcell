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
  #include "rng.h"
  #include "vector.h"
  #include "strfunc.h"
  #include "mcell_structs.h"
  #include "mem_util.h"
  #include "sym_table.h"
  #include "diffuse_util.h"
  #include "mdlparse_util.h"
  #include "mdlparse.h"
  #include "util.h"
  #include "react_output.h"

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


  #define min(x,y) ((x)<(y)) ? (x): (y)
  
  #define min3(x,y,z) ((x)<(y)) ? (((x)<(z))?(x):(z)) : (((y)<(z))?(y):(z))


  int imin3(int f1, int f2, int f3)
  {
    return (min(f1,min(f2,f3)));
  }

  long long llmin3(long long f1,long long f2,long long f3)
  {
    return min(f1,min(f2,f3));
  }


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
struct output_expression *cnt;
struct release_evaluator *rev;
} 


%{
  #include "mdllex.flex.c"
%}


%pure_parser

%name-prefix="mdl"
%output="mdlparse.bison.c"


%token <tok> ABS
%token <tok> ABSORPTIVE
%token <tok> ACCURATE_3D_REACTIONS
%token <tok> ACOS
%token <tok> ALL_DATA
%token <tok> ALL_CROSSINGS
%token <tok> ALL_ELEMENTS
%token <tok> ALL_ENCLOSED
%token <tok> ALL_EVENTS
%token <tok> ALL_HITS
%token <tok> ALL_ITERATIONS
%token <tok> ALL_MESHES
%token <tok> ALL_MOLECULES
%token <tok> ALL_NOTIFICATIONS
%token <tok> ALL_TIMES
%token <tok> ALL_WARNINGS
%token <tok> AREA_OCCUPIED
%token <tok> ASCII
%token <tok> ASIN
%token <tok> ASPECT_RATIO
%token <tok> ATAN
%token <tok> BACK
%token <tok> BACK_CROSSINGS
%token <tok> BACK_HITS
%token <tok> BINDING_POLE
%token <tok> BOTH_POLES
%token <tok> BOTTOM
%token <tok> BOX
%token <tok> BOX_TRIANGULATION_REPORT
%token <tok> BRIEF
%token <tok> CEIL
%token <tok> CENTER_MOLECULES_ON_GRID
%token <tok> CHARGE
%token <tok> CHECKPOINT_INFILE
%token <tok> CHECKPOINT_OUTFILE
%token <tok> CHECKPOINT_ITERATIONS
%token <tok> CLAMP_CONCENTRATION
%token <tok> CLOSE_PARTITION_SPACING
%token <tok> COLOR
%token <tok> COLOR_SIDE
%token <tok> COLOR_EFFECTOR
%token <tok> COMPARTMENT
%token <tok> CONCENTRATION
%token <tok> CORNERS
%token <tok> COS
%token <tok> COUNT
%token <tok> CUSTOM_SPACE_STEP
%token <tok> CUSTOM_RK
%token <tok> CUSTOM_TIME_STEP
%token <tok> CUBIC
%token <tok> CUBIC_RELEASE_SITE
%token <tok> CUMULATE_FOR_EACH_TIME_STEP
%token <tok> DIFFUSION_CONSTANT_2D
%token <tok> DIFFUSION_CONSTANT_3D
%token <tok> DEFINE
%token <tok> DEFINE_EFFECTOR_SITE_POSITIONS
%token <tok> DEFINE_MOLECULE
%token <tok> DEFINE_MOLECULES
%token <tok> DEFINE_REACTIONS
%token <tok> DEFINE_RELEASE_PATTERN
%token <tok> DEFINE_SURFACE_REGIONS
%token <tok> DEFINE_SURFACE_CLASS
%token <tok> DEFINE_SURFACE_CLASSES
%token <tok> DEGENERATE_POLYGONS
%token <tok> DELAY
%token <tok> DENSITY
%token <tok> DIFFUSION_CONSTANT_REPORT
%token <tok> DX
%token <tok> DREAMM_V3
%token <tok> DREAMM_V3_GROUPED
%token <tok> EFFECTOR
%token <tok> EFFECTOR_GRID_DENSITY
%token <tok> EFFECTOR_POSITIONS
%token <tok> EFFECTOR_STATE 
%token <tok> EFFECTOR_STATES
%token <tok> EITHER_POLE
%token <tok> ELEMENT
%token <tok> ELEMENT_CONNECTIONS
%token <tok> ELLIPTIC
%token <tok> ELLIPTIC_RELEASE_SITE
%token <tok> EOF_TOK
%token <tok> ERROR
%token <tok> ESTIMATE_CONCENTRATION
%token <tok> EXCLUDE_ELEMENTS
%token <tok> EXCLUDE_PATCH
%token <tok> EXCLUDE_REGION
%token <tok> EXP
%token <tok> EXPRESSION
%token <tok> FALSE
%token <tok> FCLOSE
%token <tok> FILENAME
%token <tok> FILE_OUTPUT_REPORT
%token <tok> FINAL_SUMMARY
%token <tok> FLOOR
%token <tok> FOPEN
%token <tok> FOR_EACH_EFFECTOR
%token <tok> FOR_EACH_MOLECULE
%token <tok> FOR_EACH_TIME_STEP
%token <tok> FPRINT_TIME
%token <tok> FPRINTF
%token <tok> FRAME_DATA
%token <tok> FRONT
%token <tok> FRONT_CROSSINGS
%token <tok> FRONT_HITS
%token <tok> GAUSSIAN_RELEASE_NUMBER
%token <tok> GEOMETRY
%token <tok> HEADER
%token <tok> HIGH_PROBABILITY_THRESHOLD
%token <tok> HIGH_REACTION_PROBABILITY
%token <tok> IGNORED
%token <tok> INCLUDE_ELEMENTS
%token <tok> INCLUDE_FILE
%token <tok> INCLUDE_MOLECULES
%token <tok> INCLUDE_OBJECT
%token <tok> INCLUDE_PATCH
%token <tok> INCLUDE_REGION
%token <tok> INITIAL_EVENTS
%token <tok> INPUT_FILE
%token <tok> INSTANTIATE
%token <tok> INTEGER
%token <tok> INTERACTION_RADIUS
%token <tok> INTERIM_EVENTS
%token <tok> IRIT
%token <tok> ITERATION_FRAME_DATA
%token <tok> ITERATION_LIST
%token <tok> ITERATION_NUMBERS
%token <tok> ITERATION_REPORT
%token <tok> ITERATIONS
%token <tok> FULLY_RANDOM
%token <tok> LEFT
%token <tok> LIFETIME_TOO_SHORT
%token <tok> LIFETIME_THRESHOLD
%token <tok> LIST
%token <tok> LOCATION
%token <tok> LOG
%token <tok> LOG10
%token <tok> MAX_TOK
%token <tok> MCELL_GENERIC
%token <tok> MEAN_DIAMETER
%token <tok> MEAN_NUMBER
%token <tok> MESHES
%token <tok> MICRO_REV_SURF_RX
%token <tok> MICRO_REV_VOL_RX
%token <tok> MIN_TOK
%token <tok> MISSED_REACTIONS
%token <tok> MISSED_REACTION_THRESHOLD
%token <tok> MISSING_SURFACE_ORIENTATION
%token <tok> MOD
%token <tok> MODE
%token <tok> MODIFY_SURFACE_REGIONS
%token <tok> MOLECULE
%token <tok> MOLECULES
%token <tok> MOLECULE_DENSITY
%token <tok> MOLECULE_NUMBER
%token <tok> MOLECULE_POSITIONS
%token <tok> MOLECULE_STATES
%token <tok> MOLECULE_POSITIONS_STATES
%token <tok> MOLECULE_FILE_PREFIX
%token <tok> NAME
%token <tok> NAME_LIST
%token <tok> NEGATIVE_POLE
%token <tok> NEGATIVE_DIFFUSION_CONSTANT
%token <tok> NEGATIVE_REACTION_RATE
%token <tok> NO
%token <tok> NONE
%token <tok> NORMAL
%token <tok> NOTIFICATIONS
%token <tok> NO_SPECIES
%token <tok> NUMBER
%token <tok> NUMBER_BOUND
%token <tok> NUMBER_OF_TRAINS
%token <tok> NUMBER_TO_RELEASE
%token <tok> OBJECT
%token <tok> OBJECT_FILE_PREFIXES
%token <tok> OFF
%token <tok> ON
%token <tok> ORIENTATION
%token <tok> ORIENTATIONS
%token <tok> OUTPUT_BUFFER_SIZE
%token <tok> OVERWRITTEN_OUTPUT_FILE
%token <tok> POSITIONS
%token <tok> PARALLEL_PARTITION
%token <tok> PART
%token <tok> PARTITION_LOCATION_REPORT
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
%token <tok> PROBABILITY_REPORT
%token <tok> PROBABILITY_REPORT_THRESHOLD
%token <tok> PROGRESS_REPORT
%token <tok> RADIAL_DIRECTIONS
%token <tok> RADIAL_SUBDIVISIONS
%token <tok> RADIANCE
%token <tok> RAYSHADE
%token <tok> RAND_UNIFORM
%token <tok> RAND_GAUSSIAN
%token <tok> REACTION_DATA_OUTPUT
%token <tok> REACTION_GROUP
%token <tok> REAL
%token <tok> RECTANGULAR_RELEASE_SITE
%token <tok> RECTANGULAR_TOKEN
%token <tok> REFERENCE_STATE
%token <tok> REFLECTIVE
%token <tok> REFERENCE_DIFFUSION_CONSTANT
%token <tok> REGION_DATA
%token <tok> RELEASE_EVENT_REPORT
%token <tok> RELEASE_INTERVAL
%token <tok> RELEASE_PATTERN
%token <tok> RELEASE_PROBABILITY
%token <tok> RELEASE_SITE
%token <tok> REMOVE_ELEMENTS
%token <tok> RENDERMAN
%token <tok> RIGHT
%token <tok> ROTATE
%token <tok> ROUND_OFF
%token <tok> SCALE
%token <tok> SEED
%token <tok> SHAPE
%token <tok> SIN
%token <tok> SITE_DIAMETER
%token <tok> SITE_RADIUS
%token <tok> SPACE_STEP
%token <tok> SPECIFIED_EFFECTORS
%token <tok> SPECIFIED_MOLECULES
%token <tok> SPHERICAL
%token <tok> SPHERICAL_RELEASE_SITE
%token <tok> SPHERICAL_SHELL
%token <tok> SPHERICAL_SHELL_SITE
%token <tok> SPRINTF
%token <tok> SQRT
%token <tok> STANDARD_DEVIATION
%token <tok> STATE
%token <tok> STATE_VALUES
%token <tok> STEP
%token <tok> STR_VALUE
%token <tok> STRING_TO_NUM
%token <tok> SUMMATION_OPERATOR
%token <tok> SUM_OVER_ALL_EFFECTORS
%token <tok> SUM_OVER_ALL_MOLECULES
%token <tok> SUM_OVER_ALL_TIME_STEPS
%token <tok> SURFACE_CLASS
%token <tok> SURFACE_MOLECULE_AREA
%token <tok> SURFACE_POSITIONS
%token <tok> SURFACE_STATES
%token <tok> TAN
%token <tok> TARGET_ONLY
%token <tok> TET_ELEMENT_CONNECTIONS
%token <tok> TIME_LIST
%token <tok> TIME_POINTS
%token <tok> TIME_STEP
%token <tok> TIME_STEP_MAX
%token <tok> TO
%token <tok> TOP
%token <tok> TRAIN_DURATION
%token <tok> TRAIN_INTERVAL
%token <tok> TRANSFORM
%token <tok> TRANSLATE
%token <tok> TRANSPARENT
%token <tok> TRIGGER
%token <tok> TRUE
%token <tok> UNLIMITED
%token <tok> USELESS_VOLUME_ORIENTATION
%token <tok> VAR
%token <tok> VACANCY_SEARCH_DISTANCE
%token <tok> VARYING_PROBABILITY_REPORT
%token <tok> VERTEX_LIST
%token <tok> VIZ_DATA_OUTPUT
%token <tok> VIZ_OUTPUT
%token <tok> VIZ_VALUE
%token <tok> VOLUME_DEPENDENT_RELEASE_NUMBER
%token <tok> VOXEL_IMAGE_MODE
%token <tok> VOXEL_LIST
%token <tok> VOXEL_VOLUME_MODE
%token <tok> WARNING
%token <tok> WARNINGS
%token <tok> WILDCARD_VAR
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
%type <tok> space_def
%type <tok> iteration_def
%type <tok> grid_density_def
%type <tok> interact_radius_def
%type <tok> radial_directions_def
%type <tok> radial_subdivisions_def
%type <tok> optional_flag_def
%type <tok> optional_flag_randomgrid_def
%type <tok> optional_flag_expandlist_def
%type <tok> optional_vacancy_search_def
%type <tok> assignment_stmt
%type <tok> equals_or_to
%type <tok> partition_def
%type <tok> molecules_def
%type <tok> define_one_molecule
%type <tok> define_multiple_molecules
%type <tok> surface_classes_def
%type <tok> define_one_surface_class
%type <tok> define_multiple_surface_classes
%type <tok> surface_rxn_type
%type <tok> chkpt_stmt
%type <tok> release_pattern_def
%type <tok> site_size_cmd
%type <tok> physical_object_def
%type <tok> existing_obj_define_surface_regions
%type <tok> mod_surface_regions
%type <tok> partition_dimension
%type <tok> instance_def
%type <tok> include_stmt
%type <tok> viz_data_output_def
%type <tok> viz_output_def 
%type <tok> output_def
%type <tok> io_stmt
%type <tok> fopen_stmt
%type <tok> fclose_stmt
%type <tok> printf_stmt
%type <tok> fprintf_stmt
%type <tok> print_time_stmt
%type <tok> fprint_time_stmt
%type <tok> sprintf_stmt
%type <tok> notification_def
%type <tok> notification_list
%type <tok> notification_item_def
%type <tok> warnings_def
%type <tok> warning_list
%type <tok> warning_item_def
%type <tok> reversibility_def
%type <tok> surface_reversibility_def
%type <tok> volume_reversibility_def
%type <tok> end_of_mdl_file


%type <tok> rx_net_def
%type <tok> rx_stmt
%type <tok> list_rx_stmts
%type <tok> rxn
%type <tok> list_rxns
%type <tok> rx_group_def
%type <tok> product
%type <tok> list_products

%type <tok> notify_level
%type <tok> warning_level
%type <tok> boolean
%type <tok> release_site_geom_old
%type <tok> side_name
%type <tok> patch_type
%type <tok> prev_region_type
%type <tok> remove_side

%type <tok> iteration_frame_data_item
%type <tok> viz_molecules_one_item
%type <tok> viz_meshes_one_item
%type <tok> custom_header
%type <tok> file_arrow

%type <tok> hit_spec
%type <tok> opt_hit_spec

%type <sym> assign_var
%type <sym> existing_var_only
%type <sym> existing_str_var
%type <sym> existing_num_var
%type <sym> existing_array
%type <sym> existing_num_or_array

%type <dbl> num_value
%type <dbl> num_expr
%type <dbl> num_expr_only
%type <dbl> intOrReal 
%type <dbl> diffusion_def
%type <dbl> reference_diffusion_def
%type <dbl> mol_timestep_def
%type <dbl> target_def
%type <dbl> atomic_rate
%type <dbl> output_buffer_size_def

%type <str> str_value
%type <str> str_expr
%type <str> str_expr_only
%type <str> file_name
%type <str> opt_custom_header

%type <sym> new_molecule
%type <sym> existing_molecule
%type <sym> existing_surface_molecule
%type <sym> existing_molecule_opt_orient
%type <sym> new_surface_class
%type <sym> existing_surface_class
%type <sym> new_release_pattern
/* %type <sym> existing_release_pattern */
%type <sym> existing_release_pattern_xor_rxpn
%type <sym> object_def
%type <sym> new_object
%type <sym> existing_object
%type <sym> existing_one_or_multiple_objects
%type <sym> existing_one_or_multiple_molecules
%type <sym> object_ref
%type <sym> existing_object_ref
%type <sym> meta_object_def
%type <sym> release_site_def
%type <sym> release_site_def_new
%type <sym> release_site_def_old
%type <sym> box_def
%type <sym> polygon_list_def
%type <sym> voxel_list_def
%type <sym> new_region
%type <sym> existing_region
%type <sym> existing_logicalOrPhysical
%type <sym> new_rxn_pathname
%type <sym> existing_rxpn_or_molecule
%type <sym> count_location_specifier

%type <sym> new_file_stream
%type <sym> existing_file_stream

%type <str> file_mode
%type <str> format_string

%type <nel> array_value
%type <nel> array_expr_only

%type <vec3> point

%type <cnt> count_expr
%type <cnt> count_value
%type <cnt> count_syntax

%type <rev> release_region_expr

/**********************

%type <tok> two_way_unimolecular_rxn
%type <tok> two_way_bimolecular_rxn

%type <tok> r_spec
%type <tok> event_spec
%type <tok> orientation_cmd
%type <tok> wall_prop
%type <tok> orientation
%type <tok> lig_spec
%type <tok> t_spec
%type <tok> viz_mode_def
%type <tok> viz_output_def
%type <tok> list_viz_output_cmds
%type <tok> viz_output_cmd
%type <tok> viz_iteration_def 
%type <tok> viz_output_block_def
%type <tok> viz_time_def 
%type <tok> viz_object_prefixes_def
%type <tok> list_viz_object_prefixes
%type <tok> viz_object_prefix
%type <tok> viz_molecule_prefix_def
%type <tok> viz_state_values_def
%type <tok> list_viz_state_values
%type <tok> viz_state_value
%type <tok> voxel_image_mode_def
%type <tok> voxel_volume_mode_def
%type <tok> effector_site_def
%type <tok> add_effector
%type <tok> wall_prop_cmd
%type <tok> polarity
%type <tok> polarity_spec
%type <tok> viz_frame_data_def
%type <tok> list_frame_data_specs
%type <tok> frame_data_spec
%type <tok> parallel_partition_def
%type <tok> partition_plane

%type <dbl> rate
%type <dbl> effector_quantity_cmd

%type <sym> transition
%type <sym> existing_molecule_or_reaction_state
%type <sym> reaction_state
%type <sym> existing_reaction_state



%type <evnt> event 


***************************************/



%right '='
%left '&'
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
	| space_def
	| iteration_def
	| optional_flag_def
	| grid_density_def
        | interact_radius_def
	| radial_directions_def
	| radial_subdivisions_def
	| partition_def
	| assignment_stmt
	| molecules_def
	| surface_classes_def
	| rx_net_def
	| chkpt_stmt
	| release_pattern_def
	| physical_object_def
	| instance_def
	| existing_obj_define_surface_regions
	| mod_surface_regions
	| io_stmt
	| viz_data_output_def
	| viz_output_def
	| output_def
        | notification_def
        | warnings_def
        | reversibility_def

/*	
	| parallel_partition_def
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
             if (prepare_reactions(mdlpvp))
	     {
	       mdlerror("Failed to initialize reactions",mdlpvp);
	       return(1);
	     }
/*
            partition_volume(volume);
            build_ligand_table();
            if (volp->procnum == 0) print_rx();
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
  mdlpvp->line_num[mdlpvp->include_stack_ptr] = 1;
  volp->curr_file=my_strdup(mdlpvp->a_str);
  if(volp->curr_file == NULL){
    sprintf(mdlpvp->mdl_err_msg,"%s %s","Out of memory while switching to include file:",mdlpvp->a_str);
    mdlerror(mdlpvp->mdl_err_msg,mdlpvp);
    return(1);
  }
  yyless(0);
  yyclearin;
  mdl_switch_to_buffer(mdl_create_buffer(mdlin,YY_BUF_SIZE));
  no_printf("Now parsing MDL file: %s\n",volp->curr_file);
  fflush(stderr);
  mdlpvp->include_flag = 0;
};

reversibility_def: surface_reversibility_def | volume_reversibility_def;

surface_reversibility_def: MICRO_REV_SURF_RX '=' boolean
{
  mdlpvp->vol->surface_reversibility=$<tok>3;
};

volume_reversibility_def: MICRO_REV_VOL_RX '=' boolean
{
  mdlpvp->vol->volume_reversibility=$<tok>3;
};

notification_def: NOTIFICATIONS '{' notification_list '}';

notification_list: notification_item_def
        | notification_list notification_item_def;

notification_item_def:
        ALL_NOTIFICATIONS '=' boolean
{
  byte notify_value;
  if ($<tok>3==1) notify_value = NOTIFY_FULL;
  else notify_value = NOTIFY_NONE;
  
  mdlpvp->vol->notify->progress_report = notify_value;
  mdlpvp->vol->notify->diffusion_constants = notify_value;
  mdlpvp->vol->notify->reaction_probabilities = notify_value;
  mdlpvp->vol->notify->time_varying_reactions = notify_value;
  mdlpvp->vol->notify->partition_location = notify_value;
  mdlpvp->vol->notify->box_triangulation = notify_value;
  mdlpvp->vol->notify->custom_iterations = notify_value;
  mdlpvp->vol->notify->release_events = notify_value;
  mdlpvp->vol->notify->file_writes = notify_value;
  mdlpvp->vol->notify->final_summary = notify_value;
}
        | PROGRESS_REPORT '=' boolean
{
  if ($<tok>3==1) mdlpvp->vol->notify->progress_report = NOTIFY_FULL;
  else mdlpvp->vol->notify->progress_report = NOTIFY_NONE;
}
        | DIFFUSION_CONSTANT_REPORT '=' notify_level
{
  mdlpvp->vol->notify->diffusion_constants = $<tok>3;
}
        | PROBABILITY_REPORT '=' boolean
{
  if ($<tok>3==1) mdlpvp->vol->notify->reaction_probabilities = NOTIFY_FULL;
  else mdlpvp->vol->notify->reaction_probabilities = NOTIFY_NONE;
}
        | VARYING_PROBABILITY_REPORT '=' boolean
{
  if ($<tok>3==1) mdlpvp->vol->notify->time_varying_reactions = NOTIFY_FULL;
  else mdlpvp->vol->notify->time_varying_reactions = NOTIFY_NONE;
}
        | PROBABILITY_REPORT_THRESHOLD '=' num_expr
{
  mdlpvp->vol->notify->reaction_prob_notify = $<dbl>3;
}
        | PARTITION_LOCATION_REPORT '=' boolean
{
  if ($<tok>3==1) mdlpvp->vol->notify->partition_location = NOTIFY_FULL;
  else mdlpvp->vol->notify->partition_location = NOTIFY_NONE;
}
        | BOX_TRIANGULATION_REPORT '=' boolean
{
  if ($<tok>3==1) mdlpvp->vol->notify->box_triangulation = NOTIFY_FULL;
  else mdlpvp->vol->notify->box_triangulation = NOTIFY_NONE;
}
        | RELEASE_EVENT_REPORT '=' boolean
{
  if ($<tok>3==1) mdlpvp->vol->notify->release_events = NOTIFY_FULL;
  else mdlpvp->vol->notify->release_events = NOTIFY_NONE;
}
        | FILE_OUTPUT_REPORT '=' boolean
{
  if ($<tok>3==1) mdlpvp->vol->notify->file_writes = NOTIFY_FULL;
  else mdlpvp->vol->notify->file_writes = NOTIFY_NONE;
}
        | FINAL_SUMMARY '=' boolean
{
  if ($<tok>3==1) mdlpvp->vol->notify->final_summary = NOTIFY_FULL;
  else mdlpvp->vol->notify->final_summary = NOTIFY_NONE;
}
        | ITERATION_REPORT '=' boolean
{
  if (mdlpvp->vol->log_freq == -1) /* Not set on command line */
  {
    if ($<tok>3==1) mdlpvp->vol->notify->custom_iterations = NOTIFY_FULL;
    else mdlpvp->vol->notify->custom_iterations = NOTIFY_NONE;
  }
}
        | ITERATION_REPORT '=' num_expr
{
  if (mdlpvp->vol->log_freq == -1) /* Not set on command line */
  {
    mdlpvp->vol->notify->custom_iterations = NOTIFY_CUSTOM;
    if ($<dbl>3 < 1.0)
    {
      mdlerror("Invalid iteration number reporting interval: use value >= 1");
      return 1;
    }
    else mdlpvp->vol->notify->custom_iteration_value = (long long)$<dbl>3;
  }
};

notify_level: boolean
{
  if ($<tok>1==1) $$=NOTIFY_FULL;
  else $$=NOTIFY_NONE;
}
        | BRIEF
{
  $$=NOTIFY_BRIEF;
};


warnings_def: WARNINGS '{' warning_list '}';

warning_list: warning_item_def
        | warning_list warning_item_def;

warning_item_def:
        ALL_WARNINGS '=' warning_level
{
  byte warn_value = (byte)$<tok>3;
  mdlpvp->vol->notify->neg_diffusion = warn_value;
  mdlpvp->vol->notify->neg_reaction = warn_value;
  mdlpvp->vol->notify->high_reaction_prob = warn_value;
  mdlpvp->vol->notify->close_partitions = warn_value;
  mdlpvp->vol->notify->degenerate_polys = warn_value;
  mdlpvp->vol->notify->overwritten_file = warn_value;
  
  if (warn_value==WARN_ERROR) warn_value=WARN_WARN;
  mdlpvp->vol->notify->short_lifetime = warn_value;
  mdlpvp->vol->notify->missed_reactions = warn_value;
  mdlpvp->vol->notify->missed_surf_orient = warn_value;
  mdlpvp->vol->notify->useless_vol_orient = warn_value;
}
        | NEGATIVE_DIFFUSION_CONSTANT '=' warning_level
{
  mdlpvp->vol->notify->neg_diffusion = (byte)$<tok>3;
}
        | NEGATIVE_REACTION_RATE '=' warning_level
{
  mdlpvp->vol->notify->neg_reaction = (byte)$<tok>3;
}
        | HIGH_REACTION_PROBABILITY '=' warning_level
{
  mdlpvp->vol->notify->high_reaction_prob = (byte)$<tok>3;
}
        |  HIGH_PROBABILITY_THRESHOLD '=' num_expr
{
  mdlpvp->vol->notify->reaction_prob_warn = $<dbl>3;
}
        | CLOSE_PARTITION_SPACING '=' warning_level
{
  mdlpvp->vol->notify->close_partitions = (byte)$<tok>3;
}
        | DEGENERATE_POLYGONS '=' warning_level
{
  mdlpvp->vol->notify->degenerate_polys = (byte)$<tok>3;
}
        | OVERWRITTEN_OUTPUT_FILE '=' warning_level
{
  mdlpvp->vol->notify->overwritten_file = (byte)$<tok>3;
}
        | LIFETIME_TOO_SHORT '=' warning_level
{
  mdlpvp->vol->notify->short_lifetime = (byte)$<tok>3;
}
        | LIFETIME_THRESHOLD '=' num_expr
{
  double lifetime = $<dbl>3;
  if (lifetime < 0.0)
  {
    mdlerror("Molecule lifetimes are measured in iterations and cannot be negative");
    return 1;
  }
  mdlpvp->vol->notify->short_lifetime_value = (long long)lifetime;
}
        | MISSED_REACTIONS '=' warning_level
{
  mdlpvp->vol->notify->missed_reactions = (byte)$<tok>3;
}
        | MISSED_REACTION_THRESHOLD '=' num_expr
{
  double rxfrac = $<dbl>3;
  if (rxfrac < 0.0 || rxfrac > 1.0)
  {
    mdlerror("Values for fraction of reactions missed should be between 0 and 1");
    return 1;
  }
  mdlpvp->vol->notify->missed_reaction_value = rxfrac;
}
	| MISSING_SURFACE_ORIENTATION '=' warning_level
{
  mdlpvp->vol->notify->missed_surf_orient = (byte)$<tok>3;
}
	| USELESS_VOLUME_ORIENTATION '=' warning_level
{
  mdlpvp->vol->notify->useless_vol_orient = (byte)$<tok>3;
};

warning_level:
        IGNORED
{
  $$=WARN_COPE;
}
        | WARNING
{
  $$=WARN_WARN;
}
        | ERROR
{
  $$=WARN_ERROR;
}


assignment_stmt: assign_var '=' num_expr_only
{
  mdlpvp->gp=$<sym>1;
  if ((mdlpvp->gp->value=(void *)malloc(sizeof(double)))==NULL) {
    sprintf(mdlpvp->mdl_err_msg,"%s %s","Out of memory while assigning variable:",mdlpvp->gp->name);
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
  if(mdlpvp->gp->value == NULL){
    sprintf(mdlpvp->mdl_err_msg,"%s %s","Out of memory while assigning string:",($<str>3));
    mdlerror(mdlpvp->mdl_err_msg,mdlpvp);
    return(1);
  }
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
      sprintf(mdlpvp->mdl_err_msg,"%s %s","Out of memory while assigning variable:",mdlpvp->gp->name);
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
    if(mdlpvp->gp->value == NULL){
    	sprintf(mdlpvp->mdl_err_msg,"%s %s",
          "Out of memory while assigning string:",(char *)mdlpvp->tp->value);
    	mdlerror(mdlpvp->mdl_err_msg,mdlpvp);
    	return(1);
    }
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


array_expr_only: '['
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
};


existing_num_or_array: VAR 
{
  if (mdlpvp->cval_2!=NULL) {
    mdlpvp->sym_name=mdlpvp->cval_2;
  }
  else {
    mdlpvp->sym_name=mdlpvp->cval;
  }
  if ((mdlpvp->gp=retrieve_sym(mdlpvp->sym_name,ARRAY,volp->main_sym_table))==NULL) {
    if ((mdlpvp->gp=retrieve_sym(mdlpvp->sym_name,DBL,volp->main_sym_table))==NULL) {
      sprintf(mdlpvp->mdl_err_msg,"%s %s","Undefined numeric or array variable:",mdlpvp->sym_name);
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
          sprintf(mdlpvp->mdl_err_msg,"%s %s","Out of memory while assigning variable:",mdlpvp->sym_name);
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
	| RAND_UNIFORM
{
  $$=rng_dbl(mdlpvp->vol->rng);
}
	| RAND_GAUSSIAN
{
  $$=rng_gauss(mdlpvp->vol->rng);
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
	| RAND_UNIFORM
{
  $$=rng_dbl(mdlpvp->vol->rng);
}
	| RAND_GAUSSIAN
{
  $$=rng_gauss(mdlpvp->vol->rng);
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
  if(my_strcat($<str>1,$<str>3) != NULL){
      $$=my_strcat($<str>1,$<str>3);
  }else{
    sprintf(mdlpvp->mdl_err_msg,"%s ","Out of memory while parsing string");
    mdlerror(mdlpvp->mdl_err_msg);
    return (1);
  }
};


str_value: STR_VALUE
{
  if(strip_quotes(mdlpvp->strval) != NULL){
  	$$=strip_quotes(mdlpvp->strval);
  }else{
    sprintf(mdlpvp->mdl_err_msg,"%s ","Out of memory while parsing string");
    mdlerror(mdlpvp->mdl_err_msg);
    return (1);
  }
  free(mdlpvp->strval);
}
        | WILDCARD_VAR
{
  if (mdlpvp->cval_2!=NULL) {  
    mdlpvp->strval=mdlpvp->cval_2;
  }   
  else {
    mdlpvp->strval=mdlpvp->cval;
  } 

  if(strip_quotes(mdlpvp->strval) != NULL){
      $$ = strip_quotes(mdlpvp->strval);
  }else{
      sprintf(mdlpvp->mdl_err_msg,"%s ","Out of memory while parsing string");
      mdlerror(mdlpvp->mdl_err_msg);
      free((void *)mdlpvp->strval);
    
      return (1);
  }
   if (mdlpvp->strval==mdlpvp->cval) {
       mdlpvp->cval=NULL;
   }
   else {
       mdlpvp->cval_2=NULL;
   }
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
  if(strip_quotes(mdlpvp->strval) != NULL){
      $$=strip_quotes(mdlpvp->strval);
  }else{
    sprintf(mdlpvp->mdl_err_msg,"%s ","Out of memory while parsing string");
    mdlerror(mdlpvp->mdl_err_msg);
    return (1);
  }
  free(mdlpvp->strval);
}
	| WILDCARD_VAR
{
  if (mdlpvp->cval_2!=NULL) {  
    mdlpvp->strval=mdlpvp->cval_2;
  }   
  else {
    mdlpvp->strval=mdlpvp->cval;
  } 

  if(strip_quotes(mdlpvp->strval) != NULL){
      $$ = strip_quotes(mdlpvp->strval);
  }else{
      sprintf(mdlpvp->mdl_err_msg,"%s ","Out of memory while parsing string");
      mdlerror(mdlpvp->mdl_err_msg);
      free((void *)mdlpvp->strval);
    
      return (1);
  }
   if (mdlpvp->strval==mdlpvp->cval) {
       mdlpvp->cval=NULL;
   }
   else {
       mdlpvp->cval_2=NULL;
   }
   free(mdlpvp->strval);

}
	| INPUT_FILE
{
  $$=volp->mdl_infile_name;
}
	| str_expr '&' str_expr
{
      $$=my_strcat($<str>1,$<str>3);
  if(my_strcat($<str>1,$<str>3) != NULL){
      $$=my_strcat($<str>1,$<str>3);
  }else{
    sprintf(mdlpvp->mdl_err_msg,"%s ","Out of memory while parsing string");
    mdlerror(mdlpvp->mdl_err_msg);
    return (1);
  }
};


array_expr:  num_expr
{
  if ((mdlpvp->elp=(struct num_expr_list *)malloc
      (sizeof(struct num_expr_list)))==NULL) {
    mdlerror("Out of memory while creating numerical array");
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
    mdlerror("Out of memory while creating numerical array");
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

space_def: SPACE_STEP '=' num_expr
{
volp->space_step=$<dbl>3;
if (volp->space_step<0) {
  sprintf(mdlpvp->mdl_err_msg,"Space step = %g\n\tSetting to %g\n",volp->space_step,-volp->space_step);
  mdl_warning(mdlpvp);
  volp->space_step = -volp->space_step;
}
no_printf("Space step = %g\n",volp->space_step);
volp->space_step *= 0.5*sqrt(MY_PI) / volp->length_unit; /* Use internal units, convert from mean to characterstic length */
fflush(stderr);
};

time_max_def: TIME_STEP_MAX '=' num_expr
{
volp->time_step_max = $<dbl>3;
if (volp->time_step_max<0) {
  sprintf(mdlpvp->mdl_err_msg,"Maximum time step = %g\n\tSetting to %g\n",volp->time_step_max,-volp->time_step_max);
  mdl_warning(mdlpvp);
  volp->time_step_max=-volp->time_step_max;
}
no_printf("Maximum time step = %g\n",volp->time_step_max);
fflush(stderr);
};

iteration_def: ITERATIONS '=' num_expr
{
if (volp->iterations==0) {
  volp->iterations=(long long) $<dbl>3;
}
no_printf("Iterations = %lld\n",volp->iterations);
fflush(stderr);
};

optional_flag_def: optional_flag_randomgrid_def 
	| optional_flag_expandlist_def 
	| optional_vacancy_search_def;

optional_flag_randomgrid_def: CENTER_MOLECULES_ON_GRID '=' boolean
{
  mdlpvp->vol->randomize_gmol_pos = !($<tok>3);
};

optional_flag_expandlist_def: ACCURATE_3D_REACTIONS '=' boolean
{
  mdlpvp->vol->use_expanded_list = $<tok>3;
};

optional_vacancy_search_def: VACANCY_SEARCH_DISTANCE '=' num_expr
{
  mdlpvp->vol->vacancy_search_dist2 = $<dbl>3;
  if (mdlpvp->vol->vacancy_search_dist2 < 0) mdlpvp->vol->vacancy_search_dist2=0;
};


radial_directions_def: RADIAL_DIRECTIONS '=' num_expr
{
  volp->radial_directions = (int) $<dbl>3;
  volp->num_directions=0;
  if (volp->d_step!=NULL) {
    free(volp->d_step);
  }
  if ((volp->d_step=init_d_step(volp->radial_directions,&volp->num_directions))==NULL) {
    mdlerror("Out of memory while creating d_step data for molecule");
    return(1);
  }

  /* Mask must contain at least every direction */
  for (volp->directions_mask = 1 ; volp->directions_mask < volp->num_directions ; volp->directions_mask <<= 1) {}
  if (volp->directions_mask > (1<<18))
  {
    mdlerror("Too many RADIAL_DIRECTIONS requested (max 131072).\n");
    return(1);
  }
  volp->directions_mask -= 1;

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
  if (volp->radial_subdivisions <= 0)
  {
    mdlerror("Must choose a positive number of radial subdivisions.");
    return 1;
  }

  if (volp->r_step!=NULL) free(volp->r_step);
  if (volp->r_step_surface!=NULL) free(volp->r_step_surface);

  volp->r_step = init_r_step(volp->radial_subdivisions);
  volp->r_step_surface = init_r_step_surface(volp->radial_subdivisions);
  
  if (volp->r_step==NULL || volp->r_step_surface==NULL)
  {
    mdlerror("Out of memory while creating r_step data for molecule");
    return(1);
  }
  
  no_printf("radial subdivisions = %d\n",volp->radial_subdivisions);
};


grid_density_def: EFFECTOR_GRID_DENSITY '=' num_expr
{
  volp->grid_density=$<dbl>3;
  no_printf("Max density = %f\n",volp->grid_density);
  
  volp->space_step*=volp->length_unit;
  volp->length_unit=1.0/sqrt(volp->grid_density);
  volp->space_step/=volp->length_unit;
  
  no_printf("Length unit = %f\n",volp->length_unit);
  mdlpvp->mc_factor=1.0e11*volp->grid_density*sqrt(MY_PI*volp->time_unit)/N_AV;
  mdlpvp->transport_mc_factor=6.2415e18*mdlpvp->mc_factor;
  fflush(stderr);
};

interact_radius_def: INTERACTION_RADIUS '=' num_expr
{
  volp->rx_radius_3d = $<dbl>3;
  no_printf("Molecule-molecule interaction radius = %f\n",volp->rx_radius_3d);
};


partition_def: partition_dimension '='
{
  mdlpvp->num_pos=0;
  mdlpvp->el_head=NULL;
  mdlpvp->el_tail=NULL;
}
        '[' list_range_specs ']' 
{
  int i;
  if ((mdlpvp->dblp=(double *)malloc((mdlpvp->num_pos+2)*sizeof(double)))==NULL) {
    mdlerror("Out of memory while creating volume partition data");
    return(1);
  }
  sort_num_expr_list(mdlpvp->el_head);
  i=1;
  mdlpvp->elp=mdlpvp->el_head;
  while(mdlpvp->elp!=NULL) {
    mdlpvp->dblp[i++]=mdlpvp->elp->value/volp->length_unit;
    mdlpvp->elp=mdlpvp->elp->next;
  }
  mdlpvp->dblp[0]=-GIGANTIC;
  mdlpvp->dblp[mdlpvp->num_pos+1]=GIGANTIC;
  switch ($<tok>1) {
  case X_PARTS:
    if (volp->x_partitions!=NULL) {
      free(volp->x_partitions);
    }
    volp->nx_parts=mdlpvp->num_pos+2;
    volp->x_partitions=mdlpvp->dblp;
    break;
  case Y_PARTS:
    if (volp->y_partitions!=NULL) {
      free(volp->y_partitions);
    }
    volp->ny_parts=mdlpvp->num_pos+2;
    volp->y_partitions=mdlpvp->dblp;
    break;
  case Z_PARTS:
    if (volp->z_partitions!=NULL) {
      free(volp->z_partitions);
    }
    volp->nz_parts=mdlpvp->num_pos+2;
    volp->z_partitions=mdlpvp->dblp;
    break;
  }
};

partition_dimension: PARTITION_X {$$=X_PARTS;} 
	| PARTITION_Y {$$=Y_PARTS;}
	| PARTITION_Z {$$=Z_PARTS;}
;

molecules_def: define_one_molecule
	| define_multiple_molecules
;


define_one_molecule: DEFINE_MOLECULE
{
  if (volp->notify->diffusion_constants==NOTIFY_BRIEF) fprintf(volp->log_file,"Defining molecule with the following diffusion constant:\n");
}
        molecule_stmt
{
  if (volp->notify->diffusion_constants==NOTIFY_BRIEF) fprintf(volp->log_file,"\n");
};


define_multiple_molecules: DEFINE_MOLECULES
{
  if (volp->notify->diffusion_constants==NOTIFY_BRIEF) fprintf(volp->log_file,"Defining molecules with the following diffusion constants:\n");
}
        '{'
	list_molecule_stmts
	'}'
{
  if (volp->notify->diffusion_constants==NOTIFY_BRIEF) fprintf(volp->log_file,"\n");
};


list_molecule_stmts: molecule_stmt
	| list_molecule_stmts molecule_stmt
;


molecule_stmt: new_molecule '{'
{
  mdlpvp->gp=$<sym>1;
  mdlpvp->specp=(struct species *)mdlpvp->gp->value;
  mdlpvp->specp->sym=mdlpvp->gp;
}
	reference_diffusion_def
	diffusion_def
        mol_timestep_def
	target_def
	'}'
{
  mdlpvp->specp->D_ref=$<dbl>4;
  mdlpvp->specp->D=$<dbl>5;
  mdlpvp->specp->time_step=$<dbl>6;
  if (volp->time_unit==0) {
    sprintf(mdlpvp->mdl_err_msg,"%s %s","TIME_STEP not yet specified.  Cannot define molecule:",mdlpvp->specp->sym->name);
    mdlerror(mdlpvp->mdl_err_msg,mdlpvp);
    return(1);
  }
  if (mdlpvp->specp->D_ref==0) {
    mdlpvp->specp->D_ref=mdlpvp->specp->D;
  }
  mdlpvp->mc_factor=1.0e11*volp->grid_density*sqrt(MY_PI*volp->time_unit)/N_AV;
  mdlpvp->transport_mc_factor=6.2415e18*mdlpvp->mc_factor;

  if (mdlpvp->specp->D==0) /* Immobile (boring) */
  {
    mdlpvp->specp->space_step=0.0;
    mdlpvp->specp->time_step=1.0;
  }
  else if (mdlpvp->specp->time_step != 0.0) /* Custom timestep */
  {
    if (mdlpvp->specp->time_step < 0) /* Hack--negative value means space step */
    {
      mdlpvp->specp->space_step = -mdlpvp->specp->time_step;
      mdlpvp->specp->time_step = (mdlpvp->specp->space_step*mdlpvp->specp->space_step)*MY_PI/(16.0 * 1.0e8 * mdlpvp->specp->D)/volp->time_unit;
      mdlpvp->specp->space_step /= volp->length_unit;
    }
    else
    {
      mdlpvp->specp->space_step = sqrt( 4.0 * 1.0e8 * mdlpvp->specp->D * mdlpvp->specp->time_step ) / volp->length_unit;
      mdlpvp->specp->time_step /= volp->time_unit;
    }
  }
  else if (volp->space_step==0) /* Global timestep */
  {
    mdlpvp->specp->space_step=sqrt(4.0*1.0e8*mdlpvp->specp->D*volp->time_unit)/volp->length_unit;
    mdlpvp->specp->time_step=1.0;
  }
  else /* Global spacestep */
  {
    mdlpvp->specp->space_step = volp->space_step;
    mdlpvp->specp->time_step = (volp->space_step*volp->space_step*volp->length_unit*volp->length_unit)*MY_PI/(16.0 * 1.0e8 * mdlpvp->specp->D)/volp->time_unit;
  }
  
  if ($<dbl>7 != 0.0)
  {
    mdlpvp->specp->flags |= CANT_INITIATE;
  }
  
  if (volp->r_step==NULL) {
    if ((volp->r_step=init_r_step(volp->radial_subdivisions))==NULL) {
      mdlerror("Out of memory while creating r_step data for molecule");
      return(1);
    }
  }
  if (volp->r_step_surface==NULL)
  {
    volp->r_step_surface = init_r_step_surface(volp->radial_subdivisions);
    if (volp->r_step_surface==NULL)
    {
      mdlerror("Cannot store r_step_surface data.");
      return 1;
    }
  }
  if (volp->d_step==NULL) {
    if ((volp->d_step=init_d_step(volp->radial_directions,&volp->num_directions))==NULL) {
      mdlerror("Out of memory while creating d_step data for molecule");
      return(1);
    }
    for (volp->directions_mask = 1 ; volp->directions_mask < volp->num_directions ; volp->directions_mask <<= 1) {}
    if (volp->directions_mask > (1<<18))
    {
      mdlerror("Internal error: bad number of default RADIAL_DIRECTIONS (max 131072).\n");
      return(1);
    }
    volp->directions_mask -= 1;
    
  }
  
  if (mdlpvp->specp->time_step==1.0)
  {
    mdlpvp->l_perp_bar=sqrt(4*1.0e8*mdlpvp->specp->D*volp->time_unit/MY_PI);
    mdlpvp->l_perp_rms=sqrt(2*1.0e8*mdlpvp->specp->D*volp->time_unit);
    mdlpvp->l_r_bar=2*mdlpvp->l_perp_bar;
    mdlpvp->l_r_rms=sqrt(6*1.0e8*mdlpvp->specp->D*volp->time_unit);
    if (volp->notify->diffusion_constants==NOTIFY_FULL)
    {
      if (volp->procnum == 0) {
        fprintf(volp->log_file,"\nMCell: Theoretical average diffusion distances for molecule %s:\n",mdlpvp->specp->sym->name);
        fprintf(volp->log_file,"\tl_r_bar = %.9g microns\n",mdlpvp->l_r_bar);
        fprintf(volp->log_file,"\tl_r_rms = %.9g microns\n",mdlpvp->l_r_rms);
        fprintf(volp->log_file,"\tl_perp_bar = %.9g microns\n",mdlpvp->l_perp_bar);
        fprintf(volp->log_file,"\tl_perp_rms = %.9g microns\n\n",mdlpvp->l_perp_rms);
      }
    }
    else if (volp->notify->diffusion_constants==NOTIFY_BRIEF)
    {
      if (volp->procnum==0) fprintf(volp->log_file,"  l_r_bar=%.9g um for %s\n",mdlpvp->l_r_bar,mdlpvp->specp->sym->name);
    }
  }
  else
  {
    if (volp->notify->diffusion_constants==NOTIFY_FULL)
    {
      if (volp->procnum == 0)
      {
        fprintf(volp->log_file,"\nMCell: Theoretical average diffusion time for molecule %s:\n",mdlpvp->specp->sym->name);
        fprintf(volp->log_file,"\tl_r_bar fixed at %.9g microns\n",volp->length_unit*mdlpvp->specp->space_step*2.0/sqrt(MY_PI));
        fprintf(volp->log_file,"\tPosition update every %.3e seconds (%.3g timesteps)\n\n",
                mdlpvp->specp->time_step*volp->time_unit,mdlpvp->specp->time_step);
      }
    }
    else if (volp->notify->diffusion_constants==NOTIFY_BRIEF)
    {
      if (volp->procnum==0) fprintf(volp->log_file,"  delta t=%.3g timesteps for %s\n",mdlpvp->specp->time_step,mdlpvp->specp->sym->name);
    }
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
    if ((mdlpvp->gp=retrieve_sym(mdlpvp->sym_name,RXPN,volp->main_sym_table))==NULL) {
      if ((mdlpvp->gp=store_sym(mdlpvp->sym_name,MOL,volp->main_sym_table))==NULL) {
        sprintf(mdlpvp->mdl_err_msg,"%s %s","Out of memory while creating molecule:",mdlpvp->sym_name);
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
      sprintf(mdlpvp->mdl_err_msg,"%s %s","Molecule already defined as a named reaction pathway:",mdlpvp->sym_name);
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


existing_one_or_multiple_molecules: VAR
{

  if (mdlpvp->cval_2!=NULL) {  
    mdlpvp->sym_name=mdlpvp->cval_2;
  }   
  else {
    mdlpvp->sym_name=mdlpvp->cval;
  } 

    /* here is just one molecule */
  
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

}
    	| WILDCARD_VAR
{
  int i;
  struct sym_table_list *stl;
  struct sym_table *sym_t;
  char *wildcard_string;

  if (mdlpvp->cval_2!=NULL) {  
    mdlpvp->sym_name=mdlpvp->cval_2;
  }   
  else {
    mdlpvp->sym_name=mdlpvp->cval;
  } 

  if(strip_quotes(mdlpvp->sym_name) != NULL){
      wildcard_string = strip_quotes(mdlpvp->sym_name);
  }else{
      sprintf(mdlpvp->mdl_err_msg,"%s ","Out of memory while parsing wildcard variable");
      mdlerror(mdlpvp->mdl_err_msg);
      free((void *)mdlpvp->sym_name);
    
      return (1);
  }

    /* here is a wildcard molecule name */
    /* do a full sym_table scan comparing each key */

    mdlpvp->sym_table_list_head = NULL;

    for(i = 0; i < SYM_HASHSIZE; i++)
    {
       for(sym_t = volp->main_sym_table[i]; sym_t != NULL; sym_t = sym_t->next)
       {
         if(is_wildcard_match(wildcard_string, sym_t->name)){ 

           if(sym_t->sym_type == MOL)
           {
               stl = (struct sym_table_list *)mem_get(mdlpvp->sym_list_mem);
                if(stl == NULL){
                  mdlerror("Out of memory while parsing wildcard variable", mdlpvp);
                  return 1;
                }
                stl->node = sym_t;
                stl->next = NULL;
                if(mdlpvp->sym_table_list_head == NULL){
                   mdlpvp->sym_table_list_head = stl;
                }else{
                   stl->next = mdlpvp->sym_table_list_head;
                   mdlpvp->sym_table_list_head = stl;
                }
            }
          }

         } /* end for */
    
    } /* end for */
   
  if (mdlpvp->sym_name==mdlpvp->cval) {
       mdlpvp->cval=NULL;
   }
   else {
       mdlpvp->cval_2=NULL;
   }
   free((void *)mdlpvp->sym_name);
   free(wildcard_string);
   
   mdlpvp->sym_table_list_head = sort_sym_list_by_name(mdlpvp->sym_table_list_head);
   $$=NULL;
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


diffusion_def: DIFFUSION_CONSTANT_3D '=' num_expr
{
  mdlpvp->specp->flags -= (mdlpvp->specp->flags & ON_GRID);
  if (volp->notify->neg_diffusion==WARN_COPE)
  {
    if ($<dbl>3 < 0) $$ = 0;
    else $$=$<dbl>3;
  }
  else if (volp->notify->neg_diffusion==WARN_WARN)
  {
    if ($<dbl>3 < 0)
    {
      mdlerror("Negative diffusion constant found, setting to zero and continuing.");
      $$ = 0;
    }
    else $$=$<dbl>3;
  }
  else
  {
    if ($<dbl>3 < 0)
    {
      mdlerror("Error: diffusion constants should be zero or positive.");
      return 1;
    }
    else $$=$<dbl>3;
  }
}
	| DIFFUSION_CONSTANT_2D '=' num_expr
{
  mdlpvp->specp->flags |= ON_GRID;
  if (volp->notify->neg_diffusion==WARN_COPE)
  {
    if ($<dbl>3 < 0) $$ = 0;
    else $$=$<dbl>3;
  }
  else if (volp->notify->neg_diffusion==WARN_WARN)
  {
    if ($<dbl>3 < 0)
    {
      mdlerror("Negative diffusion constant found, setting to zero and continuing.");
      $$ = 0;
    }
    else $$=$<dbl>3;
  }
  else
  {
    if ($<dbl>3 < 0)
    {
      mdlerror("Error: diffusion constants should be zero or positive.");
      return 1;
    }
    else $$=$<dbl>3;
  }
};


reference_diffusion_def: /* empty */
{
  $$=0;
}
	| REFERENCE_DIFFUSION_CONSTANT '=' num_expr
{
  $$=$<dbl>3;
};


mol_timestep_def: /* empty */
{
  $$=0.0;
}
	| CUSTOM_TIME_STEP '=' num_expr
{
  if ($<dbl>3 > 0)
  {
    $$=$<dbl>3;
  }
  else
  {
    mdlerror("Zero or negative custom time step is disallowed (ignoring).");
    $$ = 0.0;
  }
}
	| CUSTOM_SPACE_STEP '=' num_expr
{
  if ($<dbl>3 > 0)
  {
    $$ = - $<dbl>3;
  }
  else
  {
    mdlerror("Zero or negative custom space step is disallowed (ignoring).");
    $$ = 0.0;
  }
};

target_def: /* empty */
{
  $$ = 0.0;
}
	| TARGET_ONLY
{
  $$=1.0;
};



surface_classes_def: define_one_surface_class
	| define_multiple_surface_classes
;


define_one_surface_class: DEFINE_SURFACE_CLASS surface_class_stmt
;


define_multiple_surface_classes: DEFINE_SURFACE_CLASSES '{'
	list_surface_class_stmts
	'}'
;


list_surface_class_stmts: surface_class_stmt
	| list_surface_class_stmts surface_class_stmt
;


surface_class_stmt: new_surface_class '{'
{
  mdlpvp->stp1=$<sym>1;
  mdlpvp->specp=(struct species *)mdlpvp->gp->value;
  mdlpvp->eff_dat_head=NULL;
}
	list_surface_prop_stmts
	'}'
;


new_surface_class: new_molecule
{
  mdlpvp->gp=$<sym>1;
  mdlpvp->specp=(struct species *)mdlpvp->gp->value;
  mdlpvp->specp->sym=mdlpvp->gp;
  mdlpvp->specp->flags=IS_SURFACE;
  $$=mdlpvp->gp;
};


existing_surface_class: existing_molecule
{
  mdlpvp->gp=$<sym>1;
  mdlpvp->specp=(struct species *)mdlpvp->gp->value;
  mdlpvp->sym_name=mdlpvp->gp->name;
  if (mdlpvp->specp->flags!=IS_SURFACE) {
    sprintf(mdlpvp->mdl_err_msg,"%s %s","Undefined surface type:",mdlpvp->sym_name);
    mdlerror(mdlpvp->mdl_err_msg,mdlpvp);
    return(1);
  }
  $$=mdlpvp->gp;
};


list_surface_prop_stmts: surface_prop_stmt
	| list_surface_prop_stmts surface_prop_stmt
;


surface_prop_stmt: surface_rxn_stmt
	| surface_class_mol_stmt
	| surface_class_viz_value_stmt 
;


surface_rxn_stmt: surface_rxn_type equals_or_to existing_molecule_opt_orient
{
  mdlpvp->stp2=$<sym>3;
  mdlpvp->specp=(struct species *)mdlpvp->stp2->value;
  if (mdlpvp->specp->flags==IS_SURFACE) {
    sprintf(mdlpvp->mdl_err_msg,"%s %s -%s-> ...",
      "Illegal reaction between two surfaces in surface reaction:",
      mdlpvp->stp2->name,mdlpvp->stp1->name);
    mdlerror(mdlpvp->mdl_err_msg,mdlpvp);
    return(1);
  }
  mdlpvp->sym_name=concat_rx_name(mdlpvp->stp1->name,mdlpvp->stp2->name);
  if(mdlpvp->sym_name == NULL) {
    sprintf(mdlpvp->mdl_err_msg,"%s %s -%s-> ...",
      "Out of memory while parsing surface reaction:",mdlpvp->stp1->name,mdlpvp->stp2->name);
    mdlerror(mdlpvp->mdl_err_msg,mdlpvp);
    return(1);
  }
  if ((mdlpvp->stp3=retrieve_sym(mdlpvp->sym_name,RX,volp->main_sym_table))
      !=NULL) {
  }
  else if ((mdlpvp->stp3=store_sym(mdlpvp->sym_name,RX,volp->main_sym_table))
      ==NULL) {
    sprintf(mdlpvp->mdl_err_msg,"%s %s -%s-> ...",
      "Out of memory while creating surface reaction:",mdlpvp->stp2->name,mdlpvp->stp1->name);
    mdlerror(mdlpvp->mdl_err_msg,mdlpvp);
    return(1);
  }
  if ((mdlpvp->pathp=(struct pathway *)mem_get(mdlpvp->path_mem))==NULL) {
    sprintf(mdlpvp->mdl_err_msg,"%s %s -%s-> ...",
      "Out of memory while creating surface reaction:",mdlpvp->stp2->name,mdlpvp->stp1->name);
    mdlerror(mdlpvp->mdl_err_msg,mdlpvp);
    return(1);
  }
  mdlpvp->rxnp=(struct rxn *)mdlpvp->stp3->value;
  mdlpvp->rxnp->n_reactants=2;
  mdlpvp->rxnp->n_pathways++;
  mdlpvp->pathp->pathname=NULL;
  mdlpvp->pathp->reactant1=(struct species *)mdlpvp->stp1->value;
  mdlpvp->pathp->reactant2=(struct species *)mdlpvp->stp2->value;
  mdlpvp->pathp->reactant3=NULL;
  mdlpvp->pathp->km=GIGANTIC;
  mdlpvp->pathp->kcat=0;
  mdlpvp->pathp->km_filename=NULL;

  if (mdlpvp->orient_class==0)
  {
    mdlpvp->pathp->orientation1=0;
    mdlpvp->pathp->orientation2=1;
    mdlpvp->pathp->orientation3=0;
  }
  else
  {
    mdlpvp->pathp->orientation1=1;
    mdlpvp->pathp->orientation2=(mdlpvp->orient_class<0) ? -1 : 1;
    mdlpvp->pathp->orientation3=0;
  }

  switch ($<tok>1) {
    case RFLCT:
      if ((mdlpvp->prodp=(struct product *)mem_get(mdlpvp->prod_mem))==NULL) {
        sprintf(mdlpvp->mdl_err_msg,"%s %s -%s-> ...",
          "Out of memory while creating surface reaction:",mdlpvp->stp2->name,mdlpvp->stp1->name);
        mdlerror(mdlpvp->mdl_err_msg,mdlpvp);
        return(1);
      }
      mdlpvp->pathp->kcat = KCAT_RATE_REFLECTIVE;
      mdlpvp->prodp->prod=mdlpvp->pathp->reactant2;
      mdlpvp->prodp->orientation=1;
      mdlpvp->prodp->next=NULL;
      mdlpvp->pathp->product_head=mdlpvp->prodp;
      break;
    case TRANSP:
      if ((mdlpvp->prodp=(struct product *)mem_get(mdlpvp->prod_mem))==NULL) {
        sprintf(mdlpvp->mdl_err_msg,"%s %s -%s-> ...",
          "Out of memory while creating surface reaction:",mdlpvp->stp2->name,mdlpvp->stp1->name);
        mdlerror(mdlpvp->mdl_err_msg,mdlpvp);
        return(1);
      }
      mdlpvp->pathp->kcat = KCAT_RATE_TRANSPARENT;
      mdlpvp->prodp->prod=mdlpvp->pathp->reactant2;
      mdlpvp->prodp->orientation=-1;
      mdlpvp->prodp->next=NULL;
      mdlpvp->pathp->product_head=mdlpvp->prodp;
      break;
    case SINK:
      mdlpvp->pathp->product_head=NULL;
      break;
    default:
      mdlerror("Unknown special surface type.");
      return 1;
      break;
  }

  mdlpvp->pathp->next=mdlpvp->rxnp->pathway_head;
  mdlpvp->rxnp->pathway_head=mdlpvp->pathp;

#ifdef DEBUG
  no_printf("Surface reaction defined:\n");
  no_printf("  %s[%d] -%s[%d]->",
  mdlpvp->rxnp->pathway_head->reactant2->sym->name,
  mdlpvp->rxnp->pathway_head->orientation2,
  mdlpvp->rxnp->pathway_head->reactant1->sym->name,
  mdlpvp->rxnp->pathway_head->orientation1);
  for (mdlpvp->prodp=mdlpvp->rxnp->pathway_head->product_head;
       mdlpvp->prodp!=NULL;mdlpvp->prodp=mdlpvp->prodp->next) {
    if (mdlpvp->prodp!=mdlpvp->rxnp->pathway_head->product_head) {
      no_printf(" +");
    }
    no_printf(" %s[%d]",mdlpvp->prodp->prod->sym->name,mdlpvp->prodp->orientation);
  }
  no_printf(" [%.9g,%.9g]\n",mdlpvp->rxnp->pathway_head->km,mdlpvp->rxnp->pathway_head->kcat);
#endif
}
        | CLAMP_CONCENTRATION existing_molecule_opt_orient '=' num_expr
{
  mdlpvp->stp2=$<sym>2;
  mdlpvp->specp=(struct species *)mdlpvp->stp2->value;
  if (mdlpvp->specp->flags==IS_SURFACE) {
    sprintf(mdlpvp->mdl_err_msg,"%s %s -%s-> ...",
      "Illegal reaction between two surfaces in surface reaction:",
      mdlpvp->stp2->name,mdlpvp->stp1->name);
    mdlerror(mdlpvp->mdl_err_msg,mdlpvp);
    return(1);
  }
  if (mdlpvp->specp->flags&ON_GRID)
  {
    mdlerror("Concentration clamp does not work on surface molecules.");
    return 1;
  }
  if (mdlpvp->specp->flags&NOT_FREE || mdlpvp->specp->D <= 0.0)
  {
    mdlerror("Concentration clamp must be applied to molecule diffusing in 3D");
    return 1;
   }
   if (($<dbl>4)<0)
   {
     mdlerror("Concentration can only be clamped to positive values.");
     return 1;
   }
    
  mdlpvp->sym_name=concat_rx_name(mdlpvp->stp1->name,mdlpvp->stp2->name);
  if(mdlpvp->sym_name == NULL) {
    sprintf(mdlpvp->mdl_err_msg,"%s %s -%s-> ...",
      "Memory allocation error:",mdlpvp->stp1->name,mdlpvp->stp2->name);
    mdlerror(mdlpvp->mdl_err_msg,mdlpvp);
    return(1);
  }
  if ((mdlpvp->stp3=retrieve_sym(mdlpvp->sym_name,RX,volp->main_sym_table))
      !=NULL) {
  }
  else if ((mdlpvp->stp3=store_sym(mdlpvp->sym_name,RX,volp->main_sym_table))
      ==NULL) {
    sprintf(mdlpvp->mdl_err_msg,"%s %s -%s-> ...",
      "Cannot store surface reaction:",mdlpvp->stp2->name,mdlpvp->stp1->name);
    mdlerror(mdlpvp->mdl_err_msg,mdlpvp);
    return(1);
  }
  if ((mdlpvp->pathp=(struct pathway *)mem_get(mdlpvp->path_mem))==NULL) {
    sprintf(mdlpvp->mdl_err_msg,"%s %s -%s-> ...",
      "Cannot store surface reaction:",mdlpvp->stp2->name,mdlpvp->stp1->name);
    mdlerror(mdlpvp->mdl_err_msg,mdlpvp);
    return(1);
  }
  mdlpvp->rxnp=(struct rxn *)mdlpvp->stp3->value;
  mdlpvp->rxnp->n_reactants=2;
  mdlpvp->rxnp->n_pathways++;
  mdlpvp->pathp->pathname=NULL;
  mdlpvp->pathp->reactant1=(struct species *)mdlpvp->stp1->value;
  mdlpvp->pathp->reactant2=(struct species *)mdlpvp->stp2->value;
  mdlpvp->pathp->reactant3=NULL;
  
  mdlpvp->pathp->km=GIGANTIC;
  
  /* Hack: kcat doesn't make sense for surfaces, so store CCN here! */
  mdlpvp->pathp->kcat=$<dbl>4;
  mdlpvp->pathp->km_filename=NULL;
  
  if (mdlpvp->orient_class==0)
  {
    mdlpvp->pathp->orientation1=0;
    mdlpvp->pathp->orientation2=1;
    mdlpvp->pathp->orientation3=0;
  }
  else
  {
    mdlpvp->pathp->orientation1=1;
    mdlpvp->pathp->orientation2=(mdlpvp->orient_class<0) ? -1 : 1;
    mdlpvp->pathp->orientation3=0;
  }
  mdlpvp->pathp->product_head=NULL;

  mdlpvp->pathp->next=mdlpvp->rxnp->pathway_head;
  mdlpvp->rxnp->pathway_head=mdlpvp->pathp;
};


surface_rxn_type: REFLECTIVE {$$=RFLCT;}
	| TRANSPARENT {$$=TRANSP;}
	| ABSORPTIVE {$$=SINK;}
;

equals_or_to: '=' {$$=TO} | TO;


surface_class_mol_stmt: surface_mol_stmt
{
  mdlpvp->specp=(struct species *)mdlpvp->stp1->value;
  mdlpvp->specp->eff_dat_head=mdlpvp->eff_dat_head;
};


surface_mol_stmt: mol_quant_type '{'
	list_surface_mol_quant
	'}'
;


mol_quant_type: MOLECULE_DENSITY
{
  mdlpvp->mol_quant_type=EFFDENS;
}
	| MOLECULE_NUMBER
{
  mdlpvp->mol_quant_type=EFFNUM;
};


list_surface_mol_quant: surface_mol_quant
	| list_surface_mol_quant surface_mol_quant
;


surface_mol_quant: existing_surface_molecule '=' num_expr
{
  mdlpvp->stp2=$<sym>1;
  mdlpvp->specp=(struct species *)mdlpvp->stp2->value;
  if ((mdlpvp->effdp=(struct eff_dat *)malloc(sizeof(struct eff_dat)))==NULL) {
    sprintf(mdlpvp->mdl_err_msg,"%s %s","Out of memory while storing data for surface molecule:",mdlpvp->stp2->name);
    return(1);
  }
  mdlpvp->effdp->next=mdlpvp->eff_dat_head;
  mdlpvp->eff_dat_head=mdlpvp->effdp;
  mdlpvp->effdp->eff=mdlpvp->specp;
  mdlpvp->effdp->quantity_type=mdlpvp->mol_quant_type;
  mdlpvp->effdp->quantity=$<dbl>3;
  mdlpvp->effdp->orientation=mdlpvp->orient_class;
  if (mdlpvp->orient_specified==0)
  {
    if (mdlpvp->vol->notify->missed_surf_orient==WARN_ERROR)
    {
      mdlerror("Error: surface orientation not specified for released surface molecule\n  (use ; or ', or ,' for random orientation)");
      return 1;
    }
    else if (mdlpvp->vol->notify->missed_surf_orient==WARN_WARN)
    {
      mdlerror("Warning: surface orientation not specified for released surface molecule\n  (use ; or ', or ,' for random orientation)");
    }
  }
};


surface_class_viz_value_stmt: VIZ_VALUE '=' num_expr
{
  mdlpvp->specp->region_viz_value = (int)$<dbl>3;
};


existing_surface_molecule: existing_molecule
{
  mdlpvp->orient_specified=0;
  mdlpvp->orient_class=0;
}
	orientation_class
{
  mdlpvp->stp2=$<sym>1;
  mdlpvp->specp=(struct species *)mdlpvp->stp2->value;
  if ((mdlpvp->specp->flags & ON_GRID) == 0)
  {
    sprintf(mdlpvp->mdl_err_msg,"%s %s","Invalid surface molecule specified:",mdlpvp->stp2->name);
    mdlerror(mdlpvp->mdl_err_msg,mdlpvp);
    return(1);
  }
  $$=mdlpvp->stp2;
};

existing_molecule_opt_orient: existing_molecule
{
  mdlpvp->orient_specified=0;
  mdlpvp->orient_class=0;
}
	orientation_class
{
  $$=$<sym>1;
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
  volp->chkpt_iterations=(long long) $<dbl>3;
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
  if (mdlpvp->rpatp->release_interval<=0)
  {
    mdlerror("Release interval must be set to a positive number.");
    return 1;
  }
  if (mdlpvp->rpatp->train_interval<=0)
  {
    mdlerror("Train interval must be set to a positive number.");
    return 1;
  }
  if (mdlpvp->rpatp->train_duration > mdlpvp->rpatp->train_interval)
  {
    mdlerror("Train duration must not be longer than the train interval.");
    return 1;
  }
  if (mdlpvp->rpatp->train_duration <= 0)
  {
    mdlerror("Train duration must be set to a positive number.");
  }
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
      sprintf(mdlpvp->mdl_err_msg,"%s %s","Out of memory while creating release pattern:",mdlpvp->sym_name);
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

/*
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
*/

existing_release_pattern_xor_rxpn: VAR
{
  struct sym_table *rpat;
  struct sym_table *rxpn;
  
  if (mdlpvp->cval_2!=NULL) mdlpvp->sym_name=mdlpvp->cval_2;
  else mdlpvp->sym_name=mdlpvp->cval;
  
  rpat=retrieve_sym(mdlpvp->sym_name,RPAT,volp->main_sym_table);
  rxpn=retrieve_sym(mdlpvp->sym_name,RXPN,volp->main_sym_table);
  
  if (mdlpvp->sym_name==mdlpvp->cval) mdlpvp->cval=NULL;
  else mdlpvp->cval_2=NULL;
  free(mdlpvp->sym_name);
  mdlpvp->sym_name=NULL;

  if (rpat==NULL && rxpn==NULL)
  {
    mdlerror("Cannot find named release pattern.");
    return 1;
  }
  else if (rpat!=NULL && rxpn!=NULL)
  {
    mdlerror("Named release pattern might be a pattern or a reaction pathway.  Please change one name.");
    return 1;
  }
  else
  {
    if (rpat!=NULL) mdlpvp->gp=rpat;
    else mdlpvp->gp=rxpn;
  }
  
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
        | box_def
{
        $$=$<sym>1;
}
        | polygon_list_def
{
        $$=$<sym>1;
}
        | voxel_list_def
{
        $$=$<sym>1;
};


meta_object_def:
        new_object OBJECT '{'
{
	mdlpvp->gp=$<sym>1;
        mdlpvp->objp=(struct object *)mdlpvp->gp->value;
        mdlpvp->objp->object_type=META_OBJ;
        mdlpvp->objp->parent=mdlpvp->curr_obj;
        mdlpvp->curr_obj=mdlpvp->objp;
}
        list_objects
        list_opt_object_cmds
        '}'
{
        mdlpvp->curr_obj->parent->n_walls+=mdlpvp->curr_obj->n_walls;
        mdlpvp->curr_obj->parent->n_walls_actual+=mdlpvp->curr_obj->n_walls_actual;
        mdlpvp->curr_obj->parent->n_verts+=mdlpvp->curr_obj->n_verts;
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
      sprintf(mdlpvp->mdl_err_msg,"%s %s","Out of memory while creating object:",mdlpvp->sym_name);
      mdlerror(mdlpvp->mdl_err_msg,mdlpvp);
      return(1);
    }
    mdlpvp->object_name_list->name=NULL;
    mdlpvp->object_name_list->prev=NULL;
    mdlpvp->object_name_list->next=NULL;
    mdlpvp->object_name_list_end=mdlpvp->object_name_list;
  }
  if ((mdlpvp->object_name_list_end=concat_obj_name(mdlpvp->object_name_list_end,mdlpvp->sym_name))==NULL) {
    sprintf(mdlpvp->mdl_err_msg,"%s %s","Out of memory while creating object:",mdlpvp->sym_name);
    return(1);
  }
  mdlpvp->obj_name=mdlpvp->object_name_list_end->name;
  if ((mdlpvp->gp=retrieve_sym(mdlpvp->obj_name,OBJ,volp->main_sym_table))==NULL) {
    if ((mdlpvp->gp=store_sym(mdlpvp->obj_name,OBJ,volp->main_sym_table))==NULL) {
      sprintf(mdlpvp->mdl_err_msg,"%s %s","Out of memory while creating object:",mdlpvp->obj_name);
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
  if(mdlpvp->prefix_name == NULL){
    sprintf(mdlpvp->mdl_err_msg,"%s %s","Out of memory while parsing object:",mdlpvp->sym_name);
    mdlerror(mdlpvp->mdl_err_msg,mdlpvp);
    return (1);
  }
  no_printf("found existing object %s\n",mdlpvp->objp->sym->name);
  no_printf("first name of existing object %s is %s\n",mdlpvp->sym_name,get_first_name(mdlpvp->sym_name));
  no_printf("prefix name of existing object %s is %s\n",mdlpvp->sym_name,mdlpvp->prefix_name);
  fflush(stderr);
  fflush(stderr);
#ifdef KELP
  mdlpvp->objp->sym->ref_count++;
  no_printf("ref_count: %d\n",mdlpvp->objp->sym->ref_count);
#endif
  $$=mdlpvp->objp->sym;
};

existing_one_or_multiple_objects: VAR
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
     if(mdlpvp->prefix_name == NULL){
        sprintf(mdlpvp->mdl_err_msg,"%s %s","Out of memory while parsing object:",mdlpvp->sym_name);
        mdlerror(mdlpvp->mdl_err_msg,mdlpvp);
        return (1);
     }
     no_printf("found existing object %s\n",mdlpvp->objp->sym->name);
     no_printf("first name of existing object %s is %s\n",mdlpvp->sym_name,get_first_name(mdlpvp->sym_name));
     no_printf("prefix name of existing object %s is %s\n",mdlpvp->sym_name,mdlpvp->prefix_name);
     fflush(stderr);
#ifdef KELP
     mdlpvp->objp->sym->ref_count++;
     no_printf("ref_count: %d\n",mdlpvp->objp->sym->ref_count);
#endif
     $$=mdlpvp->objp->sym;

     if (mdlpvp->sym_name==mdlpvp->cval) {
        mdlpvp->cval=NULL;
     }
     else {
       mdlpvp->cval_2=NULL;
     }
     fflush(stderr);
}
     | WILDCARD_VAR
{


  int i;
  struct sym_table *sym_t;
  struct sym_table_list *stl;
  char *wildcard_string;

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

  if(strip_quotes(mdlpvp->sym_name) != NULL){
      wildcard_string = strip_quotes(mdlpvp->sym_name);
  }else{
      sprintf(mdlpvp->mdl_err_msg,"%s ","Out of memory while parsing wildcard variable");
      mdlerror(mdlpvp->mdl_err_msg);
      free((void *)mdlpvp->sym_name);
    
      return (1);
  }

    /* here is a wildcard object */
    /* do a full sym_table scan comparing each key */

    mdlpvp->sym_table_list_head = NULL;

    for(i = 0; i < SYM_HASHSIZE; i++)
    {
       for(sym_t = volp->main_sym_table[i]; sym_t != NULL; sym_t = sym_t->next)
       {
          if(is_wildcard_match(wildcard_string, sym_t->name)){ 

           if(sym_t->sym_type == OBJ)
           {
               stl = (struct sym_table_list *)mem_get(mdlpvp->sym_list_mem);
                if(stl == NULL){
                  mdlerror("Out of memory while parsing wildcard variable", mdlpvp);
                  return 1;
                }
                stl->node = sym_t;
                stl->next = NULL;
                if(mdlpvp->sym_table_list_head == NULL){
                   mdlpvp->sym_table_list_head = stl;
                }else{
                   stl->next = mdlpvp->sym_table_list_head;
                   mdlpvp->sym_table_list_head = stl;
                }
            }
          }

         } /* end for */
    
    } /* end for */

     if (mdlpvp->sym_name==mdlpvp->cval) {
        mdlpvp->cval=NULL;
     }
     else {
       mdlpvp->cval_2=NULL;
     }
     free((void *)mdlpvp->sym_name);
     free(wildcard_string);
     fflush(stderr);
     
     mdlpvp->sym_table_list_head = sort_sym_list_by_name(mdlpvp->sym_table_list_head);
     
     $$=NULL; /* Means we need to look in list */
};

list_opt_object_cmds: /* empty */
	| list_opt_object_cmds opt_object_cmd
;


opt_object_cmd: transformation
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
  mdlpvp->curr_obj->parent->n_walls+=mdlpvp->curr_obj->n_walls;
  mdlpvp->curr_obj->parent->n_walls_actual+=mdlpvp->curr_obj->n_walls_actual;
  mdlpvp->curr_obj->parent->n_verts+=mdlpvp->curr_obj->n_verts;
  mdlpvp->curr_obj=mdlpvp->curr_obj->parent;
  if (mdlpvp->object_name_list_end->prev!=NULL) {
    mdlpvp->object_name_list_end=mdlpvp->object_name_list_end->prev;
  }
  else {
    mdlpvp->object_name_list_end->name=NULL;
  }
  $$=$<sym>1;
};


release_site_def: release_site_def_old
{
  $$ = $<sym>1;
}
	| release_site_def_new
{
  $$ = $<sym>1;
};


release_site_def_new: new_object RELEASE_SITE '{'
{
  mdlpvp->gp=$<sym>1;
  mdlpvp->objp=(struct object *)mdlpvp->gp->value;
  if ((mdlpvp->rsop=(struct release_site_obj *)malloc
              (sizeof(struct release_site_obj)))==NULL) {
    mdlerror("Out of memory while creating release site");
    return(1);
  }
  mdlpvp->rsop->location=NULL;
  mdlpvp->rsop->mol_type=NULL;
  mdlpvp->rsop->release_number_method=CONSTNUM;
  mdlpvp->rsop->release_shape = SHAPE_UNDEFINED;
  mdlpvp->rsop->orientation=0;
  mdlpvp->rsop->release_number=0;
  mdlpvp->rsop->mean_diameter=0;
  mdlpvp->rsop->concentration=0;
  mdlpvp->rsop->standard_deviation=0;
  mdlpvp->rsop->diameter=NULL;
  mdlpvp->rsop->region_data=NULL;
  mdlpvp->rsop->mol_list=NULL;
  mdlpvp->rsop->release_prob=1.0;
  mdlpvp->rsop->pattern=volp->default_release_pattern;

  mdlpvp->objp->object_type=REL_SITE_OBJ;
  mdlpvp->objp->contents=mdlpvp->rsop;
  mdlpvp->objp->parent=mdlpvp->curr_obj;

  mdlpvp->curr_obj=mdlpvp->objp;
}
	release_site_geom
	list_release_site_cmds
	list_opt_object_cmds
	'}'
{
  no_printf("Release site %s defined:\n",mdlpvp->curr_obj->sym->name);
  if (mdlpvp->rsop->release_shape!=SHAPE_LIST && mdlpvp->rsop->mol_type==NULL)
  {
    mdlerror("Must specify molecule to release using MOLECULE=molecule_name.");
    return 1;
  }
  if (mdlpvp->rsop->release_number_method==CCNNUM)
  {
    if ((mdlpvp->rsop->mol_type->flags&NOT_FREE)==0 && mdlpvp->rsop->release_number != -3)
    {
      mdlerror("CONCENTRATION must be used with molecules that can diffuse in 3D");
      if ((mdlpvp->rsop->mol_type->flags&NOT_FREE)==ON_GRID)
      {
	mdlerror("  (Use DENSITY for molecules diffusing in 2D.)");
      }
      return 1;
    }
    else if ((mdlpvp->rsop->mol_type->flags&NOT_FREE)==ON_GRID && mdlpvp->rsop->release_number != -2)
    {
      mdlerror("DENSITY must be used with molecules that can diffuse in 2D");
      if ((mdlpvp->rsop->mol_type->flags&NOT_FREE)==0)
      {
	mdlerror("  (Use CONCENTRATION for molecules diffusing in 3D.)");
      }
      return 1;
    }
  }
  if (mdlpvp->rsop->release_shape!=SHAPE_REGION)
  {
    if (mdlpvp->rsop->location==NULL)
    {
      if (mdlpvp->rsop->release_shape!=SHAPE_LIST || mdlpvp->rsop->mol_list==NULL)
      {
        mdlerror("Release site is missing location.\n");
        return 1;
      }
      else
      {
	mdlpvp->rsop->location = (struct vector3*)malloc(sizeof(struct vector3));
	if (mdlpvp->rsop->location==NULL)
	{
	  mdlerror("Out of memory storing region location");
	  return 1;
	}
	mdlpvp->rsop->location->x = 0;
	mdlpvp->rsop->location->y = 0;
	mdlpvp->rsop->location->z = 0;	
      }
    }
    no_printf("\tLocation = [%f,%f,%f]\n",mdlpvp->rsop->location->x,mdlpvp->rsop->location->y,mdlpvp->rsop->location->z);
    mdlpvp->rsop->location->x/=mdlpvp->vol->length_unit;
    mdlpvp->rsop->location->y/=mdlpvp->vol->length_unit;
    mdlpvp->rsop->location->z/=mdlpvp->vol->length_unit;
  }
  mdlpvp->curr_obj=mdlpvp->curr_obj->parent;
  if (mdlpvp->object_name_list_end->prev!=NULL) {
    mdlpvp->object_name_list_end=mdlpvp->object_name_list_end->prev;
  }
  else {
    mdlpvp->object_name_list_end->name=NULL;
  }
  $$=$<sym>1;
};

release_site_geom: SHAPE '=' release_region_expr
{
  struct release_evaluator *re = $<rev>3;
  struct release_region_data *rrd;
  
  mdlpvp->objp = mdlpvp->curr_obj;
  
  mdlpvp->rsop->release_shape = SHAPE_REGION;
  mdlpvp->vol->releases_on_regions_flag = 1;
  
  rrd = (struct release_region_data*)malloc(sizeof(struct release_region_data));
  if (rrd==NULL)
  {
    mdlerror("Out of memory while trying to create release site on region");
    return 1;
  }
  
  rrd->n_walls_included = -1; /* Indicates uninitialized state */
  rrd->cum_area_list = NULL;
  rrd->wall_index = NULL;
  rrd->obj_index = NULL;
  rrd->n_objects = -1;
  rrd->owners = NULL;
  rrd->in_release = NULL;
  rrd->self = mdlpvp->curr_obj;
  
  rrd->expression = re;
  mdlpvp->rsop->region_data = rrd;
  
  if (check_release_regions(re,mdlpvp->curr_obj,mdlpvp->vol->root_instance))
  {
    mdlerror("Trying to release on a region that the release site cannot see!\n  Try grouping the release site and the corresponding geometry with an OBJECT.");
    return 1;
  }
}
	| SHAPE '=' SPHERICAL
{
  mdlpvp->rsop->release_shape = SHAPE_SPHERICAL;
}
	| SHAPE '=' CUBIC
{
  mdlpvp->rsop->release_shape = SHAPE_CUBIC;
}
	| SHAPE '=' ELLIPTIC
{
  mdlpvp->rsop->release_shape = SHAPE_ELLIPTIC;
}
	| SHAPE '=' RECTANGULAR_TOKEN
{
  mdlpvp->rsop->release_shape = SHAPE_RECTANGULAR
}
	| SHAPE '=' SPHERICAL_SHELL
{
  mdlpvp->rsop->release_shape = SHAPE_SPHERICAL_SHELL;
}
	| SHAPE '=' LIST
{
  mdlpvp->rsop->release_shape = SHAPE_LIST;
  mdlpvp->rsop->release_number_method = CONSTNUM;
};

release_region_expr:
	existing_region
{
  struct sym_table *my_sym = $<sym>1;
  struct release_evaluator *re;
  
  mdlpvp->objp = mdlpvp->curr_obj;  /* Fix up vars from "existing_region" */
  
  re = (struct release_evaluator*)malloc(sizeof(struct release_evaluator));
  if (re==NULL)
  {
    mdlerror("Out of memory while trying to create release site on region");
    return 1;
  }
  
  re->op = REXP_NO_OP | REXP_LEFT_REGION;
  re->left = my_sym->value;
  re->right = NULL;
  
  ((struct region*)re->left)->flags |= COUNT_CONTENTS;
  
  $$=re;
}
	| '(' release_region_expr ')'
{
  $$ = $<rev>2;
}
	| release_region_expr '+' release_region_expr
{
  struct release_evaluator *re;
  re = pack_release_expr($<rev>1,$<rev>3,REXP_UNION);
  if (re==NULL)
  {
    mdlerror("Out of memory while trying to create release site on region");
    return 1;
  }
  $$ = re;
}
	| release_region_expr '-' release_region_expr
{
  struct release_evaluator *re;
  re = pack_release_expr($<rev>1,$<rev>3,REXP_SUBTRACTION);
  if (re==NULL)
  {
    mdlerror("Out of memory while trying to create release site on region");
    return 1;
  }
  $$ = re;
}
	| release_region_expr '*' release_region_expr
{
  struct release_evaluator *re;
  re = pack_release_expr($<rev>1,$<rev>3,REXP_INTERSECTION);
  if (re==NULL)
  {
    mdlerror("Out of memory while trying to create release site on region");
    return 1;
  }
  $$ = re;
};


release_site_def_old: new_object release_site_geom_old '{'
{
  mdlpvp->gp=$<sym>1;
  mdlpvp->objp=(struct object *)mdlpvp->gp->value;
  if ((mdlpvp->rsop=(struct release_site_obj *)malloc
              (sizeof(struct release_site_obj)))==NULL) {
    mdlerror("Out of memory while creating release site");
    return(1);
  }
  mdlpvp->rsop->location=NULL;
  mdlpvp->rsop->mol_type=NULL;
  mdlpvp->rsop->release_number_method=CONSTNUM;
  mdlpvp->rsop->release_shape = $<tok>2;
  mdlpvp->rsop->release_number=0;
  mdlpvp->rsop->mean_diameter=0;
  mdlpvp->rsop->concentration=0;
  mdlpvp->rsop->standard_deviation=0;
  mdlpvp->rsop->diameter=NULL;
  mdlpvp->rsop->region_data=NULL;
  mdlpvp->rsop->mol_list=NULL;
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
  mdlpvp->rsop->location->x/=mdlpvp->vol->length_unit;
  mdlpvp->rsop->location->y/=mdlpvp->vol->length_unit;
  mdlpvp->rsop->location->z/=mdlpvp->vol->length_unit;

  mdlpvp->curr_obj=mdlpvp->curr_obj->parent;
  if (mdlpvp->object_name_list_end->prev!=NULL) {
    mdlpvp->object_name_list_end=mdlpvp->object_name_list_end->prev;
  }
  else {
    mdlpvp->object_name_list_end->name=NULL;
  }
  $$=$<sym>1;
};


release_site_geom_old: SPHERICAL_RELEASE_SITE
{
  $$=SHAPE_SPHERICAL;
}
	| CUBIC_RELEASE_SITE
{
  $$=SHAPE_CUBIC;
}
	| ELLIPTIC_RELEASE_SITE
{
  $$=SHAPE_ELLIPTIC;
}
	| RECTANGULAR_RELEASE_SITE
{
  $$=SHAPE_RECTANGULAR;
}
	| SPHERICAL_SHELL_SITE
{
  $$=SHAPE_SPHERICAL_SHELL;
};


list_release_site_cmds: release_site_cmd
	| list_release_site_cmds release_site_cmd
;


release_site_cmd:
	LOCATION '=' point
{
  mdlpvp->rsop->location=$<vec3>3;
}
	| MOLECULE '=' existing_molecule_opt_orient
{
  mdlpvp->gp=$<sym>3;
  mdlpvp->rsop->mol_type=(struct species *)mdlpvp->gp->value;
  if ((mdlpvp->rsop->mol_type->flags&NOT_FREE)==0 &&
      mdlpvp->rsop->release_shape==SHAPE_REGION)
  {
    mdlpvp->vol->place_waypoints_flag=1;
  }
  mdlpvp->rsop->orientation = mdlpvp->orient_class;
  if (mdlpvp->rsop->mol_type->flags&ON_GRID)
  {
    if (!mdlpvp->orient_specified)
    {
      if (mdlpvp->vol->notify->missed_surf_orient==WARN_ERROR)
      {
	mdlerror("Error: surface orientation not specified for released surface molecule\n  (use ; or ', or ,' for random orientation)");
	return 1;
      }
      else if (mdlpvp->vol->notify->missed_surf_orient==WARN_WARN)
      {
	mdlerror("Warning: surface orientation not specified for released surface molecule\n  (use ; or ', or ,' for random orientation)");
      }
    }
  }
  else if ((mdlpvp->rsop->mol_type->flags&NOT_FREE)==0)
  {
    if (mdlpvp->orient_specified)
    {
      if (mdlpvp->vol->notify->useless_vol_orient==WARN_ERROR)
      {
	mdlerror("Error: orientation not used for released volume molecule");
	return 1;
      }
      else if (mdlpvp->vol->notify->useless_vol_orient==WARN_WARN)
      {
	mdlerror("Warning: orientation not used for released volume molecule");
      }
    }
  }
  else
  {
    mdlerror("Error: cannot release a surface class instead of a molecule.");
    return 1;
  }
}
	| release_number_cmd
{
  if (mdlpvp->rsop->release_shape==SHAPE_LIST)
  {
    mdlerror("Molecules are already specified in a list--cannot set number or density.");
    return 1;
  }
}
	| site_size_cmd '=' num_expr_only
{
  double scaling_factor = 1.0/volp->length_unit;
  if ($<tok>1==SITE_RADIUS) scaling_factor*=2;
  
  if ((mdlpvp->rsop->diameter=(struct vector3 *)malloc(sizeof(struct vector3)))==NULL) {
    mdlerror("Out of memory while storing release site diameter");
    return(1);
  }
  mdlpvp->rsop->diameter->x = $<dbl>3 * scaling_factor;
  mdlpvp->rsop->diameter->y = $<dbl>3 * scaling_factor;
  mdlpvp->rsop->diameter->z = $<dbl>3 * scaling_factor;
}
	| site_size_cmd '=' array_expr_only
{
  double scaling_factor = 1.0/volp->length_unit;
  if ($<tok>1==SITE_RADIUS) scaling_factor*=2;
  
  if (mdlpvp->rsop->release_shape==SHAPE_LIST)
  {
    mdlerror("Release list diameters must be single valued.");
    return 1;
  }
  mdlpvp->el_head=$<nel>3;
  if ((mdlpvp->rsop->diameter=(struct vector3 *)malloc(sizeof(struct vector3)))==NULL) {
    mdlerror("Out of memory while storing release site diameter");
    return(1);
  }
  mdlpvp->elp=mdlpvp->el_head;
  if (mdlpvp->elp!=NULL) {
    mdlpvp->rsop->diameter->x=mdlpvp->elp->value * scaling_factor;
    mdlpvp->elp=mdlpvp->elp->next;
  }
  else {
    mdlerror("Three dimensional value required");
    return(1);
  }
  if (mdlpvp->elp!=NULL) {
    mdlpvp->rsop->diameter->y=mdlpvp->elp->value * scaling_factor;
    mdlpvp->elp=mdlpvp->elp->next;
  }
  else {
    mdlerror("Three dimensional value required");
    return(1);
  }
  if (mdlpvp->elp!=NULL) {
    mdlpvp->rsop->diameter->z=mdlpvp->elp->value * scaling_factor;
    mdlpvp->elp=mdlpvp->elp->next;
  }
  else {
    mdlerror("Three dimensional value required");
    return(1);
  }
  if (mdlpvp->elp!=NULL) {
    mdlerror("Three dimensional value required");
    return(1);
  }
}
	| site_size_cmd '=' existing_num_or_array
{
  double scaling_factor = 1.0/volp->length_unit;
  if ($<tok>1==SITE_RADIUS) scaling_factor *= 2;
  
  mdlpvp->gp=$<sym>3;
  if ((mdlpvp->rsop->diameter=(struct vector3 *)malloc(sizeof(struct vector3)))==NULL) {
    mdlerror("Out of memory while storing release site diameter");
    return(1);
  }
  switch (mdlpvp->gp->sym_type) {
  case DBL:
    mdlpvp->tmp_dbl = *(double *)mdlpvp->gp->value;
    mdlpvp->rsop->diameter->x = mdlpvp->tmp_dbl * scaling_factor;
    mdlpvp->rsop->diameter->y = mdlpvp->tmp_dbl * scaling_factor;
    mdlpvp->rsop->diameter->z = mdlpvp->tmp_dbl * scaling_factor;
    break;
  case ARRAY:
    if (mdlpvp->rsop->release_shape==SHAPE_LIST)
    {
      mdlerror("Release list diameters must be single valued.");
      return 1;
    }
    mdlpvp->el_head=(struct num_expr_list *)mdlpvp->gp->value;
    mdlpvp->elp=mdlpvp->el_head;
    if (mdlpvp->elp!=NULL) {
      mdlpvp->rsop->diameter->x=mdlpvp->elp->value * scaling_factor;
      mdlpvp->elp=mdlpvp->elp->next;
    }
    else {
      mdlerror("Three dimensional value required");
      return(1);
    }
    if (mdlpvp->elp!=NULL) {
      mdlpvp->rsop->diameter->y=mdlpvp->elp->value * scaling_factor;
      mdlpvp->elp=mdlpvp->elp->next;
    }
    else {
      mdlerror("Three dimensional value required");
      return(1);
    }
    if (mdlpvp->elp!=NULL) {
      mdlpvp->rsop->diameter->z=mdlpvp->elp->value * scaling_factor;
      mdlpvp->elp=mdlpvp->elp->next;
    }
    else {
      mdlerror("Three dimensional value required");
      return(1);
    }
    if (mdlpvp->elp!=NULL) {
      mdlerror("Three dimensional value required");
      return(1);
    }
    break;
  default:
    mdlerror("Diameter must either be a number or a 3-valued vector.");
    return 1;
  }
}
	| RELEASE_PROBABILITY '=' num_expr
{
  if (mdlpvp->rsop->release_prob==MAGIC_PATTERN_PROBABILITY)
  {
    mdlerror("Ignoring release probability for reaction-triggered releases.");
  }
  else
  {
    mdlpvp->rsop->release_prob=$<dbl>3;
    if (mdlpvp->rsop->release_prob<0)
    {
      mdlerror("Release probability cannot be less than 0.");
      return 1;
    }
    if (mdlpvp->rsop->release_prob>1)
    {
      mdlerror("Release probability cannot be greater than 1.");
      return 1;
    }
  }
}
	| RELEASE_PATTERN '=' existing_release_pattern_xor_rxpn
{
  mdlpvp->gp=$<sym>3;
  mdlpvp->rsop->pattern=(struct release_pattern *)mdlpvp->gp->value;
  
  if (mdlpvp->gp->sym_type==RXPN) /* Careful!  We've put a rxn_pathname into the "pattern" pointer! */
  {
    if (mdlpvp->rsop->release_prob!=1.0)
    {
      mdlerror("Ignoring release probability for reaction-triggered releases.");
    }
    mdlpvp->rsop->release_prob = MAGIC_PATTERN_PROBABILITY;   /* Magic number indicating a reaction-triggered release */
  }
}
	| MOLECULE_POSITIONS
{
  if (mdlpvp->rsop->release_shape != SHAPE_LIST)
  {
    mdlerror("You must use the LIST shape to specify molecule positions in a release.");
    return 1;
  }
}
	'{' molecule_release_pos_list '}'
{
  struct release_single_molecule *rsm,*old_head;
  
  old_head = mdlpvp->rsop->mol_list;
  mdlpvp->rsop->mol_list = NULL;
  mdlpvp->rsop->release_number = 0;
  while (old_head!=NULL)  /* Flip list from stack order to queue order */
  {
    rsm = old_head;
    old_head = old_head->next;
    
    rsm->next = mdlpvp->rsop->mol_list;
    mdlpvp->rsop->mol_list = rsm;
    mdlpvp->rsop->release_number++;
  }
};

site_size_cmd:
	SITE_DIAMETER {$$=SITE_DIAMETER;}
	| SITE_RADIUS {$$=SITE_RADIUS;};


release_number_cmd: constant_release_number_cmd
	| gaussian_release_number_cmd
	| volume_dependent_number_cmd
	| concentration_dependent_release_cmd
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
  mdlpvp->rsop->release_number=(int) $<dbl>5;
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

concentration_dependent_release_cmd: CONCENTRATION '=' num_expr
{
  mdlpvp->rsop->release_number_method=CCNNUM;
  mdlpvp->rsop->release_number=-3; /* Expect 3D molecule */
  mdlpvp->rsop->concentration = $<dbl>3;
}
	| DENSITY '=' num_expr
{
  mdlpvp->rsop->release_number_method=CCNNUM;
  mdlpvp->rsop->release_number=-2; /* Expect 2D molecule */
  mdlpvp->rsop->concentration = $<dbl>3;
};

molecule_release_pos_list: molecule_release_pos
	| molecule_release_pos_list molecule_release_pos;
	
molecule_release_pos:
	existing_molecule_opt_orient point
{
  struct release_single_molecule *rsm;
  struct vector3 temp_v3;
  
  memcpy(&temp_v3,$<vec3>2,sizeof(struct vector3));
  free($<vec3>2);
  
  rsm = (struct release_single_molecule*)malloc(sizeof(struct release_single_molecule));
  
  if (rsm==NULL)
  {
    mdlerror("Out of memory reading molecule positions");
    return 1;
  }
  
  rsm->orient = mdlpvp->orient_class;
  rsm->loc.x = temp_v3.x / mdlpvp->vol->length_unit;
  rsm->loc.y = temp_v3.y / mdlpvp->vol->length_unit;
  rsm->loc.z = temp_v3.z / mdlpvp->vol->length_unit;
  rsm->mol_type = (struct species*)( ($<sym>1)->value );
  
  if (rsm->mol_type->flags&ON_GRID)
  {
    if (!mdlpvp->orient_specified)
    {
      if (mdlpvp->vol->notify->missed_surf_orient==WARN_ERROR)
      {
	mdlerror("Error: surface orientation not specified for released surface molecule\n  (use ; or ', or ,' for random orientation)");
	return 1;
      }
      else if (mdlpvp->vol->notify->missed_surf_orient==WARN_WARN)
      {
	mdlerror("Warning: surface orientation not specified for released surface molecule\n  (use ; or ', or ,' for random orientation)");
      }
    }
  }
  else if ((rsm->mol_type->flags&NOT_FREE)==0)
  {
    if (mdlpvp->orient_specified)
    {
      if (mdlpvp->vol->notify->useless_vol_orient==WARN_ERROR)
      {
	mdlerror("Error: orientation not used for released volume molecule");
	return 1;
      }
      else if (mdlpvp->vol->notify->useless_vol_orient==WARN_WARN)
      {
	mdlerror("Warning: orientation not used for released volume molecule");
      }
    }
  }
  else
  {
    mdlerror("Error: cannot release a surface class instead of a molecule.");
    return 1;
  }

  rsm->next = mdlpvp->rsop->mol_list;
  mdlpvp->rsop->mol_list = rsm;        /* Acquire list in reverse order */
};

point: array_value
{
  mdlpvp->el_head=$<nel>1;
  if ((mdlpvp->pntp1=(struct vector3 *)malloc(sizeof(struct vector3)))==NULL) {
    mdlerror("Out of memory while creating points");
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
  if (mdlpvp->elp!=NULL) {
    mdlerror("Three dimensional value required");
    return(1);
  }
  $$=mdlpvp->pntp1;
};


polygon_list_def: new_object POLYGON_LIST '{'
{
  mdlpvp->gp=$<sym>1;
  mdlpvp->objp=(struct object *)mdlpvp->gp->value;
  if ((mdlpvp->pop=(struct polygon_object *)malloc
              (sizeof(struct polygon_object)))==NULL) {
    mdlerror("Out of memory while creating polygon list object");
    return(1);
  }
  mdlpvp->pop->polygon_data=NULL;
  mdlpvp->pop->n_walls=0;
  mdlpvp->pop->n_verts=0;
  mdlpvp->pop->fully_closed=0;
  mdlpvp->pop->surf_class=NULL;
  mdlpvp->pop->side_removed=NULL;

  if ((mdlpvp->opp=(struct ordered_poly *)malloc
              (sizeof(struct ordered_poly)))==NULL) {
    mdlerror("Out of memory while creating polygon list object");
    return(1);
  }
  mdlpvp->opp->vertex=NULL;
  mdlpvp->opp->normal=NULL;
  mdlpvp->opp->element=NULL;
  mdlpvp->opp->n_verts=0;
  mdlpvp->opp->n_walls=0;
  mdlpvp->pop->polygon_data=(void *)mdlpvp->opp;

  mdlpvp->objp->object_type=POLY_OBJ;
  mdlpvp->objp->contents=mdlpvp->pop;
  mdlpvp->objp->parent=mdlpvp->curr_obj;
  mdlpvp->obj_name=mdlpvp->objp->sym->name;
  mdlpvp->curr_obj=mdlpvp->objp;
  mdlpvp->region_list_head=mdlpvp->objp->regions;
}
	vertex_list_cmd
        element_connection_cmd
{
  u_int i,j;

  mdlpvp->pop->n_walls=mdlpvp->n_walls;
  mdlpvp->pop->n_verts=mdlpvp->n_verts;
  mdlpvp->opp->n_walls=mdlpvp->n_walls;
  mdlpvp->opp->n_verts=mdlpvp->n_verts;
  if ((mdlpvp->pop->surf_class=(struct species **)malloc
              (mdlpvp->pop->n_walls*sizeof(struct species *)))==NULL) {
    mdlerror("Out of memory while creating polygon list object");
    return(1);
  }
  for (i=0;i<mdlpvp->pop->n_walls;i++) mdlpvp->pop->surf_class[i]=volp->g_surf;
  mdlpvp->pop->side_removed = new_bit_array(mdlpvp->pop->n_walls);
  if (mdlpvp->pop->side_removed==NULL) {
    mdlerror("Out of memory while creating polygon list object");
    return(1);
  }
  set_all_bits(mdlpvp->pop->side_removed,0);

  if ((mdlpvp->opp->vertex=(struct vector3 *)malloc
              (mdlpvp->opp->n_verts*sizeof(struct vector3)))==NULL) {
    mdlerror("Out of memory while creating polygon list object");
    return(1);
  }

  mdlpvp->vlp=mdlpvp->vertex_head;
  if (mdlpvp->vlp->normal!=NULL) {
    if ((mdlpvp->opp->normal=(struct vector3 *)malloc
                (mdlpvp->opp->n_verts*sizeof(struct vector3)))==NULL) {
      mdlerror("Out of memory while creating polygon list object");
      return(1);
    }
  }
  for (i=0;i<mdlpvp->opp->n_verts;i++) {
    mdlpvp->opp->vertex[i].x=mdlpvp->vlp->vertex->x;
    mdlpvp->opp->vertex[i].y=mdlpvp->vlp->vertex->y;
    mdlpvp->opp->vertex[i].z=mdlpvp->vlp->vertex->z;
    mdlpvp->vlp_temp=mdlpvp->vlp;
    free(mdlpvp->vlp_temp->vertex);
    if (mdlpvp->opp->normal!=NULL) {
      mdlpvp->opp->normal[i].x=mdlpvp->vlp->normal->x;
      mdlpvp->opp->normal[i].y=mdlpvp->vlp->normal->y;
      mdlpvp->opp->normal[i].z=mdlpvp->vlp->normal->z;
      free(mdlpvp->vlp_temp->normal);
    }
    mdlpvp->vlp=mdlpvp->vlp->next;
    free(mdlpvp->vlp_temp);
  }
  if ((mdlpvp->edp=(struct element_data *)malloc
              (mdlpvp->opp->n_walls*sizeof(struct element_data)))==NULL) {
    mdlerror("Out of memory while creating polygon list object");
    return(1);
  }
  mdlpvp->opp->element=mdlpvp->edp;
  mdlpvp->eclp=mdlpvp->connection_head;
  for (i=0;i<mdlpvp->opp->n_walls;i++) {
    if (mdlpvp->eclp->n_verts != 3)
    {
      mdlerror("All polygons must have three vertices.");
      return(1);
    }
    mdlpvp->elp=mdlpvp->eclp->connection_list;
    for (j=0;j<mdlpvp->eclp->n_verts;j++) {
      mdlpvp->edp[i].vertex_index[j]=(int)mdlpvp->elp->value;
      mdlpvp->elp_temp=mdlpvp->elp;
      mdlpvp->elp=mdlpvp->elp->next;
      free(mdlpvp->elp_temp);
    }
    mdlpvp->eclp_temp=mdlpvp->eclp;
    mdlpvp->eclp=mdlpvp->eclp->next;
    free(mdlpvp->eclp_temp);
  }

  /* Create object default region on polygon list object: */
  if ((mdlpvp->rp=make_new_region(volp,mdlpvp->obj_name,"ALL",mdlpvp->mdl_err_msg))==NULL) {
    mdlerror(mdlpvp->mdl_err_msg,mdlpvp);
    return(1);
  }
  if ((mdlpvp->rlp=(struct region_list *)malloc(sizeof(struct region_list)))==NULL) {
    sprintf(mdlpvp->mdl_err_msg,"%s %s",
      "Out of memory while creating object default region name:",mdlpvp->rp->sym->name);
    mdlerror(mdlpvp->mdl_err_msg,mdlpvp);
    return(1);
  }
  mdlpvp->rp->region_last_name="ALL";
  mdlpvp->rp->parent=mdlpvp->curr_obj;
  mdlpvp->curr_obj->num_regions++;

  if ((mdlpvp->elmlp=(struct element_list *)malloc
             (sizeof(struct element_list)))==NULL) {
    mdlerror("Out of memory while creating element list for object default region");
    return(1);
  }
  mdlpvp->elmlp->begin=0;
  mdlpvp->elmlp->end=mdlpvp->pop->n_walls-1;
  mdlpvp->elmlp->special=NULL;
  mdlpvp->elmlp->next=mdlpvp->rp->element_list_head;
  mdlpvp->rp->element_list_head=mdlpvp->elmlp;
  mdlpvp->curr_obj->n_walls=mdlpvp->pop->n_walls;
  mdlpvp->curr_obj->n_verts=mdlpvp->pop->n_verts;
  mdlpvp->n_walls_actual = mdlpvp->pop->n_walls;
  normalize_elements(mdlpvp->rp,0);

  mdlpvp->rp->surf_class=NULL;
  mdlpvp->rlp->reg=mdlpvp->rp;
  mdlpvp->rlp->next=mdlpvp->region_list_head;
  mdlpvp->region_list_head=mdlpvp->rlp;
  no_printf("Creating object default region: %s\n",mdlpvp->rp->sym->name);
}
	list_opt_polygon_object_cmds
	list_opt_object_cmds
	'}'
{
  mdlpvp->curr_obj->regions=mdlpvp->region_list_head;
  remove_gaps_from_regions(mdlpvp->curr_obj);
  mdlpvp->n_walls_actual = mdlpvp->curr_obj->n_walls_actual;
  no_printf("Polygon list %s defined:\n",mdlpvp->curr_obj->sym->name);
  no_printf(" n_verts = %d\n",mdlpvp->pop->n_verts);
  no_printf(" n_walls = %d\n",mdlpvp->pop->n_walls);
  mdlpvp->curr_obj->parent->n_walls+=mdlpvp->curr_obj->n_walls;
  mdlpvp->curr_obj->parent->n_walls_actual+=mdlpvp->curr_obj->n_walls_actual;
  mdlpvp->curr_obj->parent->n_verts+=mdlpvp->curr_obj->n_verts;
  mdlpvp->curr_obj=mdlpvp->curr_obj->parent;
  if (mdlpvp->object_name_list_end->prev!=NULL) {
    mdlpvp->object_name_list_end=mdlpvp->object_name_list_end->prev;
  }
  else {
    mdlpvp->object_name_list_end->name=NULL;
  }
  $$=$<sym>1;
};


voxel_list_def: new_object VOXEL_LIST '{'
{
  mdlpvp->gp=$<sym>1;
  mdlpvp->objp=(struct object *)mdlpvp->gp->value;
  if ((mdlpvp->vop=(struct voxel_object *)malloc
              (sizeof(struct voxel_object)))==NULL) {
    mdlerror("Out of memory while creating voxel list object");
    return(1);
  }
  mdlpvp->vop->voxel_data=NULL;
  mdlpvp->vop->n_voxels=0;
  mdlpvp->vop->n_verts=0;
  mdlpvp->vop->fully_closed=0;

  
  if ((mdlpvp->ovp=(struct ordered_voxel *)malloc
              (sizeof(struct ordered_voxel)))==NULL) {
    mdlerror("Out of memory while creating voxel list object");
    return(1);
  }
  mdlpvp->ovp->vertex=NULL;
  mdlpvp->ovp->element=NULL;
  mdlpvp->ovp->neighbor=NULL;
  mdlpvp->ovp->n_verts=0;
  mdlpvp->ovp->n_voxels=0;
  mdlpvp->vop->voxel_data=(void *)mdlpvp->ovp;

  mdlpvp->objp->object_type=VOXEL_OBJ;
  mdlpvp->objp->contents=mdlpvp->vop;
  mdlpvp->objp->parent=mdlpvp->curr_obj;
  mdlpvp->obj_name=mdlpvp->objp->sym->name;
  mdlpvp->curr_obj=mdlpvp->objp;
}
	vertex_list_cmd
        tet_element_connection_cmd 
	'}'
{
  u_int i,j;

  mdlpvp->vop->n_voxels=mdlpvp->n_voxels;
  mdlpvp->vop->n_verts=mdlpvp->n_verts;
  mdlpvp->ovp->n_voxels=mdlpvp->n_voxels;
  mdlpvp->ovp->n_verts=mdlpvp->n_verts;

  if ((mdlpvp->ovp->vertex=(struct vector3 *)malloc
              (mdlpvp->ovp->n_verts*sizeof(struct vector3)))==NULL) {
    mdlerror("Out of memory while creating voxel list object");
    return(1);
  }

  mdlpvp->vlp=mdlpvp->vertex_head;
  for (i=0;i<mdlpvp->ovp->n_verts;i++) {
    mdlpvp->ovp->vertex[i].x=mdlpvp->vlp->vertex->x;
    mdlpvp->ovp->vertex[i].y=mdlpvp->vlp->vertex->y;
    mdlpvp->ovp->vertex[i].z=mdlpvp->vlp->vertex->z;
    mdlpvp->vlp_temp=mdlpvp->vlp;
    free(mdlpvp->vlp_temp->vertex);
    mdlpvp->vlp=mdlpvp->vlp->next;
    free(mdlpvp->vlp_temp);
  }
  if ((mdlpvp->tedp=(struct tet_element_data *)malloc
              (mdlpvp->ovp->n_voxels*sizeof(struct tet_element_data)))==NULL) {
    mdlerror("Out of memory while creating voxel list object");
    return(1);
  }
  mdlpvp->ovp->element=mdlpvp->tedp;
  mdlpvp->eclp=mdlpvp->connection_head;
  for (i=0;i<mdlpvp->ovp->n_voxels;i++) {
    if (mdlpvp->eclp->n_verts != 4)
    {
      mdlerror("All voxels must have four vertices.");
      return(1);
    }
    mdlpvp->tedp[i].n_verts=mdlpvp->eclp->n_verts;
    mdlpvp->elp=mdlpvp->eclp->connection_list;
    for (j=0;j<mdlpvp->eclp->n_verts;j++) {
      mdlpvp->tedp[i].vertex_index[j]=(int)mdlpvp->elp->value;
      mdlpvp->elp_temp=mdlpvp->elp;
      mdlpvp->elp=mdlpvp->elp->next;
      free(mdlpvp->elp_temp);
    }
    mdlpvp->eclp_temp=mdlpvp->eclp;
    mdlpvp->eclp=mdlpvp->eclp->next;
    free(mdlpvp->eclp_temp);
  }
  no_printf("Voxel list %s defined:\n",mdlpvp->curr_obj->sym->name);
  no_printf(" n_verts = %d\n",mdlpvp->vop->n_verts);
  /*mdlpvp->curr_obj->parent->n_verts+=mdlpvp->curr_obj->n_verts; 
  mdlpvp->curr_obj=mdlpvp->curr_obj->parent;
  if (mdlpvp->object_name_list_end->prev!=NULL) {
    mdlpvp->object_name_list_end=mdlpvp->object_name_list_end->prev;
  }
  else {
    mdlpvp->object_name_list_end->name=NULL;
  }
  */

  $$=$<sym>1;

};

vertex_list_cmd: VERTEX_LIST '{' 
{
  mdlpvp->n_verts=0;
  mdlpvp->vertex_head=NULL;
  mdlpvp->vertex_tail=NULL;
}
	list_points 
        '}'
;


list_points: point
{
  mdlpvp->pntp1=$<vec3>1;
  if ((mdlpvp->vlp=(struct vertex_list *)malloc(sizeof(struct vertex_list)))==NULL) {
    mdlerror("Out of memory while creating vertices");
    return(1);
  }
  mdlpvp->n_verts++;
  mdlpvp->vlp->vertex=mdlpvp->pntp1;
  mdlpvp->vlp->normal=NULL;
  if (mdlpvp->vertex_tail==NULL) {
    mdlpvp->vertex_tail=mdlpvp->vlp;
  }
  mdlpvp->vertex_tail->next=mdlpvp->vlp;
  mdlpvp->vlp->next=NULL;
  mdlpvp->vertex_tail=mdlpvp->vlp;
  if (mdlpvp->vertex_head==NULL) {
    mdlpvp->vertex_head=mdlpvp->vlp;
  }
}
        | point NORMAL point
{
  mdlpvp->pntp1=$<vec3>1;
  mdlpvp->pntp2=$<vec3>3;
  if ((mdlpvp->vlp=(struct vertex_list *)malloc(sizeof(struct vertex_list)))==NULL) {
    mdlerror("Out of memory while creating normals");
    return(1);
  }
  mdlpvp->n_verts++;
  mdlpvp->vlp->vertex=mdlpvp->pntp1;
  mdlpvp->vlp->normal=mdlpvp->pntp2;
  if (mdlpvp->vertex_tail==NULL) {
    mdlpvp->vertex_tail=mdlpvp->vlp;
  }
  mdlpvp->vertex_tail->next=mdlpvp->vlp;
  mdlpvp->vlp->next=NULL;
  mdlpvp->vertex_tail=mdlpvp->vlp;
  if (mdlpvp->vertex_head==NULL) {
    mdlpvp->vertex_head=mdlpvp->vlp;
  }
}
	| list_points point
{
  mdlpvp->pntp1=$<vec3>2;
  if ((mdlpvp->vlp=(struct vertex_list *)malloc(sizeof(struct vertex_list)))==NULL) {
    mdlerror("Out of memory while creating vertices");
    return(1);
  }
  mdlpvp->n_verts++;
  mdlpvp->vlp->vertex=mdlpvp->pntp1;
  mdlpvp->vlp->normal=NULL;
  if (mdlpvp->vertex_tail==NULL) {
    mdlpvp->vertex_tail=mdlpvp->vlp;
  }
  mdlpvp->vertex_tail->next=mdlpvp->vlp;
  mdlpvp->vlp->next=NULL;
  mdlpvp->vertex_tail=mdlpvp->vlp;
  if (mdlpvp->vertex_head==NULL) {
    mdlpvp->vertex_head=mdlpvp->vlp;
  }
}
	| list_points point NORMAL point
{
  mdlpvp->pntp1=$<vec3>2;
  mdlpvp->pntp2=$<vec3>4;
  if ((mdlpvp->vlp=(struct vertex_list *)malloc(sizeof(struct vertex_list)))==NULL) {
    mdlerror("Out of memory while creating normals");
    return(1);
  }
  mdlpvp->n_verts++;
  mdlpvp->vlp->vertex=mdlpvp->pntp1;
  mdlpvp->vlp->normal=mdlpvp->pntp2;
  if (mdlpvp->vertex_tail==NULL) {
    mdlpvp->vertex_tail=mdlpvp->vlp;
  }
  mdlpvp->vertex_tail->next=mdlpvp->vlp;
  mdlpvp->vlp->next=NULL;
  mdlpvp->vertex_tail=mdlpvp->vlp;
  if (mdlpvp->vertex_head==NULL) {
    mdlpvp->vertex_head=mdlpvp->vlp;
  }
};


element_connection_cmd: ELEMENT_CONNECTIONS '{'
{
  mdlpvp->n_walls=0;
  mdlpvp->n_walls_actual=0;
  mdlpvp->connection_head=NULL;
  mdlpvp->connection_tail=NULL;
}
        list_arrays
        '}'
{
  mdlpvp->n_walls_actual=mdlpvp->n_walls;
};

tet_element_connection_cmd: TET_ELEMENT_CONNECTIONS '{'
{
  mdlpvp->n_voxels=0;
  mdlpvp->connection_head=NULL;
  mdlpvp->connection_tail=NULL;
}
        list_tet_arrays
        '}'
;

list_arrays: array_value
{
  mdlpvp->elp=$<nel>1;
  if ((mdlpvp->eclp=(struct element_connection_list *)malloc
              (sizeof(struct element_connection_list)))==NULL) {
    mdlerror("Out of memory while creating element connections");
    return(1);
  }
  mdlpvp->n_walls++;
  mdlpvp->eclp->connection_list=mdlpvp->elp;
  mdlpvp->eclp->n_verts=0;
  while (mdlpvp->elp!=NULL) {
    mdlpvp->eclp->n_verts++;
    mdlpvp->elp=mdlpvp->elp->next;
  }
  if (mdlpvp->eclp->n_verts!=3) {
    mdlerror("Non-triangular element found in polygon list object");
    return(1);
  }
  if (mdlpvp->connection_tail==NULL) {
    mdlpvp->connection_tail=mdlpvp->eclp;
  }
  mdlpvp->connection_tail->next=mdlpvp->eclp;
  mdlpvp->eclp->next=NULL;
  mdlpvp->connection_tail=mdlpvp->eclp;
  if (mdlpvp->connection_head==NULL) {
    mdlpvp->connection_head=mdlpvp->eclp;
  }
}
        | list_arrays array_value
{
  mdlpvp->elp=$<nel>2;
  if ((mdlpvp->eclp=(struct element_connection_list *)malloc
              (sizeof(struct element_connection_list)))==NULL) {
    mdlerror("Out of memory while creating element connections");
    return(1);
  }
  mdlpvp->n_walls++;
  mdlpvp->eclp->connection_list=mdlpvp->elp;
  mdlpvp->eclp->n_verts=0;
  while (mdlpvp->elp!=NULL) {
    mdlpvp->eclp->n_verts++;
    mdlpvp->elp=mdlpvp->elp->next;
  }
  if (mdlpvp->connection_tail==NULL) {
    mdlpvp->connection_tail=mdlpvp->eclp;
  }
  mdlpvp->connection_tail->next=mdlpvp->eclp;
  mdlpvp->eclp->next=NULL;
  mdlpvp->connection_tail=mdlpvp->eclp;
  if (mdlpvp->connection_head==NULL) {
    mdlpvp->connection_head=mdlpvp->eclp;
  }
};

list_tet_arrays: array_value
{
  mdlpvp->elp=$<nel>1;
  if ((mdlpvp->eclp=(struct element_connection_list *)malloc
              (sizeof(struct element_connection_list)))==NULL) {
    mdlerror("Out of memory while creating element connections");
    return(1);
  }
  mdlpvp->n_voxels++;
  mdlpvp->eclp->connection_list=mdlpvp->elp;
  mdlpvp->eclp->n_verts=0;
  while (mdlpvp->elp!=NULL) {
    mdlpvp->eclp->n_verts++;
    mdlpvp->elp=mdlpvp->elp->next;
  }
  if (mdlpvp->eclp->n_verts!=4) {
    mdlerror("Non-tetrahedron element found in voxel list object");
    return(1);
  }
  if (mdlpvp->connection_tail==NULL) {
    mdlpvp->connection_tail=mdlpvp->eclp;
  }
  mdlpvp->connection_tail->next=mdlpvp->eclp;
  mdlpvp->eclp->next=NULL;
  mdlpvp->connection_tail=mdlpvp->eclp;
  if (mdlpvp->connection_head==NULL) {
    mdlpvp->connection_head=mdlpvp->eclp;
  }
}
        | list_tet_arrays array_value
{
  mdlpvp->elp=$<nel>2;
  if ((mdlpvp->eclp=(struct element_connection_list *)malloc
              (sizeof(struct element_connection_list)))==NULL) {
    mdlerror("Out of memory while creating element connections");
    return(1);
  }
  mdlpvp->n_voxels++;
  mdlpvp->eclp->connection_list=mdlpvp->elp;
  mdlpvp->eclp->n_verts=0;
  while (mdlpvp->elp!=NULL) {
    mdlpvp->eclp->n_verts++;
    mdlpvp->elp=mdlpvp->elp->next;
  }
  if (mdlpvp->connection_tail==NULL) {
    mdlpvp->connection_tail=mdlpvp->eclp;
  }
  mdlpvp->connection_tail->next=mdlpvp->eclp;
  mdlpvp->eclp->next=NULL;
  mdlpvp->connection_tail=mdlpvp->eclp;
  if (mdlpvp->connection_head==NULL) {
    mdlpvp->connection_head=mdlpvp->eclp;
  }
};

boolean: TRUE { $$=1; }
	| FALSE { $$=0; }
	| YES { $$=1; }
	| NO { $$=0; }
        | ON { $$=1; }
        | OFF { $$=0; }
;


box_def: new_object BOX '{'
{
  mdlpvp->gp=$<sym>1;
  mdlpvp->objp=(struct object *)mdlpvp->gp->value;
  if ((mdlpvp->pop=(struct polygon_object *)malloc
              (sizeof(struct polygon_object)))==NULL) {
    mdlerror("Out of memory while creating box object");
    return(1);
  }
  mdlpvp->pop->polygon_data=NULL;
  mdlpvp->pop->sb=NULL;
  mdlpvp->pop->n_walls=0;
  mdlpvp->pop->fully_closed=1;
  mdlpvp->pop->surf_class=NULL;
  mdlpvp->pop->side_removed=NULL;

  if ((mdlpvp->opp=(struct ordered_poly *)malloc
              (sizeof(struct ordered_poly)))==NULL) {
    mdlerror("Out of memory while creating box object");
    return(1);
  }
  mdlpvp->opp->vertex=NULL;
  mdlpvp->opp->normal=NULL;
  mdlpvp->opp->element=NULL;
  mdlpvp->opp->n_verts=0;
  mdlpvp->opp->n_walls=0;
  mdlpvp->pop->polygon_data=(void *)mdlpvp->opp;

  mdlpvp->llf=NULL;
  mdlpvp->urb=NULL;

  mdlpvp->objp->object_type=BOX_OBJ;
  mdlpvp->objp->contents=mdlpvp->pop;
  mdlpvp->objp->parent=mdlpvp->curr_obj;
  mdlpvp->obj_name=mdlpvp->objp->sym->name;
  mdlpvp->curr_obj=mdlpvp->objp;
  mdlpvp->region_list_head=mdlpvp->objp->regions;

  /* Create object default region on box object: */
  if ((mdlpvp->rp=make_new_region(volp,mdlpvp->obj_name,"ALL",mdlpvp->mdl_err_msg))==NULL) {
    mdlerror(mdlpvp->mdl_err_msg,mdlpvp);
    return(1);
  }
  if ((mdlpvp->rlp=(struct region_list *)malloc(sizeof(struct region_list)))==NULL) {
    sprintf(mdlpvp->mdl_err_msg,"%s %s",
      "Out of memory while creating object default region name:",mdlpvp->rp->sym->name);
    mdlerror(mdlpvp->mdl_err_msg,mdlpvp);
    return(1);
  }
  mdlpvp->rp->region_last_name="ALL";
  mdlpvp->rp->parent=mdlpvp->curr_obj;
  mdlpvp->curr_obj->num_regions++;

  if ((mdlpvp->elmlp=(struct element_list *)malloc
             (sizeof(struct element_list)))==NULL) {
    mdlerror("Out of memory while creating element list for object default region");
    return(1);
  }
  mdlpvp->elmlp->begin=ALL_SIDES;
  mdlpvp->elmlp->end=ALL_SIDES;
  mdlpvp->elmlp->special=NULL;
  mdlpvp->elmlp->next=mdlpvp->rp->element_list_head;
  mdlpvp->rp->element_list_head=mdlpvp->elmlp;

  mdlpvp->rp->surf_class=NULL;
  mdlpvp->rlp->reg=mdlpvp->rp;
  mdlpvp->rlp->next=mdlpvp->region_list_head;
  mdlpvp->region_list_head=mdlpvp->rlp;
  no_printf("Creating object default region: %s\n",mdlpvp->rp->sym->name);
}
	CORNERS '=' point ',' point
{
  mdlpvp->llf=$<vec3>7;
  mdlpvp->urb=$<vec3>9;
  mdlpvp->pop->sb = init_cuboid(mdlpvp->llf,mdlpvp->urb);
  free(mdlpvp->llf);
  free(mdlpvp->urb);
  mdlpvp->llf = mdlpvp->urb = NULL;
  if (mdlpvp->pop->sb == NULL)
  {
    mdlerror("Out of memory while storing box corners");
    return 1;
  }
  mdlpvp->box_aspect_ratio = 0.0;
}
	opt_aspect_ratio_def
	list_opt_polygon_object_cmds
{
  int i;

  if (mdlpvp->box_aspect_ratio >= 2.0)
  {
    i = reaspect_cuboid(mdlpvp->pop->sb,mdlpvp->box_aspect_ratio);
    if (i)
    {
      mdlerror("Error setting up box geometry");
      return 1;
    }
  }
  for (mdlpvp->rlp=mdlpvp->region_list_head ; mdlpvp->rlp!=NULL; mdlpvp->rlp=mdlpvp->rlp->next)
  {
    i = normalize_elements(mdlpvp->rlp->reg,0);
    if (i)
    {
      mdlerror("Error setting up elements in box regions");
      return 1;
    }
  }
  i = polygonalize_cuboid(mdlpvp->opp,mdlpvp->pop->sb);
  if (i)
  {
    mdlerror("Could not turn box object into polygons");
    return 1;
  }
  else if (volp->notify->box_triangulation==NOTIFY_FULL)
  {
    fprintf(volp->log_file,"Box object %s converted into %d polygons\n",mdlpvp->curr_obj->sym->name,mdlpvp->opp->n_walls);
  }
  
  mdlpvp->pop->n_walls = mdlpvp->opp->n_walls;
  mdlpvp->pop->n_verts = mdlpvp->opp->n_verts;
  mdlpvp->n_walls_actual = mdlpvp->opp->n_walls;
  /*mdlpvp->pop->n_walls = 12;*/
  /*mdlpvp->pop->n_verts = 8;*/
  
  if ((mdlpvp->pop->surf_class=(struct species **)malloc
              (mdlpvp->pop->n_walls*sizeof(struct species *)))==NULL) {
    mdlerror("Out of memory while creating box object");
    return(1);
  }
  for (i=0;i<mdlpvp->pop->n_walls;i++) mdlpvp->pop->surf_class[i]=volp->g_surf;
  mdlpvp->pop->side_removed = new_bit_array(mdlpvp->pop->n_walls);
  if (mdlpvp->pop->side_removed==NULL) {
    mdlerror("Out of memory while creating box object");
    return(1);
  }
  set_all_bits(mdlpvp->pop->side_removed,0);
}	  
	list_opt_object_cmds
	'}'
{
  mdlpvp->curr_obj->regions=mdlpvp->region_list_head;
  remove_gaps_from_regions(mdlpvp->curr_obj);
  mdlpvp->n_walls_actual = mdlpvp->curr_obj->n_walls_actual;
#ifdef DEBUG
  no_printf("Box %s defined:\n",mdlpvp->curr_obj->sym->name);
  no_printf("    LLF Corner = [ %f, %f, %f ]\n",mdlpvp->llf->x,mdlpvp->llf->y,mdlpvp->llf->z);
  no_printf("    URB Corner = [ %f, %f, %f ]\n",mdlpvp->urb->x,mdlpvp->urb->y,mdlpvp->urb->z);
#endif
  mdlpvp->curr_obj->n_walls=mdlpvp->pop->n_walls;
  mdlpvp->curr_obj->n_verts=mdlpvp->pop->n_verts;
  mdlpvp->curr_obj->parent->n_walls+=mdlpvp->curr_obj->n_walls;
  mdlpvp->curr_obj->parent->n_walls_actual+=mdlpvp->curr_obj->n_walls_actual;
  mdlpvp->curr_obj->parent->n_verts+=mdlpvp->curr_obj->n_verts;
  mdlpvp->curr_obj=mdlpvp->curr_obj->parent;
  if (mdlpvp->object_name_list_end->prev!=NULL) {
    mdlpvp->object_name_list_end=mdlpvp->object_name_list_end->prev;
  }
  else {
    mdlpvp->object_name_list_end->name=NULL;
  }
  $$=$<sym>1;
};

opt_aspect_ratio_def:
	| aspect_ratio_def
;

aspect_ratio_def: ASPECT_RATIO '=' num_expr
{
  mdlpvp->box_aspect_ratio = $<dbl>3;
  if (!(mdlpvp->box_aspect_ratio >= 2.0))
  {
    mdlerror("Invalid aspect ratio selected (must be greater than or equal to 2.0)");
    return 1;
  }
}

list_opt_polygon_object_cmds: /* empty */
	| list_opt_polygon_object_cmds opt_polygon_object_cmd
;


opt_polygon_object_cmd:
	remove_side
	| in_obj_define_surface_regions
;


remove_side: REMOVE_ELEMENTS '{'
{
  mdlpvp->rp = retrieve_old_region(volp,mdlpvp->obj_name,"REMOVED",mdlpvp->mdl_err_msg);
  if (mdlpvp->rp==NULL)
  {
    if ((mdlpvp->rp=make_new_region(volp,mdlpvp->obj_name,"REMOVED",mdlpvp->mdl_err_msg))==NULL) {
      mdlerror(mdlpvp->mdl_err_msg,mdlpvp);
      return(1);
    }
    if ((mdlpvp->rlp=(struct region_list *)malloc(sizeof(struct region_list)))==NULL) {
      sprintf(mdlpvp->mdl_err_msg,"%s %s",
	"Out of memory while creating object default region name:",mdlpvp->rp->sym->name);
      mdlerror(mdlpvp->mdl_err_msg,mdlpvp);
      return(1);
    }
    mdlpvp->rp->region_last_name="REMOVED";
    mdlpvp->rp->parent=mdlpvp->curr_obj;
    mdlpvp->rp->surf_class=(struct species*)&mdlpvp->rp->surf_class;
    mdlpvp->rlp->reg=mdlpvp->rp;
    mdlpvp->rlp->next = mdlpvp->region_list_head;
    mdlpvp->region_list_head = mdlpvp->rlp;
    mdlpvp->objp->num_regions++;
    mdlpvp->element_list_head = NULL;
  }
  else
  {
    mdlpvp->element_list_head = mdlpvp->rp->element_list_head;
  }
  mdlpvp->exclude_me = 0;  /* Default is named regions are in removal list */
}
	remove_element_specifier_list '}'
{
  mdlpvp->rp->element_list_head = mdlpvp->element_list_head;
  if (mdlpvp->objp->object_type==POLY_OBJ)
  {
    if (normalize_elements(mdlpvp->rp,0))
    {
      mdlerror("Improper element specification: out of range or out of memory");
      return 1;
    }
  }
};

remove_element_specifier_list:
	element_specifier_list
	| just_an_element_list
;


side_name: TOP {$$=Z_POS;}
	| BOTTOM {$$=Z_NEG;}
	| FRONT {$$=Y_NEG;}
	| BACK {$$=Y_POS;}
	| LEFT {$$=X_NEG;}
	| RIGHT {$$=X_POS;}
	| ALL_ELEMENTS {$$=ALL_SIDES;}
;

element_specifier_list: element_specifier
	| element_specifier_list ',' element_specifier
;

element_specifier:
	incl_element_list_stmt
	| excl_element_list_stmt
	| prev_region_stmt
	| patch_statement
;

incl_element_list_stmt: INCLUDE_ELEMENTS 
{
  mdlpvp->exclude_me = 0;
}
	'=' '[' list_element_specs ']'
;

excl_element_list_stmt: EXCLUDE_ELEMENTS
{
  mdlpvp->exclude_me = 1;
}
	'=' '[' list_element_specs ']'
;

just_an_element_list:
	list_element_specs
;

list_element_specs: element_spec
	| list_element_specs ',' element_spec
;


element_spec: num_expr
{
  if ((mdlpvp->elmlp=(struct element_list *)malloc
             (sizeof(struct element_list)))==NULL) {
    mdlerror("Out of memory while creating element list item");
    return(1);
  }
  mdlpvp->elmlp->begin=(unsigned int) ($<dbl>1);
  mdlpvp->elmlp->end=mdlpvp->elmlp->begin;
  mdlpvp->elmlp->next=mdlpvp->element_list_head;
  mdlpvp->element_list_head=mdlpvp->elmlp;
  if (mdlpvp->exclude_me) mdlpvp->elmlp->special = (struct element_special*)mdlpvp->elmlp;
  else mdlpvp->elmlp->special=NULL;
}
	| num_expr TO num_expr
{
  if ((mdlpvp->elmlp=(struct element_list *)malloc
             (sizeof(struct element_list)))==NULL) {
    mdlerror("Out of memory while creating element list item");
    return(1);
  }
  mdlpvp->elmlp->begin=(unsigned int) ($<dbl>1);
  mdlpvp->elmlp->end=(unsigned int) ($<dbl>3);
  mdlpvp->elmlp->next=mdlpvp->element_list_head;
  mdlpvp->element_list_head=mdlpvp->elmlp;
  if (mdlpvp->exclude_me) mdlpvp->elmlp->special = (struct element_special*)mdlpvp->elmlp;
  else mdlpvp->elmlp->special=NULL;
}
	| side_name
{
  if ((mdlpvp->elmlp=(struct element_list *)malloc
             (sizeof(struct element_list)))==NULL) {
    mdlerror("Out of memory while creating element list item");
    return(1);
  }
  mdlpvp->elmlp->next=mdlpvp->element_list_head;
  mdlpvp->element_list_head=mdlpvp->elmlp;
  if (mdlpvp->exclude_me) mdlpvp->elmlp->special = (struct element_special*)mdlpvp->elmlp;
  else mdlpvp->elmlp->special=NULL;

  if ($<tok>1==ALL_SIDES && mdlpvp->objp->object_type==POLY_OBJ) {
    mdlpvp->elmlp->begin=0;
    mdlpvp->elmlp->end=mdlpvp->pop->n_walls-1;
  }
  else if (mdlpvp->objp->object_type==POLY_OBJ) {
    mdlerror("Illegal reference to polygon list element by side-name");
    return(1);
  }
  else {
    mdlpvp->elmlp->begin=$<tok>1;
    mdlpvp->elmlp->end=$<tok>1;
  }

};

prev_region_stmt: prev_region_type '=' VAR
{
  struct sym_table *stp;
  char *reg_name;
  char *full_reg_name;
  
  mdlpvp->elmlp = (struct element_list*)malloc(sizeof(struct element_list));
  if (mdlpvp->elmlp==NULL)
  {
    mdlerror("Out of memory while creating element list item");
    return 1;
  }
  mdlpvp->elmlp->special = (struct element_special*)malloc(sizeof(struct element_special));
  if (mdlpvp->elmlp->special==NULL)
  {
    mdlerror("Out of memory while creating element list item");
    return 1;
  }
  mdlpvp->elmlp->begin = 0;
  mdlpvp->elmlp->end = 0;
  mdlpvp->elmlp->special->exclude = (byte)$<tok>1;
  
  if (mdlpvp->cval_2!=NULL)
  {
    reg_name = mdlpvp->cval_2;
    mdlpvp->cval_2 = NULL;
  }
  else
  {
    reg_name = mdlpvp->cval;
    mdlpvp->cval=NULL;
  }
  
  full_reg_name = (char*)malloc( strlen(reg_name) + strlen(mdlpvp->curr_obj->sym->name) + 2 );
  if (full_reg_name==NULL)
  {
    mdlerror("Out of memory while creating element list item");
    return 1;
  }
  
  strcpy(full_reg_name,mdlpvp->curr_obj->sym->name);
  strcat(full_reg_name,",");
  strcat(full_reg_name,reg_name);

  stp=retrieve_sym(full_reg_name,REG,volp->main_sym_table);

  if (reg_name==mdlpvp->cval) mdlpvp->cval=NULL;
  else mdlpvp->cval_2 = NULL;
  
  if (stp==NULL)
  {
    sprintf(mdlpvp->mdl_err_msg,"%s %s","Undefined region:",full_reg_name);
    mdlerror(mdlpvp->mdl_err_msg,mdlpvp);
    return 1;
  }
   
  mdlpvp->elmlp->special->referent = (struct region*)stp->value;
  
  if (mdlpvp->elmlp->special->referent==mdlpvp->rp)
  {
    mdlerror("Self-referential region include.  No paradoxes, please.");
    return 1;
  }
  free(full_reg_name);
  free(reg_name);
  mdlpvp->elmlp->next=mdlpvp->element_list_head;
  mdlpvp->element_list_head=mdlpvp->elmlp;
};

prev_region_type: INCLUDE_REGION { $$=0; }
	| EXCLUDE_REGION { $$=1; }
;

patch_statement: patch_type '=' point ',' point
{
  struct vector3 *temp_llf;
  struct vector3 *temp_urb;
  
  if (mdlpvp->objp->object_type!=BOX_OBJ)
  {
    mdlerror("Must use PATCH specifier on a BOX object only.");
    return 1;
  }
  
  mdlpvp->elmlp = (struct element_list*)malloc(sizeof(struct element_list));
  if (mdlpvp->elmlp==NULL)
  {
    mdlerror("Out of memory while creating element list item");
    return 1;
  }
  mdlpvp->elmlp->special = (struct element_special*)malloc(sizeof(struct element_special));
  if (mdlpvp->elmlp->special==NULL)
  {
    mdlerror("Out of memory while creating element list item");
    return 1;
  }
  mdlpvp->elmlp->begin = 0;
  mdlpvp->elmlp->end = 0;
  mdlpvp->elmlp->special->referent = NULL;
  mdlpvp->elmlp->special->exclude = (byte)$<tok>1;
  temp_llf = $<vec3>3;
  temp_urb = $<vec3>5;
  memcpy(&(mdlpvp->elmlp->special->corner1),temp_llf,sizeof(struct vector3));
  memcpy(&(mdlpvp->elmlp->special->corner2),temp_urb,sizeof(struct vector3));
  if (refine_cuboid(temp_llf,temp_urb,mdlpvp->pop->sb,mdlpvp->vol->grid_density))
  {
    mdlerror("Could not refine box to include new patch");
    return 1;
  }
  free(temp_llf);
  free(temp_urb);
  mdlpvp->elmlp->next=mdlpvp->element_list_head;
  mdlpvp->element_list_head=mdlpvp->elmlp;
};

patch_type: INCLUDE_PATCH { $$=0; }
	| EXCLUDE_PATCH { $$=1; }
;


in_obj_define_surface_regions: DEFINE_SURFACE_REGIONS '{'
	list_in_obj_surface_region_defs
	'}'
;


existing_obj_define_surface_regions: DEFINE_SURFACE_REGIONS '{'
	list_existing_obj_surface_region_defs
	'}'
;


list_existing_obj_surface_region_defs: existing_obj_surface_region_def
	| list_existing_obj_surface_region_defs existing_obj_surface_region_def
;


existing_obj_surface_region_def: existing_object
{
  mdlpvp->gp=$<sym>1;
  mdlpvp->objp=(struct object *)mdlpvp->gp->value;
  mdlpvp->obj_name=mdlpvp->objp->sym->name;
  if (mdlpvp->objp->object_type!=BOX_OBJ && mdlpvp->objp->object_type!=POLY_OBJ) {
    sprintf(mdlpvp->mdl_err_msg,"%s %s","Cannot define region on non-surface object:",mdlpvp->objp->sym->name);
    mdlerror(mdlpvp->mdl_err_msg,mdlpvp);
    return(1);
  }
  mdlpvp->pop=mdlpvp->objp->contents;
  mdlpvp->region_list_head=mdlpvp->objp->regions;
}
	'[' new_region ']' '{'
{
  mdlpvp->gp=$<sym>4;
  mdlpvp->rp=(struct region *)mdlpvp->gp->value;
  mdlpvp->rp->parent = mdlpvp->objp;  /* Was set to root of world, fixing */
  mdlpvp->eff_dat_head=mdlpvp->rp->eff_dat_head;
  mdlpvp->element_list_head = NULL;
}
	element_specifier_list
{
  mdlpvp->rp->element_list_head = mdlpvp->element_list_head;
  if (normalize_elements(mdlpvp->rp,1))
  {
    mdlerror("Improper element specification: out of range, out of memory,\n  or using patch specifier in already-created box");
    return 1;
  }
}
	list_opt_surface_region_stmts
        list_opt_surface_region_viz_values_stmts  
	'}'
{
  mdlpvp->rp->eff_dat_head=mdlpvp->eff_dat_head;
  mdlpvp->objp->regions=mdlpvp->region_list_head;
  mdlpvp->objp->num_regions++;
};


list_opt_surface_region_stmts: /* empty */
	| list_opt_surface_region_stmts opt_surface_region_stmt
;


list_opt_surface_region_viz_values_stmts: /* empty */ 
                                          |  surface_region_viz_value_stmt
;


opt_surface_region_stmt: set_surface_class_stmt
	| surface_mol_stmt
;


set_surface_class_stmt:	SURFACE_CLASS '=' existing_surface_class 
{
  mdlpvp->gp=$<sym>3;
  mdlpvp->rp->surf_class=(struct species *)mdlpvp->gp->value;

  if(mdlpvp->rp->surf_class->region_viz_value > 0){
     mdlpvp->rp->region_viz_value = mdlpvp->rp->surf_class->region_viz_value;
  }
};

surface_region_viz_value_stmt: VIZ_VALUE '=' num_expr
{
  /* if surface_class->region_viz_value is already defined print warning */
  if(mdlpvp->rp->surf_class != NULL){
     if(mdlpvp->rp->surf_class->region_viz_value > 0){
         sprintf(mdlpvp->mdl_err_msg,"ATTENTION: region_viz_value defined both through SURFACE_CLASS and VIZ_VALUE statements.\n");
         mdlerror(mdlpvp->mdl_err_msg,mdlpvp);
     }
   }
   /*this statement overwrites the results of the 'set_surface_class_stmt' */
   mdlpvp->rp->region_viz_value = (int)$<dbl>3;
};

list_in_obj_surface_region_defs: in_obj_surface_region_def
	| list_in_obj_surface_region_defs in_obj_surface_region_def
;


in_obj_surface_region_def: new_region '{'
{

  mdlpvp->gp=$<sym>1;
  mdlpvp->rp=(struct region *)mdlpvp->gp->value;
  mdlpvp->eff_dat_head=mdlpvp->rp->eff_dat_head;
  mdlpvp->element_list_head = NULL;
}
	element_specifier_list
{
  mdlpvp->rp->element_list_head = mdlpvp->element_list_head;
  if (mdlpvp->objp->object_type==POLY_OBJ)
  {
    if (normalize_elements(mdlpvp->rp,0))
    {
      mdlerror("Improper element specification: out of range or out of memory");
      return 1;
    }
  }
}
	list_opt_surface_region_stmts
        list_opt_surface_region_viz_values_stmts   
	'}'
{
  mdlpvp->pop=(struct polygon_object *)mdlpvp->objp->contents;
  mdlpvp->elmlp=mdlpvp->element_list_head;
  mdlpvp->objp->num_regions++;
  mdlpvp->rp->eff_dat_head=mdlpvp->eff_dat_head;
  if((mdlpvp->rp->surf_class != NULL) && (mdlpvp->rp->surf_class->region_viz_value > 0)){
      mdlpvp->rp->region_viz_value = mdlpvp->rp->surf_class->region_viz_value;
  }
};


new_region: VAR
{
  if (mdlpvp->cval_2!=NULL) {
    mdlpvp->sym_name=mdlpvp->cval_2;
  }
  else {
    mdlpvp->sym_name=mdlpvp->cval;
  }
  if ((mdlpvp->rp=make_new_region(volp,mdlpvp->obj_name,mdlpvp->sym_name,mdlpvp->mdl_err_msg))==NULL) {
    mdlerror(mdlpvp->mdl_err_msg,mdlpvp);
    return(1);
  }
  if ((mdlpvp->rlp=(struct region_list *)malloc(sizeof(struct region_list)))==NULL) {
    sprintf(mdlpvp->mdl_err_msg,"%s %s","Out of memory while creating region name:",mdlpvp->rp->sym->name);
    mdlerror(mdlpvp->mdl_err_msg,mdlpvp);
    return(1);
  }
  if (mdlpvp->sym_name==mdlpvp->cval) {
    mdlpvp->cval=NULL;
  }
  else {
    mdlpvp->cval_2=NULL;
  }
  mdlpvp->rp->region_last_name=mdlpvp->sym_name;
  mdlpvp->rp->parent=mdlpvp->curr_obj;
  mdlpvp->rp->surf_class=NULL;
  mdlpvp->rlp->reg=mdlpvp->rp;
  mdlpvp->rlp->next=mdlpvp->region_list_head;
  mdlpvp->region_list_head=mdlpvp->rlp;
  no_printf("Creating new region: %s\n",mdlpvp->rp->sym->name);
  fflush(stderr);
  $$=mdlpvp->rp->sym;
};


mod_surface_regions: MODIFY_SURFACE_REGIONS '{'
	list_existing_surface_region_refs
	'}'
;


list_existing_surface_region_refs: existing_surface_region_ref
	| list_existing_surface_region_refs existing_surface_region_ref
;


existing_surface_region_ref: existing_region '{'
{
  mdlpvp->gp=$<sym>1;
  mdlpvp->rp=(struct region *)mdlpvp->gp->value;
  mdlpvp->eff_dat_head=mdlpvp->rp->eff_dat_head;
}
	opt_surface_region_stmt
	list_opt_surface_region_stmts
        list_opt_surface_region_viz_values_stmts  
	'}'
{
  mdlpvp->rp->eff_dat_head=mdlpvp->eff_dat_head;
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
  mdlpvp->obj_name=mdlpvp->objp->sym->name;
  strncpy(mdlpvp->temp_str,"",1024);
  strncpy(mdlpvp->temp_str,mdlpvp->obj_name,1022);
  strcat(mdlpvp->temp_str,",");   
  mdlpvp->region_name=my_strcat(mdlpvp->temp_str,mdlpvp->sym_name);
  if(mdlpvp->region_name == NULL){
    sprintf(mdlpvp->mdl_err_msg,"%s ","Out of memory while parsing region");
    mdlerror(mdlpvp->mdl_err_msg);
    return(1);
  }
  if ((mdlpvp->gp=retrieve_sym(mdlpvp->region_name,REG,volp->main_sym_table))==NULL) {
    sprintf(mdlpvp->mdl_err_msg,"%s %s","Undefined region:",mdlpvp->region_name);
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
  free((void *)mdlpvp->region_name);
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
  mdlpvp->curr_obj->parent->n_walls+=mdlpvp->curr_obj->n_walls;
  mdlpvp->curr_obj->parent->n_walls_actual+=mdlpvp->curr_obj->n_walls_actual;
  mdlpvp->curr_obj->parent->n_verts+=mdlpvp->curr_obj->n_verts;
  mdlpvp->curr_obj=volp->root_object;
  if (mdlpvp->object_name_list_end->prev!=NULL) {
    mdlpvp->object_name_list_end=mdlpvp->object_name_list_end->prev;
  }
  else {
    mdlpvp->object_name_list_end->name=NULL;
  }
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


list_dashes: '-'
	| list_dashes '-'
;


right_arrow: list_dashes '>';
left_arrow: '<' list_dashes;
double_arrow: left_arrow '>';

right_cat_arrow: list_dashes existing_molecule_opt_orient right_arrow
{
  mdlpvp->pathp->reactant3 = (struct species*)( ($<sym>2)->value );
  if (mdlpvp->orient_specified) mdlpvp->pathp->orientation3 = mdlpvp->orient_class;
  else mdlpvp->pathp->orientation3 = ORIENT_NOT_SET;
};

double_cat_arrow: left_arrow existing_molecule_opt_orient right_arrow
{
  mdlpvp->pathp->reactant3 = (struct species*)( ($<sym>2)->value );
  if (mdlpvp->orient_specified) mdlpvp->pathp->orientation3 = mdlpvp->orient_class;
  else mdlpvp->pathp->orientation3 = ORIENT_NOT_SET;
};

unidir_arrow: right_arrow
	| right_cat_arrow;

bidir_arrow: double_arrow
	| double_cat_arrow;

reaction_arrow: unidir_arrow
{
  mdlpvp->bidirectional_arrow=0;
}
	| bidir_arrow
{
  mdlpvp->bidirectional_arrow=1;
};


new_rxn_pathname: /* empty */
{
  $$ = NULL;
}
	| ':' VAR
{
  if (mdlpvp->cval_2!=NULL) {
    mdlpvp->sym_name=mdlpvp->cval_2;
  }
  else {
    mdlpvp->sym_name=mdlpvp->cval;
  }
  if ((mdlpvp->gp=retrieve_sym(mdlpvp->sym_name,RXPN,volp->main_sym_table))==NULL) {
    if ((mdlpvp->gp=retrieve_sym(mdlpvp->sym_name,MOL,volp->main_sym_table))==NULL) {
      if ((mdlpvp->gp=store_sym(mdlpvp->sym_name,RXPN,volp->main_sym_table))==NULL) {
        sprintf(mdlpvp->mdl_err_msg,"%s %s","Out of memory while creating reaction name:",mdlpvp->sym_name);
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
      sprintf(mdlpvp->mdl_err_msg,"%s %s","Named reaction pathway already defined as a molecule:",mdlpvp->sym_name);
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
    sprintf(mdlpvp->mdl_err_msg,"%s %s","Named reaction pathway already defined:",mdlpvp->sym_name);
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


existing_rxpn_or_molecule: VAR
{

  if (mdlpvp->cval_2!=NULL) {  
    mdlpvp->sym_name=mdlpvp->cval_2;
  }   
  else {
    mdlpvp->sym_name=mdlpvp->cval;
  } 

  if ((mdlpvp->gp=retrieve_sym(mdlpvp->sym_name,RXPN,volp->main_sym_table))==NULL) {
    if ((mdlpvp->gp=retrieve_sym(mdlpvp->sym_name,MOL,volp->main_sym_table))==NULL) {
      sprintf(mdlpvp->mdl_err_msg,"%s %s","Undefined reaction or molecule name:",mdlpvp->sym_name);
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

existing_many_rxpns_or_molecules: WILDCARD_VAR
{
  int i;
  struct sym_table_list *stl;
  struct sym_table *sym_t;
  char *wildcard_string;

  if (mdlpvp->cval_2!=NULL) {  
    mdlpvp->sym_name=mdlpvp->cval_2;
  }   
  else {
    mdlpvp->sym_name=mdlpvp->cval;
  } 

  if(strip_quotes(mdlpvp->sym_name) != NULL){
      wildcard_string = strip_quotes(mdlpvp->sym_name);
  }else{
      sprintf(mdlpvp->mdl_err_msg,"%s ","Memory allocation error\n");
      mdlerror(mdlpvp->mdl_err_msg);
      free((void *)mdlpvp->sym_name);
    
      return (1);
  }

    /* here is a wildcard molecule/reaction name */
    /* do a full sym_table scan comparing each key */

    mdlpvp->sym_table_list_head = NULL;

    for(i = 0; i < SYM_HASHSIZE; i++)
    {
       for(sym_t = volp->main_sym_table[i]; sym_t != NULL; sym_t = sym_t->next)
       {
         if(is_wildcard_match(wildcard_string, sym_t->name)){ 

           if((sym_t->sym_type == MOL) || (sym_t->sym_type == RXPN))
           {
               stl = (struct sym_table_list *)mem_get(mdlpvp->sym_list_mem);
                if(stl == NULL){
                  mdlerror("Memory allocation error.\n", mdlpvp);
                  return 1;
                }
                stl->node = sym_t;
                stl->next = NULL;
                if(mdlpvp->sym_table_list_head == NULL){
                   mdlpvp->sym_table_list_head = stl;
                }else{
                   stl->next = mdlpvp->sym_table_list_head;
                   mdlpvp->sym_table_list_head = stl;
                }
            }
          }

         } /* end for */ 
    
    } /* end for */ 
   
  if (mdlpvp->sym_name==mdlpvp->cval) {
       mdlpvp->cval=NULL;
   }
   else {
       mdlpvp->cval_2=NULL;
   }
   free((void *)mdlpvp->sym_name);
   free(wildcard_string);
   
   mdlpvp->sym_table_list_head = sort_sym_list_by_name(mdlpvp->sym_table_list_head);

   /* set a flag about wildcard presence */
   mdlpvp->count_flags |= WILDCARD_PRESENT;

};

existing_one_or_many_rxpns_or_mols: existing_rxpn_or_molecule
{
  mdlpvp->sym_table_list_head = (struct sym_table_list*)mem_get(mdlpvp->sym_list_mem);
  if (mdlpvp->sym_table_list_head==NULL)
  {
    mdlerror("Out of memory trying to retrieve a molecule or reaction pathway name");
    return 1;
  }
  
  mdlpvp->sym_table_list_head->next=NULL;
  mdlpvp->sym_table_list_head->node=$<sym>1;
}
        | existing_many_rxpns_or_molecules;


rxn:
{
  mdlpvp->pathp = (struct pathway*)mem_get(mdlpvp->path_mem);
  if (mdlpvp->pathp==NULL)
  {
    mdlerror("Out of memory while creating reaction.");
    return 1;
  }
  mdlpvp->pathp->pathname = NULL;
  mdlpvp->pathp->reactant1=NULL;
  mdlpvp->pathp->reactant2=NULL;
  mdlpvp->pathp->reactant3=NULL;
  mdlpvp->pathp->km=0;
  mdlpvp->pathp->kcat=0;
  mdlpvp->pathp->km_filename=NULL;
}
  reactant_list reaction_arrow
{
  char *rx_name;
  int num_surfaces;
  int num_vol_mols;
  int oriented_count;
  
  rx_name = create_rx_name(mdlpvp->pathp);
  if (rx_name==NULL)
  {
    mdlerror("Out of memory while creating reaction.");
    return 1;
  }
  
  mdlpvp->gp = retrieve_sym(rx_name,RX,volp->main_sym_table);
  if (mdlpvp->gp==NULL)
  {
    mdlpvp->gp = store_sym(rx_name,RX,volp->main_sym_table);
    if (mdlpvp->gp==NULL)
    {
      mdlerror("Out of memory while creating reaction.");
      return 1;
    }
  }
  
  mdlpvp->rxnp = (struct rxn*)mdlpvp->gp->value;
  mdlpvp->rxnp->n_pathways++;
  mdlpvp->rxnp->n_reactants=1;
  if (mdlpvp->pathp->reactant2!=NULL) mdlpvp->rxnp->n_reactants++;
  if (mdlpvp->pathp->reactant3!=NULL) mdlpvp->rxnp->n_reactants++;
  
  mdlpvp->prod_all_3d=1;
  num_surfaces = 0;
  num_vol_mols = 0;  

  oriented_count = 0;
  if ((mdlpvp->pathp->reactant1->flags&NOT_FREE)!=0)
  {
    if (mdlpvp->pathp->reactant1->flags&IS_SURFACE) num_surfaces++;
    mdlpvp->prod_all_3d=0;
  }else{
    num_vol_mols++;
  }
  if (mdlpvp->pathp->orientation1 == ORIENT_NOT_SET) mdlpvp->pathp->orientation1 = 0;
  else oriented_count++;
  if (mdlpvp->pathp->reactant2!=NULL)
  {
    if ((mdlpvp->pathp->reactant2->flags&NOT_FREE)!=0)
    {
      if (mdlpvp->pathp->reactant2->flags&IS_SURFACE) num_surfaces++;
      mdlpvp->prod_all_3d=0;
    }else{
       num_vol_mols++;
    }
    if (mdlpvp->pathp->orientation2 == ORIENT_NOT_SET) mdlpvp->pathp->orientation2 = 0;
    else oriented_count++;      
  }
  if (mdlpvp->pathp->reactant3!=NULL)
  {
    if ((mdlpvp->pathp->reactant3->flags&NOT_FREE)!=0)
    {
      if (mdlpvp->pathp->reactant3->flags&IS_SURFACE) num_surfaces++;
      mdlpvp->prod_all_3d=0;
    }else{
      num_vol_mols++;
    }
    if (mdlpvp->pathp->orientation3 == ORIENT_NOT_SET) mdlpvp->pathp->orientation3 = 0;
    else oriented_count++;      
  }
  
  /* Check for invalid reaction specifications */
  if (num_surfaces>1)
  {
    mdlerror("Too many surfaces--reactions can take place on at most one surface.");
    return 1;
  }
  if (num_surfaces==mdlpvp->rxnp->n_reactants)
  {
    mdlerror("Reactants cannot consist entirely of surfaces.\n  Use a surface release site instead!");
    return 1;
  }
  if((num_vol_mols == 2) && (num_surfaces == 1)){
    mdlerror("Reaction between two volume molecules and a surface is not defined.\n");
    return 1;
  }
  if (mdlpvp->prod_all_3d)
  {
    if (oriented_count!=0)
    {
      if (mdlpvp->vol->notify->useless_vol_orient==WARN_ERROR)
      {
	mdlerror("Error: orientation specified for molecule in reaction in volume");
	return 1;
      }
      else if (mdlpvp->vol->notify->useless_vol_orient==WARN_WARN)
      {
	mdlerror("Warning: orientation specified for molecule in reaction in volume");
      }
    }
  }
  else
  {
    if (mdlpvp->rxnp->n_reactants != oriented_count)
    {
      if (mdlpvp->vol->notify->missed_surf_orient==WARN_ERROR)
      {
	mdlerror("Error: orientation not specified for molecule in reaction at surface\n  (use ; or ', or ,' for random orientation)");
	return 1;
      }
      else if (mdlpvp->vol->notify->missed_surf_orient==WARN_WARN)
      {
	mdlerror("Warning: orientation not specified for molecule in reaction at surface\n  (use ; or ', or ,' for random orientation)");
      }      
    }
  }
  if ((mdlpvp->rxnp->n_reactants==3) && (num_surfaces == 0))
  {
       mdlerror("Sorry, only reactions with two reactants and one surface are implemented yet.");
       return 1;
  }
 
  /* Copy catalyst to products */ 
  if ((mdlpvp->pathp->reactant3!=NULL) && (mdlpvp->pathp->reactant2 == NULL))
  {
    mdlpvp->prodp = (struct product*)mem_get(mdlpvp->prod_mem);
    if (mdlpvp->prodp==NULL)
    {
      mdlerror("Out of memory while creating reaction.");
      return 1;
    }
    
    mdlpvp->prodp->prod = mdlpvp->pathp->reactant3;
    if (mdlpvp->prod_all_3d) mdlpvp->prodp->orientation=0;
    else mdlpvp->prodp->orientation = mdlpvp->pathp->orientation3;
    mdlpvp->prodp->next = mdlpvp->pathp->product_head;
    mdlpvp->pathp->product_head = mdlpvp->prodp;
    
    if (mdlpvp->pathp->reactant2==NULL)
    {
      mdlpvp->pathp->reactant2 = mdlpvp->pathp->reactant3;
      mdlpvp->pathp->orientation2 = mdlpvp->pathp->orientation3;
      mdlpvp->pathp->reactant3 = NULL;
      mdlpvp->pathp->orientation3 = 0;
    }
  }
}
	list_products rx_rate_syntax new_rxn_pathname
{
  mdlpvp->gp=$<sym>7;
  if (mdlpvp->gp!=NULL)
  {
    mdlpvp->rxpnp=(struct rxn_pathname *)mdlpvp->gp->value;
    mdlpvp->rxpnp->rx=mdlpvp->rxnp;
    mdlpvp->pathp->pathname=mdlpvp->rxpnp;
  }
  mdlpvp->pathp->km=mdlpvp->fwd_km;
  mdlpvp->pathp->kcat=mdlpvp->fwd_kcat;
  if (mdlpvp->fwd_rate_filename != NULL)
  {
    struct pathway *tpp;
    
    mdlpvp->pathp->km_filename = mdlpvp->fwd_rate_filename;
    mdlpvp->fwd_rate_filename = NULL;
    
    if (mdlpvp->rxnp->pathway_head == NULL)
    {
      mdlpvp->rxnp->pathway_head = mdlpvp->pathp;
      mdlpvp->pathp->next = NULL;
    }
    else  /* Move varying reactions to the end of the list */
    {
      for ( tpp = mdlpvp->rxnp->pathway_head ; 
            tpp->next != NULL && tpp->next->km_filename==NULL ; 
            tpp = tpp->next ) {}
      mdlpvp->pathp->next = tpp->next;
      tpp->next = mdlpvp->pathp;
    }
  }
  else
  {
    mdlpvp->pathp->next=mdlpvp->rxnp->pathway_head;
    mdlpvp->rxnp->pathway_head=mdlpvp->pathp;
  }
  
  /* Create reverse reaction if we need to */
  if (mdlpvp->bidirectional_arrow)
  {
    if (mdlpvp->bidirectional_arrow==1) /* Hack to notice if we got both rates */
    {
      mdlerror("Reversible reaction indicated but no reverse rate supplied.");
      return 1;
    }
    if (invert_current_reaction_pathway(mdlpvp))
    {
      mdlerror("Error creating reverse reaction.");
      return 1;
    }
  }
  
};

reactant_list: reactant |
	reactant_list '+' reactant;

reactant: existing_molecule
{
  mdlpvp->orient_specified=0;
  mdlpvp->orient_class=0;
}
	orientation_class
{
  if (mdlpvp->pathp->reactant1==NULL)
  {
    mdlpvp->pathp->reactant1 = (struct species*)( ($<sym>1)->value );
    if (mdlpvp->orient_specified) mdlpvp->pathp->orientation1 = mdlpvp->orient_class;
    else mdlpvp->pathp->orientation1 = ORIENT_NOT_SET;
  }
  else if (mdlpvp->pathp->reactant2==NULL)
  {
    mdlpvp->pathp->reactant2 = (struct species*)( ($<sym>1)->value );
    if (mdlpvp->orient_specified) mdlpvp->pathp->orientation2 = mdlpvp->orient_class;
    else mdlpvp->pathp->orientation2 = ORIENT_NOT_SET;    
  }
  else if (mdlpvp->pathp->reactant3==NULL)
  {
    mdlpvp->pathp->reactant3 = (struct species*)( ($<sym>1)->value );
    if (mdlpvp->orient_specified) mdlpvp->pathp->orientation3 = mdlpvp->orient_class;
    else mdlpvp->pathp->orientation3 = ORIENT_NOT_SET;    
  }
  else
  {
    mdlerror("Too many reactants--maximum number is three.");
    return 1;
  }
};


list_products: product
	| list_products '+' product
;


product: existing_molecule
{
  mdlpvp->orient_specified=0;
  mdlpvp->orient_class=0;
}
	orientation_class
{
  mdlpvp->gp=$<sym>1;
  if ((mdlpvp->prodp=(struct product *)mem_get(mdlpvp->prod_mem))==NULL) {
    sprintf(mdlpvp->mdl_err_msg,"%s %s %s","Out of memory while creating reaction:",
      mdlpvp->rxnp->sym->name," -> ... ");
    mdlerror(mdlpvp->mdl_err_msg,mdlpvp);
    return(1);
  }
  mdlpvp->prodp->prod=(struct species *)mdlpvp->gp->value;
  if (mdlpvp->prod_all_3d) mdlpvp->prodp->orientation=0;
  else mdlpvp->prodp->orientation=mdlpvp->orient_class;

  mdlpvp->prodp->next=mdlpvp->pathp->product_head;
  mdlpvp->pathp->product_head=mdlpvp->prodp;

  if ((mdlpvp->prodp->prod->flags&IS_SURFACE)==0) /* Surfaces can be omitted, so ignore them */
  {
    if (mdlpvp->prod_all_3d==0)
    {
      if (!mdlpvp->orient_specified)
      {
	if (mdlpvp->vol->notify->missed_surf_orient==WARN_ERROR)
	{
	  mdlerror("Error: product orientation not specified in reaction with orientation\n  (use ; or ', or ,' for random orientation)");
	  return 1;
	}
	else if (mdlpvp->vol->notify->missed_surf_orient==WARN_WARN)
	{
	  mdlerror("Warning: product orientation not specified for molecule in reaction at surface\n  (use ; or ', or ,' for random orientation)");
	}
      }
    }
    else
    {
      if ((mdlpvp->prodp->prod->flags&NOT_FREE)!=0)
      {
	mdlerror("Reaction has only volume reactants but is trying to create a surface product");
	return 1;
      }
      if (mdlpvp->orient_specified)
      {
	if (mdlpvp->vol->notify->useless_vol_orient==WARN_ERROR)
	{
	  mdlerror("Error: orientation specified for molecule in reaction in volume");
	  return 1;
	}
	else if (mdlpvp->vol->notify->useless_vol_orient==WARN_WARN)
	{
	  mdlerror("Warning: orientation specified for molecule in reaction at surface");
	}
      }
    }
  }
}
	| NO_SPECIES
{
};


orientation_class: /* empty */
	| list_orient_marks
	| orient_class_number
	| ';'
{
  mdlpvp->orient_class=0;
  mdlpvp->orient_specified=1;
}
;

list_orient_marks:
        head_mark
	| tail_mark
        | list_orient_marks head_mark
	| list_orient_marks tail_mark;


head_mark: '\''
{
  mdlpvp->orient_class++;
  mdlpvp->orient_specified=1;
};


tail_mark: ','
{
  mdlpvp->orient_class--;
  mdlpvp->orient_specified=1;
};

orient_class_number: '{' num_expr '}'
{
  /* printf("!!oc = %d\n",(int)$<dbl>2); */
  mdlpvp->orient_class=(short)$<dbl>2;
  mdlpvp->orient_specified=1;
};


rx_rate_syntax:
{
  mdlpvp->fwd_km = mdlpvp->bkw_km = 0;
  mdlpvp->fwd_kcat = mdlpvp->bkw_kcat = 0;
  mdlpvp->fwd_rate_filename = mdlpvp->bkw_rate_filename = NULL;
}
	rx_rate1or2;
	
rx_rate1or2: rx_rate1
	| rx_rate2;

rx_rate1: '[' atomic_rate ']'
{
  if ( (mdlpvp->bkw_km!=0.0 || mdlpvp->bkw_rate_filename!=NULL) &&
       (mdlpvp->fwd_km==0.0 && mdlpvp->fwd_rate_filename==NULL) )
  {
    mdlerror("Invalid reaction rate specification: must specify a forward rate.");
    return 1;
  }
  if (mdlpvp->fwd_km<0)
  {
    mdlpvp->fwd_km=0.0;
    if (volp->notify->neg_reaction==WARN_ERROR)
    {
      mdlerror("Error: reaction rates should be zero or positive.");
      return 1;
    }
    else if (volp->notify->neg_reaction==WARN_WARN) mdlerror("Warning: negative reaction rate; setting to zero and continuing.");
  }
}

rx_rate2: '[' atomic_rate ',' atomic_rate ']'
{
  if (mdlpvp->fwd_km<0 || mdlpvp->bkw_km<0)
  {
    if (mdlpvp->fwd_km<0) mdlpvp->fwd_km=0.0;
    if (mdlpvp->bkw_km<0) mdlpvp->bkw_km=0.0;
    
    if (volp->notify->neg_reaction==WARN_ERROR)
    {
      mdlerror("Error: reaction rates should be zero or positive.");
      return 1;
    }
    else if (volp->notify->neg_reaction==WARN_WARN) mdlerror("Warning: negative reaction rate; setting to zero and continuing.");
  }
}

atomic_rate: num_expr_only
{
  if (mdlpvp->fwd_km!=0 || mdlpvp->fwd_rate_filename!=NULL)
  {
    mdlerror("Forward reaction rate already specified--can't specify again.");
    return 1;
  }
  mdlpvp->fwd_km = $<dbl>1;
}
	| '>' num_expr_only
{
  if (mdlpvp->fwd_km!=0 || mdlpvp->fwd_rate_filename!=NULL)
  {
    mdlerror("Forward reaction rate already specified--can't specify again.");
    return 1;
  }
  mdlpvp->fwd_km = $<dbl>2;
}
	| '<' num_expr_only
{
  if (mdlpvp->bkw_km!=0 || mdlpvp->bkw_rate_filename!=NULL)
  {
    mdlerror("Reverse reaction rate already specified--can't specify again.");
    return 1;
  }
  if (!mdlpvp->bidirectional_arrow)
  {
    mdlerror("Reverse rate specified but the reaction isn't reversible.");
    return 1;
  }
  mdlpvp->bidirectional_arrow++; /* Hack to indicate we've read a rate */
  mdlpvp->bkw_km = $<dbl>2;
}
	| str_expr_only
{
  if (mdlpvp->fwd_km!=0 || mdlpvp->fwd_rate_filename!=NULL)
  {
    mdlerror("Forward reaction rate already specified--can't specify again.");
    return 1;
  }
  mdlpvp->fwd_rate_filename = $<str>1;
}
	| '>' str_expr_only
{
  if (mdlpvp->fwd_km!=0 || mdlpvp->fwd_rate_filename!=NULL)
  {
    mdlerror("Forward reaction rate already specified--can't specify again.");
    return 1;
  }
  mdlpvp->fwd_rate_filename = $<str>2;
}
	| '<' str_expr_only
{
  if (mdlpvp->bkw_km!=0 || mdlpvp->bkw_rate_filename!=NULL)
  {
    mdlerror("Reverse reaction rate already specified--can't specify again.");
    return 1;
  }
  if (!mdlpvp->bidirectional_arrow)
  {
    mdlerror("Reverse rate specified but the reaction isn't reversible.");
    return 1;
  }
  mdlpvp->bidirectional_arrow++; /* Hack to indicate we've read a rate */
  mdlpvp->bkw_rate_filename = $<str>2;
}
	| existing_var_only
{
  struct sym_table *gp = $<sym>1;
  if (mdlpvp->fwd_km!=0 || mdlpvp->fwd_rate_filename!=NULL)
  {
    mdlerror("Forward reaction rate already specified--can't specify again.");
    return 1;
  }
  switch (gp->sym_type)
  {
    case DBL:
      mdlpvp->fwd_km = *(double *)gp->value;
      break;
    case STR:
      mdlpvp->fwd_rate_filename = my_strdup((char*)gp->value);
      if(mdlpvp->fwd_rate_filename == NULL)
      {
	mdlerror("Out of memory while storing reaction rate filename");
	return 1;
      }
      break;
    default:
      mdlerror("Invalid variable used for rates: must be number or filename");
      return(1);
      break;
  }
}
	| '>' existing_var_only
{
  struct sym_table *gp = $<sym>2;
  if (mdlpvp->fwd_km!=0 || mdlpvp->fwd_rate_filename!=NULL)
  {
    mdlerror("Forward reaction rate already specified--can't specify again.");
    return 1;
  }
  switch (gp->sym_type)
  {
    case DBL:
      mdlpvp->fwd_km = *(double *)gp->value;
      break;
    case STR:
      mdlpvp->fwd_rate_filename = my_strdup((char*)gp->value);
      if(mdlpvp->fwd_rate_filename == NULL)
      {
	mdlerror("Out of memory while storing reaction rate filename");
	return 1;
      }
      break;
    default:
      mdlerror("Invalid variable used for rates: must be number or filename");
      return(1);
      break;
  }
}	| '<' existing_var_only
{
  struct sym_table *gp = $<sym>2;
  if (mdlpvp->bkw_km!=0 || mdlpvp->bkw_rate_filename!=NULL)
  {
    mdlerror("Reverse reaction rate already specified--can't specify again.");
    return 1;
  }
  if (!mdlpvp->bidirectional_arrow)
  {
    mdlerror("Reverse rate specified but the reaction isn't reversible.");
    return 1;
  }
  mdlpvp->bidirectional_arrow++; /* Hack to indicate we've read a rate */
  switch (gp->sym_type)
  {
    case DBL:
      mdlpvp->bkw_km = *(double *)gp->value;
      break;
    case STR:
      mdlpvp->bkw_rate_filename = my_strdup((char*)gp->value);
      if(mdlpvp->bkw_rate_filename == NULL)
      {
	mdlerror("Out of memory while storing reaction rate filename");
	return 1;
      }
      break;
    default:
      mdlerror("Invalid variable used for rates: must be number or filename");
      return(1);
      break;
  }
};
  

viz_output_def: VIZ_OUTPUT '{'
	list_viz_output_cmds
	'}'
;

list_viz_output_cmds: viz_output_cmd
	| list_viz_output_cmds viz_output_cmd
;

viz_output_cmd: viz_mode_def
    	| viz_filename_prefix_def 
	| viz_molecules_block_def
	| viz_meshes_block_def 
	| viz_output_block_def
	| viz_iteration_frame_data_def
	| viz_molecule_prefix_def
	| viz_object_prefixes_def
	| viz_state_values_def
;

viz_filename_prefix_def: FILENAME '=' str_expr
{
  volp->file_prefix_name = $<str>3;

  if(volp->viz_mode == DX_MODE) {
      volp->molecule_prefix_name = $<str>3;
  }
};

viz_molecules_block_def: MOLECULES '{'
	list_viz_molecules_block_cmds
	'}'
;

list_viz_molecules_block_cmds: viz_molecules_block_cmd
	| list_viz_molecules_block_cmds viz_molecules_block_cmd
;

viz_molecules_block_cmd: 
		viz_molecules_name_list_cmd
        	| viz_molecules_time_points_def
        	| viz_molecules_iteration_numbers_def 
;

viz_molecules_name_list_cmd: NAME_LIST '{' 
		                viz_include_mols_cmd
                             '}'
;

viz_include_mols_cmd: list_mol_name_specs 
                      | list_mol_name_specs_state_values
                      | list_all_mols_specs
;

list_mol_name_specs: mol_name_spec 
	| list_mol_name_specs  mol_name_spec 
;

mol_name_spec: existing_one_or_multiple_molecules
{
  if(volp->viz_mode == DX_MODE){
    mdlerror("In DX MODE the state value for the molecule should be specified.\n");
    return(1);
  }

  struct sym_table_list *stl; 
 
  if( ($<sym>1) != NULL )
  {
     /* here is just one molecule name */
  
     mdlpvp->gp=$<sym>1;
     mdlpvp->specp=(struct species *)mdlpvp->gp->value;
     mdlpvp->specp->viz_state = INCLUDE_OBJ;

  }else{
      /* here are several molecules names because of using wildcards */ 
      
      stl = mdlpvp->sym_table_list_head;
      if (stl==NULL)
      {
	mdlerror("No molecules matching wildcard found.");
	return 1;
      }
      while(stl != NULL)
      {    
         mdlpvp->gp = stl->node;
         mdlpvp->specp=(struct species *)mdlpvp->gp->value;
         mdlpvp->specp->viz_state = INCLUDE_OBJ;

         stl = stl->next;
      }
  
      /* free allocated memory  */
      while (mdlpvp->sym_table_list_head != NULL)
      {
	stl=mdlpvp->sym_table_list_head->next;
	mem_put(mdlpvp->sym_list_mem,mdlpvp->sym_table_list_head);
	mdlpvp->sym_table_list_head=stl;
      }

  }

};

list_mol_name_specs_state_values: mol_name_specs_state_value
	| list_mol_name_specs_state_values  mol_name_specs_state_value
;

mol_name_specs_state_value: existing_molecule '=' num_expr
{
  mdlpvp->gp=$<sym>1;
  mdlpvp->specp=(struct species *)mdlpvp->gp->value;

  /* set the 'viz_state' value */ 
  mdlpvp->viz_state=(int) $<dbl>3;
  if (mdlpvp->gp->sym_type == MOL) {
      mdlpvp->existing_state=0;  
      mdlpvp->slp=mdlpvp->mol_state_head;
      while (mdlpvp->slp!=NULL && !mdlpvp->existing_state) {
        mdlpvp->existing_state=mdlpvp->existing_state||(mdlpvp->viz_state==mdlpvp->slp->state);
        mdlpvp->slp=mdlpvp->slp->next;
      }
      if (!mdlpvp->existing_state) {
        if ((mdlpvp->slp=(struct state_list *)malloc 
            (sizeof(struct state_list)))==NULL) {
          mdlerror("Out of memory while creating state list");
          return(1);
        }
        mdlpvp->slp->state=mdlpvp->viz_state;
        mdlpvp->slp->name=mdlpvp->gp->name;
        mdlpvp->slp->next=mdlpvp->mol_state_head;
        mdlpvp->mol_state_head=mdlpvp->slp;
      }
      mdlpvp->specp->viz_state=mdlpvp->viz_state;

  }
  if((volp->viz_output_flag & VIZ_MOLECULES_STATES) == 0){
      volp->viz_output_flag |= VIZ_MOLECULES_STATES;
  }
};

list_all_mols_specs: ALL_MOLECULES
{
  if(volp->viz_mode == DX_MODE){
    mdlerror("In DX MODE the state value for the molecule should be specified.\n");
    return(1);
  }
  volp->viz_output_flag |= VIZ_ALL_MOLECULES;
};

viz_molecules_time_points_def: TIME_POINTS '{'
                               viz_molecules_time_points_cmds
                          '}'
;

viz_molecules_time_points_cmds: viz_molecules_time_points_one_cmd
    | viz_molecules_time_points_cmds viz_molecules_time_points_one_cmd
;

viz_molecules_time_points_one_cmd: molecules_time_points_range_cmd 
                                   | molecules_time_points_all_times_cmd    
;

molecules_time_points_range_cmd: viz_molecules_one_item '@' '[' 
{
  mdlpvp->num_pos=0;
  mdlpvp->el_head=NULL;
  mdlpvp->el_tail=NULL;
}
     list_range_specs ']'
{
  sort_num_expr_list(mdlpvp->el_head);
  int temp = $<tok>1;

  if((volp->viz_mode == DREAMM_V3_GROUPED_MODE) || (volp->viz_mode == DREAMM_V3_MODE)){
     if ((mdlpvp->fdlp=(struct frame_data_list *)malloc
           (sizeof(struct frame_data_list)))==NULL) {
                mdlerror("Out of memory while creating time points");
                return(1);
      }
      mdlpvp->fdlp->list_type=OUTPUT_BY_TIME_LIST;
      mdlpvp->fdlp->type=$<tok>1;
      mdlpvp->fdlp->viz_iterationll=-1;
      mdlpvp->fdlp->n_viz_iterations=0;
      mdlpvp->fdlp->iteration_list=mdlpvp->el_head;
      mdlpvp->fdlp->curr_viz_iteration=mdlpvp->el_head;
      mdlpvp->fdlp->next=volp->frame_data_head;
      volp->frame_data_head = mdlpvp->fdlp;  
  }else if(volp->viz_mode == DX_MODE){ 
     if((temp == MOL_POS) || (temp == ALL_MOL_DATA)){
        /* create four frames */
        if ((mdlpvp->fdlp=(struct frame_data_list *)malloc
           (sizeof(struct frame_data_list)))==NULL) {
                mdlerror("Out of memory while creating time points");
                return(1);
         }
         mdlpvp->fdlp->list_type=OUTPUT_BY_TIME_LIST;
         mdlpvp->fdlp->type=EFF_POS;
         mdlpvp->fdlp->viz_iterationll=-1;
         mdlpvp->fdlp->n_viz_iterations=0;
         mdlpvp->fdlp->iteration_list=mdlpvp->el_head;
         mdlpvp->fdlp->curr_viz_iteration=mdlpvp->el_head;
         mdlpvp->fdlp->next=volp->frame_data_head;
         volp->frame_data_head = mdlpvp->fdlp;  

         if ((mdlpvp->fdlp=(struct frame_data_list *)malloc
           (sizeof(struct frame_data_list)))==NULL) {
                mdlerror("Out of memory while creating time points");
                return(1);
          }
          mdlpvp->fdlp->list_type=OUTPUT_BY_TIME_LIST;
          mdlpvp->fdlp->type=EFF_STATES;
          mdlpvp->fdlp->viz_iterationll=-1;
          mdlpvp->fdlp->n_viz_iterations=0;
          mdlpvp->fdlp->iteration_list=mdlpvp->el_head;
          mdlpvp->fdlp->curr_viz_iteration=mdlpvp->el_head;
          mdlpvp->fdlp->next=volp->frame_data_head;
          volp->frame_data_head = mdlpvp->fdlp;  

          if ((mdlpvp->fdlp=(struct frame_data_list *)malloc
                (sizeof(struct frame_data_list)))==NULL) {
                   mdlerror("Out of memory while creating time points");
                   return(1);
          }
          mdlpvp->fdlp->list_type=OUTPUT_BY_TIME_LIST;
          mdlpvp->fdlp->type=MOL_POS;
          mdlpvp->fdlp->viz_iterationll=-1;
          mdlpvp->fdlp->n_viz_iterations=0;
          mdlpvp->fdlp->iteration_list=mdlpvp->el_head;
          mdlpvp->fdlp->curr_viz_iteration=mdlpvp->el_head;
          mdlpvp->fdlp->next=volp->frame_data_head;
          volp->frame_data_head = mdlpvp->fdlp;  
     
          if ((mdlpvp->fdlp=(struct frame_data_list *)malloc
               (sizeof(struct frame_data_list)))==NULL) {
                    mdlerror("Out of memory while creating time points");
                    return(1);
          }
          mdlpvp->fdlp->list_type=OUTPUT_BY_TIME_LIST;
          mdlpvp->fdlp->type=MOL_STATES;
          mdlpvp->fdlp->viz_iterationll=-1;
          mdlpvp->fdlp->n_viz_iterations=0;
          mdlpvp->fdlp->iteration_list=mdlpvp->el_head;
          mdlpvp->fdlp->curr_viz_iteration=mdlpvp->el_head;
          mdlpvp->fdlp->next=volp->frame_data_head;
          volp->frame_data_head = mdlpvp->fdlp;  
     }else if(temp == MOL_ORIENT){
         /* do nothing */
     }
  }
};

molecules_time_points_all_times_cmd: viz_molecules_one_item '@' ALL_TIMES 
{
  mdlpvp->num_pos=0;
  mdlpvp->el_head=NULL;
  mdlpvp->el_tail=NULL;
  long long step;
  int temp;

  for(step = 0; step <= volp->iterations; step++)
  {
  	if ((mdlpvp->elp=(struct num_expr_list *)malloc
           (sizeof(struct num_expr_list)))==NULL) {
    		mdlerror("Out of memory while creating time points");
                return(1);
        }
 
        mdlpvp->elp->value=(double)(step*(volp->time_unit));
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
  sort_num_expr_list(mdlpvp->el_head); 
  temp = $<tok>1;

  if((volp->viz_mode == DREAMM_V3_GROUPED_MODE) || (volp->viz_mode == DREAMM_V3_MODE))
  {
     if ((mdlpvp->fdlp=(struct frame_data_list *)malloc
           (sizeof(struct frame_data_list)))==NULL) {
              mdlerror("Out of memory while creating time points");
              return(1);
     }
     mdlpvp->fdlp->list_type=OUTPUT_BY_TIME_LIST;
     mdlpvp->fdlp->type=$<tok>1;
     mdlpvp->fdlp->viz_iterationll=-1;
     mdlpvp->fdlp->n_viz_iterations=0;
     mdlpvp->fdlp->iteration_list=mdlpvp->el_head;
     mdlpvp->fdlp->curr_viz_iteration=mdlpvp->el_head;
     mdlpvp->fdlp->next=volp->frame_data_head;
     volp->frame_data_head = mdlpvp->fdlp;  
  }else if(volp->viz_mode == DX_MODE){ 
     if((temp == MOL_POS) || (temp == ALL_MOL_DATA)){
        /* create four frames */
        if ((mdlpvp->fdlp=(struct frame_data_list *)malloc
           (sizeof(struct frame_data_list)))==NULL) {
                mdlerror("Out of memory while creating time points");
                return(1);
         }
         mdlpvp->fdlp->list_type=OUTPUT_BY_TIME_LIST;
         mdlpvp->fdlp->type=EFF_POS;
         mdlpvp->fdlp->viz_iterationll=-1;
         mdlpvp->fdlp->n_viz_iterations=0;
         mdlpvp->fdlp->iteration_list=mdlpvp->el_head;
         mdlpvp->fdlp->curr_viz_iteration=mdlpvp->el_head;
         mdlpvp->fdlp->next=volp->frame_data_head;
         volp->frame_data_head = mdlpvp->fdlp;  

         if ((mdlpvp->fdlp=(struct frame_data_list *)malloc
           (sizeof(struct frame_data_list)))==NULL) {
                mdlerror("Out of memory while creating time points");
                return(1);
          }
          mdlpvp->fdlp->list_type=OUTPUT_BY_TIME_LIST;
          mdlpvp->fdlp->type=EFF_STATES;
          mdlpvp->fdlp->viz_iterationll=-1;
          mdlpvp->fdlp->n_viz_iterations=0;
          mdlpvp->fdlp->iteration_list=mdlpvp->el_head;
          mdlpvp->fdlp->curr_viz_iteration=mdlpvp->el_head;
          mdlpvp->fdlp->next=volp->frame_data_head;
          volp->frame_data_head = mdlpvp->fdlp;  

          if ((mdlpvp->fdlp=(struct frame_data_list *)malloc
                (sizeof(struct frame_data_list)))==NULL) {
                   mdlerror("Out of memory while creating time points");
                   return(1);
          }
          mdlpvp->fdlp->list_type=OUTPUT_BY_TIME_LIST;
          mdlpvp->fdlp->type=MOL_POS;
          mdlpvp->fdlp->viz_iterationll=-1;
          mdlpvp->fdlp->n_viz_iterations=0;
          mdlpvp->fdlp->iteration_list=mdlpvp->el_head;
          mdlpvp->fdlp->curr_viz_iteration=mdlpvp->el_head;
          mdlpvp->fdlp->next=volp->frame_data_head;
          volp->frame_data_head = mdlpvp->fdlp;  
     
          if ((mdlpvp->fdlp=(struct frame_data_list *)malloc
               (sizeof(struct frame_data_list)))==NULL) {
                    mdlerror("Out of memory while creating time points");
                    return(1);
          }
          mdlpvp->fdlp->list_type=OUTPUT_BY_TIME_LIST;
          mdlpvp->fdlp->type=MOL_STATES;
          mdlpvp->fdlp->viz_iterationll=-1;
          mdlpvp->fdlp->n_viz_iterations=0;
          mdlpvp->fdlp->iteration_list=mdlpvp->el_head;
          mdlpvp->fdlp->curr_viz_iteration=mdlpvp->el_head;
          mdlpvp->fdlp->next=volp->frame_data_head;
          volp->frame_data_head = mdlpvp->fdlp;  
     }else if(temp == MOL_ORIENT){
         /* do nothing */
     }
  }
};

viz_molecules_iteration_numbers_def: ITERATION_NUMBERS '{'
                                  viz_molecules_iteration_numbers_cmds
                                 '}'
;

viz_molecules_iteration_numbers_cmds: viz_molecules_iteration_numbers_one_cmd
  | viz_molecules_iteration_numbers_cmds viz_molecules_iteration_numbers_one_cmd
;

viz_molecules_iteration_numbers_one_cmd: molecules_iteration_numbers_range_cmd 
                            | molecules_iteration_numbers_all_iterations_cmd    
;

molecules_iteration_numbers_range_cmd: viz_molecules_one_item '@' '[' 
{
  mdlpvp->num_pos=0;
  mdlpvp->el_head=NULL;
  mdlpvp->el_tail=NULL;
}
     list_range_specs ']'
{
  sort_num_expr_list(mdlpvp->el_head);

  int temp = $<tok>1;

  if((volp->viz_mode == DREAMM_V3_GROUPED_MODE) || (volp->viz_mode == DREAMM_V3_MODE))
  {
     if ((mdlpvp->fdlp=(struct frame_data_list *)malloc
           (sizeof(struct frame_data_list)))==NULL) {
              mdlerror("Out of memory while creating iteration numbers");
              return(1);
      }
      mdlpvp->fdlp->list_type=OUTPUT_BY_ITERATION_LIST;
      mdlpvp->fdlp->type=$<tok>1;
      mdlpvp->fdlp->viz_iterationll=-1;
      mdlpvp->fdlp->n_viz_iterations=0;
      mdlpvp->fdlp->iteration_list=mdlpvp->el_head;
      mdlpvp->fdlp->curr_viz_iteration=mdlpvp->el_head;
      mdlpvp->fdlp->next=volp->frame_data_head;
      volp->frame_data_head = mdlpvp->fdlp;  
  }else if(volp->viz_mode == DX_MODE){ 
     if((temp == MOL_POS) || (temp == ALL_MOL_DATA)){
        /* create four frames */
        if ((mdlpvp->fdlp=(struct frame_data_list *)malloc
           (sizeof(struct frame_data_list)))==NULL) {
                mdlerror("Out of memory while creating time points");
                return(1);
         }
         mdlpvp->fdlp->list_type=OUTPUT_BY_ITERATION_LIST;
         mdlpvp->fdlp->type=EFF_POS;
         mdlpvp->fdlp->viz_iterationll=-1;
         mdlpvp->fdlp->n_viz_iterations=0;
         mdlpvp->fdlp->iteration_list=mdlpvp->el_head;
         mdlpvp->fdlp->curr_viz_iteration=mdlpvp->el_head;
         mdlpvp->fdlp->next=volp->frame_data_head;
         volp->frame_data_head = mdlpvp->fdlp;  

         if ((mdlpvp->fdlp=(struct frame_data_list *)malloc
           (sizeof(struct frame_data_list)))==NULL) {
                mdlerror("Out of memory while creating time points");
                return(1);
          }
          mdlpvp->fdlp->list_type=OUTPUT_BY_ITERATION_LIST;
          mdlpvp->fdlp->type=EFF_STATES;
          mdlpvp->fdlp->viz_iterationll=-1;
          mdlpvp->fdlp->n_viz_iterations=0;
          mdlpvp->fdlp->iteration_list=mdlpvp->el_head;
          mdlpvp->fdlp->curr_viz_iteration=mdlpvp->el_head;
          mdlpvp->fdlp->next=volp->frame_data_head;
          volp->frame_data_head = mdlpvp->fdlp;  

          if ((mdlpvp->fdlp=(struct frame_data_list *)malloc
                (sizeof(struct frame_data_list)))==NULL) {
                   mdlerror("Out of memory while creating time points");
                   return(1);
          }
          mdlpvp->fdlp->list_type=OUTPUT_BY_ITERATION_LIST;
          mdlpvp->fdlp->type=MOL_POS;
          mdlpvp->fdlp->viz_iterationll=-1;
          mdlpvp->fdlp->n_viz_iterations=0;
          mdlpvp->fdlp->iteration_list=mdlpvp->el_head;
          mdlpvp->fdlp->curr_viz_iteration=mdlpvp->el_head;
          mdlpvp->fdlp->next=volp->frame_data_head;
          volp->frame_data_head = mdlpvp->fdlp;  
     
          if ((mdlpvp->fdlp=(struct frame_data_list *)malloc
               (sizeof(struct frame_data_list)))==NULL) {
                    mdlerror("Out of memory while creating time points");
                    return(1);
          }
          mdlpvp->fdlp->list_type=OUTPUT_BY_ITERATION_LIST;
          mdlpvp->fdlp->type=MOL_STATES;
          mdlpvp->fdlp->viz_iterationll=-1;
          mdlpvp->fdlp->n_viz_iterations=0;
          mdlpvp->fdlp->iteration_list=mdlpvp->el_head;
          mdlpvp->fdlp->curr_viz_iteration=mdlpvp->el_head;
          mdlpvp->fdlp->next=volp->frame_data_head;
          volp->frame_data_head = mdlpvp->fdlp;  
     }else if(temp == MOL_ORIENT){
         /* do nothing */
     }
  }

};

molecules_iteration_numbers_all_iterations_cmd: viz_molecules_one_item '@' 
                                                 ALL_ITERATIONS 
{
  mdlpvp->num_pos=0;
  mdlpvp->el_head=NULL;
  mdlpvp->el_tail=NULL;
  long long step;
  int temp;

  for(step = 0; step <= volp->iterations; step++)
  {
  	if ((mdlpvp->elp=(struct num_expr_list *)malloc
           (sizeof(struct num_expr_list)))==NULL) {
    		mdlerror("Out of memory while creating iteration numbers");
                return(1);
        }
 
        mdlpvp->elp->value=(double)(step);
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
  sort_num_expr_list(mdlpvp->el_head); 
  temp = $<tok>1;

  if((volp->viz_mode == DREAMM_V3_GROUPED_MODE) || (volp->viz_mode == DREAMM_V3_MODE))
  {
     if ((mdlpvp->fdlp=(struct frame_data_list *)malloc
           (sizeof(struct frame_data_list)))==NULL) {
              mdlerror("Out of memory while creating iteration numbers");
              return(1);
     }
     mdlpvp->fdlp->list_type=OUTPUT_BY_ITERATION_LIST;
     mdlpvp->fdlp->type=$<tok>1;
     mdlpvp->fdlp->viz_iterationll=-1;
     mdlpvp->fdlp->n_viz_iterations=0;
     mdlpvp->fdlp->iteration_list=mdlpvp->el_head;
     mdlpvp->fdlp->curr_viz_iteration=mdlpvp->el_head;
     mdlpvp->fdlp->next=volp->frame_data_head;
     volp->frame_data_head = mdlpvp->fdlp;  
  }else if(volp->viz_mode == DX_MODE){ 
     if((temp == MOL_POS) || (temp == ALL_MOL_DATA)){
        /* create four frames */
        if ((mdlpvp->fdlp=(struct frame_data_list *)malloc
           (sizeof(struct frame_data_list)))==NULL) {
                mdlerror("Out of memory while creating time points");
                return(1);
         }
         mdlpvp->fdlp->list_type=OUTPUT_BY_ITERATION_LIST;
         mdlpvp->fdlp->type=EFF_POS;
         mdlpvp->fdlp->viz_iterationll=-1;
         mdlpvp->fdlp->n_viz_iterations=0;
         mdlpvp->fdlp->iteration_list=mdlpvp->el_head;
         mdlpvp->fdlp->curr_viz_iteration=mdlpvp->el_head;
         mdlpvp->fdlp->next=volp->frame_data_head;
         volp->frame_data_head = mdlpvp->fdlp;  

         if ((mdlpvp->fdlp=(struct frame_data_list *)malloc
           (sizeof(struct frame_data_list)))==NULL) {
                mdlerror("Out of memory while creating time points");
                return(1);
          }
          mdlpvp->fdlp->list_type=OUTPUT_BY_ITERATION_LIST;
          mdlpvp->fdlp->type=EFF_STATES;
          mdlpvp->fdlp->viz_iterationll=-1;
          mdlpvp->fdlp->n_viz_iterations=0;
          mdlpvp->fdlp->iteration_list=mdlpvp->el_head;
          mdlpvp->fdlp->curr_viz_iteration=mdlpvp->el_head;
          mdlpvp->fdlp->next=volp->frame_data_head;
          volp->frame_data_head = mdlpvp->fdlp;  

          if ((mdlpvp->fdlp=(struct frame_data_list *)malloc
                (sizeof(struct frame_data_list)))==NULL) {
                   mdlerror("Out of memory while creating time points");
                   return(1);
          }
          mdlpvp->fdlp->list_type=OUTPUT_BY_ITERATION_LIST;
          mdlpvp->fdlp->type=MOL_POS;
          mdlpvp->fdlp->viz_iterationll=-1;
          mdlpvp->fdlp->n_viz_iterations=0;
          mdlpvp->fdlp->iteration_list=mdlpvp->el_head;
          mdlpvp->fdlp->curr_viz_iteration=mdlpvp->el_head;
          mdlpvp->fdlp->next=volp->frame_data_head;
          volp->frame_data_head = mdlpvp->fdlp;  
     
          if ((mdlpvp->fdlp=(struct frame_data_list *)malloc
               (sizeof(struct frame_data_list)))==NULL) {
                    mdlerror("Out of memory while creating time points");
                    return(1);
          }
          mdlpvp->fdlp->list_type=OUTPUT_BY_ITERATION_LIST;
          mdlpvp->fdlp->type=MOL_STATES;
          mdlpvp->fdlp->viz_iterationll=-1;
          mdlpvp->fdlp->n_viz_iterations=0;
          mdlpvp->fdlp->iteration_list=mdlpvp->el_head;
          mdlpvp->fdlp->curr_viz_iteration=mdlpvp->el_head;
          mdlpvp->fdlp->next=volp->frame_data_head;
          volp->frame_data_head = mdlpvp->fdlp;  
     }else if(temp == MOL_ORIENT){
         /* do nothing */
     }
  }

};

viz_molecules_one_item: ALL_DATA
{
  $$=ALL_MOL_DATA;
}
	| POSITIONS
{
  $$=MOL_POS;
}
	| ORIENTATIONS
{
  $$=MOL_ORIENT;
};


viz_meshes_block_def: MESHES '{'
	list_viz_meshes_block_cmds
	'}'
;

list_viz_meshes_block_cmds: viz_meshes_block_cmd
	| list_viz_meshes_block_cmds viz_meshes_block_cmd
;

viz_meshes_block_cmd: 
		viz_meshes_name_list_cmd
        	| viz_meshes_time_points_def 
         	| viz_meshes_iteration_numbers_def  
;

viz_meshes_name_list_cmd: NAME_LIST '{' 
		                viz_include_meshes_cmd
                             '}'
;

viz_include_meshes_cmd: list_meshes_name_specs 
                      | list_meshes_name_specs_state_values 
                      | list_all_meshes_specs   
;

list_meshes_name_specs: mesh_one_name_spec
           | list_meshes_name_specs mesh_one_name_spec


mesh_one_name_spec: existing_one_or_multiple_objects
{
  if(volp->viz_mode == DX_MODE){
    mdlerror("In DX MODE the state value for the object should be specified.\n");
    return(1);
  }
  if(volp->file_prefix_name == NULL){
    mdlerror("The keyword FILENAME should be specified.\n");
    return(1);
  }
  
  u_int i;
  struct sym_table_list *stl; 
 
  if(($<sym>1) != NULL)
  {
     /* here is just one mesh object name */

     mdlpvp->gp=$<sym>1;
     // create viz_obj object 
     mdlpvp->objp=(struct object *)mdlpvp->gp->value;
     mdlpvp->sym_name=mdlpvp->gp->name;
     if ((mdlpvp->vizp = (struct viz_obj *)malloc
        (sizeof(struct viz_obj)))==NULL) {
           mdlerror("Out of memory while creating viz object");
           return(1);
     }
  
     mdlpvp->objp->viz_obj = mdlpvp->vizp;
     mdlpvp->vizp->name = my_strdup(volp->file_prefix_name);
     if(mdlpvp->vizp->name == NULL){
         sprintf(mdlpvp->mdl_err_msg,"%s %s","Out of memory while parsing object:",
          (char *)mdlpvp->sym_name);
         mdlerror(mdlpvp->mdl_err_msg,mdlpvp);
         return(1);
     }
  
     mdlpvp->vizp->full_name = my_strdup(mdlpvp->full_name);
     if(mdlpvp->vizp->full_name == NULL){
        sprintf(mdlpvp->mdl_err_msg,"%s %s","Out of memory while parsing object:",
          (char *)mdlpvp->sym_name);
        mdlerror(mdlpvp->mdl_err_msg,mdlpvp);
        return(1);
     }
  
     mdlpvp->vizp->viz_child_head = NULL;
     mdlpvp->vizp->obj=mdlpvp->objp;
     mdlpvp->vizp->next=volp->viz_obj_head;
     volp->viz_obj_head = mdlpvp->vizp;
 
     // set viz_state value of INCLUDE_OBJ for the object 
     if (mdlpvp->gp->sym_type == OBJ) {
         switch (mdlpvp->objp->object_type) {
            case META_OBJ:
              if (set_viz_state_value(mdlpvp->objp,INCLUDE_OBJ)) {
	         sprintf(mdlpvp->mdl_err_msg,"Out of memory while creating viz state value for meta object %s",mdlpvp->gp->name);
                 mdlerror(mdlpvp->mdl_err_msg,mdlpvp);
                 return(1);
              }
              break;
            case BOX_OBJ:
                mdlpvp->pop=(struct polygon_object *)mdlpvp->objp->contents;
                if (mdlpvp->objp->viz_state==NULL) {
                    if ((mdlpvp->objp->viz_state=(int *)malloc
                         (mdlpvp->pop->n_walls*sizeof(int)))==NULL) {
	                 sprintf(mdlpvp->mdl_err_msg,"Out of memory while creating viz state value for box object %s",mdlpvp->gp->name);
                         mdlerror(mdlpvp->mdl_err_msg,mdlpvp);
                         return(1);
                    }
                 }
                 for (i=0;i<mdlpvp->pop->n_walls;i++) {
                     mdlpvp->objp->viz_state[i]=INCLUDE_OBJ;
                 }
                 break;
             case POLY_OBJ:
                 mdlpvp->pop=(struct polygon_object *)mdlpvp->objp->contents;
                 if (mdlpvp->objp->viz_state==NULL) {
                     if ((mdlpvp->objp->viz_state=(int *)malloc
                        (mdlpvp->pop->n_walls*sizeof(int)))==NULL) {
	                   sprintf(mdlpvp->mdl_err_msg,"Out of memory while creating viz state value for polygon list object %s",mdlpvp->gp->name);
                           mdlerror(mdlpvp->mdl_err_msg,mdlpvp);
                           return(1);
                     }
                  }
                  for (i=0;i<mdlpvp->pop->n_walls;i++) {
                      mdlpvp->objp->viz_state[i]=INCLUDE_OBJ;
                  }
                  break;
             default:
                  mdlerror("Cannot set viz state value of this type of object");
                  return(1);
                  break;

         } /* end switch */
      }

   }else{

      /* here are several mesh objects names because of using wildcards */ 
      
      stl = mdlpvp->sym_table_list_head;
      if (stl==NULL)
      {
	mdlerror("No objects matching wildcard found.");
	return 1;
      }
      while(stl != NULL)
      {    

         mdlpvp->gp = stl->node;
         /* create viz_obj object */
         mdlpvp->objp=(struct object *)mdlpvp->gp->value;
         if((mdlpvp->objp->object_type == REL_SITE_OBJ) || 
             (mdlpvp->objp->object_type == META_OBJ)){ 
                  stl = stl->next;
                  continue;
         }  
         mdlpvp->sym_name=mdlpvp->gp->name;

         if ((mdlpvp->vizp = (struct viz_obj *)malloc
            (sizeof(struct viz_obj)))==NULL) {
                mdlerror("Out of memory while creating viz object");
                return(1);
         } 

         mdlpvp->objp->viz_obj = mdlpvp->vizp;
         mdlpvp->vizp->name = my_strdup(volp->file_prefix_name);
         if(mdlpvp->vizp->name == NULL){
             sprintf(mdlpvp->mdl_err_msg,"%s %s","Out of memory while parsing object:",
                (char *)mdlpvp->sym_name);
                mdlerror(mdlpvp->mdl_err_msg,mdlpvp);
                return(1);
         }
  
         mdlpvp->vizp->full_name = my_strdup(mdlpvp->sym_name); 
         if(mdlpvp->vizp->full_name == NULL){
             sprintf(mdlpvp->mdl_err_msg,"%s %s","Out of memory while parsing object:",
                  (char *)mdlpvp->sym_name);
                  mdlerror(mdlpvp->mdl_err_msg,mdlpvp);
                  return(1);
         }
  
         mdlpvp->vizp->viz_child_head = NULL;
         mdlpvp->vizp->obj=mdlpvp->objp;
         mdlpvp->vizp->next=volp->viz_obj_head;
         volp->viz_obj_head = mdlpvp->vizp;
 
         /* set viz_state value of INCLUDE_OBJ for the object */
         switch (mdlpvp->objp->object_type) {
               case BOX_OBJ:
                  mdlpvp->pop=(struct polygon_object *)mdlpvp->objp->contents;
                  if (mdlpvp->objp->viz_state==NULL) {
                    if ((mdlpvp->objp->viz_state=(int *)malloc
                         (mdlpvp->pop->n_walls*sizeof(int)))==NULL) {
	                 sprintf(mdlpvp->mdl_err_msg,"Out of memory while creating viz state value for box object %s",mdlpvp->gp->name);
                         mdlerror(mdlpvp->mdl_err_msg,mdlpvp);
                         return(1);
                    }
                   }
                   for (i=0;i<mdlpvp->pop->n_walls;i++) {
                      mdlpvp->objp->viz_state[i]=INCLUDE_OBJ;
                   }
                   break;
                case POLY_OBJ:
                   mdlpvp->pop=(struct polygon_object *)mdlpvp->objp->contents;
                   if (mdlpvp->objp->viz_state==NULL) {
                     if ((mdlpvp->objp->viz_state=(int *)malloc
                        (mdlpvp->pop->n_walls*sizeof(int)))==NULL) {
	                   sprintf(mdlpvp->mdl_err_msg,"Out of memory while creating viz state value for polygon list object %s",mdlpvp->gp->name);
                           mdlerror(mdlpvp->mdl_err_msg,mdlpvp);
                           return(1);
                     }
                   }
                   for (i=0;i<mdlpvp->pop->n_walls;i++) {
                      mdlpvp->objp->viz_state[i]=INCLUDE_OBJ;
                   }
                   break;
                default:
                    mdlerror("Cannot set viz state value of this type of object");
                    return(1);
                    break;

         } /* end switch */

         stl = stl->next;
      } /* end while(stl != NULL) */

      /* free allocated memory  */
      while (mdlpvp->sym_table_list_head != NULL)
      {
	stl=mdlpvp->sym_table_list_head->next;
	mem_put(mdlpvp->sym_list_mem,mdlpvp->sym_table_list_head);
	mdlpvp->sym_table_list_head=stl;
      }

   } /* end else */

};


list_meshes_name_specs_state_values: mesh_one_name_spec_state_value
    | list_meshes_name_specs_state_values mesh_one_name_spec_state_value


mesh_one_name_spec_state_value: existing_object '=' num_expr
{

  u_int i;

  mdlpvp->gp=$<sym>1;
  mdlpvp->viz_state=(int) $<dbl>3;
  /* create viz_obj object */
  mdlpvp->objp=(struct object *)mdlpvp->gp->value;
  mdlpvp->sym_name=mdlpvp->gp->name;
  if ((mdlpvp->vizp = (struct viz_obj *)malloc
        (sizeof(struct viz_obj)))==NULL) {
    mdlerror("Out of memory while creating viz object");
    return(1);
  }
  
  mdlpvp->objp->viz_obj = mdlpvp->vizp;
  mdlpvp->vizp->name = my_strdup(volp->file_prefix_name);
  if(mdlpvp->vizp->name == NULL){
    sprintf(mdlpvp->mdl_err_msg,"%s %s","Out of memory while parsing object:",
          (char *)mdlpvp->sym_name);
    mdlerror(mdlpvp->mdl_err_msg,mdlpvp);
    return(1);
  }
  
  mdlpvp->vizp->full_name = my_strdup(mdlpvp->full_name);
  if(mdlpvp->vizp->full_name == NULL){
    sprintf(mdlpvp->mdl_err_msg,"%s %s","Out of memory while parsing object:",
          (char *)mdlpvp->sym_name);
    mdlerror(mdlpvp->mdl_err_msg,mdlpvp);
    return(1);
  }
  
  mdlpvp->vizp->viz_child_head = NULL;
  mdlpvp->vizp->obj=mdlpvp->objp;
  mdlpvp->vizp->next=volp->viz_obj_head;
  volp->viz_obj_head = mdlpvp->vizp;

  if((volp->viz_output_flag & VIZ_SURFACE_STATES) == 0){
      volp->viz_output_flag |= VIZ_SURFACE_STATES;
  }
 
  /* set viz_state value for the object */
  if (mdlpvp->gp->sym_type == OBJ) {
      mdlpvp->existing_state=0;  
      mdlpvp->slp=mdlpvp->surf_state_head;
      while (mdlpvp->slp!=NULL && !mdlpvp->existing_state) {
        mdlpvp->existing_state=mdlpvp->existing_state||(mdlpvp->viz_state==mdlpvp->slp->state);
        mdlpvp->slp=mdlpvp->slp->next;
      }
      if (!mdlpvp->existing_state) {
        if ((mdlpvp->slp=(struct state_list *)malloc 
            (sizeof(struct state_list)))==NULL) {
          mdlerror("Out of memory while creating viz state list");
          return(1);
        }
        mdlpvp->slp->state=mdlpvp->viz_state;
        mdlpvp->slp->name=mdlpvp->gp->name;
        mdlpvp->slp->next=mdlpvp->surf_state_head;
        mdlpvp->surf_state_head=mdlpvp->slp;
      }
      switch (mdlpvp->objp->object_type) {
        case META_OBJ:
          if (set_viz_state_value(mdlpvp->objp,mdlpvp->viz_state)) {
	    sprintf(mdlpvp->mdl_err_msg,"Out of memory while creating viz state value for meta object %s",mdlpvp->gp->name);
            mdlerror(mdlpvp->mdl_err_msg,mdlpvp);
            return(1);
          }
          break;
        case BOX_OBJ:
          mdlpvp->pop=(struct polygon_object *)mdlpvp->objp->contents;
          if (mdlpvp->objp->viz_state==NULL) {
            if ((mdlpvp->objp->viz_state=(int *)malloc
                 (mdlpvp->pop->n_walls*sizeof(int)))==NULL) {
	      sprintf(mdlpvp->mdl_err_msg,"Out of memory while creating viz state value for box object %s",mdlpvp->gp->name);
              mdlerror(mdlpvp->mdl_err_msg,mdlpvp);
              return(1);
            }
          }
          for (i=0;i<mdlpvp->pop->n_walls;i++) {
            mdlpvp->objp->viz_state[i]=mdlpvp->viz_state;
          }
          break;
        case POLY_OBJ:
          mdlpvp->pop=(struct polygon_object *)mdlpvp->objp->contents;
          if (mdlpvp->objp->viz_state==NULL) {
            if ((mdlpvp->objp->viz_state=(int *)malloc
                 (mdlpvp->pop->n_walls*sizeof(int)))==NULL) {
	      sprintf(mdlpvp->mdl_err_msg,"Out of memory while creating viz state value for polygon list object %s",mdlpvp->gp->name);
              mdlerror(mdlpvp->mdl_err_msg,mdlpvp);
              return(1);
            }
          }
          for (i=0;i<mdlpvp->pop->n_walls;i++) {
            mdlpvp->objp->viz_state[i]=mdlpvp->viz_state;
          }
          break;
        default:
          mdlerror("Cannot set viz state value of this type of object");
          return(1);
          break;
      }
    }

};

list_all_meshes_specs: ALL_MESHES
{
  if(volp->viz_mode == DX_MODE){
    mdlerror("The keyword ALL_MESHES cannot be used in DX MODE.\n");
    return(1);
  }
  if(volp->file_prefix_name == NULL){
    mdlerror("The keyword FILENAME should be specified.\n");
    return(1);
  }

  u_int i;
  struct object *o;
  for(o = volp->root_instance->first_child; o != NULL; o = o->next)
  {
     mdlpvp->objp = o;
     mdlpvp->sym_name=o->sym->name;

     if ((mdlpvp->vizp = (struct viz_obj *)malloc
        (sizeof(struct viz_obj)))==NULL) {
        mdlerror("Out of memory while creating viz object");
      return(1);
     }
  
     mdlpvp->objp->viz_obj = mdlpvp->vizp;
     mdlpvp->vizp->name = my_strdup(volp->file_prefix_name);
     if(mdlpvp->vizp->name == NULL){
        sprintf(mdlpvp->mdl_err_msg,"%s %s","Out of memory while parsing object:",
          (char *)mdlpvp->sym_name);
        mdlerror(mdlpvp->mdl_err_msg,mdlpvp);
       return(1);
     }
  
     mdlpvp->vizp->full_name = my_strdup(o->sym->name);
     if(mdlpvp->vizp->full_name == NULL){
        sprintf(mdlpvp->mdl_err_msg,"%s %s","Out of memory while parsing object:",
          (char *)mdlpvp->sym_name);
        mdlerror(mdlpvp->mdl_err_msg,mdlpvp);
        return(1);
     }
  
     mdlpvp->vizp->viz_child_head = NULL;
     mdlpvp->vizp->obj=mdlpvp->objp;
     mdlpvp->vizp->next=volp->viz_obj_head;
     volp->viz_obj_head = mdlpvp->vizp;
  
     /* set viz_state value of INCLUDE_OBJ for the object */
         switch (mdlpvp->objp->object_type) {
           case META_OBJ:

             if (set_viz_state_value(mdlpvp->objp,INCLUDE_OBJ)) {
	       sprintf(mdlpvp->mdl_err_msg,"Out of memory while creating viz state value for meta object %s",mdlpvp->gp->name);
               mdlerror(mdlpvp->mdl_err_msg,mdlpvp);
               return(1);
             }
             break;
           case BOX_OBJ:

             mdlpvp->pop=(struct polygon_object *)mdlpvp->objp->contents;
             if (mdlpvp->objp->viz_state==NULL) {
                if ((mdlpvp->objp->viz_state=(int *)malloc
                    (mdlpvp->pop->n_walls*sizeof(int)))==NULL) {
	          sprintf(mdlpvp->mdl_err_msg,"Out of memory while creating viz state value for box object %s",mdlpvp->gp->name);
                  mdlerror(mdlpvp->mdl_err_msg,mdlpvp);
                  return(1);
                }
             }
             for (i=0;i<mdlpvp->pop->n_walls;i++) {
                mdlpvp->objp->viz_state[i]=INCLUDE_OBJ;
             }
             break;
           case POLY_OBJ:

             mdlpvp->pop=(struct polygon_object *)mdlpvp->objp->contents;
             if (mdlpvp->objp->viz_state==NULL) {
                if ((mdlpvp->objp->viz_state=(int *)malloc
                      (mdlpvp->pop->n_walls*sizeof(int)))==NULL) {
	        sprintf(mdlpvp->mdl_err_msg,"Out of memory while creating viz state value for polygon list object %s",mdlpvp->gp->name);
                mdlerror(mdlpvp->mdl_err_msg,mdlpvp);
                return(1);
              }
            }
            for (i=0;i<mdlpvp->pop->n_walls;i++) {
               mdlpvp->objp->viz_state[i]=INCLUDE_OBJ;
            }
            break;
          default:
             mdlerror("Cannot set viz state value of this type of object");
             return(1);
             break;
         }

  }


};

viz_meshes_time_points_def: TIME_POINTS '{'
                               viz_meshes_time_points_cmds
                          '}'
;

viz_meshes_time_points_cmds: viz_meshes_time_points_one_cmd
    | viz_meshes_time_points_cmds viz_meshes_time_points_one_cmd
;

viz_meshes_time_points_one_cmd: meshes_time_points_range_cmd 
                                   | meshes_time_points_all_times_cmd    
;

meshes_time_points_range_cmd: viz_meshes_one_item '@' '[' 
{
  mdlpvp->num_pos=0;
  mdlpvp->el_head=NULL;
  mdlpvp->el_tail=NULL;
}
     list_range_specs ']'
{
  sort_num_expr_list(mdlpvp->el_head);

  int temp = $<tok>1;
  if((volp->viz_mode == DREAMM_V3_GROUPED_MODE) || (volp->viz_mode == DREAMM_V3_MODE)){
     if ((mdlpvp->fdlp=(struct frame_data_list *)malloc
           (sizeof(struct frame_data_list)))==NULL) {
              mdlerror("Out of memory while creating time points");
              return(1);
     }
     mdlpvp->fdlp->list_type=OUTPUT_BY_TIME_LIST;
     mdlpvp->fdlp->type=$<tok>1;
     mdlpvp->fdlp->viz_iterationll=-1;
     mdlpvp->fdlp->n_viz_iterations=0;
     mdlpvp->fdlp->iteration_list=mdlpvp->el_head;
     mdlpvp->fdlp->curr_viz_iteration=mdlpvp->el_head;
     mdlpvp->fdlp->next=volp->frame_data_head;
     volp->frame_data_head = mdlpvp->fdlp; 
  }else if((volp->viz_mode == DX_MODE) && (temp == REG_DATA)){
      /* do nothing */
      fprintf(stderr, "REGION_DATA cannot be displayed in DX_MODE, please use \
                 DREAMM_V3_GROUPED (or DREAMM_V3) mode.\n");
  }else if(volp->viz_mode == DX_MODE){ 
     if((temp == MESH_GEOMETRY) || (temp == ALL_MESH_DATA)){
        /* create two frames - SURF_POS and SURF_STATES */
        if ((mdlpvp->fdlp=(struct frame_data_list *)malloc
              (sizeof(struct frame_data_list)))==NULL) {
                 mdlerror("Out of memory while creating time points");
                 return(1);
        }
        mdlpvp->fdlp->list_type=OUTPUT_BY_TIME_LIST;
        mdlpvp->fdlp->type=SURF_POS;
        mdlpvp->fdlp->viz_iterationll=-1;
        mdlpvp->fdlp->n_viz_iterations=0;
        mdlpvp->fdlp->iteration_list=mdlpvp->el_head;
        mdlpvp->fdlp->curr_viz_iteration=mdlpvp->el_head;
        mdlpvp->fdlp->next=volp->frame_data_head;
        volp->frame_data_head = mdlpvp->fdlp; 

        if ((mdlpvp->fdlp=(struct frame_data_list *)malloc
           (sizeof(struct frame_data_list)))==NULL) {
              mdlerror("Out of memory while creating time points");
              return(1);
        }
        mdlpvp->fdlp->list_type=OUTPUT_BY_TIME_LIST;
        mdlpvp->fdlp->type=SURF_STATES;
        mdlpvp->fdlp->viz_iterationll=-1;
        mdlpvp->fdlp->n_viz_iterations=0;
        mdlpvp->fdlp->iteration_list=mdlpvp->el_head;
        mdlpvp->fdlp->curr_viz_iteration=mdlpvp->el_head;
        mdlpvp->fdlp->next=volp->frame_data_head;
        volp->frame_data_head = mdlpvp->fdlp; 
     }
  }
};

meshes_time_points_all_times_cmd: viz_meshes_one_item '@' ALL_TIMES 
{
  mdlpvp->num_pos=0;
  mdlpvp->el_head=NULL;
  mdlpvp->el_tail=NULL;
  long long step;
  for(step = 0; step <= volp->iterations; step++)
  {
  	if ((mdlpvp->elp=(struct num_expr_list *)malloc
           (sizeof(struct num_expr_list)))==NULL) {
    		mdlerror("Out of memory while creating time points");
                return(1);
        }
 
        mdlpvp->elp->value=(double)(step*(volp->time_unit));
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
  sort_num_expr_list(mdlpvp->el_head); 
  
  int temp = $<tok>1;
  if((volp->viz_mode == DREAMM_V3_GROUPED_MODE) || (volp->viz_mode == DREAMM_V3_MODE)){
     
     if ((mdlpvp->fdlp=(struct frame_data_list *)malloc
        (sizeof(struct frame_data_list)))==NULL) {
           mdlerror("Out of memory while creating time points");
           return(1);
     }
     mdlpvp->fdlp->list_type=OUTPUT_BY_TIME_LIST;
     mdlpvp->fdlp->type=$<tok>1;
     mdlpvp->fdlp->viz_iterationll=-1;
     mdlpvp->fdlp->n_viz_iterations=0;
     mdlpvp->fdlp->iteration_list=mdlpvp->el_head;
     mdlpvp->fdlp->curr_viz_iteration=mdlpvp->el_head;
     mdlpvp->fdlp->next=volp->frame_data_head;
     volp->frame_data_head = mdlpvp->fdlp;  
  }else if((volp->viz_mode == DX_MODE) && (temp == REG_DATA)){
      /* do nothing */
      fprintf(stderr, "REGION_DATA cannot be dispalyed in DX_MODE, please use \
               DREAMM_V3_GROUPED (or DREAMM_V3) mode.\n");
  }else if(volp->viz_mode == DX_MODE){
     if((temp == MESH_GEOMETRY) || (temp == ALL_MESH_DATA))
     {
        /* create two frames - SURF_POS and SURF_STATES */
        if ((mdlpvp->fdlp=(struct frame_data_list *)malloc
           (sizeof(struct frame_data_list)))==NULL) {
              mdlerror("Out of memory while creating time points");
              return(1);
        }
        mdlpvp->fdlp->list_type=OUTPUT_BY_TIME_LIST;
        mdlpvp->fdlp->type=SURF_POS;
        mdlpvp->fdlp->viz_iterationll=-1;
        mdlpvp->fdlp->n_viz_iterations=0;
        mdlpvp->fdlp->iteration_list=mdlpvp->el_head;
        mdlpvp->fdlp->curr_viz_iteration=mdlpvp->el_head;
        mdlpvp->fdlp->next=volp->frame_data_head;
        volp->frame_data_head = mdlpvp->fdlp;  

        if ((mdlpvp->fdlp=(struct frame_data_list *)malloc
           (sizeof(struct frame_data_list)))==NULL) {
              mdlerror("Out of memory while creating time points");
              return(1);
        }
        mdlpvp->fdlp->list_type=OUTPUT_BY_TIME_LIST;
        mdlpvp->fdlp->type=SURF_STATES;
        mdlpvp->fdlp->viz_iterationll=-1;
        mdlpvp->fdlp->n_viz_iterations=0;
        mdlpvp->fdlp->iteration_list=mdlpvp->el_head;
        mdlpvp->fdlp->curr_viz_iteration=mdlpvp->el_head;
        mdlpvp->fdlp->next=volp->frame_data_head;
        volp->frame_data_head = mdlpvp->fdlp;  
    }
  }
};

viz_meshes_iteration_numbers_def: ITERATION_NUMBERS '{'
                                  viz_meshes_iteration_numbers_cmds
                                 '}'
;

viz_meshes_iteration_numbers_cmds: viz_meshes_iteration_numbers_one_cmd
  | viz_meshes_iteration_numbers_cmds viz_meshes_iteration_numbers_one_cmd
;

viz_meshes_iteration_numbers_one_cmd: meshes_iteration_numbers_range_cmd 
                            | meshes_iteration_numbers_all_iterations_cmd    
;

meshes_iteration_numbers_range_cmd: viz_meshes_one_item '@' '[' 
{
  mdlpvp->num_pos=0;
  mdlpvp->el_head=NULL;
  mdlpvp->el_tail=NULL;
}
     list_range_specs ']'
{
  sort_num_expr_list(mdlpvp->el_head);

  int temp = $<tok>1;
  if((volp->viz_mode == DREAMM_V3_GROUPED_MODE) || (volp->viz_mode == DREAMM_V3_MODE)){
     if ((mdlpvp->fdlp=(struct frame_data_list *)malloc
           (sizeof(struct frame_data_list)))==NULL) {
                mdlerror("Out of memory while creating iteration numbers");
                return(1);
     }
     mdlpvp->fdlp->list_type=OUTPUT_BY_ITERATION_LIST;
     mdlpvp->fdlp->type=$<tok>1;
     mdlpvp->fdlp->viz_iterationll=-1;
     mdlpvp->fdlp->n_viz_iterations=0;
     mdlpvp->fdlp->iteration_list=mdlpvp->el_head;
     mdlpvp->fdlp->curr_viz_iteration=mdlpvp->el_head;
     mdlpvp->fdlp->next=volp->frame_data_head;
     volp->frame_data_head = mdlpvp->fdlp;  
  }else if((volp->viz_mode == DX_MODE) && (temp == REG_DATA)){
      /* do nothing */
      fprintf(stderr, "REGION_DATA cannot be displayed in DX_MODE, please use \
       DREAMM_V3_GROUPED (or DREAMM_V3) mode.\n");   
  }else if(volp->viz_mode == DX_MODE){
     if((temp == MESH_GEOMETRY) || (temp == ALL_MESH_DATA)){
        /* create two frames - SURF_POS and SURF_STATES */
        if ((mdlpvp->fdlp=(struct frame_data_list *)malloc
           (sizeof(struct frame_data_list)))==NULL) {
                mdlerror("Out of memory while creating iteration numbers");
                return(1);
        }
        mdlpvp->fdlp->list_type=OUTPUT_BY_ITERATION_LIST;
        mdlpvp->fdlp->type=SURF_POS;
        mdlpvp->fdlp->viz_iterationll=-1;
        mdlpvp->fdlp->n_viz_iterations=0;
        mdlpvp->fdlp->iteration_list=mdlpvp->el_head;
        mdlpvp->fdlp->curr_viz_iteration=mdlpvp->el_head;
        mdlpvp->fdlp->next=volp->frame_data_head;
        volp->frame_data_head = mdlpvp->fdlp;  

        if ((mdlpvp->fdlp=(struct frame_data_list *)malloc
           (sizeof(struct frame_data_list)))==NULL) {
                mdlerror("Out of memory while creating iteration numbers");
                return(1);
        }
        mdlpvp->fdlp->list_type=OUTPUT_BY_ITERATION_LIST;
        mdlpvp->fdlp->type=SURF_STATES;
        mdlpvp->fdlp->viz_iterationll=-1;
        mdlpvp->fdlp->n_viz_iterations=0;
        mdlpvp->fdlp->iteration_list=mdlpvp->el_head;
        mdlpvp->fdlp->curr_viz_iteration=mdlpvp->el_head;
        mdlpvp->fdlp->next=volp->frame_data_head;
        volp->frame_data_head = mdlpvp->fdlp;  


     }

  }
};

meshes_iteration_numbers_all_iterations_cmd: viz_meshes_one_item '@' 
                                                 ALL_ITERATIONS 
{
  mdlpvp->num_pos=0;
  mdlpvp->el_head=NULL;
  mdlpvp->el_tail=NULL;
  long long step;
  for(step = 0; step <= volp->iterations; step++)
  {
  	if ((mdlpvp->elp=(struct num_expr_list *)malloc
           (sizeof(struct num_expr_list)))==NULL) {
    		mdlerror("Out of memory while creating iteration numbers");
                return(1);
        }
 
        mdlpvp->elp->value=(double)(step);
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
  sort_num_expr_list(mdlpvp->el_head);

  int temp = $<tok>1;
  if((volp->viz_mode == DREAMM_V3_GROUPED_MODE) || (volp->viz_mode == DREAMM_V3_MODE)){
 
     if ((mdlpvp->fdlp=(struct frame_data_list *)malloc
        (sizeof(struct frame_data_list)))==NULL) {
            mdlerror("Out of memory while creating iteration numbers");
            return(1);
     }
     mdlpvp->fdlp->list_type=OUTPUT_BY_ITERATION_LIST;
     mdlpvp->fdlp->type=$<tok>1;
     mdlpvp->fdlp->viz_iterationll=-1;
     mdlpvp->fdlp->n_viz_iterations=0;
     mdlpvp->fdlp->iteration_list=mdlpvp->el_head;
     mdlpvp->fdlp->curr_viz_iteration=mdlpvp->el_head;
     mdlpvp->fdlp->next=volp->frame_data_head;
     volp->frame_data_head = mdlpvp->fdlp;  
  }else if((volp->viz_mode == DX_MODE) && (temp == REG_DATA)){
      /* do nothing */
      fprintf(stderr, "REGION_DATA cannot be displayed in DX_MODE, please use \
         DREAMM_V3_GROUPED (or DREAMM_V3) mode.\n");
  }else if(volp->viz_mode == DX_MODE) {
     if((temp == MESH_GEOMETRY) || (temp == ALL_MESH_DATA))
     {
        /* create two frames - SURF_POS and SURF_STATES */

        if ((mdlpvp->fdlp=(struct frame_data_list *)malloc
           (sizeof(struct frame_data_list)))==NULL) {
               mdlerror("Out of memory while creating iteration numbers");
               return(1);
        }
        mdlpvp->fdlp->list_type=OUTPUT_BY_ITERATION_LIST;
        mdlpvp->fdlp->type=SURF_POS;
        mdlpvp->fdlp->viz_iterationll=-1;
        mdlpvp->fdlp->n_viz_iterations=0;
        mdlpvp->fdlp->iteration_list=mdlpvp->el_head;
        mdlpvp->fdlp->curr_viz_iteration=mdlpvp->el_head;
        mdlpvp->fdlp->next=volp->frame_data_head;
        volp->frame_data_head = mdlpvp->fdlp;  

        if ((mdlpvp->fdlp=(struct frame_data_list *)malloc
           (sizeof(struct frame_data_list)))==NULL) {
               mdlerror("Out of memory while creating iteration numbers");
               return(1);
        }
        mdlpvp->fdlp->list_type=OUTPUT_BY_ITERATION_LIST;
        mdlpvp->fdlp->type=SURF_STATES;
        mdlpvp->fdlp->viz_iterationll=-1;
        mdlpvp->fdlp->n_viz_iterations=0;
        mdlpvp->fdlp->iteration_list=mdlpvp->el_head;
        mdlpvp->fdlp->curr_viz_iteration=mdlpvp->el_head;
        mdlpvp->fdlp->next=volp->frame_data_head;
        volp->frame_data_head = mdlpvp->fdlp;  

     }
  }
};


viz_meshes_one_item: ALL_DATA
{
  $$=ALL_MESH_DATA;
}
	| GEOMETRY
{
  $$=MESH_GEOMETRY;
}
	| REGION_DATA
{
  $$=REG_DATA;
};


/* old viz_output style */
viz_data_output_def: VIZ_DATA_OUTPUT '{'
	list_viz_data_output_cmds
	'}'
;


list_viz_data_output_cmds: viz_data_output_cmd
	| list_viz_data_output_cmds viz_data_output_cmd
;


viz_data_output_cmd:
	viz_mode_def
	| voxel_image_mode_def
	| voxel_volume_mode_def
	| viz_output_block_def
	| viz_iteration_frame_data_def
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
	| MODE '=' DREAMM_V3
{
  volp->viz_mode = DREAMM_V3_MODE; 
}
	| MODE '=' DREAMM_V3_GROUPED
{
  volp->viz_mode = DREAMM_V3_GROUPED_MODE; 
}
	| MODE '=' CUSTOM_RK
{
  volp->viz_mode = RK_MODE;
  volp->rk_mode_var = NULL;
}
	| MODE '=' CUSTOM_RK '['
{
  mdlpvp->el_head=NULL;
  mdlpvp->el_tail=NULL;
  mdlpvp->num_pos=0;
}
	list_range_specs ']'
{
  struct num_expr_list *nel;
  int i;
  double *parts_array;
  int *bins_array;
  
  sort_num_expr_list(mdlpvp->el_head);
  
  volp->rk_mode_var = (struct rk_mode_data*)malloc(sizeof(struct rk_mode_data));
  parts_array = (double*)malloc(mdlpvp->num_pos * sizeof(double));
  bins_array = (int*)malloc( (mdlpvp->num_pos+1) * sizeof(int) );
  if (volp->rk_mode_var==NULL || parts_array==NULL || bins_array==NULL)
  {
    mdlerror("out of memory while setting up for visualization\n");
    return 1;
  }
  
  for (i=0,nel=mdlpvp->el_head ; nel!=NULL ; nel=nel->next,i++)
  {
    parts_array[i] = nel->value/volp->length_unit;
  }
  for ( ; mdlpvp->el_head != NULL ; mdlpvp->el_head = nel )
  {
    nel = mdlpvp->el_head->next;
    free(mdlpvp->el_head);
  }
  for (i=0;i<mdlpvp->num_pos+1;i++) bins_array[i] = 0;
  volp->viz_mode = RK_MODE;
  volp->rk_mode_var->n_bins = mdlpvp->num_pos+1;
  volp->rk_mode_var->bins = bins_array;
  volp->rk_mode_var->parts = parts_array;
  volp->rk_mode_var->n_written = 0;
}
	point
{
  volp->rk_mode_var->direction = $<vec3>9;
  if (vect_length(volp->rk_mode_var->direction)==0)
  {
    mdlerror("Cannot bin along a zero-length vector.\n");
    return 1;
  }
  normalize(volp->rk_mode_var->direction);
}
	| MODE '=' ASCII
{
  volp->viz_mode = ASCII_MODE;
};


voxel_image_mode_def: VOXEL_IMAGE_MODE '=' boolean
{
  volp->voxel_image_mode = $<tok>3;
};


voxel_volume_mode_def: VOXEL_VOLUME_MODE '=' boolean
{
  volp->voxel_volume_mode = $<tok>3;
};


viz_output_block_def: viz_iteration_def
	| viz_time_def 
;


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
    mdlerror("Out of memory while creating iteration list");
    return(1);
  }
  mdlpvp->fdlp->list_type=OUTPUT_BY_ITERATION_LIST;
  mdlpvp->fdlp->type=ALL_FRAME_DATA;
  mdlpvp->fdlp->viz_iterationll=-1;
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
    mdlerror("Out of memory while creating iteration list");
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
  for (mdlpvp->tmp_dbl=$<dbl>2;
       mdlpvp->tmp_dbl<$<dbl>4 || !distinguishable(mdlpvp->tmp_dbl,$<dbl>4,EPSILON);
       mdlpvp->tmp_dbl+=$<dbl>6)
  {
    if ((mdlpvp->elp=(struct num_expr_list *)malloc
        (sizeof(struct num_expr_list)))==NULL) {
      mdlerror("Out of memory while creating iteration list");
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
      mdlerror("Out of memory while creating iteration list");
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
    mdlerror("Out of memory while creating iteration list");
    return(1);
  }
  mdlpvp->fdlp->list_type=OUTPUT_BY_TIME_LIST;
  mdlpvp->fdlp->type=ALL_FRAME_DATA;
  mdlpvp->fdlp->viz_iterationll=-1;
  mdlpvp->fdlp->n_viz_iterations=0;
  mdlpvp->fdlp->iteration_list=mdlpvp->el_head;
  mdlpvp->fdlp->curr_viz_iteration=mdlpvp->el_head;
  mdlpvp->fdlp->next=volp->frame_data_head;
  volp->frame_data_head = mdlpvp->fdlp;  
};


viz_iteration_frame_data_def: ITERATION_FRAME_DATA '{'
	list_iteration_frame_data_specs
	'}'
;


list_iteration_frame_data_specs: iteration_frame_data_spec
	| list_iteration_frame_data_specs iteration_frame_data_spec
;


iteration_frame_data_spec: iteration_frame_data_item '='
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
    mdlerror("Out of memory while creating iteration frame data");
    return(1);
  }
  mdlpvp->fdlp->list_type=OUTPUT_BY_ITERATION_LIST;
  mdlpvp->fdlp->type=$<tok>1;
  mdlpvp->fdlp->viz_iterationll=-1;
  mdlpvp->fdlp->n_viz_iterations=0;
  mdlpvp->fdlp->iteration_list=mdlpvp->el_head;
  mdlpvp->fdlp->curr_viz_iteration=mdlpvp->el_head;
  mdlpvp->fdlp->next=volp->frame_data_head;
  volp->frame_data_head = mdlpvp->fdlp;
};


iteration_frame_data_item:  ALL_DATA
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


viz_object_prefixes_def: OBJECT_FILE_PREFIXES '{'
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
    mdlerror("Out of memory while creating viz object");
    return(1);
  }
  mdlpvp->objp->viz_obj = mdlpvp->vizp;
  mdlpvp->vizp->name = $<str>3;
  mdlpvp->vizp->full_name = my_strdup(mdlpvp->full_name);
  if(mdlpvp->vizp->full_name == NULL){
    sprintf(mdlpvp->mdl_err_msg,"%s %s","Out of memory while parsing object:",
          (char *)mdlpvp->full_name);
    mdlerror(mdlpvp->mdl_err_msg,mdlpvp);
    return(1);
  }
  mdlpvp->vizp->obj = mdlpvp->objp;
  mdlpvp->vizp->viz_child_head = NULL;
  mdlpvp->vizp->next = volp->viz_obj_head;
  volp->viz_obj_head = mdlpvp->vizp;
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
  u_int i;

  mdlpvp->gp=$<sym>1;
  mdlpvp->viz_state=(int) $<dbl>3;
  switch (mdlpvp->gp->sym_type) {
    case OBJ:
      mdlpvp->existing_state=0;  
      mdlpvp->slp=mdlpvp->surf_state_head;
      while (mdlpvp->slp!=NULL && !mdlpvp->existing_state) {
        mdlpvp->existing_state=mdlpvp->existing_state||(mdlpvp->viz_state==mdlpvp->slp->state);
        mdlpvp->slp=mdlpvp->slp->next;
      }
      if (!mdlpvp->existing_state) {
        if ((mdlpvp->slp=(struct state_list *)malloc 
            (sizeof(struct state_list)))==NULL) {
          mdlerror("Out of memory while creating state list");
          return(1);
        }
        mdlpvp->slp->state=mdlpvp->viz_state;
        mdlpvp->slp->name=mdlpvp->gp->name;
        mdlpvp->slp->next=mdlpvp->surf_state_head;
        mdlpvp->surf_state_head=mdlpvp->slp;
      }
      mdlpvp->objp=(struct object *)mdlpvp->gp->value;
      switch (mdlpvp->objp->object_type) {
        case META_OBJ:
          if (set_viz_state_value(mdlpvp->objp,mdlpvp->viz_state)) {
	    sprintf(mdlpvp->mdl_err_msg,"Out of memory while creating viz state value for meta object %s",mdlpvp->gp->name);
            mdlerror(mdlpvp->mdl_err_msg,mdlpvp);
            return(1);
          }
          break;
        case BOX_OBJ:
          mdlpvp->pop=(struct polygon_object *)mdlpvp->objp->contents;
          if (mdlpvp->objp->viz_state==NULL) {
            if ((mdlpvp->objp->viz_state=(int *)malloc
                 (mdlpvp->pop->n_walls*sizeof(int)))==NULL) {
	      sprintf(mdlpvp->mdl_err_msg,"Out of memory while creating viz state value for box object %s",mdlpvp->gp->name);
              mdlerror(mdlpvp->mdl_err_msg,mdlpvp);
              return(1);
            }
          }
          for (i=0;i<mdlpvp->pop->n_walls;i++) {
            mdlpvp->objp->viz_state[i]=mdlpvp->viz_state;
          }
          break;
        case POLY_OBJ:
          mdlpvp->pop=(struct polygon_object *)mdlpvp->objp->contents;
          if (mdlpvp->objp->viz_state==NULL) {
            if ((mdlpvp->objp->viz_state=(int *)malloc
                 (mdlpvp->pop->n_walls*sizeof(int)))==NULL) {
	      sprintf(mdlpvp->mdl_err_msg,"Out of memory while creating viz state value for polygon list object %s",mdlpvp->gp->name);
              mdlerror(mdlpvp->mdl_err_msg,mdlpvp);
              return(1);
            }
          }
          for (i=0;i<mdlpvp->pop->n_walls;i++) {
            mdlpvp->objp->viz_state[i]=mdlpvp->viz_state;
          }
          break;
        default:
          mdlerror("Cannot set viz state value of this type of object");
          return(1);
          break;
      }
      break;
    case MOL:
      mdlpvp->existing_state=0;  
      mdlpvp->slp=mdlpvp->mol_state_head;
      while (mdlpvp->slp!=NULL && !mdlpvp->existing_state) {
        mdlpvp->existing_state=mdlpvp->existing_state||(mdlpvp->viz_state==mdlpvp->slp->state);
        mdlpvp->slp=mdlpvp->slp->next;
      }
      if (!mdlpvp->existing_state) {
        if ((mdlpvp->slp=(struct state_list *)malloc 
            (sizeof(struct state_list)))==NULL) {
          mdlerror("Out of memory while creating state list");
          return(1);
        }
        mdlpvp->slp->state=mdlpvp->viz_state;
        mdlpvp->slp->name=mdlpvp->gp->name;
        mdlpvp->slp->next=mdlpvp->mol_state_head;
        mdlpvp->mol_state_head=mdlpvp->slp;
      }
      mdlpvp->specp=(struct species *)mdlpvp->gp->value;
      mdlpvp->specp->viz_state=mdlpvp->viz_state;
      break;
  }
}
	| existing_region '=' num_expr
{
  u_int i;

  mdlpvp->gp=$<sym>1;
  mdlpvp->viz_state=(int) $<dbl>3;

  mdlpvp->existing_state=0;  
  mdlpvp->slp=mdlpvp->surf_state_head;
  while (mdlpvp->slp!=NULL && !mdlpvp->existing_state) {
    mdlpvp->existing_state=mdlpvp->existing_state||(mdlpvp->viz_state==mdlpvp->slp->state);
    mdlpvp->slp=mdlpvp->slp->next;
  }
  if (!mdlpvp->existing_state) {
    if ((mdlpvp->slp=(struct state_list *)malloc 
        (sizeof(struct state_list)))==NULL) {
      mdlerror("Out of memory while creating state list");
      return(1);
    }
    mdlpvp->slp->state=mdlpvp->viz_state;
    mdlpvp->slp->name=mdlpvp->gp->name;
    mdlpvp->slp->next=mdlpvp->surf_state_head;
    mdlpvp->surf_state_head=mdlpvp->slp;
  }
  mdlpvp->objp=(struct object *)mdlpvp->gp->value;

  mdlpvp->pop=(struct polygon_object *)mdlpvp->objp->contents;
  if (mdlpvp->objp->viz_state==NULL) {
    if ((mdlpvp->objp->viz_state=(int *)malloc
        (mdlpvp->pop->n_walls*sizeof(int)))==NULL) {
      sprintf(mdlpvp->mdl_err_msg,"Out of memory while creating viz state value for elements of object %s",mdlpvp->gp->name);
      mdlerror(mdlpvp->mdl_err_msg,mdlpvp);
      return(1);
    }
    for (i=0;i<mdlpvp->pop->n_walls;i++) {
      mdlpvp->objp->viz_state[i]=EXCLUDE_OBJ;
    }
  }
  mdlpvp->elmlp=mdlpvp->element_list_head;
  while (mdlpvp->elmlp!=NULL) {
    if (mdlpvp->elmlp->begin < 0
        || mdlpvp->elmlp->end > mdlpvp->pop->n_walls-1) {
      mdlerror("Cannot set viz state value -- element out of range");
      return(1);
    }
    for (i=mdlpvp->elmlp->begin;i<=mdlpvp->elmlp->end;i++) {
      mdlpvp->objp->viz_state[i]=mdlpvp->viz_state;
    }
    mdlpvp->elmlp=mdlpvp->elmlp->next;
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


output_def: REACTION_DATA_OUTPUT '{'
{
  mdlpvp->obp = insert_new_output_block(mdlpvp);
  if (mdlpvp->obp==NULL)
  {
    mdlerror("Out of memory creating reaction data output");
    return 1;
  }

  mdlpvp->header_comment=NULL;
}
       output_buffer_size_def
{
  mdlpvp->obp->buffersize=$<dbl>4;
  mdlpvp->obp->time_array=(double*)malloc(mdlpvp->obp->buffersize*sizeof(double));
  if (mdlpvp->obp->time_array==NULL)
  {
    mdlerror("Out of memory creating output buffer.");
    return 1;
  }
}
       output_timer_def list_count_cmds '}';

output_buffer_size_def: OUTPUT_BUFFER_SIZE '=' num_expr
{
  double temp_value = $<dbl>3;
  if (!(temp_value >= 1.0 && temp_value < UINT_MAX))
  {
    mdlerror("Buffer size invalid.  Suggested range is 100-1000000.\n");
    return 1;
  }
  $$=$<dbl>3;
}
        | /* empty */
{
  $$=COUNTBUFFERSIZE;
}

output_timer_def:  step_time_def
{
  /* Pick a good buffer size */
  if (volp->chkpt_iterations) {
    mdlpvp->n_output=(long long)(volp->chkpt_iterations/mdlpvp->output_freq+1);
    mdlpvp->obp->buffersize=llmin3(volp->chkpt_iterations-volp->start_time+1,mdlpvp->n_output,mdlpvp->obp->buffersize); 
  }
  else {
    mdlpvp->n_output=(long long)(volp->iterations/mdlpvp->output_freq+1);
    mdlpvp->obp->buffersize=llmin3(volp->iterations-volp->start_time+1,mdlpvp->n_output,mdlpvp->obp->buffersize);
  }

  no_printf("Default output step time definition:\n");
  no_printf("  output step time = %g\n",mdlpvp->obp->step_time);
  no_printf("  output buffersize = %u\n",mdlpvp->obp->buffersize);
}
                | iteration_time_def 
                | real_time_def 
;


step_time_def: STEP '=' num_expr
{
  mdlpvp->obp->step_time=$<dbl>3;
  mdlpvp->obp->timer_type=OUTPUT_BY_STEP;

  mdlpvp->output_freq=mdlpvp->obp->step_time/volp->time_unit;
  
  if (mdlpvp->output_freq>volp->iterations && mdlpvp->output_freq>1)
  {
    mdlpvp->output_freq=(volp->iterations>1)?volp->iterations:1;
    mdlpvp->obp->step_time=mdlpvp->output_freq*volp->time_unit;
    sprintf(mdlpvp->mdl_err_msg,"Output step time too long\n\tSetting output step time to %g microseconds\n",mdlpvp->obp->step_time*1.0e6);
    mdl_warning(mdlpvp);
  }
  else if (mdlpvp->output_freq<1)
  {
    mdlpvp->output_freq=1;
    mdlpvp->obp->step_time=mdlpvp->output_freq*volp->time_unit;
    sprintf(mdlpvp->mdl_err_msg,"Output step time too short\n\tSetting output step time to %g microseconds\n",mdlpvp->obp->step_time*1.0e-6);
    mdl_warning(mdlpvp);
  }
}
        | /* empty */
{
  mdlpvp->output_freq=1;
  mdlpvp->obp->step_time=volp->time_unit;

  mdlpvp->obp->timer_type=OUTPUT_BY_STEP;
};


iteration_time_def: ITERATION_LIST '='
{
  mdlpvp->num_pos=0;
  mdlpvp->el_head=NULL;
  mdlpvp->el_tail=NULL;
}
	'[' list_range_specs ']'
{
  mdlpvp->n_output=mdlpvp->num_pos;
  mdlpvp->obp->timer_type=OUTPUT_BY_ITERATION_LIST;

  /* Pick a good buffer size */
  if (volp->chkpt_iterations) {
    mdlpvp->obp->buffersize=llmin3(volp->chkpt_iterations-volp->start_time+1,mdlpvp->n_output,mdlpvp->obp->buffersize);
  }
  else {
    mdlpvp->obp->buffersize=llmin3(volp->iterations-volp->start_time+1,mdlpvp->n_output,mdlpvp->obp->buffersize); 
  }

  sort_num_expr_list(mdlpvp->el_head);
  mdlpvp->obp->time_list_head=mdlpvp->el_head;
  mdlpvp->obp->time_now=NULL;
};


real_time_def: TIME_LIST '='
{
  mdlpvp->num_pos=0;
  mdlpvp->el_head=NULL;
  mdlpvp->el_tail=NULL;
}
	'[' list_range_specs ']'
{
  mdlpvp->n_output=mdlpvp->num_pos;
  mdlpvp->obp->timer_type=OUTPUT_BY_TIME_LIST;

  /* Pick a good buffer size */
  if (volp->chkpt_iterations) {
    mdlpvp->obp->buffersize=llmin3(volp->chkpt_iterations-volp->start_time+1,mdlpvp->n_output,mdlpvp->obp->buffersize);
  }
  else {
    mdlpvp->obp->buffersize=llmin3(volp->iterations-volp->start_time+1,mdlpvp->n_output,mdlpvp->obp->buffersize);
  }

  sort_num_expr_list(mdlpvp->el_head);
  mdlpvp->obp->time_list_head=mdlpvp->el_head;
  mdlpvp->obp->time_now=NULL;
};


list_count_cmds:
	count_cmd 
	| list_count_cmds count_cmd
	| custom_header
	| list_count_cmds custom_header;

count_cmd: '{'
{
  struct output_set *os;
  
  os = insert_new_output_set(mdlpvp->obp,mdlpvp->header_comment);
  if (os==NULL)
  {
    mdlerror("Out of memory while creating output list");
    return 1;
  }
  mdlpvp->count_flags = 0;
}
	list_count_exprs '}' file_arrow outfile_syntax
{
  struct output_set *os;
  struct output_column *oc,*new_head;
  int i;
  
  os = mdlpvp->obp->data_set_head;
  os->file_flags = $<tok>5;
  
  if ((mdlpvp->count_flags & (TRIGGER_PRESENT|COUNT_PRESENT)) == (TRIGGER_PRESENT|COUNT_PRESENT))
  {
    mdlerror("Cannot mix TRIGGER and COUNT statements.  Use separate files.");
    return 1;
  }
  
  i=check_reaction_output_file(os,mdlpvp->vol->err_file);
  if (i)
  {
    mdlerror("Unable to write data to output file.");
    return 1;
  }
  
  /* Reverse order so users see things in the order they asked for them */
  for (new_head=NULL ; os->column_head!=NULL ; os->column_head=oc)
  {
    oc = os->column_head->next;
    os->column_head->next = new_head;
    new_head = os->column_head;
  }
  os->column_head = new_head;
};

	
custom_header:
	HEADER '=' NONE
{
  mdlpvp->header_comment=NULL;
}
	| HEADER '=' boolean
{
  if ($<tok>3==0) mdlpvp->header_comment=NULL;
  else mdlpvp->header_comment="";
}
	| HEADER '=' str_expr
{
  mdlpvp->header_comment = $<str>3;
}


list_count_exprs:
	single_count_expr
	| list_count_exprs ',' single_count_expr
	
single_count_expr:
        count_expr opt_custom_header
{
  struct output_expression *oer,*oe;
  struct output_column *oc;
  struct output_set *os;
  int n;
  
  os = mdlpvp->obp->data_set_head;
  oer = $<cnt>1; /* Root of count expression */
  
  if (oer->oper==',' && ($<str>2)!=NULL)
  {
    mdlerror("Cannot use custom column headers with wildcard expansion");
    return 1;
  }
  
  if ($<str>2 != NULL) oer->title=$<str>2;
  
  /* If we have a list of results, go through to build column stack */
  for (oe=first_oexpr_tree(oer),n=0;oe!=NULL;oe=next_oexpr_tree(oe),n++)
  {
    oc = insert_new_output_column(os);
    if (oc==NULL)
    {
      mdlerror("Out of memory setting out output column");
      return 1;
    }
    
    oc->expr=oe;
    set_oexpr_column(oe,oc);
    
    switch (oe->expr_flags&OEXPR_TYPE_MASK)
    {
      case OEXPR_TYPE_INT:
        oc->data_type=INT;
        oc->buffer = (int*)malloc(os->block->buffersize*sizeof(int));
        break;
      case OEXPR_TYPE_DBL:
        oc->data_type=DBL;
        oc->buffer = (double*)malloc(os->block->buffersize*sizeof(double));
        break;
      case OEXPR_TYPE_TRIG:
        oc->data_type=TRIG_STRUCT;
        oc->buffer = (struct output_trigger_data*)malloc(os->block->buffersize*sizeof(struct output_trigger_data));
        break;
      default:
        mdlerror("Could not figure out what type of count data to store");
        return 1;
        break;
    }
    if (oc->buffer==NULL)
    {
      mdlerror("Out of memory creating buffer to store output data");
      return 1;
    }
  }
};

count_expr: num_value
{
  struct output_expression *oe = new_output_expr(mdlpvp->vol->oexpr_mem);
  if (oe==NULL)
  {
    mdlerror("Out of memory creating output expression");
    return 1;
  }
  oe->expr_flags=OEXPR_TYPE_DBL | OEXPR_TYPE_CONST;
  oe->value = $<dbl>1;
  oe->oper = '=';
  $$=oe;
}
        | count_value
{
  $$=$<cnt>1;
}
	| '(' count_expr ')' 
{
  $$=join_oexpr_tree(($<cnt>2),NULL,'(',mdlpvp->vol->oexpr_mem);
}
	| count_expr '+' count_expr
{
  $$=join_oexpr_tree(($<cnt>1),($<cnt>3) , '+',mdlpvp->vol->oexpr_mem);
}
        | count_expr '-' count_expr
{
  $$=join_oexpr_tree(($<cnt>1),($<cnt>3),'-',mdlpvp->vol->oexpr_mem);
}
        | count_expr '*' count_expr
{
  $$=join_oexpr_tree(($<cnt>1),($<cnt>3),'*',mdlpvp->vol->oexpr_mem);
}
        | count_expr '/' count_expr
{
  $$=join_oexpr_tree(($<cnt>1),($<cnt>3),'/',mdlpvp->vol->oexpr_mem);
}
        | '-' count_expr %prec UNARYMINUS
{
  $$=join_oexpr_tree(($<cnt>2),NULL,'_',mdlpvp->vol->oexpr_mem);
}
        | SUMMATION_OPERATOR '(' count_expr ')'
{
  oexpr_flood_convert($<cnt>3,',','+');
  eval_oexpr_tree($<cnt>3,0);
  $$=$<cnt>3;
};


count_value: COUNT 
{
  if (mdlpvp->count_flags&TRIGGER_PRESENT)
  {
    mdlerror("Can't put TRIGGER and COUNT statements in the same output file");
    return 1;
  }
  mdlpvp->count_flags |= COUNT_PRESENT;
}
        '[' count_syntax ']'
{
  $$=$<cnt>4;
}
        | EXPRESSION '[' num_expr ']'
{
  struct output_expression *oe = new_output_expr(mdlpvp->vol->oexpr_mem);
  if (oe==NULL)
  {
    mdlerror("Out of memory creating output expression");
    return 1;
  }
  oe->expr_flags=OEXPR_TYPE_DBL | OEXPR_TYPE_CONST;
  oe->value = $<dbl>3;
  oe->oper = '=';
  $$=oe;
}
	| TRIGGER
{
  if (mdlpvp->count_flags&COUNT_PRESENT)
  {
    mdlerror("Can't put TRIGGER and COUNT statements in the same output file");
    return 1;
  }
  mdlpvp->count_flags |= TRIGGER_PRESENT;
}
	'[' count_syntax ']'
{
  $$=$<cnt>4;
};



file_arrow: '>'
{
  $$=FILE_OVERWRITE;
}
	| '=' '>'
{
  $$=FILE_SUBSTITUTE;
}
	| '>' '>'
{
  $$=FILE_APPEND;
}
	| '>' '>' '>'
{
  $$=FILE_APPEND_HEADER;
}
	| '+' '>'
{
  $$=FILE_CREATE;
};

outfile_syntax: file_name
{
  mdlpvp->obp->data_set_head->outfile_name=$<str>1;
  no_printf("Counter output file set to %s\n",mdlpvp->obp->data_set_head->outfile_name);
};


count_syntax:  existing_one_or_many_rxpns_or_mols ',' count_location_specifier opt_hit_spec
{
  byte report_type,report_flags;
  struct sym_table_list *stl;
  struct output_expression *oe,*oet;
  struct output_expression *oe_head=NULL,*oe_tail=NULL;
  struct output_request *orq;
  struct output_request *or_head=NULL,*or_tail=NULL;

  if (mdlpvp->sym_table_list_head==NULL)
  {
    mdlerror("Wildcard matching found no matches for count output.");
    return 1;
  }
  
  if ($<sym>3==NULL)
  {
    report_flags=REPORT_WORLD;
    if ($<tok>4 != REPORT_NOTHING)
    {
      mdlerror("Invalid combination of WORLD with other counting options");
      return 1;
    }
  }
  else report_flags=0;
  
  if (mdlpvp->count_flags&TRIGGER_PRESENT) report_flags|=REPORT_TRIGGER;
  if ($<tok>4&REPORT_ENCLOSED) report_flags|=REPORT_ENCLOSED;
  
  for (stl=mdlpvp->sym_table_list_head ; stl!=NULL ; stl=stl->next)
  {
    if (stl->node->sym_type==MOL)
    {
      if (($<tok>4&REPORT_TYPE_MASK)==REPORT_NOTHING) report_type = REPORT_CONTENTS;
      else report_type=($<tok>4&REPORT_TYPE_MASK);
    }
    else
    {
      report_type = REPORT_RXNS;
      if (($<tok>4&REPORT_TYPE_MASK)!=REPORT_NOTHING)
      {
        char error_message[strlen(stl->node->name)+256];
        sprintf(error_message,"Invalid counting options used with reaction pathway %s",stl->node->name);
        mdlerror(error_message);
        return 1;
      }
    }
    
    orq = (struct output_request*)mem_get(mdlpvp->vol->outp_request_mem);
    oe = new_output_expr(mdlpvp->vol->oexpr_mem);
    if (orq==NULL || oe==NULL)
    {
      mdlerror("Out of memory storing requested counts");
      return 1;
    }
    orq->next=NULL;
    orq->requester=oe;
    orq->count_target=stl->node;
    orq->count_location=$<sym>3;
    orq->report_type=report_type|report_flags;
    
    oe->left=orq;
    oe->oper='#';
    oe->expr_flags=OEXPR_LEFT_REQUEST;
    if (orq->report_type&REPORT_TRIGGER) oe->expr_flags|=OEXPR_TYPE_TRIG;
    else if ((orq->report_type&REPORT_TYPE_MASK)!=REPORT_CONTENTS) oe->expr_flags|=OEXPR_TYPE_DBL;
    else oe->expr_flags|=OEXPR_TYPE_INT;
    
    if (oe_tail==NULL)
    {
      oe_head=oe_tail=oe;
      or_head=or_tail=orq;
    }
    else
    {
      or_tail->next=orq;
      or_tail=orq;
      
      oet = new_output_expr(mdlpvp->vol->oexpr_mem);
      if (oet==NULL)
      {
        mdlerror("Out of memory storing requested counts");
        return 1;
      }
      if (oe_tail->up==NULL)
      {
        oe_head=oet;
      }
      else
      {
        oet->up = oe_tail->up;
        if (oet->up->left==oe_tail) oet->up->left=oet;
        else oet->up->right=oet;
      }
      oet->left=oe_tail;
      oe_tail->up=oet;
      oet->right=oe;
      oe->up=oet;
      oet->oper=',';
      oet->expr_flags = OEXPR_LEFT_OEXPR|OEXPR_RIGHT_OEXPR|(oe->expr_flags&OEXPR_TYPE_MASK);
      oe_tail=oe;
    }
  }
  
  or_tail->next=mdlpvp->vol->output_request_head;
  mdlpvp->vol->output_request_head=or_head;
  
  /* free allocated memory */ 
  while (mdlpvp->sym_table_list_head != NULL)
  {
    stl=mdlpvp->sym_table_list_head->next;
    mem_put(mdlpvp->sym_list_mem,mdlpvp->sym_table_list_head);
    mdlpvp->sym_table_list_head=stl;
  }
  
  $$=oe_head;
};

count_location_specifier: WORLD
{
  $$=NULL;
}
        | existing_region
{
  $$=$<sym>1;
}
        | existing_object
{
  $$=$<sym>1;
};

opt_hit_spec:
{
  $$=REPORT_NOTHING;
}
        | ',' hit_spec
{
  $$=$<tok>2;
};

hit_spec: FRONT_HITS { $$ = REPORT_FRONT_HITS; }
	| BACK_HITS { $$ = REPORT_BACK_HITS; }
	| ALL_HITS { $$ = REPORT_ALL_HITS; }
	| FRONT_CROSSINGS { $$ = REPORT_FRONT_CROSSINGS; }
	| BACK_CROSSINGS { $$ = REPORT_BACK_CROSSINGS; }
	| ALL_CROSSINGS { $$ = REPORT_ALL_CROSSINGS; }
	| ESTIMATE_CONCENTRATION { $$ = REPORT_CONCENTRATION; }
	| ALL_ENCLOSED { $$ = REPORT_ENCLOSED; }
;

opt_custom_header:
{
  $$ = NULL;
}
	| ':' str_expr
{
  $$ = $<str>2;
};



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
  mdlpvp->filep=(struct file_stream *)mdlpvp->gp->value;
  mdlpvp->filep->name=$<str>4;
  mdlpvp->a_str=$<str>6;
  if ((mdlpvp->filep->stream=fopen(mdlpvp->filep->name,mdlpvp->a_str))==NULL) {
    sprintf(mdlpvp->mdl_err_msg,"%s %s","Cannot open file:",mdlpvp->filep->name);
    mdlerror(mdlpvp->mdl_err_msg,mdlpvp);
    return(1);
  }
};


new_file_stream: VAR
{
  if (mdlpvp->cval_2!=NULL) {
    mdlpvp->sym_name=mdlpvp->cval_2;
  }
  else {
    mdlpvp->sym_name=mdlpvp->cval;
  }

  if ((mdlpvp->gp=store_sym(mdlpvp->sym_name,FSTRM,volp->main_sym_table))==NULL) {
    sprintf(mdlpvp->mdl_err_msg,"%s %s","Out of memory while creating file stream:",mdlpvp->sym_name);
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


file_mode: str_expr
{
  char c;

  mdlpvp->a_str=$<str>1;
  c=mdlpvp->a_str[0];
  if (c!='r' 
      && c!='w'
      && c!='a') {
    sprintf(mdlpvp->mdl_err_msg,"%s %s","Invalid file mode:",mdlpvp->a_str);
    mdlerror(mdlpvp->mdl_err_msg,mdlpvp);
    free((void *)mdlpvp->cval);
    mdlpvp->cval=NULL;
    return(1);
  }
  $$=mdlpvp->a_str;
};


fclose_stmt: FCLOSE '(' existing_file_stream ')'
{
  mdlpvp->gp=$<sym>3;
  mdlpvp->filep=(struct file_stream *)mdlpvp->gp->value;
  if (fclose(mdlpvp->filep->stream)!=0) {
    sprintf(mdlpvp->mdl_err_msg,"%s %s","Error closing file:",mdlpvp->filep->name);
    mdlerror(mdlpvp->mdl_err_msg,mdlpvp);
    return(1);
  }
};


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
};


printf_stmt: PRINTF arg_list_init '(' format_string ',' list_args ')'
{
  mdlpvp->a_str=$<str>4;
  if (my_fprintf(stderr,mdlpvp->a_str,mdlpvp->arg_list,mdlpvp->num_args)) {
    sprintf(mdlpvp->mdl_err_msg,"%s %s","Could not print to stderr:",mdlpvp->a_str);
    mdlerror(mdlpvp->mdl_err_msg,mdlpvp);
    return(1);
  }
}
	| PRINTF arg_list_init '(' format_string ')'
{
  mdlpvp->a_str=$<str>4;
  if (my_fprintf(stderr,mdlpvp->a_str,NULL,mdlpvp->num_args)) {
    sprintf(mdlpvp->mdl_err_msg,"%s %s","Could not print to stderr:",mdlpvp->a_str);
    mdlerror(mdlpvp->mdl_err_msg,mdlpvp);
    return(1);
  }
};


arg_list_init: /* empty */
{
	mdlpvp->num_args=0;
};


format_string: str_expr
{
  char *rem_str;
  char c;
  int pos1;

  mdlpvp->a_str=$<str>1;
  rem_str=mdlpvp->a_str;
  strcpy(mdlpvp->format_str,"");
  while(rem_str!=NULL) {
    pos1=strcspn(rem_str,"\\");
    if(pos1==strlen(rem_str)) {  /* no \ found */
      strcat(mdlpvp->format_str,rem_str);
      rem_str=NULL;
    }
    else {  /* found a \ */
      strncat(mdlpvp->format_str,rem_str,pos1);
      c=rem_str[pos1+1];
      switch(c) {
      case 'n':
        strcat(mdlpvp->format_str,"\n");
        break;
      case 't':
        strcat(mdlpvp->format_str,"\t");
        break;
      case '\\':
        strcat(mdlpvp->format_str,"\\");
        break;
      case '\"':
        strcat(mdlpvp->format_str,"\"");
        break;
      }
      rem_str=rem_str+pos1+2;
    }
  }
  $$=mdlpvp->format_str;
};


list_args: num_expr_only
{
        if ((mdlpvp->arg_list[mdlpvp->num_args].arg_value=(void *)double_dup($<dbl>1))==NULL) {
          sprintf(mdlpvp->mdl_err_msg,"%s","Out of memory while creating format argument");
          mdlerror(mdlpvp->mdl_err_msg,mdlpvp);
          return(1);
        }
        mdlpvp->arg_list[mdlpvp->num_args++].arg_type=DBL;
}
	| list_args ',' num_expr_only
{
        if ((mdlpvp->arg_list[mdlpvp->num_args].arg_value=(void *)double_dup($<dbl>3))==NULL) {
          sprintf(mdlpvp->mdl_err_msg,"%s","Out of memory while creating format argument");
          mdlerror(mdlpvp->mdl_err_msg,mdlpvp);
          return(1);
        }
        mdlpvp->arg_list[mdlpvp->num_args++].arg_type=DBL;
}
	| str_expr_only
{
        mdlpvp->arg_list[mdlpvp->num_args].arg_value=(void *)$<str>1;
        mdlpvp->arg_list[mdlpvp->num_args++].arg_type=STR;
}
	| list_args ',' str_expr_only
{
        mdlpvp->arg_list[mdlpvp->num_args].arg_value=(void *)$<str>3;
        mdlpvp->arg_list[mdlpvp->num_args++].arg_type=STR;
}
	| existing_var_only
{
  mdlpvp->gp=$<sym>1;
  mdlpvp->arg_list[mdlpvp->num_args].arg_value=mdlpvp->gp->value;
  switch (mdlpvp->gp->sym_type) {
  case DBL:
    mdlpvp->arg_list[mdlpvp->num_args++].arg_type=DBL;
    break;
  case STR:
    mdlpvp->arg_list[mdlpvp->num_args++].arg_type=STR;
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
  mdlpvp->arg_list[mdlpvp->num_args].arg_value=mdlpvp->gp->value;
  switch (mdlpvp->gp->sym_type) {
  case DBL:
    mdlpvp->arg_list[mdlpvp->num_args++].arg_type=DBL;
    break;
  case STR:
    mdlpvp->arg_list[mdlpvp->num_args++].arg_type=STR;
    break;
  default:
    mdlerror("Invalid variable type referenced");
    return(1);
    break;
  }
};


fprintf_stmt: FPRINTF arg_list_init '(' existing_file_stream ',' format_string ',' list_args ')' 
{
  mdlpvp->gp=$<sym>4;
  mdlpvp->filep=(struct file_stream *)mdlpvp->gp->value;
  mdlpvp->a_str=$<str>6;
  if (my_fprintf(mdlpvp->filep->stream,mdlpvp->a_str,mdlpvp->arg_list,mdlpvp->num_args)) {
    sprintf(mdlpvp->mdl_err_msg,"%s %s","Could not print to file:",mdlpvp->filep->name);
    mdlerror(mdlpvp->mdl_err_msg,mdlpvp);
    return(1);
  }
}
	| FPRINTF arg_list_init '(' existing_file_stream ',' format_string ')'
{
  mdlpvp->gp=$<sym>4;
  mdlpvp->filep=(struct file_stream *)mdlpvp->gp->value;
  mdlpvp->a_str=$<str>6;
  if (my_fprintf(mdlpvp->filep->stream,mdlpvp->a_str,NULL,mdlpvp->num_args)) {
    sprintf(mdlpvp->mdl_err_msg,"%s %s","Could not print to file:",mdlpvp->filep->name);
    mdlerror(mdlpvp->mdl_err_msg,mdlpvp);
    return(1);
  }
};


sprintf_stmt: SPRINTF arg_list_init '(' assign_var ',' format_string ',' list_args ')'
{ 
  mdlpvp->gp=$<sym>4;
  mdlpvp->a_str=$<str>6;
  if (my_sprintf(mdlpvp->str_buf2,mdlpvp->a_str,mdlpvp->arg_list,mdlpvp->num_args)) {
    sprintf(mdlpvp->mdl_err_msg,"%s %s","Could not sprintf to variable:",mdlpvp->gp->name);
    mdlerror(mdlpvp->mdl_err_msg,mdlpvp);
    return(1);
  }
  mdlpvp->gp->sym_type=STR;
  mdlpvp->gp->value=(void *)my_strdup(mdlpvp->str_buf2);
  if(mdlpvp->gp->value == NULL){
    sprintf(mdlpvp->mdl_err_msg,"%s %s","Out of memory while parsing string:",
          (char *)mdlpvp->str_buf2);
    mdlerror(mdlpvp->mdl_err_msg,mdlpvp);
    return(1);
  }
}
        | SPRINTF arg_list_init '(' assign_var ',' format_string ')'
{
  mdlpvp->gp=$<sym>4;
  mdlpvp->a_str=$<str>6;
  if (my_sprintf(mdlpvp->str_buf,mdlpvp->a_str,NULL,mdlpvp->num_args)) {
    sprintf(mdlpvp->mdl_err_msg,"%s %s","Could not sprintf to variable:",mdlpvp->gp->name);
    mdlerror(mdlpvp->mdl_err_msg,mdlpvp);
    return(1);
  }
  mdlpvp->gp->sym_type=STR;
  mdlpvp->gp->value=(void *)my_strdup(mdlpvp->str_buf);
  if(mdlpvp->gp->value == NULL){
    sprintf(mdlpvp->mdl_err_msg,"%s %s","Out of memory while parsing string:",
          mdlpvp->str_buf);
    mdlerror(mdlpvp->mdl_err_msg,mdlpvp);
    return(1);
  }
};

print_time_stmt: PRINT_TIME '(' format_string ')'
{
  time_t the_time;

  mdlpvp->a_str=$<str>3;
  the_time=time(NULL);
  strftime(mdlpvp->time_str,128,mdlpvp->a_str,localtime(&the_time));
  if (volp->procnum == 0) fprintf(stderr,"%s",mdlpvp->time_str);
};


fprint_time_stmt: FPRINT_TIME '(' existing_file_stream ',' format_string ')'
{
  time_t the_time;

  mdlpvp->gp=$<sym>3;
  mdlpvp->filep=(struct file_stream *)mdlpvp->gp->value;
  mdlpvp->a_str=$<str>5;
  the_time=time(NULL);
  strftime(mdlpvp->time_str,128,mdlpvp->a_str,localtime(&the_time));
  if (fprintf(mdlpvp->filep->stream,"%s",mdlpvp->time_str)==EOF) {
    sprintf(mdlpvp->mdl_err_msg,"%s %s","Could not print to file:",mdlpvp->filep->name);
    mdlerror(mdlpvp->mdl_err_msg,mdlpvp);
    return(1);
  }
};





%%





/* Begin Bison Epilogue: */


struct mdlparse_vars *clunky_mpvp_for_errors;

void mdlerror(char *s,...)
{
/*  va_list ap; */
  FILE *log_file;
  struct mdlparse_vars *mpvp = clunky_mpvp_for_errors;

/*
  va_start(ap,s);
  mpvp=va_arg(ap,struct mdlparse_vars *);  
  va_end(ap);
*/

  log_file=stderr;
  if (mpvp->vol->err_file!=NULL) {
    log_file=mpvp->vol->err_file;
  }

  if (mpvp->vol->procnum == 0) {
	fprintf(log_file,"MCell: error on line: %d of file: %s\n  %s\n",\
	        mpvp->line_num[mpvp->include_stack_ptr],mpvp->vol->curr_file,s);
	
	fflush(log_file);
  }
	return;
}

void mdlerror_nested(char *s)
{
  FILE *log_file;
  struct mdlparse_vars *mpvp = clunky_mpvp_for_errors;

  log_file=stderr;
  if (mpvp->vol->err_file!=NULL) {
    log_file=mpvp->vol->err_file;
  }

  if (mpvp->vol->procnum == 0) {
    fprintf(log_file,"MCell error: %s\n",s);
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
  clunky_mpvp_for_errors = mpvp;

  mpvp->vol=vol;
  
  mpvp->cval=NULL;
  mpvp->cval_2=NULL;
  mpvp->strval=NULL;
  mpvp->ival=0;
  mpvp->rval=0;
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

  /* create default reflective surface reaction */
  if ((mpvp->path_mem=create_mem(sizeof(struct pathway),4096))==NULL)
  {
    sprintf(mpvp->mdl_err_msg,"Out of memory while creating reaction pathways\n");
    fprintf(vol->log_file,mpvp->mdl_err_msg);
    return(1);
  } 
  if ((mpvp->prod_mem=create_mem(sizeof(struct product),4096))==NULL) {
    sprintf(mpvp->mdl_err_msg,"Out of memory while creating reaction products\n");
    fprintf(vol->log_file,mpvp->mdl_err_msg);
      return(1);
  } 
  mpvp->sym_list_mem = create_mem(sizeof(struct sym_table_list),4096);
  if (mpvp->sym_list_mem==NULL)
  {
    fprintf(vol->err_file,"Out of memory while getting ready to store lists of molecules and reactions");
    return 1;
  }
#if 0
  mpvp->sym_name=concat_rx_name(vol->g_surf->sym->name,vol->g_mol->sym->name);
  if ((mpvp->gp=retrieve_sym(mpvp->sym_name,RX,vol->main_sym_table))
      !=NULL) {
  }
  else if ((mpvp->gp=store_sym(mpvp->sym_name,RX,vol->main_sym_table))
      ==NULL) {
    sprintf(mpvp->mdl_err_msg,"%s %s -%s-> ...",
      "Out of memory while creating surface reaction:",vol->g_mol->sym->name,vol->g_surf->sym->name);
    fprintf(vol->log_file,mpvp->mdl_err_msg);
    return(1);
  }
  if ((mpvp->pathp=(struct pathway *)mem_get(mpvp->path_mem))==NULL) {
    sprintf(mpvp->mdl_err_msg,"%s %s -%s-> ...",
      "Out of memory while creating surface reaction:",vol->g_mol->sym->name,vol->g_surf->sym->name);
    fprintf(vol->log_file,mpvp->mdl_err_msg);
    return(1);
  }
  mpvp->rxnp=(struct rxn *)mpvp->gp->value;
  mpvp->rxnp->n_reactants=2;
  mpvp->rxnp->n_pathways++;
  mpvp->pathp->pathname=NULL;
  mpvp->pathp->reactant1=vol->g_surf;
  mpvp->pathp->reactant2=vol->g_mol;
  mpvp->pathp->reactant3=NULL;
  mpvp->pathp->km=GIGANTIC;
  mpvp->pathp->kcat=0;
  mpvp->pathp->orientation1=0;
  mpvp->pathp->orientation2=1;
  mpvp->pathp->orientation3=0;
  mpvp->pathp->pcr=NULL;

  if ((mpvp->prodp=(struct product *)mem_get(mpvp->prod_mem))==NULL) {
    sprintf(mpvp->mdl_err_msg,"%s %s -%s-> ...",
      "Out of memory while creating surface reaction:",vol->g_mol->sym->name,vol->g_surf->sym->name);
    fprintf(vol->log_file,mpvp->mdl_err_msg);
    return(1);
  }
  mpvp->prodp->prod=mpvp->pathp->reactant2;
  mpvp->prodp->orientation=1;
  mpvp->prodp->next=NULL;
  mpvp->pathp->product_head=mpvp->prodp;
  mpvp->pathp->next=mpvp->rxnp->pathway_head;
  mpvp->rxnp->pathway_head=mpvp->pathp;
#endif

  mpvp->gp=NULL;
  mpvp->tp=NULL;
  mpvp->sym_name=NULL;

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


