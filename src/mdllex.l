%{
  #include <stdio.h>
  #include <string.h>
  #include <stdlib.h>
  #include <unistd.h>
  #include "strfunc.h"
  #include "mdlparse.h"

  #define YY_DECL int mdllex( YYSTYPE *mdllval, struct mdlparse_vars *mpvp )

  #define YY_NO_UNPUT

  #ifdef __cplusplus
  #define INPUT_FUNC yyinput
  #else
  #define INPUT_FUNC input
  #endif
%}

%option prefix="mdl"
%option outfile="mdllex.flex.c"
%option noyywrap

%p 10000
%a 10000 
%o 10000
%e 2000
%n 1000

R ([0-9]+\.)|([0-9]*\.[0-9]+)|([0-9]+[eE][+-]?[0-9]+)|([0-9]+\.[0-9]*[eE][+-]?[0-9]+)|([0-9]*\.[0-9]+[eE][+-]?[0-9]+)
I ([0-9]+)

%%


"/*"            {
	                /* handle comments and nested comments */
                        int     c1, c2;
                        int     comline, incomment = 1;
 
                        comline = mpvp->line_num[mpvp->include_stack_ptr];
 
                        do {
                          while ((c1 = INPUT_FUNC()) != '*' 
	                          && c1 != '/' && c1 != EOF) {
                            if (c1 == '\n') {
                              mpvp->line_num[mpvp->include_stack_ptr]++;
	                    }
	                  }
 
                          c2 = INPUT_FUNC();
                          if (c2 == '\n') {
                            mpvp->line_num[mpvp->include_stack_ptr]++;
	                  }
 
	                  if (c1 == '*' && c2 == '/') {
	                    incomment--;
	                  }
   
                          if (c1 == '/' && c2 == '*') {
                            incomment++;
                            comline = mpvp->line_num[mpvp->include_stack_ptr];
                          }
 
                        } while (incomment && c2 != EOF);
 
                        if (c1 == EOF || c2 == EOF) {
                          sprintf(mpvp->mdl_err_msg,
	                    "Unterminated comment started on line %d",
	                     comline);
	                  mpvp->line_num[mpvp->include_stack_ptr]--;
                          mdlerror(mpvp->mdl_err_msg);
	                  return(1);
                        }
                }
[ \t]+			;
"ABS"			{return(ABS);}
"ABSORPTIVE"		{return(ABSORPTIVE);}
"ACCURATE_3D_REACTIONS" {return(ACCURATE_3D_REACTIONS);}
"ACOS"			{return(ACOS);}
"ALL_DATA"		{return(ALL_DATA);}
"ALL_CROSSINGS"		{return(ALL_CROSSINGS);}
"ALL_ELEMENTS"		{return(ALL_ELEMENTS);}
"ALL_ENCLOSED"		{return(ALL_ENCLOSED);}
"ALL_EVENTS"		{return(ALL_EVENTS);}
"ALL_HITS"		{return(ALL_HITS);}
"ALL_ITERATIONS"        {return(ALL_ITERATIONS);}
"ALL_MESHES"            {return (ALL_MESHES);}
"ALL_MOLECULES"         {return (ALL_MOLECULES);}
"ALL_NOTIFICATIONS"     {return(ALL_NOTIFICATIONS);}
"ALL_TIMES"		{return(ALL_TIMES);}
"ALL_WARNINGS"          {return(ALL_WARNINGS);}
"AREA_OCCUPIED"		{return(AREA_OCCUPIED);}
"ASCII"			{return(ASCII);}
"ASIN"			{return(ASIN);}
"ASPECT_RATIO"		{return(ASPECT_RATIO);}
"ATAN"			{return(ATAN);}
"BACK"			{return(BACK);}
"BACK_CROSSINGS"	{return(BACK_CROSSINGS);}
"BACK_HITS"		{return(BACK_HITS);}
"BINDING_POLE"		{return(BINDING_POLE);}
"BOTH_POLES"		{return(BOTH_POLES);}
"BOTTOM"		{return(BOTTOM);}
"BOX"			{return(BOX);}
"BOX_TRIANGULATION_REPORT" {return(BOX_TRIANGULATION_REPORT);}
"BRIEF"                 {return(BRIEF);}
"CEIL"			{return(CEIL);}
"CENTER_MOLECULES_ON_GRID" {return(CENTER_MOLECULES_ON_GRID);}
"CHARGE"		{return(CHARGE);}
"CHECKPOINT_INFILE"	{return(CHECKPOINT_INFILE);}
"CHECKPOINT_OUTFILE"	{return(CHECKPOINT_OUTFILE);}
"CHECKPOINT_ITERATIONS"	{return(CHECKPOINT_ITERATIONS);}
"CLAMP" |
"CLAMP_CONC" |
"CLAMP_CONCENTRATION"	{return(CLAMP_CONCENTRATION);}
"CLOSE_PARTITION_SPACING" {return(CLOSE_PARTITION_SPACING);}
"COLOR"			{return(COLOR);}
"COLOR_EFFECTOR"	{return(COLOR_EFFECTOR);}
"COLOR_SIDE"		{return(COLOR_SIDE);}
"COMPARTMENT"		{return(COMPARTMENT);}
"CORNERS"		{return(CORNERS);}
"COS"			{return(COS);}
"CONC" |
"CONCENTRATION"		{return(CONCENTRATION);}
"COUNT"			{return(COUNT);}
"CUBIC"			{return(CUBIC);}
"CUBIC_RELEASE_SITE"	{return(CUBIC_RELEASE_SITE);}
"CUSTOM_SPACE_STEP"	{return(CUSTOM_SPACE_STEP); }
"CUSTOM_RK"             {return(CUSTOM_RK);}
"CUSTOM_TIME_STEP"	{return(CUSTOM_TIME_STEP);}
"CUMULATE_FOR_EACH_TIME_STEP" {return(CUMULATE_FOR_EACH_TIME_STEP);}
"D_3D" |
"DIFFUSION_CONSTANT" |
"DIFFUSION_CONSTANT_3D"	{return(DIFFUSION_CONSTANT_3D);}
"D_2D" |
"DIFFUSION_CONSTANT_2D"	{return(DIFFUSION_CONSTANT_2D);}
"DEFINE"		{return(DEFINE);}
"DEFINE_EFFECTOR_SITE_POSITIONS"	{return(DEFINE_EFFECTOR_SITE_POSITIONS);}
"DEFINE_MOLECULE" 	{return(DEFINE_MOLECULE);}
"DEFINE_MOLECULES" 	{return(DEFINE_MOLECULES);}
"DEFINE_REACTIONS"	{return(DEFINE_REACTIONS);}
"DEFINE_RELEASE_PATTERN"	{return(DEFINE_RELEASE_PATTERN);}
"DEFINE_SURFACE_REGIONS"	{return(DEFINE_SURFACE_REGIONS);}
"DEFINE_SURFACE_CLASS"	{return(DEFINE_SURFACE_CLASS);}
"DEFINE_SURFACE_CLASSES"	{return(DEFINE_SURFACE_CLASSES);}
"DEGENERATE_POLYGONS"   {return(DEGENERATE_POLYGONS);}
"DELAY"			{return(DELAY);}
"DENSITY"		{return(DENSITY);}
"DIFFUSION_CONSTANT_REPORT" {return(DIFFUSION_CONSTANT_REPORT);}
"DX"			{return(DX);}
"DREAMM_V3"		{return(DREAMM_V3);}
"DREAMM_V3_GROUPED"	{return(DREAMM_V3_GROUPED);}
"EFFECTOR"		{return(EFFECTOR);}
"EFFECTOR_GRID_DENSITY" |
"SURFACE_GRID_DENSITY"	{return(EFFECTOR_GRID_DENSITY);}
"EFFECTOR_POSITIONS"	{return(EFFECTOR_POSITIONS);}
"EFFECTOR_STATE"	{return(EFFECTOR_STATE);}
"EFFECTOR_STATES"	{return(EFFECTOR_STATES);}
"EITHER_POLE"		{return(EITHER_POLE);}
"ELEMENT"		{return(ELEMENT);}
"ELEMENT_CONNECTIONS"	{return(ELEMENT_CONNECTIONS);}
"ELEMENT_LIST"		{return(INCLUDE_ELEMENTS);}
"ELLIPTIC"		{return(ELLIPTIC);}
"ELLIPTIC_RELEASE_SITE" {return(ELLIPTIC_RELEASE_SITE);}
"ERROR"                 {return(ERROR);}
"ESTIMATE_CONC" |
"ESTIMATE_CONCENTRATION" {return(ESTIMATE_CONCENTRATION);}
"EXCLUDE_ELEMENTS"	{return(EXCLUDE_ELEMENTS);}
"EXCLUDE_PATCH"		{return(EXCLUDE_PATCH);}
"EXCLUDE_REGION"	{return(EXCLUDE_REGION);}
"EXP"			{return(EXP);}
"EXPRESSION"		{return(EXPRESSION);}
"FALSE"			{return(FALSE);}
"FILENAME"	        {return(FILENAME);}
"FILE_OUTPUT_REPORT"    {return(FILE_OUTPUT_REPORT);}
"FINAL_SUMMARY"         {return(FINAL_SUMMARY);}
"FLOOR"			{return(FLOOR);}
"FOR_EACH_TIME_STEP"	{return(FOR_EACH_TIME_STEP);}
"FOR_EACH_MOLECULE" |
"FOR_EACH_LIGAND"	{return(FOR_EACH_MOLECULE);}
"FOR_EACH_EFFECTOR"	{return(FOR_EACH_EFFECTOR);}
"FRAME_DATA"		{return(FRAME_DATA);}
"FRONT"			{return(FRONT);}
"FRONT_CROSSINGS"	{return(FRONT_CROSSINGS);}
"FRONT_HITS"		{return(FRONT_HITS);}
"FULLY_RANDOM"		{return(FULLY_RANDOM);}
"GAUSSIAN_RELEASE_NUMBER" {return(GAUSSIAN_RELEASE_NUMBER);}
"GEOMETRY" 		{return(GEOMETRY);}
"HEADER"		{return(HEADER);}
"HIGH_PROBABILITY_THRESHOLD" {return(HIGH_PROBABILITY_THRESHOLD);}
"HIGH_REACTION_PROBABILITY" {return(HIGH_REACTION_PROBABILITY);}
"IGNORE" |
"IGNORED"               {return(IGNORED);}
"INCLUDE_ELEMENTS"	{return(INCLUDE_ELEMENTS);}
"INCLUDE_FILE"		{return(INCLUDE_FILE);}
"INCLUDE_PATCH"		{return(INCLUDE_PATCH);}
"INCLUDE_REGION"	{return(INCLUDE_REGION);}
"INITIAL_EVENTS"	{return(INITIAL_EVENTS);}
"INPUT_FILE"		{return(INPUT_FILE);}
"INSTANTIATE"		{return(INSTANTIATE);}
"INTERACTION_RADIUS"    {return(INTERACTION_RADIUS);}
"INTERIM_EVENTS"	{return(INTERIM_EVENTS);}
"IRIT"			{return(IRIT);}
"ITERATIONS"		{return(ITERATIONS);}
"ITERATION_FRAME_DATA"	{return(ITERATION_FRAME_DATA);}
"ITERATION_LIST"	{return(ITERATION_LIST);}
"ITERATION_NUMBERS"     {return(ITERATION_NUMBERS);}
"ITERATION_REPORT"      {return(ITERATION_REPORT);}
"LEFT"			{return(LEFT);}
"LIFETIME_TOO_SHORT"    {return(LIFETIME_TOO_SHORT);}
"LIFETIME_THRESHOLD"    {return(LIFETIME_THRESHOLD);}
"LIST"			{return(LIST);}
"LOCATION"		{return(LOCATION);}
"LOG10"			{return(LOG10);}
"LOG"			{return(LOG);}
"MAX"			{return(MAX_TOK);}
"MCELL_GENERIC"		{return(MCELL_GENERIC);}
"MEAN_DIAMETER"		{return(MEAN_DIAMETER);}
"MEAN_NUMBER"		{return(MEAN_NUMBER);}
"MESHES"		{return(MESHES);}
"MICRO_REV_SURF_RX"     {return(MICRO_REV_SURF_RX);}
"MICRO_REV_VOL_RX"      {return(MICRO_REV_VOL_RX);}
"MIN"			{return(MIN_TOK);}
"MISSED_REACTIONS"      {return(MISSED_REACTIONS);}
"MISSED_REACTION_THRESHOLD" {return(MISSED_REACTION_THRESHOLD);}
"MISSING_SURFACE_ORIENTATION" {return(MISSING_SURFACE_ORIENTATION);}
"MOD"			{return(MOD);}
"MODE"			{return(MODE);}
"MODIFY_SURFACE_REGIONS" {return(MODIFY_SURFACE_REGIONS);}
"MOLECULE_DENSITY"	{return(MOLECULE_DENSITY);}
"MOLECULE_NUMBER"	{return(MOLECULE_NUMBER);}
"MOLECULE" |
"LIGAND"		{return(MOLECULE);}
"MOLECULES"             {return(MOLECULES);}
"MOLECULE_POSITIONS" |
"LIGAND_POSITIONS"	{return(MOLECULE_POSITIONS);}
"MOLECULE_STATES" |
"LIGAND_STATES"		{return(MOLECULE_STATES);}
"MOLECULE_POSITIONS_STATES"	{return(MOLECULE_POSITIONS_STATES);}
"MOLECULE_FILE_PREFIX"	{return(MOLECULE_FILE_PREFIX);}
"NAME"			{return(NAME);}
"NAME_LIST"		{return(NAME_LIST);}
"NEGATIVE_DIFFUSION_CONSTANT" {return(NEGATIVE_DIFFUSION_CONSTANT);}
"NEGATIVE_REACTION_RATE" {return(NEGATIVE_REACTION_RATE);}
"NEGATIVE_POLE"		{return(NEGATIVE_POLE);}
"NORMAL"		{return(NORMAL);}
"NO"			{return(NO);}
"NONE"			{return(NONE);}
"NOTIFICATIONS"         {return(NOTIFICATIONS);}
"NULL"			{return(NO_SPECIES);}
"NUMBER"		{return(NUMBER);}
"NUMBER_BOUND"		{return(NUMBER_BOUND);}
"NUMBER_OF_TRAINS"	{return(NUMBER_OF_TRAINS);}
"NUMBER_TO_RELEASE"	{return(NUMBER_TO_RELEASE);}
"OBJECT"		{return(OBJECT);}
"OBJECT_FILE_PREFIXES"	{return(OBJECT_FILE_PREFIXES);}
"OFF"                   {return(OFF);}
"ON"                    {return(ON);}
"ORIENTATION"		{return(ORIENTATION);}
"ORIENTATIONS"		{return(ORIENTATIONS);}
"OUTPUT_BUFFER_SIZE"    {return(OUTPUT_BUFFER_SIZE);}
"OVERWRITTEN_OUTPUT_FILE" {return(OVERWRITTEN_OUTPUT_FILE);}
"PARALLEL_PARTITION"	{return(PARALLEL_PARTITION);}
"PART"			{return(PART);}
"PARTITION_LOCATION_REPORT" {return(PARTITION_LOCATION_REPORT);}
"PARTITION_X"		{return(PARTITION_X);}
"PARTITION_Y"		{return(PARTITION_Y);}
"PARTITION_Z"		{return(PARTITION_Z);}
"PARTS"			{return(PARTS);}
"PI"			{return(PI_TOK);}
"POLE_ORIENTATION"	{return(POLE_ORIENTATION);}
"POLYGON_LIST"		{return(POLYGON_LIST);}
"POSITIONS"             {return(POSITIONS);}
"POSITIVE_BACK"		{return(POSITIVE_BACK);}
"POSITIVE_FRONT"	{return(POSITIVE_FRONT);}
"POSITIVE_POLE"		{return(POSITIVE_POLE);}
"PROBABILITY_REPORT"    {return(PROBABILITY_REPORT);}
"PROBABILITY_REPORT_THRESHOLD" {return(PROBABILITY_REPORT_THRESHOLD);}
"POVRAY"		{return(POVRAY);}
"PROGRESS_REPORT"       {return(PROGRESS_REPORT);}
"RADIAL_DIRECTIONS"	{return(RADIAL_DIRECTIONS);}
"RADIAL_SUBDIVISIONS"	{return(RADIAL_SUBDIVISIONS);}
"RADIANCE"		{return(RADIANCE);}
"RAND_UNIFORM"		{return(RAND_UNIFORM);}
"RAND_GAUSSIAN"		{return(RAND_GAUSSIAN);}
"RAYSHADE"		{return(RAYSHADE);}
"REACTION_DATA_OUTPUT"	{return(REACTION_DATA_OUTPUT);}
"REACTION_GROUP"	{return(REACTION_GROUP);}
"RECTANGULAR"		{return(RECTANGULAR_TOKEN);}
"RECTANGULAR_RELEASE_SITE" {return(RECTANGULAR_RELEASE_SITE); }
"REFERENCE_DIFFUSION_CONSTANT" {return(REFERENCE_DIFFUSION_CONSTANT);}
"REFERENCE_STATE"	{return(REFERENCE_STATE);}
"REFLECTIVE"		{return(REFLECTIVE);}
"REGION_DATA"	        {return(REGION_DATA);}
"RELEASE_EVENT_REPORT"  {return(RELEASE_EVENT_REPORT);}
"RELEASE_INTERVAL"	{return(RELEASE_INTERVAL);}
"RELEASE_PATTERN"	{return(RELEASE_PATTERN);}
"RELEASE_PROBABILITY"	{return(RELEASE_PROBABILITY);}
"RELEASE_SITE"		{return(RELEASE_SITE);}
"REMOVE_ELEMENTS"	{return(REMOVE_ELEMENTS);}
"RENDERMAN"		{return(RENDERMAN);}
"RIGHT"			{return(RIGHT);}
"ROTATE"		{return(ROTATE);}
"ROUND_OFF"		{return(ROUND_OFF);}
"SCALE"			{return(SCALE);}
"SEED"			{return(SEED);}
"SHAPE"			{return(SHAPE);}
"SHOW_EXACT_TIME"       {return(SHOW_EXACT_TIME);}
"SIN"			{return(SIN);}
"SITE_DIAMETER"		{return(SITE_DIAMETER);}
"SITE_RADIUS"		{return(SITE_RADIUS);}
"SPACE_STEP"		{return(SPACE_STEP);}
"SPECIFIED_MOLECULES" |
"SPECIFIED_LIGANDS"	{return(SPECIFIED_MOLECULES);}
"SPECIFIED_EFFECTORS"	{return(SPECIFIED_EFFECTORS);}
"SPHERICAL"		{return(SPHERICAL);}
"SPHERICAL_RELEASE_SITE" {return(SPHERICAL_RELEASE_SITE);}
"SPHERICAL_SHELL"	{return(SPHERICAL_SHELL);}
"SPHERICAL_SHELL_SITE"	{return(SPHERICAL_SHELL_SITE); }
"SQRT"			{return(SQRT);}
"STANDARD_DEVIATION"    {return(STANDARD_DEVIATION);}
"STATE_VALUES"		{return(STATE_VALUES);}
"STATE"			{return(STATE);}
"STEP"			{return(STEP);}
"STRING_TO_NUM"		{return(STRING_TO_NUM);}
"SUM"                   {return(SUMMATION_OPERATOR);}
"SUM_OVER_ALL_MOLECULES" |
"SUM_OVER_ALL_LIGANDS"	{return(SUM_OVER_ALL_MOLECULES);}
"SUM_OVER_ALL_EFFECTORS" {return(SUM_OVER_ALL_EFFECTORS);}
"SUM_OVER_ALL_TIME_STEPS" {return(SUM_OVER_ALL_TIME_STEPS);}
"SURFACE_CLASS"		{return(SURFACE_CLASS);}
"SURFACE_MOLECULE_AREA" {return(SURFACE_MOLECULE_AREA);}
"SURFACE_POSITIONS"	{return(SURFACE_POSITIONS);}
"SURFACE_STATES"	{return(SURFACE_STATES);}
"TAN"			{return(TAN);}
"TARGET_ONLY"		{return(TARGET_ONLY);}
"TET_ELEMENT_CONNECTIONS" {return(TET_ELEMENT_CONNECTIONS);}
"TIME_LIST"             {return(TIME_LIST);}
"TIME_POINTS"           {return(TIME_POINTS);}
"TIME_STEP"		{return(TIME_STEP);}
"TIME_STEP_MAX"         {return(TIME_STEP_MAX);}
"TO"			{return(TO);}
"TOP"			{return(TOP);}
"TRAIN_DURATION"	{return(TRAIN_DURATION);}
"TRAIN_INTERVAL"	{return(TRAIN_INTERVAL);}
"TRANSFORM"		{return(TRANSFORM);}
"TRANSLATE"		{return(TRANSLATE);}
"TRANSPARENT"		{return(TRANSPARENT);}
"TRIGGER"		{return(TRIGGER);}
"TRUE"			{return(TRUE);}
"UNLIMITED"		{return(UNLIMITED);}
"VACANCY_SEARCH_DISTANCE" {return(VACANCY_SEARCH_DISTANCE);}
"VARYING_PROBABILITY_REPORT" {return(VARYING_PROBABILITY_REPORT);}
"USELESS_VOLUME_ORIENTATION" {return(USELESS_VOLUME_ORIENTATION); }
"VERTEX_LIST"		{return(VERTEX_LIST);}
"VIZ_DATA_OUTPUT"	{return(VIZ_DATA_OUTPUT);}
"VIZ_OUTPUT"	        {return(VIZ_OUTPUT);}
"VIZ_VALUE"	        {return(VIZ_VALUE);}
"VOLUME_DEPENDENT_RELEASE_NUMBER" {return(VOLUME_DEPENDENT_RELEASE_NUMBER);}
"VOXEL_IMAGE_MODE"	{return(VOXEL_IMAGE_MODE);}
"VOXEL_LIST"            {return(VOXEL_LIST);}
"VOXEL_VOLUME_MODE"	{return(VOXEL_VOLUME_MODE);}
"WARNING"               {return(WARNING);}
"WARNINGS"              {return(WARNINGS);}
"WORLD"			{return(WORLD);}
"YES"			{return(YES);}
"X"			{return(X_TOK);}
"XY"			{return(XY_TOK);}
"XZ"			{return(XZ_TOK);}
"XYZ"			{return(XYZ_TOK);}
"Y"			{return(Y_TOK);}
"YZ"			{return(YZ_TOK);}
"Z"			{return(Z_TOK);}
"printf"		{return(PRINTF);}
"fprintf"		{return(FPRINTF);}
"sprintf"		{return(SPRINTF);}
"print_time"		{return(PRINT_TIME);}
"fprint_time"		{return(FPRINT_TIME);}
"=fopen"		{return(FOPEN);}
"= fopen"		{return(FOPEN);}
"fclose"		{return(FCLOSE);}
{R}		       {sscanf((const char *)mdltext,"%lf",&mpvp->rval); return(REAL);}
{I}			{sscanf((const char *)mdltext,"%d",&mpvp->ival); return(INTEGER);}
[A-Za-z][A-Za-z0-9\_\.]* {
			   if (mpvp->cval==NULL) {
			     mpvp->cval_2=NULL;
			     mpvp->cval=strdup(mdltext);
			   }
			   else {
			     mpvp->cval_2=mpvp->cval;
			     mpvp->cval=strdup(mdltext);
			   }
                           if(mpvp->cval == NULL){
                                fprintf(stderr,"File '%s', Line %ld: Out of memory while parsing variable name: %s\n", __FILE__, (long)__LINE__, mdltext);
                                exit(1);
                           }

			   return(VAR);
			 }
\"[A-Za-z][A-Za-z0-9\_\.]*;\" | 
\"[A-Za-z][A-Za-z0-9\_\.]*(\,)+\" | 
\"[A-Za-z][A-Za-z0-9\_\.]*(\')+\" {
			   if (mpvp->cval==NULL) {
			     mpvp->cval_2=NULL;
			     mpvp->cval=strdup(mdltext);
			   }
			   else {
			     mpvp->cval_2=mpvp->cval;
			     mpvp->cval=strdup(mdltext);
			   }
                           if(mpvp->cval == NULL){
                                fprintf(stderr,"File '%s', Line %ld: Out of memory while parsing variable name: %s\n", __FILE__, (long)__LINE__, mdltext);
                                exit(1);
                           }

			   return(VAR_ORIENT_IN_QUOTES);
			 }
[A-Za-z][A-Za-z0-9\_\.]*\{-?[0-9]+\} {
			   if (mpvp->cval==NULL) {
			     mpvp->cval_2=NULL;
			     mpvp->cval=strdup(mdltext);
			   }
			   else {
			     mpvp->cval_2=mpvp->cval;
			     mpvp->cval=strdup(mdltext);
			   }
                           if(mpvp->cval == NULL){
                                fprintf(stderr,"File '%s', Line %ld: Out of memory while parsing variable name: %s\n", __FILE__, (long)__LINE__, mdltext);
                                exit(1);
                           }

			   return(VAR_ORIENT_IN_BRACES);
                         }
\"[A-Za-z\[\]\-\*\?][A-Za-z0-9\_\.\[\]\-\*\?]*\" {
			   if (mpvp->cval==NULL) {
			     mpvp->cval_2=NULL;
			     mpvp->cval=strdup(mdltext);
			   }
			   else {
			     mpvp->cval_2=mpvp->cval;
			     mpvp->cval=strdup(mdltext);
			   }
                           if(mpvp->cval == NULL){
                                fprintf(stderr,"File '%s', Line %ld: Out of memory while parsing wildcard variable name: %s\n", __FILE__, (long)__LINE__, mdltext);
                                exit(1);
                           }

			   return(WILDCARD_VAR);
			 }
\"(([^\"]*)|((\\\")*))*\"  {
                             mpvp->strval=strdup(mdltext); 
			     if(mpvp->strval == NULL){
                                fprintf(stderr,"File '%s', Line %ld: Out of memory while parsing string: %s\n", __FILE__, (long)__LINE__, mdltext);
                                exit(1);
                             }
                             return(STR_VALUE);
                           }

[\'\,\(\)\/\-\+\=\^\[\]\{\}\|\<\>\*\#\~\@\:\&\;] {return(mdltext[0]);}
.			;
\n			{mpvp->line_num[mpvp->include_stack_ptr]++;}
<<EOF>>			{return(EOF_TOK);}
