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
"ACOS"			{return(ACOS);}
"ALL"			{return(ALL);}
"ALL_CROSSINGS"		{return(ALL_CROSSINGS);}
"ALL_ELEMENTS"		{return(ALL_ELEMENTS);}
"ALL_EVENTS"		{return(ALL_EVENTS);}
"ALL_HITS"		{return(ALL_HITS);}
"ASIN"			{return(ASIN);}
"ATAN"			{return(ATAN);}
"BACK"			{return(BACK);}
"BACK_CROSSINGS"	{return(BACK_CROSSINGS);}
"BACK_HITS"		{return(BACK_HITS);}
"BINDING_POLE"		{return(BINDING_POLE);}
"BOTH_POLES"		{return(BOTH_POLES);}
"BOTTOM"		{return(BOTTOM);}
"BOX"			{return(BOX);}
"CEIL"			{return(CEIL);}
"CHARGE"		{return(CHARGE);}
"CHECKPOINT_INFILE"	{return(CHECKPOINT_INFILE);}
"CHECKPOINT_OUTFILE"	{return(CHECKPOINT_OUTFILE);}
"CHECKPOINT_ITERATIONS"	{return(CHECKPOINT_ITERATIONS);}
"SURFACE_CLASS"		{return(SURFACE_CLASS);}
"COLOR"			{return(COLOR);}
"COLOR_EFFECTOR"	{return(COLOR_EFFECTOR);}
"COLOR_SIDE"		{return(COLOR_SIDE);}
"COMPARTMENT"		{return(COMPARTMENT);}
"CORNERS"		{return(CORNERS);}
"COS"			{return(COS);}
"CONCENTRATION"		{return(CONCENTRATION);}
"COUNT"			{return(COUNT);}
"CUBIC_RELEASE_SITE"	{return(CUBIC_RELEASE_SITE);}
"CUSTOM_RK"             {return(CUSTOM_RK);}
"CUMULATE_FOR_EACH_TIME_STEP" {return(CUMULATE_FOR_EACH_TIME_STEP);}
"DIFFUSION_CONSTANT" |
"DIFFUSION_CONSTANT_3D"	{return(DIFFUSION_CONSTANT_3D);}
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
"DELAY"			{return(DELAY);}
"DENSITY"		{return(DENSITY);}
"DX"			{return(DX);}
"EFFECTOR"		{return(EFFECTOR);}
"EFFECTOR_GRID_DENSITY"	{return(EFFECTOR_GRID_DENSITY);}
"EFFECTOR_POSITIONS"	{return(EFFECTOR_POSITIONS);}
"EFFECTOR_STATE"	{return(EFFECTOR_STATE);}
"EFFECTOR_STATES"	{return(EFFECTOR_STATES);}
"EITHER_POLE"		{return(EITHER_POLE);}
"ELEMENT"		{return(ELEMENT);}
"ELEMENT_CONNECTIONS"	{return(ELEMENT_CONNECTIONS);}
"ELEMENT_LIST"		{return(ELEMENT_LIST);}
"ELLIPTIC_RELEASE_SITE" {return(ELLIPTIC_RELEASE_SITE);}
"EXP"			{return(EXP);}
"EXPRESSION"		{return(EXPRESSION);}
"FALSE"			{return(FALSE);}
"FLOOR"			{return(FLOOR);}
"FOR_EACH_TIME_STEP"	{return(FOR_EACH_TIME_STEP);}
"FOR_EACH_MOLECULE" |
"FOR_EACH_LIGAND"	{return(FOR_EACH_MOLECULE);}
"FOR_EACH_EFFECTOR"	{return(FOR_EACH_EFFECTOR);}
"FRONT"			{return(FRONT);}
"FRONT_CROSSINGS"	{return(FRONT_CROSSINGS);}
"FRONT_HITS"		{return(FRONT_HITS);}
"FULLY_RANDOM"		{return(FULLY_RANDOM);}
"GAUSSIAN_RELEASE_NUMBER" {return(GAUSSIAN_RELEASE_NUMBER);}
"INCLUDE_FILE"		{return(INCLUDE_FILE);}
"INITIAL_EVENTS"	{return(INITIAL_EVENTS);}
"INPUT_FILE"		{return(INPUT_FILE);}
"INSTANTIATE"		{return(INSTANTIATE);}
"INTERACTION_RADIUS"    {return(INTERACTION_RADIUS);}
"INTERIM_EVENTS"	{return(INTERIM_EVENTS);}
"IRIT"			{return(IRIT);}
"ITERATIONS"		{return(ITERATIONS);}
"ITERATION_FRAME_DATA"	{return(ITERATION_FRAME_DATA);}
"ITERATION_LIST"	{return(ITERATION_LIST);}
"TIME_LIST"             {return(TIME_LIST);}
"LEFT"			{return(LEFT);}
"MEAN_DIAMETER"		{return(MEAN_DIAMETER);}
"MEAN_NUMBER"		{return(MEAN_NUMBER);}
"MODIFY_SURFACE_REGIONS" {return(MODIFY_SURFACE_REGIONS);}
"MOLECULE_DENSITY"	{return(MOLECULE_DENSITY);}
"MOLECULE_NUMBER"	{return(MOLECULE_NUMBER);}
"MOLECULE" |
"LIGAND"		{return(MOLECULE);}
"MOLECULE_POSITIONS" |
"LIGAND_POSITIONS"	{return(MOLECULE_POSITIONS);}
"MOLECULE_STATES" |
"LIGAND_STATES"		{return(MOLECULE_STATES);}
"LOCATION"		{return(LOCATION);}
"LOG10"			{return(LOG10);}
"LOG"			{return(LOG);}
"MAX"			{return(MAX_TOK);}
"MIN"			{return(MIN_TOK);}
"MCELL_GENERIC"		{return(MCELL_GENERIC);}
"MOD"			{return(MOD);}
"MODE"			{return(MODE);}
"MOLECULE_FILE_PREFIX"	{return(MOLECULE_FILE_PREFIX);}
"NAME"			{return(NAME);}
"NEGATIVE_POLE"		{return(NEGATIVE_POLE);}
"NORMAL"		{return(NORMAL);}
"NO"			{return(NO);}
"NONE"			{return(NONE);}
"NULL"			{return(NO_SPECIES);}
"NUMBER"		{return(NUMBER);}
"NUMBER_BOUND"		{return(NUMBER_BOUND);}
"NUMBER_OF_TRAINS"	{return(NUMBER_OF_TRAINS);}
"NUMBER_TO_RELEASE"	{return(NUMBER_TO_RELEASE);}
"OBJECT"		{return(OBJECT);}
"OBJECT_FILE_PREFIXES"	{return(OBJECT_FILE_PREFIXES);}
"ORIENTATION"		{return(ORIENTATION);}
"PARALLEL_PARTITION"	{return(PARALLEL_PARTITION);}
"PART"			{return(PART);}
"PARTITION_X"		{return(PARTITION_X);}
"PARTITION_Y"		{return(PARTITION_Y);}
"PARTITION_Z"		{return(PARTITION_Z);}
"PARTS"			{return(PARTS);}
"PI"			{return(PI_TOK);}
"POLE_ORIENTATION"	{return(POLE_ORIENTATION);}
"POLYGON_LIST"		{return(POLYGON_LIST);}
"POSITIVE_BACK"		{return(POSITIVE_BACK);}
"POSITIVE_FRONT"	{return(POSITIVE_FRONT);}
"POSITIVE_POLE"		{return(POSITIVE_POLE);}
"POVRAY"		{return(POVRAY);}
"RADIAL_DIRECTIONS"	{return(RADIAL_DIRECTIONS);}
"RADIAL_SUBDIVISIONS"	{return(RADIAL_SUBDIVISIONS);}
"RADIANCE"		{return(RADIANCE);}
"RAND_UNIFORM"		{return(RAND_UNIFORM);}
"RAND_GAUSSIAN"		{return(RAND_GAUSSIAN);}
"RAYSHADE"		{return(RAYSHADE);}
"REACTION_DATA_OUTPUT"	{return(REACTION_DATA_OUTPUT);}
"REACTION_GROUP"	{return(REACTION_GROUP);}
"RECTANGULAR_RELEASE_SITE" {return(RECTANGULAR_RELEASE_SITE); }
"REFERENCE_DIFFUSION_CONSTANT" {return(REFERENCE_DIFFUSION_CONSTANT);}
"REFERENCE_STATE"	{return(REFERENCE_STATE);}
"REFLECTIVE"		{return(REFLECTIVE);}
"RELEASE_INTERVAL"	{return(RELEASE_INTERVAL);}
"RELEASE_PATTERN"	{return(RELEASE_PATTERN);}
"RELEASE_PROBABILITY"	{return(RELEASE_PROBABILITY);}
"REMOVE_ELEMENT"	{return(REMOVE_ELEMENT);}
"RENDERMAN"		{return(RENDERMAN);}
"RIGHT"			{return(RIGHT);}
"ROTATE"		{return(ROTATE);}
"ROUND_OFF"		{return(ROUND_OFF);}
"SCALE"			{return(SCALE);}
"SEED"			{return(SEED);}
"SIN"			{return(SIN);}
"SITE_DIAMETER"		{return(SITE_DIAMETER);}
"SPECIFIED_MOLECULES" |
"SPECIFIED_LIGANDS"	{return(SPECIFIED_MOLECULES);}
"SPECIFIED_EFFECTORS"	{return(SPECIFIED_EFFECTORS);}
"SPHERICAL_RELEASE_SITE" {return(SPHERICAL_RELEASE_SITE);}
"SQRT"			{return(SQRT);}
"STANDARD_DEVIATION"    {return(STANDARD_DEVIATION);}
"STATE_VALUES"		{return(STATE_VALUES);}
"STATE"			{return(STATE);}
"STEP"			{return(STEP);}
"STRING_TO_NUM"		{return(STRING_TO_NUM);}
"SUM_OVER_ALL_MOLECULES" |
"SUM_OVER_ALL_LIGANDS"	{return(SUM_OVER_ALL_MOLECULES);}
"SUM_OVER_ALL_EFFECTORS" {return(SUM_OVER_ALL_EFFECTORS);}
"SUM_OVER_ALL_TIME_STEPS" {return(SUM_OVER_ALL_TIME_STEPS);}
"SURFACE_POSITIONS"	{return(SURFACE_POSITIONS);}
"SURFACE_STATES"	{return(SURFACE_STATES);}
"TAN"			{return(TAN);}
"TIME_STEP"		{return(TIME_STEP);}
"TIME_STEP_MAX"         {return(TIME_STEP_MAX);}
"TO"			{return(TO);}
"TOP"			{return(TOP);}
"TRAIN_DURATION"	{return(TRAIN_DURATION);}
"TRAIN_INTERVAL"	{return(TRAIN_INTERVAL);}
"TRANSFORM"		{return(TRANSFORM);}
"TRANSLATE"		{return(TRANSLATE);}
"TRANSPARENT"		{return(TRANSPARENT);}
"TRUE"			{return(TRUE);}
"UNLIMITED"		{return(UNLIMITED);}
"VERTEX_LIST"		{return(VERTEX_LIST);}
"VIZ_DATA_OUTPUT"	{return(VIZ_DATA_OUTPUT);}
"VOLUME_DEPENDENT_RELEASE_NUMBER" {return(VOLUME_DEPENDENT_RELEASE_NUMBER);}
"VOXEL_IMAGE_MODE"	{return(VOXEL_IMAGE_MODE);}
"VOXEL_VOLUME_MODE"	{return(VOXEL_VOLUME_MODE);}
"WINDOW"		{return(WINDOW);}
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
\"(([^\"]*)|((\\\")*))*\"  {mpvp->strval=my_strdup(mdltext); return(STR_VALUE);}
{R}		       {sscanf((const char *)mdltext,"%lf",&mpvp->rval); return(REAL);}
{I}			{sscanf((const char *)mdltext,"%d",&mpvp->ival); return(INTEGER);}
[A-Za-z][A-Za-z0-9\_\.]* {
			   if (mpvp->cval==NULL) {
			     mpvp->cval_2=NULL;
			     mpvp->cval=my_strdup(mdltext);
			   }
			   else {
			     mpvp->cval_2=mpvp->cval;
			     mpvp->cval=my_strdup(mdltext);
			   }
			   return(VAR);
			 }
"<"[-]+">"		{return(DBL_ARROW);}
"<"[-]+			{return(LFT_ARROW);}
[-]+">"			{return(RT_ARROW);}
[\'\,\(\)\/\-\+\=\^\[\]\{\}\|\>\*\#\~\@\:\&\;]	{return(mdltext[0]);}
.			;
\n			{mpvp->line_num[mpvp->include_stack_ptr]++;}
<<EOF>>			{return(EOF_TOK);}
