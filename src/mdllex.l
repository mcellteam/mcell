%{
  #include "config.h"
  #include <stdio.h>
  #include <stdlib.h>
  #include <string.h>
  #include <unistd.h>
  #include "create_species.h"
  #include "strfunc.h"
  #include "mdlparse_aux.h"
  #include "mdlparse.h"
  #include "mem_util.h"
  #include "logging.h"

  #define YY_DECL int mdllex( YYSTYPE *yylval, struct mdlparse_vars *parse_state, yyscan_t yyscanner )
  #define YY_NO_UNPUT

  #ifdef __cplusplus
  #define INPUT_FUNC yyinput
  #else
  #define INPUT_FUNC input
  #endif

void *mdlalloc(yy_size_t size, yyscan_t yyscanner)
{
  UNUSED(yyscanner);
  return CHECKED_MALLOC(size, "memory buffer for parser");
}

void *mdlrealloc(void *ptr, yy_size_t size, yyscan_t yyscanner)
{
  UNUSED(yyscanner);

  ptr = (void *) realloc((char *) ptr, size);
  if (ptr == NULL)
    mcell_allocfailed("Failed to resize memory buffer for parser.");
  return ptr;
}

void mdlfree (void * ptr, yyscan_t yyscanner)
{
  UNUSED(yyscanner);

  free((char *) ptr);
}

#define YY_FATAL_ERROR(msg) mcell_error("%s", msg)

%}

/* Define state for parsing comments */
%x IN_COMMENT

/* Reentrant lexer allows safer handling of include files from the parser */
%option reentrant

/* Set up function name prefixes and output file name */
%option prefix="mdl"
%option outfile="lex.yy.c"

/* Enable condition state stack */
%option stack

/* Disable yywrap feature */
%option noyywrap

/* Eliminate unused symbols to eliminate compiler warnings */
%option noinput
%option nounput
%option noyy_top_state
%option noyyalloc
%option noyyrealloc
%option noyyfree

R ([0-9]+\.)|([0-9]*\.[0-9]+)|([0-9]+[eE][+-]?[0-9]+)|([0-9]+\.[0-9]*[eE][+-]?[0-9]+)|([0-9]*\.[0-9]+[eE][+-]?[0-9]+)
I ([0-9]+)
IDENT [A-Za-z][A-Za-z0-9\_\.]*

%%

"//"[^\n]*"\n"          { parse_state->line_num[parse_state->include_stack_ptr - 1] ++; }

<INITIAL>"/*"   {
                  parse_state->comment_started = parse_state->line_num[parse_state->include_stack_ptr - 1];
                  yy_push_state(IN_COMMENT, yyscanner);
                }
<IN_COMMENT>{
[^*/\n]*        { }
"/"[^*\n]       { }
"*"+[^*/\n]*    { }
"/*"            { yy_push_state(YY_START, yyscanner); }
"*"+"/"         { yy_pop_state(yyscanner); }
<<EOF>>         {
                  parse_state->line_num[parse_state->include_stack_ptr - 1]--;
                  mdlerror_fmt(parse_state,
                               "Unterminated comment started on line %d",
                               parse_state->comment_started);
                  return 1;
                }
}

"FORMAT"                { return FORMAT; }
[ \t]+			;
"=="                    { return EQUAL; }
"!="                    { return NOT_EQUAL; }
"ABS"			{return(ABS);}
"ABSORPTIVE"		{return(ABSORPTIVE);}
"ACCURATE_3D_REACTIONS" {return(ACCURATE_3D_REACTIONS);}
"ACOS"			{return(ACOS);}
"ALL_DATA"		{return(ALL_DATA);}
"ALL_CROSSINGS"		{return(ALL_CROSSINGS);}
"ALL_ELEMENTS"		{return(ALL_ELEMENTS);}
"ALL_ENCLOSED"		{return(ALL_ENCLOSED);}
"ALL_HITS"		{return(ALL_HITS);}
"ALL_ITERATIONS"        {return(ALL_ITERATIONS);}
"ALL_MESHES"            {return (ALL_MESHES);}
"ALL_MOLECULES"         {return (ALL_MOLECULES);}
"ALL_NOTIFICATIONS"     {return(ALL_NOTIFICATIONS);}
"ALL_TIMES"		{return(ALL_TIMES);}
"ALL_WARNINGS"          {return(ALL_WARNINGS);}
"ASCII"			{return(ASCII);}
"ASIN"			{return(ASIN);}
"ASPECT_RATIO"		{return(ASPECT_RATIO);}
"ATAN"			{return(ATAN);}
"BACK"			{return(BACK);}
"BACK_CROSSINGS"	{return(BACK_CROSSINGS);}
"BACK_HITS"		{return(BACK_HITS);}
"BINARY"		{return(BINARY);}
"BOTTOM"		{return(BOTTOM);}
"BOX"			{return(BOX);}
"BOX_TRIANGULATION_REPORT" {return(BOX_TRIANGULATION_REPORT);}
"BRIEF"                 {return(BRIEF);}
"CEIL"			{return(CEIL);}
"CELLBLENDER"		{return(CELLBLENDER);}
"CENTER_MOLECULES_ON_GRID" {return(CENTER_MOLECULES_ON_GRID);}
"CHECKPOINT_INFILE"	{return(CHECKPOINT_INFILE);}
"CHECKPOINT_OUTFILE"	{return(CHECKPOINT_OUTFILE);}
"CHECKPOINT_ITERATIONS"	{return(CHECKPOINT_ITERATIONS);}
"CHECKPOINT_REALTIME"	{return(CHECKPOINT_REALTIME);}
"CHECKPOINT_REPORT"	{return(CHECKPOINT_REPORT);}
"CLAMP" |
"CLAMP_CONC" |
"CLAMP_CONCENTRATION"	{return(CLAMP_CONCENTRATION);}
"CLOSE_PARTITION_SPACING"     {return CLOSE_PARTITION_SPACING;}
"COMPLEX_PLACEMENT_ATTEMPTS"  { return COMPLEX_PLACEMENT_ATTEMPTS; }
"COMPLEX_PLACEMENT_FAILURE"   { return COMPLEX_PLACEMENT_FAILURE; }
"COMPLEX_PLACEMENT_FAILURE_THRESHOLD"   { return COMPLEX_PLACEMENT_FAILURE_THRESHOLD; }
"COMPLEX_RATE"		      { return COMPLEX_RATE; }
"CORNERS"		{return(CORNERS);}
"COS"			{return(COS);}
"CONC" |
"CONCENTRATION"		{return(CONCENTRATION);}
"COUNT"			{return(COUNT);}
"CUBIC"			{return(CUBIC);}
"CUBIC_RELEASE_SITE"	{return(CUBIC_RELEASE_SITE);}
"CUSTOM_SPACE_STEP"	{return(CUSTOM_SPACE_STEP); }
"CUSTOM_TIME_STEP"	{return(CUSTOM_TIME_STEP);}
"D_3D" |
"DIFFUSION_CONSTANT" |
"DIFFUSION_CONSTANT_3D"	{return(DIFFUSION_CONSTANT_3D);}
"D_2D" |
"DIFFUSION_CONSTANT_2D"	{return(DIFFUSION_CONSTANT_2D);}
"DEFAULT"		{ return DEFAULT; }
"DEFINE_COMPLEX_MOLECULE" { return DEFINE_COMPLEX_MOLECULE; }
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
"EFFECTOR_GRID_DENSITY" |
"SURFACE_GRID_DENSITY"	{return(EFFECTOR_GRID_DENSITY);}
"EFFECTOR_POSITIONS"	{return(EFFECTOR_POSITIONS);}
"EFFECTOR_STATES"	{return(EFFECTOR_STATES);}
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
"EXIT"                  {return(EXIT);}
"EXP"			{return(EXP);}
"EXPRESSION"		{return(EXPRESSION);}
"FALSE"			{return(FALSE);}
"FILENAME"	        {return(FILENAME);}
"FILENAME_PREFIX"       {return FILENAME_PREFIX;}
"FILE_OUTPUT_REPORT"    {return(FILE_OUTPUT_REPORT);}
"FINAL_SUMMARY"         {return(FINAL_SUMMARY);}
"FLOOR"			{return(FLOOR);}
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
"INPUT_FILE"		{return(INPUT_FILE);}
"INSTANTIATE"		{return(INSTANTIATE);}
"INTERACTION_RADIUS"    {return(INTERACTION_RADIUS);}
"INVALID_OUTPUT_STEP_TIME"   { return INVALID_OUTPUT_STEP_TIME; }
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
"MAXIMUM_STEP_LENGTH"	{return(MAXIMUM_STEP_LENGTH);}
"MEAN_DIAMETER"		{return(MEAN_DIAMETER);}
"MEAN_NUMBER"		{return(MEAN_NUMBER);}
"MEMORY_PARTITION_X"    { return MEMORY_PARTITION_X; }
"MEMORY_PARTITION_Y"    { return MEMORY_PARTITION_Y; }
"MEMORY_PARTITION_Z"    { return MEMORY_PARTITION_Z; }
"MEMORY_PARTITION_POOL" { return MEMORY_PARTITION_POOL; }
"MESHES"		{return(MESHES);}
"MICROSCOPIC_REVERSIBILITY" {return(MICROSCOPIC_REVERSIBILITY);}
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
"MOLECULE_COLLISION_REPORT" {return(MOLECULE_COLLISION_REPORT);} 
"MOLECULE_POSITIONS" |
"LIGAND_POSITIONS"	{return(MOLECULE_POSITIONS);}
"MOLECULE_STATES" |
"LIGAND_STATES"		{return(MOLECULE_STATES);}
"MOLECULE_FILE_PREFIX"	{return(MOLECULE_FILE_PREFIX);}
"MOLECULE_PLACEMENT_FAILURE"    { return MOLECULE_PLACEMENT_FAILURE; }
"NAME_LIST"		{return(NAME_LIST);}
"NEGATIVE_DIFFUSION_CONSTANT" {return(NEGATIVE_DIFFUSION_CONSTANT);}
"NEGATIVE_REACTION_RATE" {return(NEGATIVE_REACTION_RATE);}
"NO"			{return(NO);}
"NOEXIT"                {return(NOEXIT);}
"NONE"			{return(NONE);}
"NOTIFICATIONS"         {return(NOTIFICATIONS);}
"NULL"			{return(NO_SPECIES);}
"NUMBER_OF_SUBUNITS" |
"NUMBER_OF_SLOTS"      {return(NUMBER_OF_SUBUNITS);}
"NUMBER_OF_TRAINS"	{return(NUMBER_OF_TRAINS);}
"NUMBER_TO_RELEASE"	{return(NUMBER_TO_RELEASE);}
"OBJECT"		{return(OBJECT);}
"OBJECT_FILE_PREFIXES"	{return(OBJECT_FILE_PREFIXES);}
"OFF"                   {return(OFF);}
"ON"                    {return(ON);}
"ORIENTATIONS"		{return(ORIENTATIONS);}
"OUTPUT_BUFFER_SIZE"    {return(OUTPUT_BUFFER_SIZE);}
"OVERWRITTEN_OUTPUT_FILE" {return(OVERWRITTEN_OUTPUT_FILE);}
"PARTITION_LOCATION_REPORT" {return(PARTITION_LOCATION_REPORT);}
"PARTITION_X"		{return(PARTITION_X);}
"PARTITION_Y"		{return(PARTITION_Y);}
"PARTITION_Z"		{return(PARTITION_Z);}
"PI"			{return(PI_TOK);}
"POLYGON_LIST"		{return(POLYGON_LIST);}
"POSITIONS"             {return(POSITIONS);}
"PROBABILITY_REPORT"    {return(PROBABILITY_REPORT);}
"PROBABILITY_REPORT_THRESHOLD" {return(PROBABILITY_REPORT_THRESHOLD);}
"PROGRESS_REPORT"       {return(PROGRESS_REPORT);}
"RADIAL_DIRECTIONS"	{return(RADIAL_DIRECTIONS);}
"RADIAL_SUBDIVISIONS"	{return(RADIAL_SUBDIVISIONS);}
"RAND_UNIFORM"		{return(RAND_UNIFORM);}
"RAND_GAUSSIAN"		{return(RAND_GAUSSIAN);}
"RATE_RULES"		{ return RATE_RULES; }
"REACTION_DATA_OUTPUT"	{return(REACTION_DATA_OUTPUT);}
"REACTION_OUTPUT_REPORT" { return REACTION_OUTPUT_REPORT; }
"REACTION_GROUP"	{return(REACTION_GROUP);}
"RECTANGULAR"		{return(RECTANGULAR_TOKEN);}
"RECTANGULAR_RELEASE_SITE" {return(RECTANGULAR_RELEASE_SITE); }
"REFERENCE_DIFFUSION_CONSTANT" {return(REFERENCE_DIFFUSION_CONSTANT);}
"REFLECTIVE"		{return(REFLECTIVE);}
"REGION_DATA"	        {return(REGION_DATA);}
"RELEASE_EVENT_REPORT"  {return(RELEASE_EVENT_REPORT);}
"RELEASE_INTERVAL"	{return(RELEASE_INTERVAL);}
"RELEASE_PATTERN"	{return(RELEASE_PATTERN);}
"RELEASE_PROBABILITY"	{return(RELEASE_PROBABILITY);}
"RELEASE_SITE"		{return(RELEASE_SITE);}
"REMOVE_ELEMENTS"	{return(REMOVE_ELEMENTS);}
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
"SPHERICAL"		{return(SPHERICAL);}
"SPHERICAL_RELEASE_SITE" {return(SPHERICAL_RELEASE_SITE);}
"SPHERICAL_SHELL"	{return(SPHERICAL_SHELL);}
"SPHERICAL_SHELL_SITE"	{return(SPHERICAL_SHELL_SITE); }
"SQRT"			{return(SQRT);}
"STANDARD_DEVIATION"    {return(STANDARD_DEVIATION);}
"STATE_VALUES"		{return(STATE_VALUES);}
"STEP"			{return(STEP);}
"STRING_TO_NUM"		{return(STRING_TO_NUM);}
"SUBUNIT" |
"SLOT"      { return SUBUNIT; }
"SUBUNIT_RELATIONSHIPS" |
"SLOT_RELATIONSHIPS"       { return SUBUNIT_RELATIONSHIPS; }
"SUM"                   {return(SUMMATION_OPERATOR);}
"SURFACE_CLASS"		{return(SURFACE_CLASS);}
"SURFACE_ONLY"          {return(SURFACE_ONLY);}
"SURFACE_POSITIONS"	{return(SURFACE_POSITIONS);}
"SURFACE_STATES"	{return(SURFACE_STATES);}
"TAN"			{return(TAN);}
"TARGET_ONLY"		{return(TARGET_ONLY);}
"TET_ELEMENT_CONNECTIONS" {return(TET_ELEMENT_CONNECTIONS);}
"THROUGHPUT_REPORT"     {return THROUGHPUT_REPORT;}
"TIME_LIST"             {return(TIME_LIST);}
"TIME_POINTS"           {return(TIME_POINTS);}
"TIME_STEP"		{return(TIME_STEP);}
"TIME_STEP_MAX"         {return(TIME_STEP_MAX);}
"TO"			{return(TO);}
"TOP"			{return(TOP);}
"TRAIN_DURATION"	{return(TRAIN_DURATION);}
"TRAIN_INTERVAL"	{return(TRAIN_INTERVAL);}
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
"VIZ_MESH_FORMAT"	{return(VIZ_MESH_FORMAT);}
"VIZ_MOLECULE_FORMAT"	{return(VIZ_MOLECULE_FORMAT);}
"VIZ_OUTPUT"	        {return(VIZ_OUTPUT);}
"VIZ_OUTPUT_REPORT"     { return VIZ_OUTPUT_REPORT; }
"VIZ_VALUE"	        {return(VIZ_VALUE);}
"VOLUME_DATA_OUTPUT"    {return VOLUME_DATA_OUTPUT; }
"VOLUME_OUTPUT_REPORT"  { return VOLUME_OUTPUT_REPORT; }
"VOLUME_DEPENDENT_RELEASE_NUMBER" {return(VOLUME_DEPENDENT_RELEASE_NUMBER);}
"VOLUME_ONLY"           {return(VOLUME_ONLY);}
"VOXEL_COUNT"           {return VOXEL_COUNT; }
"VOXEL_LIST"            {return(VOXEL_LIST);}
"VOXEL_SIZE"            {return VOXEL_SIZE; }
"WARNING"               {return(WARNING);}
"WARNINGS"              {return(WARNINGS);}
"WORLD"			{return(WORLD);}
"YES"			{return(YES);}
"printf"		{return(PRINTF);}
"fprintf"		{return(FPRINTF);}
"sprintf"		{return(SPRINTF);}
"print_time"		{return(PRINT_TIME);}
"fprint_time"		{return(FPRINT_TIME);}
"=fopen"		{return(FOPEN);}
"= fopen"		{return(FOPEN);}
"fclose"		{return(FCLOSE);}
{R}		        {sscanf((const char *)yytext,"%lf",&yylval->dbl); return REAL;}
{I}			{sscanf((const char *)yytext,"%d",&yylval->ival); return INTEGER;}
{IDENT}                 {
                             yylval->str = CHECKED_STRDUP(yytext, "string token");
                             return VAR;
                        }
\"(([^\"]*)|((\\\")*))*\"  {
                             yylval->str = CHECKED_STRDUP(yytext, "string token"); 
                             return STR_VALUE;
                           }

[\'\,\(\)\/\-\+\=\^\[\]\{\}\|\<\>\*\#\~\@\:\&\;] {return(yytext[0]);}
.			;
<INITIAL,IN_COMMENT>\n+ {parse_state->line_num[parse_state->include_stack_ptr - 1] += yyleng; }
<<EOF>>			{ yyterminate(); }
