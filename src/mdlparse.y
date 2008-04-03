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
  #include "mdlparse_aux.h"
  #include "mdlparse.h"
  #include "util.h"
  #include "react_output.h"
  #include "macromolecule.h"


  #define volp mdlpvp->vol

  typedef void *yyscan_t;
  int mdllex_init(yyscan_t *ptr_yy_globals) ;
  int mdllex_destroy(yyscan_t yyscanner);
  void mdlrestart(FILE *infile, yyscan_t scanner);
  int mdllex(YYSTYPE *yylval, struct mdlparse_vars *mdlpvp, yyscan_t scanner);

  static int mdlparse_file(struct mdlparse_vars *mpvp, char const *name);

#if 0
  #define CHECK(a)  do { if ((a) != 0) return 1; } while (0)
  #define CHECKN(a) do { if ((a) == NULL) return 1; } while (0)
#else
  #define CHECK(a)  do { if ((a) != 0)    { fprintf(stderr, "Fail: %s:%d\n", __FILE__, __LINE__); return 1; } } while (0)
  #define CHECKN(a) do { if ((a) == NULL) { fprintf(stderr, "Fail: %s:%d\n", __FILE__, __LINE__); return 1; } } while (0)
#endif

  #undef yyerror
  #define yyerror(a, b, c) mdlerror(a, c)
%}


%union {
int ival;
int tok;
double dbl;
char *str;
struct sym_table *sym;
struct vector3 *vec3;
struct num_expr_list_head nlist;
struct release_evaluator *rev;
struct sym_table_list *symlist;
struct output_times *otimes;

/* Reaction output */
struct output_times_inlist ro_otimes;
struct output_column_list ro_cols;
struct output_set *ro_set;
struct output_set_list ro_sets;
struct output_expression *cnt;

/* Viz output */
struct frame_data_list_head frame_list;

/* Region definitions */
struct element_list *elem_list_item;
struct element_list_head elem_list;
struct region *reg;

/* Diffusion constants */
struct diffusion_constant diff_const;

/* Geometry */
struct vertex_list_head vertlist;
struct vertex_list *vertlistitem;
struct element_connection_list_head ecl;
struct element_connection_list *elem_conn;
struct object *obj;
struct object_list obj_list;
struct voxel_object *voxel;

/* Molecule species */
struct species_opt_orient mol_type;
struct species_opt_orient_list mol_type_list;
struct species *mol_spec;
struct species_list species_lst;
struct species_list_item *species_lst_item;
struct eff_dat *effector;
struct eff_dat_list effector_list;

/* Reactions */
struct reaction_arrow react_arrow;
struct reaction_rate  react_rate;
struct reaction_rates react_rates;

/* Release sites/patterns */
struct release_pattern rpat;
struct release_single_molecule *rsm;
struct release_single_molecule_list rsm_list;

/* printf arguments */
struct arg *printfarg;
struct arg_list printfargs;

/* Macromolecules */
struct macro_topology *mmol_topo;
struct macro_subunit_spec *mmol_su_comp;
struct macro_subunit_assignment *mmol_su_assign;
struct macro_geometry *mmol_geom;
struct macro_relationship *mmol_su_rel;
struct macro_rate_ruleset *mmol_rate_ruleset;
struct macro_rate_rule *mmol_rate_rule;
struct macro_rate_clause *mmol_rate_clause;

/* Macromolecules (after the species is built) */
struct macro_relation_state *relation_state;
}

%pure_parser

%lex-param {struct mdlparse_vars *mdlpvp}
%lex-param {yyscan_t scanner}
%parse-param {struct mdlparse_vars *mdlpvp}
%parse-param {yyscan_t scanner}
%name-prefix="mdl"

%token       ABS
%token       ABSORPTIVE
%token       ACCURATE_3D_REACTIONS
%token       ACOS
%token       ALL_CROSSINGS
%token       ALL_DATA
%token       ALL_ELEMENTS
%token       ALL_ENCLOSED
%token       ALL_HITS
%token       ALL_ITERATIONS
%token       ALL_MESHES
%token       ALL_MOLECULES
%token       ALL_NOTIFICATIONS
%token       ALL_TIMES
%token       ALL_WARNINGS
%token       ASCII
%token       ASIN
%token       ASPECT_RATIO
%token       ATAN
%token       BACK
%token       BACK_CROSSINGS
%token       BACK_HITS
%token       BINARY
%token       BOTTOM
%token       BOX
%token       BOX_TRIANGULATION_REPORT
%token       BRIEF
%token       CEIL
%token       CENTER_MOLECULES_ON_GRID
%token       CHECKPOINT_INFILE
%token       CHECKPOINT_ITERATIONS
%token       CHECKPOINT_OUTFILE
%token       CHECKPOINT_REALTIME
%token       CLAMP_CONCENTRATION
%token       CLOSE_PARTITION_SPACING
%token       COMPLEX_RATE
%token       CONCENTRATION
%token       CORNERS
%token       COS
%token       COUNT
%token       CUBIC
%token       CUBIC_RELEASE_SITE
%token       CUSTOM_RK
%token       CUSTOM_SPACE_STEP
%token       CUSTOM_TIME_STEP
%token       DEFAULT
%token       DEFINE_COMPLEX_MOLECULE
%token       DEFINE_MOLECULE
%token       DEFINE_MOLECULES
%token       DEFINE_REACTIONS
%token       DEFINE_RELEASE_PATTERN
%token       DEFINE_SURFACE_CLASS
%token       DEFINE_SURFACE_CLASSES
%token       DEFINE_SURFACE_REGIONS
%token       DEGENERATE_POLYGONS
%token       DELAY
%token       DENSITY
%token       DIFFUSION_CONSTANT_2D
%token       DIFFUSION_CONSTANT_3D
%token       DIFFUSION_CONSTANT_REPORT
%token       DREAMM_V3
%token       DREAMM_V3_GROUPED
%token       DX
%token       EFFECTOR_GRID_DENSITY
%token       EFFECTOR_POSITIONS
%token       EFFECTOR_STATES
%token       ELEMENT_CONNECTIONS
%token       ELLIPTIC
%token       ELLIPTIC_RELEASE_SITE
%token       EQUAL
%token       ERROR
%token       ESTIMATE_CONCENTRATION
%token       EXCLUDE_ELEMENTS
%token       EXCLUDE_PATCH
%token       EXCLUDE_REGION
%token       EXIT
%token       EXP
%token       EXPRESSION
%token       FALSE
%token       FCLOSE
%token       FILENAME
%token       FILENAME_PREFIX
%token       FILE_OUTPUT_REPORT
%token       FINAL_SUMMARY
%token       FLOOR
%token       FOPEN
%token       FPRINTF
%token       FPRINT_TIME
%token       FRONT
%token       FRONT_CROSSINGS
%token       FRONT_HITS
%token       GAUSSIAN_RELEASE_NUMBER
%token       GEOMETRY
%token       HEADER
%token       HIGH_PROBABILITY_THRESHOLD
%token       HIGH_REACTION_PROBABILITY
%token       IGNORED
%token       INCLUDE_ELEMENTS
%token       INCLUDE_FILE
%token       INCLUDE_PATCH
%token       INCLUDE_REGION
%token       INPUT_FILE
%token       INSTANTIATE
%token <ival> INTEGER
%token       FULLY_RANDOM
%token       INTERACTION_RADIUS
%token       ITERATION_FRAME_DATA
%token       ITERATION_LIST
%token       ITERATION_NUMBERS
%token       ITERATION_REPORT
%token       ITERATIONS
%token       LEFT
%token       LIFETIME_THRESHOLD
%token       LIFETIME_TOO_SHORT
%token       LIST
%token       LOCATION
%token       LOG
%token       LOG10
%token       MAX_TOK
%token       MEAN_DIAMETER
%token       MEAN_NUMBER
%token       MESHES
%token       MICROSCOPIC_REVERSIBILITY
%token       MIN_TOK
%token       MISSED_REACTIONS
%token       MISSED_REACTION_THRESHOLD
%token       MISSING_SURFACE_ORIENTATION
%token       MOD
%token       MODE
%token       MODIFY_SURFACE_REGIONS
%token       MOLECULE
%token       MOLECULE_DENSITY
%token       MOLECULE_FILE_PREFIX
%token       MOLECULE_NUMBER
%token       MOLECULE_POSITIONS
%token       MOLECULES
%token       MOLECULE_STATES
%token       NAME_LIST
%token       NEGATIVE_DIFFUSION_CONSTANT
%token       NEGATIVE_REACTION_RATE
%token       NO
%token       NOEXIT
%token       NONE
%token       NORMAL
%token       NO_SPECIES
%token       NOT_EQUAL
%token       NOTIFICATIONS
%token       NUMBER_OF_SUBUNITS
%token       NUMBER_OF_TRAINS
%token       NUMBER_TO_RELEASE
%token       OBJECT
%token       OBJECT_FILE_PREFIXES
%token       OFF
%token       ON
%token       ORIENTATIONS
%token       OUTPUT_BUFFER_SIZE
%token       OVERWRITTEN_OUTPUT_FILE
%token       PARTITION_LOCATION_REPORT
%token       PARTITION_X
%token       PARTITION_Y
%token       PARTITION_Z
%token       PI_TOK
%token       POLYGON_LIST
%token       POSITIONS
%token       PRINTF
%token       PRINT_TIME
%token       PROBABILITY_REPORT
%token       PROBABILITY_REPORT_THRESHOLD
%token       PROGRESS_REPORT
%token       RADIAL_DIRECTIONS
%token       RADIAL_SUBDIVISIONS
%token       RAND_GAUSSIAN
%token       RAND_UNIFORM
%token       RATE_RULES
%token       REACTION_DATA_OUTPUT
%token       REACTION_GROUP
%token <dbl> REAL
%token       RECTANGULAR_RELEASE_SITE
%token       RECTANGULAR_TOKEN
%token       REFERENCE_DIFFUSION_CONSTANT
%token       REFLECTIVE
%token       REGION_DATA
%token       RELEASE_EVENT_REPORT
%token       RELEASE_INTERVAL
%token       RELEASE_PATTERN
%token       RELEASE_PROBABILITY
%token       RELEASE_SITE
%token       REMOVE_ELEMENTS
%token       RIGHT
%token       ROTATE
%token       ROUND_OFF
%token       SCALE
%token       SEED
%token       SHAPE
%token       SHOW_EXACT_TIME
%token       SIN
%token       SITE_DIAMETER
%token       SITE_RADIUS
%token       SPACE_STEP
%token       SPHERICAL
%token       SPHERICAL_RELEASE_SITE
%token       SPHERICAL_SHELL
%token       SPHERICAL_SHELL_SITE
%token       SPRINTF
%token       SQRT
%token       STANDARD_DEVIATION
%token       STATE_VALUES
%token       STEP
%token       STRING_TO_NUM
%token <str> STR_VALUE
%token       SUBUNIT
%token       SUBUNIT_RELATIONSHIPS
%token       SUMMATION_OPERATOR
%token       SURFACE_CLASS
%token       SURFACE_ONLY
%token       SURFACE_POSITIONS
%token       SURFACE_STATES
%token       TAN
%token       TARGET_ONLY
%token       TET_ELEMENT_CONNECTIONS
%token       THROUGHPUT_REPORT
%token       TIME_LIST
%token       TIME_POINTS
%token       TIME_STEP
%token       TIME_STEP_MAX
%token       TO
%token       TOP
%token       TRAIN_DURATION
%token       TRAIN_INTERVAL
%token       TRANSLATE
%token       TRANSPARENT
%token       TRIGGER
%token       TRUE
%token       UNLIMITED
%token       USELESS_VOLUME_ORIENTATION
%token       VACANCY_SEARCH_DISTANCE
%token <str> VAR
%token       VARYING_PROBABILITY_REPORT
%token       VERTEX_LIST
%token       VIZ_DATA_OUTPUT
%token       VIZ_MESH_FORMAT
%token       VIZ_MOLECULE_FORMAT
%token       VIZ_OUTPUT
%token       VIZ_VALUE
%token       VOLUME_DATA_OUTPUT
%token       VOLUME_DEPENDENT_RELEASE_NUMBER
%token       VOLUME_ONLY
%token       VOXEL_COUNT
%token       VOXEL_IMAGE_MODE
%token       VOXEL_LIST
%token       VOXEL_SIZE
%token       VOXEL_VOLUME_MODE
%token       WARNING
%token       WARNINGS
%token       WORLD
%token       YES

/* Utility non-terminals */
%type <str> file_name
%type <sym> existing_object
%type <symlist> mesh_object_or_wildcard
%type <sym> existing_region
%type <vec3> point
%type <vec3> point_or_num
%type <tok> boolean
%type <mol_type> orientation_class
%type <mol_type> list_orient_marks
%type <mol_type> orient_class_number
%type <nlist> list_range_specs
%type <nlist> range_spec

/* Assignment/expression non-terminals */
%type <sym> assign_var
%type <sym> existing_var_only
%type <nlist> array_value
%type <nlist> array_expr_only
%type <nlist> array_expr
%type <sym> existing_array
%type <dbl> num_expr
%type <dbl> num_value
%type <dbl> intOrReal
%type <dbl> num_expr_only
%type <sym> existing_num_var
%type <dbl> arith_expr
%type <str> str_expr
%type <sym> existing_str_var
%type <str> str_expr_only

/* I/O non-terminals */
%type <sym> new_file_stream
%type <sym> existing_file_stream
%type <str> file_mode
%type <str> format_string
%type <printfargs> list_args
%type <printfarg> list_arg

/* Notification non-terminals */
%type <tok> notify_bilevel
%type <tok> notify_level

/* Warning non-terminals */
%type <tok> warning_level

/* Checkpoint non-terminals */
%type <tok> exit_or_no
%type <dbl> time_expr

/* Partition non-terminals */
%type <tok> partition_dimension

/* Molecule definition non-terminals */
%type <species_lst> list_molecule_stmts
%type <mol_spec> molecule_stmt
%type <sym> new_molecule
%type <dbl> reference_diffusion_def
%type <diff_const> diffusion_def
%type <dbl> mol_timestep_def
%type <ival> target_def

/* Complex molecule definition non-terminals */
%type <str> complex_mol_name
%type <mmol_topo> complex_mol_topology
%type <mmol_su_comp> complex_mol_subunit_component
%type <mmol_su_comp> complex_mol_subunit_spec
%type <mmol_su_assign> complex_mol_subunit_assignment
%type <mmol_su_assign> complex_mol_subunits
%type <mmol_geom> complex_mol_geometry
%type <mmol_geom> complex_mol_subunit_locations
%type <mmol_geom> complex_mol_subunit_location
%type <mmol_su_rel> complex_mol_relationships
%type <mmol_su_rel> complex_mol_relationship_list
%type <mmol_su_rel> complex_mol_relationship
%type <mmol_rate_ruleset> complex_mol_rates
%type <mmol_rate_ruleset> complex_mol_rate_list
%type <mmol_rate_ruleset> complex_mol_rate
%type <mmol_rate_rule> complex_mol_rate_rules
%type <mmol_rate_rule> complex_mol_rate_rule
%type <mmol_rate_clause> complex_mol_rate_clause_list
%type <mmol_rate_clause> complex_mol_rate_clauses
%type <mmol_rate_clause> complex_mol_rate_clause
%type <ival> equal_or_not

/* Molecule utility non-terminals */
%type <sym> existing_molecule
%type <mol_type> existing_surface_molecule
%type <mol_type> existing_molecule_opt_orient
%type <sym> existing_macromolecule

/* Surface class non-terminals */
%type <sym> existing_surface_class
%type <tok> surface_rxn_type
%type <effector_list> surface_mol_stmt
%type <effector_list> list_surface_mol_density
%type <effector_list> list_surface_mol_num
%type <effector> surface_mol_quant

/* Reaction network non-terminals */
%type <react_arrow> right_cat_arrow
%type <react_arrow> double_cat_arrow
%type <react_arrow> reaction_arrow
%type <sym> new_rxn_pathname
%type <mol_type_list> reactant_list
%type <mol_type> reactant
%type <mol_type> existing_molecule_or_subunit
%type <mol_type> opt_reactant_surface_class
%type <mol_type> reactant_surface_class
%type <mol_type_list> product_list
%type <mol_type> product
%type <react_rates> rx_rate_syntax rx_dir_rate rx_rate1 rx_rate2
%type <react_rate> atomic_rate

/* Release pattern non-terminals */
%type <sym> new_release_pattern
%type <sym> existing_release_pattern_xor_rxpn
%type <rpat> list_req_release_pattern_cmds
%type <ival> train_count

/* Object type non-terminals */
%type <obj> object_def
%type <obj> meta_object_def
%type <obj> release_site_def release_site_def_new release_site_def_old
%type <obj> box_def
%type <obj> polygon_list_def
%type <obj> voxel_list_def
%type <sym> new_object
%type <obj_list> list_objects
%type <obj> object_ref
%type <obj> existing_object_ref

/* Release-site non-terminals */
%type <rev> release_region_expr
%type <tok> release_site_geom_old
%type <sym> existing_num_or_array
%type <tok> site_size_cmd
%type <rsm_list> molecule_release_pos_list
%type <rsm> molecule_release_pos

/* Polygon/voxel non-terminals */
%type <vertlist> vertex_list_cmd list_points
%type <vertlistitem> single_vertex
%type <ecl> element_connection_cmd tet_element_connection_cmd
%type <elem_conn> element_connection element_connection_tet
%type <ecl> list_element_connections list_tet_arrays

/* Region specification non-terminals */
%type <elem_list> remove_element_specifier_list
%type <tok> side_name
%type <elem_list> element_specifier_list element_specifier
%type <elem_list> incl_element_list_stmt excl_element_list_stmt
%type <elem_list> just_an_element_list
%type <elem_list> list_element_specs
%type <elem_list_item> element_spec
%type <elem_list_item> prev_region_stmt patch_statement
%type <tok> prev_region_type patch_type
%type <reg> new_region

/* Box non-terminals */
%type <dbl> opt_aspect_ratio_def

/* Reaction output non-terminals */
%type <dbl> output_buffer_size_def
%type <ro_otimes> output_timer_def
%type <ro_otimes> real_time_def iteration_time_def step_time_def
%type <ro_sets> list_count_cmds
%type <ro_set> count_cmd count_stmt
%type <str> custom_header_value
%type <ro_cols> single_count_expr list_count_exprs
%type <cnt> count_expr count_value
%type <tok> file_arrow
%type <str> outfile_syntax
%type <sym> existing_rxpn_or_molecule
%type <mol_type> existing_molecule_required_orient_braces
%type <cnt> count_syntax count_syntax_1 count_syntax_2 count_syntax_3
%type <cnt> count_syntax_macromol count_syntax_macromol_subunit
%type <relation_state> opt_macromol_relation_states
%type <relation_state> macromol_relation_state_list macromol_relation_state
%type <ival> macromol_relation_name
%type <sym> count_location_specifier
%type <tok> opt_hit_spec hit_spec
%type <str> opt_custom_header

/* Viz output non-terminals */
%type <ival> viz_mode_def
%type <ival> viz_mesh_format_def
%type <ival> viz_molecule_format_def
%type <frame_list> viz_frames_def
%type <tok> viz_meshes_one_item viz_molecules_one_item
%type <tok> iteration_frame_data_item
%type <frame_list> viz_molecules_block_def
%type <frame_list> list_viz_molecules_block_cmds viz_molecules_block_cmd
%type <symlist> existing_one_or_multiple_molecules
%type <nlist> viz_time_spec
%type <frame_list> viz_molecules_time_points_def viz_molecules_time_points_cmds
%type <frame_list> viz_molecules_time_points_one_cmd
%type <nlist> viz_iteration_spec
%type <frame_list> viz_molecules_iteration_numbers_def
%type <frame_list> viz_molecules_iteration_numbers_cmds
%type <frame_list> viz_molecules_iteration_numbers_one_cmd
%type <frame_list> viz_meshes_block_def
%type <frame_list> list_viz_meshes_block_cmds viz_meshes_block_cmd
%type <frame_list> viz_meshes_time_points_def viz_meshes_time_points_cmds
%type <frame_list> viz_meshes_time_points_one_cmd
%type <frame_list> viz_meshes_iteration_numbers_def
%type <frame_list> viz_meshes_iteration_numbers_cmds
%type <frame_list> viz_meshes_iteration_numbers_one_cmd
%type <frame_list> viz_frames_def_old viz_output_block_def viz_iteration_def
%type <frame_list> viz_time_def
%type <frame_list> viz_iteration_frame_data_def list_iteration_frame_data_specs
%type <frame_list> iteration_frame_data_spec
%type <sym> existing_logicalOrPhysical

/* Volume output non-terminals */
%type <str> volume_output_filename_prefix
%type <species_lst> volume_output_molecule_list volume_output_molecule_decl
%type <species_lst_item> volume_output_molecule
%type <species_lst>  volume_output_molecules
%type <vec3> volume_output_location volume_output_voxel_size
%type <vec3> volume_output_voxel_count
%type <otimes> volume_output_times_def

/* Operator associativities and precendences */
%right '='
%left '&' ':'
%left '+' '-'
%left '*' '/'
%left '^'
%left UNARYMINUS

%%

/* Start */
mdl_format:
        mdl_stmt_list
;

mdl_stmt_list:
        mdl_stmt
      | mdl_stmt_list mdl_stmt
;


mdl_stmt:
        include_stmt
      | assignment_stmt
      | io_stmt
      | notification_def
      | warnings_def
      | chkpt_stmt
      | parameter_def
      | partition_def
      | molecules_def
      | surface_classes_def
      | rx_net_def
      | release_pattern_def
      | physical_object_def
      | instance_def
      | existing_obj_define_surface_regions
      | mod_surface_regions
      | output_def
      | viz_output_def
      | viz_data_output_def
      | volume_output_def
;

/* =================================================================== */
/* Utility definitions */
file_name: str_expr
;

existing_object: VAR                                  { CHECKN($$ = mdl_existing_object(mdlpvp, $1)); }
;

mesh_object_or_wildcard: existing_object              { CHECKN($$ = mdl_singleton_symbol_list(mdlpvp, $1)); }
                       | STR_VALUE                    {
                                                        char *stripped;
                                                        CHECKN(stripped = mdl_strip_quotes(mdlpvp, $1));
                                                        CHECKN($$ = mdl_meshes_by_wildcard(mdlpvp, stripped));
                                                      }
;

existing_region: existing_object '[' VAR ']'          { CHECKN($$ = mdl_existing_region(mdlpvp, $1, $3)); }
;

point: array_value                                    {
                                                        if ($1.value_count != 3)
                                                        {
                                                          mdlerror(mdlpvp, "Three dimensional value required");
                                                          return 1;
                                                        }

                                                        if (($$ = (struct vector3 *)malloc(sizeof(struct vector3)))==NULL)
                                                        {
                                                          mdlerror(mdlpvp, "Out of memory while creating points");
                                                          return 1;
                                                        }

                                                        $$->x = $1.value_head->value;
                                                        $$->y = $1.value_head->next->value;
                                                        $$->z = $1.value_tail->value;
                                                      }
;

point_or_num: point
            | num_expr_only                           {
                                                        struct vector3 *vec;
                                                        if ((vec = (struct vector3 *)malloc(sizeof(struct vector3))) == NULL)
                                                        {
                                                          mdlerror(mdlpvp, "Out of memory while creating point");
                                                          return 1;
                                                        }
                                                        vec->x = $1;
                                                        vec->y = $1;
                                                        vec->z = $1;
                                                        $$ = vec;
                                                      }
;

boolean: TRUE                                         { $$ = 1; }
       | FALSE                                        { $$ = 0; }
       | YES                                          { $$ = 1; }
       | NO                                           { $$ = 0; }
       | ON                                           { $$ = 1; }
       | OFF                                          { $$ = 0; }
;

orientation_class: /* empty */                        { $$.orient_set = 0; }
                 | list_orient_marks
                 | orient_class_number
                 | ';'                                { $$.orient_set = 1; $$.orient = 0; }
;

list_orient_marks:
          head_mark                                   { $$.orient = 1; $$.orient_set = 1; }
        | tail_mark                                   { $$.orient = -1; $$.orient_set = 1; }
        | list_orient_marks head_mark                 { $$ = $1; ++ $$.orient; }
        | list_orient_marks tail_mark                 { $$ = $1; -- $$.orient; }
;

head_mark: '\''
;

tail_mark: ','
;

orient_class_number: '{' num_expr '}'                 {
                                                        $$.orient = (int) $2;
                                                        $$.orient_set = 1;
                                                        if ($$.orient != $2)
                                                        {
                                                          mdlerror(mdlpvp, "Error: Molecule orientation specified inside braces must be an integer between -32768 and 32767\n");
                                                          return 1;
                                                        }
                                                      }
;

list_range_specs:
          range_spec
        | list_range_specs ',' range_spec             {
                                                        if ($1.value_tail)
                                                        {
                                                          $$ = $1;
                                                          $$.value_count += $3.value_count;
                                                          $$.value_tail->next = $3.value_head;
                                                          $$.value_tail = $3.value_tail;
                                                        }
                                                        else
                                                          $$ = $3;
                                                      }
;

range_spec: num_expr                                  {
                                                        struct num_expr_list *nel = (struct num_expr_list *) malloc(sizeof(struct num_expr_list));
                                                        if (nel == NULL)
                                                        {
                                                          mdlerror(mdlpvp, "Out of memory while creating iteration list");
                                                          return 1;
                                                        }
                                                        $$.value_head = $$.value_tail = nel;
                                                        $$.value_count = 1;
                                                        $$.shared = 0;
                                                        $$.value_head->value = $1;
                                                        $$.value_head->next = NULL;
                                                      }
        | '[' num_expr TO num_expr STEP num_expr ']'  { CHECK(mdl_generate_range(mdlpvp, &$$, $2, $4, $6)); }
;


/* =================================================================== */
/* Include files */
include_stmt: INCLUDE_FILE '=' str_expr               {
                                                        char *include_path = mdl_find_include_file(mdlpvp, $3, volp->curr_file);
                                                        if (mdlparse_file(mdlpvp, include_path))
                                                        {
                                                          free(include_path);
                                                          return 1;
                                                        }
                                                        free(include_path);
                                                      }
;

/* =================================================================== */
/* Expressions */

assignment_stmt:
        assign_var '=' num_expr_only                  { CHECK(mdl_assign_variable_double(mdlpvp, $1, $3)); }
      | assign_var '=' str_expr_only                  { CHECK(mdl_assign_variable_string(mdlpvp, $1, $3)); }
      | assign_var '=' existing_var_only              { CHECK(mdl_assign_variable(mdlpvp, $1, $3)); }
      | assign_var '=' '[' array_expr ']'             { CHECK(mdl_assign_variable_array(mdlpvp, $1, $4.value_head)); }
;

assign_var: VAR                                       { CHECKN($$ = mdl_get_or_create_variable(mdlpvp, $1)); }
;

existing_var_only: VAR                                { CHECKN($$ = mdl_existing_variable(mdlpvp, $1)); }
;

array_value: array_expr_only
           | existing_array                           {
                                                        struct num_expr_list *elp;
                                                        $$.value_head = (struct num_expr_list *) $1->value;
                                                        $$.value_count = 1;
                                                        for (elp = $$.value_head; elp->next != NULL; elp = elp->next)
                                                          ++ $$.value_count;
                                                        $$.value_tail = elp;
                                                        $$.shared = 1;
                                                      }
;

array_expr_only: '[' array_expr ']'                   { mdl_debug_dump_array($2.head); $$ = $2; }
;

array_expr:
        num_expr                                      {
                                                        if (($$.value_tail = $$.value_head = (struct num_expr_list *)malloc
                                                            (sizeof(struct num_expr_list)))==NULL) {
                                                          mdlerror(mdlpvp, "Out of memory while creating numerical array");
                                                          return 1;
                                                        }
                                                        $$.value_tail->value = $1;
                                                        $$.value_tail->next = NULL;
                                                        $$.value_count = 1;
                                                        $$.shared = 0;
                                                      }
      | array_expr ',' num_expr                       {
                                                        $$ = $1;
                                                        if (($$.value_tail = $$.value_tail->next = (struct num_expr_list *)malloc
                                                            (sizeof(struct num_expr_list)))==NULL) {
                                                          mdlerror(mdlpvp, "Out of memory while creating numerical array");
                                                          return 1;
                                                        }
                                                        $$.value_tail->value = $3;
                                                        $$.value_tail->next = NULL;
                                                        ++ $$.value_count;
                                                        $$.shared = 0;
                                                      }
;

existing_array: VAR                                   { CHECKN($$ = mdl_existing_array(mdlpvp, $1)); }
;

num_expr: num_value
        | arith_expr
;

num_value: intOrReal
         | existing_num_var                           { $$ = *(double *) $1->value; }
;

intOrReal: INTEGER                                    { $$ = $1; }
         | REAL
;

num_expr_only: intOrReal
             | arith_expr
;

existing_num_var: VAR                                 { CHECKN($$ = mdl_existing_double(mdlpvp, $1)); }
;

arith_expr:
        '(' num_expr ')'                              { $$ = $2; }
      | EXP '(' num_expr ')'                          { $$ = exp($3); }
      | LOG '(' num_expr ')'                          { $$ = log($3); }
      | LOG10 '(' num_expr ')'                        { $$ = log10($3); }
      | MAX_TOK '(' num_expr ',' num_expr ')'         { $$ = max2d($3, $5); }
      | MIN_TOK '(' num_expr ',' num_expr ')'         { $$ = min2d($3, $5); }
      | ROUND_OFF '(' num_expr ',' num_expr ')'       {
                                                        char fmt_string[1024];
                                                        fmt_string[0] = '\0';
                                                        snprintf(fmt_string, 1024, "%.*g", (int) $3, $5);
                                                        $$=strtod(fmt_string, (char **)NULL);
                                                      }
      |  FLOOR '(' num_expr ')'                       { $$ = floor($3); }
      |  CEIL '(' num_expr ')'                        { $$ = ceil($3); }
      | SIN '(' num_expr ')'                          { $$ = sin($3); }
      | COS '(' num_expr ')'                          { $$ = cos($3); }
      | TAN '(' num_expr ')'                          { $$ = tan($3); }
      | ASIN '(' num_expr ')'                         { $$ = asin($3); }
      | ACOS '(' num_expr ')'                         { $$ = acos($3); }
      | ATAN '(' num_expr ')'                         { $$ = atan($3); }
      | SQRT '(' num_expr ')'                         { $$ = sqrt($3); }
      | ABS '(' num_expr ')'                          { $$ = fabs($3); }
      | MOD '(' num_expr ',' num_expr ')'             { $$ = fmod($3, $5); }
      | PI_TOK                                        { $$ = MY_PI; }
      | RAND_UNIFORM                                  {
                                                        if (mdlpvp->vol->notify->final_summary == NOTIFY_FULL) {
                                                          mdlpvp->vol->random_number_use++;
                                                        }
                                                        $$ = rng_dbl(mdlpvp->vol->rng);
                                                      }
      | RAND_GAUSSIAN                                 { $$ = rng_gauss(mdlpvp->vol->rng); }
      | SEED                                          { $$ = volp->seed_seq; }
      | STRING_TO_NUM '(' str_expr ')'                {
                                                        $$ = strtod($3, (char **)NULL);
                                                        if (errno==ERANGE)
                                                        {
                                                          mdlerror_fmt(mdlpvp, "Error converting string to number: %s", $3);
                                                          free($3);
                                                          return 1;
                                                        }
                                                        free($3);
                                                      }
      | num_expr '+' num_expr                         { $$ = $1 + $3; }
      | num_expr '-' num_expr                         { $$ = $1 - $3; }
      | num_expr '*' num_expr                         { $$ = $1 * $3; }
      | num_expr '/' num_expr                         { $$ = $1 / $3; }
      | num_expr '^' num_expr                         { $$ = pow($1, $3); }
      | '-' num_expr %prec UNARYMINUS                 { $$ = -$2; }
      | '+' num_expr %prec UNARYMINUS                 { $$ = $2; }
;

str_expr:
        str_expr_only
      | existing_str_var                              { CHECKN($$ = mdl_strdup(mdlpvp, (char const *) $1->value)); }
;

str_expr_only:
        STR_VALUE                                     { CHECKN($$ = mdl_strip_quotes(mdlpvp, $1)); }
      | INPUT_FILE                                    { CHECKN($$ = mdl_strdup(mdlpvp, volp->mdl_infile_name)); }
      | str_expr '&' str_expr                         { CHECKN($$ = mdl_strcat(mdlpvp, $1, $3)); }
;

existing_str_var: VAR                                 { CHECKN($$ = mdl_existing_string(mdlpvp, $1)); }
;

/* =================================================================== */
/* I/O */

io_stmt: fopen_stmt
       | fclose_stmt
       | printf_stmt
       | fprintf_stmt
       | sprintf_stmt
       | print_time_stmt
       | fprint_time_stmt
;

fopen_stmt: new_file_stream FOPEN
            '(' file_name ',' file_mode ')'           {
                                                        struct file_stream *filep = (struct file_stream *) $1->value;
                                                        filep->name=$4;
                                                        if ((filep->stream=fopen(filep->name, $6)) == NULL)
                                                        {
                                                          mdlerror_fmt(mdlpvp, "Cannot open file: %s", filep->name);
                                                          /* XXX: Free filep? */
                                                          return 1;
                                                        }
                                                      }
;

new_file_stream: VAR                                  {
                                                        if (($$ = store_sym($1, FSTRM, volp->main_sym_table, NULL)) == NULL)
                                                        {
                                                          mdlerror_fmt(mdlpvp, "Out of memory while creating file stream: %s", $1);
                                                          free($1);
                                                          return 1;
                                                        }
                                                        free($1);
                                                      }
;

file_mode: str_expr                                   { $$ = $1; CHECK(mdl_valid_file_mode(mdlpvp, $1)); }
;

fclose_stmt: FCLOSE '(' existing_file_stream ')'      {
                                                        struct file_stream *filep=(struct file_stream *) $3->value;
                                                        if (fclose(filep->stream)!=0) {
                                                          mdlerror_fmt(mdlpvp, "Error closing file: %s", filep->name);
                                                          return 1;
                                                        }
                                                      }
;

existing_file_stream: VAR                             { CHECKN($$ = mdl_existing_file_stream(mdlpvp, $1)); }
;

format_string: str_expr                               { CHECKN($$ = mdl_expand_string_escapes(mdlpvp, $1)); }
;

list_args: /* empty */                                { $$.arg_head = $$.arg_tail = NULL; }
         | list_args ',' list_arg                     {
                                                        $$ = $1;
                                                        if ($$.arg_tail)
                                                          $$.arg_tail = $$.arg_tail->next = $3;
                                                        else
                                                          $$.arg_tail = $$.arg_head = $3;
                                                        $3->next = NULL;
                                                      }
;

list_arg: num_expr_only                               { CHECKN($$ = mdl_new_printf_arg_double(mdlpvp, $1)); }
        | str_expr_only                               { CHECKN($$ = mdl_new_printf_arg_string(mdlpvp, $1)); }
        | existing_var_only                           {
                                                        switch ($1->sym_type)
                                                        {
                                                          case DBL: CHECKN($$ = mdl_new_printf_arg_double(mdlpvp, *(double *) $1->value)); break;
                                                          case STR: CHECKN($$ = mdl_new_printf_arg_string(mdlpvp, (char *) $1->value)); break;
                                                          default:
                                                            mdlerror(mdlpvp, "Invalid variable type referenced");
                                                            return 1;
                                                        }
                                                      }
;

printf_stmt: PRINTF '(' format_string list_args ')'   { CHECK(mdl_printf(mdlpvp, $3, $4.arg_head)); }
;

fprintf_stmt:
          FPRINTF '('
            existing_file_stream ','
            format_string list_args ')'               { CHECK(mdl_fprintf(mdlpvp, (struct file_stream *) $3->value, $5, $6.arg_head)); }
;

sprintf_stmt:
          SPRINTF '('
            assign_var ','
            format_string list_args ')'               { CHECK(mdl_sprintf(mdlpvp, $3, $5, $6.arg_head)); }
;

print_time_stmt: PRINT_TIME '(' format_string ')'     { mdl_print_time(mdlpvp, $3); }
;

fprint_time_stmt:
          FPRINT_TIME '('
            existing_file_stream ','
            format_string ')'                         { CHECK(mdl_fprint_time(mdlpvp, $3, $5)); }
;

/* =================================================================== */
/* Notifications */

notification_def:
        NOTIFICATIONS '{' notification_list '}'
;

notification_list:
        notification_item_def
      | notification_list notification_item_def
;

notification_item_def:
        ALL_NOTIFICATIONS '=' notify_bilevel          { mdl_set_all_notifications(mdlpvp->vol, $3); }
      | PROGRESS_REPORT '=' notify_bilevel            { mdlpvp->vol->notify->progress_report        = $3; }
      | DIFFUSION_CONSTANT_REPORT '=' notify_level    { mdlpvp->vol->notify->diffusion_constants    = $3; }
      | PROBABILITY_REPORT '=' notify_bilevel         { mdlpvp->vol->notify->reaction_probabilities = $3; }
      | VARYING_PROBABILITY_REPORT '=' notify_bilevel { mdlpvp->vol->notify->time_varying_reactions = $3; }
      | PROBABILITY_REPORT_THRESHOLD '=' num_expr     { mdlpvp->vol->notify->reaction_prob_notify   = $3; }
      | PARTITION_LOCATION_REPORT '=' notify_bilevel  { mdlpvp->vol->notify->partition_location     = $3; }
      | BOX_TRIANGULATION_REPORT '=' notify_bilevel   { mdlpvp->vol->notify->box_triangulation      = $3; }
      | RELEASE_EVENT_REPORT '=' notify_bilevel       { mdlpvp->vol->notify->release_events         = $3; }
      | FILE_OUTPUT_REPORT '=' notify_bilevel         { mdlpvp->vol->notify->file_writes            = $3; }
      | FINAL_SUMMARY '=' notify_bilevel              { mdlpvp->vol->notify->final_summary          = $3; }
      | THROUGHPUT_REPORT '=' notify_bilevel          { mdlpvp->vol->notify->throughput_report      = $3; }
      | ITERATION_REPORT '=' notify_bilevel           {
                                                        /* Only if not set on command line */
                                                        if (mdlpvp->vol->log_freq == -1)
                                                        {
                                                          mdlpvp->vol->notify->custom_iterations = $3;
                                                        }
                                                      }
      | ITERATION_REPORT '=' num_expr                 { CHECK(mdl_set_iteration_report_freq(mdlpvp, (long long) $3)); }
;

notify_bilevel:
        boolean                                       { $$ = ($1 ? NOTIFY_FULL : NOTIFY_NONE); }
;

notify_level:
        boolean                                       { $$ = ($1 ? NOTIFY_FULL : NOTIFY_NONE); }
      | BRIEF                                         { $$ = NOTIFY_BRIEF; }
;

/* =================================================================== */
/* Warnings */

warnings_def:
        WARNINGS '{' warning_list '}'
;

warning_list:
        warning_item_def
      | warning_list warning_item_def
;

warning_item_def:
        ALL_WARNINGS '=' warning_level                { mdl_set_all_warnings(mdlpvp->vol, (byte) $3); }
      | NEGATIVE_DIFFUSION_CONSTANT '=' warning_level { mdlpvp->vol->notify->neg_diffusion = (byte)$3; }
      | NEGATIVE_REACTION_RATE '=' warning_level      { mdlpvp->vol->notify->neg_reaction = (byte)$3; }
      | HIGH_REACTION_PROBABILITY '=' warning_level   { mdlpvp->vol->notify->high_reaction_prob = (byte)$3; }
      | HIGH_PROBABILITY_THRESHOLD '=' num_expr       { mdlpvp->vol->notify->reaction_prob_warn = $3; }
      | CLOSE_PARTITION_SPACING '=' warning_level     { mdlpvp->vol->notify->close_partitions = (byte)$3; }
      | DEGENERATE_POLYGONS '=' warning_level         { mdlpvp->vol->notify->degenerate_polys = (byte)$3; }
      | OVERWRITTEN_OUTPUT_FILE '=' warning_level     { mdlpvp->vol->notify->overwritten_file = (byte)$3; }
      | LIFETIME_TOO_SHORT '=' warning_level          { mdlpvp->vol->notify->short_lifetime = (byte)$3; }
      | LIFETIME_THRESHOLD '=' num_expr               { CHECK(mdl_set_lifetime_warning_threshold(mdlpvp, (long long) $3)); }
      | MISSED_REACTIONS '=' warning_level            { mdlpvp->vol->notify->missed_reactions = (byte)$3; }
      | MISSED_REACTION_THRESHOLD '=' num_expr        { CHECK(mdl_set_missed_reaction_warning_threshold(mdlpvp, $3)); }
      | MISSING_SURFACE_ORIENTATION '=' warning_level { mdlpvp->vol->notify->missed_surf_orient = (byte)$3; }
      | USELESS_VOLUME_ORIENTATION '=' warning_level  { mdlpvp->vol->notify->useless_vol_orient = (byte)$3; }
;

warning_level:
        IGNORED                                       { $$ = WARN_COPE;  }
      | WARNING                                       { $$ = WARN_WARN;  }
      | ERROR                                         { $$ = WARN_ERROR; }
;

/* =================================================================== */
/* Checkpoint configuration */

chkpt_stmt: CHECKPOINT_INFILE '=' file_name           {
                                                        FILE *file;
                                                        if (volp->chkpt_infile == NULL)
                                                        {
                                                          volp->chkpt_infile=$3;
                                                          if ((file = fopen(volp->chkpt_infile,"r")) == NULL)
                                                          {
                                                            volp->chkpt_init=1;
                                                          }
                                                          else
                                                          {
                                                            volp->chkpt_init=0;
                                                            fclose(file);
                                                          }
                                                          volp->chkpt_flag = 1;
                                                        }
                                                      }
        | CHECKPOINT_OUTFILE '=' file_name            { volp->chkpt_outfile=$3; volp->chkpt_flag = 1; }
        | CHECKPOINT_ITERATIONS '=' num_expr          {
                                                        volp->chkpt_iterations = (long long) $3;
                                                        if(volp->chkpt_iterations <= 0)
                                                        {
                                                           mdlerror(mdlpvp, "Error: CHECKPOINT_ITERATIONS must be a positive integer\n");
                                                           return 1;
                                                        }
                                                        volp->chkpt_flag = 1;
                                                      }
        | CHECKPOINT_REALTIME '='
          time_expr exit_or_no                        { mdl_set_realtime_checkpoint(mdlpvp, (long) $3, $4); }
;

exit_or_no:                                           { $$ = 0; }
          | NOEXIT                                    { $$ = 1; }
          | EXIT                                      { $$ = 0; }
;

time_expr:
          num_expr                                    { /* seconds */     $$ = $1; }
        | num_expr ':' num_expr                       { /* mm:ss */       $$ = $1 * 60 + $3; }
        | num_expr ':' num_expr ':' num_expr          { /* hh:mm:ss */    $$ = $1 * 3600 + $3 * 60 + $5; }
        | num_expr ':' num_expr ':' num_expr
          ':' num_expr                                { /* dd:hh:mm:ss */ $$ = $1 * 86400 + $3 * 3600 + $5 * 60 + $7; }
;

/* =================================================================== */
/* Simulation parameters */

parameter_def:
          TIME_STEP '=' num_expr                      { mdl_set_time_step(mdlpvp, $3); }
        | SPACE_STEP '=' num_expr                     { mdl_set_space_step(mdlpvp, $3); }
        | TIME_STEP_MAX '=' num_expr                  { mdl_set_max_time_step(mdlpvp, $3); }
        | ITERATIONS '=' num_expr                     { CHECK(mdl_set_num_iterations(mdlpvp, (long long) $3)); }
        | CENTER_MOLECULES_ON_GRID '=' boolean        { mdlpvp->vol->randomize_gmol_pos = !($3); }
        | ACCURATE_3D_REACTIONS '=' boolean           { mdlpvp->vol->use_expanded_list = $3; }
        | VACANCY_SEARCH_DISTANCE '=' num_expr        { mdlpvp->vol->vacancy_search_dist2 = max2d($3, 0.0); }
        | RADIAL_DIRECTIONS '=' num_expr              { CHECK(mdl_set_num_radial_directions(mdlpvp, (int) $3)); }
        | RADIAL_DIRECTIONS '=' FULLY_RANDOM          { volp->fully_random = 1; }
        | RADIAL_SUBDIVISIONS '=' num_expr            { CHECK(mdl_set_num_radial_directions(mdlpvp, (int) $3)); }
        | EFFECTOR_GRID_DENSITY '=' num_expr          { CHECK(mdl_set_grid_density(mdlpvp, $3)); }
        | INTERACTION_RADIUS '=' num_expr             { volp->rx_radius_3d = $3; }
        | MICROSCOPIC_REVERSIBILITY '=' boolean       { mdlpvp->vol->surface_reversibility=$3; mdlpvp->vol->volume_reversibility=$3; }
        | MICROSCOPIC_REVERSIBILITY '=' SURFACE_ONLY  { mdlpvp->vol->surface_reversibility=1;  mdlpvp->vol->volume_reversibility=0;  }
        | MICROSCOPIC_REVERSIBILITY '=' VOLUME_ONLY   { mdlpvp->vol->surface_reversibility=0;  mdlpvp->vol->volume_reversibility=1;  }
;

/* =================================================================== */
/* Partitions */

partition_def:
          partition_dimension '='
          '[' list_range_specs ']'                    { CHECK(mdl_set_partition(mdlpvp, $1, $4.value_head, $4.value_count)); }
;

partition_dimension:
          PARTITION_X                                 { $$ = X_PARTS; }
        | PARTITION_Y                                 { $$ = Y_PARTS; }
        | PARTITION_Z                                 { $$ = Z_PARTS; }
;

/* =================================================================== */
/* Molecule definitions */

molecules_def:
          define_one_molecule
        | define_multiple_molecules
        | define_complex_molecule
;

define_one_molecule: DEFINE_MOLECULE molecule_stmt    {
                                                        if (volp->procnum == 0)
                                                        {
                                                          if (volp->notify->diffusion_constants == NOTIFY_BRIEF)
                                                            fprintf(volp->log_file,"Defining molecule with the following diffusion constant:\n");
                                                          mdl_report_diffusion_distances(volp->log_file, $2, volp->time_unit, volp->length_unit, volp->notify->diffusion_constants);
                                                          no_printf("Molecule %s defined with D = %g\n", $2->sym->name, $2->D);
                                                          if (volp->notify->diffusion_constants == NOTIFY_BRIEF)
                                                            fprintf(volp->log_file, "\n");
                                                        }
                                                      }
;

define_multiple_molecules: DEFINE_MOLECULES
                          '{'
                              list_molecule_stmts
                          '}'                         {
                                                        if (volp->procnum == 0)
                                                        {
                                                          struct species_list_item *ptrl;
                                                          if (volp->notify->diffusion_constants == NOTIFY_BRIEF)
                                                            fprintf(volp->log_file,"Defining molecules with the following theoretical average diffusion distances:\n");
                                                          for (ptrl = $3.species_head; ptrl != NULL; ptrl = ptrl->next)
                                                          {
                                                            struct species *spec = (struct species *) ptrl->spec;
                                                            mdl_report_diffusion_distances(volp->log_file, spec, volp->time_unit, volp->length_unit, volp->notify->diffusion_constants);
                                                            no_printf("Molecule %s defined with D = %g\n", spec->sym->name, spec->D);
                                                          }
                                                          if (volp->notify->diffusion_constants == NOTIFY_BRIEF)
                                                            fprintf(volp->log_file,"\n");
                                                        }
                                                        mem_put_list(mdlpvp->species_list_mem, $3.species_head);
                                                      }
;

list_molecule_stmts:
          molecule_stmt                               {
                                                        struct species_list_item *ptrl;
                                                        ptrl = (struct species_list_item *) mem_get(mdlpvp->species_list_mem);
                                                        if (ptrl == NULL) {
                                                          mdlerror_fmt(mdlpvp, "Out of memory while parsing molecules");
                                                          return 1;
                                                        }
                                                        ptrl->spec = $1;
                                                        ptrl->next = NULL;
                                                        $$.species_tail = $$.species_head = ptrl;
                                                        $$.species_count = 1;
                                                      }
        | list_molecule_stmts molecule_stmt           {
                                                        struct species_list_item *ptrl;
                                                        ptrl = (struct species_list_item *) mem_get(mdlpvp->species_list_mem);
                                                        if (ptrl == NULL) {
                                                          mdlerror_fmt(mdlpvp, "Out of memory while parsing molecules");
                                                          return 1;
                                                        }
                                                        ptrl->spec = $2;
                                                        ptrl->next = NULL;
                                                        $$ = $1;
                                                        $$.species_tail = $$.species_tail->next = ptrl;
                                                        ++ $$.species_count;
                                                      }
;

molecule_stmt:
          new_molecule '{'
              reference_diffusion_def
              diffusion_def
              mol_timestep_def
              target_def
          '}'                                         { CHECKN($$ = mdl_assemble_mol_species(mdlpvp, $1, $3, $4.D, $4.is_2d, $5, $6)); }
;

new_molecule: VAR                                     { CHECKN($$ = mdl_new_molecule(mdlpvp, $1)); }
;

reference_diffusion_def:
          /* empty */                                 { $$ = 0; }
        | REFERENCE_DIFFUSION_CONSTANT '=' num_expr   { $$ = $3; }
;

diffusion_def:
          DIFFUSION_CONSTANT_3D '=' num_expr          { $$.is_2d = 0; $$.D = $3; CHECK(mdl_check_diffusion_constant(mdlpvp, & $$.D)); }
        | DIFFUSION_CONSTANT_2D '=' num_expr          { $$.is_2d = 1; $$.D = $3; CHECK(mdl_check_diffusion_constant(mdlpvp, & $$.D)); }
;

mol_timestep_def:
          /* empty */                                 { $$ = 0.0; }
        | CUSTOM_TIME_STEP '=' num_expr               {
                                                        if ($3 > 0)
                                                          $$ = $3;
                                                        else
                                                        {
                                                          mdlerror(mdlpvp, "Zero or negative custom time step is disallowed (ignoring).");
                                                          $$ = 0.0;
                                                        }
                                                      }
        | CUSTOM_SPACE_STEP '=' num_expr              {
                                                        if ($3 > 0)
                                                          $$ = - $3;
                                                        else
                                                        {
                                                          mdlerror(mdlpvp, "Zero or negative custom space step is disallowed (ignoring).");
                                                          $$ = 0.0;
                                                        }
                                                      }
;

target_def: /* empty */                               { $$ = 0; }
          | TARGET_ONLY                               { $$ = 1; }
;

define_complex_molecule:
          DEFINE_COMPLEX_MOLECULE complex_mol_name    { mdlpvp->complex_name = $2; mdlpvp->complex_type = 0; }
          '{'
              complex_mol_topology                    { mdlpvp->complex_topo = $5; }
              complex_mol_subunits
              complex_mol_geometry
              complex_mol_relationships               { mdlpvp->complex_relations = $9; }
              complex_mol_rates
          '}'                                         { CHECK(mdl_assemble_complex_species(mdlpvp, $2, $5, $7, $8, $9, $11)); }
;

complex_mol_name: VAR                                 /* XXX: need check for unique name here? */
;

complex_mol_topology:
          NUMBER_OF_SUBUNITS '=' array_value          { CHECKN($$ = mdl_assemble_topology(mdlpvp, &$3)); }
;

complex_mol_subunits:
          complex_mol_subunit_assignment
        | complex_mol_subunits
          complex_mol_subunit_assignment              { $2->next = $1; $$ = $2; }
;

complex_mol_subunit_assignment:
          SUBUNIT '[' complex_mol_subunit_spec ']'
          '=' existing_molecule_opt_orient            { CHECKN($$ = mdl_assemble_complex_subunit_assignment(mdlpvp, $3, & $6)); }
;

complex_mol_subunit_spec:
          complex_mol_subunit_component
        | complex_mol_subunit_spec ','
          complex_mol_subunit_component               { $3->next = $1; $$ = $3; }
;

complex_mol_subunit_component: num_expr               { CHECKN($$ = mdl_assemble_subunit_spec_component(mdlpvp, $1, $1)); }
                             | num_expr ':' num_expr  { CHECKN($$ = mdl_assemble_subunit_spec_component(mdlpvp, $1, $3)); }
;

complex_mol_geometry:
          SHAPE '{' complex_mol_subunit_locations '}' { $$ = $3; CHECK(mdl_validate_complex_geometry(mdlpvp, mdlpvp->complex_topo, $3)); }
;

complex_mol_subunit_locations:
          complex_mol_subunit_location
        | complex_mol_subunit_locations
          complex_mol_subunit_location                { $2->next = $1; $$ = $2; }
;

complex_mol_subunit_location:
          SUBUNIT '[' array_expr ']' '=' point        { CHECKN($$ = mdl_assemble_complex_geometry(mdlpvp, mdlpvp->complex_topo, &$3, $6)); }
;

complex_mol_relationships:
          SUBUNIT_RELATIONSHIPS
          '{' complex_mol_relationship_list '}'       { $$ = $3; CHECK(mdl_validate_complex_relationships(mdlpvp, mdlpvp->complex_topo, $3)); }
;

complex_mol_relationship_list:
          /* empty */                                 { $$ = NULL; }
        | complex_mol_relationship_list
          complex_mol_relationship                    { $2->next = $1; $$ = $2; }
;

complex_mol_relationship: VAR '=' array_value         { CHECKN($$ = mdl_assemble_complex_relationship(mdlpvp, mdlpvp->complex_topo, $1, &$3)); }
;

complex_mol_rates:
          RATE_RULES '{' complex_mol_rate_list '}'    { $$ = $3; CHECK(mdl_validate_complex_rates(mdlpvp, $3)); }
;

complex_mol_rate_list:
          /* empty */                                 { $$ = NULL; }
        | complex_mol_rate_list complex_mol_rate      { $2->next = $1; $$ = $2; }
;

complex_mol_rate: VAR '{' complex_mol_rate_rules '}'  { CHECKN($$ = mdl_assemble_complex_ruleset(mdlpvp, $1, $3)); }
;

complex_mol_rate_rules:
          complex_mol_rate_rule
        | complex_mol_rate_rules
          complex_mol_rate_rule                       { $2->next = $1; $$ = $2; }
;

complex_mol_rate_rule:
          complex_mol_rate_clauses ':' num_expr       { CHECKN($$ = mdl_assemble_complex_rate_rule(mdlpvp, $1, $3)); }
;

complex_mol_rate_clauses:
          complex_mol_rate_clause_list
        | DEFAULT                                     { $$ = NULL; }
;

complex_mol_rate_clause_list:
          complex_mol_rate_clause
        | complex_mol_rate_clause_list '&'
          complex_mol_rate_clause                     { $3->next = $1; $$ = $3; }
;

complex_mol_rate_clause:
     VAR equal_or_not existing_molecule_opt_orient    { CHECKN($$ = mdl_assemble_complex_rate_rule_clause(mdlpvp, mdlpvp->complex_relations, $1, $2, &$3)); }
;

equal_or_not: EQUAL                                   { $$ = 0; }
            | NOT_EQUAL                               { $$ = 1; }
;

existing_molecule: VAR                                { CHECKN($$ = mdl_existing_molecule(mdlpvp, $1)); }
;

existing_surface_molecule:
          existing_molecule_opt_orient                {
                                                        $$ = $1;
                                                        struct species *specp=(struct species *) $$.mol_type->value;
                                                        if ((specp->flags & ON_GRID) == 0)
                                                        {
                                                          mdlerror_fmt(mdlpvp, "Invalid surface molecule specified: %s", specp->sym->name);
                                                          return 1;
                                                        }
                                                      }
;

existing_molecule_opt_orient:
          existing_molecule orientation_class         {
                                                        $$ = $2;
                                                        if (! $$.orient_set)
                                                          $$.orient = 0;
                                                        $$.mol_type = $1;
                                                      }
;

existing_macromolecule: VAR                           { CHECKN($$ = mdl_existing_macromolecule(mdlpvp, $1)); }
;

/* =================================================================== */
/* Surface class definitions */

surface_classes_def:
          define_one_surface_class
        | define_multiple_surface_classes
;


define_one_surface_class:
          DEFINE_SURFACE_CLASS surface_class_stmt
;

define_multiple_surface_classes:
          DEFINE_SURFACE_CLASSES '{'
            list_surface_class_stmts
          '}'
;

list_surface_class_stmts:
          surface_class_stmt
        | list_surface_class_stmts
          surface_class_stmt
;

surface_class_stmt:
          new_molecule '{'                            {
                                                        struct species *specp = (struct species *) $1->value;
                                                        specp->flags = IS_SURFACE;
                                                        mdlpvp->current_surface_class = specp;
                                                      }
            list_surface_prop_stmts
          '}'                                         { mdlpvp->current_surface_class = NULL; }
;

existing_surface_class: VAR                           { CHECKN($$ = mdl_existing_surface_class(mdlpvp, $1)); }
;

list_surface_prop_stmts:
          /* empty */
        | list_surface_prop_stmts
          surface_prop_stmt
;

surface_prop_stmt:
          surface_rxn_stmt
        | surface_class_mol_stmt
        | surface_class_viz_value_stmt
;

surface_rxn_stmt:
          surface_rxn_type
            equals_or_to
            existing_molecule_opt_orient              { CHECKN(mdl_assemble_surface_reaction(mdlpvp, $1, mdlpvp->current_surface_class, $3.mol_type, $3.orient)); }
        | CLAMP_CONCENTRATION
            existing_molecule_opt_orient '='
            num_expr                                  { CHECKN(mdl_assemble_concentration_clamp_reaction(mdlpvp, mdlpvp->current_surface_class, $2.mol_type, $2.orient, $4)); }
;

surface_rxn_type: REFLECTIVE                          { $$ = RFLCT; }
                | TRANSPARENT                         { $$ = TRANSP; }
                | ABSORPTIVE                          { $$ = SINK; }
;

equals_or_to: '='
            | TO
;

surface_class_mol_stmt: surface_mol_stmt              { mdlpvp->current_surface_class->eff_dat_head = $1.eff_head; }
;

surface_mol_stmt:
          MOLECULE_DENSITY
          '{'
              list_surface_mol_density
          '}'                                         { $$ = $3; }
        | MOLECULE_NUMBER
          '{'
              list_surface_mol_num
          '}'                                         { $$ = $3; }
;

list_surface_mol_density:
          surface_mol_quant                           {
                                                        $1->quantity_type = EFFDENS;
                                                        $$.eff_tail = $$.eff_head = $1;
                                                      }
        | list_surface_mol_density
          surface_mol_quant                           {
                                                        $$ = $1;
                                                        $2->quantity_type = EFFDENS;
                                                        $$.eff_tail = $$.eff_tail->next = $2;
                                                      }
;

list_surface_mol_num:
          surface_mol_quant                           {
                                                        $1->quantity_type = EFFNUM;
                                                        $$.eff_tail = $$.eff_head = $1;
                                                      }
        | list_surface_mol_num
          surface_mol_quant                           {
                                                        $$ = $1;
                                                        $2->quantity_type = EFFNUM;
                                                        $$.eff_tail = $$.eff_tail->next = $2;
                                                      }
;

surface_mol_quant:
          existing_surface_molecule '=' num_expr      { CHECKN($$ = mdl_new_effector_data(mdlpvp, &$1, $3)); }
;

surface_class_viz_value_stmt:
          VIZ_VALUE '=' num_expr                      { mdlpvp->current_surface_class->region_viz_value = (int) $3; }
;

/* =================================================================== */
/* Reaction network definitions */

rx_net_def:
          DEFINE_REACTIONS '{'
            list_rx_stmts
          '}'
;

list_rx_stmts: rx_stmt
             | list_rx_stmts rx_stmt
;

rx_stmt: rx_group_def
       | rxn
;

rx_group_def:
          REACTION_GROUP reaction_group_name
          '{' list_rxns '}'
;

reaction_group_name: VAR                              { free($1); }
;

list_rxns: rxn
         | list_rxns rxn
;

list_dashes: '-'
           | list_dashes '-'
;

right_arrow:  list_dashes '>';
left_arrow:   '<' list_dashes;
double_arrow: left_arrow '>';

right_cat_arrow:
          list_dashes existing_molecule_opt_orient
          right_arrow                                 { $$.catalyst = $2; $$.flags = ARROW_CATALYTIC; }
;

double_cat_arrow:
          left_arrow existing_molecule_opt_orient
          right_arrow                                 { $$.catalyst = $2; $$.flags = ARROW_CATALYTIC | ARROW_BIDIRECTIONAL; }
;


reaction_arrow:
          right_arrow                                 { $$.catalyst.mol_type = NULL; $$.flags = 0; }
        | right_cat_arrow
        | double_arrow                                { $$.catalyst.mol_type = NULL; $$.flags = ARROW_BIDIRECTIONAL; }
        | double_cat_arrow
;

new_rxn_pathname: /* empty */                         { $$ = NULL; }
                | ':' VAR                             { CHECKN($$ = mdl_new_rxn_pathname(mdlpvp, $2)); }
;

rxn:
          reactant_list opt_reactant_surface_class
          reaction_arrow product_list rx_rate_syntax
          new_rxn_pathname                            { CHECKN(mdl_assemble_reaction(mdlpvp, $1.mol_type_head, &$2, &$3, $4.mol_type_head, &$5, $6)); }
;

reactant_list: reactant                               { CHECKN($$.mol_type_head = $$.mol_type_tail = mdl_new_reaction_player(mdlpvp, &$1)); }
             | reactant_list '+' reactant             { $$ = $1; CHECKN($$.mol_type_tail = $$.mol_type_tail->next = mdl_new_reaction_player(mdlpvp, &$3)); }
;

reactant: existing_molecule_or_subunit
;

existing_molecule_or_subunit:
          existing_molecule_opt_orient                { $$ = $1; $$.is_subunit = 0; }
        | '(' existing_molecule_opt_orient ')'        { $$ = $2; $$.is_subunit = 1; }
;

opt_reactant_surface_class:
          /* empty */                                 { $$.mol_type = NULL; }
        | '@' reactant_surface_class                  { $$ = $2; }
;

reactant_surface_class:
          existing_surface_class orientation_class    { $$ = $2; $$.mol_type = $1; }
;

product_list: product                                 { CHECKN($$.mol_type_head = $$.mol_type_tail = mdl_new_reaction_player(mdlpvp, &$1)); }
             | product_list '+' product               { $$ = $1; CHECKN($$.mol_type_tail = $$.mol_type_tail->next = mdl_new_reaction_player(mdlpvp, &$3)); }
;

product: NO_SPECIES                                   { $$.mol_type = NULL; $$.orient_set = 0; }
       | existing_molecule_or_subunit
;

rx_rate_syntax:
          rx_rate1
        | rx_rate2
;

rx_rate1: '[' rx_dir_rate ']'                         {
                                                        if ($2.forward_rate.rate_type == RATE_UNSET)
                                                        {
                                                          mdlerror(mdlpvp, "Invalid reaction rate specification: must specify a forward rate.");
                                                          return 1;
                                                        }

                                                        $$ = $2;
                                                      }
;

rx_rate2: '[' rx_dir_rate ',' rx_dir_rate ']'         {
                                                        if (($2.forward_rate.rate_type  != RATE_UNSET && $4.forward_rate.rate_type  != RATE_UNSET)  ||
                                                            ($2.backward_rate.rate_type != RATE_UNSET && $4.backward_rate.rate_type != RATE_UNSET))
                                                        {
                                                          mdlerror_fmt(mdlpvp, "Error: When two reaction rates are specified, one must be a forward rate, and one must be a reverse rate");
                                                          return 1;
                                                        }

                                                        $$ = $2;
                                                        if ($4.forward_rate.rate_type != RATE_UNSET)
                                                          $$.forward_rate = $4.forward_rate;
                                                        else
                                                          $$.backward_rate = $4.backward_rate;
                                                      }
;

rx_dir_rate:
          atomic_rate                                 { $$.forward_rate = $1; $$.backward_rate.rate_type = RATE_UNSET; CHECK(mdl_valid_rate(mdlpvp, &$1)); }
        | '>' atomic_rate                             { $$.forward_rate = $2; $$.backward_rate.rate_type = RATE_UNSET; CHECK(mdl_valid_rate(mdlpvp, &$2)); }
        | '<' atomic_rate                             { $$.backward_rate = $2; $$.forward_rate.rate_type = RATE_UNSET; CHECK(mdl_valid_rate(mdlpvp, &$2)); }
;

atomic_rate:
          num_expr_only                               { $$.rate_type = RATE_CONSTANT; $$.v.rate_constant = $1; }
        | str_expr_only                               { $$.rate_type = RATE_FILE; $$.v.rate_file = $1; }
        | existing_var_only                           {
                                                        struct sym_table *rate_sym = $1;
                                                        switch (rate_sym->sym_type)
                                                        {
                                                          case DBL:
                                                            $$.rate_type = RATE_CONSTANT;
                                                            $$.v.rate_constant = *(double *) rate_sym->value;
                                                            break;

                                                          case STR:
                                                            $$.rate_type = RATE_FILE;
                                                            if (($$.v.rate_file = mdl_strdup(mdlpvp, (char *) rate_sym->value)) == NULL)
                                                              return 1;
                                                            break;

                                                          default:
                                                            mdlerror(mdlpvp, "Invalid variable used for rates: must be number or filename");
                                                            return 1;
                                                        }
                                                      }
        | COMPLEX_RATE existing_macromolecule VAR     {
                                                        struct species *complex_species = (struct species *) $2->value;
                                                        if (! (complex_species->flags & IS_COMPLEX))
                                                        {
                                                          mdlerror_fmt(mdlpvp,
                                                                       "The molecule '%s' specified in the complex reaction rate is not a complex.",
                                                                       $2->name);
                                                          free($3);
                                                          return 1;
                                                        }

                                                        $$.rate_type = RATE_COMPLEX;
                                                        if (($$.v.rate_complex = macro_lookup_ruleset((struct complex_species *) complex_species, $3)) == NULL)
                                                        {
                                                          mdlerror_fmt(mdlpvp,
                                                                       "The complex '%s' has no rate rule named '%s'.",
                                                                       $2->name, $3);
                                                          free($3);
                                                          return 1;
                                                        }
                                                        free($3);
                                                      }
;

/* =================================================================== */
/* Release pattern definitions */

release_pattern_def:
          DEFINE_RELEASE_PATTERN
          new_release_pattern
          '{'
              list_req_release_pattern_cmds
          '}'                                         { CHECK(mdl_set_release_pattern(mdlpvp, $2, &$4)); }
;

new_release_pattern: VAR                              { CHECKN($$ = mdl_new_release_pattern(mdlpvp, $1)); }
;

existing_release_pattern_xor_rxpn: VAR                { CHECKN($$ = mdl_existing_release_pattern_or_rxn_pathname(mdlpvp, $1)); }
;

list_req_release_pattern_cmds:
          /* empty */                                 { memset(&$$, 0, sizeof($$)); }
        | list_req_release_pattern_cmds
          DELAY '=' num_expr                          { $$ = $1; $$.delay = $4 / mdlpvp->vol->time_unit; }
        | list_req_release_pattern_cmds
          RELEASE_INTERVAL '=' num_expr               { $$ = $1; $$.release_interval = $4 / mdlpvp->vol->time_unit; }
        | list_req_release_pattern_cmds
          TRAIN_INTERVAL '=' num_expr                 { $$ = $1; $$.train_interval = $4 / volp->time_unit; }
        | list_req_release_pattern_cmds
          TRAIN_DURATION '=' num_expr                 { $$ = $1; $$.train_duration = $4 / volp->time_unit; }
        | list_req_release_pattern_cmds
          NUMBER_OF_TRAINS '=' train_count            { $$ = $1; $$.number_of_trains = $4; }
;

train_count: num_expr                                 { $$ = (int) $1; }
           | UNLIMITED                                { $$ = INT_MAX; }
;

/* =================================================================== */
/* Instance definitions */

instance_def:
          INSTANTIATE                                 { mdlpvp->current_object = mdlpvp->vol->root_instance; }
          new_object OBJECT
          start_object
            list_objects                              {
                                                        struct object *objp = (struct object *) $3->value;
                                                        mdl_add_child_objects(volp->root_instance, objp, objp);
                                                      }
            list_opt_object_cmds
          end_object                                  {
                                                        struct object *new_object = (struct object *) $3->value;
                                                        mdl_add_child_objects(new_object, $6.obj_head, $6.obj_tail);
                                                        new_object->object_type = META_OBJ;
                                                        mdlpvp->current_object = mdlpvp->vol->root_object;
                                                      }
;

/* =================================================================== */
/* Object type definitions */

physical_object_def: object_def                       { mdl_add_child_objects(mdlpvp->vol->root_object, $1, $1); }
;

object_def: meta_object_def
          | release_site_def
          | box_def
          | polygon_list_def
          | voxel_list_def
;

/* This non-terminal has dangerous side effects, modifying the global parser
 * state.  Notably, it updates current_object and adds a name to the object
 * name-list.  Each occurrence of this non-terminal in the grammar MUST have a
 * corresponding end_object (or must explicitly call mdl_finish_object(mdlpvp))
 * when the object scope is closed.
 */
new_object: VAR                                       { CHECKN($$ = mdl_start_object(mdlpvp, $1)); }
;

start_object: '{'
;

end_object: '}'                                       { mdl_finish_object(mdlpvp); }
;

list_opt_object_cmds:
          /* empty */
        | list_opt_object_cmds opt_object_cmd
;


opt_object_cmd: transformation
;

transformation:
          TRANSLATE '=' point                         { mdl_transform_translate(mdlpvp, mdlpvp->current_object->t_matrix, $3); }
        | SCALE '=' point_or_num                      { mdl_transform_scale(mdlpvp, mdlpvp->current_object->t_matrix, $3); }
        | ROTATE '=' point ',' num_expr               { CHECK(mdl_transform_rotate(mdlpvp, mdlpvp->current_object->t_matrix, $3, $5)); }
;

/* Object type: Meta-objects */
meta_object_def:
        new_object OBJECT
        start_object
          list_objects                                { mdlpvp->current_object = (struct object *) $1->value; }
          list_opt_object_cmds
        end_object                                    {
                                                        struct object *new_object = (struct object *) $1->value;
                                                        new_object->object_type=META_OBJ;
                                                        mdl_add_child_objects(new_object, $4.obj_head, $4.obj_tail);
                                                        $$ = new_object;
                                                      }
;

list_objects:
        object_ref                                    { $1->next = NULL; $$.obj_tail = $$.obj_head = $1; }
      | list_objects object_ref                       { $2->next = NULL; $$ = $1; $$.obj_tail = $$.obj_tail->next = $2; }
;

object_ref: existing_object_ref
          | object_def
;

existing_object_ref:
        new_object OBJECT existing_object
        start_object                                  { CHECK(mdl_deep_copy_object(mdlpvp, (struct object *) $1->value, (struct object *) $3->value)); }
          list_opt_object_cmds
        end_object                                    { $$ = (struct object *) $1->value; }
;

/* Object type: Release sites */
release_site_def: release_site_def_old
                | release_site_def_new
;

release_site_def_new:
          new_object RELEASE_SITE
          start_object                                {
                                                        struct object *objp;
                                                        objp = (struct object *) $1->value;
                                                        objp->object_type = REL_SITE_OBJ;
                                                        CHECKN(objp->contents = mdlpvp->current_release_site = mdl_new_release_site(mdlpvp, $1->name));
                                                      }
            release_site_geom
            list_release_site_cmds
            list_opt_object_cmds
          end_object                                  {
                                                        struct object *objp_new = (struct object *) $1->value;
                                                        no_printf("Release site %s defined:\n",$1->name);
                                                        CHECK(mdl_is_release_site_valid(mdlpvp, mdlpvp->current_release_site));
                                                        mdlpvp->current_release_site = NULL;
                                                        $$ = objp_new;
                                                      }
;

release_site_geom: SHAPE '=' release_region_expr      { CHECK(mdl_set_release_site_geometry_region(mdlpvp, mdlpvp->current_release_site, mdlpvp->current_object, $3)); }
                 | SHAPE '=' existing_object          { CHECK(mdl_set_release_site_geometry_object(mdlpvp, mdlpvp->current_release_site, (struct object *) $3->value)); }
                 | SHAPE '=' SPHERICAL                { mdlpvp->current_release_site->release_shape = SHAPE_SPHERICAL; }
                 | SHAPE '=' CUBIC                    { mdlpvp->current_release_site->release_shape = SHAPE_CUBIC; }
                 | SHAPE '=' ELLIPTIC                 { mdlpvp->current_release_site->release_shape = SHAPE_ELLIPTIC; }
                 | SHAPE '=' RECTANGULAR_TOKEN        { mdlpvp->current_release_site->release_shape = SHAPE_RECTANGULAR; }
                 | SHAPE '=' SPHERICAL_SHELL          { mdlpvp->current_release_site->release_shape = SHAPE_SPHERICAL_SHELL; }
                 | SHAPE '=' LIST                     {
                                                        mdlpvp->current_release_site->release_shape = SHAPE_LIST;
                                                        mdlpvp->current_release_site->release_number_method = CONSTNUM;
                                                      }
;

release_region_expr:
       existing_region                                { CHECKN($$ = mdl_new_release_region_expr_term(mdlpvp, $1)); }
     | '(' release_region_expr ')'                    { $$ = $2; }
     | release_region_expr '+' release_region_expr    { CHECKN($$ = mdl_new_release_region_expr_binary(mdlpvp, $1, $3, REXP_UNION)); }
     | release_region_expr '-' release_region_expr    { CHECKN($$ = mdl_new_release_region_expr_binary(mdlpvp, $1, $3, REXP_SUBTRACTION)); }
     | release_region_expr '*' release_region_expr    { CHECKN($$ = mdl_new_release_region_expr_binary(mdlpvp, $1, $3, REXP_INTERSECTION)); }
     | release_region_expr '&' release_region_expr    { CHECKN($$ = mdl_new_release_region_expr_binary(mdlpvp, $1, $3, REXP_INCLUSION)); }
;

release_site_def_old:
          new_object release_site_geom_old
          start_object                                {
                                                        struct object *objp=(struct object *) $1->value;
                                                        objp->object_type = REL_SITE_OBJ;
                                                        CHECKN(objp->contents = mdlpvp->current_release_site = mdl_new_release_site(mdlpvp, $1->name));
                                                        mdlpvp->current_release_site->release_shape = $2;
                                                      }
            list_release_site_cmds
            list_opt_object_cmds
          end_object                                  {
                                                        struct object *objp_new = (struct object *) $1->value;
                                                        no_printf("Release site %s defined:\n", $1->name);
                                                        no_printf("\tLocation = [%f,%f,%f]\n",
                                                                  mdlpvp->current_release_site->location->x,
                                                                  mdlpvp->current_release_site->location->y,
                                                                  mdlpvp->current_release_site->location->z);
                                                        mdlpvp->current_release_site = NULL;
                                                        $$ = objp_new;
                                                      }
;

release_site_geom_old: SPHERICAL_RELEASE_SITE         { $$ = SHAPE_SPHERICAL; }
                     | CUBIC_RELEASE_SITE             { $$ = SHAPE_CUBIC; }
                     | ELLIPTIC_RELEASE_SITE          { $$ = SHAPE_ELLIPTIC; }
                     | RECTANGULAR_RELEASE_SITE       { $$ = SHAPE_RECTANGULAR; }
                     | SPHERICAL_SHELL_SITE           { $$ = SHAPE_SPHERICAL_SHELL; }
;

list_release_site_cmds:
          release_site_cmd
        | list_release_site_cmds release_site_cmd
;

existing_num_or_array: VAR                            { CHECKN($$ = mdl_existing_num_or_array(mdlpvp, $1)); }
;

release_site_cmd:
          LOCATION '=' point                          {
                                                        mdlpvp->current_release_site->location = $3;
                                                        mdlpvp->current_release_site->location->x *= mdlpvp->vol->r_length_unit;
                                                        mdlpvp->current_release_site->location->y *= mdlpvp->vol->r_length_unit;
                                                        mdlpvp->current_release_site->location->z *= mdlpvp->vol->r_length_unit;
                                                      }
        | MOLECULE '=' existing_molecule_opt_orient   { CHECK(mdl_set_release_site_molecule(mdlpvp, mdlpvp->current_release_site, & $3)); }
        | release_number_cmd                          {
                                                        if (mdlpvp->current_release_site->release_shape == SHAPE_LIST)
                                                        {
                                                          mdlerror(mdlpvp, "Molecules are already specified in a list--cannot set number or density.");
                                                          return 1;
                                                        }
                                                      }
        | site_size_cmd '=' num_expr_only             { CHECK(mdl_set_release_site_diameter(mdlpvp, mdlpvp->current_release_site, $3 * (($1 == SITE_RADIUS) ? 2.0 : 1.0))); }
        | site_size_cmd '=' array_expr_only           { CHECK(mdl_set_release_site_diameter_array(mdlpvp, mdlpvp->current_release_site, $3.value_count, $3.value_head, ($1 == SITE_RADIUS) ? 2.0 : 1.0)); }
        | site_size_cmd '=' existing_num_or_array     {
                                                        double scaling_factor = ($1 == SITE_RADIUS) ? 2.0 : 1.0;
                                                        struct num_expr_list *elp;
                                                        int count = 0;
                                                        if ((mdlpvp->current_release_site->diameter=(struct vector3 *)malloc(sizeof(struct vector3)))==NULL)
                                                        {
                                                          mdlerror(mdlpvp, "Out of memory while storing release site diameter");
                                                          return 1;
                                                        }
                                                        switch ($3->sym_type) {
                                                          case DBL:
                                                            CHECK(mdl_set_release_site_diameter(mdlpvp, mdlpvp->current_release_site, *(double *) $3->value * scaling_factor));
                                                            break;
                                                          case ARRAY:
                                                            for (elp = (struct num_expr_list *) $3->value; elp != NULL && count < 4; ++ count, elp = elp->next)
                                                              ;
                                                            CHECK(mdl_set_release_site_diameter_array(mdlpvp, mdlpvp->current_release_site, count, elp, scaling_factor));
                                                            break;
                                                          default:
                                                            mdlerror(mdlpvp, "Diameter must either be a number or a 3-valued vector.");
                                                            return 1;
                                                        }
                                                      }
        | RELEASE_PROBABILITY '=' num_expr            { CHECK(mdl_set_release_site_probability(mdlpvp, mdlpvp->current_release_site, $3)); }
        | RELEASE_PATTERN '='
          existing_release_pattern_xor_rxpn           { CHECK(mdl_set_release_site_pattern(mdlpvp, mdlpvp->current_release_site, $3)); }
        | MOLECULE_POSITIONS
          '{' molecule_release_pos_list '}'           {
                                                        if (mdlpvp->current_release_site->release_shape != SHAPE_LIST)
                                                        {
                                                          mdlerror(mdlpvp, "You must use the LIST shape to specify molecule positions in a release.");
                                                          return 1;
                                                        }

                                                        struct release_single_molecule *rsm;
                                                        if (mdlpvp->current_release_site->mol_list == NULL)
                                                          mdlpvp->current_release_site->mol_list = $3.rsm_head;
                                                        else
                                                        {
                                                          for (rsm = mdlpvp->current_release_site->mol_list; rsm->next != NULL; rsm = rsm->next)
                                                            ;
                                                          rsm->next = $3.rsm_head;
                                                        }
                                                        mdlpvp->current_release_site->release_number += $3.rsm_count;
                                                      }
;

site_size_cmd:
          SITE_DIAMETER                               { $$ = SITE_DIAMETER; }
        | SITE_RADIUS                                 { $$ = SITE_RADIUS; }
;

release_number_cmd:
          constant_release_number_cmd
        | gaussian_release_number_cmd
        | volume_dependent_number_cmd
        | concentration_dependent_release_cmd
;


constant_release_number_cmd:
          NUMBER_TO_RELEASE '=' num_expr              { mdlpvp->current_release_site->release_number_method=CONSTNUM; mdlpvp->current_release_site->release_number= $3; }
        | GAUSSIAN_RELEASE_NUMBER '{'
          MEAN_NUMBER '=' num_expr
          '}'                                         { mdlpvp->current_release_site->release_number_method=CONSTNUM; mdlpvp->current_release_site->release_number= $5; }
;

gaussian_release_number_cmd:
          GAUSSIAN_RELEASE_NUMBER '{'
            MEAN_NUMBER '=' num_expr
            STANDARD_DEVIATION '=' num_expr
          '}'                                         {
                                                        mdlpvp->current_release_site->release_number_method=GAUSSNUM;
                                                        mdlpvp->current_release_site->release_number= $5;
                                                        mdlpvp->current_release_site->standard_deviation=$8;
                                                      }
;

volume_dependent_number_cmd:
          VOLUME_DEPENDENT_RELEASE_NUMBER '{'
            MEAN_DIAMETER '=' num_expr
            STANDARD_DEVIATION '=' num_expr
            CONCENTRATION '=' num_expr
            '}'                                       {
                                                        mdlpvp->current_release_site->release_number_method = VOLNUM;
                                                        mdlpvp->current_release_site->mean_diameter = $5;
                                                        mdlpvp->current_release_site->standard_deviation = $8;
                                                        mdlpvp->current_release_site->concentration = $11;
                                                      }
;

concentration_dependent_release_cmd:
          CONCENTRATION '=' num_expr                  {
                                                        struct release_site_obj *rsop = mdlpvp->current_release_site;
                                                        if (rsop->release_shape == SHAPE_SPHERICAL_SHELL)
                                                        {
                                                          mdlerror_fmt(mdlpvp,
                                                                       "Release site '%s' is a spherical shell; concentration-based release is not supported on a spherical shell\n",
                                                                       rsop->name);
                                                          return 1;
                                                        }
                                                        rsop->release_number_method = CCNNUM;
                                                        rsop->release_number = -3; /* Expect 3D molecule */
                                                        rsop->concentration = $3;
                                                      }
        | DENSITY '=' num_expr                        {
                                                        mdlpvp->current_release_site->release_number_method=CCNNUM;
                                                        mdlpvp->current_release_site->release_number=-2; /* Expect 2D molecule */
                                                        mdlpvp->current_release_site->concentration = $3;
                                                      }
;

molecule_release_pos_list:
          molecule_release_pos                        { $$.rsm_tail = $$.rsm_head = $1; $$.rsm_count = 1; }
        | molecule_release_pos_list
          molecule_release_pos                        { $$ = $1; $$.rsm_tail = $$.rsm_tail->next = $2; ++ $$.rsm_count; }
;

molecule_release_pos:
          existing_molecule_opt_orient point          { CHECKN($$ = mdl_new_release_single_molecule(mdlpvp, &$1, $2)); }
;

/* Object type: Polygons */
polygon_list_def:
          new_object POLYGON_LIST
          start_object
            vertex_list_cmd
            element_connection_cmd                    {
                                                        CHECKN(mdlpvp->current_polygon = mdl_new_polygon_list(mdlpvp, $1,
                                                                                                 $4.vertex_count, $4.vertex_head,
                                                                                                 $5.connection_count, $5.connection_head));
                                                      }
            list_opt_polygon_object_cmds
            list_opt_object_cmds
          end_object                                  {
                                                        struct object *objp = (struct object *) $1->value;
                                                        mdl_remove_gaps_from_regions(objp);
                                                        no_printf("Polygon list %s defined:\n",$1->name);
                                                        no_printf(" n_verts = %d\n",mdlpvp->current_polygon->n_verts);
                                                        no_printf(" n_walls = %d\n",mdlpvp->current_polygon->n_walls);
                                                        CHECK(mdl_check_degenerate_polygon_list(mdlpvp, objp));
                                                        mdlpvp->current_polygon = NULL;
                                                        $$ = objp;
                                                      }
;

vertex_list_cmd: VERTEX_LIST '{' list_points '}'      { $$ = $3; }
;

single_vertex: point                                  {
                                                        struct vertex_list *vlp;
                                                        if ((vlp=(struct vertex_list *)malloc(sizeof(struct vertex_list)))==NULL)
                                                        {
                                                          mdlerror(mdlpvp, "Out of memory while creating vertices");
                                                          return 1;
                                                        }
                                                        vlp->vertex = $1;
                                                        vlp->normal = NULL;
                                                        vlp->next = NULL;
                                                        $$ = vlp;
                                                      }
             | point NORMAL point                     {
                                                        struct vertex_list *vlp;
                                                        if ((vlp=(struct vertex_list *)malloc(sizeof(struct vertex_list)))==NULL)
                                                        {
                                                          mdlerror(mdlpvp, "Out of memory while creating normals");
                                                          return 1;
                                                        }
                                                        vlp->vertex = $1;
                                                        vlp->normal = $3;
                                                        vlp->next = NULL;
                                                        $$ = vlp;
                                                      }
;

list_points: single_vertex                            {
                                                        $$.vertex_tail = $$.vertex_head = $1;
                                                        $$.vertex_count = 1;
                                                      }
           | list_points single_vertex                {
                                                        $$ = $1;
                                                        $$.vertex_tail = $$.vertex_tail->next = $2;
                                                        ++ $$.vertex_count;
                                                      }
;

element_connection_cmd:
          ELEMENT_CONNECTIONS
          '{' list_element_connections '}'            { $$ = $3; }
;

list_element_connections:
          element_connection                          {
                                                        $$.connection_head = $$.connection_tail = $1;
                                                        $$.connection_count = 1;
                                                      }
        | list_element_connections
          element_connection                          {
                                                        $$ = $1;
                                                        $$.connection_tail = $$.connection_tail->next = $2;
                                                        ++ $$.connection_count;
                                                      }
;

element_connection: array_value                       {
                                                        struct element_connection_list *eclp = (struct element_connection_list *) malloc(sizeof(struct element_connection_list));
                                                        if (eclp == NULL)
                                                        {
                                                          mdlerror(mdlpvp, "Out of memory while creating element connections");
                                                          return 1;
                                                        }

                                                        eclp->connection_list = $1.value_head;
                                                        eclp->n_verts = $1.value_count;
                                                        eclp->next = NULL;

                                                        if (eclp->n_verts != 3)
                                                        {
                                                          mdlerror(mdlpvp, "Non-triangular element found in polygon list object");
                                                          return 1;
                                                        }
                                                        $$ = eclp;
                                                      }
;

list_opt_polygon_object_cmds:
          /* empty */
        | list_opt_polygon_object_cmds
          opt_polygon_object_cmd
;


opt_polygon_object_cmd:
          remove_side
        | in_obj_define_surface_regions
;

remove_side:
          REMOVE_ELEMENTS '{'                         { CHECKN(mdlpvp->current_region = mdl_get_region(mdlpvp, mdlpvp->current_object, "REMOVED")); }
            remove_element_specifier_list
          '}'                                         {
                                                        mdlpvp->current_region->element_list_head = $4.elml_head;
                                                        if (mdlpvp->current_object->object_type == POLY_OBJ)
                                                        {
                                                          if (mdl_normalize_elements(mdlpvp, mdlpvp->current_region,0))
                                                          {
                                                            mdlerror(mdlpvp, "Improper element specification: out of range or out of memory");
                                                            return 1;
                                                          }
                                                        }
                                                      }
;

remove_element_specifier_list:
          element_specifier_list
        | just_an_element_list
;

side_name: TOP                                        { $$ = Z_POS; }
         | BOTTOM                                     { $$ = Z_NEG; }
         | FRONT                                      { $$ = Y_NEG; }
         | BACK                                       { $$ = Y_POS; }
         | LEFT                                       { $$ = X_NEG; }
         | RIGHT                                      { $$ = X_POS; }
         | ALL_ELEMENTS                               { $$ = ALL_SIDES; }
;

element_specifier_list:
          element_specifier
        | element_specifier_list ','
          element_specifier                           { $$ = $1; $$.elml_tail->next = $3.elml_head; $$.elml_tail = $3.elml_tail; }
;

element_specifier:
          incl_element_list_stmt
        | excl_element_list_stmt
        | prev_region_stmt                            { $$.elml_tail = $$.elml_head = $1; }
        | patch_statement                             { $$.elml_tail = $$.elml_head = $1; }
;

incl_element_list_stmt:
          INCLUDE_ELEMENTS '='
          '[' list_element_specs ']'                  { $$ = $4; }
;

excl_element_list_stmt:
          EXCLUDE_ELEMENTS '='
          '[' list_element_specs ']'                  {
                                                        struct element_list *elmlp;
                                                        $$ = $4;
                                                        for (elmlp = $$.elml_head; elmlp != NULL; elmlp = elmlp->next)
                                                          elmlp->special = (struct element_special *) elmlp;
                                                      }
;

just_an_element_list: list_element_specs
;

list_element_specs:
          element_spec                                { $$.elml_tail = $$.elml_head = $1; }
        | list_element_specs ',' element_spec         { $$ = $1; $$.elml_tail = $$.elml_tail->next = $3; }
;

element_spec: num_expr                                { CHECKN($$ = mdl_new_element_list(mdlpvp, (unsigned int) $1, (unsigned int) $1)); }
            | num_expr TO num_expr                    { CHECKN($$ = mdl_new_element_list(mdlpvp, (unsigned int) $1, (unsigned int) $3)); }
            | side_name                               {
                                                        unsigned int begin, end;
                                                        if ($1 == ALL_SIDES  &&  mdlpvp->current_object->object_type == POLY_OBJ)
                                                        {
                                                          begin = 0;
                                                          end = mdlpvp->current_polygon->n_walls-1;
                                                        }
                                                        else if (mdlpvp->current_object->object_type==POLY_OBJ)
                                                        {
                                                          mdlerror(mdlpvp, "Illegal reference to polygon list element by side-name");
                                                          return 1;
                                                        }
                                                        else {
                                                          begin = $1;
                                                          end = $1;
                                                        }
                                                        CHECKN($$ = mdl_new_element_list(mdlpvp, begin, end));
                                                      }
;

prev_region_stmt: prev_region_type '=' VAR            { CHECKN($$ = mdl_new_element_previous_region(mdlpvp, mdlpvp->current_object, mdlpvp->current_region, $3, $1)); }
;

prev_region_type: INCLUDE_REGION                      { $$ = 0; }
                | EXCLUDE_REGION                      { $$ = 1; }
;

patch_statement: patch_type '=' point ',' point       {
                                                        if (mdlpvp->current_object->object_type!=BOX_OBJ)
                                                        {
                                                          mdlerror(mdlpvp, "Must use PATCH specifier on a BOX object only.");
                                                          return 1;
                                                        }

                                                        CHECKN($$ = mdl_new_element_patch(mdlpvp, mdlpvp->current_polygon->sb, $3, $5, $1));
                                                      }
;

patch_type: INCLUDE_PATCH                             { $$ = 0; }
          | EXCLUDE_PATCH                             { $$ = 1; }
;

in_obj_define_surface_regions:
          DEFINE_SURFACE_REGIONS '{'
            list_in_obj_surface_region_defs
          '}'
;

list_in_obj_surface_region_defs:
          in_obj_surface_region_def
        | list_in_obj_surface_region_defs
          in_obj_surface_region_def
;

in_obj_surface_region_def:
          new_region '{'                              { mdlpvp->current_region = $1; }
            element_specifier_list                    {
                                                        $1->element_list_head = $4.elml_head;
                                                        if (mdlpvp->current_object->object_type == POLY_OBJ)
                                                        {
                                                          if (mdl_normalize_elements(mdlpvp, mdlpvp->current_region,0))
                                                          {
                                                            mdlerror(mdlpvp, "Improper element specification: out of range or out of memory");
                                                            return 1;
                                                          }
                                                        }
                                                      }
            list_opt_surface_region_stmts
          '}'                                         {
                                                        if (mdlpvp->current_region->surf_class != NULL  &&
                                                            mdlpvp->current_region->surf_class->region_viz_value > 0)
                                                          mdlpvp->current_region->region_viz_value = mdlpvp->current_region->surf_class->region_viz_value;
                                                          mdlpvp->current_region = NULL;
                                                      }
;

/* Object type: Voxel list */
voxel_list_def:
          new_object VOXEL_LIST
          start_object
            vertex_list_cmd
            tet_element_connection_cmd                {
                                                        struct voxel_object *vop;
                                                        CHECKN(vop = mdl_new_voxel_list(mdlpvp, $1,
                                                                                        $4.vertex_count, $4.vertex_head,
                                                                                        $5.connection_count, $5.connection_head));
                                                        no_printf("Voxel list %s defined:\n",$1->name);
                                                        no_printf(" n_verts = %d\n",vop->n_verts);
                                                      }
            list_opt_object_cmds
          end_object                                  { $$ = (struct object *) $1->value; }
;

tet_element_connection_cmd:
          TET_ELEMENT_CONNECTIONS
          '{' list_tet_arrays '}'                     { $$ = $3; }
;

element_connection_tet: array_value                   {
                                                        struct element_connection_list *eclp = (struct element_connection_list *) malloc(sizeof(struct element_connection_list));
                                                        if (eclp == NULL)
                                                        {
                                                          mdlerror(mdlpvp, "Out of memory while creating element connections");
                                                          return 1;
                                                        }
                                                        eclp->connection_list = $1.value_head;
                                                        eclp->n_verts = $1.value_count;
                                                        eclp->next = NULL;
                                                        if (eclp->n_verts != 4)
                                                        {
                                                          mdlerror(mdlpvp, "Non-tetrahedron element found in voxel list object");
                                                          return 1;
                                                        }
                                                        $$ = eclp;
                                                      }
;

list_tet_arrays:
          element_connection_tet                      {
                                                        $$.connection_head = $$.connection_tail = $1;
                                                        $$.connection_count = 1;
                                                      }
        | list_tet_arrays element_connection_tet      {
                                                        $$ = $1;
                                                        $$.connection_tail = $$.connection_tail->next = $2;
                                                        ++ $$.connection_count;
                                                      }
;

/* Object type: Boxes */
box_def: new_object BOX
          start_object
            CORNERS '=' point ',' point
            opt_aspect_ratio_def                      { CHECKN(mdlpvp->current_polygon = mdl_new_box_object(mdlpvp, $1, $6, $8)); }
            list_opt_polygon_object_cmds              { CHECK(mdl_triangulate_box_object(mdlpvp, $1, mdlpvp->current_polygon, $9)); }
            list_opt_object_cmds
          end_object                                  {
                                                        struct object *objp = (struct object *) $1->value;
                                                        mdl_remove_gaps_from_regions(objp);
                                                        objp->n_walls = mdlpvp->current_polygon->n_walls;
                                                        objp->n_verts = mdlpvp->current_polygon->n_verts;
                                                        CHECK(mdl_check_degenerate_polygon_list(mdlpvp, objp));
                                                        mdlpvp->current_polygon = NULL;
                                                        $$ = objp;
                                                      }
;

opt_aspect_ratio_def: /* empty */                     { $$ = 0.0; }
                    | ASPECT_RATIO '=' num_expr       {
                                                        $$ = $3;
                                                        if (!($$ >= 2.0))
                                                        {
                                                          mdlerror(mdlpvp, "Invalid aspect ratio selected (must be greater than or equal to 2.0)");
                                                          return 1;
                                                        }
                                                      }
;

/* =================================================================== */
/* Surface region definitions */

existing_obj_define_surface_regions:
    DEFINE_SURFACE_REGIONS '{'
      list_existing_obj_surface_region_defs
    '}'
;

list_existing_obj_surface_region_defs:
      existing_obj_surface_region_def
    | list_existing_obj_surface_region_defs
      existing_obj_surface_region_def
;

existing_obj_surface_region_def:
          existing_object                             {
                                                        struct object *objp = (struct object *) $1->value;
                                                        if (objp->object_type != BOX_OBJ && objp->object_type != POLY_OBJ)
                                                        {
                                                          mdlerror_fmt(mdlpvp, "Cannot define region on non-surface object: %s", $1->name);
                                                          return 1;
                                                        }
                                                        mdlpvp->current_polygon = objp->contents;
                                                        mdlpvp->current_object = objp;
                                                      }
            '[' new_region ']'                        { mdlpvp->current_region = $4; }
            '{' element_specifier_list                {
                                                        mdlpvp->current_region->element_list_head = $8.elml_head;
                                                        if (mdl_normalize_elements(mdlpvp, mdlpvp->current_region, 1))
                                                        {
                                                          mdlerror_fmt(mdlpvp, "Improper element specification: out of range, out of memory,\n  or using patch specifier in already-created box");
                                                          return 1;
                                                        }
                                                      }
            list_opt_surface_region_stmts
          '}'                                         {
                                                        mdlpvp->current_region = NULL;
                                                        mdlpvp->current_polygon = NULL;
                                                        mdlpvp->current_object = NULL;
                                                      }
;

new_region: VAR                                       { CHECKN($$ = mdl_create_region(mdlpvp, mdlpvp->current_object, $1)); }
;

list_opt_surface_region_stmts:
          /* empty */
        | list_opt_surface_region_stmts
          opt_surface_region_stmt
;

opt_surface_region_stmt:
          set_surface_class_stmt
        | surface_mol_stmt                            { $1.eff_tail->next = mdlpvp->current_region->eff_dat_head; mdlpvp->current_region->eff_dat_head = $1.eff_head; }
        | surface_region_viz_value_stmt
;

set_surface_class_stmt:
          SURFACE_CLASS '='
          existing_surface_class                      {
                                                        mdlpvp->current_region->surf_class = (struct species *) $3->value;
                                                        if (mdlpvp->current_region->surf_class->region_viz_value > 0)
                                                           mdlpvp->current_region->region_viz_value = mdlpvp->current_region->surf_class->region_viz_value;
                                                      }
;

surface_region_viz_value_stmt:
          VIZ_VALUE '=' num_expr                      {
                                                        /* if surface_class->region_viz_value is already defined print warning */
                                                        if (mdlpvp->current_region->surf_class != NULL  &&
                                                            mdlpvp->current_region->surf_class->region_viz_value > 0)
                                                          mdlerror(mdlpvp, "ATTENTION: region_viz_value defined both through SURFACE_CLASS and VIZ_VALUE statements.\n");

                                                        mdlpvp->current_region->region_viz_value = (int) $3;
                                                      }
;

/* =================================================================== */
/* Surface region modifications */

mod_surface_regions:
          MODIFY_SURFACE_REGIONS '{'
            list_existing_surface_region_refs
          '}'
;

list_existing_surface_region_refs:
          existing_surface_region_ref
        | list_existing_surface_region_refs
          existing_surface_region_ref
;

existing_surface_region_ref:
          existing_region '{'                         { mdlpvp->current_region = (struct region *) $1->value; }
            list_opt_surface_region_stmts
          '}'                                         { mdlpvp->current_region = NULL; }
;

/* =================================================================== */
/* Reaction output definitions */

output_def:
          REACTION_DATA_OUTPUT '{'
            output_buffer_size_def                    {
                                                        mdlpvp->header_comment = NULL;  /* No header by default */
                                                        mdlpvp->exact_time_flag = 1;    /* Print exact_time column in TRIGGER output by default */
                                                      }
            output_timer_def
            list_count_cmds
          '}'                                         {
                                                        struct output_block *obp;
                                                        struct output_set *os;
                                                        CHECKN(obp = mdl_new_output_block(mdlpvp, (int) $3));
                                                        if ($5.type == OUTPUT_BY_STEP)
                                                          mdl_set_reaction_output_timer_step(mdlpvp, obp, $5.step);
                                                        else if ($5.type == OUTPUT_BY_ITERATION_LIST)
                                                          mdl_set_reaction_output_timer_iterations(mdlpvp, obp, $5.values.value_count, $5.values.value_head);
                                                        else if ($5.type == OUTPUT_BY_TIME_LIST)
                                                          mdl_set_reaction_output_timer_times(mdlpvp, obp, $5.values.value_count, $5.values.value_head);
                                                        else
                                                        {
                                                          mdlerror_fmt(mdlpvp, "Internal error: Invalid output timer def (%d)\n", $5.type);
                                                          return 1;
                                                        }
                                                        obp->data_set_head = $6.set_head;
                                                        for (os = obp->data_set_head; os != NULL; os = os->next)
                                                          os->block = obp;
                                                        CHECK(mdl_output_block_finalize(mdlpvp, obp));
                                                        obp->next = mdlpvp->vol->output_block_head;
                                                        mdlpvp->vol->output_block_head = obp;
                                                      }
;

output_buffer_size_def:
          /* empty */                                 { $$ = COUNTBUFFERSIZE; }
        | OUTPUT_BUFFER_SIZE '=' num_expr             {
                                                        double temp_value = $3;
                                                        if (!(temp_value >= 1.0 && temp_value < UINT_MAX))
                                                        {
                                                          mdlerror(mdlpvp, "Buffer size invalid.  Suggested range is 100-1000000.\n");
                                                          return 1;
                                                        }
                                                        $$ = $3;
                                                      }
;

output_timer_def: step_time_def
                | iteration_time_def
                | real_time_def
;

step_time_def: STEP '=' num_expr                      { $$.type = OUTPUT_BY_STEP; $$.step = $3; }
             | /* empty */                            { $$.type = OUTPUT_BY_STEP; $$.step = 1.0; }
;

iteration_time_def:
          ITERATION_LIST '=' '[' list_range_specs ']' {
                                                        $$.type = OUTPUT_BY_ITERATION_LIST;
                                                        $$.values = $4;
                                                      }
;


real_time_def:
          TIME_LIST '=' '[' list_range_specs ']'      {
                                                        $$.type = OUTPUT_BY_TIME_LIST;
                                                        $$.values = $4;
                                                      }
;

list_count_cmds:
          count_cmd                                   { $$.set_head = $$.set_tail = $1; }
        | list_count_cmds
          count_cmd                                   {
                                                        $$ = $1;
                                                        if ($2 != NULL)
                                                        {
                                                          if ($$.set_tail != NULL)
                                                            $$.set_tail = $$.set_tail->next = $2;
                                                        }
                                                      }
;

count_cmd:
          count_stmt
        | custom_header                               { $$ = NULL; }
        | exact_time_toggle                           { $$ = NULL; }
;

count_stmt:
          '{'                                         {
                                                        struct output_set *os;
                                                        CHECKN(os = mdl_new_output_set(mdlpvp, mdlpvp->header_comment));
                                                        os->exact_time_flag = mdlpvp->exact_time_flag;
                                                        mdlpvp->count_flags = 0;
                                                        $<ro_set>$ = os;
                                                      }
            list_count_exprs
          '}' file_arrow outfile_syntax               {
                                                        struct output_set *os = $<ro_set>2;
                                                        struct output_column *oc;
                                                        os->column_head = $3.column_head;
                                                        os->file_flags = $5;
                                                        os->outfile_name = $6;
                                                        no_printf("Counter output file set to %s\n", os->outfile_name);

                                                        if ((mdlpvp->count_flags & (TRIGGER_PRESENT|COUNT_PRESENT)) == (TRIGGER_PRESENT|COUNT_PRESENT))
                                                        {
                                                          mdlerror(mdlpvp, "Cannot mix TRIGGER and COUNT statements.  Use separate files.");
                                                          return 1;
                                                        }

                                                        for (oc = os->column_head; oc != NULL; oc = oc->next)
                                                          oc->set = os;

                                                        CHECK(mdl_check_reaction_output_file(mdlpvp, os));
                                                        $$ = os;
                                                      }
;

custom_header_value:
          NONE                                        { $$ = NULL; }
        | boolean                                     { $$ = ($1 ? "" : NULL); }
        | str_expr                                    { $$ = $1; }
;

custom_header:
          HEADER '=' custom_header_value              { mdlpvp->header_comment = $3; }
;

exact_time_toggle:
          SHOW_EXACT_TIME '=' boolean                 { mdlpvp->exact_time_flag = $3; }
;

list_count_exprs:
          single_count_expr
        | list_count_exprs ','
          single_count_expr                           {
                                                        $$ = $1;
                                                        $$.column_tail->next = $3.column_head;
                                                        $$.column_tail = $3.column_tail;
                                                      }
;

single_count_expr:
        count_expr opt_custom_header                  {
                                                        struct output_expression *oer,*oe;
                                                        struct output_column *oc;

                                                        $$.column_head = NULL;
                                                        $$.column_tail = NULL;

                                                        oer = $1; /* Root of count expression */

                                                        if (oer->oper == ',' && $2 != NULL)
                                                        {
                                                          mdlerror(mdlpvp, "Cannot use custom column headers with wildcard expansion");
                                                          return 1;
                                                        }

                                                        if ($2 != NULL) oer->title=$2;

                                                        /* If we have a list of results, go through to build column stack */
                                                        for (oe=first_oexpr_tree(oer);oe!=NULL;oe=next_oexpr_tree(oe))
                                                        {
                                                          CHECKN(oc = mdl_new_output_column(mdlpvp));
                                                          if (! $$.column_head)
                                                            $$.column_head = $$.column_tail = oc;
                                                          else
                                                            $$.column_tail = $$.column_tail->next = oc;

                                                          oc->expr=oe;
                                                          set_oexpr_column(oe, oc);
                                                        }
                                                      }
;

count_expr:
          num_value                                   { CHECKN($$ = mdl_new_oexpr_constant(mdlpvp, $1)); }
        | count_value
        | '(' count_expr ')'                          { CHECKN($$ = mdl_join_oexpr_tree(mdlpvp, $2, NULL, '(')); }
        | count_expr '+' count_expr                   { CHECKN($$ = mdl_join_oexpr_tree(mdlpvp, $1,   $3, '+')); }
        | count_expr '-' count_expr                   { CHECKN($$ = mdl_join_oexpr_tree(mdlpvp, $1,   $3, '-')); }
        | count_expr '*' count_expr                   { CHECKN($$ = mdl_join_oexpr_tree(mdlpvp, $1,   $3, '*')); }
        | count_expr '/' count_expr                   { CHECKN($$ = mdl_join_oexpr_tree(mdlpvp, $1,   $3, '/')); }
        | '-' count_expr %prec UNARYMINUS             { CHECKN($$ = mdl_join_oexpr_tree(mdlpvp, $2, NULL, '_')); }
        | SUMMATION_OPERATOR '(' count_expr ')'       { CHECKN($$ = mdl_sum_oexpr($3)); }
;


count_value:
          COUNT                                       { mdlpvp->count_flags |= COUNT_PRESENT; }
          '[' count_syntax ']'                        { $$ = $4; }
        | EXPRESSION '[' num_expr ']'                 { CHECKN($$ = mdl_new_oexpr_constant(mdlpvp, $3)); }
        | TRIGGER                                     { mdlpvp->count_flags |= TRIGGER_PRESENT; }
          '[' count_syntax ']'                        { $$ = $4; }
;

file_arrow: '>'                                       { $$ = FILE_OVERWRITE; }
          | '=' '>'                                   { $$ = FILE_SUBSTITUTE; }
          | '>' '>'                                   { $$ = FILE_APPEND; }
          | '>' '>' '>'                               { $$ = FILE_APPEND_HEADER; }
          | '+' '>'                                   { $$ = FILE_CREATE; }
;

outfile_syntax: file_name
;

existing_rxpn_or_molecule: VAR                        { CHECKN($$ = mdl_existing_rxn_pathname_or_molecule(mdlpvp, $1)); }
;

existing_molecule_required_orient_braces:
          VAR orient_class_number                     {
                                                        $$ = $2;
                                                        if ($$.orient > 0)
                                                          $$.orient = 1;
                                                        else if ($$.orient < 0)
                                                          $$.orient = -1;
                                                        CHECKN($$.mol_type = mdl_existing_molecule(mdlpvp, $1));
                                                      }
;

count_syntax: count_syntax_1
            | count_syntax_2
            | count_syntax_3
            | count_syntax_macromol
;

count_syntax_1:
    existing_rxpn_or_molecule ','
    count_location_specifier opt_hit_spec             { CHECKN($$ = mdl_count_syntax_1(mdlpvp, $1, $3, $4, mdlpvp->count_flags)); }
;

count_syntax_2:
    existing_molecule_required_orient_braces ','
    count_location_specifier opt_hit_spec             { CHECKN($$ = mdl_count_syntax_2(mdlpvp, $1.mol_type, $1.orient, $3, $4, mdlpvp->count_flags)); }
;

count_syntax_3:
    STR_VALUE  ','
    count_location_specifier opt_hit_spec             { CHECKN($$ = mdl_count_syntax_3(mdlpvp, $1, $3, $4, mdlpvp->count_flags)); }
;

count_syntax_macromol:
    count_syntax_macromol_subunit
;

count_syntax_macromol_subunit:
          SUBUNIT '{'
            existing_macromolecule                    { mdlpvp->current_complex = (struct complex_species *) $3->value; }
            orientation_class
            ':' existing_molecule_opt_orient
            opt_macromol_relation_states
          '}' ',' count_location_specifier            {
                                                        mdlpvp->current_complex = NULL;
                                                        struct complex_species *macromol = (struct complex_species *) $3->value;
                                                        struct species_opt_orient master_orientation = $5;
                                                        struct species_opt_orient subunit = $7;
                                                        struct macro_relation_state *relation_states = $8;
                                                        struct sym_table *location = $11;
                                                        CHECKN($$ = mdl_count_syntax_macromol_subunit(mdlpvp, macromol, &master_orientation, & subunit, relation_states, location));
                                                      }
;

opt_macromol_relation_states:
          /* empty */                                 { $$ = NULL; }
        | '[' macromol_relation_state_list ']'        { $$ = $2; }
;

macromol_relation_state_list:
          macromol_relation_state
        | macromol_relation_state_list '&'
          macromol_relation_state                     { $3->next = $1; $$ = $3; }
;

macromol_relation_state:
          macromol_relation_name
            equal_or_not
            existing_molecule_opt_orient              { CHECKN($$ = mdl_assemble_complex_relation_state(mdlpvp, $1, $2, & $3)); }
;

macromol_relation_name: VAR                           {
                                                        int rel_idx = macro_lookup_relation(mdlpvp->current_complex, $1);
                                                        if (rel_idx == -1)
                                                        {
                                                          mdlerror_fmt(mdlpvp,
                                                                       "In subunit specification for COUNT statement, relation '%s' does not exist within the complex '%s'",
                                                                       $1,
                                                                       mdlpvp->current_complex->base.sym->name);
                                                          return 1;
                                                        }

                                                        $$ = rel_idx;
                                                      }
;

count_location_specifier: WORLD                       { $$ = NULL; }
                        | existing_region             { $$ = $1; }
                        | existing_object             { $$ = $1; }
;

opt_hit_spec: /* empty */                             { $$ = REPORT_NOTHING; }
            | ',' hit_spec                            { $$ = $2; }
;

hit_spec: FRONT_HITS                                  { $$ = REPORT_FRONT_HITS; }
        | BACK_HITS                                   { $$ = REPORT_BACK_HITS; }
        | ALL_HITS                                    { $$ = REPORT_ALL_HITS; }
        | FRONT_CROSSINGS                             { $$ = REPORT_FRONT_CROSSINGS; }
        | BACK_CROSSINGS                              { $$ = REPORT_BACK_CROSSINGS; }
        | ALL_CROSSINGS                               { $$ = REPORT_ALL_CROSSINGS; }
        | ESTIMATE_CONCENTRATION                      { $$ = REPORT_CONCENTRATION; }
        | ALL_ENCLOSED                                { $$ = REPORT_ENCLOSED; }
;

opt_custom_header: /* empty */                        { $$ = NULL; }
                 | ':' str_expr                       { $$ = $2; }
;

/* =================================================================== */
/* New-style viz output definitions */

viz_output_def:
          VIZ_OUTPUT '{'
            viz_output_maybe_mode_cmd
            viz_mesh_format_maybe_cmd
            viz_molecule_format_maybe_cmd
            list_viz_output_cmds
          '}'                                         {
                                                        if (volp->viz_mode == DREAMM_V3_MODE  ||
                                                            volp->viz_mode == DREAMM_V3_GROUPED_MODE)
                                                        {
                                                          if (volp->file_prefix_name == NULL)
                                                          {
                                                            mdlerror(mdlpvp, "Inside VIZ_OUTPUT block the required keyword FILENAME is missing.\n");
                                                            return 1;
                                                          }
                                                        }
                                                      }
;

list_viz_output_cmds:
          viz_output_cmd
        | list_viz_output_cmds
          viz_output_cmd
;

viz_output_maybe_mode_cmd:
          /* empty */                                 {
                                                        if (volp->viz_mode == -1)
                                                          volp->viz_mode = DREAMM_V3_MODE;
                                                      }
                         | viz_mode_def               {
                                                        if (volp->viz_mode != -1  &&  volp->viz_mode != $1)
                                                        {
                                                          mdlerror_fmt(mdlpvp, "Only one visualization mode is allowed in a given MDL file");
                                                          return 1;
                                                        }
                                                        volp->viz_mode = $1;
                                                      }
;

viz_mode_def: MODE '=' NONE                           { $$ = NO_VIZ_MODE; }
            | MODE '=' DX                             { $$ = DX_MODE; }
            | MODE '=' DREAMM_V3                      { $$ = DREAMM_V3_MODE; }
            | MODE '=' DREAMM_V3_GROUPED              { $$ = DREAMM_V3_GROUPED_MODE; }
            | MODE '=' CUSTOM_RK                      {
                                                        /* If we already had a CUSTOM_RK block, we're going to overwrite the volp->rk_mode_var data... */
                                                        if (volp->viz_mode == RK_MODE)
                                                        {
                                                          mdlerror_fmt(mdlpvp, "Only one CUSTOM_RK visualization block is allowed in an MDL file");
                                                          return 1;
                                                        }
                                                        volp->rk_mode_var = NULL;
                                                        $$ = RK_MODE;
                                                      }
            | MODE '=' CUSTOM_RK
              '[' list_range_specs ']' point          {
                                                        /* If we already had a CUSTOM_RK block, we're going to overwrite the volp->rk_mode_var data... */
                                                        if (volp->viz_mode == RK_MODE)
                                                        {
                                                          mdlerror_fmt(mdlpvp, "Only one CUSTOM_RK visualization block is allowed in an MDL file");
                                                          return 1;
                                                        }

                                                        CHECKN(volp->rk_mode_var = mdl_new_rk_mode_var(mdlpvp, $5.value_count, $5.value_head, $7));
                                                        $$ = RK_MODE;
                                                      }
            | MODE '=' ASCII                          { $$ = ASCII_MODE; }
;

viz_mesh_format_maybe_cmd: /* empty */                {
                                                        if (volp->viz_mode == DREAMM_V3_MODE)
                                                          volp->viz_output_flag |= VIZ_MESH_FORMAT_BINARY;
                                                      }
                         | viz_mesh_format_def        {
                                                        if (volp->viz_mode != DREAMM_V3_MODE)
                                                        {
                                                          mdlerror_fmt(mdlpvp, "VIZ_MESH_FORMAT command is allowed only in DREAMM_V3 mode.\n");
                                                          return 1;
                                                        }
                                                        volp->viz_output_flag |= $1;
                                                        if ((volp->viz_output_flag & VIZ_MESH_FORMAT_ASCII) &&
                                                            (volp->viz_output_flag & VIZ_MESH_FORMAT_BINARY))
                                                        {
                                                          mdlerror_fmt(mdlpvp, "BINARY and ASCII options for the VIZ_MESH_FORMAT command are mutually exclusive.\n");
                                                          return 1;
                                                        }
                                                      }
;

viz_mesh_format_def: VIZ_MESH_FORMAT '=' BINARY       { $$ = VIZ_MESH_FORMAT_BINARY; }
                   | VIZ_MESH_FORMAT '=' ASCII        { $$ = VIZ_MESH_FORMAT_ASCII; }
;

viz_molecule_format_maybe_cmd:
          /* empty */                                 {
                                                        if (volp->viz_mode == DREAMM_V3_MODE)
                                                          volp->viz_output_flag |= VIZ_MOLECULE_FORMAT_BINARY;
                                                      }
        | viz_molecule_format_def                     {
                                                        if (volp->viz_mode != DREAMM_V3_MODE)
                                                        {
                                                          mdlerror_fmt(mdlpvp, "VIZ_MOLECULE_FORMAT command is allowed only in DREAMM_V3 mode.\n");
                                                          return 1;
                                                        }
                                                        volp->viz_output_flag |= $1;
                                                        if ((volp->viz_output_flag & VIZ_MOLECULE_FORMAT_ASCII) &&
                                                            (volp->viz_output_flag & VIZ_MOLECULE_FORMAT_BINARY))
                                                        {
                                                          mdlerror_fmt(mdlpvp, "BINARY and ASCII options for the VIZ_MOLECULE_FORMAT command are mutually exclusive.\n");
                                                          return 1;
                                                        }
                                                      }
;

viz_molecule_format_def:
          VIZ_MOLECULE_FORMAT '=' BINARY              { $$ = VIZ_MOLECULE_FORMAT_BINARY; }
        | VIZ_MOLECULE_FORMAT '=' ASCII               { $$ = VIZ_MOLECULE_FORMAT_ASCII; }
;

viz_output_cmd:
          viz_filename_prefix_def
        | viz_frames_def                              {
                                                        if ($1.frame_head)
                                                        {
                                                          $1.frame_tail->next = mdlpvp->vol->frame_data_head;
                                                          mdlpvp->vol->frame_data_head = $1.frame_head;
                                                        }
                                                      }
        | viz_molecule_prefix_def
        | viz_object_prefixes_def
        | viz_state_values_def
;

viz_frames_def:
          viz_molecules_block_def
        | viz_meshes_block_def
        | viz_output_block_def
        | viz_iteration_frame_data_def
;

viz_filename_prefix_def: FILENAME '=' str_expr        {
                                                        volp->file_prefix_name = $3;
                                                        if (volp->viz_mode == DX_MODE)
                                                        {
                                                          volp->molecule_prefix_name = $3;
                                                        }
                                                      }
;

viz_molecules_block_def:
          MOLECULES '{'
            list_viz_molecules_block_cmds
          '}'                                         { $$ = $3; }
;

list_viz_molecules_block_cmds:
          viz_molecules_block_cmd
        | list_viz_molecules_block_cmds
          viz_molecules_block_cmd                     {
                                                        $$ = $1;
                                                        if ($$.frame_tail)
                                                        {
                                                          $$.frame_tail->next = $2.frame_head;
                                                          if ($2.frame_tail)
                                                            $$.frame_tail = $2.frame_tail;
                                                        }
                                                        else
                                                          $$ = $2;
                                                      }
;

viz_molecules_block_cmd:
          viz_molecules_name_list_cmd                 { $$.frame_head = $$.frame_tail = NULL; }
        | viz_molecules_time_points_def
        | viz_molecules_iteration_numbers_def
;

viz_molecules_name_list_cmd:
          NAME_LIST '{'
            viz_include_mols_cmd
          '}'
;

viz_include_mols_cmd:
          list_mol_name_specs
        | list_mol_name_specs_state_values
        | list_all_mols_specs
;

list_mol_name_specs:
          mol_name_spec
        | list_mol_name_specs mol_name_spec
;

existing_one_or_multiple_molecules:
          VAR                                         { CHECKN($$ = mdl_existing_molecule_list(mdlpvp, $1)); }
        | STR_VALUE                                   { CHECKN($$ = mdl_existing_molecules_wildcard(mdlpvp, $1)); }
;

mol_name_spec: existing_one_or_multiple_molecules     {
                                                        if (volp->viz_mode == DX_MODE)
                                                        {
                                                          mem_put_list(mdlpvp->sym_list_mem, $1);
                                                          mdlerror(mdlpvp, "In DX MODE the state value for the molecule should be specified.\n");
                                                          return 1;
                                                        }

                                                        /* Mark all specified molecules */
                                                        struct sym_table_list *stl;
                                                        for (stl = $1; stl != NULL; stl = stl->next)
                                                        {
                                                          struct species *specp = (struct species *) stl->node->value;
                                                          specp->viz_state = INCLUDE_OBJ;
                                                        }

                                                        /* free allocated memory  */
                                                        mem_put_list(mdlpvp->sym_list_mem, $1);
                                                      }
;

list_mol_name_specs_state_values:
          mol_name_specs_state_value
        | list_mol_name_specs_state_values
          mol_name_specs_state_value
;

mol_name_specs_state_value:
          existing_molecule '=' num_expr              {
                                                        struct species *specp = (struct species *) $1->value;

                                                        /* set the 'viz_state' value */
                                                        if ($1->sym_type == MOL)
                                                            specp->viz_state = (int) $3;
                                                        volp->viz_output_flag |= VIZ_MOLECULES_STATES;
                                                      }
;

list_all_mols_specs: ALL_MOLECULES                    {
                                                        if (volp->viz_mode == DX_MODE)
                                                        {
                                                          mdlerror(mdlpvp, "In DX MODE the state value for the molecule should be specified.\n");
                                                          return 1;
                                                        }
                                                        volp->viz_output_flag |= VIZ_ALL_MOLECULES;
                                                      }
;


viz_time_spec:
          ALL_TIMES                                   { CHECK(mdl_new_viz_all_times(mdlpvp, & $$)); }
        | '[' list_range_specs ']'                    { $$ = $2; }
;
viz_molecules_time_points_def:
          TIME_POINTS '{'
            viz_molecules_time_points_cmds
          '}'                                         { $$ = $3; }
;

viz_molecules_time_points_cmds:
          viz_molecules_time_points_one_cmd
        | viz_molecules_time_points_cmds
          viz_molecules_time_points_one_cmd           { $$ = $1; $$.frame_tail->next = $2.frame_head; $$.frame_tail = $2.frame_tail; }
;

viz_molecules_time_points_one_cmd:
          viz_molecules_one_item '@'
          viz_time_spec                               {
                                                        struct frame_data_list *fdlp;
                                                        CHECKN(fdlp = mdl_create_viz_mol_frames(mdlpvp,
                                                                                                OUTPUT_BY_TIME_LIST,
                                                                                                $1,
                                                                                                volp->viz_mode,
                                                                                                $3.value_head));
                                                        $$.frame_head = fdlp;
                                                        while (fdlp->next != NULL)
                                                          fdlp = fdlp->next;
                                                        $$.frame_tail = fdlp;
                                                      }
;

viz_iteration_spec:
          ALL_ITERATIONS                              { CHECK(mdl_new_viz_all_iterations(mdlpvp, & $$)); }
        | '[' list_range_specs ']'                    { $$ = $2; }
;

viz_molecules_iteration_numbers_def:
          ITERATION_NUMBERS '{'
            viz_molecules_iteration_numbers_cmds
          '}'                                         { $$ = $3; }
;

viz_molecules_iteration_numbers_cmds:
          viz_molecules_iteration_numbers_one_cmd
        | viz_molecules_iteration_numbers_cmds
          viz_molecules_iteration_numbers_one_cmd     { $$ = $1; $$.frame_tail->next = $2.frame_head; $$.frame_tail = $2.frame_tail; }
;

viz_molecules_iteration_numbers_one_cmd:
          viz_molecules_one_item '@'
          viz_iteration_spec                          {
                                                        struct frame_data_list *fdlp;
                                                        CHECKN(fdlp = mdl_create_viz_mol_frames(mdlpvp,
                                                                                                OUTPUT_BY_ITERATION_LIST,
                                                                                                $1,
                                                                                                volp->viz_mode,
                                                                                                $3.value_head));
                                                        $$.frame_head = fdlp;
                                                        while (fdlp->next != NULL)
                                                          fdlp = fdlp->next;
                                                        $$.frame_tail = fdlp;
                                                      }
;

viz_molecules_one_item: ALL_DATA                      { $$ = ALL_MOL_DATA; }
                      | POSITIONS                     { $$ = MOL_POS; }
                      | ORIENTATIONS                  { $$ = MOL_ORIENT; }
;

viz_meshes_block_def:
          MESHES '{'
            list_viz_meshes_block_cmds
          '}'                                         { $$ = $3; }
;

list_viz_meshes_block_cmds:
          viz_meshes_block_cmd
        | list_viz_meshes_block_cmds
          viz_meshes_block_cmd                        {
                                                        $$ = $1;
                                                        if ($$.frame_tail)
                                                          $$.frame_tail->next = $2.frame_head;
                                                        if ($2.frame_tail)
                                                          $$.frame_tail = $2.frame_tail;
                                                      }
;

viz_meshes_block_cmd:
          viz_meshes_name_list_cmd                    { $$.frame_head = $$.frame_tail = NULL; }
        | viz_meshes_time_points_def
        | viz_meshes_iteration_numbers_def
;

viz_meshes_name_list_cmd:
          NAME_LIST '{'
            viz_include_meshes_cmd
          '}'
;

viz_include_meshes_cmd:
          list_meshes_name_specs
        | list_meshes_name_specs_state_values
        | list_all_meshes_specs
;

list_meshes_name_specs:
          mesh_one_name_spec
        | list_meshes_name_specs
          mesh_one_name_spec
;

mesh_one_name_spec:
          mesh_object_or_wildcard                     {
                                                        struct sym_table_list *stl;
                                                        if (volp->viz_mode == DX_MODE)
                                                        {
                                                          mdlerror(mdlpvp, "In DX MODE the state value for the object should be specified.\n");
                                                          return 1;
                                                        }

                                                        for (stl = $1; stl != NULL; stl = stl->next)
                                                        {
                                                          /* create viz_obj object */
                                                          struct object *objp = (struct object *) stl->node->value;
                                                          if((objp->object_type == REL_SITE_OBJ))
                                                            continue;
                                                          CHECK(mdl_add_viz_object(mdlpvp, stl->node, INCLUDE_OBJ));
                                                        }
                                                        mem_put_list(mdlpvp->sym_list_mem, $1);
                                                      }
;

list_meshes_name_specs_state_values:
          mesh_one_name_spec_state_value
        | list_meshes_name_specs_state_values
          mesh_one_name_spec_state_value
;

mesh_one_name_spec_state_value:
          existing_object '=' num_expr                {
                                                        volp->viz_output_flag |= VIZ_SURFACE_STATES;
                                                        CHECK(mdl_add_viz_object(mdlpvp, $1, (int) $3));
                                                      }
;

list_all_meshes_specs: ALL_MESHES                     {
                                                        struct object *o;
                                                        if (volp->viz_mode == DX_MODE)
                                                        {
                                                          mdlerror(mdlpvp, "The keyword ALL_MESHES cannot be used in DX MODE.\n");
                                                          return 1;
                                                        }

                                                        for (o = volp->root_instance->first_child; o != NULL; o = o->next)
                                                          CHECK(mdl_add_viz_object(mdlpvp, o->sym, INCLUDE_OBJ));
                                                      }
;

viz_meshes_time_points_def:
          TIME_POINTS '{'
            viz_meshes_time_points_cmds
          '}'                                         { $$ = $3; }
;

viz_meshes_time_points_cmds:
          viz_meshes_time_points_one_cmd
        | viz_meshes_time_points_cmds
          viz_meshes_time_points_one_cmd              { $$ = $1; $$.frame_tail->next = $2.frame_head; $$.frame_tail = $2.frame_tail; }
;

viz_meshes_time_points_one_cmd:
          viz_meshes_one_item '@'
          viz_time_spec                               {
                                                        struct frame_data_list *fdlp;
                                                        CHECKN(fdlp = mdl_create_viz_mesh_frames(mdlpvp,
                                                                                                 OUTPUT_BY_TIME_LIST,
                                                                                                 $1,
                                                                                                 volp->viz_mode,
                                                                                                 $3.value_head));
                                                        $$.frame_head = fdlp;
                                                        while (fdlp->next != NULL)
                                                          fdlp = fdlp->next;
                                                        $$.frame_tail = fdlp;
                                                      }
;

viz_meshes_iteration_numbers_def:
          ITERATION_NUMBERS '{'
            viz_meshes_iteration_numbers_cmds
          '}'                                         { $$ = $3; }
;

viz_meshes_iteration_numbers_cmds:
          viz_meshes_iteration_numbers_one_cmd
        | viz_meshes_iteration_numbers_cmds
          viz_meshes_iteration_numbers_one_cmd        { $$ = $1; $$.frame_tail->next = $2.frame_head; $$.frame_tail = $2.frame_tail; }
;

viz_meshes_iteration_numbers_one_cmd:
          viz_meshes_one_item '@'
          viz_iteration_spec                          {
                                                        struct frame_data_list *fdlp;
                                                        CHECKN(fdlp = mdl_create_viz_mesh_frames(mdlpvp,
                                                                                                 OUTPUT_BY_ITERATION_LIST,
                                                                                                 $1,
                                                                                                 volp->viz_mode,
                                                                                                 $3.value_head));
                                                        $$.frame_head = fdlp;
                                                        while (fdlp->next != NULL)
                                                          fdlp = fdlp->next;
                                                        $$.frame_tail = fdlp;
                                                      }
;

viz_meshes_one_item: ALL_DATA                         { $$ = ALL_MESH_DATA; }
                   | GEOMETRY                         { $$ = MESH_GEOMETRY; }
                   | REGION_DATA                      { $$ = REG_DATA; }
;

/* =================================================================== */
/* Old-style viz output definitions */

viz_data_output_def:
          VIZ_DATA_OUTPUT '{'
            viz_output_maybe_mode_cmd
            list_viz_data_output_cmds
          '}'
;


list_viz_data_output_cmds:
          /* empty */
        | list_viz_data_output_cmds
          viz_data_output_cmd
;


viz_data_output_cmd:
          viz_frames_def_old                          {
                                                        if ($1.frame_head)
                                                        {
                                                          $1.frame_tail->next = mdlpvp->vol->frame_data_head;
                                                          mdlpvp->vol->frame_data_head = $1.frame_head;
                                                        }
                                                      }
        | viz_molecule_prefix_def
        | viz_object_prefixes_def
        | viz_state_values_def
        | voxel_image_mode_def
        | voxel_volume_mode_def
;

viz_frames_def_old:
          viz_output_block_def
        | viz_iteration_frame_data_def
;

voxel_image_mode_def:
          VOXEL_IMAGE_MODE '=' boolean                { volp->voxel_image_mode = $3; }
;

voxel_volume_mode_def:
          VOXEL_VOLUME_MODE '=' boolean               { volp->voxel_volume_mode = $3; }
;

viz_output_block_def:
          viz_iteration_def
        | viz_time_def
;

viz_iteration_def:
          ITERATION_LIST '='
          '[' list_range_specs ']'                    {
                                                        struct frame_data_list *fdlp;
                                                        mdl_sort_numeric_list($4.value_head);
                                                        CHECKN(fdlp = mdl_create_viz_frame(mdlpvp,
                                                                                           OUTPUT_BY_ITERATION_LIST,
                                                                                           ALL_FRAME_DATA,
                                                                                           $4.value_head));
                                                        $$.frame_tail = $$.frame_head = fdlp;
                                                      }
;

viz_time_def:
          TIME_LIST '='
          '[' list_range_specs ']'                    {
                                                        struct frame_data_list *fdlp;
                                                        mdl_sort_numeric_list($4.value_head);
                                                        CHECKN(fdlp = mdl_create_viz_frame(mdlpvp,
                                                                                           OUTPUT_BY_TIME_LIST,
                                                                                           ALL_FRAME_DATA,
                                                                                           $4.value_head));
                                                        $$.frame_tail = $$.frame_head = fdlp;
                                                      }
;

viz_iteration_frame_data_def:
          ITERATION_FRAME_DATA '{'
            list_iteration_frame_data_specs
          '}'                                         { $$ = $3; }
;

list_iteration_frame_data_specs:
          iteration_frame_data_spec
        | list_iteration_frame_data_specs
          iteration_frame_data_spec                   { $$ = $1; $$.frame_tail->next = $2.frame_head; $$.frame_tail = $2.frame_tail; }
;

iteration_frame_data_spec:
          iteration_frame_data_item '='
          '[' list_range_specs ']'                    {
                                                        struct frame_data_list *fdlp;
                                                        mdl_sort_numeric_list($4.value_head);
                                                        CHECKN(fdlp = mdl_create_viz_frame(mdlpvp,
                                                                                           OUTPUT_BY_ITERATION_LIST,
                                                                                           $1,
                                                                                           $4.value_head));
                                                        $$.frame_tail = $$.frame_head = fdlp;
                                                      }
;

iteration_frame_data_item:
          ALL_DATA                                    { $$ = ALL_FRAME_DATA; }
        | EFFECTOR_POSITIONS                          { $$ = EFF_POS;        }
        | EFFECTOR_STATES                             { $$ = EFF_STATES;     }
        | MOLECULE_POSITIONS                          { $$ = MOL_POS;        }
        | MOLECULE_STATES                             { $$ = MOL_STATES;     }
        | SURFACE_POSITIONS                           { $$ = SURF_POS;       }
        | SURFACE_STATES                              { $$ = SURF_STATES;    }
;

viz_molecule_prefix_def:
          MOLECULE_FILE_PREFIX '=' str_expr           {
                                                        if (volp->viz_mode == DREAMM_V3_MODE  ||
                                                            volp->viz_mode == DREAMM_V3_GROUPED_MODE)
                                                        {
                                                          mdlerror(mdlpvp, "MOLECULE_FILE_PREFIX canot be used with the DREAMM viz output modes");
                                                          return 1;
                                                        }
                                                        volp->molecule_prefix_name = $3;
                                                      }
;

viz_object_prefixes_def:
          OBJECT_FILE_PREFIXES '{'
            list_viz_object_prefixes
          '}'                                         {
                                                        if (volp->viz_mode == DREAMM_V3_MODE  ||
                                                            volp->viz_mode == DREAMM_V3_GROUPED_MODE)
                                                        {
                                                          mdlerror(mdlpvp, "OBJECT_FILE_PREFIXES canot be used with the DREAMM viz output modes");
                                                          return 1;
                                                        }
                                                      }
;

list_viz_object_prefixes:
          viz_object_prefix
        | list_viz_object_prefixes
          viz_object_prefix
;

viz_object_prefix: existing_object '=' str_expr       {
                                                        struct object *objp;
                                                        struct viz_obj *vizp;
                                                        objp = (struct object *) $1->value;
                                                        if ((vizp = (struct viz_obj *) malloc(sizeof(struct viz_obj))) == NULL)
                                                        {
                                                          mdlerror(mdlpvp, "Out of memory while creating viz object");
                                                          return 1;
                                                        }
                                                        objp->viz_obj = vizp;
                                                        vizp->name = $3;
                                                        vizp->full_name = mdl_strdup(mdlpvp, $1->name);
                                                        if (vizp->full_name == NULL)
                                                          return 1;
                                                        vizp->obj = objp;
                                                        vizp->viz_child_head = NULL;
                                                        vizp->next = volp->viz_obj_head;
                                                        volp->viz_obj_head = vizp;
                                                      }
;

viz_state_values_def:
          STATE_VALUES '{'
            list_viz_state_values
          '}'
;

list_viz_state_values:
          viz_state_value
        | list_viz_state_values
          viz_state_value
;

existing_logicalOrPhysical: VAR                       { CHECKN($$ = mdl_existing_molecule_or_object(mdlpvp, $1)); }
;

viz_state_value:
          existing_logicalOrPhysical '='
          num_expr                                    {
                                                        struct species *specp;
                                                        int viz_state = (int) $3;
                                                        switch ($1->sym_type)
                                                        {
                                                          case OBJ:
                                                            CHECK(mdl_set_object_viz_state(mdlpvp, $1, viz_state));
                                                            break;

                                                          case MOL:
                                                            specp = (struct species *) $1->value;
                                                            specp->viz_state = viz_state;
                                                            break;
                                                        }
                                                      }
        | existing_region '=' num_expr                { CHECK(mdl_set_region_viz_state(mdlpvp, (struct region *) $1->value, (int) $3)); }
;

/* =================================================================== */
/* Volume data output mode */

volume_output_def:
          VOLUME_DATA_OUTPUT '{'
            volume_output_filename_prefix
            volume_output_molecule_list
            volume_output_location
            volume_output_voxel_size
            volume_output_voxel_count
            volume_output_times_def
          '}'                                         {
                                                        struct volume_output_item *vo;
                                                        CHECKN(vo = mdl_new_volume_output_item(mdlpvp, $3, & $4, $5, $6, $7, $8));
                                                        vo->next = mdlpvp->vol->volume_output_head;
                                                        mdlpvp->vol->volume_output_head = vo;
                                                      }
;

volume_output_filename_prefix:
          FILENAME_PREFIX '=' str_expr                { $$ = $3; }
;

volume_output_molecule_list:
          volume_output_molecule_decl
        | volume_output_molecule_list
          volume_output_molecule_decl                 {
                                                        $$ = $1;
                                                        $$.species_count += $2.species_count;
                                                        $$.species_tail->next = $2.species_head;
                                                        $$.species_tail = $2.species_tail;
                                                      }
;

volume_output_molecule_decl:
          MOLECULES '=' volume_output_molecules       { $$ = $3; }
;

volume_output_molecule: VAR                           {
                                                        struct sym_table *sp;
                                                        struct species_list_item *ptrl;
                                                        CHECKN(sp = mdl_existing_molecule(mdlpvp, $1));

                                                        ptrl = (struct species_list_item *) mem_get(mdlpvp->species_list_mem);
                                                        if (ptrl == NULL)
                                                        {
                                                          mdlerror_fmt(mdlpvp, "Out of memory while parsing molecule list");
                                                          return 1;
                                                        }
                                                        ptrl->spec = (struct species *) sp->value;
                                                        ptrl->next = NULL;
                                                        $$ = ptrl;
                                                      }
;

volume_output_molecules:
          volume_output_molecule                      { $$.species_tail = $$.species_head = $1; $$.species_count = 1; }
        | volume_output_molecules '+'
          volume_output_molecule                      { $$ = $1; $$.species_tail = $$.species_tail->next = $3; ++ $$.species_count; }
;

volume_output_location:
          LOCATION '=' point                          { $$ = $3; }
;

volume_output_voxel_size:
          VOXEL_SIZE '=' point_or_num                 { $$ = $3; }
;

volume_output_voxel_count:
          VOXEL_COUNT '=' point_or_num                {
                                                        if ($3->x < 1.0)
                                                        {
                                                          mdl_warning(mdlpvp, "Voxel count (x dimension) too small.  Setting x count to 1.");
                                                          $3->x = 1.0;
                                                        }
                                                        if ($3->y < 1.0)
                                                        {
                                                          mdl_warning(mdlpvp, "Voxel count (y dimension) too small.  Setting y count to 1.");
                                                          $3->y = 1.0;
                                                        }
                                                        if ($3->z < 1.0)
                                                        {
                                                          mdl_warning(mdlpvp, "Voxel count (z dimension) too small.  Setting z count to 1.");
                                                          $3->z = 1.0;
                                                        }
                                                        $$ = $3;
                                                      }
;

volume_output_times_def:
          /* empty */                                 { CHECKN($$ = mdl_new_output_times_default(mdlpvp)); }
        | STEP '=' num_expr                           { CHECKN($$ = mdl_new_output_times_step(mdlpvp, $3)); }
        | ITERATION_LIST '=' '[' list_range_specs ']' { CHECKN($$ = mdl_new_output_times_iterations(mdlpvp, & $4)); }
        | TIME_LIST '=' '[' list_range_specs ']'      { CHECKN($$ = mdl_new_output_times_time(mdlpvp, & $4)); }
;

%%





/* Begin Bison Epilogue: */

/* mdlerror: Standard error callback from parser.
 *
 *   mpvp: the parser state variables
 *   str:  the error message to display
 */
void mdlerror(struct mdlparse_vars *mpvp, char const *str)
{
  mdlerror_fmt(mpvp, "%s", str);
}

/* mdlerror_fmt: Print a formatted error message regarding an error in the MDL
 *               file.
 *
 *   mpvp: the parser state variables
 *   fmt:  the printf-style format string
 */
void mdlerror_fmt(struct mdlparse_vars *mpvp, char const *fmt, ...)
{
  va_list arglist;
  FILE *log_file = stderr;

  if (mpvp->vol->procnum != 0)
    return;

  if (mpvp->vol->err_file != NULL)
    log_file = mpvp->vol->err_file;

  /* print error location */
  if (mpvp->include_stack_ptr == 0)
    fprintf(log_file,
            "MCell: error occurred after parsing file: %s\n  ",
            mpvp->vol->curr_file);
  else
    fprintf(log_file,
            "MCell: error on line: %d of file: %s\n  ",
            mpvp->line_num[mpvp->include_stack_ptr - 1], mpvp->vol->curr_file);

  /* format error message */
  va_start(arglist, fmt);
  vfprintf(log_file, fmt, arglist);
  va_end(arglist);

  /* terminate error message and flush */
  fprintf(log_file, "\n");
  fflush(log_file);
}

/* mdlerror_file: Open and parse an MDL file.
 *
 *   mpvp: the parser state variables
 *   name: the path to the MDL file
 */
static int mdlparse_file(struct mdlparse_vars *mpvp, char const *name)
{
  int failure;
  int cur_stack = mpvp->include_stack_ptr ++;
  FILE *infile;
  yyscan_t scanner;
  char const *prev_file;

  /* Put filename and line number on stack */
  if (cur_stack >= MAX_INCLUDE_DEPTH)
  {
    mdlerror_fmt(mpvp, "Includes nested too deeply at file %s, included from %s:%d\n", name, mpvp->include_filename[cur_stack-1], mpvp->line_num[cur_stack-1]);
    return 1;
  }
  mpvp->line_num[cur_stack] = 1;
  mpvp->include_filename[cur_stack] = name;

  /* Open file, or know the reason why */
  no_printf("Opening file %s\n", name);
  fflush(mpvp->vol->err_file);
  if ((infile = fopen(name,"r")) == NULL)
  {
    int err = errno;
    mdlerror_fmt(mpvp, "Couldn't open file %s, included from %s:%d: %s\n", name, mpvp->include_filename[cur_stack-1], mpvp->line_num[cur_stack-1], strerror(err));
    -- mpvp->include_stack_ptr;
    return 1;
  }

  /* Create and initialize a lexer */
  if (mdllex_init(&scanner))
  {
    int err = errno;
    if (err == ENOMEM)
      mdlerror_fmt(mpvp, "Couldn't initialize lexer for file %s, included from %s:%d: out of memory\n", name, mpvp->include_filename[cur_stack-1], mpvp->line_num[cur_stack-1]);
    else if (err == EINVAL)
      mdlerror_fmt(mpvp, "Couldn't initialize lexer for file %s, included from %s:%d: internal error (invalid argument)\n", name, mpvp->include_filename[cur_stack-1], mpvp->line_num[cur_stack-1]);
    else
      mdlerror_fmt(mpvp, "Couldn't initialize lexer for file %s, included from %s:%d: internal error\n", name, mpvp->include_filename[cur_stack-1], mpvp->line_num[cur_stack-1]);
    fclose(infile);
    -- mpvp->include_stack_ptr;
    return 1;
  }
  mdlrestart(infile, scanner);

  /* Parse this file */
  prev_file = mpvp->vol->curr_file;
  mpvp->vol->curr_file = name;
  failure = mdlparse(mpvp, scanner);
  mpvp->vol->curr_file = prev_file;
  -- mpvp->include_stack_ptr;

  /* Clean up! */
  fclose(infile);
  mdllex_destroy(scanner);

  return failure;
}

/* mdlerror_init: Set up and parse the top-level MDL file.
 *
 *   vol: the world to populate
 */
int mdlparse_init(struct volume *vol)
{
  int failure;
  struct mdlparse_vars mpv;

  vol->initialization_state = "parsing";
  memset(&mpv, 0, sizeof(struct mdlparse_vars));

  mpv.vol=vol;
  mpv.include_stack_ptr=0;
  mpv.current_object = vol->root_object;
  mpv.vol->macro_count_request_head = NULL;

  /* Create memory pools for parsing */
  if ((mpv.path_mem=create_mem(sizeof(struct pathway),4096))==NULL)
  {
    fprintf(vol->log_file, "Out of memory while creating reaction pathways\n");
    return 1;
  }
  if ((mpv.prod_mem=create_mem(sizeof(struct product),4096))==NULL) {
    fprintf(vol->log_file, "Out of memory while creating reaction products\n");
    return 1;
  }
  mpv.sym_list_mem = create_mem(sizeof(struct sym_table_list),4096);
  if (mpv.sym_list_mem==NULL)
  {
    fprintf(vol->err_file,"Out of memory while getting ready to store lists of molecules and reactions");
    return 1;
  }
  mpv.species_list_mem = create_mem(sizeof(struct species_list_item), 1024);
  if (mpv.species_list_mem == NULL)
  {
    fprintf(vol->err_file,"Out of memory while allocating temporary space for species list parsing");
    return 1;
  }
  mpv.mol_data_list_mem = create_mem(sizeof(struct species_opt_orient), 1024);
  if (mpv.mol_data_list_mem == NULL)
  {
    fprintf(vol->err_file,"Out of memory while allocating temporary space for species list parsing");
    return 1;
  }
  mpv.output_times_mem = create_mem(sizeof(struct output_times), 1024);
  if (mpv.output_times_mem == NULL)
  {
    fprintf(vol->err_file,"Out of memory while allocating temporary space for volume output parsing");
    return 1;
  }

  /* Start parsing at the top-level file */
  vol->curr_file = vol->mdl_infile_name;
  failure = mdlparse_file(&mpv, vol->mdl_infile_name);

  /* If we succeeded, prepare the reactions */
  if (! failure  &&  prepare_reactions(&mpv))
  {
    mdlerror(&mpv, "Failed to initialize reactions");
    failure = 1;
  }

  /* Free leftover object names */
  while (mpv.object_name_list)
  {
    struct name_list *l = mpv.object_name_list->next;
    free(mpv.object_name_list);
    mpv.object_name_list = l;
  }

  /* Destroy memory pools */
  delete_mem(mpv.species_list_mem);
  delete_mem(mpv.mol_data_list_mem);
  delete_mem(mpv.output_times_mem);
  delete_mem(mpv.sym_list_mem);
  delete_mem(mpv.prod_mem);
  delete_mem(mpv.path_mem);

  vol->initialization_state = "initializing";

  return failure;
}
