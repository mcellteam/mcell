/* A Bison parser, made by GNU Bison 3.0.4.  */

/* Bison interface for Yacc-like parsers in C

   Copyright (C) 1984, 1989-1990, 2000-2015 Free Software Foundation, Inc.

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.  */

/* As a special exception, you may create a larger work that contains
   part or all of the Bison parser skeleton and distribute that work
   under terms of your choice, so long as that work isn't itself a
   parser generator using the skeleton or a modified version thereof
   as a parser skeleton.  Alternatively, if you modify or redistribute
   the parser skeleton itself, you may (at your option) remove this
   special exception, which will cause the skeleton and the resulting
   Bison output files to be licensed under the GNU General Public
   License without this special exception.

   This special exception was added by the Free Software Foundation in
   version 2.2 of Bison.  */

#ifndef YY_MDL_HOME_JCZECH_MCELL_BUILD_DEPS_MDLPARSE_H_INCLUDED
# define YY_MDL_HOME_JCZECH_MCELL_BUILD_DEPS_MDLPARSE_H_INCLUDED
/* Debug traces.  */
#ifndef YYDEBUG
# define YYDEBUG 0
#endif
#if YYDEBUG
extern int mdldebug;
#endif

/* Token type.  */
#ifndef YYTOKENTYPE
# define YYTOKENTYPE
  enum yytokentype
  {
    ABS = 258,
    ABSORPTIVE = 259,
    ACCURATE_3D_REACTIONS = 260,
    ACOS = 261,
    ALL_CROSSINGS = 262,
    ALL_DATA = 263,
    ALL_ELEMENTS = 264,
    ALL_ENCLOSED = 265,
    ALL_HITS = 266,
    ALL_ITERATIONS = 267,
    ALL_MOLECULES = 268,
    ALL_NOTIFICATIONS = 269,
    ALL_TIMES = 270,
    ALL_WARNINGS = 271,
    ASCII = 272,
    ASIN = 273,
    ASPECT_RATIO = 274,
    ATAN = 275,
    BACK = 276,
    BACK_CROSSINGS = 277,
    BACK_HITS = 278,
    BOTTOM = 279,
    BOX = 280,
    BOX_TRIANGULATION_REPORT = 281,
    BRIEF = 282,
    CEIL = 283,
    CELLBLENDER = 284,
    CENTER_MOLECULES_ON_GRID = 285,
    CHECKPOINT_INFILE = 286,
    CHECKPOINT_ITERATIONS = 287,
    CHECKPOINT_OUTFILE = 288,
    CHECKPOINT_REALTIME = 289,
    CHECKPOINT_REPORT = 290,
    CLAMP_CONCENTRATION = 291,
    CLOSE_PARTITION_SPACING = 292,
    CONCENTRATION = 293,
    CORNERS = 294,
    COS = 295,
    COUNT = 296,
    CUBIC = 297,
    CUBIC_RELEASE_SITE = 298,
    CUSTOM_SPACE_STEP = 299,
    CUSTOM_TIME_STEP = 300,
    DEFINE_MOLECULE = 301,
    DEFINE_MOLECULES = 302,
    DEFINE_REACTIONS = 303,
    DEFINE_RELEASE_PATTERN = 304,
    DEFINE_SURFACE_CLASS = 305,
    DEFINE_SURFACE_CLASSES = 306,
    DEFINE_SURFACE_REGIONS = 307,
    DEGENERATE_POLYGONS = 308,
    DELAY = 309,
    DENSITY = 310,
    DIFFUSION_CONSTANT_2D = 311,
    DIFFUSION_CONSTANT_3D = 312,
    DIFFUSION_CONSTANT_REPORT = 313,
    DYNAMIC_GEOMETRY = 314,
    DYNAMIC_GEOMETRY_MOLECULE_PLACEMENT = 315,
    EFFECTOR_GRID_DENSITY = 316,
    ELEMENT_CONNECTIONS = 317,
    ELLIPTIC = 318,
    ELLIPTIC_RELEASE_SITE = 319,
    EQUAL = 320,
    ERROR = 321,
    ESTIMATE_CONCENTRATION = 322,
    EXCLUDE_ELEMENTS = 323,
    EXCLUDE_PATCH = 324,
    EXCLUDE_REGION = 325,
    EXIT = 326,
    EXP = 327,
    EXPRESSION = 328,
    EXTERN = 329,
    FALSE = 330,
    FCLOSE = 331,
    FILENAME = 332,
    FILENAME_PREFIX = 333,
    FILE_OUTPUT_REPORT = 334,
    FINAL_SUMMARY = 335,
    FLOOR = 336,
    FOPEN = 337,
    FORMAT = 338,
    FPRINTF = 339,
    FPRINT_TIME = 340,
    FRONT = 341,
    FRONT_CROSSINGS = 342,
    FRONT_HITS = 343,
    GAUSSIAN_RELEASE_NUMBER = 344,
    GEOMETRY = 345,
    GRAPH_PATTERN = 346,
    HEADER = 347,
    HIGH_PROBABILITY_THRESHOLD = 348,
    HIGH_REACTION_PROBABILITY = 349,
    IGNORED = 350,
    INCLUDE_ELEMENTS = 351,
    INCLUDE_FILE = 352,
    INCLUDE_PATCH = 353,
    INCLUDE_REGION = 354,
    INPUT_FILE = 355,
    INSTANTIATE = 356,
    LLINTEGER = 357,
    FULLY_RANDOM = 358,
    INTERACTION_RADIUS = 359,
    ITERATION_LIST = 360,
    ITERATION_NUMBERS = 361,
    ITERATION_REPORT = 362,
    ITERATIONS = 363,
    KEEP_CHECKPOINT_FILES = 364,
    LEFT = 365,
    LIFETIME_THRESHOLD = 366,
    LIFETIME_TOO_SHORT = 367,
    LIST = 368,
    LOCATION = 369,
    LOG = 370,
    LOG10 = 371,
    MAX_TOK = 372,
    MAXIMUM_STEP_LENGTH = 373,
    MEAN_DIAMETER = 374,
    MEAN_NUMBER = 375,
    MEMORY_PARTITION_X = 376,
    MEMORY_PARTITION_Y = 377,
    MEMORY_PARTITION_Z = 378,
    MEMORY_PARTITION_POOL = 379,
    MICROSCOPIC_REVERSIBILITY = 380,
    MIN_TOK = 381,
    MISSED_REACTIONS = 382,
    MISSED_REACTION_THRESHOLD = 383,
    MISSING_SURFACE_ORIENTATION = 384,
    MOD = 385,
    MODE = 386,
    MODIFY_SURFACE_REGIONS = 387,
    MOLECULE = 388,
    MOLECULE_COLLISION_REPORT = 389,
    MOLECULE_DENSITY = 390,
    MOLECULE_NUMBER = 391,
    MOLECULE_POSITIONS = 392,
    MOLECULES = 393,
    MOLECULE_PLACEMENT_FAILURE = 394,
    NAME_LIST = 395,
    NEAREST_POINT = 396,
    NEAREST_TRIANGLE = 397,
    NEGATIVE_DIFFUSION_CONSTANT = 398,
    NEGATIVE_REACTION_RATE = 399,
    NO = 400,
    NOEXIT = 401,
    NONE = 402,
    NO_SPECIES = 403,
    NOT_EQUAL = 404,
    NOTIFICATIONS = 405,
    NUMBER_OF_SUBUNITS = 406,
    NUMBER_OF_TRAINS = 407,
    NUMBER_TO_RELEASE = 408,
    OBJECT = 409,
    OFF = 410,
    ON = 411,
    ORIENTATIONS = 412,
    OUTPUT_BUFFER_SIZE = 413,
    INVALID_OUTPUT_STEP_TIME = 414,
    LARGE_MOLECULAR_DISPLACEMENT = 415,
    ADD_REMOVE_MESH = 416,
    OVERWRITTEN_OUTPUT_FILE = 417,
    PARTITION_LOCATION_REPORT = 418,
    PARTITION_X = 419,
    PARTITION_Y = 420,
    PARTITION_Z = 421,
    PERIODIC_BOX = 422,
    PERIODIC_X = 423,
    PERIODIC_Y = 424,
    PERIODIC_Z = 425,
    PERIODIC_TRADITIONAL = 426,
    PI_TOK = 427,
    POLYGON_LIST = 428,
    POSITIONS = 429,
    PRINTF = 430,
    PRINT_TIME = 431,
    PROBABILITY_REPORT = 432,
    PROBABILITY_REPORT_THRESHOLD = 433,
    PROGRESS_REPORT = 434,
    RADIAL_DIRECTIONS = 435,
    RADIAL_SUBDIVISIONS = 436,
    RAND_GAUSSIAN = 437,
    RAND_UNIFORM = 438,
    REACTION_DATA_OUTPUT = 439,
    REACTION_OUTPUT_REPORT = 440,
    REAL = 441,
    RECTANGULAR_RELEASE_SITE = 442,
    RECTANGULAR_TOKEN = 443,
    REFLECTIVE = 444,
    RELEASE_EVENT_REPORT = 445,
    RELEASE_INTERVAL = 446,
    RELEASE_PATTERN = 447,
    RELEASE_PROBABILITY = 448,
    RELEASE_SITE = 449,
    REMOVE_ELEMENTS = 450,
    RIGHT = 451,
    ROTATE = 452,
    ROUND_OFF = 453,
    SCALE = 454,
    SEED = 455,
    SHAPE = 456,
    SHOW_EXACT_TIME = 457,
    SIN = 458,
    SITE_DIAMETER = 459,
    SITE_RADIUS = 460,
    SPACE_STEP = 461,
    SPHERICAL = 462,
    SPHERICAL_RELEASE_SITE = 463,
    SPHERICAL_SHELL = 464,
    SPHERICAL_SHELL_SITE = 465,
    SPRINTF = 466,
    SQRT = 467,
    STANDARD_DEVIATION = 468,
    PERIODIC_BOX_INITIAL = 469,
    STEP = 470,
    STRING_TO_NUM = 471,
    STR_VALUE = 472,
    SUBUNIT = 473,
    SUBUNIT_RELATIONSHIPS = 474,
    SUMMATION_OPERATOR = 475,
    SURFACE_CLASS = 476,
    SURFACE_ONLY = 477,
    TAN = 478,
    TARGET_ONLY = 479,
    TET_ELEMENT_CONNECTIONS = 480,
    THROUGHPUT_REPORT = 481,
    TIME_LIST = 482,
    TIME_POINTS = 483,
    TIME_STEP = 484,
    TIME_STEP_MAX = 485,
    TO = 486,
    TOP = 487,
    TRAIN_DURATION = 488,
    TRAIN_INTERVAL = 489,
    TRANSLATE = 490,
    TRANSPARENT = 491,
    TRIGGER = 492,
    TRUE = 493,
    UNLIMITED = 494,
    USELESS_VOLUME_ORIENTATION = 495,
    VACANCY_SEARCH_DISTANCE = 496,
    VAR = 497,
    VARYING_PROBABILITY_REPORT = 498,
    VERTEX_LIST = 499,
    VIZ_OUTPUT = 500,
    VIZ_OUTPUT_REPORT = 501,
    VIZ_VALUE = 502,
    VOLUME_DATA_OUTPUT = 503,
    VOLUME_OUTPUT_REPORT = 504,
    VOLUME_DEPENDENT_RELEASE_NUMBER = 505,
    VOLUME_ONLY = 506,
    VOXEL_COUNT = 507,
    VOXEL_LIST = 508,
    VOXEL_SIZE = 509,
    WARNING = 510,
    WARNINGS = 511,
    WORLD = 512,
    YES = 513,
    UNARYMINUS = 514
  };
#endif

/* Value type.  */
#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED

union YYSTYPE
{
#line 67 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1909  */

int ival;
int tok;
double dbl;
long long llival;
char *str;
struct sym_entry *sym;
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
struct mcell_species mol_type;
struct mcell_species_list mol_type_list;
struct mcell_species_spec *mcell_mol_spec;
struct parse_mcell_species_list mcell_species_lst;
struct sm_dat *surf_mol_dat;
struct sm_dat_list surf_mol_dat_list;
struct species_list species_lst;
struct species_list_item *species_lst_item;

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


#line 380 "/home/jczech/mcell/build/deps/mdlparse.h" /* yacc.c:1909  */
};

typedef union YYSTYPE YYSTYPE;
# define YYSTYPE_IS_TRIVIAL 1
# define YYSTYPE_IS_DECLARED 1
#endif



int mdlparse (struct mdlparse_vars *parse_state, yyscan_t scanner);

#endif /* !YY_MDL_HOME_JCZECH_MCELL_BUILD_DEPS_MDLPARSE_H_INCLUDED  */
