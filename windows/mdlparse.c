/* A Bison parser, made by GNU Bison 3.0.4.  */

/* Bison implementation for Yacc-like parsers in C

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

/* C LALR(1) parser skeleton written by Richard Stallman, by
   simplifying the original so-called "semantic" parser.  */

/* All symbols defined below should begin with yy or YY, to avoid
   infringing on user name space.  This should be done even for local
   variables, as they might otherwise be expanded by user macros.
   There are some unavoidable exceptions within include files to
   define necessary library symbols; they are noted "INFRINGES ON
   USER NAME SPACE" below.  */

/* Identify Bison output.  */
#define YYBISON 1

/* Bison version.  */
#define YYBISON_VERSION "3.0.4"

/* Skeleton name.  */
#define YYSKELETON_NAME "yacc.c"

/* Pure parsers.  */
#define YYPURE 1

/* Push parsers.  */
#define YYPUSH 0

/* Pull parsers.  */
#define YYPULL 1


/* Substitute the variable and function names.  */
#define yyparse         mdlparse
#define yylex           mdllex
#define yyerror         mdlerror
#define yydebug         mdldebug
#define yynerrs         mdlnerrs


/* Copy the first part of user declarations.  */
#line 1 "../src/../src/mdlparse.y" /* yacc.c:339  */

  #include "config.h"
  #include <stdio.h>
  #include <stdlib.h>
  #include <stdarg.h>
  #include <string.h>
  #include <time.h>
  #include <math.h>
  #include <float.h>
  #include <limits.h>
  #include <errno.h>

  #include "rng.h"
  #include "logging.h"
  #include "vector.h"
  #include "strfunc.h"
  #include "mem_util.h"
  #include "sym_table.h"
  #include "diffuse_util.h"
  #include "mdlparse_util.h"
  #include "mdlparse_aux.h"
  #include "util.h"
  #include "react_output.h"
  #include "macromolecule.h"

  #include "mcell_misc.h"
  #include "mcell_structs.h"
  #include "mcell_viz.h"
  #include "mcell_release.h"
  #include "mcell_objects.h"

  /* make sure to declare yyscan_t before including mdlparse.h */
  typedef void *yyscan_t;
  #include "mdlparse.h"

  int mdllex_init(yyscan_t *ptr_yy_globals) ;
  int mdllex_destroy(yyscan_t yyscanner);
  void mdlrestart(FILE *infile, yyscan_t scanner);
  int mdllex(YYSTYPE *yylval, struct mdlparse_vars *parse_state, yyscan_t scanner);

  static int mdlparse_file(struct mdlparse_vars *parse_state, char const *name);


#ifdef DEBUG_MDL_PARSER
  #define FAILCHECK(t) do { mcell_error_nodie("Parser fail: %s:%d (%s)\n", __FILE__, __LINE__, t); return 1; } while(0)
#else
  #define FAILCHECK(t) return 1
#endif
#define CHECK(a)  do { if ((a) != 0) FAILCHECK("non-zero"); } while (0)
#define CHECKN(a) do { if ((a) == NULL) FAILCHECK("NULL"); } while (0)
#define CHECKF(a)  do {                                               \
                        if (isnan(a))                                 \
                        {                                             \
                          mdlerror(parse_state, "Expression result is not a number"); \
                          FAILCHECK("NaN");                           \
                        }                                             \
                        else if (isinf(a))                            \
                        {                                             \
                          mdlerror(parse_state, "Expression result is infinite"); \
                          FAILCHECK("Infinite");                      \
                        }                                             \
                      } while(0)

  #undef yyerror
  #define yyerror(a, b, c) mdlerror(a, c)

#line 139 "mdlparse.c" /* yacc.c:339  */

# ifndef YY_NULLPTR
#  if defined __cplusplus && 201103L <= __cplusplus
#   define YY_NULLPTR nullptr
#  else
#   define YY_NULLPTR 0
#  endif
# endif

/* Enabling verbose error messages.  */
#ifdef YYERROR_VERBOSE
# undef YYERROR_VERBOSE
# define YYERROR_VERBOSE 1
#else
# define YYERROR_VERBOSE 0
#endif

/* In a future release of Bison, this section will be replaced
   by #include "y.tab.h".  */
#ifndef YY_MDL_MDLPARSE_H_INCLUDED
# define YY_MDL_MDLPARSE_H_INCLUDED
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
    ALL_MESHES = 268,
    ALL_MOLECULES = 269,
    ALL_NOTIFICATIONS = 270,
    ALL_TIMES = 271,
    ALL_WARNINGS = 272,
    ASCII = 273,
    ASIN = 274,
    ASPECT_RATIO = 275,
    ATAN = 276,
    BACK = 277,
    BACK_CROSSINGS = 278,
    BACK_HITS = 279,
    BINARY = 280,
    BOTTOM = 281,
    BOX = 282,
    BOX_TRIANGULATION_REPORT = 283,
    BRIEF = 284,
    CEIL = 285,
    CELLBLENDER = 286,
    CENTER_MOLECULES_ON_GRID = 287,
    CHECKPOINT_INFILE = 288,
    CHECKPOINT_ITERATIONS = 289,
    CHECKPOINT_OUTFILE = 290,
    CHECKPOINT_REALTIME = 291,
    CHECKPOINT_REPORT = 292,
    CLAMP_CONCENTRATION = 293,
    CLOSE_PARTITION_SPACING = 294,
    COMPLEX_PLACEMENT_ATTEMPTS = 295,
    COMPLEX_PLACEMENT_FAILURE = 296,
    COMPLEX_PLACEMENT_FAILURE_THRESHOLD = 297,
    COMPLEX_RATE = 298,
    CONCENTRATION = 299,
    CORNERS = 300,
    COS = 301,
    COUNT = 302,
    CUBIC = 303,
    CUBIC_RELEASE_SITE = 304,
    CUSTOM_SPACE_STEP = 305,
    CUSTOM_TIME_STEP = 306,
    DEFAULT = 307,
    DEFINE_COMPLEX_MOLECULE = 308,
    DEFINE_MOLECULE = 309,
    DEFINE_MOLECULES = 310,
    DEFINE_REACTIONS = 311,
    DEFINE_RELEASE_PATTERN = 312,
    DEFINE_SURFACE_CLASS = 313,
    DEFINE_SURFACE_CLASSES = 314,
    DEFINE_SURFACE_REGIONS = 315,
    DEGENERATE_POLYGONS = 316,
    DELAY = 317,
    DENSITY = 318,
    DIFFUSION_CONSTANT_2D = 319,
    DIFFUSION_CONSTANT_3D = 320,
    DIFFUSION_CONSTANT_REPORT = 321,
    DREAMM_V3 = 322,
    DREAMM_V3_GROUPED = 323,
    EFFECTOR_GRID_DENSITY = 324,
    ELEMENT_CONNECTIONS = 325,
    ELLIPTIC = 326,
    ELLIPTIC_RELEASE_SITE = 327,
    EQUAL = 328,
    ERROR = 329,
    ESTIMATE_CONCENTRATION = 330,
    EXCLUDE_ELEMENTS = 331,
    EXCLUDE_PATCH = 332,
    EXCLUDE_REGION = 333,
    EXIT = 334,
    EXP = 335,
    EXPRESSION = 336,
    FALSE = 337,
    FCLOSE = 338,
    FILENAME = 339,
    FILENAME_PREFIX = 340,
    FILE_OUTPUT_REPORT = 341,
    FINAL_SUMMARY = 342,
    FLOOR = 343,
    FOPEN = 344,
    FORMAT = 345,
    FPRINTF = 346,
    FPRINT_TIME = 347,
    FRONT = 348,
    FRONT_CROSSINGS = 349,
    FRONT_HITS = 350,
    GAUSSIAN_RELEASE_NUMBER = 351,
    GEOMETRY = 352,
    HEADER = 353,
    HIGH_PROBABILITY_THRESHOLD = 354,
    HIGH_REACTION_PROBABILITY = 355,
    IGNORED = 356,
    INCLUDE_ELEMENTS = 357,
    INCLUDE_FILE = 358,
    INCLUDE_PATCH = 359,
    INCLUDE_REGION = 360,
    INPUT_FILE = 361,
    INSTANTIATE = 362,
    LLINTEGER = 363,
    FULLY_RANDOM = 364,
    INTERACTION_RADIUS = 365,
    ITERATION_LIST = 366,
    ITERATION_NUMBERS = 367,
    ITERATION_REPORT = 368,
    ITERATIONS = 369,
    KEEP_CHECKPOINT_FILES = 370,
    LEFT = 371,
    LIFETIME_THRESHOLD = 372,
    LIFETIME_TOO_SHORT = 373,
    LIST = 374,
    LOCATION = 375,
    LOG = 376,
    LOG10 = 377,
    MAX_TOK = 378,
    MAXIMUM_STEP_LENGTH = 379,
    MEAN_DIAMETER = 380,
    MEAN_NUMBER = 381,
    MEMORY_PARTITION_X = 382,
    MEMORY_PARTITION_Y = 383,
    MEMORY_PARTITION_Z = 384,
    MEMORY_PARTITION_POOL = 385,
    MESHES = 386,
    MICROSCOPIC_REVERSIBILITY = 387,
    MIN_TOK = 388,
    MISSED_REACTIONS = 389,
    MISSED_REACTION_THRESHOLD = 390,
    MISSING_SURFACE_ORIENTATION = 391,
    MOD = 392,
    MODE = 393,
    MODIFY_SURFACE_REGIONS = 394,
    MOLECULE = 395,
    MOLECULE_COLLISION_REPORT = 396,
    MOLECULE_DENSITY = 397,
    MOLECULE_NUMBER = 398,
    MOLECULE_POSITIONS = 399,
    MOLECULES = 400,
    MOLECULE_PLACEMENT_FAILURE = 401,
    NAME_LIST = 402,
    NEGATIVE_DIFFUSION_CONSTANT = 403,
    NEGATIVE_REACTION_RATE = 404,
    NO = 405,
    NOEXIT = 406,
    NONE = 407,
    NO_SPECIES = 408,
    NOT_EQUAL = 409,
    NOTIFICATIONS = 410,
    NUMBER_OF_SUBUNITS = 411,
    NUMBER_OF_TRAINS = 412,
    NUMBER_TO_RELEASE = 413,
    OBJECT = 414,
    OFF = 415,
    ON = 416,
    ORIENTATIONS = 417,
    OUTPUT_BUFFER_SIZE = 418,
    INVALID_OUTPUT_STEP_TIME = 419,
    OVERWRITTEN_OUTPUT_FILE = 420,
    PARTITION_LOCATION_REPORT = 421,
    PARTITION_X = 422,
    PARTITION_Y = 423,
    PARTITION_Z = 424,
    PI_TOK = 425,
    POLYGON_LIST = 426,
    POSITIONS = 427,
    PRINTF = 428,
    PRINT_TIME = 429,
    PROBABILITY_REPORT = 430,
    PROBABILITY_REPORT_THRESHOLD = 431,
    PROGRESS_REPORT = 432,
    RADIAL_DIRECTIONS = 433,
    RADIAL_SUBDIVISIONS = 434,
    RAND_GAUSSIAN = 435,
    RAND_UNIFORM = 436,
    RATE_RULES = 437,
    REACTION_DATA_OUTPUT = 438,
    REACTION_OUTPUT_REPORT = 439,
    REACTION_GROUP = 440,
    REAL = 441,
    RECTANGULAR_RELEASE_SITE = 442,
    RECTANGULAR_TOKEN = 443,
    REFLECTIVE = 444,
    REGION_DATA = 445,
    RELEASE_EVENT_REPORT = 446,
    RELEASE_INTERVAL = 447,
    RELEASE_PATTERN = 448,
    RELEASE_PROBABILITY = 449,
    RELEASE_SITE = 450,
    REMOVE_ELEMENTS = 451,
    RIGHT = 452,
    ROTATE = 453,
    ROUND_OFF = 454,
    SCALE = 455,
    SEED = 456,
    SHAPE = 457,
    SHOW_EXACT_TIME = 458,
    SIN = 459,
    SITE_DIAMETER = 460,
    SITE_RADIUS = 461,
    SPACE_STEP = 462,
    SPHERICAL = 463,
    SPHERICAL_RELEASE_SITE = 464,
    SPHERICAL_SHELL = 465,
    SPHERICAL_SHELL_SITE = 466,
    SPRINTF = 467,
    SQRT = 468,
    STANDARD_DEVIATION = 469,
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
    VIZ_MESH_FORMAT = 500,
    VIZ_MOLECULE_FORMAT = 501,
    VIZ_OUTPUT = 502,
    VIZ_OUTPUT_REPORT = 503,
    VIZ_VALUE = 504,
    VOLUME_DATA_OUTPUT = 505,
    VOLUME_OUTPUT_REPORT = 506,
    VOLUME_DEPENDENT_RELEASE_NUMBER = 507,
    VOLUME_ONLY = 508,
    VOXEL_COUNT = 509,
    VOXEL_LIST = 510,
    VOXEL_SIZE = 511,
    WARNING = 512,
    WARNINGS = 513,
    WORLD = 514,
    YES = 515,
    UNARYMINUS = 516
  };
#endif
/* Tokens.  */
#define ABS 258
#define ABSORPTIVE 259
#define ACCURATE_3D_REACTIONS 260
#define ACOS 261
#define ALL_CROSSINGS 262
#define ALL_DATA 263
#define ALL_ELEMENTS 264
#define ALL_ENCLOSED 265
#define ALL_HITS 266
#define ALL_ITERATIONS 267
#define ALL_MESHES 268
#define ALL_MOLECULES 269
#define ALL_NOTIFICATIONS 270
#define ALL_TIMES 271
#define ALL_WARNINGS 272
#define ASCII 273
#define ASIN 274
#define ASPECT_RATIO 275
#define ATAN 276
#define BACK 277
#define BACK_CROSSINGS 278
#define BACK_HITS 279
#define BINARY 280
#define BOTTOM 281
#define BOX 282
#define BOX_TRIANGULATION_REPORT 283
#define BRIEF 284
#define CEIL 285
#define CELLBLENDER 286
#define CENTER_MOLECULES_ON_GRID 287
#define CHECKPOINT_INFILE 288
#define CHECKPOINT_ITERATIONS 289
#define CHECKPOINT_OUTFILE 290
#define CHECKPOINT_REALTIME 291
#define CHECKPOINT_REPORT 292
#define CLAMP_CONCENTRATION 293
#define CLOSE_PARTITION_SPACING 294
#define COMPLEX_PLACEMENT_ATTEMPTS 295
#define COMPLEX_PLACEMENT_FAILURE 296
#define COMPLEX_PLACEMENT_FAILURE_THRESHOLD 297
#define COMPLEX_RATE 298
#define CONCENTRATION 299
#define CORNERS 300
#define COS 301
#define COUNT 302
#define CUBIC 303
#define CUBIC_RELEASE_SITE 304
#define CUSTOM_SPACE_STEP 305
#define CUSTOM_TIME_STEP 306
#define DEFAULT 307
#define DEFINE_COMPLEX_MOLECULE 308
#define DEFINE_MOLECULE 309
#define DEFINE_MOLECULES 310
#define DEFINE_REACTIONS 311
#define DEFINE_RELEASE_PATTERN 312
#define DEFINE_SURFACE_CLASS 313
#define DEFINE_SURFACE_CLASSES 314
#define DEFINE_SURFACE_REGIONS 315
#define DEGENERATE_POLYGONS 316
#define DELAY 317
#define DENSITY 318
#define DIFFUSION_CONSTANT_2D 319
#define DIFFUSION_CONSTANT_3D 320
#define DIFFUSION_CONSTANT_REPORT 321
#define DREAMM_V3 322
#define DREAMM_V3_GROUPED 323
#define EFFECTOR_GRID_DENSITY 324
#define ELEMENT_CONNECTIONS 325
#define ELLIPTIC 326
#define ELLIPTIC_RELEASE_SITE 327
#define EQUAL 328
#define ERROR 329
#define ESTIMATE_CONCENTRATION 330
#define EXCLUDE_ELEMENTS 331
#define EXCLUDE_PATCH 332
#define EXCLUDE_REGION 333
#define EXIT 334
#define EXP 335
#define EXPRESSION 336
#define FALSE 337
#define FCLOSE 338
#define FILENAME 339
#define FILENAME_PREFIX 340
#define FILE_OUTPUT_REPORT 341
#define FINAL_SUMMARY 342
#define FLOOR 343
#define FOPEN 344
#define FORMAT 345
#define FPRINTF 346
#define FPRINT_TIME 347
#define FRONT 348
#define FRONT_CROSSINGS 349
#define FRONT_HITS 350
#define GAUSSIAN_RELEASE_NUMBER 351
#define GEOMETRY 352
#define HEADER 353
#define HIGH_PROBABILITY_THRESHOLD 354
#define HIGH_REACTION_PROBABILITY 355
#define IGNORED 356
#define INCLUDE_ELEMENTS 357
#define INCLUDE_FILE 358
#define INCLUDE_PATCH 359
#define INCLUDE_REGION 360
#define INPUT_FILE 361
#define INSTANTIATE 362
#define LLINTEGER 363
#define FULLY_RANDOM 364
#define INTERACTION_RADIUS 365
#define ITERATION_LIST 366
#define ITERATION_NUMBERS 367
#define ITERATION_REPORT 368
#define ITERATIONS 369
#define KEEP_CHECKPOINT_FILES 370
#define LEFT 371
#define LIFETIME_THRESHOLD 372
#define LIFETIME_TOO_SHORT 373
#define LIST 374
#define LOCATION 375
#define LOG 376
#define LOG10 377
#define MAX_TOK 378
#define MAXIMUM_STEP_LENGTH 379
#define MEAN_DIAMETER 380
#define MEAN_NUMBER 381
#define MEMORY_PARTITION_X 382
#define MEMORY_PARTITION_Y 383
#define MEMORY_PARTITION_Z 384
#define MEMORY_PARTITION_POOL 385
#define MESHES 386
#define MICROSCOPIC_REVERSIBILITY 387
#define MIN_TOK 388
#define MISSED_REACTIONS 389
#define MISSED_REACTION_THRESHOLD 390
#define MISSING_SURFACE_ORIENTATION 391
#define MOD 392
#define MODE 393
#define MODIFY_SURFACE_REGIONS 394
#define MOLECULE 395
#define MOLECULE_COLLISION_REPORT 396
#define MOLECULE_DENSITY 397
#define MOLECULE_NUMBER 398
#define MOLECULE_POSITIONS 399
#define MOLECULES 400
#define MOLECULE_PLACEMENT_FAILURE 401
#define NAME_LIST 402
#define NEGATIVE_DIFFUSION_CONSTANT 403
#define NEGATIVE_REACTION_RATE 404
#define NO 405
#define NOEXIT 406
#define NONE 407
#define NO_SPECIES 408
#define NOT_EQUAL 409
#define NOTIFICATIONS 410
#define NUMBER_OF_SUBUNITS 411
#define NUMBER_OF_TRAINS 412
#define NUMBER_TO_RELEASE 413
#define OBJECT 414
#define OFF 415
#define ON 416
#define ORIENTATIONS 417
#define OUTPUT_BUFFER_SIZE 418
#define INVALID_OUTPUT_STEP_TIME 419
#define OVERWRITTEN_OUTPUT_FILE 420
#define PARTITION_LOCATION_REPORT 421
#define PARTITION_X 422
#define PARTITION_Y 423
#define PARTITION_Z 424
#define PI_TOK 425
#define POLYGON_LIST 426
#define POSITIONS 427
#define PRINTF 428
#define PRINT_TIME 429
#define PROBABILITY_REPORT 430
#define PROBABILITY_REPORT_THRESHOLD 431
#define PROGRESS_REPORT 432
#define RADIAL_DIRECTIONS 433
#define RADIAL_SUBDIVISIONS 434
#define RAND_GAUSSIAN 435
#define RAND_UNIFORM 436
#define RATE_RULES 437
#define REACTION_DATA_OUTPUT 438
#define REACTION_OUTPUT_REPORT 439
#define REACTION_GROUP 440
#define REAL 441
#define RECTANGULAR_RELEASE_SITE 442
#define RECTANGULAR_TOKEN 443
#define REFLECTIVE 444
#define REGION_DATA 445
#define RELEASE_EVENT_REPORT 446
#define RELEASE_INTERVAL 447
#define RELEASE_PATTERN 448
#define RELEASE_PROBABILITY 449
#define RELEASE_SITE 450
#define REMOVE_ELEMENTS 451
#define RIGHT 452
#define ROTATE 453
#define ROUND_OFF 454
#define SCALE 455
#define SEED 456
#define SHAPE 457
#define SHOW_EXACT_TIME 458
#define SIN 459
#define SITE_DIAMETER 460
#define SITE_RADIUS 461
#define SPACE_STEP 462
#define SPHERICAL 463
#define SPHERICAL_RELEASE_SITE 464
#define SPHERICAL_SHELL 465
#define SPHERICAL_SHELL_SITE 466
#define SPRINTF 467
#define SQRT 468
#define STANDARD_DEVIATION 469
#define STEP 470
#define STRING_TO_NUM 471
#define STR_VALUE 472
#define SUBUNIT 473
#define SUBUNIT_RELATIONSHIPS 474
#define SUMMATION_OPERATOR 475
#define SURFACE_CLASS 476
#define SURFACE_ONLY 477
#define TAN 478
#define TARGET_ONLY 479
#define TET_ELEMENT_CONNECTIONS 480
#define THROUGHPUT_REPORT 481
#define TIME_LIST 482
#define TIME_POINTS 483
#define TIME_STEP 484
#define TIME_STEP_MAX 485
#define TO 486
#define TOP 487
#define TRAIN_DURATION 488
#define TRAIN_INTERVAL 489
#define TRANSLATE 490
#define TRANSPARENT 491
#define TRIGGER 492
#define TRUE 493
#define UNLIMITED 494
#define USELESS_VOLUME_ORIENTATION 495
#define VACANCY_SEARCH_DISTANCE 496
#define VAR 497
#define VARYING_PROBABILITY_REPORT 498
#define VERTEX_LIST 499
#define VIZ_MESH_FORMAT 500
#define VIZ_MOLECULE_FORMAT 501
#define VIZ_OUTPUT 502
#define VIZ_OUTPUT_REPORT 503
#define VIZ_VALUE 504
#define VOLUME_DATA_OUTPUT 505
#define VOLUME_OUTPUT_REPORT 506
#define VOLUME_DEPENDENT_RELEASE_NUMBER 507
#define VOLUME_ONLY 508
#define VOXEL_COUNT 509
#define VOXEL_LIST 510
#define VOXEL_SIZE 511
#define WARNING 512
#define WARNINGS 513
#define WORLD 514
#define YES 515
#define UNARYMINUS 516

/* Value type.  */
#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED

union YYSTYPE
{
#line 69 "../src/../src/mdlparse.y" /* yacc.c:355  */

int ival;
int tok;
double dbl;
long long llival;
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

/* Macromolecules */
struct macro_subunit_assignment_list mmol_subunits;
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

#line 780 "mdlparse.c" /* yacc.c:355  */
};

typedef union YYSTYPE YYSTYPE;
# define YYSTYPE_IS_TRIVIAL 1
# define YYSTYPE_IS_DECLARED 1
#endif



int mdlparse (struct mdlparse_vars *parse_state, yyscan_t scanner);

#endif /* !YY_MDL_MDLPARSE_H_INCLUDED  */

/* Copy the second part of user declarations.  */

#line 796 "mdlparse.c" /* yacc.c:358  */

#ifdef short
# undef short
#endif

#ifdef YYTYPE_UINT8
typedef YYTYPE_UINT8 yytype_uint8;
#else
typedef unsigned char yytype_uint8;
#endif

#ifdef YYTYPE_INT8
typedef YYTYPE_INT8 yytype_int8;
#else
typedef signed char yytype_int8;
#endif

#ifdef YYTYPE_UINT16
typedef YYTYPE_UINT16 yytype_uint16;
#else
typedef unsigned short int yytype_uint16;
#endif

#ifdef YYTYPE_INT16
typedef YYTYPE_INT16 yytype_int16;
#else
typedef short int yytype_int16;
#endif

#ifndef YYSIZE_T
# ifdef __SIZE_TYPE__
#  define YYSIZE_T __SIZE_TYPE__
# elif defined size_t
#  define YYSIZE_T size_t
# elif ! defined YYSIZE_T
#  include <stddef.h> /* INFRINGES ON USER NAME SPACE */
#  define YYSIZE_T size_t
# else
#  define YYSIZE_T unsigned int
# endif
#endif

#define YYSIZE_MAXIMUM ((YYSIZE_T) -1)

#ifndef YY_
# if defined YYENABLE_NLS && YYENABLE_NLS
#  if ENABLE_NLS
#   include <libintl.h> /* INFRINGES ON USER NAME SPACE */
#   define YY_(Msgid) dgettext ("bison-runtime", Msgid)
#  endif
# endif
# ifndef YY_
#  define YY_(Msgid) Msgid
# endif
#endif

#ifndef YY_ATTRIBUTE
# if (defined __GNUC__                                               \
      && (2 < __GNUC__ || (__GNUC__ == 2 && 96 <= __GNUC_MINOR__)))  \
     || defined __SUNPRO_C && 0x5110 <= __SUNPRO_C
#  define YY_ATTRIBUTE(Spec) __attribute__(Spec)
# else
#  define YY_ATTRIBUTE(Spec) /* empty */
# endif
#endif

#ifndef YY_ATTRIBUTE_PURE
# define YY_ATTRIBUTE_PURE   YY_ATTRIBUTE ((__pure__))
#endif

#ifndef YY_ATTRIBUTE_UNUSED
# define YY_ATTRIBUTE_UNUSED YY_ATTRIBUTE ((__unused__))
#endif

#if !defined _Noreturn \
     && (!defined __STDC_VERSION__ || __STDC_VERSION__ < 201112)
# if defined _MSC_VER && 1200 <= _MSC_VER
#  define _Noreturn __declspec (noreturn)
# else
#  define _Noreturn YY_ATTRIBUTE ((__noreturn__))
# endif
#endif

/* Suppress unused-variable warnings by "using" E.  */
#if ! defined lint || defined __GNUC__
# define YYUSE(E) ((void) (E))
#else
# define YYUSE(E) /* empty */
#endif

#if defined __GNUC__ && 407 <= __GNUC__ * 100 + __GNUC_MINOR__
/* Suppress an incorrect diagnostic about yylval being uninitialized.  */
# define YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN \
    _Pragma ("GCC diagnostic push") \
    _Pragma ("GCC diagnostic ignored \"-Wuninitialized\"")\
    _Pragma ("GCC diagnostic ignored \"-Wmaybe-uninitialized\"")
# define YY_IGNORE_MAYBE_UNINITIALIZED_END \
    _Pragma ("GCC diagnostic pop")
#else
# define YY_INITIAL_VALUE(Value) Value
#endif
#ifndef YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
# define YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
# define YY_IGNORE_MAYBE_UNINITIALIZED_END
#endif
#ifndef YY_INITIAL_VALUE
# define YY_INITIAL_VALUE(Value) /* Nothing. */
#endif


#if ! defined yyoverflow || YYERROR_VERBOSE

/* The parser invokes alloca or malloc; define the necessary symbols.  */

# ifdef YYSTACK_USE_ALLOCA
#  if YYSTACK_USE_ALLOCA
#   ifdef __GNUC__
#    define YYSTACK_ALLOC __builtin_alloca
#   elif defined __BUILTIN_VA_ARG_INCR
#    include <alloca.h> /* INFRINGES ON USER NAME SPACE */
#   elif defined _AIX
#    define YYSTACK_ALLOC __alloca
#   elif defined _MSC_VER
#    include <malloc.h> /* INFRINGES ON USER NAME SPACE */
#    define alloca _alloca
#   else
#    define YYSTACK_ALLOC alloca
#    if ! defined _ALLOCA_H && ! defined EXIT_SUCCESS
#     include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
      /* Use EXIT_SUCCESS as a witness for stdlib.h.  */
#     ifndef EXIT_SUCCESS
#      define EXIT_SUCCESS 0
#     endif
#    endif
#   endif
#  endif
# endif

# ifdef YYSTACK_ALLOC
   /* Pacify GCC's 'empty if-body' warning.  */
#  define YYSTACK_FREE(Ptr) do { /* empty */; } while (0)
#  ifndef YYSTACK_ALLOC_MAXIMUM
    /* The OS might guarantee only one guard page at the bottom of the stack,
       and a page size can be as small as 4096 bytes.  So we cannot safely
       invoke alloca (N) if N exceeds 4096.  Use a slightly smaller number
       to allow for a few compiler-allocated temporary stack slots.  */
#   define YYSTACK_ALLOC_MAXIMUM 4032 /* reasonable circa 2006 */
#  endif
# else
#  define YYSTACK_ALLOC YYMALLOC
#  define YYSTACK_FREE YYFREE
#  ifndef YYSTACK_ALLOC_MAXIMUM
#   define YYSTACK_ALLOC_MAXIMUM YYSIZE_MAXIMUM
#  endif
#  if (defined __cplusplus && ! defined EXIT_SUCCESS \
       && ! ((defined YYMALLOC || defined malloc) \
             && (defined YYFREE || defined free)))
#   include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
#   ifndef EXIT_SUCCESS
#    define EXIT_SUCCESS 0
#   endif
#  endif
#  ifndef YYMALLOC
#   define YYMALLOC malloc
#   if ! defined malloc && ! defined EXIT_SUCCESS
void *malloc (YYSIZE_T); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
#  ifndef YYFREE
#   define YYFREE free
#   if ! defined free && ! defined EXIT_SUCCESS
void free (void *); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
# endif
#endif /* ! defined yyoverflow || YYERROR_VERBOSE */


#if (! defined yyoverflow \
     && (! defined __cplusplus \
         || (defined YYSTYPE_IS_TRIVIAL && YYSTYPE_IS_TRIVIAL)))

/* A type that is properly aligned for any stack member.  */
union yyalloc
{
  yytype_int16 yyss_alloc;
  YYSTYPE yyvs_alloc;
};

/* The size of the maximum gap between one aligned stack and the next.  */
# define YYSTACK_GAP_MAXIMUM (sizeof (union yyalloc) - 1)

/* The size of an array large to enough to hold all stacks, each with
   N elements.  */
# define YYSTACK_BYTES(N) \
     ((N) * (sizeof (yytype_int16) + sizeof (YYSTYPE)) \
      + YYSTACK_GAP_MAXIMUM)

# define YYCOPY_NEEDED 1

/* Relocate STACK from its old location to the new one.  The
   local variables YYSIZE and YYSTACKSIZE give the old and new number of
   elements in the stack, and YYPTR gives the new location of the
   stack.  Advance YYPTR to a properly aligned location for the next
   stack.  */
# define YYSTACK_RELOCATE(Stack_alloc, Stack)                           \
    do                                                                  \
      {                                                                 \
        YYSIZE_T yynewbytes;                                            \
        YYCOPY (&yyptr->Stack_alloc, Stack, yysize);                    \
        Stack = &yyptr->Stack_alloc;                                    \
        yynewbytes = yystacksize * sizeof (*Stack) + YYSTACK_GAP_MAXIMUM; \
        yyptr += yynewbytes / sizeof (*yyptr);                          \
      }                                                                 \
    while (0)

#endif

#if defined YYCOPY_NEEDED && YYCOPY_NEEDED
/* Copy COUNT objects from SRC to DST.  The source and destination do
   not overlap.  */
# ifndef YYCOPY
#  if defined __GNUC__ && 1 < __GNUC__
#   define YYCOPY(Dst, Src, Count) \
      __builtin_memcpy (Dst, Src, (Count) * sizeof (*(Src)))
#  else
#   define YYCOPY(Dst, Src, Count)              \
      do                                        \
        {                                       \
          YYSIZE_T yyi;                         \
          for (yyi = 0; yyi < (Count); yyi++)   \
            (Dst)[yyi] = (Src)[yyi];            \
        }                                       \
      while (0)
#  endif
# endif
#endif /* !YYCOPY_NEEDED */

/* YYFINAL -- State number of the termination state.  */
#define YYFINAL  148
/* YYLAST -- Last index in YYTABLE.  */
#define YYLAST   2861

/* YYNTOKENS -- Number of terminals.  */
#define YYNTOKENS  282
/* YYNNTS -- Number of nonterminals.  */
#define YYNNTS  338
/* YYNRULES -- Number of rules.  */
#define YYNRULES  697
/* YYNSTATES -- Number of states.  */
#define YYNSTATES  1366

/* YYTRANSLATE[YYX] -- Symbol number corresponding to YYX as returned
   by yylex, with out-of-bounds checking.  */
#define YYUNDEFTOK  2
#define YYMAXUTOK   516

#define YYTRANSLATE(YYX)                                                \
  ((unsigned int) (YYX) <= YYMAXUTOK ? yytranslate[YYX] : YYUNDEFTOK)

/* YYTRANSLATE[TOKEN-NUM] -- Symbol number corresponding to TOKEN-NUM
   as returned by yylex, without out-of-bounds checking.  */
static const yytype_uint16 yytranslate[] =
{
       0,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,   262,   273,
     277,   278,   266,   264,   274,   265,     2,   267,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,   263,   272,
     280,   261,   279,     2,   281,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,   270,     2,   271,   268,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,   275,     2,   276,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     1,     2,     3,     4,
       5,     6,     7,     8,     9,    10,    11,    12,    13,    14,
      15,    16,    17,    18,    19,    20,    21,    22,    23,    24,
      25,    26,    27,    28,    29,    30,    31,    32,    33,    34,
      35,    36,    37,    38,    39,    40,    41,    42,    43,    44,
      45,    46,    47,    48,    49,    50,    51,    52,    53,    54,
      55,    56,    57,    58,    59,    60,    61,    62,    63,    64,
      65,    66,    67,    68,    69,    70,    71,    72,    73,    74,
      75,    76,    77,    78,    79,    80,    81,    82,    83,    84,
      85,    86,    87,    88,    89,    90,    91,    92,    93,    94,
      95,    96,    97,    98,    99,   100,   101,   102,   103,   104,
     105,   106,   107,   108,   109,   110,   111,   112,   113,   114,
     115,   116,   117,   118,   119,   120,   121,   122,   123,   124,
     125,   126,   127,   128,   129,   130,   131,   132,   133,   134,
     135,   136,   137,   138,   139,   140,   141,   142,   143,   144,
     145,   146,   147,   148,   149,   150,   151,   152,   153,   154,
     155,   156,   157,   158,   159,   160,   161,   162,   163,   164,
     165,   166,   167,   168,   169,   170,   171,   172,   173,   174,
     175,   176,   177,   178,   179,   180,   181,   182,   183,   184,
     185,   186,   187,   188,   189,   190,   191,   192,   193,   194,
     195,   196,   197,   198,   199,   200,   201,   202,   203,   204,
     205,   206,   207,   208,   209,   210,   211,   212,   213,   214,
     215,   216,   217,   218,   219,   220,   221,   222,   223,   224,
     225,   226,   227,   228,   229,   230,   231,   232,   233,   234,
     235,   236,   237,   238,   239,   240,   241,   242,   243,   244,
     245,   246,   247,   248,   249,   250,   251,   252,   253,   254,
     255,   256,   257,   258,   259,   260,   269
};

#if YYDEBUG
  /* YYRLINE[YYN] -- Source line where rule number YYN was defined.  */
static const yytype_uint16 yyrline[] =
{
       0,   650,   650,   654,   655,   660,   661,   662,   663,   664,
     665,   666,   667,   668,   669,   670,   671,   672,   673,   674,
     675,   676,   677,   678,   679,   684,   687,   690,   693,   696,
     697,   700,   703,   706,   707,   710,   711,   712,   713,   714,
     715,   718,   719,   720,   721,   725,   726,   727,   737,   749,
     752,   755,   767,   768,   781,   782,   788,   811,   812,   813,
     814,   817,   820,   823,   824,   835,   838,   841,   842,   845,
     846,   849,   850,   853,   854,   857,   861,   862,   863,   864,
     865,   866,   867,   868,   869,   870,   871,   872,   873,   874,
     875,   876,   877,   878,   879,   880,   881,   882,   883,   884,
     885,   886,   887,   888,   889,   890,   894,   895,   899,   900,
     901,   902,   905,   911,   912,   913,   914,   915,   916,   917,
     920,   924,   927,   930,   933,   936,   939,   940,   950,   951,
     952,   964,   968,   974,   979,   983,   992,   996,   997,  1001,
    1002,  1003,  1004,  1005,  1006,  1007,  1008,  1009,  1010,  1011,
    1012,  1013,  1014,  1015,  1016,  1017,  1021,  1022,  1026,  1030,
    1031,  1038,  1042,  1043,  1047,  1048,  1049,  1050,  1051,  1052,
    1053,  1054,  1055,  1056,  1057,  1058,  1059,  1060,  1061,  1062,
    1063,  1064,  1068,  1069,  1070,  1076,  1077,  1078,  1079,  1080,
    1084,  1085,  1086,  1090,  1091,  1092,  1093,  1101,  1102,  1103,
    1104,  1105,  1106,  1107,  1108,  1109,  1110,  1111,  1112,  1113,
    1114,  1115,  1116,  1123,  1124,  1125,  1126,  1130,  1134,  1135,
    1136,  1143,  1144,  1145,  1148,  1152,  1156,  1157,  1161,  1169,
    1172,  1176,  1177,  1181,  1182,  1191,  1202,  1203,  1207,  1208,
    1219,  1221,  1224,  1219,  1229,  1233,  1237,  1238,  1243,  1248,
    1249,  1253,  1254,  1258,  1262,  1263,  1268,  1269,  1273,  1277,
    1282,  1283,  1287,  1291,  1295,  1296,  1299,  1303,  1304,  1309,
    1313,  1314,  1318,  1319,  1324,  1327,  1328,  1331,  1335,  1339,
    1347,  1354,  1355,  1360,  1364,  1370,  1371,  1376,  1376,  1381,
    1384,  1386,  1391,  1392,  1393,  1397,  1400,  1406,  1411,  1412,
    1413,  1416,  1417,  1420,  1424,  1428,  1435,  1439,  1448,  1452,
    1461,  1465,  1472,  1477,  1478,  1481,  1482,  1486,  1490,  1493,
    1494,  1497,  1498,  1501,  1502,  1503,  1506,  1511,  1517,  1518,
    1519,  1520,  1523,  1524,  1528,  1533,  1534,  1537,  1541,  1542,
    1546,  1547,  1551,  1554,  1555,  1558,  1559,  1563,  1564,  1567,
    1578,  1595,  1596,  1597,  1601,  1602,  1603,  1604,  1611,  1618,
    1621,  1625,  1632,  1634,  1636,  1638,  1640,  1644,  1645,  1652,
    1652,  1663,  1666,  1667,  1668,  1669,  1670,  1679,  1682,  1685,
    1688,  1690,  1694,  1698,  1699,  1700,  1705,  1718,  1719,  1722,
    1723,  1728,  1727,  1734,  1735,  1740,  1739,  1747,  1748,  1749,
    1750,  1751,  1752,  1753,  1754,  1761,  1762,  1763,  1764,  1765,
    1770,  1769,  1776,  1777,  1778,  1779,  1780,  1784,  1785,  1788,
    1792,  1793,  1794,  1801,  1802,  1803,  1804,  1805,  1807,  1812,
    1813,  1817,  1818,  1819,  1820,  1825,  1826,  1832,  1839,  1847,
    1848,  1852,  1853,  1858,  1861,  1869,  1866,  1884,  1887,  1890,
    1891,  1895,  1900,  1901,  1905,  1908,  1910,  1916,  1917,  1921,
    1921,  1933,  1934,  1937,  1938,  1939,  1940,  1941,  1942,  1943,
    1947,  1948,  1953,  1954,  1955,  1956,  1960,  1965,  1969,  1973,
    1974,  1977,  1978,  1979,  1982,  1985,  1986,  1989,  1992,  1993,
    1997,  2003,  2004,  2009,  2010,  2009,  2020,  2017,  2030,  2034,
    2038,  2042,  2053,  2054,  2050,  2062,  2063,  2077,  2083,  2084,
    2089,  2090,  2092,  2089,  2101,  2104,  2106,  2111,  2112,  2113,
    2117,  2121,  2128,  2134,  2135,  2140,  2140,  2150,  2149,  2160,
    2161,  2172,  2173,  2174,  2177,  2181,  2189,  2196,  2197,  2211,
    2212,  2213,  2217,  2217,  2223,  2224,  2225,  2229,  2233,  2237,
    2238,  2247,  2251,  2252,  2253,  2254,  2255,  2256,  2257,  2258,
    2259,  2264,  2264,  2266,  2267,  2267,  2271,  2272,  2273,  2274,
    2275,  2278,  2281,  2285,  2295,  2296,  2297,  2298,  2302,  2307,
    2312,  2317,  2320,  2331,  2330,  2347,  2348,  2352,  2353,  2358,
    2363,  2378,  2379,  2380,  2383,  2384,  2387,  2388,  2389,  2390,
    2391,  2392,  2393,  2394,  2397,  2398,  2405,  2405,  2414,  2415,
    2419,  2420,  2423,  2424,  2425,  2426,  2427,  2430,  2434,  2437,
    2438,  2442,  2446,  2450,  2451,  2455,  2456,  2466,  2467,  2470,
    2474,  2480,  2481,  2496,  2497,  2498,  2502,  2508,  2509,  2513,
    2514,  2518,  2520,  2524,  2525,  2529,  2530,  2533,  2539,  2540,
    2557,  2562,  2563,  2567,  2573,  2574,  2591,  2595,  2596,  2597,
    2601,  2607,  2608,  2623,  2624,  2625,  2629,  2635,  2636,  2640,
    2641,  2642,  2646,  2652,  2653,  2670,  2675,  2681,  2682,  2699,
    2703,  2704,  2705,  2712,  2728,  2732,  2733,  2743,  2746,  2764,
    2765,  2774,  2778,  2782,  2803,  2804,  2805,  2806
};
#endif

#if YYDEBUG || YYERROR_VERBOSE || 0
/* YYTNAME[SYMBOL-NUM] -- String name of the symbol SYMBOL-NUM.
   First, the terminals, then, starting at YYNTOKENS, nonterminals.  */
static const char *const yytname[] =
{
  "$end", "error", "$undefined", "ABS", "ABSORPTIVE",
  "ACCURATE_3D_REACTIONS", "ACOS", "ALL_CROSSINGS", "ALL_DATA",
  "ALL_ELEMENTS", "ALL_ENCLOSED", "ALL_HITS", "ALL_ITERATIONS",
  "ALL_MESHES", "ALL_MOLECULES", "ALL_NOTIFICATIONS", "ALL_TIMES",
  "ALL_WARNINGS", "ASCII", "ASIN", "ASPECT_RATIO", "ATAN", "BACK",
  "BACK_CROSSINGS", "BACK_HITS", "BINARY", "BOTTOM", "BOX",
  "BOX_TRIANGULATION_REPORT", "BRIEF", "CEIL", "CELLBLENDER",
  "CENTER_MOLECULES_ON_GRID", "CHECKPOINT_INFILE", "CHECKPOINT_ITERATIONS",
  "CHECKPOINT_OUTFILE", "CHECKPOINT_REALTIME", "CHECKPOINT_REPORT",
  "CLAMP_CONCENTRATION", "CLOSE_PARTITION_SPACING",
  "COMPLEX_PLACEMENT_ATTEMPTS", "COMPLEX_PLACEMENT_FAILURE",
  "COMPLEX_PLACEMENT_FAILURE_THRESHOLD", "COMPLEX_RATE", "CONCENTRATION",
  "CORNERS", "COS", "COUNT", "CUBIC", "CUBIC_RELEASE_SITE",
  "CUSTOM_SPACE_STEP", "CUSTOM_TIME_STEP", "DEFAULT",
  "DEFINE_COMPLEX_MOLECULE", "DEFINE_MOLECULE", "DEFINE_MOLECULES",
  "DEFINE_REACTIONS", "DEFINE_RELEASE_PATTERN", "DEFINE_SURFACE_CLASS",
  "DEFINE_SURFACE_CLASSES", "DEFINE_SURFACE_REGIONS",
  "DEGENERATE_POLYGONS", "DELAY", "DENSITY", "DIFFUSION_CONSTANT_2D",
  "DIFFUSION_CONSTANT_3D", "DIFFUSION_CONSTANT_REPORT", "DREAMM_V3",
  "DREAMM_V3_GROUPED", "EFFECTOR_GRID_DENSITY", "ELEMENT_CONNECTIONS",
  "ELLIPTIC", "ELLIPTIC_RELEASE_SITE", "EQUAL", "ERROR",
  "ESTIMATE_CONCENTRATION", "EXCLUDE_ELEMENTS", "EXCLUDE_PATCH",
  "EXCLUDE_REGION", "EXIT", "EXP", "EXPRESSION", "FALSE", "FCLOSE",
  "FILENAME", "FILENAME_PREFIX", "FILE_OUTPUT_REPORT", "FINAL_SUMMARY",
  "FLOOR", "FOPEN", "FORMAT", "FPRINTF", "FPRINT_TIME", "FRONT",
  "FRONT_CROSSINGS", "FRONT_HITS", "GAUSSIAN_RELEASE_NUMBER", "GEOMETRY",
  "HEADER", "HIGH_PROBABILITY_THRESHOLD", "HIGH_REACTION_PROBABILITY",
  "IGNORED", "INCLUDE_ELEMENTS", "INCLUDE_FILE", "INCLUDE_PATCH",
  "INCLUDE_REGION", "INPUT_FILE", "INSTANTIATE", "LLINTEGER",
  "FULLY_RANDOM", "INTERACTION_RADIUS", "ITERATION_LIST",
  "ITERATION_NUMBERS", "ITERATION_REPORT", "ITERATIONS",
  "KEEP_CHECKPOINT_FILES", "LEFT", "LIFETIME_THRESHOLD",
  "LIFETIME_TOO_SHORT", "LIST", "LOCATION", "LOG", "LOG10", "MAX_TOK",
  "MAXIMUM_STEP_LENGTH", "MEAN_DIAMETER", "MEAN_NUMBER",
  "MEMORY_PARTITION_X", "MEMORY_PARTITION_Y", "MEMORY_PARTITION_Z",
  "MEMORY_PARTITION_POOL", "MESHES", "MICROSCOPIC_REVERSIBILITY",
  "MIN_TOK", "MISSED_REACTIONS", "MISSED_REACTION_THRESHOLD",
  "MISSING_SURFACE_ORIENTATION", "MOD", "MODE", "MODIFY_SURFACE_REGIONS",
  "MOLECULE", "MOLECULE_COLLISION_REPORT", "MOLECULE_DENSITY",
  "MOLECULE_NUMBER", "MOLECULE_POSITIONS", "MOLECULES",
  "MOLECULE_PLACEMENT_FAILURE", "NAME_LIST", "NEGATIVE_DIFFUSION_CONSTANT",
  "NEGATIVE_REACTION_RATE", "NO", "NOEXIT", "NONE", "NO_SPECIES",
  "NOT_EQUAL", "NOTIFICATIONS", "NUMBER_OF_SUBUNITS", "NUMBER_OF_TRAINS",
  "NUMBER_TO_RELEASE", "OBJECT", "OFF", "ON", "ORIENTATIONS",
  "OUTPUT_BUFFER_SIZE", "INVALID_OUTPUT_STEP_TIME",
  "OVERWRITTEN_OUTPUT_FILE", "PARTITION_LOCATION_REPORT", "PARTITION_X",
  "PARTITION_Y", "PARTITION_Z", "PI_TOK", "POLYGON_LIST", "POSITIONS",
  "PRINTF", "PRINT_TIME", "PROBABILITY_REPORT",
  "PROBABILITY_REPORT_THRESHOLD", "PROGRESS_REPORT", "RADIAL_DIRECTIONS",
  "RADIAL_SUBDIVISIONS", "RAND_GAUSSIAN", "RAND_UNIFORM", "RATE_RULES",
  "REACTION_DATA_OUTPUT", "REACTION_OUTPUT_REPORT", "REACTION_GROUP",
  "REAL", "RECTANGULAR_RELEASE_SITE", "RECTANGULAR_TOKEN", "REFLECTIVE",
  "REGION_DATA", "RELEASE_EVENT_REPORT", "RELEASE_INTERVAL",
  "RELEASE_PATTERN", "RELEASE_PROBABILITY", "RELEASE_SITE",
  "REMOVE_ELEMENTS", "RIGHT", "ROTATE", "ROUND_OFF", "SCALE", "SEED",
  "SHAPE", "SHOW_EXACT_TIME", "SIN", "SITE_DIAMETER", "SITE_RADIUS",
  "SPACE_STEP", "SPHERICAL", "SPHERICAL_RELEASE_SITE", "SPHERICAL_SHELL",
  "SPHERICAL_SHELL_SITE", "SPRINTF", "SQRT", "STANDARD_DEVIATION", "STEP",
  "STRING_TO_NUM", "STR_VALUE", "SUBUNIT", "SUBUNIT_RELATIONSHIPS",
  "SUMMATION_OPERATOR", "SURFACE_CLASS", "SURFACE_ONLY", "TAN",
  "TARGET_ONLY", "TET_ELEMENT_CONNECTIONS", "THROUGHPUT_REPORT",
  "TIME_LIST", "TIME_POINTS", "TIME_STEP", "TIME_STEP_MAX", "TO", "TOP",
  "TRAIN_DURATION", "TRAIN_INTERVAL", "TRANSLATE", "TRANSPARENT",
  "TRIGGER", "TRUE", "UNLIMITED", "USELESS_VOLUME_ORIENTATION",
  "VACANCY_SEARCH_DISTANCE", "VAR", "VARYING_PROBABILITY_REPORT",
  "VERTEX_LIST", "VIZ_MESH_FORMAT", "VIZ_MOLECULE_FORMAT", "VIZ_OUTPUT",
  "VIZ_OUTPUT_REPORT", "VIZ_VALUE", "VOLUME_DATA_OUTPUT",
  "VOLUME_OUTPUT_REPORT", "VOLUME_DEPENDENT_RELEASE_NUMBER", "VOLUME_ONLY",
  "VOXEL_COUNT", "VOXEL_LIST", "VOXEL_SIZE", "WARNING", "WARNINGS",
  "WORLD", "YES", "'='", "'&'", "':'", "'+'", "'-'", "'*'", "'/'", "'^'",
  "UNARYMINUS", "'['", "']'", "';'", "'\\''", "','", "'{'", "'}'", "'('",
  "')'", "'>'", "'<'", "'@'", "$accept", "mdl_format", "mdl_stmt_list",
  "mdl_stmt", "str_value", "var", "file_name", "existing_object",
  "mesh_object_or_wildcard", "existing_region", "point", "point_or_num",
  "boolean", "orientation_class", "list_orient_marks", "head_mark",
  "tail_mark", "orient_class_number", "list_range_specs", "range_spec",
  "include_stmt", "assignment_stmt", "assign_var", "existing_var_only",
  "array_value", "array_expr_only", "existing_array", "num_expr",
  "num_value", "intOrReal", "num_expr_only", "existing_num_var",
  "arith_expr", "str_expr", "str_expr_only", "existing_str_var", "io_stmt",
  "fopen_stmt", "new_file_stream", "file_mode", "fclose_stmt",
  "existing_file_stream", "format_string", "list_args", "list_arg",
  "printf_stmt", "fprintf_stmt", "sprintf_stmt", "print_time_stmt",
  "fprint_time_stmt", "notification_def", "notification_list",
  "notification_item_def", "notify_bilevel", "notify_level",
  "warnings_def", "warning_list", "warning_item_def", "warning_level",
  "chkpt_stmt", "exit_or_no", "time_expr", "parameter_def",
  "memory_partition_def", "partition_def", "partition_dimension",
  "molecules_def", "define_one_molecule", "define_multiple_molecules",
  "list_molecule_stmts", "molecule_stmt", "molecule_name", "new_molecule",
  "diffusion_def", "mol_timestep_def", "target_def",
  "maximum_step_length_def", "define_complex_molecule", "$@1", "$@2",
  "$@3", "complex_mol_name", "complex_mol_topology",
  "complex_mol_subunits", "complex_mol_subunit_assignment",
  "complex_mol_subunit_spec", "complex_mol_subunit_component",
  "complex_mol_geometry", "complex_mol_subunit_locations", "subunit_coord",
  "complex_mol_subunit_location", "complex_mol_relationships",
  "complex_mol_relationship_list", "complex_mol_relationship",
  "complex_mol_rates", "complex_mol_rate_list", "complex_mol_rate",
  "complex_mol_rate_rules", "complex_mol_rate_rule",
  "complex_mol_rate_clauses", "complex_mol_rate_clause_list",
  "complex_mol_rate_clause", "equal_or_not", "existing_molecule",
  "existing_surface_molecule", "existing_molecule_opt_orient",
  "existing_macromolecule", "surface_classes_def",
  "define_one_surface_class", "define_multiple_surface_classes",
  "list_surface_class_stmts", "surface_class_stmt", "$@4",
  "existing_surface_class", "list_surface_prop_stmts", "surface_prop_stmt",
  "surface_rxn_stmt", "surface_rxn_type", "equals_or_to",
  "surface_class_mol_stmt", "surface_mol_stmt", "list_surface_mol_density",
  "list_surface_mol_num", "surface_mol_quant",
  "surface_class_viz_value_stmt", "rx_net_def", "list_rx_stmts", "rx_stmt",
  "rx_group_def", "reaction_group_name", "list_rxns", "list_dashes",
  "right_arrow", "left_arrow", "double_arrow", "right_cat_arrow",
  "double_cat_arrow", "reaction_arrow", "new_rxn_pathname", "rxn",
  "reactant_list", "reactant", "existing_molecule_or_subunit",
  "opt_reactant_surface_class", "reactant_surface_class", "product_list",
  "product", "rx_rate_syntax", "rx_rate1", "rx_rate2", "rx_dir_rate",
  "atomic_rate", "release_pattern_def", "new_release_pattern",
  "existing_release_pattern_xor_rxpn", "list_req_release_pattern_cmds",
  "train_count", "instance_def", "$@5", "physical_object_def",
  "object_def", "new_object", "start_object", "end_object",
  "list_opt_object_cmds", "opt_object_cmd", "transformation",
  "meta_object_def", "list_objects", "object_ref", "existing_object_ref",
  "$@6", "release_site_def", "release_site_def_new", "$@7",
  "release_site_geom", "release_region_expr", "release_site_def_old",
  "$@8", "release_site_geom_old", "list_release_site_cmds",
  "existing_num_or_array", "release_site_cmd", "site_size_cmd",
  "release_number_cmd", "constant_release_number_cmd",
  "gaussian_release_number_cmd", "volume_dependent_number_cmd",
  "concentration_dependent_release_cmd", "molecule_release_pos_list",
  "molecule_release_pos", "new_object_name", "polygon_list_def", "@9",
  "vertex_list_cmd", "single_vertex", "list_points",
  "element_connection_cmd", "list_element_connections",
  "element_connection", "list_opt_polygon_object_cmds",
  "opt_polygon_object_cmd", "remove_side", "$@10",
  "remove_element_specifier_list", "side_name", "element_specifier_list",
  "element_specifier", "incl_element_list_stmt", "excl_element_list_stmt",
  "just_an_element_list", "list_element_specs", "element_spec",
  "prev_region_stmt", "prev_region_type", "patch_statement", "patch_type",
  "in_obj_define_surface_regions", "list_in_obj_surface_region_defs",
  "in_obj_surface_region_def", "$@11", "$@12", "voxel_list_def", "$@13",
  "tet_element_connection_cmd", "element_connection_tet",
  "list_tet_arrays", "box_def", "$@14", "$@15", "opt_aspect_ratio_def",
  "existing_obj_define_surface_regions",
  "list_existing_obj_surface_region_defs",
  "existing_obj_surface_region_def", "$@16", "$@17", "$@18", "new_region",
  "list_opt_surface_region_stmts", "opt_surface_region_stmt",
  "set_surface_class_stmt", "surface_region_viz_value_stmt",
  "mod_surface_regions", "list_existing_surface_region_refs",
  "existing_surface_region_ref", "$@19", "output_def", "$@20",
  "output_buffer_size_def", "output_timer_def", "step_time_def",
  "iteration_time_def", "real_time_def", "list_count_cmds", "count_cmd",
  "count_stmt", "$@21", "custom_header_value", "custom_header",
  "exact_time_toggle", "list_count_exprs", "single_count_expr",
  "count_expr", "count_value", "$@22", "$@23", "file_arrow",
  "outfile_syntax", "existing_rxpn_or_molecule",
  "existing_molecule_required_orient_braces", "count_syntax",
  "count_syntax_1", "count_syntax_2", "count_syntax_3",
  "count_syntax_macromol", "subunit_molecule",
  "count_syntax_macromol_subunit", "$@24", "opt_macromol_relation_states",
  "macromol_relation_state_list", "macromol_relation_state",
  "macromol_relation_name", "count_location_specifier", "opt_hit_spec",
  "hit_spec", "opt_custom_header", "viz_output_def", "$@25",
  "list_viz_output_cmds", "viz_output_maybe_mode_cmd", "viz_mode_def",
  "viz_mesh_format_maybe_cmd", "viz_mesh_format_def",
  "viz_molecule_format_maybe_cmd", "viz_molecule_format_def",
  "viz_output_cmd", "viz_frames_def", "viz_filename_prefix_def",
  "viz_molecules_block_def", "list_viz_molecules_block_cmds",
  "viz_molecules_block_cmd", "viz_molecules_name_list_cmd",
  "optional_state", "viz_include_mols_cmd_list", "viz_include_mols_cmd",
  "existing_one_or_multiple_molecules", "viz_time_spec",
  "viz_molecules_time_points_def", "viz_molecules_time_points_cmds",
  "viz_molecules_time_points_one_cmd", "viz_iteration_spec",
  "viz_molecules_iteration_numbers_def",
  "viz_molecules_iteration_numbers_cmds",
  "viz_molecules_iteration_numbers_one_cmd", "viz_molecules_one_item",
  "viz_meshes_block_def", "list_viz_meshes_block_cmds",
  "viz_meshes_block_cmd", "viz_meshes_name_list_cmd",
  "viz_include_meshes_cmd_list", "viz_include_meshes_cmd",
  "viz_meshes_time_points_def", "viz_meshes_time_points_cmds",
  "viz_meshes_time_points_one_cmd", "viz_meshes_iteration_numbers_def",
  "viz_meshes_iteration_numbers_cmds",
  "viz_meshes_iteration_numbers_one_cmd", "viz_meshes_one_item",
  "volume_output_def", "volume_output_filename_prefix",
  "volume_output_molecule_list", "volume_output_molecule_decl",
  "volume_output_molecule", "volume_output_molecules",
  "volume_output_location", "volume_output_voxel_size",
  "volume_output_voxel_count", "volume_output_times_def", YY_NULLPTR
};
#endif

# ifdef YYPRINT
/* YYTOKNUM[NUM] -- (External) token number corresponding to the
   (internal) symbol number NUM (which must be that of a token).  */
static const yytype_uint16 yytoknum[] =
{
       0,   256,   257,   258,   259,   260,   261,   262,   263,   264,
     265,   266,   267,   268,   269,   270,   271,   272,   273,   274,
     275,   276,   277,   278,   279,   280,   281,   282,   283,   284,
     285,   286,   287,   288,   289,   290,   291,   292,   293,   294,
     295,   296,   297,   298,   299,   300,   301,   302,   303,   304,
     305,   306,   307,   308,   309,   310,   311,   312,   313,   314,
     315,   316,   317,   318,   319,   320,   321,   322,   323,   324,
     325,   326,   327,   328,   329,   330,   331,   332,   333,   334,
     335,   336,   337,   338,   339,   340,   341,   342,   343,   344,
     345,   346,   347,   348,   349,   350,   351,   352,   353,   354,
     355,   356,   357,   358,   359,   360,   361,   362,   363,   364,
     365,   366,   367,   368,   369,   370,   371,   372,   373,   374,
     375,   376,   377,   378,   379,   380,   381,   382,   383,   384,
     385,   386,   387,   388,   389,   390,   391,   392,   393,   394,
     395,   396,   397,   398,   399,   400,   401,   402,   403,   404,
     405,   406,   407,   408,   409,   410,   411,   412,   413,   414,
     415,   416,   417,   418,   419,   420,   421,   422,   423,   424,
     425,   426,   427,   428,   429,   430,   431,   432,   433,   434,
     435,   436,   437,   438,   439,   440,   441,   442,   443,   444,
     445,   446,   447,   448,   449,   450,   451,   452,   453,   454,
     455,   456,   457,   458,   459,   460,   461,   462,   463,   464,
     465,   466,   467,   468,   469,   470,   471,   472,   473,   474,
     475,   476,   477,   478,   479,   480,   481,   482,   483,   484,
     485,   486,   487,   488,   489,   490,   491,   492,   493,   494,
     495,   496,   497,   498,   499,   500,   501,   502,   503,   504,
     505,   506,   507,   508,   509,   510,   511,   512,   513,   514,
     515,    61,    38,    58,    43,    45,    42,    47,    94,   516,
      91,    93,    59,    39,    44,   123,   125,    40,    41,    62,
      60,    64
};
# endif

#define YYPACT_NINF -1157

#define yypact_value_is_default(Yystate) \
  (!!((Yystate) == (-1157)))

#define YYTABLE_NINF -445

#define yytable_value_is_error(Yytable_value) \
  0

  /* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
     STATE-NUM.  */
static const yytype_int16 yypact[] =
{
    2603,  -165,  -153,   -75,   -33,    45,    55,   114,    34,    34,
      13,    78,    34,    34,   119,   128,   161,   184,   199,   258,
     238, -1157,   284,   292,   304,   306,   318,   343,   351,   383,
     275,   282, -1157, -1157, -1157,   341,   373,   391,   403,   392,
     413,   398,   416,   418,   419, -1157,   408,   411,   412,   692,
    2603, -1157,   -46, -1157, -1157,   433, -1157, -1157,   606, -1157,
   -1157, -1157, -1157, -1157, -1157, -1157, -1157, -1157, -1157, -1157,
   -1157,   436, -1157, -1157, -1157, -1157, -1157, -1157, -1157, -1157,
   -1157, -1157, -1157, -1157,   124, -1157, -1157, -1157, -1157,   527,
   -1157, -1157, -1157, -1157, -1157, -1157, -1157, -1157,   350,   350,
     -37,  2403,   -37,  2403,  2403, -1157, -1157, -1157, -1157,   438,
      34,   232, -1157,   445, -1157,   446, -1157,    34,    34,  2403,
      34,    34,    34,   -37,    34,  2403,  2403,   350,  2403,  2403,
    2403,  2403,   588,    34,   400,   -37,   -37,  2078,  2403,   545,
    2403,    34,  2403,  2403,  2403, -1157,   625,  1366, -1157, -1157,
    1769,   437,  -158,   447, -1157, -1157,   447, -1157,   447, -1157,
   -1157,   447,   447,   447, -1157, -1157, -1157, -1157, -1157, -1157,
   -1157, -1157,   448, -1157, -1157, -1157, -1157, -1157,   461, -1157,
   -1157,   451,   453,   456,   458,   459,   462,   464,   465, -1157,
     470,   474,   478,   479,   481, -1157, -1157, -1157, -1157,   482,
   -1157,   483,   486,   487,   488,  2403,  2403,  2403, -1157,   186,
   -1157, -1157, -1157, -1157, -1157,  1294,    -6,   672,   477,    23,
    -128, -1157,    34,    34, -1157,   334, -1157,   193, -1157, -1157,
   -1157,  -183, -1157, -1157, -1157, -1157,    62, -1157, -1157, -1157,
      84, -1157,   672, -1157,   492,   498,   506,   461, -1157,   607,
   -1157,   672,   672, -1157,   672,   672,   672,   672, -1157, -1157,
   -1157,   514,   511,    92, -1157,   526,   530,   531,   533,   539,
     541,   542,   546,   547,   552,   553,   554,   558,   559,   562,
     563,   564,   571,   709, -1157,   461, -1157,   528, -1157,   672,
     672,   590, -1157,   672, -1157,   556,   672,   672,   672,   701,
     605,   725,   610,   611,   626,   629,   631,   633,   635,   636,
     638,   646,   648,   653,   655,   656,   658,   660,   663,   665,
    1065, -1157,  2177,   578, -1157, -1157,   672,   819, -1157,  1008,
     461,   668,   -37, -1157, -1157, -1157, -1157, -1157,   853,    34,
   -1157,   688, -1157,   688,   -37,   -37,  2403,  2403,  2403,  2403,
    2403,  2403,  2403,  2403,  2403,  2403,  2403,  2403,  2403,  2403,
    2403,  2403,   -37,  2403, -1157, -1157,   700, -1157, -1157,  2403,
    2403,  2403,  2403,  2403, -1157,  2403, -1157,   777,   673,   681,
     475, -1157, -1157, -1157,   669,   676, -1157, -1157, -1157,  2403,
   -1157,   326, -1157, -1157, -1157, -1157, -1157,  -147,    34,  -101,
     -31, -1157, -1157, -1157,   689, -1157, -1157, -1157,   -37,   -37,
      34, -1157, -1157, -1157,   350,   350,   350,    17,   350,   350,
    1624,   350,   350,   350,  2403,   350,    17,   350,   350,   350,
      17,    17, -1157, -1157,  -119, -1157,  2403,   -39,   -37,   697,
     717, -1157,   -37,   708,   374, -1157,   -41,   -41,   -41,  2403,
     -41,  2403,   -41,  2403,   -41,   -41,  2403,   -41,   -41,   -41,
     -41,   -41,   -41,   -41, -1157, -1157,  2403,   233, -1157,   672,
     696,   710,   801, -1157,   452,    34, -1157, -1157,   778,   707,
     773,   538,   930, -1157, -1157,   729,   774,   812,   878,   893,
     902,   924,   929,   954,   970,   781,  1104,  1252,  1268,   993,
    1043,   -18,  1048, -1157,   363,   363,   731,   731, -1157,  1312,
     740, -1157,  2403,  2403,   742,   743,   786,  -147, -1157,   886,
   -1157, -1157, -1157, -1157,   334, -1157, -1157,   752,  -136, -1157,
      63, -1157, -1157, -1157,   -70,   762,   763,   764,   765,   767,
   -1157,    20,    34, -1157,   751,   759, -1157, -1157, -1157, -1157,
   -1157, -1157, -1157, -1157, -1157, -1157,   672, -1157, -1157, -1157,
   -1157,   672, -1157, -1157, -1157, -1157, -1157, -1157, -1157,  2039,
   -1157,   672,   770,   772,   775,   -34, -1157, -1157, -1157, -1157,
     395,   776,   797, -1157,   461,    34,   793, -1157,   802, -1157,
   -1157, -1157, -1157, -1157, -1157,   672, -1157,   672, -1157,   672,
   -1157, -1157,   672, -1157, -1157, -1157, -1157, -1157, -1157, -1157,
     637, -1157,  2177,   -37,  -158,   183,   -91, -1157,   803,   538,
    -158,   790, -1157,   806,   809,   814,   811,   820,   816,   831,
     838,   839, -1157, -1157,   827,   538, -1157,   844, -1157, -1157,
   -1157, -1157, -1157,   833, -1157,   203, -1157, -1157, -1157, -1157,
   -1157, -1157, -1157, -1157, -1157, -1157,  2403,  2403,  2403,  2403,
   -1157, -1157, -1157, -1157,  2403,  -158,   891,   672,   672,  2403,
    2403, -1157,   987,  -137, -1157, -1157, -1157,   848, -1157, -1157,
     752, -1157,   752, -1157, -1157,   276, -1157,  2403,  2265,  2403,
    2403,  2403, -1157,    34,   841,   843, -1157, -1157,   858, -1157,
   -1157, -1157,  -134, -1157, -1157, -1157, -1157,   850,   259, -1157,
   -1157,   222, -1157, -1157,   668, -1157,  -158,  2403,  -158,   862,
     863, -1157,    15, -1157, -1157, -1157, -1157,   290, -1157, -1157,
   -1157, -1157, -1157,   439,   864,   289, -1157, -1157, -1157,   865,
    -158,   866,   874,  2403, -1157,   461,   854,   857, -1157,   447,
     872,   873,   894, -1157, -1157, -1157, -1157,   328,   538, -1157,
   -1157,    73,  -158, -1157,  2403,  2403,  1013,  -158,    34,    34,
    2403,    34,  2403,  1022,   -91, -1157,  2364,  -158, -1157, -1157,
    1145,  1153,  1171,  1189,  1373, -1157,   903,   129, -1157,   672,
     672,   911,   887, -1157, -1157,   253, -1157, -1157,   -70,   279,
     913, -1157, -1157,   672, -1157,   672, -1157,   672,   672,   672,
     916,    34,    34,  2403, -1157, -1157,    33, -1157, -1157,   917,
     918, -1157, -1157, -1157, -1157, -1157, -1157,   672, -1157,   703,
     350,   154, -1157, -1157, -1157, -1157, -1157,   490,   920,   899,
     909,   -30, -1157, -1157, -1157, -1157, -1157,    34, -1157,  2364,
     925,   103,   357, -1157,  -158, -1157,  -158,  2364,  -158, -1157,
   -1157, -1157, -1157, -1157, -1157,    79,   514, -1157,   440,   -91,
   -1157, -1157, -1157, -1157,    87,   -91,   672,   672,   926, -1157,
   -1157,  -158,    94, -1157,   672, -1157, -1157,   672,   937, -1157,
    1343, -1157, -1157, -1157, -1157,   151, -1157,    -1, -1157, -1157,
   -1157, -1157,  2403,  2403,   928, -1157,   985,  2403, -1157, -1157,
      34,  1824,  1824, -1157, -1157,   668,   281, -1157,    34, -1157,
    2403,   334,   944,   206, -1157,   230, -1157,   672,   334, -1157,
     933,    34,  2403, -1157, -1157,   461, -1157, -1157, -1157,   936,
     932, -1157,   154,   154, -1157,    72, -1157,  1350, -1157, -1157,
   -1157,   -37,   239,   277, -1157, -1157, -1157,  1343, -1157, -1157,
   -1157,  2364,   949,   951,   955,   939,  2403,  1197, -1157,   965,
   -1157, -1157,   296,    79,    79,    79, -1157, -1157, -1157, -1157,
    2403, -1157, -1157, -1157,  2403, -1157, -1157,   950,   953,    51,
   -1157, -1157, -1157,   672,  1392,   312, -1157,  1023,   967, -1157,
     672, -1157,    34, -1157, -1157, -1157,   279, -1157,   672, -1157,
    2403, -1157, -1157, -1157, -1157, -1157,   851, -1157,   672,   974,
    2403,   154,   975, -1157,   390,   154,  -175,   -37,   154,   154,
     154,   154, -1157,   461,   972,   976,   977,   -57, -1157, -1157,
   -1157, -1157,   979,   981,   988,   -44, -1157, -1157, -1157, -1157,
   -1157,  -158,  2403,  -158, -1157,  1257,  1001, -1157,   -91,  2403,
   -1157,   983,   983, -1157,   435,   569,    34, -1157, -1157,  2403,
    1003,  2403,   997,  -157, -1157, -1157,  1058, -1157,   998,   672,
    1007, -1157, -1157,  1032, -1157, -1157,   851, -1157, -1157, -1157,
   -1157,  1034, -1157,  1035,   272,  1281,   613,   272, -1157, -1157,
    1018,  1020,  1039,   -37,   461,   174,   174, -1157, -1157,    49,
   -1157,    49, -1157, -1157,    85, -1157,    85, -1157, -1157, -1157,
     672, -1157, -1157,  2403, -1157, -1157,   672,  1061, -1157,  1062,
     246, -1157,  1044,  1370,   672,    34, -1157,  2403, -1157, -1157,
     251,  1049,  1057, -1157,  1055,  1064, -1157, -1157,    34,  -158,
    1060,  1066,  1068,  1070,  1071,  1067, -1157, -1157, -1157, -1157,
   -1157, -1157, -1157,  1076, -1157, -1157,  1069, -1157, -1157, -1157,
   -1157, -1157,    -3, -1157,  1072,    25,    41, -1157,  1073, -1157,
   -1157, -1157,    21, -1157,  1075,    16,    27, -1157,  1077,   672,
      -1,  2403,  2403, -1157, -1157, -1157, -1157, -1157, -1157, -1157,
   -1157, -1157, -1157,   682,  1074, -1157,   851, -1157,  1083, -1157,
   -1157,   672,   316, -1157,  1088, -1157, -1157, -1157,   191,   191,
     263, -1157,  1086,    34,   -85, -1157,   -85,   -85, -1157, -1157,
   -1157, -1157, -1157,    32,  1102, -1157, -1157,   514,  1102,  1102,
   -1157, -1157, -1157,    10, -1157, -1157,    32,  1102, -1157, -1157,
   -1157, -1157,  1102, -1157, -1157,    10, -1157,  1063,   116,   851,
    2403, -1157,   191,  1103,  2403,  -158,   255,   321,   345, -1157,
    -158, -1157, -1157,   514, -1157,  1091,  1091,  1091, -1157, -1157,
   -1157,  2403, -1157, -1157, -1157, -1157, -1157, -1157, -1157, -1157,
   -1157, -1157,   -91, -1157,  1105,   851,   672, -1157,  -158,   672,
   -1157, -1157,  1099, -1157, -1157, -1157, -1157,   334,   794, -1157,
   -1157, -1157,   672, -1157,  2403, -1157, -1157,    -4,  1112, -1157,
   -1157, -1157, -1157, -1157, -1157, -1157, -1157, -1157,  1230,   302,
   -1157,    12,    11, -1157,  1114,  1119, -1157,    34, -1157, -1157,
   -1157, -1157,    34, -1157, -1157,  2403,    34, -1157,  1115, -1157,
     672, -1157,    34,  1110, -1157,   -12, -1157,    12,  1116,    34,
   -1157,    34,   -85, -1157, -1157, -1157
};

  /* YYDEFACT[STATE-NUM] -- Default reduction number in state STATE-NUM.
     Performed when YYTABLE does not specify something else to do.  Zero
     means the default is an error.  */
static const yytype_uint16 yydefact[] =
{
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,   369,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,   218,   219,   220,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,    26,     0,     0,     0,     0,
       2,     3,   377,     5,     6,     0,     7,   113,     0,   114,
     115,   116,   117,   118,   119,     8,     9,    10,    11,    13,
      12,     0,    14,   221,   222,   223,    15,   281,   282,    16,
      17,    19,    18,   371,     0,   372,   373,   394,   393,     0,
     375,   376,   374,    20,    21,    22,    23,    24,     0,     0,
       0,     0,     0,     0,     0,   244,   240,   229,   224,     0,
       0,     0,   359,     0,   230,     0,   283,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,   529,
       0,     0,     0,     0,     0,   606,     0,     0,     1,     4,
       0,     0,     0,     0,   413,   414,     0,   415,     0,   412,
     416,     0,     0,     0,    36,    38,    40,    39,    35,    37,
     202,   201,     0,   109,    25,   108,   112,   185,    27,   106,
     107,     0,     0,     0,     0,     0,     0,     0,     0,    71,
       0,     0,     0,     0,     0,    94,    96,    95,    72,     0,
      97,     0,     0,     0,     0,     0,     0,     0,    75,   190,
      67,    69,    70,    68,   186,   193,   190,   212,     0,     0,
       0,   226,     0,     0,   277,    41,   338,     0,   313,   315,
     316,   340,   335,   337,   361,   287,     0,   285,    28,   510,
       0,   508,   207,   124,     0,     0,     0,    56,   377,     0,
     370,   208,   200,   188,   213,   214,   215,   216,   210,   211,
     209,     0,     0,     0,   523,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,   137,   125,   126,     0,   205,   204,
     206,     0,   527,   198,    61,     0,   197,   199,   203,   610,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,   162,     0,    62,    59,    60,     0,    73,    57,    74,
       0,    58,     0,    66,   217,    63,    64,   378,     0,     0,
     395,     0,   410,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,   105,   104,     0,   192,   191,     0,
       0,     0,     0,     0,   187,     0,   189,     0,     0,     0,
     233,   225,   227,   318,     0,     0,    44,    49,    50,     0,
     279,    42,    45,    46,    43,   312,   314,     0,     0,     0,
       0,   290,   284,   286,     0,   507,   509,   123,     0,     0,
       0,   525,   522,   524,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,   136,   138,     0,   134,     0,     0,     0,     0,
     617,   611,     0,     0,     0,   685,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,   161,   163,     0,     0,    52,    54,
       0,     0,   377,   390,     0,   380,   387,   389,     0,     0,
       0,     0,     0,   126,   110,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,    76,    99,   100,   101,   102,   103,   194,
       0,   241,     0,     0,     0,     0,   236,     0,   339,     0,
      47,    48,   336,   289,    41,   341,   321,     0,     0,   328,
       0,   330,   329,   331,     0,     0,     0,     0,     0,     0,
     358,     0,     0,   126,     0,     0,   515,   158,   139,   146,
     154,   160,   159,   141,   148,   149,   156,   155,   157,   145,
     142,   144,   140,   151,   147,   150,   143,   153,   152,     0,
     131,   530,     0,     0,     0,     0,   531,   532,   533,   126,
       0,     0,   621,   618,   684,     0,     0,   686,     0,   184,
     182,   183,   164,   169,   178,   179,   170,   168,   167,   173,
     172,   174,   175,   176,   180,   165,   166,   181,   171,   177,
       0,    65,     0,     0,     0,     0,     0,   388,     0,     0,
       0,     0,   496,     0,     0,     0,     0,     0,     0,     0,
       0,     0,   429,   430,     0,   380,   417,     0,   422,   431,
     432,   433,   434,     0,   445,     0,    92,    89,    88,    90,
      84,    86,    77,    83,    78,    79,     0,     0,     0,     0,
      85,    91,    98,    87,     0,     0,     0,   232,   231,     0,
       0,   237,   238,     0,   319,    51,   342,   324,   322,   323,
       0,   325,     0,   345,   346,     0,   343,     0,     0,     0,
       0,     0,   300,     0,     0,     0,   298,   299,     0,   288,
     291,   292,     0,   293,   303,   294,   514,     0,     0,   135,
      31,     0,   130,   128,   129,   127,     0,     0,     0,     0,
       0,   542,     0,   537,   539,   540,   541,     0,   615,   616,
     613,   614,   612,     0,     0,     0,   622,   688,   689,   687,
       0,     0,     0,     0,    53,   122,     0,     0,    32,     0,
       0,     0,     0,   379,   386,   381,   382,     0,   380,   448,
     449,     0,     0,   380,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,   418,     0,     0,   455,   111,
       0,     0,     0,     0,   195,   245,     0,     0,   246,   235,
     234,     0,     0,   317,   320,     0,   326,   327,     0,     0,
     332,   347,   348,   362,   368,   367,   366,   363,   365,   364,
       0,     0,     0,     0,   302,   301,     0,   511,   132,     0,
       0,   526,   518,   516,   517,   519,   535,   534,   536,     0,
       0,     0,   528,   538,   133,   620,   619,     0,     0,     0,
       0,     0,   608,   626,   625,   627,   628,     0,   691,     0,
       0,   694,     0,   120,     0,   391,     0,     0,     0,   400,
     401,   404,   402,   399,   403,     0,   398,   405,   397,     0,
     447,   450,   499,   500,     0,     0,   439,   440,     0,   420,
     421,     0,     0,   441,   435,   360,   427,   426,     0,   411,
     419,   424,   423,   425,   454,     0,   452,   380,    80,    81,
      93,    82,     0,     0,     0,   247,     0,     0,   228,   344,
       0,     0,     0,   356,   354,   355,     0,   351,     0,   334,
       0,    41,     0,     0,   306,     0,   308,   311,    41,   295,
       0,     0,     0,   544,   545,   546,   547,   548,   561,     0,
       0,   564,     0,     0,   552,     0,   549,   604,   553,   624,
     623,     0,     0,     0,   607,   609,   690,    66,    33,   692,
      34,     0,     0,     0,     0,     0,     0,   505,   380,     0,
     384,   383,     0,     0,     0,     0,   396,   498,   501,   497,
       0,   443,   428,   442,     0,   451,   453,     0,     0,     0,
     456,   457,   458,   196,   251,     0,   249,     0,     0,   242,
     239,   280,     0,   352,   353,   349,     0,   333,   297,   278,
       0,   304,   307,   305,   309,   296,     0,   520,   521,     0,
       0,     0,     0,   559,     0,     0,     0,     0,     0,     0,
       0,     0,   551,   629,     0,     0,     0,     0,   661,   663,
     664,   665,     0,     0,     0,     0,   631,   633,   634,   635,
     693,     0,     0,     0,   683,     0,     0,   502,     0,     0,
     406,   407,   408,   409,     0,     0,     0,   459,   446,     0,
       0,     0,     0,     0,   254,   260,     0,   357,     0,   310,
       0,   489,   486,     0,   488,   485,   512,   470,   472,   473,
     474,     0,   475,     0,     0,     0,     0,     0,   554,   550,
       0,     0,   566,     0,   605,   555,   556,   557,   558,     0,
     668,     0,   660,   662,     0,   640,     0,   630,   632,   696,
     695,   697,    55,     0,   455,   392,   385,     0,   436,     0,
       0,   491,     0,     0,   252,     0,   250,     0,   253,   255,
       0,     0,     0,   350,     0,     0,   471,   515,     0,     0,
       0,     0,   572,     0,     0,     0,   574,   575,   576,   577,
     581,   563,   560,     0,   567,   570,   568,   571,   543,   680,
     681,   682,     0,   677,     0,     0,     0,   673,     0,   657,
     659,   658,     0,   654,     0,     0,     0,   648,     0,   506,
     503,     0,     0,   490,   492,   493,   469,   466,   464,   465,
     467,   468,   463,   481,     0,   483,   461,   462,   478,   479,
     248,   256,     0,   259,     0,   261,   264,   243,     0,     0,
       0,   484,     0,     0,     0,   573,     0,     0,   562,   565,
     569,   676,   678,     0,   638,   666,    30,    29,   638,   638,
     667,   672,   674,     0,   653,   655,     0,   638,   636,   644,
     643,   639,   638,   647,   649,     0,   380,     0,     0,     0,
       0,   460,     0,     0,     0,     0,     0,     0,     0,   513,
       0,   583,   591,   593,   592,   594,   594,   594,   651,   652,
     679,     0,   671,   670,   669,   645,   646,   675,   656,   642,
     641,   650,     0,   437,     0,   494,   482,   480,     0,   257,
     262,   263,     0,   265,   477,   476,   487,    41,     0,   580,
     578,   579,   637,   504,     0,   515,   258,     0,     0,   601,
     603,   598,   600,   597,   602,   599,   596,   595,     0,     0,
     271,     0,     0,   267,     0,   270,   272,     0,   438,   495,
     275,   276,     0,   266,   268,     0,     0,   582,   585,   274,
     269,   273,     0,     0,   590,     0,   587,     0,     0,     0,
     586,     0,     0,   588,   589,   584
};

  /* YYPGOTO[NTERM-NUM].  */
static const yytype_int16 yypgoto[] =
{
   -1157, -1157, -1157,  1286, -1017,     0,   -98,  -112, -1157,  -131,
    -586,  -779,   -76,  -510, -1157,  1002,  1004,   236, -1157,   782,
   -1157, -1157,  1258,  -139,  -145,  -149, -1157,   931,  -340,  -130,
    -140, -1157,  -123,   -84,  -133, -1157, -1157, -1157, -1157, -1157,
   -1157,   516,  -115,  -409, -1157, -1157, -1157, -1157, -1157, -1157,
   -1157, -1157,  1118,   593,   167, -1157, -1157,  1078,   529, -1157,
    1186, -1157, -1157, -1157, -1157, -1157, -1157, -1157, -1157, -1157,
     -28, -1157, -1157, -1157, -1157, -1157, -1157, -1157, -1157, -1157,
   -1157, -1157, -1157, -1157,   616, -1157,   333, -1157, -1157, -1157,
     342, -1157, -1157, -1157, -1157, -1157, -1157, -1157,    74, -1157,
   -1157,    68,    65, -1157, -1157,  -220,   201, -1157, -1157, -1157,
   -1157,   -25, -1157,   494, -1157, -1157, -1157, -1157, -1157, -1157,
     885, -1157, -1157,  -659, -1157, -1157, -1157,  1201, -1157, -1157,
   -1157,  -323,   -47, -1157, -1157, -1157, -1157, -1157, -1157,  -465,
   -1157,  1033,  -509, -1157, -1157, -1157,   634, -1157, -1157, -1157,
     423,  -270, -1157, -1157, -1157, -1157, -1157, -1157, -1157, -1157,
    -277,  -109,  -121,  -738,  -616, -1157, -1157,  1309, -1157,   959,
   -1157, -1157, -1157, -1157, -1157, -1157,  -749, -1157, -1157, -1157,
     821, -1157,  -568, -1157, -1157, -1157, -1157, -1157, -1157, -1157,
     560, -1157, -1157, -1157,  1098,   684, -1157, -1157, -1157,   557,
     327, -1157, -1157, -1157, -1157, -1157, -1068, -1041, -1157, -1157,
   -1157,  -559,   197, -1157, -1157, -1157, -1157, -1157, -1157,   330,
   -1157, -1157, -1157, -1157, -1157,   587, -1157, -1157, -1157, -1157,
   -1157, -1157, -1157,  1222, -1157, -1157, -1157,   922, -1072, -1157,
   -1157, -1157, -1157, -1157,  1205, -1157, -1157, -1157, -1157, -1157,
   -1157, -1157, -1157, -1157,   747, -1157, -1157, -1157, -1157, -1157,
   -1157,   454,  -169, -1157, -1157, -1157, -1157, -1157, -1157, -1157,
     376, -1157, -1157, -1157, -1157,   109, -1157, -1157, -1157, -1157,
     112, -1157, -1156,  -605, -1157, -1157, -1157, -1157, -1157, -1157,
   -1157, -1157, -1157, -1157, -1157,   639, -1157, -1157, -1157, -1157,
     431, -1157,  -658, -1157, -1157, -1157,   226, -1157, -1157,   291,
     241, -1157, -1157,   300, -1037, -1157, -1157,   467, -1157, -1157,
   -1157, -1157, -1157,   309, -1157, -1157,   317, -1020, -1157, -1157,
   -1157,  1046,   641, -1157, -1157, -1157, -1157, -1157
};

  /* YYDEFGOTO[NTERM-NUM].  */
static const yytype_int16 yydefgoto[] =
{
      -1,    49,    50,    51,   175,   208,   177,   261,  1238,   867,
     958,   959,   547,   390,   391,   392,   393,   394,   467,   468,
      53,    54,    55,   913,   748,   335,   336,   326,   210,   211,
     914,   212,   213,   285,   179,   180,    56,    57,    58,   746,
      59,   244,   286,   434,   715,    60,    61,    62,    63,    64,
      65,   283,   284,   548,   553,    66,   320,   321,   592,    67,
     374,   216,    68,    69,    70,    71,    72,    73,    74,   220,
     108,   109,   115,   380,   516,   672,   792,    75,   218,   666,
    1076,   106,   511,   787,   788,   995,   996,   906,  1073,  1212,
    1074,   999,  1140,  1215,  1142,  1266,  1303,  1332,  1333,  1334,
    1335,  1336,  1342,   225,   922,   226,  1002,    76,    77,    78,
     236,   116,   401,   524,   541,   700,   701,   702,   816,   703,
     822,   923,   925,   924,   705,    79,   227,   228,   229,   384,
     673,   795,   529,   530,   531,   532,   533,   534,   919,   230,
     231,   232,   233,   399,   525,   685,   686,   800,   801,   802,
     916,   917,    80,   113,   886,   400,   806,    81,   124,    82,
      83,    84,   339,   754,   616,   755,   756,    85,   475,   476,
     477,   968,    86,    87,   478,   619,   868,    88,   481,   162,
     635,   893,   636,   637,   638,   639,   640,   641,   642,   882,
     883,    89,    90,   778,   480,   760,   761,   644,   895,   896,
     897,   990,   991,  1133,  1204,  1205,  1086,  1087,  1088,  1089,
    1207,  1208,  1209,  1090,  1091,  1092,  1093,   992,  1130,  1131,
    1259,  1315,    91,   763,   622,   873,   874,    92,  1124,  1256,
    1057,    93,   240,   241,   404,   930,  1147,  1132,   711,   823,
     824,   825,    94,   263,   264,   546,    95,   437,   292,   575,
     576,   577,   578,   722,   723,   724,   831,   936,   725,   726,
     945,   946,   947,   948,  1019,  1022,  1103,  1168,  1153,  1154,
    1155,  1156,  1157,  1158,  1159,  1348,  1160,  1307,  1353,  1355,
    1356,  1357,  1275,  1309,  1327,  1032,    96,   299,   841,   440,
     441,   582,   583,   735,   736,   842,   843,   844,   845,  1045,
    1046,  1047,  1282,  1185,  1251,  1252,  1287,  1048,  1186,  1187,
    1280,  1049,  1182,  1183,  1184,   846,  1037,  1038,  1039,  1175,
    1240,  1040,  1176,  1177,  1041,  1172,  1173,  1174,    97,   301,
     444,   445,   738,   739,   588,   742,   851,   965
};

  /* YYTABLE[YYPACT[STATE-NUM]] -- What to do in state STATE-NUM.  If
     positive, shift that token.  If negative, reduce the rule whose
     number is the opposite.  If YYTABLE_NINF, syntax error.  */
static const yytype_int16 yytable[] =
{
      52,   325,   262,   385,   214,  1169,   239,   334,   105,   107,
     328,   324,   112,   114,   676,   249,   178,   331,   178,   774,
     327,   287,   170,   171,   692,   684,  1285,   329,   747,  1179,
    1247,   535,   338,   589,   759,  1179,   889,   340,  1234,   247,
     341,   342,   343,  -121,  1278,  1146,   551,   928,  1330,  1169,
      52,   253,   674,   172,   838,  1034,   260,  1169,   693,   987,
     590,  1072,   473,  1330,   719,  1206,   330,   775,  1042,   173,
    1276,  1277,   572,   367,   645,  1220,   528,  1151,   970,  1188,
    1151,   397,   221,   683,    45,  1340,  1100,   378,   379,  1101,
    1035,  1178,   237,  1179,  1170,    45,    98,   814,   398,   164,
     176,   839,   176,  1043,  1102,    45,    45,   750,    99,   751,
     107,   224,   322,   719,    45,   840,   972,   114,   238,  1138,
     243,   243,   243,   176,   248,  -444,   536,   815,   239,   678,
     223,   976,   262,   238,   708,   176,   176,   979,  1170,   793,
     223,   294,   869,   679,   752,   368,  1170,   875,   381,  1188,
     323,   153,   333,   926,   848,   569,  1178,    45,  1236,   570,
    1294,   537,   694,   695,   526,  1146,  1341,   165,  1249,   720,
     727,  1036,    45,   154,  1272,   759,   573,   166,   167,   527,
     174,   879,  1050,  1180,  1044,   753,   100,  1171,   574,  1180,
     775,  1295,   382,  1181,   181,   988,   155,   182,   473,  1181,
    1196,   938,   538,   539,   677,    45,  1365,   223,   794,   696,
     183,   403,   184,  1197,   962,   -61,   591,  1198,   720,  1112,
     107,   185,   383,   224,  1061,  1062,  1063,   224,   101,   483,
     474,  1171,  1117,   174,   470,   939,   114,   186,    45,  1171,
     238,   721,   174,  1329,   345,   540,   954,  1180,   178,   750,
    1359,   751,    45,    45,  1146,   168,   697,  1181,    45,  1360,
     662,   484,   189,   238,  1012,   367,  1014,    45,   967,   698,
     969,   187,   971,  1231,    45,    45,    45,   169,   501,   188,
     322,   989,   181,   156,  1199,   182,   752,  1343,   110,   684,
     721,   832,  1248,   543,   544,   981,   699,  1244,   183,   189,
     184,  1235,   322,  1253,    45,    45,   102,  1200,   680,   185,
     682,   157,   190,   191,   192,    45,   103,  1241,   963,   158,
    1125,    45,   910,   579,   193,   186,    45,  1068,   194,    45,
     964,   904,   176,   159,    45,   160,    45,   368,   402,   472,
     198,   552,   681,   322,   176,   176,  1025,   786,  1026,   870,
     552,  1034,  1058,   111,   552,   552,   865,   322,   584,   187,
     405,   195,   176,   977,   694,   695,   474,   188,   412,   172,
     982,   196,   197,   838,   940,   104,   859,   198,   222,   161,
     369,   370,   371,   372,   373,   173,  1035,   189,  1201,  1042,
     199,   941,   200,    45,   117,   201,    45,   224,   523,   860,
     190,   191,   192,   118,   202,   694,   695,   203,   176,   176,
     545,  1009,   193,   728,   204,   265,   194,   222,  1015,   942,
     839,   322,   119,  1202,  1043,    45,   729,   985,   266,   713,
     712,   943,   164,    45,   840,    45,   714,   267,   176,   327,
    1030,  1031,   176,   819,   694,   695,   329,   861,    45,   195,
     369,   370,   371,   372,   373,   205,   206,   835,   337,   196,
     197,   120,   730,   731,   836,   198,   268,  1036,   207,   395,
     223,   820,    45,   810,    45,   472,   121,   569,   199,   153,
     200,   779,  1011,   201,   819,   330,   269,   270,    45,   174,
    1150,   944,   202,    45,   586,   203,   174,    45,   821,   123,
     165,   154,   204,   749,   611,  1044,  1013,   612,   949,   223,
     166,   167,   820,   271,    45,   950,   862,   224,   678,   443,
     785,    45,  1193,   819,   155,   514,   515,  1213,   224,   745,
     224,  1301,   679,   569,   224,   122,   863,   818,   864,  1269,
     798,   272,   706,   205,   206,   125,   799,   732,   880,   881,
     133,   820,  1005,   126,  1313,  1006,   207,   134,   911,   912,
     973,   974,   975,  1222,   569,   127,   273,   128,   834,   323,
      45,   826,   966,   828,  1060,   274,   275,   276,  1339,   129,
    1283,  1284,   623,  1070,   277,   737,  1071,  1263,   168,  1289,
    1264,   278,  1304,   563,  1290,  1262,   929,   567,   568,   387,
     388,   624,   944,   944,   130,   865,   386,   387,   388,   389,
     169,   615,   131,   176,   333,   238,  1305,   872,   135,  1262,
     333,   369,   370,   371,   372,   373,   279,   891,   855,   371,
     372,   373,   894,   796,   625,   797,   892,   245,   246,   157,
    1292,  1003,  1004,   280,   132,   866,   327,   158,   281,  1127,
     136,   282,   137,   329,  1028,  1029,  1030,  1031,   626,  1267,
    1268,   159,   881,   160,   138,   333,   915,   139,  1098,   327,
     164,  1310,  1311,   224,   140,   141,   329,   142,   627,   143,
     144,   944,   628,   145,  1306,   944,   146,   147,   944,   944,
     944,   944,   148,   224,   150,   151,   629,   152,   163,   369,
     370,   371,   372,   373,   973,   974,   975,   161,   291,   960,
     300,  1128,  1316,   219,   332,   330,   333,   960,   333,   327,
     234,   235,   337,   345,   265,   344,   329,   327,   346,   872,
     347,   630,   631,   348,   329,   349,   350,   266,   165,   351,
     333,   352,   353,   632,   633,   935,   267,   354,   166,   167,
     894,   355,   377,   934,   937,   356,   357,   238,   358,   359,
     360,   333,   333,   361,   362,   363,   156,   333,   224,   224,
     407,   885,   408,  1023,  1024,   268,   890,   333,   915,   915,
     409,   327,   327,  1129,   410,   164,   411,   414,   329,   329,
     634,   415,   416,   172,   417,   269,   270,  1318,   224,   323,
     418,  1319,   419,   420,  1320,  1321,   435,   421,   422,   173,
     258,   921,   921,   423,   424,   425,   224,  1322,  1323,   426,
     427,   960,   271,   428,   429,   430,   168,   330,   330,   176,
     438,   327,   431,   369,   370,   371,   372,   373,   329,   439,
    -112,   259,   -75,   -75,   -75,   -75,   -75,   737,   169,   957,
     272,   436,  1096,   165,   333,   933,   333,   957,   333,  1105,
    1106,  1107,  1108,   166,   167,   238,   442,  1033,   743,  1324,
     443,   446,   447,   915,   333,   273,   327,  1028,  1029,  1030,
    1031,   333,   224,   329,   274,   275,   276,   448,  1325,  1326,
     449,  1162,   450,   277,   451,   333,   452,   453,   471,   454,
     278,   369,   370,   371,   372,   373,  1119,   455,  1121,   456,
    1001,   323,   323,  1260,   457,  1210,   458,   459,  1007,   460,
     174,   461,   330,   921,   462,   921,   463,  1080,  1081,  1082,
    -106,   523,   479,   510,   512,   279,   369,   370,   371,   372,
     373,   168,   513,  1104,   517,    45,   369,   370,   371,   372,
     373,   176,   280,  1083,   518,  1084,  1085,   281,   580,   542,
     282,   957,   581,   169,   369,   370,   371,   372,   373,   585,
     613,   614,  -444,   238,   238,   238,   593,   594,   503,   596,
     618,   598,   620,   600,   601,   432,   603,   604,   605,   606,
     607,   608,   609,   369,   370,   371,   372,   373,   621,   373,
     643,   665,  1077,   669,   670,  1167,   323,   646,   549,   550,
     671,   554,   555,   557,   558,   559,   560,   526,   562,   178,
     564,   565,   566,   687,   688,   689,   690,   176,   691,   709,
     710,   716,   209,   717,   215,   217,   718,   733,   369,   370,
     371,   372,   373,   734,  1239,   369,   370,   371,   372,   373,
     242,   333,   647,   333,   740,   656,   251,   252,   741,   254,
     255,   256,   257,  1237,   757,   762,   706,   764,   289,   290,
     765,   293,   767,   296,   297,   298,   369,   370,   371,   372,
     373,   768,   302,   -69,   -69,   -69,   -69,   -69,  1279,   766,
     648,   769,   770,  1274,  1152,  1274,  1274,  1152,  1286,   771,
     772,  1279,   773,   176,   303,   776,   304,   305,   777,   786,
    1286,   791,  1273,   678,  1273,  1273,   811,  1347,   812,   813,
    1300,   817,  1349,   829,   830,   837,   306,   849,   850,   847,
     706,   854,   853,   856,   857,   224,   364,   365,   366,   878,
    1214,  1347,   369,   370,   371,   372,   373,   888,  1221,   333,
     369,   370,   371,   372,   373,   858,   649,   369,   370,   371,
     372,   373,   675,   908,   307,   308,   369,   370,   371,   372,
     373,   650,   907,   903,   952,   238,   918,   920,   931,   932,
     651,   951,   309,   310,   953,  1250,   961,   980,   369,   370,
     371,   372,   373,   369,   370,   371,   372,   373,   984,   311,
     312,   313,   652,   997,   998,  1010,  1020,   653,  1016,  1021,
    1051,   314,  1052,   315,   316,  1054,  1053,  1056,   369,   370,
     371,   372,   373,  1001,   238,  1066,   238,   238,  1067,   317,
     318,  1274,   654,   333,   369,   370,   371,   372,   373,  1059,
    1141,  1072,  1075,   333,  1094,  1097,   333,  1109,   655,   975,
    1273,  1110,  1111,   469,  1114,   333,  1115,   369,   370,   371,
     372,   373,  1123,  1116,  1135,   333,  1302,  1137,  1144,  1143,
     333,   660,   -68,   -68,   -68,   -68,   -68,   485,   486,   487,
     488,   489,   490,   491,   492,   493,   494,   495,   496,   497,
     498,   499,   500,  1145,   502,  1148,  1149,  1164,   333,  1165,
     504,   505,   506,   507,   508,   319,   509,   369,   370,   371,
     372,   373,   369,   370,   371,   372,   373,  1331,  1166,  1195,
     519,   661,  1191,  1192,  1216,  1218,   663,   369,   370,   371,
     372,   373,  1331,  1217,  1219,  1223,   149,   224,  1228,  1293,
    1224,   464,   224,   389,  1226,  1227,  1331,  1229,  1230,  1265,
    1261,   556,  1354,  1233,  1243,   561,  1246,  1262,  1255,  1354,
    1270,   224,   238,  1281,  1298,  1308,  1314,   571,   369,   370,
     371,   372,   373,   181,  1317,  1337,   182,  1345,   657,  1196,
     595,  1346,   597,   302,   599,  1352,  1358,   602,  1225,   183,
    1362,   184,  1197,   520,   744,   521,  1198,   610,   465,   295,
     185,   433,   376,   905,  1136,   303,  1344,   304,   305,   369,
     370,   371,   372,   373,  1351,  1139,   186,   369,   370,   371,
     372,   373,  1361,   898,  1271,  1017,   704,   306,   396,  1078,
     522,   899,   909,   250,   617,   369,   370,   371,   372,   373,
     758,   482,   983,   667,   668,   871,  1080,  1081,  1082,   900,
     187,  1190,   986,   369,   370,   371,   372,   373,   188,  1297,
    1194,   978,   406,  1199,   707,   307,   308,   901,   413,   833,
    1364,  1363,  1083,  1163,  1084,  1085,  1118,  1254,   189,  1099,
     955,  1291,  1245,   309,   310,  1242,  1200,  1288,   956,  1232,
     587,   190,   191,   192,   369,   370,   371,   372,   373,     0,
     311,   312,   313,   193,  1113,     0,  1338,   194,     0,     0,
       0,     0,   314,     0,   315,   316,   369,   370,   371,   372,
     373,   369,   370,   371,   372,   373,   658,     0,  1122,     0,
     317,   318,   369,   370,   371,   372,   373,     0,     0,     0,
     195,     0,   659,   469,     0,   369,   370,   371,   372,   373,
     196,   197,  1161,     0,     0,     0,   198,   375,   369,   370,
     371,   372,   373,     0,     0,     0,     0,  1201,     0,   199,
       0,   200,     0,     0,   201,   664,   369,   370,   371,   372,
     373,     0,     0,   202,     0,     0,   203,   780,   781,   782,
     783,     0,     0,   204,     0,   784,     0,     0,     0,     0,
     789,   790,  1202,     0,     0,     0,   319,   -75,   -75,   -75,
     -75,   -75,    45,  1027,  1028,  1029,  1030,  1031,   803,   805,
     807,   808,   809,     0,     0,     0,     0,   181,     0,     0,
     182,     0,     0,     0,   205,   206,   902,   369,   370,   371,
     372,   373,     0,   183,     0,   184,     0,   207,   827,     0,
       0,     0,     0,     0,   185,  1069,   369,   370,   371,   372,
     373,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     186,     0,     0,     0,   852,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,   876,   877,     0,     0,     0,
       0,   884,     0,   887,   187,     0,   164,     0,     0,     0,
       0,     0,   188,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,   189,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,   927,   190,   191,   192,     0,     0,
       0,     0,     0,     0,     0,     0,     0,   193,     0,     0,
       0,   194,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,   181,     0,   165,   182,     0,     0,     0,     0,
       0,     0,     0,     0,   166,   167,     0,     0,   183,     0,
     184,     0,     0,     0,   195,     0,     0,     0,     0,   185,
       0,     0,     0,     0,   196,   197,     0,     0,     0,     0,
     198,     0,     0,     0,     0,   186,     0,     0,     0,     0,
       0,     0,     0,   199,     0,   200,     0,   181,   201,     0,
     182,     0,     0,   993,   994,     0,     0,   202,  1000,     0,
     203,     0,     0,   183,     0,   184,     0,   204,     0,   187,
       0,  1008,     0,     0,   185,     0,     0,   188,     0,   172,
       0,     0,   168,  1018,     0,     0,    45,   910,     0,     0,
     186,     0,     0,     0,     0,   173,     0,   189,     0,     0,
       0,     0,     0,     0,   169,     0,     0,     0,   205,   206,
     190,   191,   192,     0,     0,     0,     0,  1055,     0,     0,
       0,   207,   193,     0,   187,     0,   194,     0,     0,     0,
       0,  1064,   188,     0,   172,  1065,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     173,     0,   189,     0,     0,     0,     0,     0,     0,   195,
       0,  1079,     0,     0,     0,   190,   191,   192,     0,   196,
     197,  1095,     0,     0,     0,   198,     0,   193,     0,     0,
       0,   194,     0,     0,     0,     0,     0,     0,   199,     0,
     200,     0,     0,   201,     0,     0,     0,     0,     0,     0,
       0,     0,   202,  1120,     0,   203,   174,     0,     0,     0,
    1126,     0,   204,     0,   195,     0,     0,     0,     0,     0,
    1134,     0,   994,     0,   196,   197,     0,     0,     0,     0,
     198,    45,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,   199,     0,   200,     0,     0,   201,     0,
       0,     0,     0,   205,   206,     0,     0,   202,     0,   322,
     203,   174,   181,     0,     0,   182,   207,   204,     0,     0,
       0,     0,     0,     0,  1189,     0,     0,     0,   183,     0,
     184,     0,     0,     0,  1203,     0,    45,     0,  1211,   185,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,   181,     0,     0,   182,   186,     0,     0,   205,   206,
       0,     0,     0,     0,     0,     0,     0,   183,     0,   184,
       0,   207,     0,     0,     0,     0,     0,     0,   185,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,   187,
       0,     0,  1257,  1258,   186,     0,     0,   188,     0,   172,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,   173,     0,   189,     0,  1203,
    1203,     0,     0,     0,     0,     0,     0,     0,   187,     0,
     190,   191,   192,     0,     0,     0,   188,     0,     0,     0,
       0,     0,   193,     0,     0,     0,   194,     0,     0,     0,
     181,     0,     0,   182,     0,     0,   189,   288,     0,     0,
       0,  1296,     0,  1203,     0,  1299,   183,     0,   184,   190,
     191,   192,     0,     0,     0,     0,     0,   185,     0,   195,
       0,   193,  1312,     0,     0,   194,     0,     0,     0,   196,
     197,     0,     0,   186,     0,   198,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,   199,     0,
     200,     0,     0,   201,     0,  1328,     0,     0,   195,     0,
       0,     0,   202,     0,     0,   203,   174,   187,   196,   197,
       0,     0,   204,     0,   198,   188,     0,     0,   181,     0,
       0,   182,     0,     0,     0,     0,  1350,   199,     0,   200,
       0,    45,   201,     0,   183,   189,   184,     0,     0,     0,
       0,   202,     0,     0,   203,   185,     0,     0,   190,   191,
     192,   204,     0,   205,   206,     0,     0,     0,     0,     0,
     193,   186,     0,     0,   194,     0,   207,     0,     0,     0,
      45,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,   205,   206,     0,   187,     0,   195,     0,     0,
       0,     0,     0,   188,     0,   207,     0,   196,   197,     0,
       0,     0,     0,   198,     0,     0,     0,   181,     0,     0,
     182,     0,     0,   189,     0,     0,   199,     0,   200,     0,
       0,   201,     0,   183,     0,   184,   190,   191,   192,     0,
     202,     0,     0,   203,   185,     0,     0,     0,   193,     0,
     204,     0,   194,     0,     0,     0,   181,     0,     0,   182,
     186,     0,     0,     0,     0,     0,     0,     0,     0,    45,
       0,     0,   183,     0,   184,     0,     0,     0,     0,     0,
       0,     0,     0,   185,     0,   195,     0,     0,     0,     0,
       0,   205,   206,     0,   187,   196,   197,   466,     0,   186,
       0,   198,   188,     0,   207,     0,     0,     0,     0,     0,
       0,     0,     0,     0,   199,     0,   200,     0,     0,   201,
       0,     0,   189,     0,     0,     0,     0,     0,   202,     0,
       0,   203,     0,   187,     0,   190,   191,   192,   204,     0,
       0,   188,     0,     0,     0,     0,     0,   193,     0,     0,
       0,   194,     0,     0,   804,     0,     0,    45,     0,     0,
       0,   189,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,   190,   191,   192,     0,     0,   205,
     206,     0,     0,     0,   195,     0,   193,     0,     0,     0,
     194,     0,   207,     0,   196,   197,     0,     0,     0,     0,
     198,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,   199,     0,   200,     0,     0,   201,     0,
       0,     0,     0,   195,     0,     0,     0,   202,     0,     0,
     203,     0,     0,   196,   197,     0,     0,   204,     0,   198,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,   199,     0,   200,     0,    45,   201,     1,     0,
       0,     0,     0,     0,     0,     0,   202,     0,     0,   203,
       0,     0,     0,     0,     0,     0,   204,     0,   205,   206,
       0,     0,     0,     0,   322,     2,     3,     4,     5,     6,
       0,   207,     0,     7,     0,    45,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     8,     9,    10,    11,
      12,    13,    14,    15,     0,     0,     0,   205,   206,     0,
       0,     0,    16,     0,     0,     0,     0,     0,     0,     0,
     207,     0,     0,     0,     0,     0,    17,     0,     0,     0,
       0,     0,     0,     0,    18,    19,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,    20,     0,     0,     0,
      21,     0,     0,    22,     0,     0,     0,    23,    24,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
      25,    26,    27,    28,     0,    29,     0,     0,     0,     0,
       0,     0,    30,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,    31,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
      32,    33,    34,     0,     0,     0,    35,    36,     0,     0,
       0,    37,    38,     0,     0,     0,    39,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
      40,     0,     0,     0,     0,    41,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,    42,    43,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,    44,    45,     0,     0,     0,     0,
      46,     0,     0,    47,     0,     0,     0,     0,     0,     0,
       0,    48
};

static const yytype_int16 yycheck[] =
{
       0,   150,   133,   223,   102,     8,   118,   152,     8,     9,
     150,   150,    12,    13,   524,   124,   100,   150,   102,   635,
     150,   136,    98,    99,     4,   534,    16,   150,   614,     8,
      14,    62,   153,    74,   620,     8,   774,   158,    13,   123,
     161,   162,   163,    89,    12,  1086,    29,    14,    52,     8,
      50,   127,   517,    90,    84,   112,   132,     8,    38,    60,
     101,   218,   339,    52,    98,  1133,   150,   635,   112,   106,
    1226,  1227,   111,    79,   483,  1147,   399,  1094,   857,  1116,
    1097,   264,   110,   153,   242,    73,   261,    64,    65,   264,
     147,  1111,   117,     8,    97,   242,   261,   231,   281,    82,
     100,   131,   102,   147,   279,   242,   242,   198,   261,   200,
     110,   111,   270,    98,   242,   145,   865,   117,   118,   276,
     120,   121,   122,   123,   124,   171,   157,   261,   240,   265,
     277,   869,   263,   133,   543,   135,   136,   875,    97,   276,
     277,   141,   758,   279,   235,   151,    97,   763,   276,  1186,
     150,    27,   152,   812,   740,   274,  1176,   242,  1175,   278,
      44,   192,   142,   143,   265,  1206,   154,   150,  1185,   203,
     579,   228,   242,    49,   259,   761,   215,   160,   161,   280,
     217,   767,   961,   162,   228,   276,   261,   190,   227,   162,
     758,  1259,   220,   172,     3,   196,    72,     6,   475,   172,
       9,    47,   233,   234,   527,   242,  1362,   277,   673,   189,
      19,   236,    21,    22,   111,   261,   257,    26,   203,   276,
     220,    30,   222,   223,   973,   974,   975,   227,   261,   344,
     339,   190,   276,   217,   332,    81,   236,    46,   242,   190,
     240,   275,   217,  1315,   262,   276,   276,   162,   332,   198,
     262,   200,   242,   242,  1295,   238,   236,   172,   242,   271,
     278,   345,   108,   263,   923,    79,   925,   242,   854,   249,
     856,    80,   858,   276,   242,   242,   242,   260,   362,    88,
     270,   897,     3,   159,    93,     6,   235,   276,   275,   798,
     275,   276,   276,   408,   409,   881,   276,   276,    19,   108,
      21,   276,   270,   276,   242,   242,   261,   116,   528,    30,
     530,   187,   121,   122,   123,   242,   261,   276,   215,   195,
    1058,   242,    43,   438,   133,    46,   242,   276,   137,   242,
     227,   202,   332,   209,   242,   211,   242,   151,   276,   339,
     186,   417,   279,   270,   344,   345,   274,   218,   276,   276,
     426,   112,   968,   275,   430,   431,   277,   270,   442,    80,
     276,   170,   362,   276,   142,   143,   475,    88,   276,    90,
     276,   180,   181,    84,   220,   261,    48,   186,   185,   255,
     264,   265,   266,   267,   268,   106,   147,   108,   197,   112,
     199,   237,   201,   242,   275,   204,   242,   397,   398,    71,
     121,   122,   123,   275,   213,   142,   143,   216,   408,   409,
     410,   921,   133,    18,   223,    15,   137,   185,   928,   265,
     131,   270,   261,   232,   147,   242,    31,   276,    28,   569,
     569,   277,    82,   242,   145,   242,   569,    37,   438,   569,
     266,   267,   442,   221,   142,   143,   569,   119,   242,   170,
     264,   265,   266,   267,   268,   264,   265,    18,   275,   180,
     181,   277,    67,    68,    25,   186,    66,   228,   277,   276,
     277,   249,   242,   693,   242,   475,   277,   274,   199,    27,
     201,   278,   276,   204,   221,   569,    86,    87,   242,   217,
     218,   831,   213,   242,   120,   216,   217,   242,   276,   261,
     150,    49,   223,   615,   271,   228,   276,   274,    18,   277,
     160,   161,   249,   113,   242,    25,   188,   517,   265,   145,
     665,   242,   276,   221,    72,    50,    51,   276,   528,   613,
     530,   276,   279,   274,   534,   277,   208,   278,   210,   276,
     264,   141,   542,   264,   265,   261,   270,   152,   768,   769,
     275,   249,   271,   261,  1292,   274,   277,   275,   279,   280,
     264,   265,   266,  1149,   274,   261,   166,   261,   278,   569,
     242,   716,   215,   718,   278,   175,   176,   177,   276,   261,
    1238,  1239,    44,   271,   184,   585,   274,   271,   238,  1247,
     274,   191,   271,   426,  1252,   274,   816,   430,   431,   273,
     274,    63,   942,   943,   261,   277,   272,   273,   274,   275,
     260,   159,   261,   613,   614,   615,   271,   762,   277,   274,
     620,   264,   265,   266,   267,   268,   226,   776,   749,   266,
     267,   268,   777,   680,    96,   682,   776,   121,   122,   187,
    1256,   911,   912,   243,   261,   757,   776,   195,   248,   214,
     277,   251,   261,   776,   264,   265,   266,   267,   120,  1218,
    1219,   209,   882,   211,   261,   665,   799,   275,   278,   799,
      82,  1276,  1277,   673,   261,   277,   799,   261,   140,   261,
     261,  1021,   144,   275,  1270,  1025,   275,   275,  1028,  1029,
    1030,  1031,     0,   693,   261,    89,   158,   261,   171,   264,
     265,   266,   267,   268,   264,   265,   266,   255,   163,   849,
      85,   276,  1298,   275,   277,   799,   716,   857,   718,   849,
     275,   275,   275,   262,    15,   277,   849,   857,   277,   874,
     277,   193,   194,   277,   857,   277,   277,    28,   150,   277,
     740,   277,   277,   205,   206,   829,    37,   277,   160,   161,
     895,   277,   275,   829,   830,   277,   277,   757,   277,   277,
     277,   761,   762,   277,   277,   277,   159,   767,   768,   769,
     278,   771,   274,   942,   943,    66,   776,   777,   911,   912,
     274,   911,   912,   214,   270,    82,   275,   261,   911,   912,
     252,   261,   261,    90,   261,    86,    87,  1307,   798,   799,
     261,     7,   261,   261,    10,    11,   278,   261,   261,   106,
     222,   811,   812,   261,   261,   261,   816,    23,    24,   261,
     261,   961,   113,   261,   261,   261,   238,   911,   912,   829,
     274,   961,   261,   264,   265,   266,   267,   268,   961,   138,
     262,   253,   264,   265,   266,   267,   268,   847,   260,   849,
     141,   261,  1021,   150,   854,   152,   856,   857,   858,  1028,
    1029,  1030,  1031,   160,   161,   865,   261,   951,   231,    75,
     145,   261,   261,  1006,   874,   166,  1006,   264,   265,   266,
     267,   881,   882,  1006,   175,   176,   177,   261,    94,    95,
     261,   278,   261,   184,   261,   895,   261,   261,    45,   261,
     191,   264,   265,   266,   267,   268,  1051,   261,  1053,   261,
     910,   911,   912,   231,   261,  1135,   261,   261,   918,   261,
     217,   261,  1006,   923,   261,   925,   261,    76,    77,    78,
     262,   931,   244,   156,   261,   226,   264,   265,   266,   267,
     268,   238,   261,  1027,   275,   242,   264,   265,   266,   267,
     268,   951,   243,   102,   278,   104,   105,   248,   261,   270,
     251,   961,   245,   260,   264,   265,   266,   267,   268,   261,
     274,   261,   171,   973,   974,   975,   447,   448,   278,   450,
     202,   452,   275,   454,   455,   276,   457,   458,   459,   460,
     461,   462,   463,   264,   265,   266,   267,   268,   225,   268,
      70,   261,  1002,   261,   261,  1103,  1006,   278,   415,   416,
     224,   418,   419,   420,   421,   422,   423,   265,   425,  1103,
     427,   428,   429,   261,   261,   261,   261,  1027,   261,   278,
     271,   261,   101,   261,   103,   104,   261,   261,   264,   265,
     266,   267,   268,   246,  1175,   264,   265,   266,   267,   268,
     119,  1051,   278,  1053,   261,   274,   125,   126,   256,   128,
     129,   130,   131,  1175,   261,   275,  1066,   261,   137,   138,
     261,   140,   261,   142,   143,   144,   264,   265,   266,   267,
     268,   261,    17,   264,   265,   266,   267,   268,  1233,   275,
     278,   275,   261,  1224,  1094,  1226,  1227,  1097,  1243,   261,
     261,  1246,   275,  1103,    39,   261,    41,    42,   275,   218,
    1255,   124,  1224,   265,  1226,  1227,   275,  1337,   275,   261,
    1265,   271,  1342,   261,   261,   261,    61,   261,   254,   264,
    1130,   274,   278,   261,   261,  1135,   205,   206,   207,   126,
    1140,  1361,   264,   265,   266,   267,   268,   125,  1148,  1149,
     264,   265,   266,   267,   268,   261,   278,   264,   265,   266,
     267,   268,   276,   276,    99,   100,   264,   265,   266,   267,
     268,   278,   261,   270,   275,  1175,   263,   261,   261,   261,
     278,   261,   117,   118,   275,  1185,   261,   261,   264,   265,
     266,   267,   268,   264,   265,   266,   267,   268,   261,   134,
     135,   136,   278,   275,   219,   261,   270,   278,   275,   277,
     261,   146,   261,   148,   149,   276,   261,    20,   264,   265,
     266,   267,   268,  1223,  1224,   275,  1226,  1227,   275,   164,
     165,  1362,   278,  1233,   264,   265,   266,   267,   268,   274,
     182,   218,   275,  1243,   270,   270,  1246,   275,   278,   266,
    1362,   275,   275,   322,   275,  1255,   275,   264,   265,   266,
     267,   268,   261,   275,   261,  1265,  1266,   270,   261,   271,
    1270,   278,   264,   265,   266,   267,   268,   346,   347,   348,
     349,   350,   351,   352,   353,   354,   355,   356,   357,   358,
     359,   360,   361,   261,   363,   261,   261,   279,  1298,   279,
     369,   370,   371,   372,   373,   240,   375,   264,   265,   266,
     267,   268,   264,   265,   266,   267,   268,  1317,   279,   275,
     389,   278,   261,   261,   275,   270,   278,   264,   265,   266,
     267,   268,  1332,   276,   270,   275,    50,  1337,   271,   276,
     274,   276,  1342,   275,   274,   274,  1346,   271,   279,   261,
     276,   420,  1352,   281,   281,   424,   281,   274,   281,  1359,
     274,  1361,  1362,   261,   261,   274,   261,   436,   264,   265,
     266,   267,   268,     3,   275,   263,     6,   263,   274,     9,
     449,   262,   451,    17,   453,   270,   276,   456,  1152,    19,
     274,    21,    22,   391,   612,   391,    26,   466,   320,   141,
      30,   283,   216,   787,  1071,    39,  1332,    41,    42,   264,
     265,   266,   267,   268,  1346,  1073,    46,   264,   265,   266,
     267,   268,  1357,   278,  1223,   931,   541,    61,   227,  1006,
     397,   278,   798,   124,   475,   264,   265,   266,   267,   268,
     619,   343,   882,   512,   513,   761,    76,    77,    78,   278,
      80,  1124,   895,   264,   265,   266,   267,   268,    88,  1262,
    1130,   874,   240,    93,   542,    99,   100,   278,   263,   722,
    1361,  1359,   102,  1097,   104,   105,  1045,  1186,   108,  1025,
     841,  1255,  1182,   117,   118,  1176,   116,  1246,   847,  1172,
     444,   121,   122,   123,   264,   265,   266,   267,   268,    -1,
     134,   135,   136,   133,  1037,    -1,   276,   137,    -1,    -1,
      -1,    -1,   146,    -1,   148,   149,   264,   265,   266,   267,
     268,   264,   265,   266,   267,   268,   274,    -1,   271,    -1,
     164,   165,   264,   265,   266,   267,   268,    -1,    -1,    -1,
     170,    -1,   274,   612,    -1,   264,   265,   266,   267,   268,
     180,   181,   271,    -1,    -1,    -1,   186,   263,   264,   265,
     266,   267,   268,    -1,    -1,    -1,    -1,   197,    -1,   199,
      -1,   201,    -1,    -1,   204,   263,   264,   265,   266,   267,
     268,    -1,    -1,   213,    -1,    -1,   216,   656,   657,   658,
     659,    -1,    -1,   223,    -1,   664,    -1,    -1,    -1,    -1,
     669,   670,   232,    -1,    -1,    -1,   240,   264,   265,   266,
     267,   268,   242,   263,   264,   265,   266,   267,   687,   688,
     689,   690,   691,    -1,    -1,    -1,    -1,     3,    -1,    -1,
       6,    -1,    -1,    -1,   264,   265,   263,   264,   265,   266,
     267,   268,    -1,    19,    -1,    21,    -1,   277,   717,    -1,
      -1,    -1,    -1,    -1,    30,   263,   264,   265,   266,   267,
     268,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      46,    -1,    -1,    -1,   743,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,   764,   765,    -1,    -1,    -1,
      -1,   770,    -1,   772,    80,    -1,    82,    -1,    -1,    -1,
      -1,    -1,    88,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,   108,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,   813,   121,   122,   123,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,   133,    -1,    -1,
      -1,   137,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,     3,    -1,   150,     6,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,   160,   161,    -1,    -1,    19,    -1,
      21,    -1,    -1,    -1,   170,    -1,    -1,    -1,    -1,    30,
      -1,    -1,    -1,    -1,   180,   181,    -1,    -1,    -1,    -1,
     186,    -1,    -1,    -1,    -1,    46,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,   199,    -1,   201,    -1,     3,   204,    -1,
       6,    -1,    -1,   902,   903,    -1,    -1,   213,   907,    -1,
     216,    -1,    -1,    19,    -1,    21,    -1,   223,    -1,    80,
      -1,   920,    -1,    -1,    30,    -1,    -1,    88,    -1,    90,
      -1,    -1,   238,   932,    -1,    -1,   242,    43,    -1,    -1,
      46,    -1,    -1,    -1,    -1,   106,    -1,   108,    -1,    -1,
      -1,    -1,    -1,    -1,   260,    -1,    -1,    -1,   264,   265,
     121,   122,   123,    -1,    -1,    -1,    -1,   966,    -1,    -1,
      -1,   277,   133,    -1,    80,    -1,   137,    -1,    -1,    -1,
      -1,   980,    88,    -1,    90,   984,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
     106,    -1,   108,    -1,    -1,    -1,    -1,    -1,    -1,   170,
      -1,  1010,    -1,    -1,    -1,   121,   122,   123,    -1,   180,
     181,  1020,    -1,    -1,    -1,   186,    -1,   133,    -1,    -1,
      -1,   137,    -1,    -1,    -1,    -1,    -1,    -1,   199,    -1,
     201,    -1,    -1,   204,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,   213,  1052,    -1,   216,   217,    -1,    -1,    -1,
    1059,    -1,   223,    -1,   170,    -1,    -1,    -1,    -1,    -1,
    1069,    -1,  1071,    -1,   180,   181,    -1,    -1,    -1,    -1,
     186,   242,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,   199,    -1,   201,    -1,    -1,   204,    -1,
      -1,    -1,    -1,   264,   265,    -1,    -1,   213,    -1,   270,
     216,   217,     3,    -1,    -1,     6,   277,   223,    -1,    -1,
      -1,    -1,    -1,    -1,  1123,    -1,    -1,    -1,    19,    -1,
      21,    -1,    -1,    -1,  1133,    -1,   242,    -1,  1137,    30,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,     3,    -1,    -1,     6,    46,    -1,    -1,   264,   265,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    19,    -1,    21,
      -1,   277,    -1,    -1,    -1,    -1,    -1,    -1,    30,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    80,
      -1,    -1,  1191,  1192,    46,    -1,    -1,    88,    -1,    90,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,   106,    -1,   108,    -1,  1218,
    1219,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    80,    -1,
     121,   122,   123,    -1,    -1,    -1,    88,    -1,    -1,    -1,
      -1,    -1,   133,    -1,    -1,    -1,   137,    -1,    -1,    -1,
       3,    -1,    -1,     6,    -1,    -1,   108,   109,    -1,    -1,
      -1,  1260,    -1,  1262,    -1,  1264,    19,    -1,    21,   121,
     122,   123,    -1,    -1,    -1,    -1,    -1,    30,    -1,   170,
      -1,   133,  1281,    -1,    -1,   137,    -1,    -1,    -1,   180,
     181,    -1,    -1,    46,    -1,   186,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   199,    -1,
     201,    -1,    -1,   204,    -1,  1314,    -1,    -1,   170,    -1,
      -1,    -1,   213,    -1,    -1,   216,   217,    80,   180,   181,
      -1,    -1,   223,    -1,   186,    88,    -1,    -1,     3,    -1,
      -1,     6,    -1,    -1,    -1,    -1,  1345,   199,    -1,   201,
      -1,   242,   204,    -1,    19,   108,    21,    -1,    -1,    -1,
      -1,   213,    -1,    -1,   216,    30,    -1,    -1,   121,   122,
     123,   223,    -1,   264,   265,    -1,    -1,    -1,    -1,    -1,
     133,    46,    -1,    -1,   137,    -1,   277,    -1,    -1,    -1,
     242,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,   264,   265,    -1,    80,    -1,   170,    -1,    -1,
      -1,    -1,    -1,    88,    -1,   277,    -1,   180,   181,    -1,
      -1,    -1,    -1,   186,    -1,    -1,    -1,     3,    -1,    -1,
       6,    -1,    -1,   108,    -1,    -1,   199,    -1,   201,    -1,
      -1,   204,    -1,    19,    -1,    21,   121,   122,   123,    -1,
     213,    -1,    -1,   216,    30,    -1,    -1,    -1,   133,    -1,
     223,    -1,   137,    -1,    -1,    -1,     3,    -1,    -1,     6,
      46,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   242,
      -1,    -1,    19,    -1,    21,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    30,    -1,   170,    -1,    -1,    -1,    -1,
      -1,   264,   265,    -1,    80,   180,   181,   270,    -1,    46,
      -1,   186,    88,    -1,   277,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,   199,    -1,   201,    -1,    -1,   204,
      -1,    -1,   108,    -1,    -1,    -1,    -1,    -1,   213,    -1,
      -1,   216,    -1,    80,    -1,   121,   122,   123,   223,    -1,
      -1,    88,    -1,    -1,    -1,    -1,    -1,   133,    -1,    -1,
      -1,   137,    -1,    -1,   239,    -1,    -1,   242,    -1,    -1,
      -1,   108,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,   121,   122,   123,    -1,    -1,   264,
     265,    -1,    -1,    -1,   170,    -1,   133,    -1,    -1,    -1,
     137,    -1,   277,    -1,   180,   181,    -1,    -1,    -1,    -1,
     186,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,   199,    -1,   201,    -1,    -1,   204,    -1,
      -1,    -1,    -1,   170,    -1,    -1,    -1,   213,    -1,    -1,
     216,    -1,    -1,   180,   181,    -1,    -1,   223,    -1,   186,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,   199,    -1,   201,    -1,   242,   204,     5,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,   213,    -1,    -1,   216,
      -1,    -1,    -1,    -1,    -1,    -1,   223,    -1,   264,   265,
      -1,    -1,    -1,    -1,   270,    32,    33,    34,    35,    36,
      -1,   277,    -1,    40,    -1,   242,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    53,    54,    55,    56,
      57,    58,    59,    60,    -1,    -1,    -1,   264,   265,    -1,
      -1,    -1,    69,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
     277,    -1,    -1,    -1,    -1,    -1,    83,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    91,    92,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,   103,    -1,    -1,    -1,
     107,    -1,    -1,   110,    -1,    -1,    -1,   114,   115,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
     127,   128,   129,   130,    -1,   132,    -1,    -1,    -1,    -1,
      -1,    -1,   139,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   155,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
     167,   168,   169,    -1,    -1,    -1,   173,   174,    -1,    -1,
      -1,   178,   179,    -1,    -1,    -1,   183,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
     207,    -1,    -1,    -1,    -1,   212,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,   229,   230,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,   241,   242,    -1,    -1,    -1,    -1,
     247,    -1,    -1,   250,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,   258
};

  /* YYSTOS[STATE-NUM] -- The (internal number of the) accessing
     symbol of state STATE-NUM.  */
static const yytype_uint16 yystos[] =
{
       0,     5,    32,    33,    34,    35,    36,    40,    53,    54,
      55,    56,    57,    58,    59,    60,    69,    83,    91,    92,
     103,   107,   110,   114,   115,   127,   128,   129,   130,   132,
     139,   155,   167,   168,   169,   173,   174,   178,   179,   183,
     207,   212,   229,   230,   241,   242,   247,   250,   258,   283,
     284,   285,   287,   302,   303,   304,   318,   319,   320,   322,
     327,   328,   329,   330,   331,   332,   337,   341,   344,   345,
     346,   347,   348,   349,   350,   359,   389,   390,   391,   407,
     434,   439,   441,   442,   443,   449,   454,   455,   459,   473,
     474,   504,   509,   513,   524,   528,   568,   610,   261,   261,
     261,   261,   261,   261,   261,   287,   363,   287,   352,   353,
     275,   275,   287,   435,   287,   354,   393,   275,   275,   261,
     277,   277,   277,   261,   440,   261,   261,   261,   261,   261,
     261,   261,   261,   275,   275,   277,   277,   261,   261,   275,
     261,   277,   261,   261,   261,   275,   275,   275,     0,   285,
     261,    89,   261,    27,    49,    72,   159,   187,   195,   209,
     211,   255,   461,   171,    82,   150,   160,   161,   238,   260,
     294,   294,    90,   106,   217,   286,   287,   288,   315,   316,
     317,     3,     6,    19,    21,    30,    46,    80,    88,   108,
     121,   122,   123,   133,   137,   170,   180,   181,   186,   199,
     201,   204,   213,   216,   223,   264,   265,   277,   287,   309,
     310,   311,   313,   314,   288,   309,   343,   309,   360,   275,
     351,   352,   185,   277,   287,   385,   387,   408,   409,   410,
     421,   422,   423,   424,   275,   275,   392,   393,   287,   289,
     514,   515,   309,   287,   323,   323,   323,   315,   287,   443,
     449,   309,   309,   294,   309,   309,   309,   309,   222,   253,
     294,   289,   291,   525,   526,    15,    28,    37,    66,    86,
      87,   113,   141,   166,   175,   176,   177,   184,   191,   226,
     243,   248,   251,   333,   334,   315,   324,   324,   109,   309,
     309,   163,   530,   309,   287,   304,   309,   309,   309,   569,
      85,   611,    17,    39,    41,    42,    61,    99,   100,   117,
     118,   134,   135,   136,   146,   148,   149,   164,   165,   240,
     338,   339,   270,   287,   305,   307,   309,   311,   312,   314,
     315,   316,   277,   287,   306,   307,   308,   275,   444,   444,
     444,   444,   444,   444,   277,   262,   277,   277,   277,   277,
     277,   277,   277,   277,   277,   277,   277,   277,   277,   277,
     277,   277,   277,   277,   309,   309,   309,    79,   151,   264,
     265,   266,   267,   268,   342,   263,   342,   275,    64,    65,
     355,   276,   352,   287,   411,   387,   272,   273,   274,   275,
     295,   296,   297,   298,   299,   276,   409,   264,   281,   425,
     437,   394,   276,   393,   516,   276,   515,   278,   274,   274,
     270,   275,   276,   526,   261,   261,   261,   261,   261,   261,
     261,   261,   261,   261,   261,   261,   261,   261,   261,   261,
     261,   261,   276,   334,   325,   278,   261,   529,   274,   138,
     571,   572,   261,   145,   612,   613,   261,   261,   261,   261,
     261,   261,   261,   261,   261,   261,   261,   261,   261,   261,
     261,   261,   261,   261,   276,   339,   270,   300,   301,   309,
     288,    45,   287,   442,   443,   450,   451,   452,   456,   244,
     476,   460,   476,   324,   315,   309,   309,   309,   309,   309,
     309,   309,   309,   309,   309,   309,   309,   309,   309,   309,
     309,   315,   309,   278,   309,   309,   309,   309,   309,   309,
     156,   364,   261,   261,    50,    51,   356,   275,   278,   309,
     297,   298,   423,   287,   395,   426,   265,   280,   413,   414,
     415,   416,   417,   418,   419,    62,   157,   192,   233,   234,
     276,   396,   270,   324,   324,   287,   527,   294,   335,   335,
     335,    29,   294,   336,   335,   335,   309,   335,   335,   335,
     335,   309,   335,   336,   335,   335,   335,   336,   336,   274,
     278,   309,   111,   215,   227,   531,   532,   533,   534,   324,
     261,   245,   573,   574,   315,   261,   120,   613,   616,    74,
     101,   257,   340,   340,   340,   309,   340,   309,   340,   309,
     340,   340,   309,   340,   340,   340,   340,   340,   340,   340,
     309,   271,   274,   274,   261,   159,   446,   451,   202,   457,
     275,   225,   506,    44,    63,    96,   120,   140,   144,   158,
     193,   194,   205,   206,   252,   462,   464,   465,   466,   467,
     468,   469,   470,    70,   479,   325,   278,   278,   278,   278,
     278,   278,   278,   278,   278,   278,   274,   274,   274,   274,
     278,   278,   278,   278,   263,   261,   361,   309,   309,   261,
     261,   224,   357,   412,   421,   276,   295,   413,   265,   279,
     387,   279,   387,   153,   424,   427,   428,   261,   261,   261,
     261,   261,     4,    38,   142,   143,   189,   236,   249,   276,
     397,   398,   399,   401,   402,   406,   287,   519,   325,   278,
     271,   520,   305,   312,   316,   326,   261,   261,   261,    98,
     203,   275,   535,   536,   537,   540,   541,   325,    18,    31,
      67,    68,   152,   261,   246,   575,   576,   287,   614,   615,
     261,   256,   617,   231,   301,   315,   321,   292,   306,   289,
     198,   200,   235,   276,   445,   447,   448,   261,   462,   292,
     477,   478,   275,   505,   261,   261,   275,   261,   261,   275,
     261,   261,   261,   275,   446,   464,   261,   275,   475,   278,
     309,   309,   309,   309,   309,   306,   218,   365,   366,   309,
     309,   124,   358,   276,   421,   413,   414,   414,   264,   270,
     429,   430,   431,   309,   239,   309,   438,   309,   309,   309,
     387,   275,   275,   261,   231,   261,   400,   271,   278,   221,
     249,   276,   402,   521,   522,   523,   306,   309,   306,   261,
     261,   538,   276,   536,   278,    18,    25,   261,    84,   131,
     145,   570,   577,   578,   579,   580,   597,   264,   292,   261,
     254,   618,   309,   278,   274,   444,   261,   261,   261,    48,
      71,   119,   188,   208,   210,   277,   289,   291,   458,   446,
     276,   477,   306,   507,   508,   446,   309,   309,   126,   292,
     387,   387,   471,   472,   309,   287,   436,   309,   125,   445,
     287,   307,   312,   463,   306,   480,   481,   482,   278,   278,
     278,   278,   263,   270,   202,   366,   369,   261,   276,   428,
      43,   279,   280,   305,   312,   316,   432,   433,   263,   420,
     261,   287,   386,   403,   405,   404,   405,   309,    14,   387,
     517,   261,   261,   152,   294,   315,   539,   294,    47,    81,
     220,   237,   265,   277,   310,   542,   543,   544,   545,    18,
      25,   261,   275,   275,   276,   577,   614,   287,   292,   293,
     312,   261,   111,   215,   227,   619,   215,   292,   453,   292,
     293,   292,   458,   264,   265,   266,   445,   276,   507,   445,
     261,   292,   276,   472,   261,   276,   481,    60,   196,   446,
     483,   484,   499,   309,   309,   367,   368,   275,   219,   373,
     309,   287,   388,   433,   433,   271,   274,   287,   309,   295,
     261,   276,   405,   276,   405,   295,   275,   395,   309,   546,
     270,   277,   547,   544,   544,   274,   276,   263,   264,   265,
     266,   267,   567,   315,   112,   147,   228,   598,   599,   600,
     603,   606,   112,   147,   228,   581,   582,   583,   589,   593,
     293,   261,   261,   261,   276,   309,    20,   512,   446,   274,
     278,   458,   458,   458,   309,   309,   275,   275,   276,   263,
     271,   274,   218,   370,   372,   275,   362,   287,   432,   309,
      76,    77,    78,   102,   104,   105,   488,   489,   490,   491,
     495,   496,   497,   498,   270,   309,   544,   270,   278,   543,
     261,   264,   279,   548,   315,   544,   544,   544,   544,   275,
     275,   275,   276,   599,   275,   275,   275,   276,   582,   306,
     309,   306,   271,   261,   510,   445,   309,   214,   276,   214,
     500,   501,   519,   485,   309,   261,   368,   270,   276,   372,
     374,   182,   376,   271,   261,   261,   489,   518,   261,   261,
     218,   286,   287,   550,   551,   552,   553,   554,   555,   556,
     558,   271,   278,   552,   279,   279,   279,   288,   549,     8,
      97,   190,   607,   608,   609,   601,   604,   605,   609,     8,
     162,   172,   594,   595,   596,   585,   590,   591,   596,   309,
     482,   261,   261,   276,   501,   275,     9,    22,    26,    93,
     116,   197,   232,   309,   486,   487,   488,   492,   493,   494,
     387,   309,   371,   276,   287,   375,   275,   276,   270,   270,
     520,   287,   292,   275,   274,   299,   274,   274,   271,   271,
     279,   276,   608,   281,    13,   276,   286,   289,   290,   291,
     602,   276,   605,   281,   276,   595,   281,    14,   276,   286,
     287,   586,   587,   276,   591,   281,   511,   309,   309,   502,
     231,   276,   274,   271,   274,   261,   377,   493,   493,   276,
     274,   388,   259,   289,   291,   564,   564,   564,    12,   306,
     592,   261,   584,   584,   584,    16,   306,   588,   592,   584,
     584,   588,   446,   276,    44,   488,   309,   494,   261,   309,
     306,   276,   287,   378,   271,   271,   292,   559,   274,   565,
     565,   565,   309,   445,   261,   503,   292,   275,   295,     7,
      10,    11,    23,    24,    75,    94,    95,   566,   309,   520,
      52,   287,   379,   380,   381,   382,   383,   263,   276,   276,
      73,   154,   384,   276,   380,   263,   262,   387,   557,   387,
     309,   383,   270,   560,   287,   561,   562,   563,   276,   262,
     271,   384,   274,   562,   557,   564
};

  /* YYR1[YYN] -- Symbol number of symbol that rule YYN derives.  */
static const yytype_uint16 yyr1[] =
{
       0,   282,   283,   284,   284,   285,   285,   285,   285,   285,
     285,   285,   285,   285,   285,   285,   285,   285,   285,   285,
     285,   285,   285,   285,   285,   286,   287,   288,   289,   290,
     290,   291,   292,   293,   293,   294,   294,   294,   294,   294,
     294,   295,   295,   295,   295,   296,   296,   296,   296,   297,
     298,   299,   300,   300,   301,   301,   302,   303,   303,   303,
     303,   304,   305,   306,   306,   307,   308,   309,   309,   310,
     310,   311,   311,   312,   312,   313,   314,   314,   314,   314,
     314,   314,   314,   314,   314,   314,   314,   314,   314,   314,
     314,   314,   314,   314,   314,   314,   314,   314,   314,   314,
     314,   314,   314,   314,   314,   314,   315,   315,   316,   316,
     316,   316,   317,   318,   318,   318,   318,   318,   318,   318,
     319,   320,   321,   322,   323,   324,   325,   325,   326,   326,
     326,   327,   328,   329,   330,   331,   332,   333,   333,   334,
     334,   334,   334,   334,   334,   334,   334,   334,   334,   334,
     334,   334,   334,   334,   334,   334,   334,   334,   335,   336,
     336,   337,   338,   338,   339,   339,   339,   339,   339,   339,
     339,   339,   339,   339,   339,   339,   339,   339,   339,   339,
     339,   339,   340,   340,   340,   341,   341,   341,   341,   341,
     342,   342,   342,   343,   343,   343,   343,   344,   344,   344,
     344,   344,   344,   344,   344,   344,   344,   344,   344,   344,
     344,   344,   344,   345,   345,   345,   345,   346,   347,   347,
     347,   348,   348,   348,   349,   350,   351,   351,   352,   353,
     354,   355,   355,   356,   356,   356,   357,   357,   358,   358,
     360,   361,   362,   359,   363,   364,   365,   365,   366,   367,
     367,   368,   368,   369,   370,   370,   371,   371,   372,   373,
     374,   374,   375,   376,   377,   377,   378,   379,   379,   380,
     381,   381,   382,   382,   383,   384,   384,   385,   386,   387,
     388,   389,   389,   390,   391,   392,   392,   394,   393,   395,
     396,   396,   397,   397,   397,   398,   398,   398,   399,   399,
     399,   400,   400,   401,   402,   402,   403,   403,   404,   404,
     405,   406,   407,   408,   408,   409,   409,   410,   411,   412,
     412,   413,   413,   414,   415,   416,   417,   418,   419,   419,
     419,   419,   420,   420,   421,   422,   422,   423,   424,   424,
     425,   425,   426,   427,   427,   428,   428,   429,   429,   430,
     431,   432,   432,   432,   433,   433,   433,   433,   434,   435,
     436,   437,   437,   437,   437,   437,   437,   438,   438,   440,
     439,   441,   442,   442,   442,   442,   442,   443,   444,   445,
     446,   446,   447,   448,   448,   448,   449,   450,   450,   451,
     451,   453,   452,   454,   454,   456,   455,   457,   457,   457,
     457,   457,   457,   457,   457,   458,   458,   458,   458,   458,
     460,   459,   461,   461,   461,   461,   461,   462,   462,   463,
     464,   464,   464,   464,   464,   464,   464,   464,   464,   465,
     465,   466,   466,   466,   466,   467,   467,   468,   469,   470,
     470,   471,   471,   472,   473,   475,   474,   476,   477,   478,
     478,   479,   480,   480,   481,   482,   482,   483,   483,   485,
     484,   486,   486,   487,   487,   487,   487,   487,   487,   487,
     488,   488,   489,   489,   489,   489,   490,   491,   492,   493,
     493,   494,   494,   494,   495,   496,   496,   497,   498,   498,
     499,   500,   500,   502,   503,   501,   505,   504,   506,   507,
     508,   508,   510,   511,   509,   512,   512,   513,   514,   514,
     516,   517,   518,   515,   519,   520,   520,   521,   521,   521,
     522,   523,   524,   525,   525,   527,   526,   529,   528,   530,
     530,   531,   531,   531,   532,   533,   534,   535,   535,   536,
     536,   536,   538,   537,   539,   539,   539,   540,   541,   542,
     542,   543,   544,   544,   544,   544,   544,   544,   544,   544,
     544,   546,   545,   545,   547,   545,   548,   548,   548,   548,
     548,   549,   550,   551,   552,   552,   552,   552,   553,   554,
     555,   556,   557,   559,   558,   560,   560,   561,   561,   562,
     563,   564,   564,   564,   565,   565,   566,   566,   566,   566,
     566,   566,   566,   566,   567,   567,   569,   568,   570,   570,
     571,   571,   572,   572,   572,   572,   572,   573,   573,   574,
     574,   575,   575,   576,   576,   577,   577,   578,   578,   579,
     580,   581,   581,   582,   582,   582,   583,   584,   584,   585,
     585,   586,   586,   587,   587,   588,   588,   589,   590,   590,
     591,   592,   592,   593,   594,   594,   595,   596,   596,   596,
     597,   598,   598,   599,   599,   599,   600,   601,   601,   602,
     602,   602,   603,   604,   604,   605,   606,   607,   607,   608,
     609,   609,   609,   610,   611,   612,   612,   613,   614,   615,
     615,   616,   617,   618,   619,   619,   619,   619
};

  /* YYR2[YYN] -- Number of symbols on the right hand side of rule YYN.  */
static const yytype_uint8 yyr2[] =
{
       0,     2,     1,     1,     2,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     4,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     0,     1,     1,     1,     1,     1,     2,     2,     1,
       1,     3,     1,     3,     1,     7,     3,     3,     3,     3,
       3,     1,     1,     1,     1,     3,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     3,     4,     4,     4,
       6,     6,     6,     4,     4,     4,     4,     4,     4,     4,
       4,     4,     4,     6,     1,     1,     1,     1,     4,     3,
       3,     3,     3,     3,     2,     2,     1,     1,     1,     1,
       3,     5,     1,     1,     1,     1,     1,     1,     1,     1,
       7,     1,     1,     4,     1,     1,     0,     3,     1,     1,
       1,     5,     7,     7,     4,     6,     4,     1,     2,     3,
       3,     3,     3,     3,     3,     3,     3,     3,     3,     3,
       3,     3,     3,     3,     3,     3,     3,     3,     1,     1,
       1,     4,     1,     2,     3,     3,     3,     3,     3,     3,
       3,     3,     3,     3,     3,     3,     3,     3,     3,     3,
       3,     3,     1,     1,     1,     3,     3,     4,     3,     4,
       0,     1,     1,     1,     3,     5,     7,     3,     3,     3,
       3,     3,     3,     3,     3,     3,     3,     3,     3,     3,
       3,     3,     3,     3,     3,     3,     3,     3,     1,     1,
       1,     1,     1,     1,     2,     4,     1,     2,     7,     1,
       1,     3,     3,     0,     3,     3,     0,     1,     0,     3,
       0,     0,     0,    12,     1,     3,     1,     2,     6,     1,
       3,     1,     3,     4,     1,     2,     1,     3,     6,     4,
       0,     2,     3,     4,     0,     2,     4,     1,     2,     3,
       1,     1,     1,     3,     3,     1,     1,     1,     2,     2,
       1,     1,     1,     2,     4,     1,     2,     0,     5,     1,
       0,     2,     1,     1,     1,     3,     4,     4,     1,     1,
       1,     1,     1,     1,     4,     4,     1,     2,     1,     2,
       3,     3,     4,     1,     2,     1,     1,     5,     1,     1,
       2,     1,     2,     2,     2,     2,     3,     3,     1,     1,
       1,     1,     0,     2,     6,     1,     3,     1,     1,     3,
       0,     2,     2,     1,     3,     1,     1,     1,     1,     3,
       5,     1,     2,     2,     1,     1,     1,     3,     5,     1,
       1,     0,     4,     4,     4,     4,     4,     1,     1,     0,
       3,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       0,     2,     1,     3,     3,     5,     6,     1,     2,     1,
       1,     0,     7,     1,     1,     0,     8,     3,     3,     3,
       3,     3,     3,     3,     3,     1,     3,     3,     3,     3,
       0,     7,     1,     1,     1,     1,     1,     1,     2,     1,
       3,     3,     1,     3,     3,     3,     3,     3,     4,     1,
       1,     1,     1,     1,     1,     3,     6,     9,    12,     3,
       3,     1,     2,     2,     1,     0,     9,     4,     1,     1,
       2,     4,     1,     2,     1,     0,     2,     1,     1,     0,
       5,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     2,     1,     1,     1,     1,     5,     5,     1,     1,
       3,     1,     3,     1,     3,     1,     1,     5,     1,     1,
       4,     1,     2,     0,     0,     7,     0,     8,     4,     1,
       1,     2,     0,     0,    14,     0,     3,     4,     1,     2,
       0,     0,     0,    11,     1,     0,     2,     1,     1,     1,
       3,     3,     4,     1,     2,     0,     5,     0,     7,     0,
       3,     1,     1,     1,     3,     3,     3,     1,     2,     1,
       1,     1,     0,     6,     1,     1,     1,     3,     3,     1,
       3,     2,     1,     1,     3,     3,     3,     3,     3,     2,
       4,     0,     5,     4,     0,     5,     1,     2,     2,     3,
       2,     1,     1,     2,     1,     1,     1,     1,     4,     4,
       4,     1,     1,     0,    11,     0,     3,     1,     3,     3,
       1,     1,     1,     1,     0,     2,     1,     1,     1,     1,
       1,     1,     1,     1,     0,     2,     0,     8,     1,     2,
       0,     1,     3,     3,     3,     3,     3,     0,     1,     3,
       3,     0,     1,     3,     3,     1,     1,     1,     1,     3,
       4,     1,     2,     1,     1,     1,     4,     2,     0,     2,
       0,     2,     2,     1,     1,     1,     1,     4,     1,     2,
       3,     1,     1,     4,     1,     2,     3,     1,     1,     1,
       4,     1,     2,     1,     1,     1,     4,     2,     0,     2,
       2,     2,     4,     1,     2,     3,     4,     1,     2,     3,
       1,     1,     1,     9,     3,     1,     2,     3,     1,     1,
       3,     3,     3,     3,     0,     3,     3,     3
};


#define yyerrok         (yyerrstatus = 0)
#define yyclearin       (yychar = YYEMPTY)
#define YYEMPTY         (-2)
#define YYEOF           0

#define YYACCEPT        goto yyacceptlab
#define YYABORT         goto yyabortlab
#define YYERROR         goto yyerrorlab


#define YYRECOVERING()  (!!yyerrstatus)

#define YYBACKUP(Token, Value)                                  \
do                                                              \
  if (yychar == YYEMPTY)                                        \
    {                                                           \
      yychar = (Token);                                         \
      yylval = (Value);                                         \
      YYPOPSTACK (yylen);                                       \
      yystate = *yyssp;                                         \
      goto yybackup;                                            \
    }                                                           \
  else                                                          \
    {                                                           \
      yyerror (parse_state, scanner, YY_("syntax error: cannot back up")); \
      YYERROR;                                                  \
    }                                                           \
while (0)

/* Error token number */
#define YYTERROR        1
#define YYERRCODE       256



/* Enable debugging if requested.  */
#if YYDEBUG

# ifndef YYFPRINTF
#  include <stdio.h> /* INFRINGES ON USER NAME SPACE */
#  define YYFPRINTF fprintf
# endif

# define YYDPRINTF(Args)                        \
do {                                            \
  if (yydebug)                                  \
    YYFPRINTF Args;                             \
} while (0)

/* This macro is provided for backward compatibility. */
#ifndef YY_LOCATION_PRINT
# define YY_LOCATION_PRINT(File, Loc) ((void) 0)
#endif


# define YY_SYMBOL_PRINT(Title, Type, Value, Location)                    \
do {                                                                      \
  if (yydebug)                                                            \
    {                                                                     \
      YYFPRINTF (stderr, "%s ", Title);                                   \
      yy_symbol_print (stderr,                                            \
                  Type, Value, parse_state, scanner); \
      YYFPRINTF (stderr, "\n");                                           \
    }                                                                     \
} while (0)


/*----------------------------------------.
| Print this symbol's value on YYOUTPUT.  |
`----------------------------------------*/

static void
yy_symbol_value_print (FILE *yyoutput, int yytype, YYSTYPE const * const yyvaluep, struct mdlparse_vars *parse_state, yyscan_t scanner)
{
  FILE *yyo = yyoutput;
  YYUSE (yyo);
  YYUSE (parse_state);
  YYUSE (scanner);
  if (!yyvaluep)
    return;
# ifdef YYPRINT
  if (yytype < YYNTOKENS)
    YYPRINT (yyoutput, yytoknum[yytype], *yyvaluep);
# endif
  YYUSE (yytype);
}


/*--------------------------------.
| Print this symbol on YYOUTPUT.  |
`--------------------------------*/

static void
yy_symbol_print (FILE *yyoutput, int yytype, YYSTYPE const * const yyvaluep, struct mdlparse_vars *parse_state, yyscan_t scanner)
{
  YYFPRINTF (yyoutput, "%s %s (",
             yytype < YYNTOKENS ? "token" : "nterm", yytname[yytype]);

  yy_symbol_value_print (yyoutput, yytype, yyvaluep, parse_state, scanner);
  YYFPRINTF (yyoutput, ")");
}

/*------------------------------------------------------------------.
| yy_stack_print -- Print the state stack from its BOTTOM up to its |
| TOP (included).                                                   |
`------------------------------------------------------------------*/

static void
yy_stack_print (yytype_int16 *yybottom, yytype_int16 *yytop)
{
  YYFPRINTF (stderr, "Stack now");
  for (; yybottom <= yytop; yybottom++)
    {
      int yybot = *yybottom;
      YYFPRINTF (stderr, " %d", yybot);
    }
  YYFPRINTF (stderr, "\n");
}

# define YY_STACK_PRINT(Bottom, Top)                            \
do {                                                            \
  if (yydebug)                                                  \
    yy_stack_print ((Bottom), (Top));                           \
} while (0)


/*------------------------------------------------.
| Report that the YYRULE is going to be reduced.  |
`------------------------------------------------*/

static void
yy_reduce_print (yytype_int16 *yyssp, YYSTYPE *yyvsp, int yyrule, struct mdlparse_vars *parse_state, yyscan_t scanner)
{
  unsigned long int yylno = yyrline[yyrule];
  int yynrhs = yyr2[yyrule];
  int yyi;
  YYFPRINTF (stderr, "Reducing stack by rule %d (line %lu):\n",
             yyrule - 1, yylno);
  /* The symbols being reduced.  */
  for (yyi = 0; yyi < yynrhs; yyi++)
    {
      YYFPRINTF (stderr, "   $%d = ", yyi + 1);
      yy_symbol_print (stderr,
                       yystos[yyssp[yyi + 1 - yynrhs]],
                       &(yyvsp[(yyi + 1) - (yynrhs)])
                                              , parse_state, scanner);
      YYFPRINTF (stderr, "\n");
    }
}

# define YY_REDUCE_PRINT(Rule)          \
do {                                    \
  if (yydebug)                          \
    yy_reduce_print (yyssp, yyvsp, Rule, parse_state, scanner); \
} while (0)

/* Nonzero means print parse trace.  It is left uninitialized so that
   multiple parsers can coexist.  */
int yydebug;
#else /* !YYDEBUG */
# define YYDPRINTF(Args)
# define YY_SYMBOL_PRINT(Title, Type, Value, Location)
# define YY_STACK_PRINT(Bottom, Top)
# define YY_REDUCE_PRINT(Rule)
#endif /* !YYDEBUG */


/* YYINITDEPTH -- initial size of the parser's stacks.  */
#ifndef YYINITDEPTH
# define YYINITDEPTH 200
#endif

/* YYMAXDEPTH -- maximum size the stacks can grow to (effective only
   if the built-in stack extension method is used).

   Do not make this value too large; the results are undefined if
   YYSTACK_ALLOC_MAXIMUM < YYSTACK_BYTES (YYMAXDEPTH)
   evaluated with infinite-precision integer arithmetic.  */

#ifndef YYMAXDEPTH
# define YYMAXDEPTH 10000
#endif


#if YYERROR_VERBOSE

# ifndef yystrlen
#  if defined __GLIBC__ && defined _STRING_H
#   define yystrlen strlen
#  else
/* Return the length of YYSTR.  */
static YYSIZE_T
yystrlen (const char *yystr)
{
  YYSIZE_T yylen;
  for (yylen = 0; yystr[yylen]; yylen++)
    continue;
  return yylen;
}
#  endif
# endif

# ifndef yystpcpy
#  if defined __GLIBC__ && defined _STRING_H && defined _GNU_SOURCE
#   define yystpcpy stpcpy
#  else
/* Copy YYSRC to YYDEST, returning the address of the terminating '\0' in
   YYDEST.  */
static char *
yystpcpy (char *yydest, const char *yysrc)
{
  char *yyd = yydest;
  const char *yys = yysrc;

  while ((*yyd++ = *yys++) != '\0')
    continue;

  return yyd - 1;
}
#  endif
# endif

# ifndef yytnamerr
/* Copy to YYRES the contents of YYSTR after stripping away unnecessary
   quotes and backslashes, so that it's suitable for yyerror.  The
   heuristic is that double-quoting is unnecessary unless the string
   contains an apostrophe, a comma, or backslash (other than
   backslash-backslash).  YYSTR is taken from yytname.  If YYRES is
   null, do not copy; instead, return the length of what the result
   would have been.  */
static YYSIZE_T
yytnamerr (char *yyres, const char *yystr)
{
  if (*yystr == '"')
    {
      YYSIZE_T yyn = 0;
      char const *yyp = yystr;

      for (;;)
        switch (*++yyp)
          {
          case '\'':
          case ',':
            goto do_not_strip_quotes;

          case '\\':
            if (*++yyp != '\\')
              goto do_not_strip_quotes;
            /* Fall through.  */
          default:
            if (yyres)
              yyres[yyn] = *yyp;
            yyn++;
            break;

          case '"':
            if (yyres)
              yyres[yyn] = '\0';
            return yyn;
          }
    do_not_strip_quotes: ;
    }

  if (! yyres)
    return yystrlen (yystr);

  return yystpcpy (yyres, yystr) - yyres;
}
# endif

/* Copy into *YYMSG, which is of size *YYMSG_ALLOC, an error message
   about the unexpected token YYTOKEN for the state stack whose top is
   YYSSP.

   Return 0 if *YYMSG was successfully written.  Return 1 if *YYMSG is
   not large enough to hold the message.  In that case, also set
   *YYMSG_ALLOC to the required number of bytes.  Return 2 if the
   required number of bytes is too large to store.  */
static int
yysyntax_error (YYSIZE_T *yymsg_alloc, char **yymsg,
                yytype_int16 *yyssp, int yytoken)
{
  YYSIZE_T yysize0 = yytnamerr (YY_NULLPTR, yytname[yytoken]);
  YYSIZE_T yysize = yysize0;
  enum { YYERROR_VERBOSE_ARGS_MAXIMUM = 5 };
  /* Internationalized format string. */
  const char *yyformat = YY_NULLPTR;
  /* Arguments of yyformat. */
  char const *yyarg[YYERROR_VERBOSE_ARGS_MAXIMUM];
  /* Number of reported tokens (one for the "unexpected", one per
     "expected"). */
  int yycount = 0;

  /* There are many possibilities here to consider:
     - If this state is a consistent state with a default action, then
       the only way this function was invoked is if the default action
       is an error action.  In that case, don't check for expected
       tokens because there are none.
     - The only way there can be no lookahead present (in yychar) is if
       this state is a consistent state with a default action.  Thus,
       detecting the absence of a lookahead is sufficient to determine
       that there is no unexpected or expected token to report.  In that
       case, just report a simple "syntax error".
     - Don't assume there isn't a lookahead just because this state is a
       consistent state with a default action.  There might have been a
       previous inconsistent state, consistent state with a non-default
       action, or user semantic action that manipulated yychar.
     - Of course, the expected token list depends on states to have
       correct lookahead information, and it depends on the parser not
       to perform extra reductions after fetching a lookahead from the
       scanner and before detecting a syntax error.  Thus, state merging
       (from LALR or IELR) and default reductions corrupt the expected
       token list.  However, the list is correct for canonical LR with
       one exception: it will still contain any token that will not be
       accepted due to an error action in a later state.
  */
  if (yytoken != YYEMPTY)
    {
      int yyn = yypact[*yyssp];
      yyarg[yycount++] = yytname[yytoken];
      if (!yypact_value_is_default (yyn))
        {
          /* Start YYX at -YYN if negative to avoid negative indexes in
             YYCHECK.  In other words, skip the first -YYN actions for
             this state because they are default actions.  */
          int yyxbegin = yyn < 0 ? -yyn : 0;
          /* Stay within bounds of both yycheck and yytname.  */
          int yychecklim = YYLAST - yyn + 1;
          int yyxend = yychecklim < YYNTOKENS ? yychecklim : YYNTOKENS;
          int yyx;

          for (yyx = yyxbegin; yyx < yyxend; ++yyx)
            if (yycheck[yyx + yyn] == yyx && yyx != YYTERROR
                && !yytable_value_is_error (yytable[yyx + yyn]))
              {
                if (yycount == YYERROR_VERBOSE_ARGS_MAXIMUM)
                  {
                    yycount = 1;
                    yysize = yysize0;
                    break;
                  }
                yyarg[yycount++] = yytname[yyx];
                {
                  YYSIZE_T yysize1 = yysize + yytnamerr (YY_NULLPTR, yytname[yyx]);
                  if (! (yysize <= yysize1
                         && yysize1 <= YYSTACK_ALLOC_MAXIMUM))
                    return 2;
                  yysize = yysize1;
                }
              }
        }
    }

  switch (yycount)
    {
# define YYCASE_(N, S)                      \
      case N:                               \
        yyformat = S;                       \
      break
      YYCASE_(0, YY_("syntax error"));
      YYCASE_(1, YY_("syntax error, unexpected %s"));
      YYCASE_(2, YY_("syntax error, unexpected %s, expecting %s"));
      YYCASE_(3, YY_("syntax error, unexpected %s, expecting %s or %s"));
      YYCASE_(4, YY_("syntax error, unexpected %s, expecting %s or %s or %s"));
      YYCASE_(5, YY_("syntax error, unexpected %s, expecting %s or %s or %s or %s"));
# undef YYCASE_
    }

  {
    YYSIZE_T yysize1 = yysize + yystrlen (yyformat);
    if (! (yysize <= yysize1 && yysize1 <= YYSTACK_ALLOC_MAXIMUM))
      return 2;
    yysize = yysize1;
  }

  if (*yymsg_alloc < yysize)
    {
      *yymsg_alloc = 2 * yysize;
      if (! (yysize <= *yymsg_alloc
             && *yymsg_alloc <= YYSTACK_ALLOC_MAXIMUM))
        *yymsg_alloc = YYSTACK_ALLOC_MAXIMUM;
      return 1;
    }

  /* Avoid sprintf, as that infringes on the user's name space.
     Don't have undefined behavior even if the translation
     produced a string with the wrong number of "%s"s.  */
  {
    char *yyp = *yymsg;
    int yyi = 0;
    while ((*yyp = *yyformat) != '\0')
      if (*yyp == '%' && yyformat[1] == 's' && yyi < yycount)
        {
          yyp += yytnamerr (yyp, yyarg[yyi++]);
          yyformat += 2;
        }
      else
        {
          yyp++;
          yyformat++;
        }
  }
  return 0;
}
#endif /* YYERROR_VERBOSE */

/*-----------------------------------------------.
| Release the memory associated to this symbol.  |
`-----------------------------------------------*/

static void
yydestruct (const char *yymsg, int yytype, YYSTYPE *yyvaluep, struct mdlparse_vars *parse_state, yyscan_t scanner)
{
  YYUSE (yyvaluep);
  YYUSE (parse_state);
  YYUSE (scanner);
  if (!yymsg)
    yymsg = "Deleting";
  YY_SYMBOL_PRINT (yymsg, yytype, yyvaluep, yylocationp);

  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  YYUSE (yytype);
  YY_IGNORE_MAYBE_UNINITIALIZED_END
}




/*----------.
| yyparse.  |
`----------*/

int
yyparse (struct mdlparse_vars *parse_state, yyscan_t scanner)
{
/* The lookahead symbol.  */
int yychar;


/* The semantic value of the lookahead symbol.  */
/* Default value used for initialization, for pacifying older GCCs
   or non-GCC compilers.  */
YY_INITIAL_VALUE (static YYSTYPE yyval_default;)
YYSTYPE yylval YY_INITIAL_VALUE (= yyval_default);

    /* Number of syntax errors so far.  */
    int yynerrs;

    int yystate;
    /* Number of tokens to shift before error messages enabled.  */
    int yyerrstatus;

    /* The stacks and their tools:
       'yyss': related to states.
       'yyvs': related to semantic values.

       Refer to the stacks through separate pointers, to allow yyoverflow
       to reallocate them elsewhere.  */

    /* The state stack.  */
    yytype_int16 yyssa[YYINITDEPTH];
    yytype_int16 *yyss;
    yytype_int16 *yyssp;

    /* The semantic value stack.  */
    YYSTYPE yyvsa[YYINITDEPTH];
    YYSTYPE *yyvs;
    YYSTYPE *yyvsp;

    YYSIZE_T yystacksize;

  int yyn;
  int yyresult;
  /* Lookahead token as an internal (translated) token number.  */
  int yytoken = 0;
  /* The variables used to return semantic value and location from the
     action routines.  */
  YYSTYPE yyval;

#if YYERROR_VERBOSE
  /* Buffer for error messages, and its allocated size.  */
  char yymsgbuf[128];
  char *yymsg = yymsgbuf;
  YYSIZE_T yymsg_alloc = sizeof yymsgbuf;
#endif

#define YYPOPSTACK(N)   (yyvsp -= (N), yyssp -= (N))

  /* The number of symbols on the RHS of the reduced rule.
     Keep to zero when no symbol should be popped.  */
  int yylen = 0;

  yyssp = yyss = yyssa;
  yyvsp = yyvs = yyvsa;
  yystacksize = YYINITDEPTH;

  YYDPRINTF ((stderr, "Starting parse\n"));

  yystate = 0;
  yyerrstatus = 0;
  yynerrs = 0;
  yychar = YYEMPTY; /* Cause a token to be read.  */
  goto yysetstate;

/*------------------------------------------------------------.
| yynewstate -- Push a new state, which is found in yystate.  |
`------------------------------------------------------------*/
 yynewstate:
  /* In all cases, when you get here, the value and location stacks
     have just been pushed.  So pushing a state here evens the stacks.  */
  yyssp++;

 yysetstate:
  *yyssp = yystate;

  if (yyss + yystacksize - 1 <= yyssp)
    {
      /* Get the current used size of the three stacks, in elements.  */
      YYSIZE_T yysize = yyssp - yyss + 1;

#ifdef yyoverflow
      {
        /* Give user a chance to reallocate the stack.  Use copies of
           these so that the &'s don't force the real ones into
           memory.  */
        YYSTYPE *yyvs1 = yyvs;
        yytype_int16 *yyss1 = yyss;

        /* Each stack pointer address is followed by the size of the
           data in use in that stack, in bytes.  This used to be a
           conditional around just the two extra args, but that might
           be undefined if yyoverflow is a macro.  */
        yyoverflow (YY_("memory exhausted"),
                    &yyss1, yysize * sizeof (*yyssp),
                    &yyvs1, yysize * sizeof (*yyvsp),
                    &yystacksize);

        yyss = yyss1;
        yyvs = yyvs1;
      }
#else /* no yyoverflow */
# ifndef YYSTACK_RELOCATE
      goto yyexhaustedlab;
# else
      /* Extend the stack our own way.  */
      if (YYMAXDEPTH <= yystacksize)
        goto yyexhaustedlab;
      yystacksize *= 2;
      if (YYMAXDEPTH < yystacksize)
        yystacksize = YYMAXDEPTH;

      {
        yytype_int16 *yyss1 = yyss;
        union yyalloc *yyptr =
          (union yyalloc *) YYSTACK_ALLOC (YYSTACK_BYTES (yystacksize));
        if (! yyptr)
          goto yyexhaustedlab;
        YYSTACK_RELOCATE (yyss_alloc, yyss);
        YYSTACK_RELOCATE (yyvs_alloc, yyvs);
#  undef YYSTACK_RELOCATE
        if (yyss1 != yyssa)
          YYSTACK_FREE (yyss1);
      }
# endif
#endif /* no yyoverflow */

      yyssp = yyss + yysize - 1;
      yyvsp = yyvs + yysize - 1;

      YYDPRINTF ((stderr, "Stack size increased to %lu\n",
                  (unsigned long int) yystacksize));

      if (yyss + yystacksize - 1 <= yyssp)
        YYABORT;
    }

  YYDPRINTF ((stderr, "Entering state %d\n", yystate));

  if (yystate == YYFINAL)
    YYACCEPT;

  goto yybackup;

/*-----------.
| yybackup.  |
`-----------*/
yybackup:

  /* Do appropriate processing given the current state.  Read a
     lookahead token if we need one and don't already have one.  */

  /* First try to decide what to do without reference to lookahead token.  */
  yyn = yypact[yystate];
  if (yypact_value_is_default (yyn))
    goto yydefault;

  /* Not known => get a lookahead token if don't already have one.  */

  /* YYCHAR is either YYEMPTY or YYEOF or a valid lookahead symbol.  */
  if (yychar == YYEMPTY)
    {
      YYDPRINTF ((stderr, "Reading a token: "));
      yychar = yylex (&yylval, parse_state, scanner);
    }

  if (yychar <= YYEOF)
    {
      yychar = yytoken = YYEOF;
      YYDPRINTF ((stderr, "Now at end of input.\n"));
    }
  else
    {
      yytoken = YYTRANSLATE (yychar);
      YY_SYMBOL_PRINT ("Next token is", yytoken, &yylval, &yylloc);
    }

  /* If the proper action on seeing token YYTOKEN is to reduce or to
     detect an error, take that action.  */
  yyn += yytoken;
  if (yyn < 0 || YYLAST < yyn || yycheck[yyn] != yytoken)
    goto yydefault;
  yyn = yytable[yyn];
  if (yyn <= 0)
    {
      if (yytable_value_is_error (yyn))
        goto yyerrlab;
      yyn = -yyn;
      goto yyreduce;
    }

  /* Count tokens shifted since error; after three, turn off error
     status.  */
  if (yyerrstatus)
    yyerrstatus--;

  /* Shift the lookahead token.  */
  YY_SYMBOL_PRINT ("Shifting", yytoken, &yylval, &yylloc);

  /* Discard the shifted token.  */
  yychar = YYEMPTY;

  yystate = yyn;
  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  *++yyvsp = yylval;
  YY_IGNORE_MAYBE_UNINITIALIZED_END

  goto yynewstate;


/*-----------------------------------------------------------.
| yydefault -- do the default action for the current state.  |
`-----------------------------------------------------------*/
yydefault:
  yyn = yydefact[yystate];
  if (yyn == 0)
    goto yyerrlab;
  goto yyreduce;


/*-----------------------------.
| yyreduce -- Do a reduction.  |
`-----------------------------*/
yyreduce:
  /* yyn is the number of a rule to reduce with.  */
  yylen = yyr2[yyn];

  /* If YYLEN is nonzero, implement the default value of the action:
     '$$ = $1'.

     Otherwise, the following line sets YYVAL to garbage.
     This behavior is undocumented and Bison
     users should not rely upon it.  Assigning to YYVAL
     unconditionally makes the parser a bit smaller, and it avoids a
     GCC warning that YYVAL may be used uninitialized.  */
  yyval = yyvsp[1-yylen];


  YY_REDUCE_PRINT (yyn);
  switch (yyn)
    {
        case 28:
#line 693 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.sym) = mdl_existing_object(parse_state, (yyvsp[0].str))); }
#line 3362 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 29:
#line 696 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.symlist) = mdl_singleton_symbol_list(parse_state, (yyvsp[0].sym))); }
#line 3368 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 30:
#line 697 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.symlist) = mdl_existing_objects_wildcard(parse_state, (yyvsp[0].str))); }
#line 3374 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 31:
#line 700 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.sym) = mdl_existing_region(parse_state, (yyvsp[-3].sym), (yyvsp[-1].str))); }
#line 3380 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 32:
#line 703 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.vec3) = mdl_point(parse_state, &(yyvsp[0].nlist))); }
#line 3386 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 34:
#line 707 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.vec3) = mdl_point_scalar((yyvsp[0].dbl))); }
#line 3392 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 35:
#line 710 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = 1; }
#line 3398 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 36:
#line 711 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = 0; }
#line 3404 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 37:
#line 712 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = 1; }
#line 3410 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 38:
#line 713 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = 0; }
#line 3416 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 39:
#line 714 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = 1; }
#line 3422 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 40:
#line 715 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = 0; }
#line 3428 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 41:
#line 718 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.mol_type).orient_set = 0; }
#line 3434 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 44:
#line 721 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.mol_type).orient_set = 1; (yyval.mol_type).orient = 0; }
#line 3440 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 45:
#line 725 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.mol_type).orient = 1; (yyval.mol_type).orient_set = 1; }
#line 3446 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 46:
#line 726 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.mol_type).orient = -1; (yyval.mol_type).orient_set = 1; }
#line 3452 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 47:
#line 727 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                          (yyval.mol_type) = (yyvsp[-1].mol_type);
                                                          if ((yyval.mol_type).orient >= 32767)
                                                          {
                                                            /* Seriously?  Wow. */
                                                            mdlerror(parse_state, "Error: Molecule orientation must not be greater than 32767");
                                                            return 1;
                                                          }
                                                          ++ (yyval.mol_type).orient;
                                                      }
#line 3467 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 48:
#line 737 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                          (yyval.mol_type) = (yyvsp[-1].mol_type);
                                                          if ((yyval.mol_type).orient <= -32768)
                                                          {
                                                            /* Seriously?  Wow. */
                                                            mdlerror(parse_state, "Error: Molecule orientation must not be less than -32768");
                                                            return 1;
                                                          }
                                                          -- (yyval.mol_type).orient;
                                                      }
#line 3482 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 51:
#line 755 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                          (yyval.mol_type).orient = (int) (yyvsp[-1].dbl);
                                                          (yyval.mol_type).orient_set = 1;
                                                          if ((yyval.mol_type).orient != (yyvsp[-1].dbl))
                                                          {
                                                            mdlerror(parse_state, "Molecule orientation specified inside braces must be an integer between -32768 and 32767.");
                                                            return 1;
                                                          }
                                                      }
#line 3496 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 53:
#line 768 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                          if ((yyvsp[-2].nlist).value_tail)
                                                          {
                                                            (yyval.nlist) = (yyvsp[-2].nlist);
                                                            (yyval.nlist).value_count += (yyvsp[0].nlist).value_count;
                                                            (yyval.nlist).value_tail->next = (yyvsp[0].nlist).value_head;
                                                            (yyval.nlist).value_tail = (yyvsp[0].nlist).value_tail;
                                                          }
                                                          else
                                                            (yyval.nlist) = (yyvsp[0].nlist);
                                                      }
#line 3512 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 54:
#line 781 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mcell_generate_range_singleton(&(yyval.nlist), (yyvsp[0].dbl))); }
#line 3518 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 55:
#line 782 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_generate_range(parse_state, &(yyval.nlist), (yyvsp[-5].dbl), (yyvsp[-3].dbl), (yyvsp[-1].dbl))); }
#line 3524 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 56:
#line 788 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                          char *include_path = mdl_find_include_file((yyvsp[0].str), parse_state->vol->curr_file);
                                                          if (include_path == NULL)
                                                          {
                                                            mdlerror_fmt(parse_state, "Out of memory while trying to open include file '%s'", (yyvsp[0].str));
                                                            free((yyvsp[0].str));
                                                            return 1;
                                                          }
                                                          if (mdlparse_file(parse_state, include_path))
                                                          {
                                                            free(include_path);
                                                            free((yyvsp[0].str));
                                                            return 1;
                                                          }
                                                          free(include_path);
                                                          free((yyvsp[0].str));
                                                      }
#line 3546 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 57:
#line 811 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_assign_variable_double(parse_state, (yyvsp[-2].sym), (yyvsp[0].dbl))); }
#line 3552 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 58:
#line 812 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_assign_variable_string(parse_state, (yyvsp[-2].sym), (yyvsp[0].str))); }
#line 3558 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 59:
#line 813 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_assign_variable(parse_state, (yyvsp[-2].sym), (yyvsp[0].sym))); }
#line 3564 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 60:
#line 814 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_assign_variable_array(parse_state, (yyvsp[-2].sym), (yyvsp[0].nlist).value_head)); }
#line 3570 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 61:
#line 817 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.sym) = mdl_get_or_create_variable(parse_state, (yyvsp[0].str))); }
#line 3576 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 62:
#line 820 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.sym) = mdl_existing_variable(parse_state, (yyvsp[0].str))); }
#line 3582 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 64:
#line 824 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                          struct num_expr_list *elp;
                                                          (yyval.nlist).value_head = (struct num_expr_list *) (yyvsp[0].sym)->value;
                                                          (yyval.nlist).value_count = 1;
                                                          for (elp = (yyval.nlist).value_head; elp->next != NULL; elp = elp->next)
                                                            ++ (yyval.nlist).value_count;
                                                          (yyval.nlist).value_tail = elp;
                                                          (yyval.nlist).shared = 1;
                                                      }
#line 3596 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 65:
#line 835 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { mdl_debug_dump_array((yyvsp[-1].nlist).value_head); (yyval.nlist) = (yyvsp[-1].nlist); }
#line 3602 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 66:
#line 838 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.sym) = mdl_existing_array(parse_state, (yyvsp[0].str))); }
#line 3608 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 70:
#line 846 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.dbl) = *(double *) (yyvsp[0].sym)->value; }
#line 3614 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 71:
#line 849 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.dbl) = (yyvsp[0].llival); }
#line 3620 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 75:
#line 857 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.sym) = mdl_existing_double(parse_state, (yyvsp[0].str))); }
#line 3626 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 76:
#line 861 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.dbl) = (yyvsp[-1].dbl); }
#line 3632 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 77:
#line 862 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKF((yyval.dbl) = exp((yyvsp[-1].dbl))); }
#line 3638 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 78:
#line 863 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_expr_log(parse_state, (yyvsp[-1].dbl), &(yyval.dbl))); }
#line 3644 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 79:
#line 864 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_expr_log10(parse_state, (yyvsp[-1].dbl), &(yyval.dbl))); }
#line 3650 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 80:
#line 865 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.dbl) = max2d((yyvsp[-3].dbl), (yyvsp[-1].dbl)); }
#line 3656 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 81:
#line 866 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.dbl) = min2d((yyvsp[-3].dbl), (yyvsp[-1].dbl)); }
#line 3662 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 82:
#line 867 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.dbl) = mdl_expr_roundoff((yyvsp[-1].dbl), (int) (yyvsp[-3].dbl)); }
#line 3668 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 83:
#line 868 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.dbl) = floor((yyvsp[-1].dbl)); }
#line 3674 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 84:
#line 869 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.dbl) = ceil((yyvsp[-1].dbl)); }
#line 3680 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 85:
#line 870 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.dbl) = sin((yyvsp[-1].dbl)); }
#line 3686 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 86:
#line 871 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.dbl) = cos((yyvsp[-1].dbl)); }
#line 3692 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 87:
#line 872 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKF((yyval.dbl) = tan((yyvsp[-1].dbl))); }
#line 3698 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 88:
#line 873 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKF((yyval.dbl) = asin((yyvsp[-1].dbl))); }
#line 3704 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 89:
#line 874 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKF((yyval.dbl) = acos((yyvsp[-1].dbl))); }
#line 3710 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 90:
#line 875 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.dbl) = atan((yyvsp[-1].dbl)); }
#line 3716 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 91:
#line 876 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKF((yyval.dbl) = sqrt((yyvsp[-1].dbl))); }
#line 3722 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 92:
#line 877 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.dbl) = fabs((yyvsp[-1].dbl)); }
#line 3728 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 93:
#line 878 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_expr_mod(parse_state, (yyvsp[-3].dbl), (yyvsp[-1].dbl), &(yyval.dbl))); }
#line 3734 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 94:
#line 879 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.dbl) = MY_PI; }
#line 3740 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 95:
#line 880 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.dbl) = mdl_expr_rng_uniform(parse_state); }
#line 3746 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 96:
#line 881 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.dbl) = rng_gauss(parse_state->vol->rng); }
#line 3752 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 97:
#line 882 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.dbl) = parse_state->vol->seed_seq; }
#line 3758 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 98:
#line 883 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_expr_string_to_double(parse_state, (yyvsp[-1].str), &(yyval.dbl))); }
#line 3764 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 99:
#line 884 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKF((yyval.dbl) = (yyvsp[-2].dbl) + (yyvsp[0].dbl)); }
#line 3770 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 100:
#line 885 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKF((yyval.dbl) = (yyvsp[-2].dbl) - (yyvsp[0].dbl)); }
#line 3776 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 101:
#line 886 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKF((yyval.dbl) = (yyvsp[-2].dbl) * (yyvsp[0].dbl)); }
#line 3782 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 102:
#line 887 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_expr_div(parse_state, (yyvsp[-2].dbl), (yyvsp[0].dbl), &(yyval.dbl))); }
#line 3788 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 103:
#line 888 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_expr_pow(parse_state, (yyvsp[-2].dbl), (yyvsp[0].dbl), &(yyval.dbl))); }
#line 3794 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 104:
#line 889 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.dbl) = -(yyvsp[0].dbl); }
#line 3800 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 105:
#line 890 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.dbl) = (yyvsp[0].dbl); }
#line 3806 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 107:
#line 895 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.str) = mdl_strdup((char const *) (yyvsp[0].sym)->value)); }
#line 3812 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 108:
#line 899 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.str) = mdl_strip_quotes((yyvsp[0].str))); }
#line 3818 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 109:
#line 900 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.str) = mdl_strdup(parse_state->vol->mdl_infile_name)); }
#line 3824 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 110:
#line 901 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.str) = mdl_strcat((yyvsp[-2].str), (yyvsp[0].str))); }
#line 3830 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 111:
#line 902 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.str) = mdl_string_format(parse_state, (yyvsp[-2].str), (yyvsp[-1].printfargs).arg_head)); }
#line 3836 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 112:
#line 905 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.sym) = mdl_existing_string(parse_state, (yyvsp[0].str))); }
#line 3842 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 120:
#line 921 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_fopen(parse_state, (yyvsp[-6].sym), (yyvsp[-3].str), (yyvsp[-1].str))); }
#line 3848 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 121:
#line 924 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.sym) = mdl_new_filehandle(parse_state, (yyvsp[0].str))); }
#line 3854 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 122:
#line 927 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.str) = (yyvsp[0].str); CHECK(mdl_valid_file_mode(parse_state, (yyvsp[0].str))); }
#line 3860 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 123:
#line 930 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_fclose(parse_state, (yyvsp[-1].sym))); }
#line 3866 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 124:
#line 933 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.sym) = mdl_existing_file_stream(parse_state, (yyvsp[0].str))); }
#line 3872 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 125:
#line 936 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.str) = mdl_expand_string_escapes((yyvsp[0].str))); }
#line 3878 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 126:
#line 939 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.printfargs).arg_head = (yyval.printfargs).arg_tail = NULL; }
#line 3884 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 127:
#line 940 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                        (yyval.printfargs) = (yyvsp[-2].printfargs);
                                                        if ((yyval.printfargs).arg_tail)
                                                          (yyval.printfargs).arg_tail = (yyval.printfargs).arg_tail->next = (yyvsp[0].printfarg);
                                                        else
                                                          (yyval.printfargs).arg_tail = (yyval.printfargs).arg_head = (yyvsp[0].printfarg);
                                                        (yyvsp[0].printfarg)->next = NULL;
                                                      }
#line 3897 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 128:
#line 950 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.printfarg) = mdl_new_printf_arg_double((yyvsp[0].dbl))); }
#line 3903 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 129:
#line 951 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.printfarg) = mdl_new_printf_arg_string((yyvsp[0].str))); }
#line 3909 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 130:
#line 952 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                          switch ((yyvsp[0].sym)->sym_type)
                                                          {
                                                            case DBL: CHECKN((yyval.printfarg) = mdl_new_printf_arg_double(*(double *) (yyvsp[0].sym)->value)); break;
                                                            case STR: CHECKN((yyval.printfarg) = mdl_new_printf_arg_string((char *) (yyvsp[0].sym)->value)); break;
                                                            default:
                                                              mdlerror(parse_state, "Invalid variable type referenced");
                                                              return 1;
                                                          }
                                                      }
#line 3924 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 131:
#line 964 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_printf(parse_state, (yyvsp[-2].str), (yyvsp[-1].printfargs).arg_head)); }
#line 3930 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 132:
#line 970 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_fprintf(parse_state, (struct file_stream *) (yyvsp[-4].sym)->value, (yyvsp[-2].str), (yyvsp[-1].printfargs).arg_head)); }
#line 3936 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 133:
#line 976 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_sprintf(parse_state, (yyvsp[-4].sym), (yyvsp[-2].str), (yyvsp[-1].printfargs).arg_head)); }
#line 3942 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 134:
#line 979 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { mdl_print_time(parse_state, (yyvsp[-1].str)); }
#line 3948 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 135:
#line 985 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_fprint_time(parse_state, (yyvsp[-3].sym), (yyvsp[-1].str))); }
#line 3954 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 139:
#line 1001 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { if (!parse_state->vol->quiet_flag) mdl_set_all_notifications(parse_state->vol, (yyvsp[0].tok)); }
#line 3960 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 140:
#line 1002 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { if (!parse_state->vol->quiet_flag) parse_state->vol->notify->progress_report        = (yyvsp[0].tok); }
#line 3966 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 141:
#line 1003 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { if (!parse_state->vol->quiet_flag) parse_state->vol->notify->diffusion_constants    = (yyvsp[0].tok); }
#line 3972 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 142:
#line 1004 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { if (!parse_state->vol->quiet_flag) parse_state->vol->notify->reaction_probabilities = (yyvsp[0].tok); }
#line 3978 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 143:
#line 1005 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { if (!parse_state->vol->quiet_flag) parse_state->vol->notify->time_varying_reactions = (yyvsp[0].tok); }
#line 3984 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 144:
#line 1006 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { if (!parse_state->vol->quiet_flag) parse_state->vol->notify->reaction_prob_notify   = (yyvsp[0].dbl); }
#line 3990 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 145:
#line 1007 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { if (!parse_state->vol->quiet_flag) parse_state->vol->notify->partition_location     = (yyvsp[0].tok); }
#line 3996 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 146:
#line 1008 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { if (!parse_state->vol->quiet_flag) parse_state->vol->notify->box_triangulation      = (yyvsp[0].tok); }
#line 4002 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 147:
#line 1009 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { if (!parse_state->vol->quiet_flag) parse_state->vol->notify->release_events         = (yyvsp[0].tok); }
#line 4008 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 148:
#line 1010 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { if (!parse_state->vol->quiet_flag) parse_state->vol->notify->file_writes            = (yyvsp[0].tok); }
#line 4014 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 149:
#line 1011 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { if (!parse_state->vol->quiet_flag) parse_state->vol->notify->final_summary          = (yyvsp[0].tok); }
#line 4020 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 150:
#line 1012 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { if (!parse_state->vol->quiet_flag) parse_state->vol->notify->throughput_report      = (yyvsp[0].tok); }
#line 4026 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 151:
#line 1013 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { if (!parse_state->vol->quiet_flag) parse_state->vol->notify->reaction_output_report = (yyvsp[0].tok); }
#line 4032 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 152:
#line 1014 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { if (!parse_state->vol->quiet_flag) parse_state->vol->notify->volume_output_report   = (yyvsp[0].tok); }
#line 4038 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 153:
#line 1015 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { if (!parse_state->vol->quiet_flag) parse_state->vol->notify->viz_output_report      = (yyvsp[0].tok); }
#line 4044 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 154:
#line 1016 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { if (!parse_state->vol->quiet_flag) parse_state->vol->notify->checkpoint_report      = (yyvsp[0].tok); }
#line 4050 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 155:
#line 1017 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                          if (!parse_state->vol->quiet_flag && parse_state->vol->log_freq == ULONG_MAX)
                                                            parse_state->vol->notify->iteration_report = (yyvsp[0].tok);
                                                      }
#line 4059 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 156:
#line 1021 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { if (!parse_state->vol->quiet_flag) CHECK(mdl_set_iteration_report_freq(parse_state, (long long) (yyvsp[0].dbl))); }
#line 4065 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 157:
#line 1022 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { if (!parse_state->vol->quiet_flag) parse_state->vol->notify->molecule_collision_report    = (yyvsp[0].tok); }
#line 4071 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 158:
#line 1026 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = ((yyvsp[0].tok) ? NOTIFY_FULL : NOTIFY_NONE); }
#line 4077 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 159:
#line 1030 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = ((yyvsp[0].tok) ? NOTIFY_FULL : NOTIFY_NONE); }
#line 4083 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 160:
#line 1031 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = NOTIFY_BRIEF; }
#line 4089 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 164:
#line 1047 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { mdl_set_all_warnings(parse_state->vol, (byte) (yyvsp[0].tok)); }
#line 4095 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 165:
#line 1048 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->vol->notify->neg_diffusion = (byte)(yyvsp[0].tok); }
#line 4101 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 166:
#line 1049 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->vol->notify->neg_reaction = (byte)(yyvsp[0].tok); }
#line 4107 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 167:
#line 1050 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->vol->notify->high_reaction_prob = (byte)(yyvsp[0].tok); }
#line 4113 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 168:
#line 1051 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->vol->notify->reaction_prob_warn = (yyvsp[0].dbl); }
#line 4119 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 169:
#line 1052 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->vol->notify->close_partitions = (byte)(yyvsp[0].tok); }
#line 4125 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 170:
#line 1053 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->vol->notify->degenerate_polys = (byte)(yyvsp[0].tok); }
#line 4131 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 171:
#line 1054 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->vol->notify->overwritten_file = (byte)(yyvsp[0].tok); }
#line 4137 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 172:
#line 1055 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->vol->notify->short_lifetime = (byte)(yyvsp[0].tok); }
#line 4143 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 173:
#line 1056 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_lifetime_warning_threshold(parse_state, (long long) (yyvsp[0].dbl))); }
#line 4149 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 174:
#line 1057 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->vol->notify->missed_reactions = (byte)(yyvsp[0].tok); }
#line 4155 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 175:
#line 1058 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_missed_reaction_warning_threshold(parse_state, (yyvsp[0].dbl))); }
#line 4161 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 176:
#line 1059 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->vol->notify->missed_surf_orient = (byte)(yyvsp[0].tok); }
#line 4167 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 177:
#line 1060 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->vol->notify->useless_vol_orient = (byte)(yyvsp[0].tok); }
#line 4173 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 178:
#line 1061 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->vol->notify->complex_placement_failure = (byte) (yyvsp[0].tok); }
#line 4179 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 179:
#line 1062 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->vol->notify->complex_placement_failure_threshold = (long long) (yyvsp[0].dbl); }
#line 4185 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 180:
#line 1063 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->vol->notify->mol_placement_failure = (byte) (yyvsp[0].tok); }
#line 4191 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 181:
#line 1064 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->vol->notify->invalid_output_step_time = (byte) (yyvsp[0].tok); }
#line 4197 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 182:
#line 1068 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = WARN_COPE;  }
#line 4203 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 183:
#line 1069 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = WARN_WARN;  }
#line 4209 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 184:
#line 1070 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = WARN_ERROR; }
#line 4215 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 185:
#line 1076 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_checkpoint_infile(parse_state, (yyvsp[0].str))); }
#line 4221 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 186:
#line 1077 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_checkpoint_outfile(parse_state, (yyvsp[0].str))); }
#line 4227 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 187:
#line 1078 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_checkpoint_interval(parse_state, (yyvsp[-1].dbl), (yyvsp[0].tok))); }
#line 4233 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 188:
#line 1079 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_keep_checkpoint_files(parse_state, (yyvsp[0].tok))); }
#line 4239 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 189:
#line 1081 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_realtime_checkpoint(parse_state, (long) (yyvsp[-1].dbl), (yyvsp[0].tok))); }
#line 4245 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 190:
#line 1084 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = 0; }
#line 4251 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 191:
#line 1085 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = 1; }
#line 4257 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 192:
#line 1086 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = 0; }
#line 4263 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 193:
#line 1090 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { /* seconds */     (yyval.dbl) = (yyvsp[0].dbl); }
#line 4269 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 194:
#line 1091 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { /* mm:ss */       (yyval.dbl) = (yyvsp[-2].dbl) * 60 + (yyvsp[0].dbl); }
#line 4275 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 195:
#line 1092 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { /* hh:mm:ss */    (yyval.dbl) = (yyvsp[-4].dbl) * 3600 + (yyvsp[-2].dbl) * 60 + (yyvsp[0].dbl); }
#line 4281 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 196:
#line 1094 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { /* dd:hh:mm:ss */ (yyval.dbl) = (yyvsp[-6].dbl) * 86400 + (yyvsp[-4].dbl) * 3600 + (yyvsp[-2].dbl) * 60 + (yyvsp[0].dbl); }
#line 4287 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 197:
#line 1101 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_time_step(parse_state, (yyvsp[0].dbl))); }
#line 4293 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 198:
#line 1102 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_space_step(parse_state, (yyvsp[0].dbl))); }
#line 4299 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 199:
#line 1103 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_max_time_step(parse_state, (yyvsp[0].dbl))); }
#line 4305 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 200:
#line 1104 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_num_iterations(parse_state, (long long) (yyvsp[0].dbl))); }
#line 4311 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 201:
#line 1105 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->vol->randomize_smol_pos = !((yyvsp[0].tok)); }
#line 4317 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 202:
#line 1106 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->vol->use_expanded_list = (yyvsp[0].tok); }
#line 4323 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 203:
#line 1107 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->vol->vacancy_search_dist2 = max2d((yyvsp[0].dbl), 0.0); }
#line 4329 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 204:
#line 1108 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_num_radial_directions(parse_state, (int) (yyvsp[0].dbl))); }
#line 4335 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 205:
#line 1109 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->vol->fully_random = 1; }
#line 4341 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 206:
#line 1110 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_num_radial_subdivisions(parse_state, (int) (yyvsp[0].dbl))); }
#line 4347 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 207:
#line 1111 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_grid_density(parse_state, (yyvsp[0].dbl))); }
#line 4353 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 208:
#line 1112 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_interaction_radius(parse_state, (yyvsp[0].dbl))); }
#line 4359 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 209:
#line 1113 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->vol->surface_reversibility=(yyvsp[0].tok); parse_state->vol->volume_reversibility=(yyvsp[0].tok); }
#line 4365 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 210:
#line 1114 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->vol->surface_reversibility=1;  parse_state->vol->volume_reversibility=0;  }
#line 4371 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 211:
#line 1115 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->vol->surface_reversibility=0;  parse_state->vol->volume_reversibility=1;  }
#line 4377 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 212:
#line 1116 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_complex_placement_attempts(parse_state, (yyvsp[0].dbl))); }
#line 4383 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 213:
#line 1123 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->vol->mem_part_x = (int) (yyvsp[0].dbl); }
#line 4389 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 214:
#line 1124 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->vol->mem_part_y = (int) (yyvsp[0].dbl); }
#line 4395 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 215:
#line 1125 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->vol->mem_part_z = (int) (yyvsp[0].dbl); }
#line 4401 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 216:
#line 1126 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->vol->mem_part_pool = (int) (yyvsp[0].dbl); }
#line 4407 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 217:
#line 1130 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mcell_set_partition(parse_state->vol, (yyvsp[-2].tok), & (yyvsp[0].nlist))); }
#line 4413 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 218:
#line 1134 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = X_PARTS; }
#line 4419 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 219:
#line 1135 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = Y_PARTS; }
#line 4425 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 220:
#line 1136 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = Z_PARTS; }
#line 4431 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 224:
#line 1148 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { mdl_print_species_summary(parse_state->vol, (yyvsp[0].mcell_mol_spec)); }
#line 4437 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 225:
#line 1152 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { mdl_print_species_summaries(parse_state->vol, (yyvsp[-1].mcell_species_lst).species_head); }
#line 4443 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 226:
#line 1156 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.mcell_species_lst).species_count = 0; CHECK(mdl_add_to_species_list(&(yyval.mcell_species_lst), (yyvsp[0].mcell_mol_spec))); }
#line 4449 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 227:
#line 1157 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.mcell_species_lst) = (yyvsp[-1].mcell_species_lst); CHECK(mdl_add_to_species_list(&(yyval.mcell_species_lst), (yyvsp[0].mcell_mol_spec))); }
#line 4455 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 228:
#line 1166 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.mcell_mol_spec) = mdl_create_species(parse_state, (yyvsp[-6].str), (yyvsp[-4].diff_const).D, (yyvsp[-4].diff_const).is_2d, (yyvsp[-3].dbl), (yyvsp[-1].dbl), (yyvsp[-1].dbl) )); }
#line 4461 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 230:
#line 1172 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.sym) = mdl_new_mol_species(parse_state, (yyvsp[0].str))); }
#line 4467 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 231:
#line 1176 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.diff_const).is_2d = 0; (yyval.diff_const).D = (yyvsp[0].dbl); CHECK(mdl_check_diffusion_constant(parse_state, & (yyval.diff_const).D)); }
#line 4473 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 232:
#line 1177 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.diff_const).is_2d = 1; (yyval.diff_const).D = (yyvsp[0].dbl); CHECK(mdl_check_diffusion_constant(parse_state, & (yyval.diff_const).D)); }
#line 4479 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 233:
#line 1181 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.dbl) = 0.0; }
#line 4485 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 234:
#line 1182 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                          if ((yyvsp[0].dbl) <= 0)
                                                          {
                                                            mdlerror_fmt(parse_state, "Requested custom time step of %.15g; custom time step must be positive.", (yyvsp[0].dbl));
                                                            return 1;
                                                          }

                                                          (yyval.dbl) = (yyvsp[0].dbl);
                                                      }
#line 4499 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 235:
#line 1191 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                          if ((yyvsp[0].dbl) <= 0)
                                                          {
                                                            mdlerror_fmt(parse_state, "Requested custom space step of %.15g; custom space step must be positive.", (yyvsp[0].dbl));
                                                            return 1;
                                                          }

                                                          (yyval.dbl) = -(yyvsp[0].dbl);
                                                      }
#line 4513 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 236:
#line 1202 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.ival) = 0; }
#line 4519 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 237:
#line 1203 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.ival) = 1; }
#line 4525 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 238:
#line 1207 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.dbl) = 0; }
#line 4531 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 239:
#line 1208 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                        if ((yyvsp[0].dbl) <= 0)
                                                        {
                                                          mdlerror_fmt(parse_state, "Requested maximum step length of %.15g; maximum step length must be positive.", (yyvsp[0].dbl));
                                                          return 1;
                                                        }
                                                        (yyval.dbl) = (yyvsp[0].dbl);
                                                      }
#line 4544 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 240:
#line 1219 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->complex_name = (yyvsp[0].str); parse_state->complex_type = 0; }
#line 4550 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 241:
#line 1221 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->complex_topo = (yyvsp[0].mmol_topo); }
#line 4556 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 242:
#line 1224 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->complex_relations = (yyvsp[0].mmol_su_rel); }
#line 4562 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 243:
#line 1226 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_assemble_complex_species(parse_state, (yyvsp[-10].str), (yyvsp[-7].mmol_topo), (yyvsp[-5].mmol_subunits).assign_head, (yyvsp[-4].mmol_geom), (yyvsp[-3].mmol_su_rel), (yyvsp[-1].mmol_rate_ruleset))); }
#line 4568 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 244:
#line 1229 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.str) = (yyvsp[0].str); CHECK(mdl_valid_complex_name(parse_state, (yyvsp[0].str))); }
#line 4574 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 245:
#line 1233 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.mmol_topo) = mdl_assemble_topology(parse_state, &(yyvsp[0].nlist))); }
#line 4580 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 246:
#line 1237 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.mmol_subunits).assign_tail = (yyval.mmol_subunits).assign_head = (yyvsp[0].mmol_su_assign); }
#line 4586 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 247:
#line 1239 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.mmol_subunits) = (yyvsp[-1].mmol_subunits); (yyval.mmol_subunits).assign_tail = (yyval.mmol_subunits).assign_tail->next = (yyvsp[0].mmol_su_assign); }
#line 4592 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 248:
#line 1244 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.mmol_su_assign) = mdl_assemble_complex_subunit_assignment(parse_state, (yyvsp[-3].mmol_su_comp), & (yyvsp[0].mol_type))); }
#line 4598 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 250:
#line 1250 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { if ((yyvsp[0].mmol_su_comp)) (yyvsp[0].mmol_su_comp)->next = (yyvsp[-2].mmol_su_comp); (yyval.mmol_su_comp) = (yyvsp[0].mmol_su_comp); }
#line 4604 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 251:
#line 1253 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.mmol_su_comp) = mdl_assemble_subunit_spec_component((yyvsp[0].dbl), (yyvsp[0].dbl))); }
#line 4610 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 252:
#line 1254 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.mmol_su_comp) = mdl_assemble_subunit_spec_component((yyvsp[-2].dbl), (yyvsp[0].dbl))); }
#line 4616 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 253:
#line 1258 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.mmol_geom) = (yyvsp[-1].mmol_geom); CHECK(mdl_validate_complex_geometry(parse_state, parse_state->complex_topo, (yyvsp[-1].mmol_geom))); }
#line 4622 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 255:
#line 1264 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { if ((yyvsp[0].mmol_geom)) (yyvsp[0].mmol_geom)->next = (yyvsp[-1].mmol_geom); (yyval.mmol_geom) = (yyvsp[0].mmol_geom); }
#line 4628 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 256:
#line 1268 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mcell_generate_range_singleton(&(yyval.nlist), (yyvsp[0].dbl))); }
#line 4634 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 257:
#line 1269 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.nlist) = (yyvsp[-2].nlist); CHECK(mdl_add_range_value(&(yyval.nlist), (yyvsp[0].dbl))); }
#line 4640 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 258:
#line 1273 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.mmol_geom) = mdl_assemble_complex_geometry(parse_state, parse_state->complex_topo, &(yyvsp[-3].nlist), (yyvsp[0].vec3))); }
#line 4646 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 259:
#line 1278 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.mmol_su_rel) = (yyvsp[-1].mmol_su_rel); CHECK(mdl_validate_complex_relationships(parse_state, parse_state->complex_topo, (yyvsp[-1].mmol_su_rel))); }
#line 4652 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 260:
#line 1282 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.mmol_su_rel) = NULL; }
#line 4658 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 261:
#line 1284 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {  if ((yyvsp[0].mmol_su_rel)) (yyvsp[0].mmol_su_rel)->next = (yyvsp[-1].mmol_su_rel); (yyval.mmol_su_rel) = (yyvsp[0].mmol_su_rel); }
#line 4664 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 262:
#line 1287 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.mmol_su_rel) = mdl_assemble_complex_relationship(parse_state, parse_state->complex_topo, (yyvsp[-2].str), &(yyvsp[0].nlist))); }
#line 4670 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 263:
#line 1291 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.mmol_rate_ruleset) = (yyvsp[-1].mmol_rate_ruleset); CHECK(mdl_validate_complex_rates(parse_state, (yyvsp[-1].mmol_rate_ruleset))); }
#line 4676 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 264:
#line 1295 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.mmol_rate_ruleset) = NULL; }
#line 4682 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 265:
#line 1296 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { if ((yyvsp[0].mmol_rate_ruleset)) (yyvsp[0].mmol_rate_ruleset)->next = (yyvsp[-1].mmol_rate_ruleset); (yyval.mmol_rate_ruleset) = (yyvsp[0].mmol_rate_ruleset); }
#line 4688 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 266:
#line 1299 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.mmol_rate_ruleset) = mdl_assemble_complex_ruleset((yyvsp[-3].str), (yyvsp[-1].mmol_rate_rule))); }
#line 4694 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 268:
#line 1305 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { if ((yyvsp[0].mmol_rate_rule)) (yyvsp[0].mmol_rate_rule)->next = (yyvsp[-1].mmol_rate_rule); (yyval.mmol_rate_rule) = (yyvsp[0].mmol_rate_rule); }
#line 4700 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 269:
#line 1309 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.mmol_rate_rule) = mdl_assemble_complex_rate_rule((yyvsp[-2].mmol_rate_clause), (yyvsp[0].dbl))); }
#line 4706 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 271:
#line 1314 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.mmol_rate_clause) = NULL; }
#line 4712 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 273:
#line 1320 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { if ((yyvsp[0].mmol_rate_clause)) (yyvsp[0].mmol_rate_clause)->next = (yyvsp[-2].mmol_rate_clause); (yyval.mmol_rate_clause) = (yyvsp[0].mmol_rate_clause); }
#line 4718 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 274:
#line 1324 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.mmol_rate_clause) = mdl_assemble_complex_rate_rule_clause(parse_state, parse_state->complex_relations, (yyvsp[-2].str), (yyvsp[-1].ival), &(yyvsp[0].mol_type))); }
#line 4724 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 275:
#line 1327 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.ival) = 0; }
#line 4730 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 276:
#line 1328 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.ival) = 1; }
#line 4736 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 277:
#line 1331 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.sym) = mdl_existing_molecule(parse_state, (yyvsp[0].str))); }
#line 4742 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 278:
#line 1335 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.mol_type) = (yyvsp[0].mol_type); CHECKN((yyval.mol_type).mol_type = mdl_existing_surface_molecule(parse_state, (yyvsp[-1].str))); }
#line 4748 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 279:
#line 1339 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                        (yyval.mol_type) = (yyvsp[0].mol_type);
                                                        if (! (yyval.mol_type).orient_set)
                                                          (yyval.mol_type).orient = 0;
                                                        (yyval.mol_type).mol_type = (yyvsp[-1].sym);
                                                      }
#line 4759 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 280:
#line 1347 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.sym) = mdl_existing_macromolecule(parse_state, (yyvsp[0].str))); }
#line 4765 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 287:
#line 1376 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { mdl_start_surface_class(parse_state, (yyvsp[-1].sym)); }
#line 4771 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 288:
#line 1378 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { mdl_finish_surface_class(parse_state); }
#line 4777 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 289:
#line 1381 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.sym) = mdl_existing_surface_class(parse_state, (yyvsp[0].str))); }
#line 4783 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 295:
#line 1399 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN(mdl_assemble_surface_reaction(parse_state, (yyvsp[-2].tok), parse_state->current_surface_class, (yyvsp[0].mol_type).mol_type, (yyvsp[0].mol_type).orient)); }
#line 4789 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 296:
#line 1402 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
              struct sym_table *mol_sym = retrieve_sym("ALL_MOLECULES", parse_state->vol->mol_sym_table);
              if(!(yyvsp[0].mol_type).orient_set) (yyvsp[0].mol_type).orient = 0;
              CHECKN(mdl_assemble_surface_reaction(parse_state, (yyvsp[-3].tok), parse_state->current_surface_class, mol_sym, (yyvsp[0].mol_type).orient));}
#line 4798 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 297:
#line 1408 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN(mdl_assemble_concentration_clamp_reaction(parse_state, parse_state->current_surface_class, (yyvsp[-2].mol_type).mol_type, (yyvsp[-2].mol_type).orient, (yyvsp[0].dbl))); }
#line 4804 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 298:
#line 1411 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = RFLCT; }
#line 4810 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 299:
#line 1412 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = TRANSP; }
#line 4816 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 300:
#line 1413 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = SINK; }
#line 4822 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 303:
#line 1420 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->current_surface_class->sm_dat_head = (yyvsp[0].surf_mol_dat_list).sm_head; }
#line 4828 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 304:
#line 1427 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.surf_mol_dat_list) = (yyvsp[-1].surf_mol_dat_list); }
#line 4834 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 305:
#line 1431 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.surf_mol_dat_list) = (yyvsp[-1].surf_mol_dat_list); }
#line 4840 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 306:
#line 1435 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                          (yyvsp[0].surf_mol_dat)->quantity_type = SURFMOLDENS;
                                                          (yyval.surf_mol_dat_list).sm_tail = (yyval.surf_mol_dat_list).sm_head = (yyvsp[0].surf_mol_dat);
                                                      }
#line 4849 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 307:
#line 1440 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                          (yyval.surf_mol_dat_list) = (yyvsp[-1].surf_mol_dat_list);
                                                          (yyvsp[0].surf_mol_dat)->quantity_type = SURFMOLDENS;
                                                          (yyval.surf_mol_dat_list).sm_tail = (yyval.surf_mol_dat_list).sm_tail->next = (yyvsp[0].surf_mol_dat);
                                                      }
#line 4859 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 308:
#line 1448 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                          (yyvsp[0].surf_mol_dat)->quantity_type = SURFMOLNUM;
                                                          (yyval.surf_mol_dat_list).sm_tail = (yyval.surf_mol_dat_list).sm_head = (yyvsp[0].surf_mol_dat);
                                                      }
#line 4868 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 309:
#line 1453 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                          (yyval.surf_mol_dat_list) = (yyvsp[-1].surf_mol_dat_list);
                                                          (yyvsp[0].surf_mol_dat)->quantity_type = SURFMOLNUM;
                                                          (yyval.surf_mol_dat_list).sm_tail = (yyval.surf_mol_dat_list).sm_tail->next = (yyvsp[0].surf_mol_dat);
                                                      }
#line 4878 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 310:
#line 1461 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.surf_mol_dat) = mdl_new_surf_mol_data(parse_state, &(yyvsp[-2].mol_type), (yyvsp[0].dbl))); }
#line 4884 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 311:
#line 1465 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->current_surface_class->region_viz_value = (int) (yyvsp[0].dbl); }
#line 4890 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 318:
#line 1490 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { free((yyvsp[0].str)); }
#line 4896 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 326:
#line 1507 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.react_arrow).catalyst = (yyvsp[-1].mol_type); (yyval.react_arrow).flags = ARROW_CATALYTIC; }
#line 4902 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 327:
#line 1512 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.react_arrow).catalyst = (yyvsp[-1].mol_type); (yyval.react_arrow).flags = ARROW_CATALYTIC | ARROW_BIDIRECTIONAL; }
#line 4908 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 328:
#line 1517 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.react_arrow).catalyst.mol_type = NULL; (yyval.react_arrow).flags = 0; }
#line 4914 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 330:
#line 1519 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.react_arrow).catalyst.mol_type = NULL; (yyval.react_arrow).flags = ARROW_BIDIRECTIONAL; }
#line 4920 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 332:
#line 1523 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.sym) = NULL; }
#line 4926 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 333:
#line 1524 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.sym) = mdl_new_rxn_pathname(parse_state, (yyvsp[0].str))); }
#line 4932 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 334:
#line 1530 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN(mdl_assemble_reaction(parse_state, (yyvsp[-5].mol_type_list).mol_type_head, &(yyvsp[-4].mol_type), &(yyvsp[-3].react_arrow), (yyvsp[-2].mol_type_list).mol_type_head, &(yyvsp[-1].react_rates), (yyvsp[0].sym))); }
#line 4938 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 335:
#line 1533 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_reaction_player_singleton(parse_state, & (yyval.mol_type_list), & (yyvsp[0].mol_type))); }
#line 4944 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 336:
#line 1534 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.mol_type_list) = (yyvsp[-2].mol_type_list); CHECK(mdl_add_reaction_player(parse_state, & (yyval.mol_type_list), & (yyvsp[0].mol_type))); }
#line 4950 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 338:
#line 1541 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.mol_type) = (yyvsp[0].mol_type); (yyval.mol_type).is_subunit = 0; }
#line 4956 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 339:
#line 1542 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.mol_type) = (yyvsp[-1].mol_type); (yyval.mol_type).is_subunit = 1; }
#line 4962 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 340:
#line 1546 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.mol_type).mol_type = NULL; }
#line 4968 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 341:
#line 1547 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.mol_type) = (yyvsp[0].mol_type); }
#line 4974 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 342:
#line 1551 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.mol_type) = (yyvsp[0].mol_type); (yyval.mol_type).mol_type = (yyvsp[-1].sym); }
#line 4980 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 343:
#line 1554 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_reaction_player_singleton(parse_state, & (yyval.mol_type_list), & (yyvsp[0].mol_type))); }
#line 4986 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 344:
#line 1555 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.mol_type_list) = (yyvsp[-2].mol_type_list); CHECK(mdl_add_reaction_player(parse_state, & (yyval.mol_type_list), & (yyvsp[0].mol_type))); }
#line 4992 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 345:
#line 1558 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.mol_type).mol_type = NULL; (yyval.mol_type).orient_set = 0; }
#line 4998 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 349:
#line 1567 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                        if ((yyvsp[-1].react_rates).forward_rate.rate_type == RATE_UNSET)
                                                        {
                                                          mdlerror(parse_state, "Invalid reaction rate specification: must specify a forward rate.");
                                                          return 1;
                                                        }

                                                        (yyval.react_rates) = (yyvsp[-1].react_rates);
                                                      }
#line 5012 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 350:
#line 1578 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                        if (((yyvsp[-3].react_rates).forward_rate.rate_type  != RATE_UNSET && (yyvsp[-1].react_rates).forward_rate.rate_type  != RATE_UNSET)  ||
                                                            ((yyvsp[-3].react_rates).backward_rate.rate_type != RATE_UNSET && (yyvsp[-1].react_rates).backward_rate.rate_type != RATE_UNSET))
                                                        {
                                                          mdlerror_fmt(parse_state, "Error: When two reaction rates are specified, one must be a forward rate, and one must be a reverse rate");
                                                          return 1;
                                                        }

                                                        (yyval.react_rates) = (yyvsp[-3].react_rates);
                                                        if ((yyvsp[-1].react_rates).forward_rate.rate_type != RATE_UNSET)
                                                          (yyval.react_rates).forward_rate = (yyvsp[-1].react_rates).forward_rate;
                                                        else
                                                          (yyval.react_rates).backward_rate = (yyvsp[-1].react_rates).backward_rate;
                                                      }
#line 5031 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 351:
#line 1595 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.react_rates).forward_rate = (yyvsp[0].react_rate); (yyval.react_rates).backward_rate.rate_type = RATE_UNSET; CHECK(mdl_valid_rate(parse_state, &(yyvsp[0].react_rate))); }
#line 5037 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 352:
#line 1596 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.react_rates).forward_rate = (yyvsp[0].react_rate); (yyval.react_rates).backward_rate.rate_type = RATE_UNSET; CHECK(mdl_valid_rate(parse_state, &(yyvsp[0].react_rate))); }
#line 5043 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 353:
#line 1597 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.react_rates).backward_rate = (yyvsp[0].react_rate); (yyval.react_rates).forward_rate.rate_type = RATE_UNSET; CHECK(mdl_valid_rate(parse_state, &(yyvsp[0].react_rate))); }
#line 5049 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 354:
#line 1601 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.react_rate).rate_type = RATE_CONSTANT; (yyval.react_rate).v.rate_constant = (yyvsp[0].dbl); }
#line 5055 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 355:
#line 1602 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.react_rate).rate_type = RATE_FILE; (yyval.react_rate).v.rate_file = (yyvsp[0].str); }
#line 5061 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 356:
#line 1603 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_reaction_rate_from_var(parse_state, & (yyval.react_rate), (yyvsp[0].sym))); }
#line 5067 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 357:
#line 1604 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_reaction_rate_complex(parse_state, & (yyval.react_rate), (yyvsp[-1].sym), (yyvsp[0].str))); }
#line 5073 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 358:
#line 1615 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_release_pattern(parse_state, (yyvsp[-3].sym), &(yyvsp[-1].rpat))); }
#line 5079 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 359:
#line 1618 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.sym) = mdl_new_release_pattern(parse_state, (yyvsp[0].str))); }
#line 5085 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 360:
#line 1621 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.sym) = mdl_existing_release_pattern_or_rxn_pathname(parse_state, (yyvsp[0].str))); }
#line 5091 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 361:
#line 1625 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                        (yyval.rpat).delay = 0;
                                                        (yyval.rpat).release_interval = FOREVER;
                                                        (yyval.rpat).train_interval = FOREVER;
                                                        (yyval.rpat).train_duration = FOREVER;
                                                        (yyval.rpat).number_of_trains = 1;
                                                      }
#line 5103 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 362:
#line 1633 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.rpat) = (yyvsp[-3].rpat); (yyval.rpat).delay = (yyvsp[0].dbl) / parse_state->vol->time_unit; }
#line 5109 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 363:
#line 1635 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.rpat) = (yyvsp[-3].rpat); (yyval.rpat).release_interval = (yyvsp[0].dbl) / parse_state->vol->time_unit; }
#line 5115 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 364:
#line 1637 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.rpat) = (yyvsp[-3].rpat); (yyval.rpat).train_interval = (yyvsp[0].dbl) / parse_state->vol->time_unit; }
#line 5121 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 365:
#line 1639 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.rpat) = (yyvsp[-3].rpat); (yyval.rpat).train_duration = (yyvsp[0].dbl) / parse_state->vol->time_unit; }
#line 5127 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 366:
#line 1641 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.rpat) = (yyvsp[-3].rpat); (yyval.rpat).number_of_trains = (yyvsp[0].ival); }
#line 5133 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 367:
#line 1644 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.ival) = (int) (yyvsp[0].dbl); }
#line 5139 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 368:
#line 1645 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.ival) = INT_MAX; }
#line 5145 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 369:
#line 1652 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->current_object = parse_state->vol->root_instance; }
#line 5151 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 370:
#line 1653 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                        check_regions(parse_state->vol->root_instance, (yyvsp[0].obj));
                                                        add_child_objects(parse_state->vol->root_instance, (yyvsp[0].obj), (yyvsp[0].obj));
                                                        parse_state->current_object = parse_state->vol->root_object;
                                                      }
#line 5161 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 371:
#line 1663 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { add_child_objects(parse_state->vol->root_object, (yyvsp[0].obj), (yyvsp[0].obj)); }
#line 5167 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 377:
#line 1679 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.sym) = mdl_start_object(parse_state, (yyvsp[0].str))); }
#line 5173 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 379:
#line 1685 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { mdl_finish_object(parse_state); }
#line 5179 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 383:
#line 1698 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { transform_translate(parse_state->vol, parse_state->current_object->t_matrix, (yyvsp[0].vec3)); }
#line 5185 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 384:
#line 1699 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { transform_scale(parse_state->current_object->t_matrix, (yyvsp[0].vec3)); }
#line 5191 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 385:
#line 1700 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_transform_rotate(parse_state, parse_state->current_object->t_matrix, (yyvsp[-2].vec3), (yyvsp[0].dbl))); }
#line 5197 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 386:
#line 1709 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                          struct object *the_object = (struct object *) (yyvsp[-5].sym)->value;
                                                          the_object->object_type = META_OBJ;
                                                          add_child_objects(the_object, (yyvsp[-2].obj_list).obj_head, (yyvsp[-2].obj_list).obj_tail);
                                                          (yyval.obj) = the_object;
                                                      }
#line 5208 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 387:
#line 1718 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { mdl_object_list_singleton(& (yyval.obj_list), (yyvsp[0].obj)); }
#line 5214 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 388:
#line 1719 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.obj_list) = (yyvsp[-1].obj_list); mdl_add_object_to_list(& (yyval.obj_list), (yyvsp[0].obj)); }
#line 5220 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 391:
#line 1728 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_deep_copy_object(parse_state, (struct object *) (yyvsp[-3].sym)->value, (struct object *) (yyvsp[-1].sym)->value)); }
#line 5226 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 392:
#line 1730 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.obj) = (struct object *) (yyvsp[-6].sym)->value; }
#line 5232 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 395:
#line 1740 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_start_release_site(parse_state, (yyvsp[-2].sym), SHAPE_UNDEFINED)); }
#line 5238 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 396:
#line 1744 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.obj) = mdl_finish_release_site(parse_state, (yyvsp[-7].sym))); }
#line 5244 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 397:
#line 1747 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_release_site_geometry_region(parse_state, parse_state->current_release_site, parse_state->current_object, (yyvsp[0].rev))); }
#line 5250 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 398:
#line 1748 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_release_site_geometry_object(parse_state, parse_state->current_release_site, (struct object *) (yyvsp[0].sym)->value)); }
#line 5256 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 399:
#line 1749 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->current_release_site->release_shape = SHAPE_SPHERICAL; }
#line 5262 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 400:
#line 1750 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->current_release_site->release_shape = SHAPE_CUBIC; }
#line 5268 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 401:
#line 1751 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->current_release_site->release_shape = SHAPE_ELLIPTIC; }
#line 5274 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 402:
#line 1752 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->current_release_site->release_shape = SHAPE_RECTANGULAR; }
#line 5280 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 403:
#line 1753 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->current_release_site->release_shape = SHAPE_SPHERICAL_SHELL; }
#line 5286 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 404:
#line 1754 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                          parse_state->current_release_site->release_shape = SHAPE_LIST;
                                                          parse_state->current_release_site->release_number_method = CONSTNUM;
                                                      }
#line 5295 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 405:
#line 1761 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.rev) = new_release_region_expr_term((yyvsp[0].sym))); }
#line 5301 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 406:
#line 1762 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.rev) = (yyvsp[-1].rev); }
#line 5307 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 407:
#line 1763 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.rev) = new_release_region_expr_binary((yyvsp[-2].rev), (yyvsp[0].rev), REXP_UNION)); }
#line 5313 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 408:
#line 1764 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.rev) = new_release_region_expr_binary((yyvsp[-2].rev), (yyvsp[0].rev), REXP_SUBTRACTION)); }
#line 5319 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 409:
#line 1765 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.rev) = new_release_region_expr_binary((yyvsp[-2].rev), (yyvsp[0].rev), REXP_INTERSECTION)); }
#line 5325 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 410:
#line 1770 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_start_release_site(parse_state, (yyvsp[-2].sym), (yyvsp[-1].tok))); }
#line 5331 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 411:
#line 1773 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.obj) = mdl_finish_release_site(parse_state, (yyvsp[-6].sym))); }
#line 5337 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 412:
#line 1776 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = SHAPE_SPHERICAL; }
#line 5343 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 413:
#line 1777 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = SHAPE_CUBIC; }
#line 5349 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 414:
#line 1778 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = SHAPE_ELLIPTIC; }
#line 5355 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 415:
#line 1779 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = SHAPE_RECTANGULAR; }
#line 5361 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 416:
#line 1780 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = SHAPE_SPHERICAL_SHELL; }
#line 5367 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 419:
#line 1788 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.sym) = mdl_existing_num_or_array(parse_state, (yyvsp[0].str))); }
#line 5373 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 420:
#line 1792 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { set_release_site_location(parse_state->vol, parse_state->current_release_site, (yyvsp[0].vec3)); }
#line 5379 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 421:
#line 1793 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_release_site_molecule(parse_state, parse_state->current_release_site, & (yyvsp[0].mol_type))); }
#line 5385 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 422:
#line 1794 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                        if (parse_state->current_release_site->release_shape == SHAPE_LIST)
                                                        {
                                                          mdlerror(parse_state, "Molecules are already specified in a list--cannot set number or density.");
                                                          return 1;
                                                        }
                                                      }
#line 5397 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 423:
#line 1801 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_release_site_diameter(parse_state, parse_state->current_release_site, (yyvsp[0].dbl) * (((yyvsp[-2].tok) == SITE_RADIUS) ? 2.0 : 1.0))); }
#line 5403 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 424:
#line 1802 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_release_site_diameter_array(parse_state, parse_state->current_release_site, (yyvsp[0].nlist).value_count, (yyvsp[0].nlist).value_head, ((yyvsp[-2].tok) == SITE_RADIUS) ? 2.0 : 1.0)); }
#line 5409 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 425:
#line 1803 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_release_site_diameter_var(parse_state, parse_state->current_release_site, ((yyvsp[-2].tok) == SITE_RADIUS) ? 2.0 : 1.0, (yyvsp[0].sym))); }
#line 5415 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 426:
#line 1804 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_release_site_probability(parse_state, parse_state->current_release_site, (yyvsp[0].dbl))); }
#line 5421 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 427:
#line 1806 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_release_site_pattern(parse_state, parse_state->current_release_site, (yyvsp[0].sym))); }
#line 5427 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 428:
#line 1808 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_release_site_molecule_positions(parse_state, parse_state->current_release_site, & (yyvsp[-1].rsm_list))); }
#line 5433 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 429:
#line 1812 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = SITE_DIAMETER; }
#line 5439 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 430:
#line 1813 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = SITE_RADIUS; }
#line 5445 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 435:
#line 1825 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { set_release_site_constant_number(parse_state->current_release_site, (yyvsp[0].dbl)); }
#line 5451 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 436:
#line 1828 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { set_release_site_constant_number(parse_state->current_release_site, (yyvsp[-1].dbl)); }
#line 5457 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 437:
#line 1835 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { set_release_site_gaussian_number(parse_state->current_release_site, (yyvsp[-4].dbl), (yyvsp[-1].dbl)); }
#line 5463 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 438:
#line 1843 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { set_release_site_volume_dependent_number(parse_state->current_release_site, (yyvsp[-7].dbl), (yyvsp[-4].dbl), (yyvsp[-1].dbl)); }
#line 5469 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 439:
#line 1847 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_release_site_concentration(parse_state, parse_state->current_release_site, (yyvsp[0].dbl))); }
#line 5475 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 440:
#line 1848 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(set_release_site_density(parse_state->current_release_site, (yyvsp[0].dbl))); }
#line 5481 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 441:
#line 1852 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { release_single_molecule_singleton(& (yyval.rsm_list), (yyvsp[0].rsm)); }
#line 5487 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 442:
#line 1854 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.rsm_list) = (yyvsp[-1].rsm_list); add_release_single_molecule_to_list(& (yyval.rsm_list), (yyvsp[0].rsm)); }
#line 5493 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 443:
#line 1858 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.rsm) = mdl_new_release_single_molecule(parse_state, &(yyvsp[-1].mol_type), (yyvsp[0].vec3))); }
#line 5499 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 445:
#line 1869 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                        CHECKN((yyval.obj) = mdl_new_polygon_list(
                                                          parse_state, (yyvsp[-4].str), (yyvsp[-1].vertlist).vertex_count, (yyvsp[-1].vertlist).vertex_head,
                                                          (yyvsp[0].ecl).connection_count, (yyvsp[0].ecl).connection_head));
                                                      }
#line 5509 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 446:
#line 1878 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                          (yyval.obj) = (struct object *) (yyvsp[-3].obj);
                                                          CHECK(mdl_finish_polygon_list(parse_state, (yyval.obj)));
                                                      }
#line 5518 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 447:
#line 1884 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.vertlist) = (yyvsp[-1].vertlist); }
#line 5524 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 448:
#line 1887 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.vertlistitem) = mdl_new_vertex_list_item((yyvsp[0].vec3))); }
#line 5530 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 449:
#line 1890 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { mdl_vertex_list_singleton(& (yyval.vertlist), (yyvsp[0].vertlistitem)); }
#line 5536 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 450:
#line 1891 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.vertlist) = (yyvsp[-1].vertlist); mdl_add_vertex_to_list(& (yyval.vertlist), (yyvsp[0].vertlistitem)); }
#line 5542 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 451:
#line 1896 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.ecl) = (yyvsp[-1].ecl); }
#line 5548 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 452:
#line 1900 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { mdl_element_connection_list_singleton(& (yyval.ecl), (yyvsp[0].elem_conn)); }
#line 5554 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 453:
#line 1902 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.ecl) = (yyvsp[-1].ecl); mdl_add_element_connection_to_list(& (yyval.ecl), (yyvsp[0].elem_conn)); }
#line 5560 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 454:
#line 1905 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.elem_conn) = mdl_new_element_connection(parse_state, & (yyvsp[0].nlist))); }
#line 5566 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 459:
#line 1921 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN(parse_state->current_region = mdl_get_region(parse_state, parse_state->current_object, "REMOVED")); }
#line 5572 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 460:
#line 1923 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                          parse_state->current_region->element_list_head = (yyvsp[-1].elem_list).elml_head;
                                                          if (parse_state->current_object->object_type == POLY_OBJ)
                                                          {
                                                            CHECK(mdl_normalize_elements(parse_state, parse_state->current_region,0));
                                                          }
                                                      }
#line 5584 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 463:
#line 1937 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = Z_POS; }
#line 5590 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 464:
#line 1938 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = Z_NEG; }
#line 5596 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 465:
#line 1939 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = Y_NEG; }
#line 5602 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 466:
#line 1940 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = Y_POS; }
#line 5608 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 467:
#line 1941 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = X_NEG; }
#line 5614 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 468:
#line 1942 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = X_POS; }
#line 5620 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 469:
#line 1943 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = ALL_SIDES; }
#line 5626 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 471:
#line 1949 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.elem_list) = (yyvsp[-1].elem_list); mdl_add_elements_to_list(& (yyval.elem_list), (yyvsp[0].elem_list).elml_head, (yyvsp[0].elem_list).elml_tail); }
#line 5632 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 474:
#line 1955 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.elem_list).elml_tail = (yyval.elem_list).elml_head = (yyvsp[0].elem_list_item); }
#line 5638 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 475:
#line 1956 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.elem_list).elml_tail = (yyval.elem_list).elml_head = (yyvsp[0].elem_list_item); }
#line 5644 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 476:
#line 1961 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.elem_list) = (yyvsp[-1].elem_list); }
#line 5650 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 477:
#line 1966 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.elem_list) = (yyvsp[-1].elem_list); mdl_set_elements_to_exclude((yyval.elem_list).elml_head); }
#line 5656 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 479:
#line 1973 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.elem_list).elml_tail = (yyval.elem_list).elml_head = (yyvsp[0].elem_list_item); }
#line 5662 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 480:
#line 1974 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.elem_list) = (yyvsp[-2].elem_list); mdl_add_elements_to_list(& (yyval.elem_list), (yyvsp[0].elem_list_item), (yyvsp[0].elem_list_item)); }
#line 5668 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 481:
#line 1977 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.elem_list_item) = new_element_list((unsigned int) (yyvsp[0].dbl), (unsigned int) (yyvsp[0].dbl))); }
#line 5674 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 482:
#line 1978 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.elem_list_item) = new_element_list((unsigned int) (yyvsp[-2].dbl), (unsigned int) (yyvsp[0].dbl))); }
#line 5680 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 483:
#line 1979 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.elem_list_item) = mdl_new_element_side(parse_state, (yyvsp[0].tok))); }
#line 5686 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 484:
#line 1982 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.elem_list_item) = mdl_new_element_previous_region(parse_state, parse_state->current_object, parse_state->current_region, (yyvsp[0].str), (yyvsp[-2].tok))); }
#line 5692 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 485:
#line 1985 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = 0; }
#line 5698 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 486:
#line 1986 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = 1; }
#line 5704 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 487:
#line 1989 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.elem_list_item) = mdl_new_element_patch(parse_state, parse_state->current_polygon, (yyvsp[-2].vec3), (yyvsp[0].vec3), (yyvsp[-4].tok))); }
#line 5710 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 488:
#line 1992 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = 0; }
#line 5716 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 489:
#line 1993 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = 1; }
#line 5722 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 493:
#line 2009 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->current_region = (yyvsp[-1].reg); }
#line 5728 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 494:
#line 2010 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_region_elements(parse_state, (yyvsp[-3].reg), (yyvsp[0].elem_list).elml_head, (yyvsp[-3].reg)->parent->object_type == POLY_OBJ)); }
#line 5734 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 495:
#line 2012 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->current_region = NULL; }
#line 5740 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 496:
#line 2020 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                        CHECKN(mdl_new_voxel_list(parse_state, (yyvsp[-4].sym),
                                                                                  (yyvsp[-1].vertlist).vertex_count, (yyvsp[-1].vertlist).vertex_head,
                                                                                  (yyvsp[0].ecl).connection_count, (yyvsp[0].ecl).connection_head));
                                                      }
#line 5750 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 497:
#line 2026 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.obj) = (struct object *) (yyvsp[-7].sym)->value; }
#line 5756 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 498:
#line 2031 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.ecl) = (yyvsp[-1].ecl); }
#line 5762 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 499:
#line 2034 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.elem_conn) = mdl_new_tet_element_connection(parse_state, & (yyvsp[0].nlist))); }
#line 5768 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 500:
#line 2038 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                          (yyval.ecl).connection_head = (yyval.ecl).connection_tail = (yyvsp[0].elem_conn);
                                                          (yyval.ecl).connection_count = 1;
                                                      }
#line 5777 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 501:
#line 2042 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                          (yyval.ecl) = (yyvsp[-1].ecl);
                                                          (yyval.ecl).connection_tail = (yyval.ecl).connection_tail->next = (yyvsp[0].elem_conn);
                                                          ++ (yyval.ecl).connection_count;
                                                      }
#line 5787 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 502:
#line 2053 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN(mdl_new_box_object(parse_state, (yyvsp[-8].sym), (yyvsp[-3].vec3), (yyvsp[-1].vec3))); }
#line 5793 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 503:
#line 2054 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_triangulate_box_object(parse_state, (yyvsp[-10].sym), parse_state->current_polygon, (yyvsp[-2].dbl))); }
#line 5799 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 504:
#line 2056 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                          CHECK(mdl_finish_box_object(parse_state, (yyvsp[-13].sym)));
                                                          (yyval.obj) = (struct object *) (yyvsp[-13].sym)->value;
                                                      }
#line 5808 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 505:
#line 2062 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.dbl) = 0.0; }
#line 5814 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 506:
#line 2063 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                        (yyval.dbl) = (yyvsp[0].dbl);
                                                        if ((yyval.dbl) < 2.0)
                                                        {
                                                          mdlerror(parse_state, "Invalid aspect ratio requested (must be greater than or equal to 2.0)");
                                                          return 1;
                                                        }
                                                      }
#line 5827 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 510:
#line 2089 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_start_existing_obj_region_def(parse_state, (yyvsp[0].sym))); }
#line 5833 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 511:
#line 2090 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->current_region = (yyvsp[-1].reg); }
#line 5839 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 512:
#line 2092 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { mdl_set_region_elements(parse_state, (yyvsp[-4].reg), (yyvsp[0].elem_list).elml_head, 1); }
#line 5845 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 513:
#line 2094 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                          parse_state->current_region = NULL;
                                                          parse_state->current_polygon = NULL;
                                                          parse_state->current_object = parse_state->vol->root_object;
                                                      }
#line 5855 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 514:
#line 2101 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.reg) = mdl_create_region(parse_state, parse_state->current_object, (yyvsp[0].str))); }
#line 5861 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 518:
#line 2112 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { mdl_add_surf_mol_to_region(parse_state->current_region, & (yyvsp[0].surf_mol_dat_list)); }
#line 5867 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 520:
#line 2117 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { mdl_set_region_surface_class(parse_state, parse_state->current_region, (yyvsp[0].sym)); }
#line 5873 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 521:
#line 2121 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { mdl_set_region_region_viz_value(parse_state, parse_state->current_region, (int) (yyvsp[0].dbl)); }
#line 5879 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 525:
#line 2140 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->current_region = (struct region *) (yyvsp[-1].sym)->value; }
#line 5885 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 526:
#line 2142 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->current_region = NULL; }
#line 5891 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 527:
#line 2150 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                          parse_state->header_comment = NULL;  /* No header by default */
                                                          parse_state->exact_time_flag = 1;    /* Print exact_time column in TRIGGER output by default */
                                                      }
#line 5900 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 528:
#line 2156 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_add_reaction_output_block_to_world(parse_state, (int) (yyvsp[-4].dbl), & (yyvsp[-2].ro_otimes), & (yyvsp[-1].ro_sets))); }
#line 5906 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 529:
#line 2160 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.dbl) = COUNTBUFFERSIZE; }
#line 5912 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 530:
#line 2161 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                          double temp_value = (yyvsp[0].dbl);
                                                          if (!(temp_value >= 1.0 && temp_value < UINT_MAX))
                                                          {
                                                            mdlerror_fmt(parse_state, "Requested buffer size of %.15g lines is invalid.  Suggested range is 100-1000000.", temp_value);
                                                            return 1;
                                                          }
                                                          (yyval.dbl) = (yyvsp[0].dbl);
                                                      }
#line 5926 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 534:
#line 2177 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.ro_otimes).type = OUTPUT_BY_STEP; (yyval.ro_otimes).step = (yyvsp[0].dbl); }
#line 5932 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 535:
#line 2181 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                        (yyval.ro_otimes).type = OUTPUT_BY_ITERATION_LIST;
                                                        (yyval.ro_otimes).values = (yyvsp[0].nlist);
                                                      }
#line 5941 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 536:
#line 2189 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                        (yyval.ro_otimes).type = OUTPUT_BY_TIME_LIST;
                                                        (yyval.ro_otimes).values = (yyvsp[0].nlist);
                                                      }
#line 5950 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 537:
#line 2196 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.ro_sets).set_head = (yyval.ro_sets).set_tail = (yyvsp[0].ro_set); }
#line 5956 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 538:
#line 2198 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                        (yyval.ro_sets) = (yyvsp[-1].ro_sets);
                                                        if ((yyvsp[0].ro_set) != NULL)
                                                        {
                                                          if ((yyval.ro_sets).set_tail != NULL)
                                                            (yyval.ro_sets).set_tail = (yyval.ro_sets).set_tail->next = (yyvsp[0].ro_set);
                                                          else
                                                            (yyval.ro_sets).set_tail = (yyval.ro_sets).set_head = (yyvsp[0].ro_set);
                                                        }
                                                      }
#line 5971 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 540:
#line 2212 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.ro_set) = NULL; }
#line 5977 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 541:
#line 2213 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.ro_set) = NULL; }
#line 5983 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 542:
#line 2217 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {  parse_state->count_flags = 0; }
#line 5989 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 543:
#line 2219 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.ro_set) = mdl_populate_output_set(parse_state, parse_state->header_comment, parse_state->exact_time_flag, (yyvsp[-3].ro_cols).column_head, (yyvsp[-1].tok), (yyvsp[0].str))); }
#line 5995 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 544:
#line 2223 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.str) = NULL; }
#line 6001 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 545:
#line 2224 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.str) = ((yyvsp[0].tok) ? "" : NULL); }
#line 6007 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 546:
#line 2225 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.str) = (yyvsp[0].str); }
#line 6013 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 547:
#line 2229 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->header_comment = (yyvsp[0].str); }
#line 6019 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 548:
#line 2233 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->exact_time_flag = (yyvsp[0].tok); }
#line 6025 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 550:
#line 2239 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                          (yyval.ro_cols) = (yyvsp[-2].ro_cols);
                                                          (yyval.ro_cols).column_tail->next = (yyvsp[0].ro_cols).column_head;
                                                          (yyval.ro_cols).column_tail = (yyvsp[0].ro_cols).column_tail;
                                                      }
#line 6035 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 551:
#line 2247 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_single_count_expr(parse_state, & (yyval.ro_cols), (yyvsp[-1].cnt), (yyvsp[0].str))); }
#line 6041 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 552:
#line 2251 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.cnt) = mdl_new_oexpr_constant(parse_state, (yyvsp[0].dbl))); }
#line 6047 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 554:
#line 2253 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.cnt) = mdl_join_oexpr_tree(parse_state, (yyvsp[-1].cnt), NULL, '(')); }
#line 6053 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 555:
#line 2254 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.cnt) = mdl_join_oexpr_tree(parse_state, (yyvsp[-2].cnt),   (yyvsp[0].cnt), '+')); }
#line 6059 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 556:
#line 2255 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.cnt) = mdl_join_oexpr_tree(parse_state, (yyvsp[-2].cnt),   (yyvsp[0].cnt), '-')); }
#line 6065 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 557:
#line 2256 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.cnt) = mdl_join_oexpr_tree(parse_state, (yyvsp[-2].cnt),   (yyvsp[0].cnt), '*')); }
#line 6071 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 558:
#line 2257 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.cnt) = mdl_join_oexpr_tree(parse_state, (yyvsp[-2].cnt),   (yyvsp[0].cnt), '/')); }
#line 6077 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 559:
#line 2258 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.cnt) = mdl_join_oexpr_tree(parse_state, (yyvsp[0].cnt), NULL, '_')); }
#line 6083 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 560:
#line 2259 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.cnt) = mdl_sum_oexpr((yyvsp[-1].cnt))); }
#line 6089 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 561:
#line 2264 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->count_flags |= COUNT_PRESENT; }
#line 6095 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 562:
#line 2265 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.cnt) = (yyvsp[-1].cnt); }
#line 6101 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 563:
#line 2266 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.cnt) = mdl_new_oexpr_constant(parse_state, (yyvsp[-1].dbl))); }
#line 6107 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 564:
#line 2267 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->count_flags |= TRIGGER_PRESENT; }
#line 6113 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 565:
#line 2268 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.cnt) = (yyvsp[-1].cnt); }
#line 6119 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 566:
#line 2271 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = FILE_OVERWRITE; }
#line 6125 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 567:
#line 2272 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = FILE_SUBSTITUTE; }
#line 6131 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 568:
#line 2273 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = FILE_APPEND; }
#line 6137 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 569:
#line 2274 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = FILE_APPEND_HEADER; }
#line 6143 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 570:
#line 2275 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = FILE_CREATE; }
#line 6149 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 572:
#line 2281 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.sym) = mdl_existing_rxn_pathname_or_molecule(parse_state, (yyvsp[0].str))); }
#line 6155 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 573:
#line 2285 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                        (yyval.mol_type) = (yyvsp[0].mol_type);
                                                        if ((yyval.mol_type).orient > 0)
                                                          (yyval.mol_type).orient = 1;
                                                        else if ((yyval.mol_type).orient < 0)
                                                          (yyval.mol_type).orient = -1;
                                                        CHECKN((yyval.mol_type).mol_type = mdl_existing_molecule(parse_state, (yyvsp[-1].str)));
                                                      }
#line 6168 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 578:
#line 2303 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.cnt) = mdl_count_syntax_1(parse_state, (yyvsp[-3].sym), (yyvsp[-1].sym), (yyvsp[0].tok), parse_state->count_flags)); }
#line 6174 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 579:
#line 2308 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.cnt) = mdl_count_syntax_2(parse_state, (yyvsp[-3].mol_type).mol_type, (yyvsp[-3].mol_type).orient, (yyvsp[-1].sym), (yyvsp[0].tok), parse_state->count_flags)); }
#line 6180 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 580:
#line 2313 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.cnt) = mdl_count_syntax_3(parse_state, (yyvsp[-3].str), (yyvsp[-1].sym), (yyvsp[0].tok), parse_state->count_flags)); }
#line 6186 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 582:
#line 2320 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                          if (((struct species *) (yyvsp[0].mol_type).mol_type->value)->flags & IS_COMPLEX)
                                                          {
                                                            mdlerror_fmt(parse_state, "The molecule '%s' is a complex, and may not be used as a subunit type", (yyvsp[0].mol_type).mol_type->name);
                                                            return 1;
                                                          }
                                                      }
#line 6198 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 583:
#line 2331 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->current_complex = (struct complex_species *) (yyvsp[0].sym)->value; }
#line 6204 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 584:
#line 2335 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                          parse_state->current_complex = NULL;
                                                          struct complex_species *macromol = (struct complex_species *) (yyvsp[-8].sym)->value;
                                                          struct mcell_species master_orientation = (yyvsp[-6].mol_type);
                                                          struct mcell_species subunit = (yyvsp[-4].mol_type);
                                                          struct macro_relation_state *relation_states = (yyvsp[-3].relation_state);
                                                          struct sym_table *location = (yyvsp[0].sym);
                                                          CHECKN((yyval.cnt) = mdl_count_syntax_macromol_subunit(parse_state, macromol, &master_orientation, & subunit, relation_states, location));
                                                      }
#line 6218 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 585:
#line 2347 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.relation_state) = NULL; }
#line 6224 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 586:
#line 2348 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.relation_state) = (yyvsp[-1].relation_state); }
#line 6230 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 588:
#line 2354 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyvsp[0].relation_state)->next = (yyvsp[-2].relation_state); (yyval.relation_state) = (yyvsp[0].relation_state); }
#line 6236 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 589:
#line 2360 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.relation_state) = mdl_assemble_complex_relation_state(parse_state, (yyvsp[-2].ival), (yyvsp[-1].ival), & (yyvsp[0].mol_type))); }
#line 6242 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 590:
#line 2363 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                          int rel_idx = macro_lookup_relation(parse_state->current_complex, (yyvsp[0].str));
                                                          if (rel_idx == -1)
                                                          {
                                                            mdlerror_fmt(parse_state,
                                                                         "In subunit specification for COUNT statement, relation '%s' does not exist within the complex '%s'",
                                                                         (yyvsp[0].str),
                                                                         parse_state->current_complex->base.sym->name);
                                                            return 1;
                                                          }

                                                          (yyval.ival) = rel_idx;
                                                      }
#line 6260 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 591:
#line 2378 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.sym) = NULL; }
#line 6266 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 592:
#line 2379 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.sym) = (yyvsp[0].sym); }
#line 6272 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 593:
#line 2380 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.sym) = (yyvsp[0].sym); }
#line 6278 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 594:
#line 2383 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = REPORT_NOTHING; }
#line 6284 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 595:
#line 2384 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = (yyvsp[0].tok); }
#line 6290 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 596:
#line 2387 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = REPORT_FRONT_HITS; }
#line 6296 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 597:
#line 2388 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = REPORT_BACK_HITS; }
#line 6302 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 598:
#line 2389 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = REPORT_ALL_HITS; }
#line 6308 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 599:
#line 2390 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = REPORT_FRONT_CROSSINGS; }
#line 6314 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 600:
#line 2391 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = REPORT_BACK_CROSSINGS; }
#line 6320 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 601:
#line 2392 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = REPORT_ALL_CROSSINGS; }
#line 6326 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 602:
#line 2393 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = REPORT_CONCENTRATION; }
#line 6332 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 603:
#line 2394 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = REPORT_ENCLOSED; }
#line 6338 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 604:
#line 2397 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.str) = NULL; }
#line 6344 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 605:
#line 2398 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.str) = (yyvsp[0].str); }
#line 6350 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 606:
#line 2405 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_new_viz_output_block(parse_state)); }
#line 6356 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 607:
#line 2410 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_finish_viz_output_block(parse_state, parse_state->vol->viz_blocks)); }
#line 6362 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 610:
#line 2419 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_viz_mode(parse_state->vol->viz_blocks, CELLBLENDER_MODE)); }
#line 6368 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 611:
#line 2420 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_viz_mode(parse_state->vol->viz_blocks, (yyvsp[0].ival))); }
#line 6374 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 612:
#line 2423 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.ival) = NO_VIZ_MODE; }
#line 6380 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 613:
#line 2424 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.ival) = DREAMM_V3_MODE; }
#line 6386 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 614:
#line 2425 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.ival) = DREAMM_V3_GROUPED_MODE; }
#line 6392 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 615:
#line 2426 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.ival) = ASCII_MODE; }
#line 6398 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 616:
#line 2427 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.ival) = CELLBLENDER_MODE; }
#line 6404 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 617:
#line 2430 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                        if (parse_state->vol->viz_blocks->viz_mode == DREAMM_V3_MODE)
                                                          CHECK(mdl_set_viz_mesh_format(parse_state, parse_state->vol->viz_blocks, VIZ_MESH_FORMAT_BINARY));
                                                      }
#line 6413 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 618:
#line 2434 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_viz_mesh_format(parse_state, parse_state->vol->viz_blocks, (yyvsp[0].ival))); }
#line 6419 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 619:
#line 2437 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.ival) = VIZ_MESH_FORMAT_BINARY; }
#line 6425 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 620:
#line 2438 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.ival) = VIZ_MESH_FORMAT_ASCII; }
#line 6431 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 621:
#line 2442 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                        if (parse_state->vol->viz_blocks->viz_mode == DREAMM_V3_MODE)
                                                          CHECK(mdl_set_viz_molecule_format(parse_state, parse_state->vol->viz_blocks, VIZ_MOLECULE_FORMAT_BINARY));
                                                      }
#line 6440 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 622:
#line 2446 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_viz_molecule_format(parse_state, parse_state->vol->viz_blocks, (yyvsp[0].ival))); }
#line 6446 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 623:
#line 2450 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.ival) = VIZ_MOLECULE_FORMAT_BINARY; }
#line 6452 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 624:
#line 2451 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.ival) = VIZ_MOLECULE_FORMAT_ASCII; }
#line 6458 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 626:
#line 2456 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                        if ((yyvsp[0].frame_list).frame_head)
                                                        {
                                                          (yyvsp[0].frame_list).frame_tail->next = parse_state->vol->viz_blocks->frame_data_head;
                                                          parse_state->vol->viz_blocks->frame_data_head = (yyvsp[0].frame_list).frame_head;
                                                        }
                                                      }
#line 6470 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 629:
#line 2470 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_viz_filename_prefix(parse_state, parse_state->vol->viz_blocks, (yyvsp[0].str))); }
#line 6476 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 630:
#line 2476 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.frame_list) = (yyvsp[-1].frame_list); }
#line 6482 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 632:
#line 2482 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                        (yyval.frame_list) = (yyvsp[-1].frame_list);
                                                        if ((yyval.frame_list).frame_tail)
                                                        {
                                                          (yyval.frame_list).frame_tail->next = (yyvsp[0].frame_list).frame_head;
                                                          if ((yyvsp[0].frame_list).frame_tail)
                                                            (yyval.frame_list).frame_tail = (yyvsp[0].frame_list).frame_tail;
                                                        }
                                                        else
                                                          (yyval.frame_list) = (yyvsp[0].frame_list);
                                                      }
#line 6498 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 633:
#line 2496 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.frame_list).frame_head = (yyval.frame_list).frame_tail = NULL; }
#line 6504 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 637:
#line 2508 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_viz_state(parse_state, & (yyval.ival), (yyvsp[0].dbl))); }
#line 6510 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 638:
#line 2509 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.ival) = INCLUDE_OBJ; }
#line 6516 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 641:
#line 2519 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_viz_include_molecules(parse_state, parse_state->vol->viz_blocks, (yyvsp[-1].symlist), (yyvsp[0].ival))); }
#line 6522 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 642:
#line 2520 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_viz_include_all_molecules(parse_state->vol->viz_blocks, (yyvsp[0].ival))); }
#line 6528 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 643:
#line 2524 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.symlist) = mdl_existing_molecule_list(parse_state, (yyvsp[0].str))); }
#line 6534 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 644:
#line 2525 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.symlist) = mdl_existing_molecules_wildcard(parse_state, (yyvsp[0].str))); }
#line 6540 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 645:
#line 2529 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_new_viz_all_times(parse_state, & (yyval.nlist))); }
#line 6546 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 647:
#line 2535 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.frame_list) = (yyvsp[-1].frame_list); }
#line 6552 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 649:
#line 2541 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                        if ((yyvsp[-1].frame_list).frame_head != NULL)
                                                        {
                                                          (yyval.frame_list) = (yyvsp[-1].frame_list);
                                                          if ((yyvsp[0].frame_list).frame_head != NULL)
                                                          {
                                                            (yyval.frame_list).frame_tail->next = (yyvsp[0].frame_list).frame_head;
                                                            (yyval.frame_list).frame_tail = (yyvsp[0].frame_list).frame_tail;
                                                          }
                                                        }
                                                        else if ((yyvsp[0].frame_list).frame_head != NULL)
                                                          (yyval.frame_list) = (yyvsp[0].frame_list);
                                                      }
#line 6570 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 650:
#line 2558 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_new_viz_mol_frames(parse_state, parse_state->vol->viz_blocks, & (yyval.frame_list), OUTPUT_BY_TIME_LIST, (yyvsp[-2].tok), & (yyvsp[0].nlist))); }
#line 6576 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 651:
#line 2562 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_new_viz_all_iterations(parse_state, & (yyval.nlist))); }
#line 6582 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 653:
#line 2569 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.frame_list) = (yyvsp[-1].frame_list); }
#line 6588 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 655:
#line 2575 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                        if ((yyvsp[-1].frame_list).frame_head != NULL)
                                                        {
                                                          (yyval.frame_list) = (yyvsp[-1].frame_list);
                                                          if ((yyvsp[0].frame_list).frame_head != NULL)
                                                          {
                                                            (yyval.frame_list).frame_tail->next = (yyvsp[0].frame_list).frame_head;
                                                            (yyval.frame_list).frame_tail = (yyvsp[0].frame_list).frame_tail;
                                                          }
                                                        }
                                                        else if ((yyvsp[0].frame_list).frame_head != NULL)
                                                          (yyval.frame_list) = (yyvsp[0].frame_list);
                                                      }
#line 6606 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 656:
#line 2592 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_new_viz_mol_frames(parse_state, parse_state->vol->viz_blocks, & (yyval.frame_list), OUTPUT_BY_ITERATION_LIST, (yyvsp[-2].tok), & (yyvsp[0].nlist))); }
#line 6612 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 657:
#line 2595 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = ALL_MOL_DATA; }
#line 6618 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 658:
#line 2596 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = MOL_POS; }
#line 6624 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 659:
#line 2597 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = MOL_ORIENT; }
#line 6630 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 660:
#line 2603 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.frame_list) = (yyvsp[-1].frame_list); }
#line 6636 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 662:
#line 2609 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                        (yyval.frame_list) = (yyvsp[-1].frame_list);
                                                        if ((yyval.frame_list).frame_tail)
                                                        {
                                                          (yyval.frame_list).frame_tail->next = (yyvsp[0].frame_list).frame_head;
                                                          if ((yyvsp[0].frame_list).frame_tail)
                                                            (yyval.frame_list).frame_tail = (yyvsp[0].frame_list).frame_tail;
                                                        }
                                                        else
                                                          (yyval.frame_list) = (yyvsp[0].frame_list);
                                                      }
#line 6652 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 663:
#line 2623 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.frame_list).frame_head = (yyval.frame_list).frame_tail = NULL; }
#line 6658 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 669:
#line 2640 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_region_viz_state(parse_state, parse_state->vol->viz_blocks, (struct region *) (yyvsp[-1].sym)->value, (int) (yyvsp[0].ival))); }
#line 6664 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 670:
#line 2641 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_viz_include_meshes(parse_state, parse_state->vol->viz_blocks, (yyvsp[-1].symlist), (yyvsp[0].ival))); }
#line 6670 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 671:
#line 2642 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_viz_include_all_meshes(parse_state->vol->viz_blocks, (yyvsp[0].ival))); }
#line 6676 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 672:
#line 2648 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.frame_list) = (yyvsp[-1].frame_list); }
#line 6682 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 674:
#line 2654 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                        if ((yyvsp[-1].frame_list).frame_head != NULL)
                                                        {
                                                          (yyval.frame_list) = (yyvsp[-1].frame_list);
                                                          if ((yyvsp[0].frame_list).frame_head != NULL)
                                                          {
                                                            (yyval.frame_list).frame_tail->next = (yyvsp[0].frame_list).frame_head;
                                                            (yyval.frame_list).frame_tail = (yyvsp[0].frame_list).frame_tail;
                                                          }
                                                        }
                                                        else if ((yyvsp[0].frame_list).frame_head != NULL)
                                                          (yyval.frame_list) = (yyvsp[0].frame_list);
                                                      }
#line 6700 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 675:
#line 2671 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_new_viz_mesh_frames(parse_state, parse_state->vol->viz_blocks, & (yyval.frame_list), OUTPUT_BY_TIME_LIST, (yyvsp[-2].tok), & (yyvsp[0].nlist))); }
#line 6706 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 676:
#line 2677 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.frame_list) = (yyvsp[-1].frame_list); }
#line 6712 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 678:
#line 2683 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                        if ((yyvsp[-1].frame_list).frame_head != NULL)
                                                        {
                                                          (yyval.frame_list) = (yyvsp[-1].frame_list);
                                                          if ((yyvsp[0].frame_list).frame_head != NULL)
                                                          {
                                                            (yyval.frame_list).frame_tail->next = (yyvsp[0].frame_list).frame_head;
                                                            (yyval.frame_list).frame_tail = (yyvsp[0].frame_list).frame_tail;
                                                          }
                                                        }
                                                        else if ((yyvsp[0].frame_list).frame_head != NULL)
                                                          (yyval.frame_list) = (yyvsp[0].frame_list);
                                                      }
#line 6730 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 679:
#line 2700 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_new_viz_mesh_frames(parse_state, parse_state->vol->viz_blocks, & (yyval.frame_list), OUTPUT_BY_ITERATION_LIST, (yyvsp[-2].tok), & (yyvsp[0].nlist))); }
#line 6736 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 680:
#line 2703 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = ALL_MESH_DATA; }
#line 6742 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 681:
#line 2704 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = MESH_GEOMETRY; }
#line 6748 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 682:
#line 2705 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = REG_DATA; }
#line 6754 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 683:
#line 2719 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                          struct volume_output_item *vo;
                                                          CHECKN(vo = mdl_new_volume_output_item(parse_state, (yyvsp[-6].str), & (yyvsp[-5].species_lst), (yyvsp[-4].vec3), (yyvsp[-3].vec3), (yyvsp[-2].vec3), (yyvsp[-1].otimes)));
                                                          vo->next = parse_state->vol->volume_output_head;
                                                          parse_state->vol->volume_output_head = vo;
                                                      }
#line 6765 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 684:
#line 2728 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.str) = (yyvsp[0].str); }
#line 6771 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 686:
#line 2734 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                          (yyval.species_lst) = (yyvsp[-1].species_lst);
                                                          (yyval.species_lst).species_count += (yyvsp[0].species_lst).species_count;
                                                          (yyval.species_lst).species_tail->next = (yyvsp[0].species_lst).species_head;
                                                          (yyval.species_lst).species_tail = (yyvsp[0].species_lst).species_tail;
                                                      }
#line 6782 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 687:
#line 2743 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.species_lst) = (yyvsp[0].species_lst); }
#line 6788 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 688:
#line 2746 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                          struct sym_table *sp;
                                                          struct species_list_item *ptrl;
                                                          CHECKN(sp = mdl_existing_molecule(parse_state, (yyvsp[0].str)));

                                                          ptrl = (struct species_list_item *) mem_get(parse_state->species_list_mem);
                                                          if (ptrl == NULL)
                                                          {
                                                            mdlerror_fmt(parse_state, "Out of memory while parsing molecule list");
                                                            return 1;
                                                          }
                                                          ptrl->spec = (struct species *) sp->value;
                                                          ptrl->next = NULL;
                                                          (yyval.species_lst_item) = ptrl;
                                                      }
#line 6808 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 689:
#line 2764 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.species_lst).species_tail = (yyval.species_lst).species_head = (yyvsp[0].species_lst_item); (yyval.species_lst).species_count = 1; }
#line 6814 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 690:
#line 2766 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                        (yyval.species_lst) = (yyvsp[-2].species_lst);
                                                        (yyval.species_lst).species_tail = (yyval.species_lst).species_tail->next = (yyvsp[0].species_lst_item);
                                                        ++ (yyval.species_lst).species_count;
                                                      }
#line 6824 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 691:
#line 2774 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.vec3) = (yyvsp[0].vec3); }
#line 6830 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 692:
#line 2778 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.vec3) = (yyvsp[0].vec3); }
#line 6836 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 693:
#line 2782 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                          if ((yyvsp[0].vec3)->x < 1.0)
                                                          {
                                                            mdl_warning(parse_state, "Voxel count (x dimension) too small.  Setting x count to 1.");
                                                            (yyvsp[0].vec3)->x = 1.0;
                                                          }
                                                          if ((yyvsp[0].vec3)->y < 1.0)
                                                          {
                                                            mdl_warning(parse_state, "Voxel count (y dimension) too small.  Setting y count to 1.");
                                                            (yyvsp[0].vec3)->y = 1.0;
                                                          }
                                                          if ((yyvsp[0].vec3)->z < 1.0)
                                                          {
                                                            mdl_warning(parse_state, "Voxel count (z dimension) too small.  Setting z count to 1.");
                                                            (yyvsp[0].vec3)->z = 1.0;
                                                          }
                                                          (yyval.vec3) = (yyvsp[0].vec3);
                                                      }
#line 6859 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 694:
#line 2803 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.otimes) = mdl_new_output_times_default(parse_state)); }
#line 6865 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 695:
#line 2804 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.otimes) = mdl_new_output_times_step(parse_state, (yyvsp[0].dbl))); }
#line 6871 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 696:
#line 2805 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.otimes) = mdl_new_output_times_iterations(parse_state, & (yyvsp[0].nlist))); }
#line 6877 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 697:
#line 2806 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.otimes) = mdl_new_output_times_time(parse_state, & (yyvsp[0].nlist))); }
#line 6883 "mdlparse.c" /* yacc.c:1646  */
    break;


#line 6887 "mdlparse.c" /* yacc.c:1646  */
      default: break;
    }
  /* User semantic actions sometimes alter yychar, and that requires
     that yytoken be updated with the new translation.  We take the
     approach of translating immediately before every use of yytoken.
     One alternative is translating here after every semantic action,
     but that translation would be missed if the semantic action invokes
     YYABORT, YYACCEPT, or YYERROR immediately after altering yychar or
     if it invokes YYBACKUP.  In the case of YYABORT or YYACCEPT, an
     incorrect destructor might then be invoked immediately.  In the
     case of YYERROR or YYBACKUP, subsequent parser actions might lead
     to an incorrect destructor call or verbose syntax error message
     before the lookahead is translated.  */
  YY_SYMBOL_PRINT ("-> $$ =", yyr1[yyn], &yyval, &yyloc);

  YYPOPSTACK (yylen);
  yylen = 0;
  YY_STACK_PRINT (yyss, yyssp);

  *++yyvsp = yyval;

  /* Now 'shift' the result of the reduction.  Determine what state
     that goes to, based on the state we popped back to and the rule
     number reduced by.  */

  yyn = yyr1[yyn];

  yystate = yypgoto[yyn - YYNTOKENS] + *yyssp;
  if (0 <= yystate && yystate <= YYLAST && yycheck[yystate] == *yyssp)
    yystate = yytable[yystate];
  else
    yystate = yydefgoto[yyn - YYNTOKENS];

  goto yynewstate;


/*--------------------------------------.
| yyerrlab -- here on detecting error.  |
`--------------------------------------*/
yyerrlab:
  /* Make sure we have latest lookahead translation.  See comments at
     user semantic actions for why this is necessary.  */
  yytoken = yychar == YYEMPTY ? YYEMPTY : YYTRANSLATE (yychar);

  /* If not already recovering from an error, report this error.  */
  if (!yyerrstatus)
    {
      ++yynerrs;
#if ! YYERROR_VERBOSE
      yyerror (parse_state, scanner, YY_("syntax error"));
#else
# define YYSYNTAX_ERROR yysyntax_error (&yymsg_alloc, &yymsg, \
                                        yyssp, yytoken)
      {
        char const *yymsgp = YY_("syntax error");
        int yysyntax_error_status;
        yysyntax_error_status = YYSYNTAX_ERROR;
        if (yysyntax_error_status == 0)
          yymsgp = yymsg;
        else if (yysyntax_error_status == 1)
          {
            if (yymsg != yymsgbuf)
              YYSTACK_FREE (yymsg);
            yymsg = (char *) YYSTACK_ALLOC (yymsg_alloc);
            if (!yymsg)
              {
                yymsg = yymsgbuf;
                yymsg_alloc = sizeof yymsgbuf;
                yysyntax_error_status = 2;
              }
            else
              {
                yysyntax_error_status = YYSYNTAX_ERROR;
                yymsgp = yymsg;
              }
          }
        yyerror (parse_state, scanner, yymsgp);
        if (yysyntax_error_status == 2)
          goto yyexhaustedlab;
      }
# undef YYSYNTAX_ERROR
#endif
    }



  if (yyerrstatus == 3)
    {
      /* If just tried and failed to reuse lookahead token after an
         error, discard it.  */

      if (yychar <= YYEOF)
        {
          /* Return failure if at end of input.  */
          if (yychar == YYEOF)
            YYABORT;
        }
      else
        {
          yydestruct ("Error: discarding",
                      yytoken, &yylval, parse_state, scanner);
          yychar = YYEMPTY;
        }
    }

  /* Else will try to reuse lookahead token after shifting the error
     token.  */
  goto yyerrlab1;


/*---------------------------------------------------.
| yyerrorlab -- error raised explicitly by YYERROR.  |
`---------------------------------------------------*/
yyerrorlab:

  /* Pacify compilers like GCC when the user code never invokes
     YYERROR and the label yyerrorlab therefore never appears in user
     code.  */
  if (/*CONSTCOND*/ 0)
     goto yyerrorlab;

  /* Do not reclaim the symbols of the rule whose action triggered
     this YYERROR.  */
  YYPOPSTACK (yylen);
  yylen = 0;
  YY_STACK_PRINT (yyss, yyssp);
  yystate = *yyssp;
  goto yyerrlab1;


/*-------------------------------------------------------------.
| yyerrlab1 -- common code for both syntax error and YYERROR.  |
`-------------------------------------------------------------*/
yyerrlab1:
  yyerrstatus = 3;      /* Each real token shifted decrements this.  */

  for (;;)
    {
      yyn = yypact[yystate];
      if (!yypact_value_is_default (yyn))
        {
          yyn += YYTERROR;
          if (0 <= yyn && yyn <= YYLAST && yycheck[yyn] == YYTERROR)
            {
              yyn = yytable[yyn];
              if (0 < yyn)
                break;
            }
        }

      /* Pop the current state because it cannot handle the error token.  */
      if (yyssp == yyss)
        YYABORT;


      yydestruct ("Error: popping",
                  yystos[yystate], yyvsp, parse_state, scanner);
      YYPOPSTACK (1);
      yystate = *yyssp;
      YY_STACK_PRINT (yyss, yyssp);
    }

  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  *++yyvsp = yylval;
  YY_IGNORE_MAYBE_UNINITIALIZED_END


  /* Shift the error token.  */
  YY_SYMBOL_PRINT ("Shifting", yystos[yyn], yyvsp, yylsp);

  yystate = yyn;
  goto yynewstate;


/*-------------------------------------.
| yyacceptlab -- YYACCEPT comes here.  |
`-------------------------------------*/
yyacceptlab:
  yyresult = 0;
  goto yyreturn;

/*-----------------------------------.
| yyabortlab -- YYABORT comes here.  |
`-----------------------------------*/
yyabortlab:
  yyresult = 1;
  goto yyreturn;

#if !defined yyoverflow || YYERROR_VERBOSE
/*-------------------------------------------------.
| yyexhaustedlab -- memory exhaustion comes here.  |
`-------------------------------------------------*/
yyexhaustedlab:
  yyerror (parse_state, scanner, YY_("memory exhausted"));
  yyresult = 2;
  /* Fall through.  */
#endif

yyreturn:
  if (yychar != YYEMPTY)
    {
      /* Make sure we have latest lookahead translation.  See comments at
         user semantic actions for why this is necessary.  */
      yytoken = YYTRANSLATE (yychar);
      yydestruct ("Cleanup: discarding lookahead",
                  yytoken, &yylval, parse_state, scanner);
    }
  /* Do not reclaim the symbols of the rule whose action triggered
     this YYABORT or YYACCEPT.  */
  YYPOPSTACK (yylen);
  YY_STACK_PRINT (yyss, yyssp);
  while (yyssp != yyss)
    {
      yydestruct ("Cleanup: popping",
                  yystos[*yyssp], yyvsp, parse_state, scanner);
      YYPOPSTACK (1);
    }
#ifndef yyoverflow
  if (yyss != yyssa)
    YYSTACK_FREE (yyss);
#endif
#if YYERROR_VERBOSE
  if (yymsg != yymsgbuf)
    YYSTACK_FREE (yymsg);
#endif
  return yyresult;
}
#line 2809 "../src/../src/mdlparse.y" /* yacc.c:1906  */






/* Begin Bison Epilogue: */

/* mdlerror: Standard error callback from parser.
 *
 *   parse_state: the parser state variables
 *   str:  the error message to display
 */
void mdlerror(struct mdlparse_vars *parse_state, char const *str)
{
  mdlerror_fmt(parse_state, "%s", str);
}

/* mdlerror_fmt: Print a formatted error message regarding an error in the MDL
 *               file.
 *
 *   parse_state: the parser state variables
 *   fmt:  the printf-style format string
 */
void mdlerror_fmt(struct mdlparse_vars *parse_state, char const *fmt, ...)
{
  va_list arglist;
  if (parse_state->vol->procnum != 0)
    return;

  /* print error location */
  if (parse_state->include_stack_ptr == 0)
    mcell_error_raw("Fatal error: After parsing file %s\n",
                    parse_state->vol->curr_file);
  else
    mcell_error_raw("Fatal error: On line: %d of file %s\n",
                    parse_state->line_num[parse_state->include_stack_ptr - 1],
                    parse_state->vol->curr_file);

  /* format error message */
  va_start(arglist, fmt);
  mcell_errorv_raw(fmt, arglist);
  va_end(arglist);

  /* terminate error message and flush */
  mcell_error_raw("\n");
  mcell_die();
}

/* mdlerror_file: Open and parse an MDL file.
 *
 *   parse_state: the parser state variables
 *   name: the path to the MDL file
 */
static int mdlparse_file(struct mdlparse_vars *parse_state, char const *name)
{
  int failure;
  int cur_stack = parse_state->include_stack_ptr ++;
  FILE *infile;
  yyscan_t scanner;
  char const *prev_file;

  /* Put filename and line number on stack */
  if (cur_stack >= MAX_INCLUDE_DEPTH)
  {
    -- parse_state->include_stack_ptr;
    mdlerror_fmt(parse_state, "Includes nested too deeply at file %s\n  included from %s:%d",
                 name,
                 parse_state->include_filename[cur_stack-1],
                 parse_state->line_num[cur_stack-1]);
    return 1;
  }
  parse_state->line_num[cur_stack] = 1;
  parse_state->include_filename[cur_stack] = name;

  /* Open file, or know the reason why */
  no_printf("Opening file %s\n", name);
  if ((infile = fopen(name,"r")) == NULL)
  {
    char *err = mcell_strerror(errno);
    -- parse_state->include_stack_ptr;
    if (cur_stack > 0)
      mdlerror_fmt(parse_state, "Couldn't open file %s\n  included from %s:%d: %s",
                   name,
                   parse_state->include_filename[cur_stack-1],
                   parse_state->line_num[cur_stack-1],
                   err);
    else
      mdlerror_fmt(parse_state, "Couldn't open file %s: %s", name, err);
    return 1;
  }

  /* Create and initialize a lexer */
  if (mdllex_init(&scanner))
  {
    int err = errno;
    if (err == ENOMEM)
      mdlerror_fmt(parse_state, "Couldn't initialize lexer for file %s\n  included from %s:%d: out of memory",
                   name, parse_state->include_filename[cur_stack-1], parse_state->line_num[cur_stack-1]);
    else if (err == EINVAL)
      mdlerror_fmt(parse_state, "Couldn't initialize lexer for file %s\n  included from %s:%d: internal error (invalid argument)",
                   name,
                   parse_state->include_filename[cur_stack-1],
                   parse_state->line_num[cur_stack-1]);
    else
      mdlerror_fmt(parse_state, "Couldn't initialize lexer for file %s\n  included from %s:%d: internal error",
                   name,
                   parse_state->include_filename[cur_stack-1],
                   parse_state->line_num[cur_stack-1]);
    fclose(infile);
    -- parse_state->include_stack_ptr;
    return 1;
  }
  mdlrestart(infile, scanner);

  /* Parse this file */
  prev_file = parse_state->vol->curr_file;
  parse_state->vol->curr_file = name;
  failure = mdlparse(parse_state, scanner);
  parse_state->vol->curr_file = prev_file;
  -- parse_state->include_stack_ptr;

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
  if ((mpv.path_mem = create_mem(sizeof(struct pathway), 4096)) == NULL)
    mcell_allocfailed("Failed to allocate temporary memory pool for reaction pathways.");
  if ((mpv.prod_mem = create_mem(sizeof(struct product), 4096)) == NULL)
    mcell_allocfailed("Failed to allocate temporary memory pool for reaction products.");
  if ((mpv.sym_list_mem = create_mem(sizeof(struct sym_table_list),4096)) == NULL)
    mcell_allocfailed("Failed to allocate temporary memory pool for symbol lists.");
  if ((mpv.species_list_mem = create_mem(sizeof(struct species_list_item), 1024)) == NULL)
    mcell_allocfailed("Failed to allocate temporary memory pool for species lists.");
  if ((mpv.mol_data_list_mem = create_mem(sizeof(struct mcell_species), 1024)) == NULL)
    mcell_allocfailed("Failed to allocate temporary memory pool for oriented species lists.");
  if ((mpv.output_times_mem = create_mem(sizeof(struct output_times), 1024)) == NULL)
    mcell_allocfailed("Failed to allocate temporary memory pool for output times.");

  /* Start parsing at the top-level file */
  vol->curr_file = vol->mdl_infile_name;
  failure = mdlparse_file(&mpv, vol->mdl_infile_name);

  /* Close any open file streams */
  for (int i=0; i<vol->fstream_sym_table->n_bins; ++ i)
  {
    if (vol->fstream_sym_table->entries[i] != NULL)
    {
      for (struct sym_table *symp = vol->fstream_sym_table->entries[i];
           symp != NULL;
           symp = symp->next)
      {
        if (((struct file_stream *) symp->value)->stream == NULL)
          continue;
        mdl_fclose(&mpv, symp);
      }
    }
  }

  /* Check for required settings */
  if (! failure)
  {
    if (vol->time_unit == 0.0)
    {
      mdlerror(&mpv, "A valid model requires a time step to be specified using the TIME_STEP declaration");
      failure = 1;
    }
  }

  /* If we succeeded, prepare the reactions */
  if (failure)
  {
    mdlerror(&mpv, "Failed to parse input file");
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
