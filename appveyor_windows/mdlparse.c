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

  #include "mcell_misc.h"
  #include "mcell_structs.h"
  #include "mcell_viz.h"
  #include "mcell_release.h"
  #include "mcell_objects.h"
  #include "mcell_dyngeom.h"

  /* make sure to declare yyscan_t before including mdlparse.h */
  typedef void *yyscan_t;
  #include "mdlparse.h"

  int mdllex_init(yyscan_t *ptr_yy_globals) ;
  int mdllex_destroy(yyscan_t yyscanner);
  void mdlrestart(FILE *infile, yyscan_t scanner);
  int mdllex(YYSTYPE *yylval, struct mdlparse_vars *parse_state, yyscan_t scanner);


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

#line 137 "mdlparse.c" /* yacc.c:339  */

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
    FALSE = 329,
    FCLOSE = 330,
    FILENAME = 331,
    FILENAME_PREFIX = 332,
    FILE_OUTPUT_REPORT = 333,
    FINAL_SUMMARY = 334,
    FLOOR = 335,
    FOPEN = 336,
    FORMAT = 337,
    FPRINTF = 338,
    FPRINT_TIME = 339,
    FRONT = 340,
    FRONT_CROSSINGS = 341,
    FRONT_HITS = 342,
    GAUSSIAN_RELEASE_NUMBER = 343,
    HEADER = 344,
    HIGH_PROBABILITY_THRESHOLD = 345,
    HIGH_REACTION_PROBABILITY = 346,
    IGNORED = 347,
    INCLUDE_ELEMENTS = 348,
    INCLUDE_FILE = 349,
    INCLUDE_PATCH = 350,
    INCLUDE_REGION = 351,
    INPUT_FILE = 352,
    INSTANTIATE = 353,
    LLINTEGER = 354,
    FULLY_RANDOM = 355,
    INTERACTION_RADIUS = 356,
    ITERATION_LIST = 357,
    ITERATION_NUMBERS = 358,
    ITERATION_REPORT = 359,
    ITERATIONS = 360,
    KEEP_CHECKPOINT_FILES = 361,
    LEFT = 362,
    LIFETIME_THRESHOLD = 363,
    LIFETIME_TOO_SHORT = 364,
    LIST = 365,
    LOCATION = 366,
    LOG = 367,
    LOG10 = 368,
    MAX_TOK = 369,
    MAXIMUM_STEP_LENGTH = 370,
    MEAN_DIAMETER = 371,
    MEAN_NUMBER = 372,
    MEMORY_PARTITION_X = 373,
    MEMORY_PARTITION_Y = 374,
    MEMORY_PARTITION_Z = 375,
    MEMORY_PARTITION_POOL = 376,
    MICROSCOPIC_REVERSIBILITY = 377,
    MIN_TOK = 378,
    MISSED_REACTIONS = 379,
    MISSED_REACTION_THRESHOLD = 380,
    MISSING_SURFACE_ORIENTATION = 381,
    MOD = 382,
    MODE = 383,
    MODIFY_SURFACE_REGIONS = 384,
    MOLECULE = 385,
    MOLECULE_COLLISION_REPORT = 386,
    MOLECULE_DENSITY = 387,
    MOLECULE_NUMBER = 388,
    MOLECULE_POSITIONS = 389,
    MOLECULES = 390,
    MOLECULE_PLACEMENT_FAILURE = 391,
    NAME_LIST = 392,
    NEAREST_POINT = 393,
    NEAREST_TRIANGLE = 394,
    NEGATIVE_DIFFUSION_CONSTANT = 395,
    NEGATIVE_REACTION_RATE = 396,
    NO = 397,
    NOEXIT = 398,
    NONE = 399,
    NO_SPECIES = 400,
    NOT_EQUAL = 401,
    NOTIFICATIONS = 402,
    NUMBER_OF_SUBUNITS = 403,
    NUMBER_OF_TRAINS = 404,
    NUMBER_TO_RELEASE = 405,
    OBJECT = 406,
    OFF = 407,
    ON = 408,
    ORIENTATIONS = 409,
    OUTPUT_BUFFER_SIZE = 410,
    INVALID_OUTPUT_STEP_TIME = 411,
    LARGE_MOLECULAR_DISPLACEMENT = 412,
    ADD_REMOVE_MESH = 413,
    OVERWRITTEN_OUTPUT_FILE = 414,
    PARTITION_LOCATION_REPORT = 415,
    PARTITION_X = 416,
    PARTITION_Y = 417,
    PARTITION_Z = 418,
    PERIODIC_BOX = 419,
    PERIODIC_X = 420,
    PERIODIC_Y = 421,
    PERIODIC_Z = 422,
    PERIODIC_TRADITIONAL = 423,
    PI_TOK = 424,
    POLYGON_LIST = 425,
    POSITIONS = 426,
    PRINTF = 427,
    PRINT_TIME = 428,
    PROBABILITY_REPORT = 429,
    PROBABILITY_REPORT_THRESHOLD = 430,
    PROGRESS_REPORT = 431,
    RADIAL_DIRECTIONS = 432,
    RADIAL_SUBDIVISIONS = 433,
    RAND_GAUSSIAN = 434,
    RAND_UNIFORM = 435,
    REACTION_DATA_OUTPUT = 436,
    REACTION_OUTPUT_REPORT = 437,
    REAL = 438,
    RECTANGULAR_RELEASE_SITE = 439,
    RECTANGULAR_TOKEN = 440,
    REFLECTIVE = 441,
    RELEASE_EVENT_REPORT = 442,
    RELEASE_INTERVAL = 443,
    RELEASE_PATTERN = 444,
    RELEASE_PROBABILITY = 445,
    RELEASE_SITE = 446,
    REMOVE_ELEMENTS = 447,
    RIGHT = 448,
    ROTATE = 449,
    ROUND_OFF = 450,
    SCALE = 451,
    SEED = 452,
    SHAPE = 453,
    SHOW_EXACT_TIME = 454,
    SIN = 455,
    SITE_DIAMETER = 456,
    SITE_RADIUS = 457,
    SPACE_STEP = 458,
    SPHERICAL = 459,
    SPHERICAL_RELEASE_SITE = 460,
    SPHERICAL_SHELL = 461,
    SPHERICAL_SHELL_SITE = 462,
    SPRINTF = 463,
    SQRT = 464,
    STANDARD_DEVIATION = 465,
    PERIODIC_BOX_INITIAL = 466,
    STEP = 467,
    STRING_TO_NUM = 468,
    STR_VALUE = 469,
    SUBUNIT = 470,
    SUBUNIT_RELATIONSHIPS = 471,
    SUMMATION_OPERATOR = 472,
    SURFACE_CLASS = 473,
    SURFACE_ONLY = 474,
    TAN = 475,
    TARGET_ONLY = 476,
    TET_ELEMENT_CONNECTIONS = 477,
    THROUGHPUT_REPORT = 478,
    TIME_LIST = 479,
    TIME_POINTS = 480,
    TIME_STEP = 481,
    TIME_STEP_MAX = 482,
    TO = 483,
    TOP = 484,
    TRAIN_DURATION = 485,
    TRAIN_INTERVAL = 486,
    TRANSLATE = 487,
    TRANSPARENT = 488,
    TRIGGER = 489,
    TRUE = 490,
    UNLIMITED = 491,
    USELESS_VOLUME_ORIENTATION = 492,
    VACANCY_SEARCH_DISTANCE = 493,
    VAR = 494,
    VARYING_PROBABILITY_REPORT = 495,
    VERTEX_LIST = 496,
    VIZ_OUTPUT = 497,
    VIZ_OUTPUT_REPORT = 498,
    VIZ_VALUE = 499,
    VOLUME_DATA_OUTPUT = 500,
    VOLUME_OUTPUT_REPORT = 501,
    VOLUME_DEPENDENT_RELEASE_NUMBER = 502,
    VOLUME_ONLY = 503,
    VOXEL_COUNT = 504,
    VOXEL_LIST = 505,
    VOXEL_SIZE = 506,
    WARNING = 507,
    WARNINGS = 508,
    WORLD = 509,
    YES = 510,
    UNARYMINUS = 511
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
#define ALL_MOLECULES 268
#define ALL_NOTIFICATIONS 269
#define ALL_TIMES 270
#define ALL_WARNINGS 271
#define ASCII 272
#define ASIN 273
#define ASPECT_RATIO 274
#define ATAN 275
#define BACK 276
#define BACK_CROSSINGS 277
#define BACK_HITS 278
#define BOTTOM 279
#define BOX 280
#define BOX_TRIANGULATION_REPORT 281
#define BRIEF 282
#define CEIL 283
#define CELLBLENDER 284
#define CENTER_MOLECULES_ON_GRID 285
#define CHECKPOINT_INFILE 286
#define CHECKPOINT_ITERATIONS 287
#define CHECKPOINT_OUTFILE 288
#define CHECKPOINT_REALTIME 289
#define CHECKPOINT_REPORT 290
#define CLAMP_CONCENTRATION 291
#define CLOSE_PARTITION_SPACING 292
#define CONCENTRATION 293
#define CORNERS 294
#define COS 295
#define COUNT 296
#define CUBIC 297
#define CUBIC_RELEASE_SITE 298
#define CUSTOM_SPACE_STEP 299
#define CUSTOM_TIME_STEP 300
#define DEFINE_MOLECULE 301
#define DEFINE_MOLECULES 302
#define DEFINE_REACTIONS 303
#define DEFINE_RELEASE_PATTERN 304
#define DEFINE_SURFACE_CLASS 305
#define DEFINE_SURFACE_CLASSES 306
#define DEFINE_SURFACE_REGIONS 307
#define DEGENERATE_POLYGONS 308
#define DELAY 309
#define DENSITY 310
#define DIFFUSION_CONSTANT_2D 311
#define DIFFUSION_CONSTANT_3D 312
#define DIFFUSION_CONSTANT_REPORT 313
#define DYNAMIC_GEOMETRY 314
#define DYNAMIC_GEOMETRY_MOLECULE_PLACEMENT 315
#define EFFECTOR_GRID_DENSITY 316
#define ELEMENT_CONNECTIONS 317
#define ELLIPTIC 318
#define ELLIPTIC_RELEASE_SITE 319
#define EQUAL 320
#define ERROR 321
#define ESTIMATE_CONCENTRATION 322
#define EXCLUDE_ELEMENTS 323
#define EXCLUDE_PATCH 324
#define EXCLUDE_REGION 325
#define EXIT 326
#define EXP 327
#define EXPRESSION 328
#define FALSE 329
#define FCLOSE 330
#define FILENAME 331
#define FILENAME_PREFIX 332
#define FILE_OUTPUT_REPORT 333
#define FINAL_SUMMARY 334
#define FLOOR 335
#define FOPEN 336
#define FORMAT 337
#define FPRINTF 338
#define FPRINT_TIME 339
#define FRONT 340
#define FRONT_CROSSINGS 341
#define FRONT_HITS 342
#define GAUSSIAN_RELEASE_NUMBER 343
#define HEADER 344
#define HIGH_PROBABILITY_THRESHOLD 345
#define HIGH_REACTION_PROBABILITY 346
#define IGNORED 347
#define INCLUDE_ELEMENTS 348
#define INCLUDE_FILE 349
#define INCLUDE_PATCH 350
#define INCLUDE_REGION 351
#define INPUT_FILE 352
#define INSTANTIATE 353
#define LLINTEGER 354
#define FULLY_RANDOM 355
#define INTERACTION_RADIUS 356
#define ITERATION_LIST 357
#define ITERATION_NUMBERS 358
#define ITERATION_REPORT 359
#define ITERATIONS 360
#define KEEP_CHECKPOINT_FILES 361
#define LEFT 362
#define LIFETIME_THRESHOLD 363
#define LIFETIME_TOO_SHORT 364
#define LIST 365
#define LOCATION 366
#define LOG 367
#define LOG10 368
#define MAX_TOK 369
#define MAXIMUM_STEP_LENGTH 370
#define MEAN_DIAMETER 371
#define MEAN_NUMBER 372
#define MEMORY_PARTITION_X 373
#define MEMORY_PARTITION_Y 374
#define MEMORY_PARTITION_Z 375
#define MEMORY_PARTITION_POOL 376
#define MICROSCOPIC_REVERSIBILITY 377
#define MIN_TOK 378
#define MISSED_REACTIONS 379
#define MISSED_REACTION_THRESHOLD 380
#define MISSING_SURFACE_ORIENTATION 381
#define MOD 382
#define MODE 383
#define MODIFY_SURFACE_REGIONS 384
#define MOLECULE 385
#define MOLECULE_COLLISION_REPORT 386
#define MOLECULE_DENSITY 387
#define MOLECULE_NUMBER 388
#define MOLECULE_POSITIONS 389
#define MOLECULES 390
#define MOLECULE_PLACEMENT_FAILURE 391
#define NAME_LIST 392
#define NEAREST_POINT 393
#define NEAREST_TRIANGLE 394
#define NEGATIVE_DIFFUSION_CONSTANT 395
#define NEGATIVE_REACTION_RATE 396
#define NO 397
#define NOEXIT 398
#define NONE 399
#define NO_SPECIES 400
#define NOT_EQUAL 401
#define NOTIFICATIONS 402
#define NUMBER_OF_SUBUNITS 403
#define NUMBER_OF_TRAINS 404
#define NUMBER_TO_RELEASE 405
#define OBJECT 406
#define OFF 407
#define ON 408
#define ORIENTATIONS 409
#define OUTPUT_BUFFER_SIZE 410
#define INVALID_OUTPUT_STEP_TIME 411
#define LARGE_MOLECULAR_DISPLACEMENT 412
#define ADD_REMOVE_MESH 413
#define OVERWRITTEN_OUTPUT_FILE 414
#define PARTITION_LOCATION_REPORT 415
#define PARTITION_X 416
#define PARTITION_Y 417
#define PARTITION_Z 418
#define PERIODIC_BOX 419
#define PERIODIC_X 420
#define PERIODIC_Y 421
#define PERIODIC_Z 422
#define PERIODIC_TRADITIONAL 423
#define PI_TOK 424
#define POLYGON_LIST 425
#define POSITIONS 426
#define PRINTF 427
#define PRINT_TIME 428
#define PROBABILITY_REPORT 429
#define PROBABILITY_REPORT_THRESHOLD 430
#define PROGRESS_REPORT 431
#define RADIAL_DIRECTIONS 432
#define RADIAL_SUBDIVISIONS 433
#define RAND_GAUSSIAN 434
#define RAND_UNIFORM 435
#define REACTION_DATA_OUTPUT 436
#define REACTION_OUTPUT_REPORT 437
#define REAL 438
#define RECTANGULAR_RELEASE_SITE 439
#define RECTANGULAR_TOKEN 440
#define REFLECTIVE 441
#define RELEASE_EVENT_REPORT 442
#define RELEASE_INTERVAL 443
#define RELEASE_PATTERN 444
#define RELEASE_PROBABILITY 445
#define RELEASE_SITE 446
#define REMOVE_ELEMENTS 447
#define RIGHT 448
#define ROTATE 449
#define ROUND_OFF 450
#define SCALE 451
#define SEED 452
#define SHAPE 453
#define SHOW_EXACT_TIME 454
#define SIN 455
#define SITE_DIAMETER 456
#define SITE_RADIUS 457
#define SPACE_STEP 458
#define SPHERICAL 459
#define SPHERICAL_RELEASE_SITE 460
#define SPHERICAL_SHELL 461
#define SPHERICAL_SHELL_SITE 462
#define SPRINTF 463
#define SQRT 464
#define STANDARD_DEVIATION 465
#define PERIODIC_BOX_INITIAL 466
#define STEP 467
#define STRING_TO_NUM 468
#define STR_VALUE 469
#define SUBUNIT 470
#define SUBUNIT_RELATIONSHIPS 471
#define SUMMATION_OPERATOR 472
#define SURFACE_CLASS 473
#define SURFACE_ONLY 474
#define TAN 475
#define TARGET_ONLY 476
#define TET_ELEMENT_CONNECTIONS 477
#define THROUGHPUT_REPORT 478
#define TIME_LIST 479
#define TIME_POINTS 480
#define TIME_STEP 481
#define TIME_STEP_MAX 482
#define TO 483
#define TOP 484
#define TRAIN_DURATION 485
#define TRAIN_INTERVAL 486
#define TRANSLATE 487
#define TRANSPARENT 488
#define TRIGGER 489
#define TRUE 490
#define UNLIMITED 491
#define USELESS_VOLUME_ORIENTATION 492
#define VACANCY_SEARCH_DISTANCE 493
#define VAR 494
#define VARYING_PROBABILITY_REPORT 495
#define VERTEX_LIST 496
#define VIZ_OUTPUT 497
#define VIZ_OUTPUT_REPORT 498
#define VIZ_VALUE 499
#define VOLUME_DATA_OUTPUT 500
#define VOLUME_OUTPUT_REPORT 501
#define VOLUME_DEPENDENT_RELEASE_NUMBER 502
#define VOLUME_ONLY 503
#define VOXEL_COUNT 504
#define VOXEL_LIST 505
#define VOXEL_SIZE 506
#define WARNING 507
#define WARNINGS 508
#define WORLD 509
#define YES 510
#define UNARYMINUS 511

/* Value type.  */
#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED

union YYSTYPE
{
#line 67 "../src/../src/mdlparse.y" /* yacc.c:355  */

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


#line 755 "mdlparse.c" /* yacc.c:355  */
};

typedef union YYSTYPE YYSTYPE;
# define YYSTYPE_IS_TRIVIAL 1
# define YYSTYPE_IS_DECLARED 1
#endif



int mdlparse (struct mdlparse_vars *parse_state, yyscan_t scanner);

#endif /* !YY_MDL_MDLPARSE_H_INCLUDED  */

/* Copy the second part of user declarations.  */

#line 771 "mdlparse.c" /* yacc.c:358  */

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
#define YYFINAL  150
/* YYLAST -- Last index in YYTABLE.  */
#define YYLAST   2442

/* YYNTOKENS -- Number of terminals.  */
#define YYNTOKENS  277
/* YYNNTS -- Number of nonterminals.  */
#define YYNNTS  289
/* YYNRULES -- Number of rules.  */
#define YYNRULES  620
/* YYNSTATES -- Number of states.  */
#define YYNSTATES  1227

/* YYTRANSLATE[YYX] -- Symbol number corresponding to YYX as returned
   by yylex, with out-of-bounds checking.  */
#define YYUNDEFTOK  2
#define YYMAXUTOK   511

#define YYTRANSLATE(YYX)                                                \
  ((unsigned int) (YYX) <= YYMAXUTOK ? yytranslate[YYX] : YYUNDEFTOK)

/* YYTRANSLATE[TOKEN-NUM] -- Symbol number corresponding to TOKEN-NUM
   as returned by yylex, without out-of-bounds checking.  */
static const yytype_uint16 yytranslate[] =
{
       0,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,   257,   268,
     272,   273,   261,   259,   269,   260,     2,   262,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,   258,   267,
     275,   256,   274,     2,   276,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,   265,     2,   266,   264,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,   270,     2,   271,     2,     2,     2,     2,
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
     255,   263
};

#if YYDEBUG
  /* YYRLINE[YYN] -- Source line where rule number YYN was defined.  */
static const yytype_uint16 yyrline[] =
{
       0,   598,   598,   602,   603,   608,   609,   610,   611,   612,
     613,   614,   615,   616,   617,   618,   619,   620,   621,   622,
     623,   624,   625,   626,   627,   628,   633,   636,   639,   642,
     645,   648,   651,   652,   655,   656,   657,   658,   659,   660,
     663,   664,   665,   666,   670,   671,   672,   682,   694,   697,
     700,   712,   713,   726,   727,   733,   756,   757,   758,   759,
     762,   765,   768,   769,   780,   783,   786,   787,   790,   791,
     794,   795,   798,   799,   802,   806,   807,   808,   809,   810,
     811,   812,   813,   814,   815,   816,   817,   818,   819,   820,
     821,   822,   823,   824,   825,   826,   827,   828,   829,   830,
     831,   832,   833,   834,   835,   839,   840,   844,   845,   846,
     847,   850,   856,   857,   858,   859,   860,   861,   862,   865,
     869,   872,   875,   878,   881,   884,   885,   895,   896,   897,
     909,   913,   919,   924,   928,   937,   941,   942,   946,   947,
     948,   949,   950,   951,   952,   953,   954,   955,   956,   957,
     958,   959,   960,   961,   962,   966,   967,   971,   975,   976,
     983,   987,   988,   992,   993,   994,   995,   996,   997,   998,
     999,  1000,  1001,  1002,  1003,  1004,  1005,  1006,  1007,  1008,
    1009,  1013,  1014,  1015,  1021,  1022,  1023,  1024,  1025,  1029,
    1030,  1031,  1035,  1036,  1037,  1038,  1046,  1047,  1048,  1049,
    1050,  1051,  1052,  1053,  1054,  1055,  1056,  1057,  1058,  1059,
    1060,  1061,  1062,  1063,  1070,  1071,  1072,  1073,  1077,  1081,
    1082,  1083,  1090,  1091,  1094,  1098,  1102,  1103,  1107,  1115,
    1118,  1122,  1123,  1127,  1128,  1137,  1148,  1149,  1153,  1154,
    1164,  1168,  1172,  1185,  1186,  1191,  1195,  1201,  1202,  1207,
    1207,  1212,  1215,  1217,  1222,  1223,  1227,  1230,  1236,  1241,
    1242,  1243,  1246,  1247,  1250,  1254,  1258,  1265,  1269,  1278,
    1282,  1291,  1298,  1303,  1304,  1307,  1310,  1311,  1314,  1315,
    1316,  1319,  1324,  1330,  1331,  1332,  1333,  1336,  1337,  1341,
    1346,  1347,  1350,  1354,  1355,  1359,  1362,  1363,  1366,  1367,
    1371,  1372,  1375,  1386,  1403,  1404,  1405,  1409,  1410,  1411,
    1418,  1425,  1428,  1432,  1439,  1441,  1443,  1445,  1447,  1451,
    1452,  1459,  1459,  1470,  1473,  1474,  1475,  1476,  1477,  1486,
    1489,  1492,  1495,  1497,  1501,  1505,  1506,  1507,  1512,  1525,
    1526,  1529,  1530,  1535,  1534,  1541,  1542,  1547,  1546,  1554,
    1555,  1556,  1557,  1558,  1559,  1560,  1561,  1568,  1569,  1570,
    1571,  1572,  1577,  1576,  1583,  1584,  1585,  1586,  1587,  1591,
    1592,  1595,  1599,  1600,  1601,  1608,  1609,  1610,  1611,  1612,
    1613,  1615,  1620,  1621,  1625,  1626,  1627,  1628,  1633,  1634,
    1640,  1647,  1655,  1656,  1660,  1661,  1666,  1669,  1677,  1674,
    1692,  1695,  1698,  1699,  1703,  1708,  1709,  1713,  1716,  1718,
    1724,  1725,  1729,  1729,  1741,  1742,  1745,  1746,  1747,  1748,
    1749,  1750,  1751,  1755,  1756,  1761,  1762,  1763,  1764,  1768,
    1773,  1777,  1781,  1782,  1785,  1786,  1787,  1790,  1793,  1794,
    1797,  1800,  1801,  1805,  1811,  1812,  1817,  1818,  1817,  1828,
    1825,  1838,  1842,  1846,  1850,  1860,  1863,  1857,  1870,  1871,
    1867,  1880,  1881,  1885,  1886,  1890,  1891,  1895,  1896,  1899,
    1900,  1914,  1920,  1921,  1926,  1927,  1929,  1926,  1938,  1941,
    1943,  1948,  1949,  1953,  1960,  1966,  1967,  1972,  1972,  1982,
    1981,  1992,  1993,  2004,  2005,  2006,  2009,  2013,  2021,  2028,
    2029,  2043,  2044,  2045,  2049,  2049,  2055,  2056,  2057,  2061,
    2065,  2069,  2070,  2079,  2083,  2084,  2085,  2086,  2087,  2088,
    2089,  2090,  2091,  2096,  2096,  2098,  2099,  2099,  2103,  2104,
    2105,  2106,  2107,  2110,  2113,  2117,  2127,  2128,  2129,  2130,
    2131,  2132,  2136,  2141,  2146,  2151,  2156,  2161,  2165,  2166,
    2167,  2170,  2171,  2174,  2175,  2176,  2177,  2178,  2179,  2180,
    2181,  2184,  2185,  2192,  2192,  2199,  2200,  2204,  2205,  2208,
    2209,  2210,  2214,  2215,  2225,  2228,  2232,  2238,  2239,  2254,
    2255,  2256,  2260,  2266,  2267,  2271,  2272,  2276,  2278,  2282,
    2283,  2287,  2288,  2291,  2297,  2298,  2315,  2320,  2321,  2325,
    2331,  2332,  2349,  2353,  2354,  2355,  2362,  2378,  2382,  2383,
    2393,  2396,  2414,  2415,  2424,  2428,  2432,  2453,  2454,  2455,
    2456
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
  "ALL_MOLECULES", "ALL_NOTIFICATIONS", "ALL_TIMES", "ALL_WARNINGS",
  "ASCII", "ASIN", "ASPECT_RATIO", "ATAN", "BACK", "BACK_CROSSINGS",
  "BACK_HITS", "BOTTOM", "BOX", "BOX_TRIANGULATION_REPORT", "BRIEF",
  "CEIL", "CELLBLENDER", "CENTER_MOLECULES_ON_GRID", "CHECKPOINT_INFILE",
  "CHECKPOINT_ITERATIONS", "CHECKPOINT_OUTFILE", "CHECKPOINT_REALTIME",
  "CHECKPOINT_REPORT", "CLAMP_CONCENTRATION", "CLOSE_PARTITION_SPACING",
  "CONCENTRATION", "CORNERS", "COS", "COUNT", "CUBIC",
  "CUBIC_RELEASE_SITE", "CUSTOM_SPACE_STEP", "CUSTOM_TIME_STEP",
  "DEFINE_MOLECULE", "DEFINE_MOLECULES", "DEFINE_REACTIONS",
  "DEFINE_RELEASE_PATTERN", "DEFINE_SURFACE_CLASS",
  "DEFINE_SURFACE_CLASSES", "DEFINE_SURFACE_REGIONS",
  "DEGENERATE_POLYGONS", "DELAY", "DENSITY", "DIFFUSION_CONSTANT_2D",
  "DIFFUSION_CONSTANT_3D", "DIFFUSION_CONSTANT_REPORT", "DYNAMIC_GEOMETRY",
  "DYNAMIC_GEOMETRY_MOLECULE_PLACEMENT", "EFFECTOR_GRID_DENSITY",
  "ELEMENT_CONNECTIONS", "ELLIPTIC", "ELLIPTIC_RELEASE_SITE", "EQUAL",
  "ERROR", "ESTIMATE_CONCENTRATION", "EXCLUDE_ELEMENTS", "EXCLUDE_PATCH",
  "EXCLUDE_REGION", "EXIT", "EXP", "EXPRESSION", "FALSE", "FCLOSE",
  "FILENAME", "FILENAME_PREFIX", "FILE_OUTPUT_REPORT", "FINAL_SUMMARY",
  "FLOOR", "FOPEN", "FORMAT", "FPRINTF", "FPRINT_TIME", "FRONT",
  "FRONT_CROSSINGS", "FRONT_HITS", "GAUSSIAN_RELEASE_NUMBER", "HEADER",
  "HIGH_PROBABILITY_THRESHOLD", "HIGH_REACTION_PROBABILITY", "IGNORED",
  "INCLUDE_ELEMENTS", "INCLUDE_FILE", "INCLUDE_PATCH", "INCLUDE_REGION",
  "INPUT_FILE", "INSTANTIATE", "LLINTEGER", "FULLY_RANDOM",
  "INTERACTION_RADIUS", "ITERATION_LIST", "ITERATION_NUMBERS",
  "ITERATION_REPORT", "ITERATIONS", "KEEP_CHECKPOINT_FILES", "LEFT",
  "LIFETIME_THRESHOLD", "LIFETIME_TOO_SHORT", "LIST", "LOCATION", "LOG",
  "LOG10", "MAX_TOK", "MAXIMUM_STEP_LENGTH", "MEAN_DIAMETER",
  "MEAN_NUMBER", "MEMORY_PARTITION_X", "MEMORY_PARTITION_Y",
  "MEMORY_PARTITION_Z", "MEMORY_PARTITION_POOL",
  "MICROSCOPIC_REVERSIBILITY", "MIN_TOK", "MISSED_REACTIONS",
  "MISSED_REACTION_THRESHOLD", "MISSING_SURFACE_ORIENTATION", "MOD",
  "MODE", "MODIFY_SURFACE_REGIONS", "MOLECULE",
  "MOLECULE_COLLISION_REPORT", "MOLECULE_DENSITY", "MOLECULE_NUMBER",
  "MOLECULE_POSITIONS", "MOLECULES", "MOLECULE_PLACEMENT_FAILURE",
  "NAME_LIST", "NEAREST_POINT", "NEAREST_TRIANGLE",
  "NEGATIVE_DIFFUSION_CONSTANT", "NEGATIVE_REACTION_RATE", "NO", "NOEXIT",
  "NONE", "NO_SPECIES", "NOT_EQUAL", "NOTIFICATIONS", "NUMBER_OF_SUBUNITS",
  "NUMBER_OF_TRAINS", "NUMBER_TO_RELEASE", "OBJECT", "OFF", "ON",
  "ORIENTATIONS", "OUTPUT_BUFFER_SIZE", "INVALID_OUTPUT_STEP_TIME",
  "LARGE_MOLECULAR_DISPLACEMENT", "ADD_REMOVE_MESH",
  "OVERWRITTEN_OUTPUT_FILE", "PARTITION_LOCATION_REPORT", "PARTITION_X",
  "PARTITION_Y", "PARTITION_Z", "PERIODIC_BOX", "PERIODIC_X", "PERIODIC_Y",
  "PERIODIC_Z", "PERIODIC_TRADITIONAL", "PI_TOK", "POLYGON_LIST",
  "POSITIONS", "PRINTF", "PRINT_TIME", "PROBABILITY_REPORT",
  "PROBABILITY_REPORT_THRESHOLD", "PROGRESS_REPORT", "RADIAL_DIRECTIONS",
  "RADIAL_SUBDIVISIONS", "RAND_GAUSSIAN", "RAND_UNIFORM",
  "REACTION_DATA_OUTPUT", "REACTION_OUTPUT_REPORT", "REAL",
  "RECTANGULAR_RELEASE_SITE", "RECTANGULAR_TOKEN", "REFLECTIVE",
  "RELEASE_EVENT_REPORT", "RELEASE_INTERVAL", "RELEASE_PATTERN",
  "RELEASE_PROBABILITY", "RELEASE_SITE", "REMOVE_ELEMENTS", "RIGHT",
  "ROTATE", "ROUND_OFF", "SCALE", "SEED", "SHAPE", "SHOW_EXACT_TIME",
  "SIN", "SITE_DIAMETER", "SITE_RADIUS", "SPACE_STEP", "SPHERICAL",
  "SPHERICAL_RELEASE_SITE", "SPHERICAL_SHELL", "SPHERICAL_SHELL_SITE",
  "SPRINTF", "SQRT", "STANDARD_DEVIATION", "PERIODIC_BOX_INITIAL", "STEP",
  "STRING_TO_NUM", "STR_VALUE", "SUBUNIT", "SUBUNIT_RELATIONSHIPS",
  "SUMMATION_OPERATOR", "SURFACE_CLASS", "SURFACE_ONLY", "TAN",
  "TARGET_ONLY", "TET_ELEMENT_CONNECTIONS", "THROUGHPUT_REPORT",
  "TIME_LIST", "TIME_POINTS", "TIME_STEP", "TIME_STEP_MAX", "TO", "TOP",
  "TRAIN_DURATION", "TRAIN_INTERVAL", "TRANSLATE", "TRANSPARENT",
  "TRIGGER", "TRUE", "UNLIMITED", "USELESS_VOLUME_ORIENTATION",
  "VACANCY_SEARCH_DISTANCE", "VAR", "VARYING_PROBABILITY_REPORT",
  "VERTEX_LIST", "VIZ_OUTPUT", "VIZ_OUTPUT_REPORT", "VIZ_VALUE",
  "VOLUME_DATA_OUTPUT", "VOLUME_OUTPUT_REPORT",
  "VOLUME_DEPENDENT_RELEASE_NUMBER", "VOLUME_ONLY", "VOXEL_COUNT",
  "VOXEL_LIST", "VOXEL_SIZE", "WARNING", "WARNINGS", "WORLD", "YES", "'='",
  "'&'", "':'", "'+'", "'-'", "'*'", "'/'", "UNARYMINUS", "'^'", "'['",
  "']'", "';'", "'\\''", "','", "'{'", "'}'", "'('", "')'", "'>'", "'<'",
  "'@'", "$accept", "mdl_format", "mdl_stmt_list", "mdl_stmt", "str_value",
  "var", "file_name", "existing_object", "existing_region", "point",
  "point_or_num", "boolean", "orientation_class", "list_orient_marks",
  "head_mark", "tail_mark", "orient_class_number", "list_range_specs",
  "range_spec", "include_stmt", "assignment_stmt", "assign_var",
  "existing_var_only", "array_value", "array_expr_only", "existing_array",
  "num_expr", "num_value", "intOrReal", "num_expr_only",
  "existing_num_var", "arith_expr", "str_expr", "str_expr_only",
  "existing_str_var", "io_stmt", "fopen_stmt", "new_file_stream",
  "file_mode", "fclose_stmt", "existing_file_stream", "format_string",
  "list_args", "list_arg", "printf_stmt", "fprintf_stmt", "sprintf_stmt",
  "print_time_stmt", "fprint_time_stmt", "notification_def",
  "notification_list", "notification_item_def", "notify_bilevel",
  "notify_level", "warnings_def", "warning_list", "warning_item_def",
  "warning_level", "chkpt_stmt", "exit_or_no", "time_expr",
  "parameter_def", "memory_partition_def", "partition_def",
  "partition_dimension", "molecules_def", "define_one_molecule",
  "define_multiple_molecules", "list_molecule_stmts", "molecule_stmt",
  "molecule_name", "new_molecule", "diffusion_def", "mol_timestep_def",
  "target_def", "maximum_step_length_def", "existing_molecule",
  "existing_surface_molecule", "existing_molecule_opt_orient",
  "surface_classes_def", "define_one_surface_class",
  "define_multiple_surface_classes", "list_surface_class_stmts",
  "surface_class_stmt", "$@1", "existing_surface_class",
  "list_surface_prop_stmts", "surface_prop_stmt", "surface_rxn_stmt",
  "surface_rxn_type", "equals_or_to", "surface_class_mol_stmt",
  "surface_mol_stmt", "list_surface_mol_density", "list_surface_mol_num",
  "surface_mol_quant", "rx_net_def", "list_rx_stmts", "rx_stmt",
  "list_dashes", "right_arrow", "left_arrow", "double_arrow",
  "right_cat_arrow", "double_cat_arrow", "reaction_arrow",
  "new_rxn_pathname", "rxn", "reactant_list", "reactant",
  "opt_reactant_surface_class", "reactant_surface_class", "product_list",
  "product", "rx_rate_syntax", "rx_rate1", "rx_rate2", "rx_dir_rate",
  "atomic_rate", "release_pattern_def", "new_release_pattern",
  "existing_release_pattern_xor_rxpn", "list_req_release_pattern_cmds",
  "train_count", "instance_def", "$@2", "physical_object_def",
  "object_def", "new_object", "start_object", "end_object",
  "list_opt_object_cmds", "opt_object_cmd", "transformation",
  "meta_object_def", "list_objects", "object_ref", "existing_object_ref",
  "$@3", "release_site_def", "release_site_def_new", "$@4",
  "release_site_geom", "release_region_expr", "release_site_def_old",
  "$@5", "release_site_geom_old", "list_release_site_cmds",
  "existing_num_or_array", "release_site_cmd", "site_size_cmd",
  "release_number_cmd", "constant_release_number_cmd",
  "gaussian_release_number_cmd", "volume_dependent_number_cmd",
  "concentration_dependent_release_cmd", "molecule_release_pos_list",
  "molecule_release_pos", "new_object_name", "polygon_list_def", "@6",
  "vertex_list_cmd", "single_vertex", "list_points",
  "element_connection_cmd", "list_element_connections",
  "element_connection", "list_opt_polygon_object_cmds",
  "opt_polygon_object_cmd", "remove_side", "$@7",
  "remove_element_specifier_list", "side_name", "element_specifier_list",
  "element_specifier", "incl_element_list_stmt", "excl_element_list_stmt",
  "just_an_element_list", "list_element_specs", "element_spec",
  "prev_region_stmt", "prev_region_type", "patch_statement", "patch_type",
  "in_obj_define_surface_regions", "list_in_obj_surface_region_defs",
  "in_obj_surface_region_def", "$@8", "$@9", "voxel_list_def", "$@10",
  "tet_element_connection_cmd", "element_connection_tet",
  "list_tet_arrays", "periodic_box_def", "$@11", "$@12", "box_def", "$@13",
  "$@14", "periodic_x_def", "periodic_y_def", "periodic_z_def",
  "periodic_traditional", "opt_aspect_ratio_def",
  "existing_obj_define_surface_regions",
  "list_existing_obj_surface_region_defs",
  "existing_obj_surface_region_def", "$@15", "$@16", "$@17", "new_region",
  "list_opt_surface_region_stmts", "opt_surface_region_stmt",
  "set_surface_class_stmt", "mod_surface_regions",
  "list_existing_surface_region_refs", "existing_surface_region_ref",
  "$@18", "output_def", "$@19", "output_buffer_size_def",
  "output_timer_def", "step_time_def", "iteration_time_def",
  "real_time_def", "list_count_cmds", "count_cmd", "count_stmt", "$@20",
  "custom_header_value", "custom_header", "exact_time_toggle",
  "list_count_exprs", "single_count_expr", "count_expr", "count_value",
  "$@21", "$@22", "file_arrow", "outfile_syntax",
  "existing_rxpn_or_molecule", "existing_molecule_required_orient_braces",
  "count_syntax", "count_syntax_1", "count_syntax_2", "count_syntax_3",
  "count_syntax_periodic_1", "count_syntax_periodic_2",
  "count_syntax_periodic_3", "count_location_specifier", "opt_hit_spec",
  "hit_spec", "opt_custom_header", "viz_output_def", "$@23",
  "list_viz_output_cmds", "viz_output_maybe_mode_cmd", "viz_mode_def",
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
     505,   506,   507,   508,   509,   510,    61,    38,    58,    43,
      45,    42,    47,   511,    94,    91,    93,    59,    39,    44,
     123,   125,    40,    41,    62,    60,    64
};
# endif

#define YYPACT_NINF -1044

#define yypact_value_is_default(Yystate) \
  (!!((Yystate) == (-1044)))

#define YYTABLE_NINF -398

#define yytable_value_is_error(Yytable_value) \
  0

  /* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
     STATE-NUM.  */
static const yytype_int16 yypact[] =
{
    2189,   -28,    10,    20,    58,    63,   153,   178,   151,   166,
     178,   178,   181,   186,   226,   229,   234,   -55,    61,   198,
     269, -1044,   273,   279,   286,   292,   304,   306,   308,   317,
     318,   323, -1044, -1044, -1044,   333,   314,   328,   370,   374,
     363,   389,   388,   395,   406,   414, -1044,   391,   402,   403,
     675,  2189, -1044,   -25, -1044, -1044,   429, -1044, -1044,   606,
   -1044, -1044, -1044, -1044, -1044, -1044, -1044, -1044, -1044, -1044,
   -1044, -1044,   433, -1044, -1044, -1044, -1044, -1044, -1044, -1044,
   -1044, -1044, -1044, -1044,   293, -1044, -1044, -1044, -1044,   522,
   -1044, -1044, -1044, -1044, -1044, -1044, -1044, -1044, -1044,   197,
     197,   -48,  1971,   -48,  1971, -1044, -1044,   426,   178,   178,
   -1044,   427, -1044,   431, -1044,   178,   178,   -48,   -47,  1971,
     178,   178,   178,   -48,   178,  1971,  1971,   197,  1971,  1971,
    1971,  1971,   379,   178,  1226, -1044,   660,   -48,   -48,  1653,
    1971,   548,  1971,   178,  1971,  1971,  1971, -1044,   627,  1964,
   -1044, -1044,  1483,   434,    77,   333, -1044, -1044,   333, -1044,
     333, -1044, -1044,   333,   333,   333, -1044, -1044, -1044, -1044,
   -1044, -1044, -1044, -1044,   435, -1044, -1044, -1044, -1044, -1044,
     456, -1044, -1044,   443,   449,   451,   461,   463,   472,   477,
     478, -1044,   479,   489,   491,   493,   499, -1044, -1044, -1044,
   -1044,   508, -1044,   510,   511,   520,   523,  1971,  1971,  1971,
   -1044,   121, -1044, -1044, -1044, -1044, -1044,  1184,    40,    56,
    -110, -1044, -1044,   487, -1044,   -45, -1044, -1044,   135, -1044,
   -1044, -1044,   -18, -1044, -1044, -1044,    -6, -1044,   456,   537,
   -1044, -1044,  1169, -1044,   526,   465,   533,   456, -1044,   652,
   -1044,  1169,  1169, -1044,  1169,  1169,  1169,  1169, -1044, -1044,
   -1044,   546,   540,   117, -1044,   559,   561,   565,   566,   569,
     570,   571,   575,   578,   588,   591,   592,   596,   597,   600,
     602,   604,   612,   364, -1044,   614,   456, -1044,   577, -1044,
    1169,  1169,   617, -1044,  1169, -1044,   610,  1169,  1169,  1169,
     733,   619,   746,   626,   629,   630,   631,   635,   636,   640,
     641,   645,   653,   654,   655,   657,   658,   659,   662,   663,
     665,  1030, -1044,  1755,   467, -1044, -1044,  1169,  1238, -1044,
    1257,   537,   -48, -1044, -1044, -1044, -1044,   884,   178, -1044,
     684, -1044,   684,   -48,   -48,  1971,  1971,  1971,  1971,  1971,
    1971,  1971,  1971,  1971,  1971,  1971,  1971,  1971,  1971,  1971,
    1971,   -48,  1971,   670,   670,   517, -1044, -1044,  1971,  1971,
    1971,  1971,  1971, -1044,  1971, -1044,   683,   690,   462, -1044,
   -1044, -1044, -1044, -1044,  1971, -1044,   313, -1044, -1044, -1044,
   -1044, -1044,   178,   178,    62,   -19, -1044, -1044, -1044,   661,
   -1044, -1044, -1044,   -48,   -48,   178, -1044, -1044, -1044,   197,
     197,   197,    15,   197,   197,  1435,   197,   197,   197,  1971,
     197,    15,   197,   197,   197,    15,    15, -1044, -1044,    77,
      99, -1044,  1971,   -36,   -48,   691,    73, -1044,   -48,   695,
     218, -1044,   -39,   -39,   -39,  1971,   -39,  1971,   -39,   -39,
    1971,   -39,   -39,   -39,   -39,   -39,   -39,   -39,   -39,   -39,
   -1044, -1044,  1971,    94, -1044,  1169,   669,   697,   785, -1044,
     394,   178, -1044, -1044,   759,   688,   738,   682,   899, -1044,
   -1044,   527,   545,   581,   605,   644,   671,   681,   705,   729,
     737,  1087,  1100,  1118,  1130,   761,   777,   221,   840, -1044,
     235,   235,   670,   670, -1044,  1198,  1971,  1971,   706,   707,
     747,   871, -1044, -1044, -1044, -1044,   487, -1044, -1044,   710,
     -15, -1044,   -33, -1044, -1044, -1044,   -50,   736,   744,   748,
     749,   750, -1044,    41,   178, -1044,   721,   742, -1044, -1044,
   -1044, -1044, -1044, -1044, -1044, -1044, -1044, -1044,  1169, -1044,
   -1044, -1044, -1044,  1169, -1044, -1044, -1044, -1044, -1044, -1044,
   -1044,   740, -1044,  1605, -1044,  1169,   755,   757,   760,   -60,
   -1044, -1044, -1044, -1044,    43,   762,   745,   -53, -1044, -1044,
   -1044, -1044,   456,   178,   768, -1044,   766, -1044, -1044, -1044,
   -1044, -1044, -1044,  1169, -1044,  1169, -1044, -1044,  1169, -1044,
   -1044, -1044, -1044, -1044, -1044, -1044, -1044, -1044,   251, -1044,
    1755,   -48,    77,   -16,   -34, -1044,   770,   682,    77,   763,
   -1044,   776,   779,   772,   788,   789,   778,   791,   800,   803,
   -1044, -1044,   804,   793,   682, -1044,   808, -1044, -1044, -1044,
   -1044, -1044,   798, -1044,   145, -1044, -1044, -1044, -1044, -1044,
   -1044, -1044, -1044, -1044, -1044,  1971,  1971,  1971,  1971, -1044,
   -1044, -1044, -1044,  1971,  1169,  1169,  1971,  1971, -1044,   954,
   -1044, -1044,   812, -1044, -1044,   710, -1044,   710, -1044, -1044,
     108, -1044,  1971,  1803,  1971,  1971,  1971, -1044,   178,   806,
     815, -1044, -1044, -1044, -1044, -1044,  -142, -1044, -1044, -1044,
     807,   214, -1044, -1044,   -14,    77, -1044, -1044,   537, -1044,
      77,  1971,    77,   818,   826, -1044,    -2, -1044, -1044, -1044,
   -1044,   254, -1044, -1044, -1044,   -48,   -44, -1044, -1044, -1044,
   -1044,   827,    77,   832,   842,  1971, -1044,   456,   820,   825,
     333,   839,   847,   849, -1044, -1044, -1044, -1044,    34,   682,
   -1044, -1044,  -165,    77, -1044,  1971,  1971,   981,    77,   178,
     178,  1971,   178,  1971,    77,   990,   -34, -1044,  1933,    77,
   -1044, -1044,   864,   948,   963,   969,  1219,  1169,  1169,   851,
     837,   -65, -1044, -1044,   -50,  1037,   853, -1044, -1044,  1169,
   -1044,  1169, -1044,  1169,  1169,  1169,   858,   178,   178, -1044,
   -1044,    12, -1044, -1044,   859, -1044, -1044, -1044, -1044,   950,
   -1044,  1169, -1044,   594,   197,   242, -1044, -1044, -1044,   456,
     846,   852,   857,   -32, -1044, -1044, -1044, -1044,   178, -1044,
    1933,   873,   -22,   136, -1044,    77, -1044,    77,  1933,    77,
   -1044, -1044, -1044, -1044, -1044, -1044,  -176,   546, -1044,   310,
     -34, -1044, -1044, -1044, -1044,   -23,   -34,  1169,  1169,   885,
   -1044, -1044,    77,   120, -1044,  1169, -1044, -1044,  1169, -1044,
     887, -1044,  1268, -1044, -1044, -1044, -1044,   168, -1044,    -8,
   -1044, -1044, -1044, -1044,  1971,  1971, -1044, -1044,  1605,  1605,
   -1044, -1044,   537,   264, -1044,   178, -1044,  1971,   487,   889,
     131, -1044,   169, -1044,   487, -1044,   870,   178,   890, -1044,
   -1044, -1044,   456, -1044, -1044, -1044,   882,   876, -1044,   242,
     242, -1044,   -64, -1044,   294, -1044,    25, -1044,    25, -1044,
   -1044, -1044,  1268, -1044, -1044, -1044,  1933,   896,   897,   901,
     892,  1971,  1139, -1044,   898, -1044, -1044,   316,  -176,  -176,
    -176, -1044, -1044, -1044, -1044,  1971, -1044, -1044, -1044,  1971,
   -1044, -1044,   895,   902,   230, -1044, -1044, -1044,  1169,  1169,
   -1044, -1044, -1044,  1037, -1044,  1169, -1044,  1971, -1044, -1044,
   -1044, -1044, -1044,   423, -1044,   197,  1003,   904,  1971,   242,
     908, -1044,   405,   242,  -149,   -48,   242,   242,   242,   242,
   -1044, -1044, -1044, -1044,     9, -1044,   900,     8,    11, -1044,
     903, -1044,    77,  1971,    77, -1044,   377,   918, -1044,   -34,
    1971, -1044,   914,   914, -1044,   422,   205,   178, -1044, -1044,
     916,  1169,   924, -1044, -1044,   927, -1044, -1044,   423, -1044,
   -1044, -1044, -1044,   929, -1044,   934, -1044,   935,  1026,   216,
    1157,   819,   216, -1044, -1044,   920,   923,   930,   -48,   456,
     330,   330, -1044, -1044, -1044, -1044,    16,   945, -1044, -1044,
   -1044, -1044,   945, -1044, -1044,     7, -1044,  1169, -1044, -1044,
    1971, -1044, -1044,  1169,   949, -1044,   955,   177, -1044,   943,
    1175, -1044,   953,   961, -1044, -1044,   178,    77,   197,   958,
    1052,   966,   968,   970,   972,   982, -1044, -1044, -1044, -1044,
   -1044, -1044, -1044, -1044,   983, -1044, -1044,   979, -1044, -1044,
   -1044, -1044, -1044,  1971, -1044, -1044, -1044, -1044, -1044,  1169,
      -8,  1971,  1971, -1044, -1044, -1044, -1044, -1044, -1044, -1044,
   -1044, -1044, -1044,   481,   985, -1044,   423, -1044,   989, -1044,
    1313,  1313,   192, -1044,   993, -1044,   197,   998, -1044,   -40,
   -1044,   -40,   -40, -1044, -1044, -1044,  1169, -1044,  1019,   115,
     423,  1971, -1044,  1313,   275,   281, -1044,    77, -1044,   197,
     988, -1044,   546, -1044,   994,   995,   996,   -34, -1044,  1010,
     423,  1169, -1044, -1044, -1044, -1044, -1044, -1044,   148, -1044,
     148, -1044,   148, -1044, -1044,  1971, -1044, -1044, -1044, -1044,
   -1044, -1044, -1044, -1044, -1044,  1000, -1044,  1000,  1000,  1064,
     232,   450, -1044, -1044, -1044, -1044, -1044
};

  /* YYDEFACT[STATE-NUM] -- Default reduction number in state STATE-NUM.
     Performed when YYTABLE does not specify something else to do.  Zero
     means the default is an error.  */
static const yytype_uint16 yydefact[] =
{
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,   321,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,   219,   220,   221,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,    27,     0,     0,     0,
       0,     2,     3,   329,     5,     6,     0,     7,   112,     0,
     113,   114,   115,   116,   117,   118,     8,     9,    10,    11,
      14,    12,     0,    15,   222,   223,    16,   243,   244,    17,
      18,    20,    19,   323,     0,   324,   325,   346,   345,     0,
     327,   328,    13,   326,    21,    22,    23,    24,    25,     0,
       0,     0,     0,     0,     0,   229,   224,     0,     0,     0,
     311,     0,   230,     0,   245,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,   330,     0,     0,     0,     0,
       0,   491,     0,     0,     0,     0,     0,   563,     0,     0,
       1,     4,     0,     0,     0,     0,   365,   366,     0,   367,
       0,   364,   368,     0,     0,     0,    35,    37,    39,    38,
      34,    36,   201,   200,     0,   108,    26,   107,   111,   184,
      28,   105,   106,     0,     0,     0,     0,     0,     0,     0,
       0,    70,     0,     0,     0,     0,     0,    93,    95,    94,
      71,     0,    96,     0,     0,     0,     0,     0,     0,     0,
      74,   189,    66,    68,    69,    67,   185,   192,   189,     0,
       0,   226,   240,    40,   292,     0,   273,   275,   293,   290,
     313,   249,     0,   247,    29,   474,     0,   472,     0,   211,
     212,   213,   206,   123,     0,     0,     0,    55,   329,     0,
     322,   207,   199,   187,   214,   215,   216,   217,   209,   210,
     208,     0,     0,     0,   485,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,   136,     0,   124,   125,     0,   204,
     203,   205,     0,   489,   197,    60,     0,   196,   198,   202,
     567,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,   161,     0,    61,    58,    59,     0,    72,    56,
      73,    57,     0,    65,   218,    62,    63,     0,     0,   347,
       0,   362,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,   104,   103,     0,   191,   190,     0,     0,
       0,     0,     0,   186,     0,   188,     0,     0,   233,   225,
     227,    43,    48,    49,     0,   242,    41,    44,    45,    42,
     272,   274,     0,     0,     0,     0,   252,   246,   248,     0,
     471,   473,   122,     0,     0,     0,   487,   484,   486,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,   135,   137,     0,
       0,   133,     0,     0,     0,     0,     0,   568,     0,     0,
       0,   608,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     160,   162,     0,     0,    51,    53,     0,     0,   329,   342,
       0,   332,   339,   341,     0,     0,     0,     0,     0,   125,
     109,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,    75,
      98,    99,   100,   101,   102,   193,     0,     0,     0,     0,
     236,     0,    46,    47,   291,   251,    40,   294,   276,     0,
       0,   283,     0,   285,   284,   286,     0,     0,     0,     0,
       0,     0,   310,     0,     0,   125,     0,     0,   479,   157,
     138,   145,   153,   159,   158,   140,   147,   148,   155,   154,
     156,   144,   141,   143,   139,   150,   146,   149,   142,   152,
     151,     0,    31,     0,   130,   492,     0,     0,     0,     0,
     493,   494,   495,   125,     0,     0,     0,     0,   565,   573,
     572,   574,   607,     0,     0,   609,     0,   183,   181,   182,
     163,   168,   169,   167,   166,   172,   171,   173,   174,   175,
     177,   164,   165,   178,   179,   180,   170,   176,     0,    64,
       0,     0,     0,     0,     0,   340,     0,     0,     0,     0,
     449,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     382,   383,     0,     0,   332,   369,     0,   374,   384,   385,
     386,   387,     0,   398,     0,    91,    88,    87,    89,    83,
      85,    76,    82,    77,    78,     0,     0,     0,     0,    84,
      90,    97,    86,     0,   232,   231,     0,     0,   237,   238,
      50,   295,   279,   277,   278,     0,   280,     0,   298,   299,
       0,   296,     0,     0,     0,     0,     0,   261,     0,     0,
       0,   259,   260,   250,   253,   254,     0,   255,   264,   478,
       0,     0,   134,    30,     0,     0,   129,   127,   128,   126,
       0,     0,     0,     0,     0,   504,     0,   499,   501,   502,
     503,     0,   570,   571,   569,     0,     0,   564,   566,   611,
     612,   610,     0,     0,     0,     0,    52,   121,     0,     0,
       0,     0,     0,     0,   331,   338,   333,   334,     0,   332,
     401,   402,     0,     0,   332,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,   370,     0,     0,
     408,   110,     0,     0,     0,     0,   194,   235,   234,     0,
       0,     0,   281,   282,     0,     0,   287,   300,   301,   314,
     320,   319,   318,   315,   317,   316,     0,     0,     0,   263,
     262,     0,   475,   131,     0,   488,   482,   480,   481,   467,
     497,   496,   498,     0,     0,     0,   490,   500,   132,   575,
       0,     0,     0,     0,   577,   579,   580,   581,     0,   614,
       0,     0,   617,     0,   119,     0,   343,     0,     0,     0,
     352,   353,   356,   354,   351,   355,     0,   350,   357,   349,
       0,   400,   403,   452,   453,     0,     0,   392,   393,     0,
     372,   373,     0,     0,   394,   388,   312,   380,   379,   378,
       0,   363,   371,   376,   375,   377,   407,     0,   405,   332,
      79,    80,    92,    81,     0,     0,   228,   297,     0,     0,
     309,   307,   308,     0,   304,     0,   289,     0,    40,     0,
       0,   267,     0,   269,    40,   256,     0,     0,     0,   455,
     506,   507,   508,   509,   510,   523,     0,     0,   526,     0,
       0,   514,     0,   511,   561,   515,     0,   586,     0,   576,
     578,   613,    65,    32,   615,    33,     0,     0,     0,     0,
       0,     0,   469,   332,     0,   336,   335,     0,     0,     0,
       0,   348,   451,   454,   450,     0,   396,   381,   395,     0,
     404,   406,     0,     0,     0,   409,   410,   411,   195,   239,
     305,   306,   302,     0,   288,   258,   241,     0,   265,   268,
     266,   270,   257,     0,   483,     0,   461,     0,     0,     0,
       0,   521,     0,     0,     0,     0,     0,     0,     0,     0,
     513,   603,   605,   604,     0,   600,     0,     0,     0,   594,
       0,   616,     0,     0,     0,   606,     0,     0,   458,     0,
       0,   358,   359,   360,   361,     0,     0,     0,   412,   399,
       0,   271,     0,   442,   439,     0,   441,   438,   476,   423,
     425,   426,   427,     0,   428,     0,   468,     0,   463,     0,
       0,     0,     0,   516,   512,     0,     0,   528,     0,   562,
     517,   518,   519,   520,   599,   601,     0,   584,   582,   590,
     589,   585,   584,   593,   595,     0,   619,   618,   620,    54,
       0,   408,   344,   337,     0,   389,     0,     0,   444,     0,
       0,   303,     0,     0,   424,   479,     0,     0,     0,     0,
     465,     0,   534,     0,     0,     0,   536,   537,   538,   539,
     540,   541,   525,   522,     0,   529,   532,   530,   533,   505,
     597,   598,   602,     0,   588,   587,   591,   592,   596,   470,
     459,     0,     0,   443,   445,   446,   422,   419,   417,   418,
     420,   421,   416,   434,     0,   436,   414,   415,   431,   432,
       0,     0,     0,   437,     0,   462,     0,     0,   456,     0,
     535,     0,     0,   524,   527,   531,   583,   332,     0,     0,
       0,     0,   413,     0,     0,     0,   477,     0,   464,     0,
       0,   548,   550,   549,   551,   551,   551,     0,   390,     0,
     447,   435,   433,   430,   429,   440,   466,   457,     0,   544,
       0,   542,     0,   543,   460,     0,   479,   558,   560,   555,
     557,   554,   559,   556,   553,   551,   552,   551,   551,     0,
       0,     0,   547,   545,   546,   391,   448
};

  /* YYPGOTO[NTERM-NUM].  */
static const yytype_int16 yypgoto[] =
{
   -1044, -1044, -1044,  1221,  -749,     0,  -102,  -108,  -131,  -278,
    -786,   -95,  -492, -1044,   891,   905,   171, -1044,   676, -1044,
   -1044,  1132,  -143,  -138,  -139, -1044,   -61,  -340,  -122,  -140,
   -1044,  -116,   -83,  -114, -1044, -1044, -1044, -1044, -1044, -1044,
     495,   -99,  -433, -1044, -1044, -1044, -1044, -1044, -1044, -1044,
   -1044,  1002,  1721,    83, -1044, -1044,   971,   528, -1044,  1075,
   -1044, -1044, -1044, -1044, -1044, -1044, -1044, -1044, -1044,   -35,
   -1044, -1044, -1044, -1044, -1044, -1044, -1044, -1044,  -432, -1044,
   -1044, -1044, -1044,   -54, -1044,   387, -1044, -1044, -1044, -1044,
   -1044, -1044,   767, -1044, -1044,  -699, -1044, -1044,  1070,  -344,
    -400, -1044, -1044, -1044, -1044, -1044, -1044, -1044, -1044,   907,
   -1044, -1044, -1044,   519, -1044, -1044, -1044,   334,  -280, -1044,
   -1044, -1044, -1044, -1044, -1044, -1044, -1044,  -307,  -109,   -29,
    -752,  -577, -1044, -1044,  1182, -1044,   843, -1044, -1044, -1044,
   -1044, -1044, -1044,  -603, -1044, -1044, -1044,   693, -1044,  -608,
   -1044, -1044, -1044, -1044, -1044, -1044, -1044,   445, -1044, -1044,
   -1044,   973,   568, -1044, -1044, -1044,   436,   236, -1044, -1044,
   -1044, -1044, -1044, -1043,  -990, -1044, -1044, -1044,  -531,   154,
   -1044, -1044, -1044, -1044, -1044, -1044,   231, -1044, -1044, -1044,
   -1044, -1044,   466, -1044, -1044, -1044, -1044, -1044, -1044, -1044,
   -1044, -1044, -1044, -1044, -1044, -1044, -1044,  1093, -1044, -1044,
   -1044,   802, -1020, -1044, -1044, -1044, -1044,  1069, -1044, -1044,
   -1044, -1044, -1044, -1044, -1044, -1044, -1044,   622, -1044, -1044,
   -1044, -1044, -1044, -1044,   346,  -431, -1044, -1044, -1044, -1044,
   -1044, -1044, -1044,   288, -1044, -1044, -1044, -1044, -1044, -1044,
    -507,  -538, -1044, -1044, -1044, -1044, -1044, -1044, -1044,   765,
   -1044, -1044, -1044, -1044,   529, -1044,   271, -1044, -1044, -1044,
   -1044, -1044, -1044,   337, -1044, -1044, -1044,   359,  -866, -1044,
   -1044, -1044,   910,   538, -1044, -1044, -1044, -1044, -1044
};

  /* YYDEFGOTO[NTERM-NUM].  */
static const yytype_int16 yydefgoto[] =
{
      -1,    50,    51,    52,   177,   210,   179,   261,   848,   933,
     934,   539,   385,   386,   387,   388,   389,   463,   464,    54,
      55,    56,   890,   562,   335,   336,   327,   212,   213,   891,
     214,   215,   238,   181,   182,    57,    58,    59,   738,    60,
     244,   287,   430,   709,    61,    62,    63,    64,    65,    66,
     283,   284,   540,   545,    67,   321,   322,   590,    68,   373,
     218,    69,    70,    71,    72,    73,    74,    75,   220,   106,
     107,   113,   378,   510,   669,   780,   223,   899,   224,    76,
      77,    78,   232,   114,   396,   516,   533,   694,   695,   696,
     801,   697,   806,   900,   902,   901,    79,   225,   226,   781,
     521,   522,   523,   524,   525,   526,   896,   227,   228,   229,
     394,   517,   680,   681,   786,   787,   788,   893,   894,    80,
     111,   867,   395,   792,    81,   124,    82,    83,    84,   338,
     745,   614,   746,   747,    85,   471,   472,   473,   943,    86,
      87,   474,   617,   849,    88,   477,   164,   634,   875,   635,
     636,   637,   638,   639,   640,   641,   863,   864,    89,    90,
     770,   476,   751,   752,   643,   877,   878,   879,   965,   966,
    1090,  1144,  1145,  1038,  1039,  1040,  1041,  1147,  1148,  1149,
    1042,  1043,  1044,  1045,   967,  1087,  1088,  1170,  1206,    91,
     754,   620,   854,   855,    92,   986,  1180,    93,  1081,  1167,
    1048,  1100,  1158,   909,  1018,    94,   236,   237,   399,   906,
    1095,  1089,   704,   807,   808,    95,   263,   264,   538,    96,
     433,   293,   569,   570,   571,   572,   716,   717,   718,   815,
     913,   719,   720,   922,   923,   924,   925,   987,   990,  1058,
    1119,  1103,  1104,  1105,  1106,  1107,  1108,  1109,  1110,  1111,
    1184,  1199,  1216,  1000,    97,   300,   577,   436,   437,   578,
     579,   580,   581,   823,   824,   825,  1124,  1007,  1071,  1072,
    1128,   826,  1008,  1009,  1122,   827,  1004,  1005,  1006,    98,
     302,   440,   441,   730,   731,   586,   734,   832,   940
};

  /* YYTABLE[YYPACT[STATE-NUM]] -- What to do in state STATE-NUM.  If
     positive, shift that token.  If negative, reduce the rule whose
     number is the opposite.  If YYTABLE_NINF, syntax error.  */
static const yytype_int16 yytable[] =
{
      53,   216,   262,   239,   172,   173,   136,   105,   235,   325,
     110,   112,   329,   326,   871,   249,   334,  1001,   180,  1001,
     180,  1067,  1126,   575,   671,   904,   767,   587,  1120,   713,
     328,   469,   253,  1001,   174,   527,   330,   260,   331,   288,
     247,   211,   543,   217,   962,   687,   644,  1146,  1094,   175,
     520,    53,   945,   588,   286,   286,  -120,   766,   242,   820,
     722,   233,  1010,    46,   251,   252,   566,   254,   255,   256,
     257,   820,   723,   221,    46,  1152,   840,   688,   290,   291,
     937,   294,   576,   297,   298,   299,   799,   713,   675,   166,
     677,   240,   241,   821,   679,   678,   846,   841,   951,   903,
     323,   178,   701,   178,   954,   821,   851,  1055,   105,   222,
    1056,   366,   376,   377,   800,   112,   234,   178,   689,   690,
     243,   243,   243,   178,   248,  1057,   337,  1190,   235,    46,
     528,   339,   262,   234,   340,   341,   342,   178,   178,   714,
     721,   767,  1010,   295,   842,  -397,   363,   364,   365,   575,
    1011,   561,   324,  1189,   333,  1207,  1094,   167,  1208,  1209,
     741,   379,   742,  1002,   469,  1002,   176,   168,   169,   529,
    1210,  1211,   850,   689,   690,   672,   567,   856,   398,  1002,
    1003,   822,  1003,   367,   963,   380,  1220,   724,   568,    46,
     938,    46,   366,   822,    46,   673,  1003,   714,   743,    46,
    1094,   979,   939,   981,   804,   993,    46,   994,   576,   674,
     715,   530,   531,   589,  1181,  1212,    46,   120,   727,   843,
     105,    46,   176,    46,    46,   222,   390,   691,    99,   470,
     466,   -60,   112,    46,  1213,  1214,   234,   744,   844,   929,
     845,   676,   323,   947,   479,   673,    46,    46,   952,   180,
     170,    46,   532,   397,   135,    46,   796,   805,  1069,   674,
     286,   480,   465,   234,   367,   400,   100,  1082,   715,   816,
     171,   166,   323,    46,   692,   782,   101,   783,   497,  1068,
    1064,   323,  1073,   915,   481,   482,   483,   484,   485,   486,
     487,   488,   489,   490,   491,   492,   493,   494,   495,   496,
    1101,   498,   964,  1101,   535,   536,   846,   500,   501,   502,
     503,   504,   693,   505,   102,   916,    46,   544,   155,   103,
     286,   286,   518,   511,   689,   690,   544,   861,   862,   584,
     544,   544,   178,   121,   739,   573,   156,   519,   468,   167,
     750,   191,   323,   178,   178,  1022,  1023,  1024,   941,   168,
     169,   286,   679,   439,   548,   582,    46,   157,   553,    46,
     609,   178,   470,   610,   689,   690,  1019,   784,   563,   905,
      46,   565,   564,   785,   368,   369,   370,   371,   265,   372,
     368,   369,   370,   371,   593,   372,   595,    46,   407,   598,
     266,   957,   222,   515,   392,   368,   369,   370,   371,   267,
     372,   608,   978,   178,   178,   537,   976,    46,    46,   104,
     804,   393,   982,   323,   563,  1086,    46,    46,   771,   155,
     706,   108,   268,   707,   741,   200,   742,   809,  1197,   333,
     176,   862,   170,   323,   178,  1204,   109,   156,   178,   960,
     980,   328,   269,   270,   158,   664,   665,   330,  1133,   708,
     804,   115,   171,   166,   829,    46,   116,  1207,   157,   917,
    1208,  1209,   743,  1176,   368,   369,   370,   371,   271,   372,
     122,   468,  1210,  1211,   750,   921,   918,   159,   344,   735,
     860,    46,   117,   563,   160,   118,   869,   803,   991,   992,
     119,  1032,  1033,  1034,   661,   272,   370,   371,   161,   372,
     162,  1029,   919,  1226,   555,   740,   508,   509,   559,   560,
     368,   369,   370,   371,   920,   372,  1035,  1212,  1036,  1037,
     222,   167,   222,   563,   273,   123,   222,   818,   737,   125,
     972,   168,   169,   973,   699,   126,  1213,  1214,   274,   275,
     276,  1193,   127,   163,  1173,   613,   277,  1194,   128,   465,
    1173,   278,   995,   996,   997,   998,   999,   942,  1051,   944,
     129,   946,   130,   324,   131,  1060,  1061,  1062,  1063,   948,
     949,   950,   810,   132,   812,   948,   949,   950,   159,   921,
     921,   382,   383,   729,   956,   160,   137,   279,   133,  1021,
    1187,   998,   999,   134,   772,   773,   774,   775,   258,   161,
     138,   162,   776,   135,   280,   777,   778,   281,   970,   971,
     282,   178,   333,   234,   170,   853,   245,   246,   333,  1174,
    1175,   789,   791,   793,   794,   795,   139,   259,   874,   873,
     140,   876,  1084,   141,   171,   427,   368,   369,   370,   371,
     847,   372,   819,  1079,   163,   142,   328,  1201,  1203,   921,
     811,   144,   330,   921,  1185,  1186,   921,   921,   921,   921,
     143,   147,   145,   328,   996,   997,   998,   999,   166,   330,
     146,   892,   148,   149,   833,   150,   174,  1222,  1053,  1223,
    1224,   368,   369,   370,   371,   152,   372,   153,   222,   154,
     935,   175,   165,  1085,   857,   858,   219,   230,   935,   285,
     865,   231,   868,   292,   301,   333,   332,   343,   328,  1171,
     333,   836,   333,   344,   330,   345,   328,   853,   911,   914,
     621,   346,   330,   347,  -111,   178,   -74,   -74,   -74,   -74,
     912,   -74,   333,   348,   403,   349,   167,   622,   910,   876,
     368,   369,   370,   371,   350,   372,   168,   169,   234,   351,
     352,   353,   333,   333,   381,   382,   383,   384,   333,   222,
     222,   354,   866,   355,   333,   356,   328,   328,   872,   333,
     623,   357,   330,   330,   892,   892,   368,   369,   370,   371,
     358,   372,   359,   360,   222,   324,   368,   369,   370,   371,
     499,   372,   361,   624,  -105,   362,   935,   898,   898,   402,
     645,   222,   404,   158,   368,   369,   370,   371,   176,   372,
     406,   405,   625,   178,   328,   409,   626,   410,   646,  1154,
     330,   411,   412,   968,   969,   413,   414,   415,   729,   170,
     932,   416,   627,    46,   417,   333,   975,   333,   932,   333,
     368,   369,   370,   371,   418,   372,   234,   419,   420,   171,
     431,   328,   421,   422,   647,   333,   423,   330,   424,   892,
     425,   435,   333,   222,   368,   369,   370,   371,   426,   372,
     429,   628,   629,   432,  1076,   438,  1078,   333,   648,   434,
    1016,   439,   442,   630,   631,   443,   444,   445,   324,   324,
    1046,   446,   447,   632,  1025,   974,   448,   449,  1026,  1195,
     898,   450,   898,   368,   369,   370,   371,   515,   372,   451,
     452,   453,  1059,   454,   455,   456,  1031,   649,   457,   458,
    1215,   459,  1217,   467,  1218,   475,   534,  1050,  1121,   633,
     368,   369,   370,   371,   372,   372,   932,  1127,   611,   506,
     368,   369,   370,   371,   650,   372,   507,   574,   234,   234,
     234,   583,  1077,   612,   651,  -397,  1118,   616,   618,  1083,
     619,   642,   666,   667,   368,   369,   370,   371,   668,   372,
     518,   591,   592,   324,   594,   180,   596,   597,   652,   599,
     600,   601,   602,   603,   604,   605,   606,   607,   368,   369,
     370,   371,   682,   372,   702,   178,   368,   369,   370,   371,
     683,   372,   653,  1155,   684,   685,   686,  1070,   703,   705,
     654,   710,   333,   711,   333,   726,   712,   733,   725,  1129,
     368,   369,   370,   371,   732,   372,   748,   699,  1183,  1143,
    1183,  1183,   755,   753,   659,   756,   368,   369,   370,   371,
     183,   372,   757,   184,   758,   759,   303,   761,   760,  1102,
     660,  1182,  1102,  1182,  1182,   185,   762,   186,   178,   763,
     764,  1178,  1166,   765,   768,   187,   333,   304,   769,   779,
    1168,  1169,   673,   802,   813,   333,   797,   188,   996,   997,
     998,   999,   814,   305,  1196,   798,   828,   699,   830,  1143,
    1143,   831,  1113,   834,   835,   837,  1153,   333,   859,   368,
     369,   370,   371,   838,   372,   839,   870,   885,   886,   189,
    1191,   895,  1143,   662,   897,   907,   926,   190,   908,   174,
     306,   307,   927,   368,   369,   370,   371,   928,   372,   936,
     368,   369,   370,   371,   175,   372,   191,   880,   308,   309,
     983,   955,   670,   959,  1219,   977,   985,   988,   989,   192,
     193,   194,  1012,  1013,   310,   311,   312,  1014,  1017,   234,
     195,   234,   234,  1015,   196,  1027,   313,  1020,  1047,  1049,
     314,   315,  1028,  1052,  1080,   950,  1066,   333,   183,  1075,
    1092,   184,  1091,  1093,  1136,  1096,   316,   317,   318,   319,
    1097,  1098,  1099,   185,  1115,   186,  1137,  1116,   333,  1138,
     333,  1123,   333,   187,  1117,  1131,   197,   368,   369,   370,
     371,  1132,   372,  1135,  1156,   188,   198,   199,  1150,  1157,
     200,   881,   368,   369,   370,   371,  1151,   372,   368,   369,
     370,   371,   201,   372,   202,  1159,   882,   203,   384,  1161,
     265,  1162,   883,  1032,  1033,  1034,   204,   189,  1163,  1164,
     205,   176,   266,  1165,  1179,   190,  1172,   206,  1173,   744,
    1139,   267,  1177,  1198,  1200,  1202,  1205,   320,  1035,  1221,
    1036,  1037,   151,  1160,   191,   296,    46,   512,   368,   369,
     370,   371,  1140,   372,   268,   428,   736,   192,   193,   194,
    1188,   513,   461,   375,   984,   391,   207,   208,   195,   514,
     698,   460,   196,   887,   269,   270,   250,  1030,   958,   209,
     749,   888,   889,   961,   615,   478,   183,  1130,  1134,   184,
     852,   953,  1136,   368,   369,   370,   371,  1192,   372,   401,
     271,   185,   408,   186,  1137,  1225,   700,  1138,   817,  1054,
    1114,   187,   728,  1125,   197,  1074,   368,   369,   370,   371,
     585,   372,   930,   188,   198,   199,   655,   272,   200,   368,
     369,   370,   371,  1065,   372,     0,   931,     0,  1141,   656,
     201,     0,   202,     0,     0,   203,     0,   368,   369,   370,
     371,     0,   372,     0,   204,   189,   273,   657,   205,   368,
     369,   370,   371,   190,   372,   206,     0,     0,  1139,   658,
     274,   275,   276,     0,  1142,     0,     0,     0,   277,     0,
       0,     0,   191,   278,    46,     0,   368,   369,   370,   371,
    1140,   372,     0,  1112,     0,   192,   193,   194,   368,   369,
     370,   371,     0,   372,   207,   208,   195,     0,   183,     0,
     196,   184,   374,   368,   369,   370,   371,   209,   372,   279,
       0,     0,     0,   185,     0,   186,   663,   368,   369,   370,
     371,     0,   372,   187,     0,     0,   280,     0,     0,   281,
       0,     0,   282,     0,     0,   188,     0,   884,   368,   369,
     370,   371,   197,   372,     0,     0,   183,     0,     0,   184,
       0,     0,   198,   199,     0,     0,   200,   -68,   -68,   -68,
     -68,   185,   -68,   186,     0,     0,  1141,   189,   201,   166,
     202,   187,     0,   203,     0,   190,   -67,   -67,   -67,   -67,
       0,   -67,   204,   188,     0,     0,   205,   -74,   -74,   -74,
     -74,     0,   -74,   206,   191,     0,     0,     0,     0,     0,
       0,     0,  1142,     0,     0,     0,     0,   192,   193,   194,
       0,     0,    46,     0,     0,   189,     0,     0,   195,     0,
       0,     0,   196,   190,     0,   174,     0,     0,     0,     0,
       0,     0,   207,   208,     0,     0,     0,   167,     0,     0,
     175,     0,   191,     0,     0,   209,     0,   168,   169,     0,
       0,     0,     0,     0,     0,   192,   193,   194,     0,     0,
       0,     0,     0,     0,   197,     0,   195,     0,   183,     0,
     196,   184,     0,     0,   198,   199,     0,     0,   200,     0,
       0,     0,     0,   185,     0,   186,     0,     0,     0,     0,
     201,     0,   202,   187,     0,   203,     0,     0,     0,     0,
       0,     0,     0,     0,   204,   188,     0,     0,   205,     0,
       0,     0,   197,     0,     0,   206,   183,     0,     0,   184,
       0,     0,   198,   199,     0,     0,   200,     0,     0,     0,
     170,   185,     0,   186,    46,     0,     0,   189,   201,     0,
     202,   187,     0,   203,     0,   190,     0,   174,     0,     0,
     171,     0,   204,   188,   207,   208,   205,   176,     0,     0,
       0,     0,   175,   206,   191,     0,     0,   209,     0,     0,
       0,     0,     0,     0,     0,     0,     0,   192,   193,   194,
       0,     0,    46,     0,     0,   189,     0,     0,   195,     0,
       0,     0,   196,   190,     0,     0,     0,     0,     0,     0,
       0,     0,   207,   208,     0,     0,     0,     0,   323,     0,
       0,     0,   191,   289,     0,   209,     0,     0,   183,     0,
       0,   184,     0,     0,     0,   192,   193,   194,     0,     0,
       0,     0,     0,   185,   197,   186,   195,     0,     0,     0,
     196,     0,     0,   187,   198,   199,     0,     0,   200,     0,
       0,     0,     0,     0,     0,   188,     0,     0,     0,     0,
     201,     0,   202,     0,     0,   203,   183,     0,     0,   184,
       0,     0,     0,     0,   204,     0,     0,     0,   205,   176,
       0,   185,   197,   186,     0,   206,     0,   189,     0,     0,
       0,   187,   198,   199,     0,   190,   200,     0,     0,     0,
       0,     0,     0,   188,    46,     0,     0,     0,   201,     0,
     202,     0,     0,   203,   191,     0,     0,     0,     0,     0,
       0,     0,   204,     0,   207,   208,   205,   192,   193,   194,
       0,     0,     0,   206,     0,   189,     0,   209,   195,     0,
       0,     0,   196,   190,     0,     0,     0,     0,     0,     0,
       0,     0,    46,     0,     0,     0,     0,     0,     0,     0,
       0,     0,   191,     0,     0,     0,     0,     0,     0,     0,
       0,     0,   207,   208,     0,   192,   193,   194,     0,     0,
       0,     0,     0,     0,   197,   209,   195,     0,     0,     0,
     196,     0,     0,     0,   198,   199,   183,     0,   200,   184,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     201,   185,   202,   186,     0,   203,     0,     0,     0,     0,
       0,   187,     0,     0,   204,     0,     0,     0,   205,     0,
       0,     0,   197,   188,   183,   206,     0,   184,     0,     0,
     303,     0,   198,   199,     0,     0,   200,     0,     0,   185,
       0,   186,     0,     0,    46,     0,     0,     0,   201,   187,
     202,   304,     0,   203,     0,   189,     0,     0,     0,     0,
       0,   188,   204,   190,   207,   208,   205,   305,     0,     0,
     462,     0,     0,   206,     0,     0,     0,   209,     0,     0,
       0,     0,   191,     0,     0,     0,     0,     0,     0,   790,
       0,     0,    46,   189,     0,   192,   193,   194,     0,     0,
       0,   190,     0,     0,   306,   307,   195,     0,     0,     0,
     196,     0,   207,   208,     0,     0,     0,     0,     0,     0,
     191,     0,   308,   309,     0,   209,     0,     0,     0,     0,
       0,     0,     0,   192,   193,   194,     0,     0,   310,   311,
     312,     0,     0,     0,   195,     0,     0,     0,   196,     0,
     313,     0,   197,     0,   314,   315,     0,     0,     0,     0,
       0,     0,   198,   199,     0,     0,   200,     0,     0,     0,
     316,   317,   318,   319,     0,     0,     0,     0,   201,     0,
     202,   541,   542,   203,   546,   547,   549,   550,   551,   552,
     197,   554,   204,   556,   557,   558,   205,     0,     0,     0,
     198,   199,     0,   206,   200,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,   201,     0,   202,     0,
       0,   203,    46,     0,     0,     0,     0,     0,     0,     0,
     204,     0,     0,     0,   205,     0,     0,     0,     0,     0,
       0,   206,   207,   208,     1,     0,     0,     0,   323,     0,
       0,   320,     0,     0,     0,   209,     0,     0,     0,     0,
      46,     0,     0,     0,     0,     0,     0,     0,     0,     2,
       3,     4,     5,     6,     0,     0,     0,     0,     0,     0,
     207,   208,     0,     0,     0,     7,     8,     9,    10,    11,
      12,    13,     0,   209,     0,     0,     0,     0,    14,    15,
      16,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,    17,     0,     0,     0,     0,     0,
       0,     0,    18,    19,     0,     0,     0,     0,     0,     0,
       0,     0,     0,    20,     0,     0,     0,    21,     0,     0,
      22,     0,     0,     0,    23,    24,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,    25,    26,    27,
      28,    29,     0,     0,     0,     0,     0,     0,    30,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,    31,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
      32,    33,    34,    35,     0,     0,     0,     0,     0,     0,
       0,    36,    37,     0,     0,     0,    38,    39,     0,     0,
      40,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,    41,     0,     0,     0,     0,    42,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,    43,    44,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,    45,    46,     0,
       0,    47,     0,     0,    48,     0,     0,     0,     0,     0,
       0,     0,    49
};

static const yytype_int16 yycheck[] =
{
       0,   103,   133,   117,    99,   100,    35,     7,   116,   152,
      10,    11,   152,   152,   766,   124,   154,     8,   101,     8,
     103,    13,    15,    76,   516,    13,   634,    66,    12,    89,
     152,   338,   127,     8,    82,    54,   152,   132,   152,   138,
     123,   102,    27,   104,    52,     4,   479,  1090,  1038,    97,
     394,    51,   838,    92,   137,   138,    81,   634,   119,   103,
      17,   115,   928,   239,   125,   126,   102,   128,   129,   130,
     131,   103,    29,   108,   239,  1095,    42,    36,   139,   140,
     102,   142,   135,   144,   145,   146,   228,    89,   520,    74,
     522,   138,   139,   137,   526,   145,   272,    63,   850,   798,
     265,   101,   535,   103,   856,   137,   271,   256,   108,   109,
     259,    71,    56,    57,   256,   115,   116,   117,   132,   133,
     120,   121,   122,   123,   124,   274,   155,  1170,   236,   239,
     149,   160,   263,   133,   163,   164,   165,   137,   138,   199,
     573,   749,  1008,   143,   110,   170,   207,   208,   209,    76,
     936,   429,   152,    38,   154,     7,  1146,   142,    10,    11,
     194,   271,   196,   154,   471,   154,   214,   152,   153,   188,
      22,    23,   749,   132,   133,   519,   212,   754,   232,   154,
     171,   225,   171,   143,   192,   220,  1206,   144,   224,   239,
     212,   239,    71,   225,   239,   260,   171,   199,   232,   239,
    1190,   900,   224,   902,   218,   269,   239,   271,   135,   274,
     270,   230,   231,   252,   254,    67,   239,   272,   271,   185,
     220,   239,   214,   239,   239,   225,   271,   186,   256,   338,
     332,   256,   232,   239,    86,    87,   236,   271,   204,   271,
     206,   274,   265,   846,   343,   260,   239,   239,   271,   332,
     235,   239,   271,   271,   270,   239,   688,   271,  1007,   274,
     343,   344,   323,   263,   143,   271,   256,  1019,   270,   271,
     255,    74,   265,   239,   233,   675,   256,   677,   361,   271,
     271,   265,   271,    41,   345,   346,   347,   348,   349,   350,
     351,   352,   353,   354,   355,   356,   357,   358,   359,   360,
    1049,   362,   879,  1052,   403,   404,   272,   368,   369,   370,
     371,   372,   271,   374,   256,    73,   239,   412,    25,   256,
     403,   404,   260,   384,   132,   133,   421,   759,   760,   111,
     425,   426,   332,   272,   612,   434,    43,   275,   338,   142,
     618,    99,   265,   343,   344,   948,   949,   950,   212,   152,
     153,   434,   784,   135,   415,   438,   239,    64,   419,   239,
     266,   361,   471,   269,   132,   133,   943,   259,   269,   801,
     239,   432,   273,   265,   259,   260,   261,   262,    14,   264,
     259,   260,   261,   262,   445,   264,   447,   239,   271,   450,
      26,   271,   392,   393,   259,   259,   260,   261,   262,    35,
     264,   462,   271,   403,   404,   405,   898,   239,   239,   256,
     218,   276,   904,   265,   269,   210,   239,   239,   273,    25,
     563,   270,    58,   563,   194,   183,   196,   705,  1180,   429,
     214,   863,   235,   265,   434,  1187,   270,    43,   438,   271,
     271,   563,    78,    79,   151,   506,   507,   563,   271,   563,
     218,   270,   255,    74,   732,   239,   270,     7,    64,   217,
      10,    11,   232,   271,   259,   260,   261,   262,   104,   264,
     272,   471,    22,    23,   752,   815,   234,   184,   257,   228,
     758,   239,   256,   269,   191,   256,   764,   273,   919,   920,
     256,    68,    69,    70,   273,   131,   261,   262,   205,   264,
     207,   271,   260,   271,   421,   613,    44,    45,   425,   426,
     259,   260,   261,   262,   272,   264,    93,    67,    95,    96,
     520,   142,   522,   269,   160,   256,   526,   273,   611,   256,
     266,   152,   153,   269,   534,   256,    86,    87,   174,   175,
     176,   266,   256,   250,   269,   151,   182,   266,   256,   610,
     269,   187,   258,   259,   260,   261,   262,   835,   989,   837,
     256,   839,   256,   563,   256,   996,   997,   998,   999,   259,
     260,   261,   710,   256,   712,   259,   260,   261,   184,   919,
     920,   268,   269,   583,   862,   191,   272,   223,   270,   273,
    1167,   261,   262,   270,   655,   656,   657,   658,   219,   205,
     272,   207,   663,   270,   240,   666,   667,   243,   888,   889,
     246,   611,   612,   613,   235,   753,   121,   122,   618,  1150,
    1151,   682,   683,   684,   685,   686,   256,   248,   768,   768,
     256,   769,   210,   270,   255,   271,   259,   260,   261,   262,
     748,   264,   725,   266,   250,   256,   768,  1185,  1186,   989,
     711,   256,   768,   993,  1161,  1162,   996,   997,   998,   999,
     272,   270,   256,   785,   259,   260,   261,   262,    74,   785,
     256,   785,   270,   270,   735,     0,    82,  1215,   273,  1217,
    1218,   259,   260,   261,   262,   256,   264,    81,   688,   256,
     830,    97,   170,   271,   755,   756,   270,   270,   838,    39,
     761,   270,   763,   155,    77,   705,   272,   272,   830,   228,
     710,   740,   712,   257,   830,   272,   838,   855,   813,   814,
      38,   272,   838,   272,   257,   725,   259,   260,   261,   262,
     813,   264,   732,   272,   269,   272,   142,    55,   144,   877,
     259,   260,   261,   262,   272,   264,   152,   153,   748,   272,
     272,   272,   752,   753,   267,   268,   269,   270,   758,   759,
     760,   272,   762,   272,   764,   272,   888,   889,   768,   769,
      88,   272,   888,   889,   888,   889,   259,   260,   261,   262,
     272,   264,   272,   272,   784,   785,   259,   260,   261,   262,
     273,   264,   272,   111,   257,   272,   936,   797,   798,   273,
     273,   801,   269,   151,   259,   260,   261,   262,   214,   264,
     270,   265,   130,   813,   936,   256,   134,   256,   273,  1097,
     936,   256,   256,   884,   885,   256,   256,   256,   828,   235,
     830,   256,   150,   239,   256,   835,   897,   837,   838,   839,
     259,   260,   261,   262,   256,   264,   846,   256,   256,   255,
     273,   973,   256,   256,   273,   855,   256,   973,   256,   973,
     256,   128,   862,   863,   259,   260,   261,   262,   256,   264,
     256,   189,   190,   256,  1012,   256,  1014,   877,   273,   269,
     941,   135,   256,   201,   202,   256,   256,   256,   888,   889,
     985,   256,   256,   211,   955,   895,   256,   256,   959,  1177,
     900,   256,   902,   259,   260,   261,   262,   907,   264,   256,
     256,   256,   995,   256,   256,   256,   977,   273,   256,   256,
    1198,   256,  1200,    39,  1202,   241,   265,   988,  1066,   247,
     259,   260,   261,   262,   264,   264,   936,  1075,   269,   256,
     259,   260,   261,   262,   273,   264,   256,   256,   948,   949,
     950,   256,  1013,   256,   273,   170,  1058,   198,   270,  1020,
     222,    62,   256,   256,   259,   260,   261,   262,   221,   264,
     260,   443,   444,   973,   446,  1058,   448,   449,   273,   451,
     452,   453,   454,   455,   456,   457,   458,   459,   259,   260,
     261,   262,   256,   264,   273,   995,   259,   260,   261,   262,
     256,   264,   273,  1098,   256,   256,   256,  1007,   266,   269,
     273,   256,  1012,   256,  1014,   270,   256,   251,   256,  1080,
     259,   260,   261,   262,   256,   264,   256,  1027,  1159,  1090,
    1161,  1162,   256,   270,   273,   256,   259,   260,   261,   262,
       3,   264,   270,     6,   256,   256,    16,   256,   270,  1049,
     273,  1159,  1052,  1161,  1162,    18,   256,    20,  1058,   256,
     256,  1156,  1123,   270,   256,    28,  1066,    37,   270,   115,
    1131,  1132,   260,   266,   256,  1075,   270,    40,   259,   260,
     261,   262,   256,    53,  1179,   270,   259,  1087,   256,  1150,
    1151,   249,   273,   273,   269,   256,  1096,  1097,   117,   259,
     260,   261,   262,   256,   264,   256,   116,   256,   271,    72,
    1171,   258,  1173,   273,   256,   256,   270,    80,   168,    82,
      90,    91,   270,   259,   260,   261,   262,   270,   264,   256,
     259,   260,   261,   262,    97,   264,    99,   273,   108,   109,
     270,   256,   271,   256,  1205,   256,   256,   265,   272,   112,
     113,   114,   256,   256,   124,   125,   126,   256,    19,  1159,
     123,  1161,  1162,   271,   127,   270,   136,   269,   165,   265,
     140,   141,   270,   265,   256,   261,   276,  1177,     3,   276,
     256,     6,   266,   256,     9,   256,   156,   157,   158,   159,
     256,   256,   166,    18,   274,    20,    21,   274,  1198,    24,
    1200,   256,  1202,    28,   274,   256,   169,   259,   260,   261,
     262,   256,   264,   270,   256,    40,   179,   180,   265,   167,
     183,   273,   259,   260,   261,   262,   265,   264,   259,   260,
     261,   262,   195,   264,   197,   269,   273,   200,   270,   269,
      14,   269,   273,    68,    69,    70,   209,    72,   266,   266,
     213,   214,    26,   274,   256,    80,   271,   220,   269,   271,
      85,    35,   269,   269,   269,   269,   256,   237,    93,   269,
      95,    96,    51,  1102,    99,   143,   239,   386,   259,   260,
     261,   262,   107,   264,    58,   283,   610,   112,   113,   114,
     271,   386,   321,   218,   907,   225,   259,   260,   123,   392,
     533,   271,   127,   784,    78,    79,   124,   973,   863,   272,
     617,   274,   275,   877,   471,   342,     3,  1081,  1087,     6,
     752,   855,     9,   259,   260,   261,   262,  1173,   264,   236,
     104,    18,   263,    20,    21,   271,   534,    24,   716,   993,
    1052,    28,   577,  1072,   169,  1008,   259,   260,   261,   262,
     440,   264,   823,    40,   179,   180,   269,   131,   183,   259,
     260,   261,   262,  1004,   264,    -1,   828,    -1,   193,   269,
     195,    -1,   197,    -1,    -1,   200,    -1,   259,   260,   261,
     262,    -1,   264,    -1,   209,    72,   160,   269,   213,   259,
     260,   261,   262,    80,   264,   220,    -1,    -1,    85,   269,
     174,   175,   176,    -1,   229,    -1,    -1,    -1,   182,    -1,
      -1,    -1,    99,   187,   239,    -1,   259,   260,   261,   262,
     107,   264,    -1,   266,    -1,   112,   113,   114,   259,   260,
     261,   262,    -1,   264,   259,   260,   123,    -1,     3,    -1,
     127,     6,   258,   259,   260,   261,   262,   272,   264,   223,
      -1,    -1,    -1,    18,    -1,    20,   258,   259,   260,   261,
     262,    -1,   264,    28,    -1,    -1,   240,    -1,    -1,   243,
      -1,    -1,   246,    -1,    -1,    40,    -1,   258,   259,   260,
     261,   262,   169,   264,    -1,    -1,     3,    -1,    -1,     6,
      -1,    -1,   179,   180,    -1,    -1,   183,   259,   260,   261,
     262,    18,   264,    20,    -1,    -1,   193,    72,   195,    74,
     197,    28,    -1,   200,    -1,    80,   259,   260,   261,   262,
      -1,   264,   209,    40,    -1,    -1,   213,   259,   260,   261,
     262,    -1,   264,   220,    99,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,   229,    -1,    -1,    -1,    -1,   112,   113,   114,
      -1,    -1,   239,    -1,    -1,    72,    -1,    -1,   123,    -1,
      -1,    -1,   127,    80,    -1,    82,    -1,    -1,    -1,    -1,
      -1,    -1,   259,   260,    -1,    -1,    -1,   142,    -1,    -1,
      97,    -1,    99,    -1,    -1,   272,    -1,   152,   153,    -1,
      -1,    -1,    -1,    -1,    -1,   112,   113,   114,    -1,    -1,
      -1,    -1,    -1,    -1,   169,    -1,   123,    -1,     3,    -1,
     127,     6,    -1,    -1,   179,   180,    -1,    -1,   183,    -1,
      -1,    -1,    -1,    18,    -1,    20,    -1,    -1,    -1,    -1,
     195,    -1,   197,    28,    -1,   200,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,   209,    40,    -1,    -1,   213,    -1,
      -1,    -1,   169,    -1,    -1,   220,     3,    -1,    -1,     6,
      -1,    -1,   179,   180,    -1,    -1,   183,    -1,    -1,    -1,
     235,    18,    -1,    20,   239,    -1,    -1,    72,   195,    -1,
     197,    28,    -1,   200,    -1,    80,    -1,    82,    -1,    -1,
     255,    -1,   209,    40,   259,   260,   213,   214,    -1,    -1,
      -1,    -1,    97,   220,    99,    -1,    -1,   272,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,   112,   113,   114,
      -1,    -1,   239,    -1,    -1,    72,    -1,    -1,   123,    -1,
      -1,    -1,   127,    80,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,   259,   260,    -1,    -1,    -1,    -1,   265,    -1,
      -1,    -1,    99,   100,    -1,   272,    -1,    -1,     3,    -1,
      -1,     6,    -1,    -1,    -1,   112,   113,   114,    -1,    -1,
      -1,    -1,    -1,    18,   169,    20,   123,    -1,    -1,    -1,
     127,    -1,    -1,    28,   179,   180,    -1,    -1,   183,    -1,
      -1,    -1,    -1,    -1,    -1,    40,    -1,    -1,    -1,    -1,
     195,    -1,   197,    -1,    -1,   200,     3,    -1,    -1,     6,
      -1,    -1,    -1,    -1,   209,    -1,    -1,    -1,   213,   214,
      -1,    18,   169,    20,    -1,   220,    -1,    72,    -1,    -1,
      -1,    28,   179,   180,    -1,    80,   183,    -1,    -1,    -1,
      -1,    -1,    -1,    40,   239,    -1,    -1,    -1,   195,    -1,
     197,    -1,    -1,   200,    99,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,   209,    -1,   259,   260,   213,   112,   113,   114,
      -1,    -1,    -1,   220,    -1,    72,    -1,   272,   123,    -1,
      -1,    -1,   127,    80,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,   239,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    99,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,   259,   260,    -1,   112,   113,   114,    -1,    -1,
      -1,    -1,    -1,    -1,   169,   272,   123,    -1,    -1,    -1,
     127,    -1,    -1,    -1,   179,   180,     3,    -1,   183,     6,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
     195,    18,   197,    20,    -1,   200,    -1,    -1,    -1,    -1,
      -1,    28,    -1,    -1,   209,    -1,    -1,    -1,   213,    -1,
      -1,    -1,   169,    40,     3,   220,    -1,     6,    -1,    -1,
      16,    -1,   179,   180,    -1,    -1,   183,    -1,    -1,    18,
      -1,    20,    -1,    -1,   239,    -1,    -1,    -1,   195,    28,
     197,    37,    -1,   200,    -1,    72,    -1,    -1,    -1,    -1,
      -1,    40,   209,    80,   259,   260,   213,    53,    -1,    -1,
     265,    -1,    -1,   220,    -1,    -1,    -1,   272,    -1,    -1,
      -1,    -1,    99,    -1,    -1,    -1,    -1,    -1,    -1,   236,
      -1,    -1,   239,    72,    -1,   112,   113,   114,    -1,    -1,
      -1,    80,    -1,    -1,    90,    91,   123,    -1,    -1,    -1,
     127,    -1,   259,   260,    -1,    -1,    -1,    -1,    -1,    -1,
      99,    -1,   108,   109,    -1,   272,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,   112,   113,   114,    -1,    -1,   124,   125,
     126,    -1,    -1,    -1,   123,    -1,    -1,    -1,   127,    -1,
     136,    -1,   169,    -1,   140,   141,    -1,    -1,    -1,    -1,
      -1,    -1,   179,   180,    -1,    -1,   183,    -1,    -1,    -1,
     156,   157,   158,   159,    -1,    -1,    -1,    -1,   195,    -1,
     197,   410,   411,   200,   413,   414,   415,   416,   417,   418,
     169,   420,   209,   422,   423,   424,   213,    -1,    -1,    -1,
     179,   180,    -1,   220,   183,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,   195,    -1,   197,    -1,
      -1,   200,   239,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
     209,    -1,    -1,    -1,   213,    -1,    -1,    -1,    -1,    -1,
      -1,   220,   259,   260,     5,    -1,    -1,    -1,   265,    -1,
      -1,   237,    -1,    -1,    -1,   272,    -1,    -1,    -1,    -1,
     239,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    30,
      31,    32,    33,    34,    -1,    -1,    -1,    -1,    -1,    -1,
     259,   260,    -1,    -1,    -1,    46,    47,    48,    49,    50,
      51,    52,    -1,   272,    -1,    -1,    -1,    -1,    59,    60,
      61,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    75,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    83,    84,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    94,    -1,    -1,    -1,    98,    -1,    -1,
     101,    -1,    -1,    -1,   105,   106,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,   118,   119,   120,
     121,   122,    -1,    -1,    -1,    -1,    -1,    -1,   129,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,   147,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
     161,   162,   163,   164,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,   172,   173,    -1,    -1,    -1,   177,   178,    -1,    -1,
     181,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,   203,    -1,    -1,    -1,    -1,   208,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,   226,   227,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,   238,   239,    -1,
      -1,   242,    -1,    -1,   245,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,   253
};

  /* YYSTOS[STATE-NUM] -- The (internal number of the) accessing
     symbol of state STATE-NUM.  */
static const yytype_uint16 yystos[] =
{
       0,     5,    30,    31,    32,    33,    34,    46,    47,    48,
      49,    50,    51,    52,    59,    60,    61,    75,    83,    84,
      94,    98,   101,   105,   106,   118,   119,   120,   121,   122,
     129,   147,   161,   162,   163,   164,   172,   173,   177,   178,
     181,   203,   208,   226,   227,   238,   239,   242,   245,   253,
     278,   279,   280,   282,   296,   297,   298,   312,   313,   314,
     316,   321,   322,   323,   324,   325,   326,   331,   335,   338,
     339,   340,   341,   342,   343,   344,   356,   357,   358,   373,
     396,   401,   403,   404,   405,   411,   416,   417,   421,   435,
     436,   466,   471,   474,   482,   492,   496,   531,   556,   256,
     256,   256,   256,   256,   256,   282,   346,   347,   270,   270,
     282,   397,   282,   348,   360,   270,   270,   256,   256,   256,
     272,   272,   272,   256,   402,   256,   256,   256,   256,   256,
     256,   256,   256,   270,   270,   270,   406,   272,   272,   256,
     256,   270,   256,   272,   256,   256,   256,   270,   270,   270,
       0,   280,   256,    81,   256,    25,    43,    64,   151,   184,
     191,   205,   207,   250,   423,   170,    74,   142,   152,   153,
     235,   255,   288,   288,    82,    97,   214,   281,   282,   283,
     309,   310,   311,     3,     6,    18,    20,    28,    40,    72,
      80,    99,   112,   113,   114,   123,   127,   169,   179,   180,
     183,   195,   197,   200,   209,   213,   220,   259,   260,   272,
     282,   303,   304,   305,   307,   308,   283,   303,   337,   270,
     345,   346,   282,   353,   355,   374,   375,   384,   385,   386,
     270,   270,   359,   360,   282,   284,   483,   484,   309,   310,
     138,   139,   303,   282,   317,   317,   317,   309,   282,   405,
     411,   303,   303,   288,   303,   303,   303,   303,   219,   248,
     288,   284,   285,   493,   494,    14,    26,    35,    58,    78,
      79,   104,   131,   160,   174,   175,   176,   182,   187,   223,
     240,   243,   246,   327,   328,    39,   309,   318,   318,   100,
     303,   303,   155,   498,   303,   282,   298,   303,   303,   303,
     532,    77,   557,    16,    37,    53,    90,    91,   108,   109,
     124,   125,   126,   136,   140,   141,   156,   157,   158,   159,
     237,   332,   333,   265,   282,   299,   301,   303,   305,   306,
     308,   310,   272,   282,   300,   301,   302,   406,   406,   406,
     406,   406,   406,   272,   257,   272,   272,   272,   272,   272,
     272,   272,   272,   272,   272,   272,   272,   272,   272,   272,
     272,   272,   272,   303,   303,   303,    71,   143,   259,   260,
     261,   262,   264,   336,   258,   336,    56,    57,   349,   271,
     346,   267,   268,   269,   270,   289,   290,   291,   292,   293,
     271,   375,   259,   276,   387,   399,   361,   271,   360,   485,
     271,   484,   273,   269,   269,   265,   270,   271,   494,   256,
     256,   256,   256,   256,   256,   256,   256,   256,   256,   256,
     256,   256,   256,   256,   256,   256,   256,   271,   328,   256,
     319,   273,   256,   497,   269,   128,   534,   535,   256,   135,
     558,   559,   256,   256,   256,   256,   256,   256,   256,   256,
     256,   256,   256,   256,   256,   256,   256,   256,   256,   256,
     271,   333,   265,   294,   295,   303,   283,    39,   282,   404,
     405,   412,   413,   414,   418,   241,   438,   422,   438,   318,
     309,   303,   303,   303,   303,   303,   303,   303,   303,   303,
     303,   303,   303,   303,   303,   303,   303,   309,   303,   273,
     303,   303,   303,   303,   303,   303,   256,   256,    44,    45,
     350,   303,   291,   292,   386,   282,   362,   388,   260,   275,
     376,   377,   378,   379,   380,   381,   382,    54,   149,   188,
     230,   231,   271,   363,   265,   318,   318,   282,   495,   288,
     329,   329,   329,    27,   288,   330,   329,   329,   303,   329,
     329,   329,   329,   303,   329,   330,   329,   329,   329,   330,
     330,   286,   300,   269,   273,   303,   102,   212,   224,   499,
     500,   501,   502,   318,   256,    76,   135,   533,   536,   537,
     538,   539,   309,   256,   111,   559,   562,    66,    92,   252,
     334,   334,   334,   303,   334,   303,   334,   334,   303,   334,
     334,   334,   334,   334,   334,   334,   334,   334,   303,   266,
     269,   269,   256,   151,   408,   413,   198,   419,   270,   222,
     468,    38,    55,    88,   111,   130,   134,   150,   189,   190,
     201,   202,   211,   247,   424,   426,   427,   428,   429,   430,
     431,   432,    62,   441,   319,   273,   273,   273,   273,   273,
     273,   273,   273,   273,   273,   269,   269,   269,   269,   273,
     273,   273,   273,   258,   303,   303,   256,   256,   221,   351,
     271,   289,   376,   260,   274,   355,   274,   355,   145,   355,
     389,   390,   256,   256,   256,   256,   256,     4,    36,   132,
     133,   186,   233,   271,   364,   365,   366,   368,   369,   282,
     488,   319,   273,   266,   489,   269,   299,   306,   310,   320,
     256,   256,   256,    89,   199,   270,   503,   504,   505,   508,
     509,   319,    17,    29,   144,   256,   270,   271,   536,   282,
     560,   561,   256,   251,   563,   228,   295,   309,   315,   286,
     284,   194,   196,   232,   271,   407,   409,   410,   256,   424,
     286,   439,   440,   270,   467,   256,   256,   270,   256,   256,
     270,   256,   256,   256,   256,   270,   408,   426,   256,   270,
     437,   273,   303,   303,   303,   303,   303,   303,   303,   115,
     352,   376,   377,   377,   259,   265,   391,   392,   393,   303,
     236,   303,   400,   303,   303,   303,   355,   270,   270,   228,
     256,   367,   266,   273,   218,   271,   369,   490,   491,   286,
     300,   303,   300,   256,   256,   506,   271,   504,   273,   309,
     103,   137,   225,   540,   541,   542,   548,   552,   259,   286,
     256,   249,   564,   303,   273,   269,   406,   256,   256,   256,
      42,    63,   110,   185,   204,   206,   272,   284,   285,   420,
     408,   271,   439,   300,   469,   470,   408,   303,   303,   117,
     286,   355,   355,   433,   434,   303,   282,   398,   303,   286,
     116,   407,   282,   301,   306,   425,   300,   442,   443,   444,
     273,   273,   273,   273,   258,   256,   271,   390,   274,   275,
     299,   306,   310,   394,   395,   258,   383,   256,   282,   354,
     370,   372,   371,   372,    13,   355,   486,   256,   168,   480,
     144,   288,   309,   507,   288,    41,    73,   217,   234,   260,
     272,   304,   510,   511,   512,   513,   270,   270,   270,   271,
     541,   560,   282,   286,   287,   306,   256,   102,   212,   224,
     565,   212,   286,   415,   286,   287,   286,   420,   259,   260,
     261,   407,   271,   469,   407,   256,   286,   271,   434,   256,
     271,   443,    52,   192,   408,   445,   446,   461,   303,   303,
     395,   395,   266,   269,   282,   303,   289,   256,   271,   372,
     271,   372,   289,   270,   362,   256,   472,   514,   265,   272,
     515,   512,   512,   269,   271,   258,   259,   260,   261,   262,
     530,     8,   154,   171,   553,   554,   555,   544,   549,   550,
     555,   287,   256,   256,   256,   271,   303,    19,   481,   408,
     269,   273,   420,   420,   420,   303,   303,   270,   270,   271,
     394,   303,    68,    69,    70,    93,    95,    96,   450,   451,
     452,   453,   457,   458,   459,   460,   288,   165,   477,   265,
     303,   512,   265,   273,   511,   256,   259,   274,   516,   309,
     512,   512,   512,   512,   271,   554,   276,    13,   271,   281,
     282,   545,   546,   271,   550,   276,   300,   303,   300,   266,
     256,   475,   407,   303,   210,   271,   210,   462,   463,   488,
     447,   266,   256,   256,   451,   487,   256,   256,   256,   166,
     478,   281,   282,   518,   519,   520,   521,   522,   523,   524,
     525,   526,   266,   273,   520,   274,   274,   274,   283,   517,
      12,   300,   551,   256,   543,   543,    15,   300,   547,   303,
     444,   256,   256,   271,   463,   270,     9,    21,    24,    85,
     107,   193,   229,   303,   448,   449,   450,   454,   455,   456,
     265,   265,   489,   282,   286,   288,   256,   167,   479,   269,
     293,   269,   269,   266,   266,   274,   303,   476,   303,   303,
     464,   228,   271,   269,   455,   455,   271,   269,   288,   256,
     473,   254,   284,   285,   527,   527,   527,   408,   271,    38,
     450,   303,   456,   266,   266,   286,   288,   407,   269,   528,
     269,   528,   269,   528,   407,   256,   465,     7,    10,    11,
      22,    23,    67,    86,    87,   286,   529,   286,   286,   303,
     489,   269,   528,   528,   528,   271,   271
};

  /* YYR1[YYN] -- Symbol number of symbol that rule YYN derives.  */
static const yytype_uint16 yyr1[] =
{
       0,   277,   278,   279,   279,   280,   280,   280,   280,   280,
     280,   280,   280,   280,   280,   280,   280,   280,   280,   280,
     280,   280,   280,   280,   280,   280,   281,   282,   283,   284,
     285,   286,   287,   287,   288,   288,   288,   288,   288,   288,
     289,   289,   289,   289,   290,   290,   290,   290,   291,   292,
     293,   294,   294,   295,   295,   296,   297,   297,   297,   297,
     298,   299,   300,   300,   301,   302,   303,   303,   304,   304,
     305,   305,   306,   306,   307,   308,   308,   308,   308,   308,
     308,   308,   308,   308,   308,   308,   308,   308,   308,   308,
     308,   308,   308,   308,   308,   308,   308,   308,   308,   308,
     308,   308,   308,   308,   308,   309,   309,   310,   310,   310,
     310,   311,   312,   312,   312,   312,   312,   312,   312,   313,
     314,   315,   316,   317,   318,   319,   319,   320,   320,   320,
     321,   322,   323,   324,   325,   326,   327,   327,   328,   328,
     328,   328,   328,   328,   328,   328,   328,   328,   328,   328,
     328,   328,   328,   328,   328,   328,   328,   329,   330,   330,
     331,   332,   332,   333,   333,   333,   333,   333,   333,   333,
     333,   333,   333,   333,   333,   333,   333,   333,   333,   333,
     333,   334,   334,   334,   335,   335,   335,   335,   335,   336,
     336,   336,   337,   337,   337,   337,   338,   338,   338,   338,
     338,   338,   338,   338,   338,   338,   338,   338,   338,   338,
     338,   338,   338,   338,   339,   339,   339,   339,   340,   341,
     341,   341,   342,   342,   343,   344,   345,   345,   346,   347,
     348,   349,   349,   350,   350,   350,   351,   351,   352,   352,
     353,   354,   355,   356,   356,   357,   358,   359,   359,   361,
     360,   362,   363,   363,   364,   364,   365,   365,   365,   366,
     366,   366,   367,   367,   368,   369,   369,   370,   370,   371,
     371,   372,   373,   374,   374,   375,   376,   376,   377,   378,
     379,   380,   381,   382,   382,   382,   382,   383,   383,   384,
     385,   385,   386,   387,   387,   388,   389,   389,   390,   390,
     391,   391,   392,   393,   394,   394,   394,   395,   395,   395,
     396,   397,   398,   399,   399,   399,   399,   399,   399,   400,
     400,   402,   401,   403,   404,   404,   404,   404,   404,   405,
     406,   407,   408,   408,   409,   410,   410,   410,   411,   412,
     412,   413,   413,   415,   414,   416,   416,   418,   417,   419,
     419,   419,   419,   419,   419,   419,   419,   420,   420,   420,
     420,   420,   422,   421,   423,   423,   423,   423,   423,   424,
     424,   425,   426,   426,   426,   426,   426,   426,   426,   426,
     426,   426,   427,   427,   428,   428,   428,   428,   429,   429,
     430,   431,   432,   432,   433,   433,   434,   435,   437,   436,
     438,   439,   440,   440,   441,   442,   442,   443,   444,   444,
     445,   445,   447,   446,   448,   448,   449,   449,   449,   449,
     449,   449,   449,   450,   450,   451,   451,   451,   451,   452,
     453,   454,   455,   455,   456,   456,   456,   457,   458,   458,
     459,   460,   460,   461,   462,   462,   464,   465,   463,   467,
     466,   468,   469,   470,   470,   472,   473,   471,   475,   476,
     474,   477,   477,   478,   478,   479,   479,   480,   480,   481,
     481,   482,   483,   483,   485,   486,   487,   484,   488,   489,
     489,   490,   490,   491,   492,   493,   493,   495,   494,   497,
     496,   498,   498,   499,   499,   499,   500,   501,   502,   503,
     503,   504,   504,   504,   506,   505,   507,   507,   507,   508,
     509,   510,   510,   511,   512,   512,   512,   512,   512,   512,
     512,   512,   512,   514,   513,   513,   515,   513,   516,   516,
     516,   516,   516,   517,   518,   519,   520,   520,   520,   520,
     520,   520,   521,   522,   523,   524,   525,   526,   527,   527,
     527,   528,   528,   529,   529,   529,   529,   529,   529,   529,
     529,   530,   530,   532,   531,   533,   533,   534,   534,   535,
     535,   535,   536,   536,   537,   538,   539,   540,   540,   541,
     541,   541,   542,   543,   543,   544,   544,   545,   545,   546,
     546,   547,   547,   548,   549,   549,   550,   551,   551,   552,
     553,   553,   554,   555,   555,   555,   556,   557,   558,   558,
     559,   560,   561,   561,   562,   563,   564,   565,   565,   565,
     565
};

  /* YYR2[YYN] -- Number of symbols on the right hand side of rule YYN.  */
static const yytype_uint8 yyr2[] =
{
       0,     2,     1,     1,     2,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       4,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       0,     1,     1,     1,     1,     1,     2,     2,     1,     1,
       3,     1,     3,     1,     7,     3,     3,     3,     3,     3,
       1,     1,     1,     1,     3,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     3,     4,     4,     4,     6,
       6,     6,     4,     4,     4,     4,     4,     4,     4,     4,
       4,     4,     6,     1,     1,     1,     1,     4,     3,     3,
       3,     3,     3,     2,     2,     1,     1,     1,     1,     3,
       5,     1,     1,     1,     1,     1,     1,     1,     1,     7,
       1,     1,     4,     1,     1,     0,     3,     1,     1,     1,
       5,     7,     7,     4,     6,     4,     1,     2,     3,     3,
       3,     3,     3,     3,     3,     3,     3,     3,     3,     3,
       3,     3,     3,     3,     3,     3,     3,     1,     1,     1,
       4,     1,     2,     3,     3,     3,     3,     3,     3,     3,
       3,     3,     3,     3,     3,     3,     3,     3,     3,     3,
       3,     1,     1,     1,     3,     3,     4,     3,     4,     0,
       1,     1,     1,     3,     5,     7,     3,     3,     3,     3,
       3,     3,     3,     3,     3,     3,     3,     3,     3,     3,
       3,     3,     3,     3,     3,     3,     3,     3,     3,     1,
       1,     1,     1,     1,     2,     4,     1,     2,     7,     1,
       1,     3,     3,     0,     3,     3,     0,     1,     0,     3,
       1,     2,     2,     1,     1,     2,     4,     1,     2,     0,
       5,     1,     0,     2,     1,     1,     3,     4,     4,     1,
       1,     1,     1,     1,     1,     4,     4,     1,     2,     1,
       2,     3,     4,     1,     2,     1,     1,     2,     2,     2,
       2,     3,     3,     1,     1,     1,     1,     0,     2,     6,
       1,     3,     1,     0,     2,     2,     1,     3,     1,     1,
       1,     1,     3,     5,     1,     2,     2,     1,     1,     1,
       5,     1,     1,     0,     4,     4,     4,     4,     4,     1,
       1,     0,     3,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     0,     2,     1,     3,     3,     5,     6,     1,
       2,     1,     1,     0,     7,     1,     1,     0,     8,     3,
       3,     3,     3,     3,     3,     3,     3,     1,     3,     3,
       3,     3,     0,     7,     1,     1,     1,     1,     1,     1,
       2,     1,     3,     3,     1,     3,     3,     3,     3,     3,
       3,     4,     1,     1,     1,     1,     1,     1,     3,     6,
       9,    12,     3,     3,     1,     2,     2,     1,     0,     9,
       4,     1,     1,     2,     4,     1,     2,     1,     0,     2,
       1,     1,     0,     5,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     2,     1,     1,     1,     1,     5,
       5,     1,     1,     3,     1,     3,     1,     3,     1,     1,
       5,     1,     1,     4,     1,     2,     0,     0,     7,     0,
       8,     4,     1,     1,     2,     0,     0,    14,     0,     0,
      14,     0,     3,     0,     3,     0,     3,     0,     3,     0,
       3,     4,     1,     2,     0,     0,     0,    11,     1,     0,
       2,     1,     1,     3,     4,     1,     2,     0,     5,     0,
       7,     0,     3,     1,     1,     1,     3,     3,     3,     1,
       2,     1,     1,     1,     0,     6,     1,     1,     1,     3,
       3,     1,     3,     2,     1,     1,     3,     3,     3,     3,
       3,     2,     4,     0,     5,     4,     0,     5,     1,     2,
       2,     3,     2,     1,     1,     2,     1,     1,     1,     1,
       1,     1,     4,     4,     4,     6,     6,     6,     1,     1,
       1,     0,     2,     1,     1,     1,     1,     1,     1,     1,
       1,     0,     2,     0,     6,     1,     2,     0,     1,     3,
       3,     3,     1,     1,     1,     3,     4,     1,     2,     1,
       1,     1,     4,     2,     0,     2,     0,     2,     2,     1,
       1,     1,     1,     4,     1,     2,     3,     1,     1,     4,
       1,     2,     3,     1,     1,     1,     9,     3,     1,     2,
       3,     1,     1,     3,     3,     3,     3,     0,     3,     3,
       3
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
        case 29:
#line 642 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.sym) = mdl_existing_object(parse_state, (yyvsp[0].str))); }
#line 3155 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 30:
#line 645 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.sym) = mdl_existing_region(parse_state, (yyvsp[-3].sym), (yyvsp[-1].str))); }
#line 3161 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 31:
#line 648 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.vec3) = mdl_point(parse_state, &(yyvsp[0].nlist))); }
#line 3167 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 33:
#line 652 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.vec3) = mdl_point_scalar((yyvsp[0].dbl))); }
#line 3173 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 34:
#line 655 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = 1; }
#line 3179 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 35:
#line 656 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = 0; }
#line 3185 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 36:
#line 657 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = 1; }
#line 3191 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 37:
#line 658 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = 0; }
#line 3197 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 38:
#line 659 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = 1; }
#line 3203 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 39:
#line 660 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = 0; }
#line 3209 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 40:
#line 663 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.mol_type).orient_set = 0; }
#line 3215 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 43:
#line 666 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.mol_type).orient_set = 1; (yyval.mol_type).orient = 0; }
#line 3221 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 44:
#line 670 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.mol_type).orient = 1; (yyval.mol_type).orient_set = 1; }
#line 3227 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 45:
#line 671 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.mol_type).orient = -1; (yyval.mol_type).orient_set = 1; }
#line 3233 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 46:
#line 672 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                          (yyval.mol_type) = (yyvsp[-1].mol_type);
                                                          if ((yyval.mol_type).orient >= 32767)
                                                          {
                                                            /* Seriously?  Wow. */
                                                            mdlerror(parse_state, "molecule orientation must not be greater than 32767");
                                                            return 1;
                                                          }
                                                          ++ (yyval.mol_type).orient;
                                                      }
#line 3248 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 47:
#line 682 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                          (yyval.mol_type) = (yyvsp[-1].mol_type);
                                                          if ((yyval.mol_type).orient <= -32768)
                                                          {
                                                            /* Seriously?  Wow. */
                                                            mdlerror(parse_state, "molecule orientation must not be less than -32768");
                                                            return 1;
                                                          }
                                                          -- (yyval.mol_type).orient;
                                                      }
#line 3263 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 50:
#line 700 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                          (yyval.mol_type).orient = (int) (yyvsp[-1].dbl);
                                                          (yyval.mol_type).orient_set = 1;
                                                          if ((yyval.mol_type).orient != (yyvsp[-1].dbl))
                                                          {
                                                            mdlerror(parse_state, "molecule orientation specified inside braces must be an integer between -32768 and 32767.");
                                                            return 1;
                                                          }
                                                      }
#line 3277 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 52:
#line 713 "../src/../src/mdlparse.y" /* yacc.c:1646  */
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
#line 3293 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 53:
#line 726 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mcell_generate_range_singleton(&(yyval.nlist), (yyvsp[0].dbl))); }
#line 3299 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 54:
#line 727 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_generate_range(parse_state, &(yyval.nlist), (yyvsp[-5].dbl), (yyvsp[-3].dbl), (yyvsp[-1].dbl))); }
#line 3305 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 55:
#line 733 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                          char *include_path = mcell_find_include_file((yyvsp[0].str), parse_state->vol->curr_file);
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
#line 3327 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 56:
#line 756 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_assign_variable_double(parse_state, (yyvsp[-2].sym), (yyvsp[0].dbl))); }
#line 3333 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 57:
#line 757 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_assign_variable_string(parse_state, (yyvsp[-2].sym), (yyvsp[0].str))); }
#line 3339 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 58:
#line 758 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_assign_variable(parse_state, (yyvsp[-2].sym), (yyvsp[0].sym))); }
#line 3345 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 59:
#line 759 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_assign_variable_array(parse_state, (yyvsp[-2].sym), (yyvsp[0].nlist).value_head)); }
#line 3351 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 60:
#line 762 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.sym) = mdl_get_or_create_variable(parse_state, (yyvsp[0].str))); }
#line 3357 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 61:
#line 765 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.sym) = mdl_existing_variable(parse_state, (yyvsp[0].str))); }
#line 3363 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 63:
#line 769 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                          struct num_expr_list *elp;
                                                          (yyval.nlist).value_head = (struct num_expr_list *) (yyvsp[0].sym)->value;
                                                          (yyval.nlist).value_count = 1;
                                                          for (elp = (yyval.nlist).value_head; elp->next != NULL; elp = elp->next)
                                                            ++ (yyval.nlist).value_count;
                                                          (yyval.nlist).value_tail = elp;
                                                          (yyval.nlist).shared = 1;
                                                      }
#line 3377 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 64:
#line 780 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { mdl_debug_dump_array((yyvsp[-1].nlist).value_head); (yyval.nlist) = (yyvsp[-1].nlist); }
#line 3383 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 65:
#line 783 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.sym) = mdl_existing_array(parse_state, (yyvsp[0].str))); }
#line 3389 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 69:
#line 791 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.dbl) = *(double *) (yyvsp[0].sym)->value; }
#line 3395 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 70:
#line 794 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.dbl) = (yyvsp[0].llival); }
#line 3401 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 74:
#line 802 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.sym) = mdl_existing_double(parse_state, (yyvsp[0].str))); }
#line 3407 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 75:
#line 806 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.dbl) = (yyvsp[-1].dbl); }
#line 3413 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 76:
#line 807 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKF((yyval.dbl) = exp((yyvsp[-1].dbl))); }
#line 3419 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 77:
#line 808 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_expr_log(parse_state, (yyvsp[-1].dbl), &(yyval.dbl))); }
#line 3425 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 78:
#line 809 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_expr_log10(parse_state, (yyvsp[-1].dbl), &(yyval.dbl))); }
#line 3431 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 79:
#line 810 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.dbl) = max2d((yyvsp[-3].dbl), (yyvsp[-1].dbl)); }
#line 3437 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 80:
#line 811 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.dbl) = min2d((yyvsp[-3].dbl), (yyvsp[-1].dbl)); }
#line 3443 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 81:
#line 812 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.dbl) = mdl_expr_roundoff((yyvsp[-1].dbl), (int) (yyvsp[-3].dbl)); }
#line 3449 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 82:
#line 813 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.dbl) = floor((yyvsp[-1].dbl)); }
#line 3455 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 83:
#line 814 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.dbl) = ceil((yyvsp[-1].dbl)); }
#line 3461 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 84:
#line 815 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.dbl) = sin((yyvsp[-1].dbl)); }
#line 3467 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 85:
#line 816 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.dbl) = cos((yyvsp[-1].dbl)); }
#line 3473 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 86:
#line 817 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKF((yyval.dbl) = tan((yyvsp[-1].dbl))); }
#line 3479 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 87:
#line 818 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKF((yyval.dbl) = asin((yyvsp[-1].dbl))); }
#line 3485 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 88:
#line 819 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKF((yyval.dbl) = acos((yyvsp[-1].dbl))); }
#line 3491 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 89:
#line 820 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.dbl) = atan((yyvsp[-1].dbl)); }
#line 3497 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 90:
#line 821 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKF((yyval.dbl) = sqrt((yyvsp[-1].dbl))); }
#line 3503 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 91:
#line 822 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.dbl) = fabs((yyvsp[-1].dbl)); }
#line 3509 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 92:
#line 823 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_expr_mod(parse_state, (yyvsp[-3].dbl), (yyvsp[-1].dbl), &(yyval.dbl))); }
#line 3515 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 93:
#line 824 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.dbl) = MY_PI; }
#line 3521 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 94:
#line 825 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.dbl) = mdl_expr_rng_uniform(parse_state); }
#line 3527 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 95:
#line 826 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.dbl) = rng_gauss(parse_state->vol->rng); }
#line 3533 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 96:
#line 827 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.dbl) = parse_state->vol->seed_seq; }
#line 3539 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 97:
#line 828 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_expr_string_to_double(parse_state, (yyvsp[-1].str), &(yyval.dbl))); }
#line 3545 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 98:
#line 829 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKF((yyval.dbl) = (yyvsp[-2].dbl) + (yyvsp[0].dbl)); }
#line 3551 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 99:
#line 830 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKF((yyval.dbl) = (yyvsp[-2].dbl) - (yyvsp[0].dbl)); }
#line 3557 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 100:
#line 831 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKF((yyval.dbl) = (yyvsp[-2].dbl) * (yyvsp[0].dbl)); }
#line 3563 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 101:
#line 832 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_expr_div(parse_state, (yyvsp[-2].dbl), (yyvsp[0].dbl), &(yyval.dbl))); }
#line 3569 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 102:
#line 833 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_expr_pow(parse_state, (yyvsp[-2].dbl), (yyvsp[0].dbl), &(yyval.dbl))); }
#line 3575 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 103:
#line 834 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.dbl) = -(yyvsp[0].dbl); }
#line 3581 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 104:
#line 835 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.dbl) = (yyvsp[0].dbl); }
#line 3587 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 106:
#line 840 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.str) = mdl_strdup((char const *) (yyvsp[0].sym)->value)); }
#line 3593 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 107:
#line 844 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.str) = mdl_strip_quotes((yyvsp[0].str))); }
#line 3599 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 108:
#line 845 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.str) = mdl_strdup(parse_state->vol->mdl_infile_name)); }
#line 3605 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 109:
#line 846 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.str) = mdl_strcat((yyvsp[-2].str), (yyvsp[0].str))); }
#line 3611 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 110:
#line 847 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.str) = mdl_string_format(parse_state, (yyvsp[-2].str), (yyvsp[-1].printfargs).arg_head)); }
#line 3617 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 111:
#line 850 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.sym) = mdl_existing_string(parse_state, (yyvsp[0].str))); }
#line 3623 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 119:
#line 866 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_fopen(parse_state, (yyvsp[-6].sym), (yyvsp[-3].str), (yyvsp[-1].str))); }
#line 3629 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 120:
#line 869 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.sym) = mdl_new_filehandle(parse_state, (yyvsp[0].str))); }
#line 3635 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 121:
#line 872 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.str) = (yyvsp[0].str); CHECK(mdl_valid_file_mode(parse_state, (yyvsp[0].str))); }
#line 3641 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 122:
#line 875 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_fclose(parse_state, (yyvsp[-1].sym))); }
#line 3647 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 123:
#line 878 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.sym) = mdl_existing_file_stream(parse_state, (yyvsp[0].str))); }
#line 3653 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 124:
#line 881 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.str) = mdl_expand_string_escapes((yyvsp[0].str))); }
#line 3659 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 125:
#line 884 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.printfargs).arg_head = (yyval.printfargs).arg_tail = NULL; }
#line 3665 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 126:
#line 885 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                        (yyval.printfargs) = (yyvsp[-2].printfargs);
                                                        if ((yyval.printfargs).arg_tail)
                                                          (yyval.printfargs).arg_tail = (yyval.printfargs).arg_tail->next = (yyvsp[0].printfarg);
                                                        else
                                                          (yyval.printfargs).arg_tail = (yyval.printfargs).arg_head = (yyvsp[0].printfarg);
                                                        (yyvsp[0].printfarg)->next = NULL;
                                                      }
#line 3678 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 127:
#line 895 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.printfarg) = mdl_new_printf_arg_double((yyvsp[0].dbl))); }
#line 3684 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 128:
#line 896 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.printfarg) = mdl_new_printf_arg_string((yyvsp[0].str))); }
#line 3690 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 129:
#line 897 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                          switch ((yyvsp[0].sym)->sym_type)
                                                          {
                                                            case DBL: CHECKN((yyval.printfarg) = mdl_new_printf_arg_double(*(double *) (yyvsp[0].sym)->value)); break;
                                                            case STR: CHECKN((yyval.printfarg) = mdl_new_printf_arg_string((char *) (yyvsp[0].sym)->value)); break;
                                                            default:
                                                              mdlerror(parse_state, "invalid variable type referenced");
                                                              return 1;
                                                          }
                                                      }
#line 3705 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 130:
#line 909 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_printf(parse_state, (yyvsp[-2].str), (yyvsp[-1].printfargs).arg_head)); }
#line 3711 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 131:
#line 915 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_fprintf(parse_state, (struct file_stream *) (yyvsp[-4].sym)->value, (yyvsp[-2].str), (yyvsp[-1].printfargs).arg_head)); }
#line 3717 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 132:
#line 921 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_sprintf(parse_state, (yyvsp[-4].sym), (yyvsp[-2].str), (yyvsp[-1].printfargs).arg_head)); }
#line 3723 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 133:
#line 924 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { mdl_print_time(parse_state, (yyvsp[-1].str)); }
#line 3729 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 134:
#line 930 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_fprint_time(parse_state, (yyvsp[-3].sym), (yyvsp[-1].str))); }
#line 3735 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 138:
#line 946 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { if (!parse_state->vol->quiet_flag) mdl_set_all_notifications(parse_state->vol, (yyvsp[0].tok)); }
#line 3741 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 139:
#line 947 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { if (!parse_state->vol->quiet_flag) parse_state->vol->notify->progress_report        = (yyvsp[0].tok); }
#line 3747 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 140:
#line 948 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { if (!parse_state->vol->quiet_flag) parse_state->vol->notify->diffusion_constants    = (yyvsp[0].tok); }
#line 3753 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 141:
#line 949 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { if (!parse_state->vol->quiet_flag) parse_state->vol->notify->reaction_probabilities = (yyvsp[0].tok); }
#line 3759 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 142:
#line 950 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { if (!parse_state->vol->quiet_flag) parse_state->vol->notify->time_varying_reactions = (yyvsp[0].tok); }
#line 3765 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 143:
#line 951 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { if (!parse_state->vol->quiet_flag) parse_state->vol->notify->reaction_prob_notify   = (yyvsp[0].dbl); }
#line 3771 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 144:
#line 952 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { if (!parse_state->vol->quiet_flag) parse_state->vol->notify->partition_location     = (yyvsp[0].tok); }
#line 3777 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 145:
#line 953 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { if (!parse_state->vol->quiet_flag) parse_state->vol->notify->box_triangulation      = (yyvsp[0].tok); }
#line 3783 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 146:
#line 954 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { if (!parse_state->vol->quiet_flag) parse_state->vol->notify->release_events         = (yyvsp[0].tok); }
#line 3789 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 147:
#line 955 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { if (!parse_state->vol->quiet_flag) parse_state->vol->notify->file_writes            = (yyvsp[0].tok); }
#line 3795 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 148:
#line 956 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { if (!parse_state->vol->quiet_flag) parse_state->vol->notify->final_summary          = (yyvsp[0].tok); }
#line 3801 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 149:
#line 957 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { if (!parse_state->vol->quiet_flag) parse_state->vol->notify->throughput_report      = (yyvsp[0].tok); }
#line 3807 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 150:
#line 958 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { if (!parse_state->vol->quiet_flag) parse_state->vol->notify->reaction_output_report = (yyvsp[0].tok); }
#line 3813 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 151:
#line 959 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { if (!parse_state->vol->quiet_flag) parse_state->vol->notify->volume_output_report   = (yyvsp[0].tok); }
#line 3819 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 152:
#line 960 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { if (!parse_state->vol->quiet_flag) parse_state->vol->notify->viz_output_report      = (yyvsp[0].tok); }
#line 3825 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 153:
#line 961 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { if (!parse_state->vol->quiet_flag) parse_state->vol->notify->checkpoint_report      = (yyvsp[0].tok); }
#line 3831 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 154:
#line 962 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                          if (!parse_state->vol->quiet_flag && parse_state->vol->log_freq == ULONG_MAX)
                                                            parse_state->vol->notify->iteration_report = (yyvsp[0].tok);
                                                      }
#line 3840 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 155:
#line 966 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { if (!parse_state->vol->quiet_flag) CHECK(mdl_set_iteration_report_freq(parse_state, (long long) (yyvsp[0].dbl))); }
#line 3846 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 156:
#line 967 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { if (!parse_state->vol->quiet_flag) parse_state->vol->notify->molecule_collision_report    = (yyvsp[0].tok); }
#line 3852 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 157:
#line 971 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = ((yyvsp[0].tok) ? NOTIFY_FULL : NOTIFY_NONE); }
#line 3858 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 158:
#line 975 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = ((yyvsp[0].tok) ? NOTIFY_FULL : NOTIFY_NONE); }
#line 3864 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 159:
#line 976 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = NOTIFY_BRIEF; }
#line 3870 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 163:
#line 992 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { mdl_set_all_warnings(parse_state->vol, (byte) (yyvsp[0].tok)); }
#line 3876 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 164:
#line 993 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->vol->notify->neg_diffusion = (byte)(yyvsp[0].tok); }
#line 3882 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 165:
#line 994 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->vol->notify->neg_reaction = (byte)(yyvsp[0].tok); }
#line 3888 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 166:
#line 995 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->vol->notify->high_reaction_prob = (byte)(yyvsp[0].tok); }
#line 3894 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 167:
#line 996 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->vol->notify->reaction_prob_warn = (yyvsp[0].dbl); }
#line 3900 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 168:
#line 997 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->vol->notify->close_partitions = (byte)(yyvsp[0].tok); }
#line 3906 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 169:
#line 998 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->vol->notify->degenerate_polys = (byte)(yyvsp[0].tok); }
#line 3912 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 170:
#line 999 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->vol->notify->overwritten_file = (byte)(yyvsp[0].tok); }
#line 3918 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 171:
#line 1000 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->vol->notify->short_lifetime = (byte)(yyvsp[0].tok); }
#line 3924 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 172:
#line 1001 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_lifetime_warning_threshold(parse_state, (long long) (yyvsp[0].dbl))); }
#line 3930 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 173:
#line 1002 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->vol->notify->missed_reactions = (byte)(yyvsp[0].tok); }
#line 3936 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 174:
#line 1003 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_missed_reaction_warning_threshold(parse_state, (yyvsp[0].dbl))); }
#line 3942 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 175:
#line 1004 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->vol->notify->missed_surf_orient = (byte)(yyvsp[0].tok); }
#line 3948 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 176:
#line 1005 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->vol->notify->useless_vol_orient = (byte)(yyvsp[0].tok); }
#line 3954 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 177:
#line 1006 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->vol->notify->mol_placement_failure = (byte) (yyvsp[0].tok); }
#line 3960 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 178:
#line 1007 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->vol->notify->invalid_output_step_time = (byte) (yyvsp[0].tok); }
#line 3966 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 179:
#line 1008 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->vol->notify->large_molecular_displacement = (byte) (yyvsp[0].tok); }
#line 3972 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 180:
#line 1009 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->vol->notify->add_remove_mesh_warning = (byte) (yyvsp[0].tok); }
#line 3978 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 181:
#line 1013 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = WARN_COPE;  }
#line 3984 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 182:
#line 1014 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = WARN_WARN;  }
#line 3990 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 183:
#line 1015 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = WARN_ERROR; }
#line 3996 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 184:
#line 1021 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_checkpoint_infile(parse_state, (yyvsp[0].str))); }
#line 4002 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 185:
#line 1022 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_checkpoint_outfile(parse_state, (yyvsp[0].str))); }
#line 4008 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 186:
#line 1023 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_checkpoint_interval(parse_state, (yyvsp[-1].dbl), (yyvsp[0].tok))); }
#line 4014 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 187:
#line 1024 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_keep_checkpoint_files(parse_state, (yyvsp[0].tok))); }
#line 4020 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 188:
#line 1026 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_realtime_checkpoint(parse_state, (long) (yyvsp[-1].dbl), (yyvsp[0].tok))); }
#line 4026 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 189:
#line 1029 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = 0; }
#line 4032 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 190:
#line 1030 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = 1; }
#line 4038 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 191:
#line 1031 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = 0; }
#line 4044 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 192:
#line 1035 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { /* seconds */     (yyval.dbl) = (yyvsp[0].dbl); }
#line 4050 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 193:
#line 1036 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { /* mm:ss */       (yyval.dbl) = (yyvsp[-2].dbl) * 60 + (yyvsp[0].dbl); }
#line 4056 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 194:
#line 1037 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { /* hh:mm:ss */    (yyval.dbl) = (yyvsp[-4].dbl) * 3600 + (yyvsp[-2].dbl) * 60 + (yyvsp[0].dbl); }
#line 4062 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 195:
#line 1039 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { /* dd:hh:mm:ss */ (yyval.dbl) = (yyvsp[-6].dbl) * 86400 + (yyvsp[-4].dbl) * 3600 + (yyvsp[-2].dbl) * 60 + (yyvsp[0].dbl); }
#line 4068 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 196:
#line 1046 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_time_step(parse_state, (yyvsp[0].dbl))); }
#line 4074 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 197:
#line 1047 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_space_step(parse_state, (yyvsp[0].dbl))); }
#line 4080 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 198:
#line 1048 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_max_time_step(parse_state, (yyvsp[0].dbl))); }
#line 4086 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 199:
#line 1049 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_num_iterations(parse_state, (long long) (yyvsp[0].dbl))); }
#line 4092 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 200:
#line 1050 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->vol->randomize_smol_pos = !((yyvsp[0].tok)); }
#line 4098 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 201:
#line 1051 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->vol->use_expanded_list = (yyvsp[0].tok); }
#line 4104 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 202:
#line 1052 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->vol->vacancy_search_dist2 = max2d((yyvsp[0].dbl), 0.0); }
#line 4110 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 203:
#line 1053 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_num_radial_directions(parse_state, (int) (yyvsp[0].dbl))); }
#line 4116 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 204:
#line 1054 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->vol->fully_random = 1; }
#line 4122 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 205:
#line 1055 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_num_radial_subdivisions(parse_state, (int) (yyvsp[0].dbl))); }
#line 4128 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 206:
#line 1056 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_grid_density(parse_state, (yyvsp[0].dbl))); }
#line 4134 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 207:
#line 1057 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_interaction_radius(parse_state, (yyvsp[0].dbl))); }
#line 4140 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 208:
#line 1058 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->vol->surface_reversibility=(yyvsp[0].tok); parse_state->vol->volume_reversibility=(yyvsp[0].tok); }
#line 4146 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 209:
#line 1059 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->vol->surface_reversibility=1;  parse_state->vol->volume_reversibility=0;  }
#line 4152 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 210:
#line 1060 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->vol->surface_reversibility=0;  parse_state->vol->volume_reversibility=1;  }
#line 4158 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 211:
#line 1061 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mcell_add_dynamic_geometry_file((yyvsp[0].str), parse_state)); }
#line 4164 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 212:
#line 1062 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->vol->dynamic_geometry_molecule_placement = 0; }
#line 4170 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 213:
#line 1063 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->vol->dynamic_geometry_molecule_placement = 1; }
#line 4176 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 214:
#line 1070 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->vol->mem_part_x = (int) (yyvsp[0].dbl); }
#line 4182 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 215:
#line 1071 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->vol->mem_part_y = (int) (yyvsp[0].dbl); }
#line 4188 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 216:
#line 1072 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->vol->mem_part_z = (int) (yyvsp[0].dbl); }
#line 4194 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 217:
#line 1073 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->vol->mem_part_pool = (int) (yyvsp[0].dbl); }
#line 4200 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 218:
#line 1077 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mcell_set_partition(parse_state->vol, (yyvsp[-2].tok), & (yyvsp[0].nlist))); }
#line 4206 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 219:
#line 1081 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = X_PARTS; }
#line 4212 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 220:
#line 1082 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = Y_PARTS; }
#line 4218 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 221:
#line 1083 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = Z_PARTS; }
#line 4224 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 224:
#line 1094 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { mdl_print_species_summary(parse_state->vol, (yyvsp[0].mcell_mol_spec)); }
#line 4230 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 225:
#line 1098 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { mdl_print_species_summaries(parse_state->vol, (yyvsp[-1].mcell_species_lst).species_head); }
#line 4236 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 226:
#line 1102 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.mcell_species_lst).species_count = 0; CHECK(mdl_add_to_species_list(&(yyval.mcell_species_lst), (yyvsp[0].mcell_mol_spec))); }
#line 4242 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 227:
#line 1103 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.mcell_species_lst) = (yyvsp[-1].mcell_species_lst); CHECK(mdl_add_to_species_list(&(yyval.mcell_species_lst), (yyvsp[0].mcell_mol_spec))); }
#line 4248 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 228:
#line 1112 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.mcell_mol_spec) = mdl_create_species(parse_state, (yyvsp[-6].str), (yyvsp[-4].diff_const).D, (yyvsp[-4].diff_const).is_2d, (yyvsp[-3].dbl), (yyvsp[-2].ival), (yyvsp[-1].dbl) )); }
#line 4254 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 230:
#line 1118 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.sym) = mdl_new_mol_species(parse_state, (yyvsp[0].str))); }
#line 4260 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 231:
#line 1122 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.diff_const).is_2d = 0; (yyval.diff_const).D = (yyvsp[0].dbl); CHECK(mdl_check_diffusion_constant(parse_state, & (yyval.diff_const).D)); }
#line 4266 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 232:
#line 1123 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.diff_const).is_2d = 1; (yyval.diff_const).D = (yyvsp[0].dbl); CHECK(mdl_check_diffusion_constant(parse_state, & (yyval.diff_const).D)); }
#line 4272 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 233:
#line 1127 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.dbl) = 0.0; }
#line 4278 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 234:
#line 1128 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                          if ((yyvsp[0].dbl) <= 0)
                                                          {
                                                            mdlerror_fmt(parse_state, "Requested custom time step of %.15g; custom time step must be positive.", (yyvsp[0].dbl));
                                                            return 1;
                                                          }

                                                          (yyval.dbl) = (yyvsp[0].dbl);
                                                      }
#line 4292 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 235:
#line 1137 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                          if ((yyvsp[0].dbl) <= 0)
                                                          {
                                                            mdlerror_fmt(parse_state, "Requested custom space step of %.15g; custom space step must be positive.", (yyvsp[0].dbl));
                                                            return 1;
                                                          }

                                                          (yyval.dbl) = -(yyvsp[0].dbl);
                                                      }
#line 4306 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 236:
#line 1148 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.ival) = 0; }
#line 4312 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 237:
#line 1149 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.ival) = 1; }
#line 4318 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 238:
#line 1153 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.dbl) = 0; }
#line 4324 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 239:
#line 1154 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                        if ((yyvsp[0].dbl) <= 0)
                                                        {
                                                          mdlerror_fmt(parse_state, "Requested maximum step length of %.15g; maximum step length must be positive.", (yyvsp[0].dbl));
                                                          return 1;
                                                        }
                                                        (yyval.dbl) = (yyvsp[0].dbl);
                                                      }
#line 4337 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 240:
#line 1164 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.sym) = mdl_existing_molecule(parse_state, (yyvsp[0].str))); }
#line 4343 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 241:
#line 1168 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.mol_type) = (yyvsp[0].mol_type); CHECKN((yyval.mol_type).mol_type = mdl_existing_surface_molecule(parse_state, (yyvsp[-1].str))); }
#line 4349 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 242:
#line 1172 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                        (yyval.mol_type) = (yyvsp[0].mol_type);
                                                        if (! (yyval.mol_type).orient_set)
                                                          (yyval.mol_type).orient = 0;
                                                        (yyval.mol_type).mol_type = (yyvsp[-1].sym);
                                                      }
#line 4360 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 249:
#line 1207 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { mdl_start_surface_class(parse_state, (yyvsp[-1].sym)); }
#line 4366 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 250:
#line 1209 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { mdl_finish_surface_class(parse_state); }
#line 4372 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 251:
#line 1212 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.sym) = mdl_existing_surface_class(parse_state, (yyvsp[0].str))); }
#line 4378 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 256:
#line 1229 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN(mdl_assemble_surface_reaction(parse_state, (yyvsp[-2].tok), parse_state->current_surface_class, (yyvsp[0].mol_type).mol_type, (yyvsp[0].mol_type).orient)); }
#line 4384 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 257:
#line 1232 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
              struct sym_entry *mol_sym = retrieve_sym("ALL_MOLECULES", parse_state->vol->mol_sym_table);
              if(!(yyvsp[0].mol_type).orient_set) (yyvsp[0].mol_type).orient = 0;
              CHECKN(mdl_assemble_surface_reaction(parse_state, (yyvsp[-3].tok), parse_state->current_surface_class, mol_sym, (yyvsp[0].mol_type).orient));}
#line 4393 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 258:
#line 1238 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN(mdl_assemble_concentration_clamp_reaction(parse_state, parse_state->current_surface_class, (yyvsp[-2].mol_type).mol_type, (yyvsp[-2].mol_type).orient, (yyvsp[0].dbl))); }
#line 4399 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 259:
#line 1241 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = RFLCT; }
#line 4405 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 260:
#line 1242 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = TRANSP; }
#line 4411 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 261:
#line 1243 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = SINK; }
#line 4417 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 264:
#line 1250 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->current_surface_class->sm_dat_head = (yyvsp[0].surf_mol_dat_list).sm_head; }
#line 4423 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 265:
#line 1257 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.surf_mol_dat_list) = (yyvsp[-1].surf_mol_dat_list); }
#line 4429 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 266:
#line 1261 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.surf_mol_dat_list) = (yyvsp[-1].surf_mol_dat_list); }
#line 4435 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 267:
#line 1265 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                          (yyvsp[0].surf_mol_dat)->quantity_type = SURFMOLDENS;
                                                          (yyval.surf_mol_dat_list).sm_tail = (yyval.surf_mol_dat_list).sm_head = (yyvsp[0].surf_mol_dat);
                                                      }
#line 4444 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 268:
#line 1270 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                          (yyval.surf_mol_dat_list) = (yyvsp[-1].surf_mol_dat_list);
                                                          (yyvsp[0].surf_mol_dat)->quantity_type = SURFMOLDENS;
                                                          (yyval.surf_mol_dat_list).sm_tail = (yyval.surf_mol_dat_list).sm_tail->next = (yyvsp[0].surf_mol_dat);
                                                      }
#line 4454 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 269:
#line 1278 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                          (yyvsp[0].surf_mol_dat)->quantity_type = SURFMOLNUM;
                                                          (yyval.surf_mol_dat_list).sm_tail = (yyval.surf_mol_dat_list).sm_head = (yyvsp[0].surf_mol_dat);
                                                      }
#line 4463 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 270:
#line 1283 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                          (yyval.surf_mol_dat_list) = (yyvsp[-1].surf_mol_dat_list);
                                                          (yyvsp[0].surf_mol_dat)->quantity_type = SURFMOLNUM;
                                                          (yyval.surf_mol_dat_list).sm_tail = (yyval.surf_mol_dat_list).sm_tail->next = (yyvsp[0].surf_mol_dat);
                                                      }
#line 4473 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 271:
#line 1291 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.surf_mol_dat) = mdl_new_surf_mol_data(parse_state, &(yyvsp[-2].mol_type), (yyvsp[0].dbl))); }
#line 4479 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 281:
#line 1320 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.react_arrow).catalyst = (yyvsp[-1].mol_type); (yyval.react_arrow).flags = ARROW_CATALYTIC; }
#line 4485 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 282:
#line 1325 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.react_arrow).catalyst = (yyvsp[-1].mol_type); (yyval.react_arrow).flags = ARROW_CATALYTIC | ARROW_BIDIRECTIONAL; }
#line 4491 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 283:
#line 1330 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.react_arrow).catalyst.mol_type = NULL; (yyval.react_arrow).flags = 0; }
#line 4497 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 285:
#line 1332 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.react_arrow).catalyst.mol_type = NULL; (yyval.react_arrow).flags = ARROW_BIDIRECTIONAL; }
#line 4503 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 287:
#line 1336 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.sym) = NULL; }
#line 4509 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 288:
#line 1337 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.sym) = mdl_new_rxn_pathname(parse_state, (yyvsp[0].str))); }
#line 4515 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 289:
#line 1343 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN(mdl_assemble_reaction(parse_state, (yyvsp[-5].mol_type_list).mol_type_head, &(yyvsp[-4].mol_type), &(yyvsp[-3].react_arrow), (yyvsp[-2].mol_type_list).mol_type_head, &(yyvsp[-1].react_rates), (yyvsp[0].sym))); }
#line 4521 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 290:
#line 1346 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_reaction_player_singleton(parse_state, & (yyval.mol_type_list), & (yyvsp[0].mol_type))); }
#line 4527 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 291:
#line 1347 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.mol_type_list) = (yyvsp[-2].mol_type_list); CHECK(mdl_add_reaction_player(parse_state, & (yyval.mol_type_list), & (yyvsp[0].mol_type))); }
#line 4533 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 293:
#line 1354 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.mol_type).mol_type = NULL; }
#line 4539 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 294:
#line 1355 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.mol_type) = (yyvsp[0].mol_type); }
#line 4545 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 295:
#line 1359 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.mol_type) = (yyvsp[0].mol_type); (yyval.mol_type).mol_type = (yyvsp[-1].sym); }
#line 4551 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 296:
#line 1362 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_reaction_player_singleton(parse_state, & (yyval.mol_type_list), & (yyvsp[0].mol_type))); }
#line 4557 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 297:
#line 1363 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.mol_type_list) = (yyvsp[-2].mol_type_list); CHECK(mdl_add_reaction_player(parse_state, & (yyval.mol_type_list), & (yyvsp[0].mol_type))); }
#line 4563 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 298:
#line 1366 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.mol_type).mol_type = NULL; (yyval.mol_type).orient_set = 0; }
#line 4569 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 302:
#line 1375 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                        if ((yyvsp[-1].react_rates).forward_rate.rate_type == RATE_UNSET)
                                                        {
                                                          mdlerror(parse_state, "invalid reaction rate specification: must specify a forward rate.");
                                                          return 1;
                                                        }

                                                        (yyval.react_rates) = (yyvsp[-1].react_rates);
                                                      }
#line 4583 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 303:
#line 1386 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                        if (((yyvsp[-3].react_rates).forward_rate.rate_type  != RATE_UNSET && (yyvsp[-1].react_rates).forward_rate.rate_type  != RATE_UNSET)  ||
                                                            ((yyvsp[-3].react_rates).backward_rate.rate_type != RATE_UNSET && (yyvsp[-1].react_rates).backward_rate.rate_type != RATE_UNSET))
                                                        {
                                                          mdlerror_fmt(parse_state, "when two reaction rates are specified, one must be a forward rate, and one must be a reverse rate");
                                                          return 1;
                                                        }

                                                        (yyval.react_rates) = (yyvsp[-3].react_rates);
                                                        if ((yyvsp[-1].react_rates).forward_rate.rate_type != RATE_UNSET)
                                                          (yyval.react_rates).forward_rate = (yyvsp[-1].react_rates).forward_rate;
                                                        else
                                                          (yyval.react_rates).backward_rate = (yyvsp[-1].react_rates).backward_rate;
                                                      }
#line 4602 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 304:
#line 1403 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.react_rates).forward_rate = (yyvsp[0].react_rate); (yyval.react_rates).backward_rate.rate_type = RATE_UNSET; CHECK(mdl_valid_rate(parse_state, &(yyvsp[0].react_rate))); }
#line 4608 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 305:
#line 1404 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.react_rates).forward_rate = (yyvsp[0].react_rate); (yyval.react_rates).backward_rate.rate_type = RATE_UNSET; CHECK(mdl_valid_rate(parse_state, &(yyvsp[0].react_rate))); }
#line 4614 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 306:
#line 1405 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.react_rates).backward_rate = (yyvsp[0].react_rate); (yyval.react_rates).forward_rate.rate_type = RATE_UNSET; CHECK(mdl_valid_rate(parse_state, &(yyvsp[0].react_rate))); }
#line 4620 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 307:
#line 1409 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.react_rate).rate_type = RATE_CONSTANT; (yyval.react_rate).v.rate_constant = (yyvsp[0].dbl); }
#line 4626 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 308:
#line 1410 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.react_rate).rate_type = RATE_FILE; (yyval.react_rate).v.rate_file = (yyvsp[0].str); }
#line 4632 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 309:
#line 1411 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_reaction_rate_from_var(parse_state, & (yyval.react_rate), (yyvsp[0].sym))); }
#line 4638 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 310:
#line 1422 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_release_pattern(parse_state, (yyvsp[-3].sym), &(yyvsp[-1].rpat))); }
#line 4644 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 311:
#line 1425 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.sym) = mdl_new_release_pattern(parse_state, (yyvsp[0].str))); }
#line 4650 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 312:
#line 1428 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.sym) = mdl_existing_release_pattern_or_rxn_pathname(parse_state, (yyvsp[0].str))); }
#line 4656 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 313:
#line 1432 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                        (yyval.rpat).delay = 0;
                                                        (yyval.rpat).release_interval = FOREVER;
                                                        (yyval.rpat).train_interval = FOREVER;
                                                        (yyval.rpat).train_duration = FOREVER;
                                                        (yyval.rpat).number_of_trains = 1;
                                                      }
#line 4668 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 314:
#line 1440 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.rpat) = (yyvsp[-3].rpat); (yyval.rpat).delay = (yyvsp[0].dbl) / parse_state->vol->time_unit; }
#line 4674 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 315:
#line 1442 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.rpat) = (yyvsp[-3].rpat); (yyval.rpat).release_interval = (yyvsp[0].dbl) / parse_state->vol->time_unit; }
#line 4680 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 316:
#line 1444 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.rpat) = (yyvsp[-3].rpat); (yyval.rpat).train_interval = (yyvsp[0].dbl) / parse_state->vol->time_unit; }
#line 4686 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 317:
#line 1446 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.rpat) = (yyvsp[-3].rpat); (yyval.rpat).train_duration = (yyvsp[0].dbl) / parse_state->vol->time_unit; }
#line 4692 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 318:
#line 1448 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.rpat) = (yyvsp[-3].rpat); (yyval.rpat).number_of_trains = (yyvsp[0].ival); }
#line 4698 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 319:
#line 1451 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.ival) = (int) (yyvsp[0].dbl); }
#line 4704 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 320:
#line 1452 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.ival) = INT_MAX; }
#line 4710 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 321:
#line 1459 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->current_object = parse_state->vol->root_instance; }
#line 4716 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 322:
#line 1460 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                        check_regions(parse_state->vol->root_instance, (yyvsp[0].obj));
                                                        add_child_objects(parse_state->vol->root_instance, (yyvsp[0].obj), (yyvsp[0].obj));
                                                        parse_state->current_object = parse_state->vol->root_object;
                                                      }
#line 4726 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 323:
#line 1470 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { add_child_objects(parse_state->vol->root_object, (yyvsp[0].obj), (yyvsp[0].obj)); }
#line 4732 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 329:
#line 1486 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.sym) = mdl_start_object(parse_state, (yyvsp[0].str))); }
#line 4738 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 331:
#line 1492 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { mdl_finish_object(parse_state); }
#line 4744 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 335:
#line 1505 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { transform_translate(parse_state->vol, parse_state->current_object->t_matrix, (yyvsp[0].vec3)); }
#line 4750 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 336:
#line 1506 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { transform_scale(parse_state->current_object->t_matrix, (yyvsp[0].vec3)); }
#line 4756 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 337:
#line 1507 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_transform_rotate(parse_state, parse_state->current_object->t_matrix, (yyvsp[-2].vec3), (yyvsp[0].dbl))); }
#line 4762 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 338:
#line 1516 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                          struct object *the_object = (struct object *) (yyvsp[-5].sym)->value;
                                                          the_object->object_type = META_OBJ;
                                                          add_child_objects(the_object, (yyvsp[-2].obj_list).obj_head, (yyvsp[-2].obj_list).obj_tail);
                                                          (yyval.obj) = the_object;
                                                      }
#line 4773 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 339:
#line 1525 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { mdl_object_list_singleton(& (yyval.obj_list), (yyvsp[0].obj)); }
#line 4779 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 340:
#line 1526 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.obj_list) = (yyvsp[-1].obj_list); mdl_add_object_to_list(& (yyval.obj_list), (yyvsp[0].obj)); }
#line 4785 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 343:
#line 1535 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_deep_copy_object(parse_state, (struct object *) (yyvsp[-3].sym)->value, (struct object *) (yyvsp[-1].sym)->value)); }
#line 4791 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 344:
#line 1537 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.obj) = (struct object *) (yyvsp[-6].sym)->value; }
#line 4797 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 347:
#line 1547 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_start_release_site(parse_state, (yyvsp[-2].sym), SHAPE_UNDEFINED)); }
#line 4803 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 348:
#line 1551 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.obj) = mdl_finish_release_site(parse_state, (yyvsp[-7].sym))); }
#line 4809 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 349:
#line 1554 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_release_site_geometry_region(parse_state, parse_state->current_release_site, parse_state->current_object, (yyvsp[0].rev))); }
#line 4815 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 350:
#line 1555 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_release_site_geometry_object(parse_state, parse_state->current_release_site, (struct object *) (yyvsp[0].sym)->value)); }
#line 4821 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 351:
#line 1556 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->current_release_site->release_shape = SHAPE_SPHERICAL; }
#line 4827 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 352:
#line 1557 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->current_release_site->release_shape = SHAPE_CUBIC; }
#line 4833 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 353:
#line 1558 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->current_release_site->release_shape = SHAPE_ELLIPTIC; }
#line 4839 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 354:
#line 1559 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->current_release_site->release_shape = SHAPE_RECTANGULAR; }
#line 4845 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 355:
#line 1560 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->current_release_site->release_shape = SHAPE_SPHERICAL_SHELL; }
#line 4851 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 356:
#line 1561 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                          parse_state->current_release_site->release_shape = SHAPE_LIST;
                                                          parse_state->current_release_site->release_number_method = CONSTNUM;
                                                      }
#line 4860 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 357:
#line 1568 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.rev) = new_release_region_expr_term((yyvsp[0].sym))); }
#line 4866 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 358:
#line 1569 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.rev) = (yyvsp[-1].rev); }
#line 4872 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 359:
#line 1570 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.rev) = new_release_region_expr_binary((yyvsp[-2].rev), (yyvsp[0].rev), REXP_UNION)); }
#line 4878 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 360:
#line 1571 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.rev) = new_release_region_expr_binary((yyvsp[-2].rev), (yyvsp[0].rev), REXP_SUBTRACTION)); }
#line 4884 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 361:
#line 1572 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.rev) = new_release_region_expr_binary((yyvsp[-2].rev), (yyvsp[0].rev), REXP_INTERSECTION)); }
#line 4890 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 362:
#line 1577 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_start_release_site(parse_state, (yyvsp[-2].sym), (yyvsp[-1].tok))); }
#line 4896 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 363:
#line 1580 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.obj) = mdl_finish_release_site(parse_state, (yyvsp[-6].sym))); }
#line 4902 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 364:
#line 1583 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = SHAPE_SPHERICAL; }
#line 4908 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 365:
#line 1584 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = SHAPE_CUBIC; }
#line 4914 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 366:
#line 1585 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = SHAPE_ELLIPTIC; }
#line 4920 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 367:
#line 1586 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = SHAPE_RECTANGULAR; }
#line 4926 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 368:
#line 1587 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = SHAPE_SPHERICAL_SHELL; }
#line 4932 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 371:
#line 1595 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.sym) = mdl_existing_num_or_array(parse_state, (yyvsp[0].str))); }
#line 4938 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 372:
#line 1599 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { set_release_site_location(parse_state->vol, parse_state->current_release_site, (yyvsp[0].vec3)); }
#line 4944 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 373:
#line 1600 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_release_site_molecule(parse_state, parse_state->current_release_site, & (yyvsp[0].mol_type))); }
#line 4950 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 374:
#line 1601 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                        if (parse_state->current_release_site->release_shape == SHAPE_LIST)
                                                        {
                                                          mdlerror(parse_state, "molecules are already specified in a list--cannot set number or density.");
                                                          return 1;
                                                        }
                                                      }
#line 4962 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 375:
#line 1608 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_release_site_diameter(parse_state, parse_state->current_release_site, (yyvsp[0].dbl) * (((yyvsp[-2].tok) == SITE_RADIUS) ? 2.0 : 1.0))); }
#line 4968 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 376:
#line 1609 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_release_site_diameter_array(parse_state, parse_state->current_release_site, (yyvsp[0].nlist).value_count, (yyvsp[0].nlist).value_head, ((yyvsp[-2].tok) == SITE_RADIUS) ? 2.0 : 1.0)); }
#line 4974 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 377:
#line 1610 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_release_site_diameter_var(parse_state, parse_state->current_release_site, ((yyvsp[-2].tok) == SITE_RADIUS) ? 2.0 : 1.0, (yyvsp[0].sym))); }
#line 4980 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 378:
#line 1611 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_release_site_periodic_box(parse_state, parse_state->current_release_site, (yyvsp[0].vec3))); }
#line 4986 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 379:
#line 1612 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_release_site_probability(parse_state, parse_state->current_release_site, (yyvsp[0].dbl))); }
#line 4992 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 380:
#line 1614 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_release_site_pattern(parse_state, parse_state->current_release_site, (yyvsp[0].sym))); }
#line 4998 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 381:
#line 1616 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_release_site_molecule_positions(parse_state, parse_state->current_release_site, & (yyvsp[-1].rsm_list))); }
#line 5004 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 382:
#line 1620 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = SITE_DIAMETER; }
#line 5010 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 383:
#line 1621 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = SITE_RADIUS; }
#line 5016 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 388:
#line 1633 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { set_release_site_constant_number(parse_state->current_release_site, (yyvsp[0].dbl)); }
#line 5022 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 389:
#line 1636 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { set_release_site_constant_number(parse_state->current_release_site, (yyvsp[-1].dbl)); }
#line 5028 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 390:
#line 1643 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { set_release_site_gaussian_number(parse_state->current_release_site, (yyvsp[-4].dbl), (yyvsp[-1].dbl)); }
#line 5034 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 391:
#line 1651 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { set_release_site_volume_dependent_number(parse_state->current_release_site, (yyvsp[-7].dbl), (yyvsp[-4].dbl), (yyvsp[-1].dbl)); }
#line 5040 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 392:
#line 1655 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_release_site_concentration(parse_state, parse_state->current_release_site, (yyvsp[0].dbl))); }
#line 5046 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 393:
#line 1656 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(set_release_site_density(parse_state->current_release_site, (yyvsp[0].dbl))); }
#line 5052 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 394:
#line 1660 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { release_single_molecule_singleton(& (yyval.rsm_list), (yyvsp[0].rsm)); }
#line 5058 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 395:
#line 1662 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.rsm_list) = (yyvsp[-1].rsm_list); add_release_single_molecule_to_list(& (yyval.rsm_list), (yyvsp[0].rsm)); }
#line 5064 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 396:
#line 1666 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.rsm) = mdl_new_release_single_molecule(parse_state, &(yyvsp[-1].mol_type), (yyvsp[0].vec3))); }
#line 5070 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 398:
#line 1677 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                        CHECKN((yyval.obj) = mdl_new_polygon_list(
                                                          parse_state, (yyvsp[-4].str), (yyvsp[-1].vertlist).vertex_count, (yyvsp[-1].vertlist).vertex_head,
                                                          (yyvsp[0].ecl).connection_count, (yyvsp[0].ecl).connection_head));
                                                      }
#line 5080 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 399:
#line 1686 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                          (yyval.obj) = (struct object *) (yyvsp[-3].obj);
                                                          CHECK(mdl_finish_polygon_list(parse_state, (yyval.obj)));
                                                      }
#line 5089 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 400:
#line 1692 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.vertlist) = (yyvsp[-1].vertlist); }
#line 5095 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 401:
#line 1695 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.vertlistitem) = mdl_new_vertex_list_item((yyvsp[0].vec3))); }
#line 5101 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 402:
#line 1698 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { mdl_vertex_list_singleton(& (yyval.vertlist), (yyvsp[0].vertlistitem)); }
#line 5107 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 403:
#line 1699 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.vertlist) = (yyvsp[-1].vertlist); mdl_add_vertex_to_list(& (yyval.vertlist), (yyvsp[0].vertlistitem)); }
#line 5113 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 404:
#line 1704 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.ecl) = (yyvsp[-1].ecl); }
#line 5119 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 405:
#line 1708 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { mdl_element_connection_list_singleton(& (yyval.ecl), (yyvsp[0].elem_conn)); }
#line 5125 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 406:
#line 1710 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.ecl) = (yyvsp[-1].ecl); mdl_add_element_connection_to_list(& (yyval.ecl), (yyvsp[0].elem_conn)); }
#line 5131 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 407:
#line 1713 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.elem_conn) = mdl_new_element_connection(parse_state, & (yyvsp[0].nlist))); }
#line 5137 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 412:
#line 1729 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN(parse_state->current_region = mdl_get_region(parse_state, parse_state->current_object, "REMOVED")); }
#line 5143 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 413:
#line 1731 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                          parse_state->current_region->element_list_head = (yyvsp[-1].elem_list).elml_head;
                                                          if (parse_state->current_object->object_type == POLY_OBJ)
                                                          {
                                                            CHECK(mdl_normalize_elements(parse_state, parse_state->current_region,0));
                                                          }
                                                      }
#line 5155 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 416:
#line 1745 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = Z_POS; }
#line 5161 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 417:
#line 1746 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = Z_NEG; }
#line 5167 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 418:
#line 1747 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = Y_NEG; }
#line 5173 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 419:
#line 1748 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = Y_POS; }
#line 5179 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 420:
#line 1749 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = X_NEG; }
#line 5185 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 421:
#line 1750 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = X_POS; }
#line 5191 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 422:
#line 1751 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = ALL_SIDES; }
#line 5197 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 424:
#line 1757 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.elem_list) = (yyvsp[-1].elem_list); mdl_add_elements_to_list(& (yyval.elem_list), (yyvsp[0].elem_list).elml_head, (yyvsp[0].elem_list).elml_tail); }
#line 5203 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 427:
#line 1763 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.elem_list).elml_tail = (yyval.elem_list).elml_head = (yyvsp[0].elem_list_item); }
#line 5209 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 428:
#line 1764 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.elem_list).elml_tail = (yyval.elem_list).elml_head = (yyvsp[0].elem_list_item); }
#line 5215 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 429:
#line 1769 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.elem_list) = (yyvsp[-1].elem_list); }
#line 5221 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 430:
#line 1774 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.elem_list) = (yyvsp[-1].elem_list); mdl_set_elements_to_exclude((yyval.elem_list).elml_head); }
#line 5227 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 432:
#line 1781 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.elem_list).elml_tail = (yyval.elem_list).elml_head = (yyvsp[0].elem_list_item); }
#line 5233 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 433:
#line 1782 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.elem_list) = (yyvsp[-2].elem_list); mdl_add_elements_to_list(& (yyval.elem_list), (yyvsp[0].elem_list_item), (yyvsp[0].elem_list_item)); }
#line 5239 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 434:
#line 1785 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.elem_list_item) = new_element_list((unsigned int) (yyvsp[0].dbl), (unsigned int) (yyvsp[0].dbl))); }
#line 5245 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 435:
#line 1786 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.elem_list_item) = new_element_list((unsigned int) (yyvsp[-2].dbl), (unsigned int) (yyvsp[0].dbl))); }
#line 5251 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 436:
#line 1787 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.elem_list_item) = mdl_new_element_side(parse_state, (yyvsp[0].tok))); }
#line 5257 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 437:
#line 1790 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.elem_list_item) = mdl_new_element_previous_region(parse_state, parse_state->current_object, parse_state->current_region, (yyvsp[0].str), (yyvsp[-2].tok))); }
#line 5263 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 438:
#line 1793 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = 0; }
#line 5269 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 439:
#line 1794 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = 1; }
#line 5275 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 440:
#line 1797 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.elem_list_item) = mdl_new_element_patch(parse_state, parse_state->current_polygon, (yyvsp[-2].vec3), (yyvsp[0].vec3), (yyvsp[-4].tok))); }
#line 5281 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 441:
#line 1800 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = 0; }
#line 5287 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 442:
#line 1801 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = 1; }
#line 5293 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 446:
#line 1817 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->current_region = (yyvsp[-1].reg); }
#line 5299 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 447:
#line 1818 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_region_elements(parse_state, (yyvsp[-3].reg), (yyvsp[0].elem_list).elml_head, (yyvsp[-3].reg)->parent->object_type == POLY_OBJ)); }
#line 5305 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 448:
#line 1820 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->current_region = NULL; }
#line 5311 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 449:
#line 1828 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                        CHECKN(mdl_new_voxel_list(parse_state, (yyvsp[-4].sym),
                                                                                  (yyvsp[-1].vertlist).vertex_count, (yyvsp[-1].vertlist).vertex_head,
                                                                                  (yyvsp[0].ecl).connection_count, (yyvsp[0].ecl).connection_head));
                                                      }
#line 5321 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 450:
#line 1834 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.obj) = (struct object *) (yyvsp[-7].sym)->value; }
#line 5327 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 451:
#line 1839 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.ecl) = (yyvsp[-1].ecl); }
#line 5333 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 452:
#line 1842 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.elem_conn) = mdl_new_tet_element_connection(parse_state, & (yyvsp[0].nlist))); }
#line 5339 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 453:
#line 1846 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                          (yyval.ecl).connection_head = (yyval.ecl).connection_tail = (yyvsp[0].elem_conn);
                                                          (yyval.ecl).connection_count = 1;
                                                      }
#line 5348 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 454:
#line 1850 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                          (yyval.ecl) = (yyvsp[-1].ecl);
                                                          (yyval.ecl).connection_tail = (yyval.ecl).connection_tail->next = (yyvsp[0].elem_conn);
                                                          ++ (yyval.ecl).connection_count;
                                                      }
#line 5358 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 455:
#line 1860 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->vol->periodic_traditional = (yyvsp[0].tok); }
#line 5364 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 456:
#line 1863 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN(mdl_create_periodic_box(parse_state, (yyvsp[-7].vec3), (yyvsp[-5].vec3), (yyvsp[-2].tok), (yyvsp[-1].tok), (yyvsp[0].tok))); }
#line 5370 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 457:
#line 1864 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_finish_periodic_box(parse_state)); }
#line 5376 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 458:
#line 1870 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN(mdl_new_box_object(parse_state, (yyvsp[-8].sym), (yyvsp[-3].vec3), (yyvsp[-1].vec3))); }
#line 5382 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 459:
#line 1871 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_triangulate_box_object(parse_state, (yyvsp[-10].sym), parse_state->current_polygon, (yyvsp[-2].dbl))); }
#line 5388 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 460:
#line 1873 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                          CHECK(mdl_finish_box_object(parse_state, (yyvsp[-13].sym)));
                                                          (yyval.obj) = (struct object *) (yyvsp[-13].sym)->value;
                                                      }
#line 5397 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 461:
#line 1880 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = 0; }
#line 5403 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 462:
#line 1881 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = (yyvsp[0].tok); }
#line 5409 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 463:
#line 1885 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = 0; }
#line 5415 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 464:
#line 1886 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = (yyvsp[0].tok); }
#line 5421 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 465:
#line 1890 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = 0; }
#line 5427 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 466:
#line 1891 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = (yyvsp[0].tok); }
#line 5433 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 467:
#line 1895 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = 0; }
#line 5439 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 468:
#line 1896 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = (yyvsp[0].tok); }
#line 5445 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 469:
#line 1899 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.dbl) = 0.0; }
#line 5451 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 470:
#line 1900 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                        (yyval.dbl) = (yyvsp[0].dbl);
                                                        if ((yyval.dbl) < 2.0)
                                                        {
                                                          mdlerror(parse_state, "invalid aspect ratio requested (must be greater than or equal to 2.0)");
                                                          return 1;
                                                        }
                                                      }
#line 5464 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 474:
#line 1926 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_start_existing_obj_region_def(parse_state, (yyvsp[0].sym))); }
#line 5470 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 475:
#line 1927 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->current_region = (yyvsp[-1].reg); }
#line 5476 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 476:
#line 1929 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { mdl_set_region_elements(parse_state, (yyvsp[-4].reg), (yyvsp[0].elem_list).elml_head, 1); }
#line 5482 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 477:
#line 1931 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                          parse_state->current_region = NULL;
                                                          parse_state->current_polygon = NULL;
                                                          parse_state->current_object = parse_state->vol->root_object;
                                                      }
#line 5492 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 478:
#line 1938 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.reg) = mdl_create_region(parse_state, parse_state->current_object, (yyvsp[0].str))); }
#line 5498 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 482:
#line 1949 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { mdl_add_surf_mol_to_region(parse_state->current_region, & (yyvsp[0].surf_mol_dat_list)); }
#line 5504 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 483:
#line 1953 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { mdl_set_region_surface_class(parse_state, parse_state->current_region, (yyvsp[0].sym)); }
#line 5510 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 487:
#line 1972 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->current_region = (struct region *) (yyvsp[-1].sym)->value; }
#line 5516 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 488:
#line 1974 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->current_region = NULL; }
#line 5522 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 489:
#line 1982 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                          parse_state->header_comment = NULL;  /* No header by default */
                                                          parse_state->exact_time_flag = 1;    /* Print exact_time column in TRIGGER output by default */
                                                      }
#line 5531 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 490:
#line 1988 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_add_reaction_output_block_to_world(parse_state, (int) (yyvsp[-4].dbl), & (yyvsp[-2].ro_otimes), & (yyvsp[-1].ro_sets))); }
#line 5537 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 491:
#line 1992 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.dbl) = COUNTBUFFERSIZE; }
#line 5543 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 492:
#line 1993 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                          double temp_value = (yyvsp[0].dbl);
                                                          if (!(temp_value >= 1.0 && temp_value < UINT_MAX))
                                                          {
                                                            mdlerror_fmt(parse_state, "Requested buffer size of %.15g lines is invalid.  Suggested range is 100-1000000.", temp_value);
                                                            return 1;
                                                          }
                                                          (yyval.dbl) = (yyvsp[0].dbl);
                                                      }
#line 5557 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 496:
#line 2009 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.ro_otimes).type = OUTPUT_BY_STEP; (yyval.ro_otimes).step = (yyvsp[0].dbl); }
#line 5563 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 497:
#line 2013 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                        (yyval.ro_otimes).type = OUTPUT_BY_ITERATION_LIST;
                                                        (yyval.ro_otimes).values = (yyvsp[0].nlist);
                                                      }
#line 5572 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 498:
#line 2021 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                        (yyval.ro_otimes).type = OUTPUT_BY_TIME_LIST;
                                                        (yyval.ro_otimes).values = (yyvsp[0].nlist);
                                                      }
#line 5581 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 499:
#line 2028 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.ro_sets).set_head = (yyval.ro_sets).set_tail = (yyvsp[0].ro_set); }
#line 5587 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 500:
#line 2030 "../src/../src/mdlparse.y" /* yacc.c:1646  */
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
#line 5602 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 502:
#line 2044 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.ro_set) = NULL; }
#line 5608 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 503:
#line 2045 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.ro_set) = NULL; }
#line 5614 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 504:
#line 2049 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {  parse_state->count_flags = 0; }
#line 5620 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 505:
#line 2051 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.ro_set) = mdl_populate_output_set(parse_state, parse_state->header_comment, parse_state->exact_time_flag, (yyvsp[-3].ro_cols).column_head, (yyvsp[-1].tok), (yyvsp[0].str))); }
#line 5626 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 506:
#line 2055 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.str) = NULL; }
#line 5632 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 507:
#line 2056 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.str) = ((yyvsp[0].tok) ? "" : NULL); }
#line 5638 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 508:
#line 2057 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.str) = (yyvsp[0].str); }
#line 5644 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 509:
#line 2061 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->header_comment = (yyvsp[0].str); }
#line 5650 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 510:
#line 2065 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->exact_time_flag = (yyvsp[0].tok); }
#line 5656 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 512:
#line 2071 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                          (yyval.ro_cols) = (yyvsp[-2].ro_cols);
                                                          (yyval.ro_cols).column_tail->next = (yyvsp[0].ro_cols).column_head;
                                                          (yyval.ro_cols).column_tail = (yyvsp[0].ro_cols).column_tail;
                                                      }
#line 5666 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 513:
#line 2079 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_single_count_expr(parse_state, & (yyval.ro_cols), (yyvsp[-1].cnt), (yyvsp[0].str))); }
#line 5672 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 514:
#line 2083 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.cnt) = mdl_new_oexpr_constant(parse_state, (yyvsp[0].dbl))); }
#line 5678 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 516:
#line 2085 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.cnt) = mdl_join_oexpr_tree(parse_state, (yyvsp[-1].cnt), NULL, '(')); }
#line 5684 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 517:
#line 2086 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.cnt) = mdl_join_oexpr_tree(parse_state, (yyvsp[-2].cnt),   (yyvsp[0].cnt), '+')); }
#line 5690 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 518:
#line 2087 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.cnt) = mdl_join_oexpr_tree(parse_state, (yyvsp[-2].cnt),   (yyvsp[0].cnt), '-')); }
#line 5696 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 519:
#line 2088 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.cnt) = mdl_join_oexpr_tree(parse_state, (yyvsp[-2].cnt),   (yyvsp[0].cnt), '*')); }
#line 5702 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 520:
#line 2089 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.cnt) = mdl_join_oexpr_tree(parse_state, (yyvsp[-2].cnt),   (yyvsp[0].cnt), '/')); }
#line 5708 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 521:
#line 2090 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.cnt) = mdl_join_oexpr_tree(parse_state, (yyvsp[0].cnt), NULL, '_')); }
#line 5714 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 522:
#line 2091 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.cnt) = mdl_sum_oexpr((yyvsp[-1].cnt))); }
#line 5720 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 523:
#line 2096 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->count_flags |= COUNT_PRESENT; }
#line 5726 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 524:
#line 2097 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.cnt) = (yyvsp[-1].cnt); }
#line 5732 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 525:
#line 2098 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.cnt) = mdl_new_oexpr_constant(parse_state, (yyvsp[-1].dbl))); }
#line 5738 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 526:
#line 2099 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->count_flags |= TRIGGER_PRESENT; }
#line 5744 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 527:
#line 2100 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.cnt) = (yyvsp[-1].cnt); }
#line 5750 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 528:
#line 2103 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = FILE_OVERWRITE; }
#line 5756 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 529:
#line 2104 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = FILE_SUBSTITUTE; }
#line 5762 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 530:
#line 2105 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = FILE_APPEND; }
#line 5768 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 531:
#line 2106 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = FILE_APPEND_HEADER; }
#line 5774 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 532:
#line 2107 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = FILE_CREATE; }
#line 5780 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 534:
#line 2113 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.sym) = mdl_existing_rxn_pathname_or_molecule(parse_state, (yyvsp[0].str))); }
#line 5786 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 535:
#line 2117 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                        (yyval.mol_type) = (yyvsp[0].mol_type);
                                                        if ((yyval.mol_type).orient > 0)
                                                          (yyval.mol_type).orient = 1;
                                                        else if ((yyval.mol_type).orient < 0)
                                                          (yyval.mol_type).orient = -1;
                                                        CHECKN((yyval.mol_type).mol_type = mdl_existing_molecule(parse_state, (yyvsp[-1].str)));
                                                      }
#line 5799 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 542:
#line 2137 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.cnt) = mdl_count_syntax_1(parse_state, (yyvsp[-3].sym), (yyvsp[-1].sym), (yyvsp[0].tok), parse_state->count_flags)); }
#line 5805 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 543:
#line 2142 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.cnt) = mdl_count_syntax_2(parse_state, (yyvsp[-3].mol_type).mol_type, (yyvsp[-3].mol_type).orient, (yyvsp[-1].sym), (yyvsp[0].tok), parse_state->count_flags)); }
#line 5811 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 544:
#line 2147 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.cnt) = mdl_count_syntax_3(parse_state, (yyvsp[-3].str), (yyvsp[-1].sym), (yyvsp[0].tok), parse_state->count_flags)); }
#line 5817 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 545:
#line 2153 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.cnt) = mdl_count_syntax_periodic_1(parse_state, (yyvsp[-5].sym), (yyvsp[-3].sym), (yyvsp[-1].vec3), (yyvsp[0].tok), parse_state->count_flags)); }
#line 5823 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 546:
#line 2157 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.cnt) = mdl_count_syntax_periodic_2(parse_state, (yyvsp[-5].mol_type).mol_type, (yyvsp[-5].mol_type).orient, (yyvsp[-3].sym), (yyvsp[-1].vec3), (yyvsp[0].tok), parse_state->count_flags)); }
#line 5829 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 547:
#line 2162 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.cnt) = mdl_count_syntax_periodic_3(parse_state, (yyvsp[-5].str), (yyvsp[-3].sym), (yyvsp[-1].vec3), (yyvsp[0].tok), parse_state->count_flags)); }
#line 5835 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 548:
#line 2165 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.sym) = NULL; }
#line 5841 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 549:
#line 2166 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.sym) = (yyvsp[0].sym); }
#line 5847 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 550:
#line 2167 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.sym) = (yyvsp[0].sym); }
#line 5853 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 551:
#line 2170 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = REPORT_NOTHING; }
#line 5859 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 552:
#line 2171 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = (yyvsp[0].tok); }
#line 5865 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 553:
#line 2174 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = REPORT_FRONT_HITS; }
#line 5871 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 554:
#line 2175 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = REPORT_BACK_HITS; }
#line 5877 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 555:
#line 2176 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = REPORT_ALL_HITS; }
#line 5883 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 556:
#line 2177 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = REPORT_FRONT_CROSSINGS; }
#line 5889 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 557:
#line 2178 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = REPORT_BACK_CROSSINGS; }
#line 5895 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 558:
#line 2179 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = REPORT_ALL_CROSSINGS; }
#line 5901 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 559:
#line 2180 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = REPORT_CONCENTRATION; }
#line 5907 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 560:
#line 2181 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = REPORT_ENCLOSED; }
#line 5913 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 561:
#line 2184 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.str) = NULL; }
#line 5919 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 562:
#line 2185 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.str) = (yyvsp[0].str); }
#line 5925 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 563:
#line 2192 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_new_viz_output_block(parse_state)); }
#line 5931 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 564:
#line 2195 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { }
#line 5937 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 567:
#line 2204 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_viz_mode(parse_state->vol->viz_blocks, CELLBLENDER_MODE)); }
#line 5943 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 568:
#line 2205 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_viz_mode(parse_state->vol->viz_blocks, (yyvsp[0].ival))); }
#line 5949 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 569:
#line 2208 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.ival) = NO_VIZ_MODE; }
#line 5955 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 570:
#line 2209 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.ival) = ASCII_MODE; }
#line 5961 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 571:
#line 2210 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.ival) = CELLBLENDER_MODE; }
#line 5967 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 573:
#line 2215 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                        if ((yyvsp[0].frame_list).frame_head)
                                                        {
                                                          (yyvsp[0].frame_list).frame_tail->next = parse_state->vol->viz_blocks->frame_data_head;
                                                          parse_state->vol->viz_blocks->frame_data_head = (yyvsp[0].frame_list).frame_head;
                                                        }
                                                      }
#line 5979 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 575:
#line 2228 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_viz_filename_prefix(parse_state, parse_state->vol->viz_blocks, (yyvsp[0].str))); }
#line 5985 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 576:
#line 2234 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.frame_list) = (yyvsp[-1].frame_list); }
#line 5991 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 578:
#line 2240 "../src/../src/mdlparse.y" /* yacc.c:1646  */
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
#line 6007 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 579:
#line 2254 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.frame_list).frame_head = (yyval.frame_list).frame_tail = NULL; }
#line 6013 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 583:
#line 2266 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_viz_state(parse_state, & (yyval.ival), (yyvsp[0].dbl))); }
#line 6019 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 584:
#line 2267 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.ival) = INCLUDE_OBJ; }
#line 6025 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 587:
#line 2277 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_viz_include_molecules(parse_state, parse_state->vol->viz_blocks, (yyvsp[-1].symlist), (yyvsp[0].ival))); }
#line 6031 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 588:
#line 2278 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_viz_include_all_molecules(parse_state->vol->viz_blocks, (yyvsp[0].ival))); }
#line 6037 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 589:
#line 2282 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.symlist) = mdl_existing_molecule_list(parse_state, (yyvsp[0].str))); }
#line 6043 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 590:
#line 2283 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.symlist) = mdl_existing_molecules_wildcard(parse_state, (yyvsp[0].str))); }
#line 6049 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 591:
#line 2287 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_new_viz_all_times(parse_state, & (yyval.nlist))); }
#line 6055 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 593:
#line 2293 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.frame_list) = (yyvsp[-1].frame_list); }
#line 6061 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 595:
#line 2299 "../src/../src/mdlparse.y" /* yacc.c:1646  */
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
#line 6079 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 596:
#line 2316 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_new_viz_mol_frames(parse_state, parse_state->vol->viz_blocks, & (yyval.frame_list), OUTPUT_BY_TIME_LIST, (yyvsp[-2].tok), & (yyvsp[0].nlist))); }
#line 6085 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 597:
#line 2320 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_new_viz_all_iterations(parse_state, & (yyval.nlist))); }
#line 6091 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 599:
#line 2327 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.frame_list) = (yyvsp[-1].frame_list); }
#line 6097 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 601:
#line 2333 "../src/../src/mdlparse.y" /* yacc.c:1646  */
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
#line 6115 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 602:
#line 2350 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_new_viz_mol_frames(parse_state, parse_state->vol->viz_blocks, & (yyval.frame_list), OUTPUT_BY_ITERATION_LIST, (yyvsp[-2].tok), & (yyvsp[0].nlist))); }
#line 6121 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 603:
#line 2353 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = ALL_MOL_DATA; }
#line 6127 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 604:
#line 2354 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = MOL_POS; }
#line 6133 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 605:
#line 2355 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = MOL_ORIENT; }
#line 6139 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 606:
#line 2369 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                          struct volume_output_item *vo;
                                                          CHECKN(vo = mdl_new_volume_output_item(parse_state, (yyvsp[-6].str), & (yyvsp[-5].species_lst), (yyvsp[-4].vec3), (yyvsp[-3].vec3), (yyvsp[-2].vec3), (yyvsp[-1].otimes)));
                                                          vo->next = parse_state->vol->volume_output_head;
                                                          parse_state->vol->volume_output_head = vo;
                                                      }
#line 6150 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 607:
#line 2378 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.str) = (yyvsp[0].str); }
#line 6156 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 609:
#line 2384 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                          (yyval.species_lst) = (yyvsp[-1].species_lst);
                                                          (yyval.species_lst).species_count += (yyvsp[0].species_lst).species_count;
                                                          (yyval.species_lst).species_tail->next = (yyvsp[0].species_lst).species_head;
                                                          (yyval.species_lst).species_tail = (yyvsp[0].species_lst).species_tail;
                                                      }
#line 6167 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 610:
#line 2393 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.species_lst) = (yyvsp[0].species_lst); }
#line 6173 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 611:
#line 2396 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                          struct sym_entry *sp;
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
#line 6193 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 612:
#line 2414 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.species_lst).species_tail = (yyval.species_lst).species_head = (yyvsp[0].species_lst_item); (yyval.species_lst).species_count = 1; }
#line 6199 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 613:
#line 2416 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                        (yyval.species_lst) = (yyvsp[-2].species_lst);
                                                        (yyval.species_lst).species_tail = (yyval.species_lst).species_tail->next = (yyvsp[0].species_lst_item);
                                                        ++ (yyval.species_lst).species_count;
                                                      }
#line 6209 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 614:
#line 2424 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.vec3) = (yyvsp[0].vec3); }
#line 6215 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 615:
#line 2428 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.vec3) = (yyvsp[0].vec3); }
#line 6221 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 616:
#line 2432 "../src/../src/mdlparse.y" /* yacc.c:1646  */
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
#line 6244 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 617:
#line 2453 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.otimes) = mdl_new_output_times_default(parse_state)); }
#line 6250 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 618:
#line 2454 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.otimes) = mdl_new_output_times_step(parse_state, (yyvsp[0].dbl))); }
#line 6256 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 619:
#line 2455 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.otimes) = mdl_new_output_times_iterations(parse_state, & (yyvsp[0].nlist))); }
#line 6262 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 620:
#line 2456 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.otimes) = mdl_new_output_times_time(parse_state, & (yyvsp[0].nlist))); }
#line 6268 "mdlparse.c" /* yacc.c:1646  */
    break;


#line 6272 "mdlparse.c" /* yacc.c:1646  */
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
#line 2459 "../src/../src/mdlparse.y" /* yacc.c:1906  */






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
  mcell_errorv_nodie(fmt, arglist);
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
int mdlparse_file(struct mdlparse_vars *parse_state, char const *name)
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
    free(err);
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
      for (struct sym_entry *symp = vol->fstream_sym_table->entries[i];
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

  if (mpv.header_comment != 0) {
    free(mpv.header_comment); 
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
