/* A Bison parser, made by GNU Bison 3.0.2.  */

/* Bison implementation for Yacc-like parsers in C

   Copyright (C) 1984, 1989-1990, 2000-2013 Free Software Foundation, Inc.

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
#define YYBISON_VERSION "3.0.2"

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

#line 138 "mdlparse.c" /* yacc.c:339  */

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
    CONCENTRATION = 295,
    CORNERS = 296,
    COS = 297,
    COUNT = 298,
    CUBIC = 299,
    CUBIC_RELEASE_SITE = 300,
    CUSTOM_SPACE_STEP = 301,
    CUSTOM_TIME_STEP = 302,
    DEFAULT = 303,
    DEFINE_MOLECULE = 304,
    DEFINE_MOLECULES = 305,
    DEFINE_REACTIONS = 306,
    DEFINE_RELEASE_PATTERN = 307,
    DEFINE_SURFACE_CLASS = 308,
    DEFINE_SURFACE_CLASSES = 309,
    DEFINE_SURFACE_REGIONS = 310,
    DEGENERATE_POLYGONS = 311,
    DELAY = 312,
    DENSITY = 313,
    DIFFUSION_CONSTANT_2D = 314,
    DIFFUSION_CONSTANT_3D = 315,
    DIFFUSION_CONSTANT_REPORT = 316,
    EFFECTOR_GRID_DENSITY = 317,
    ELEMENT_CONNECTIONS = 318,
    ELLIPTIC = 319,
    ELLIPTIC_RELEASE_SITE = 320,
    EQUAL = 321,
    ERROR = 322,
    ESTIMATE_CONCENTRATION = 323,
    EXCLUDE_ELEMENTS = 324,
    EXCLUDE_PATCH = 325,
    EXCLUDE_REGION = 326,
    EXIT = 327,
    EXP = 328,
    EXPRESSION = 329,
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
    HEADER = 346,
    HIGH_PROBABILITY_THRESHOLD = 347,
    HIGH_REACTION_PROBABILITY = 348,
    IGNORED = 349,
    INCLUDE_ELEMENTS = 350,
    INCLUDE_FILE = 351,
    INCLUDE_PATCH = 352,
    INCLUDE_REGION = 353,
    INPUT_FILE = 354,
    INSTANTIATE = 355,
    LLINTEGER = 356,
    FULLY_RANDOM = 357,
    INTERACTION_RADIUS = 358,
    ITERATION_LIST = 359,
    ITERATION_NUMBERS = 360,
    ITERATION_REPORT = 361,
    ITERATIONS = 362,
    KEEP_CHECKPOINT_FILES = 363,
    LEFT = 364,
    LIFETIME_THRESHOLD = 365,
    LIFETIME_TOO_SHORT = 366,
    LIST = 367,
    LOCATION = 368,
    LOG = 369,
    LOG10 = 370,
    MAX_TOK = 371,
    MAXIMUM_STEP_LENGTH = 372,
    MEAN_DIAMETER = 373,
    MEAN_NUMBER = 374,
    MEMORY_PARTITION_X = 375,
    MEMORY_PARTITION_Y = 376,
    MEMORY_PARTITION_Z = 377,
    MEMORY_PARTITION_POOL = 378,
    MESHES = 379,
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
    NEGATIVE_DIFFUSION_CONSTANT = 396,
    NEGATIVE_REACTION_RATE = 397,
    NO = 398,
    NOEXIT = 399,
    NONE = 400,
    NO_SPECIES = 401,
    NOT_EQUAL = 402,
    NOTIFICATIONS = 403,
    NUMBER_OF_SUBUNITS = 404,
    NUMBER_OF_TRAINS = 405,
    NUMBER_TO_RELEASE = 406,
    OBJECT = 407,
    OFF = 408,
    ON = 409,
    ORIENTATIONS = 410,
    OUTPUT_BUFFER_SIZE = 411,
    INVALID_OUTPUT_STEP_TIME = 412,
    OVERWRITTEN_OUTPUT_FILE = 413,
    PARTITION_LOCATION_REPORT = 414,
    PARTITION_X = 415,
    PARTITION_Y = 416,
    PARTITION_Z = 417,
    PI_TOK = 418,
    POLYGON_LIST = 419,
    POSITIONS = 420,
    PRINTF = 421,
    PRINT_TIME = 422,
    PROBABILITY_REPORT = 423,
    PROBABILITY_REPORT_THRESHOLD = 424,
    PROGRESS_REPORT = 425,
    RADIAL_DIRECTIONS = 426,
    RADIAL_SUBDIVISIONS = 427,
    RAND_GAUSSIAN = 428,
    RAND_UNIFORM = 429,
    RATE_RULES = 430,
    REACTION_DATA_OUTPUT = 431,
    REACTION_OUTPUT_REPORT = 432,
    REAL = 433,
    RECTANGULAR_RELEASE_SITE = 434,
    RECTANGULAR_TOKEN = 435,
    REFLECTIVE = 436,
    REGION_DATA = 437,
    RELEASE_EVENT_REPORT = 438,
    RELEASE_INTERVAL = 439,
    RELEASE_PATTERN = 440,
    RELEASE_PROBABILITY = 441,
    RELEASE_SITE = 442,
    REMOVE_ELEMENTS = 443,
    RIGHT = 444,
    ROTATE = 445,
    ROUND_OFF = 446,
    SCALE = 447,
    SEED = 448,
    SHAPE = 449,
    SHOW_EXACT_TIME = 450,
    SIN = 451,
    SITE_DIAMETER = 452,
    SITE_RADIUS = 453,
    SPACE_STEP = 454,
    SPHERICAL = 455,
    SPHERICAL_RELEASE_SITE = 456,
    SPHERICAL_SHELL = 457,
    SPHERICAL_SHELL_SITE = 458,
    SPRINTF = 459,
    SQRT = 460,
    STANDARD_DEVIATION = 461,
    STEP = 462,
    STRING_TO_NUM = 463,
    STR_VALUE = 464,
    SUBUNIT = 465,
    SUBUNIT_RELATIONSHIPS = 466,
    SUMMATION_OPERATOR = 467,
    SURFACE_CLASS = 468,
    SURFACE_ONLY = 469,
    TAN = 470,
    TARGET_ONLY = 471,
    TET_ELEMENT_CONNECTIONS = 472,
    THROUGHPUT_REPORT = 473,
    TIME_LIST = 474,
    TIME_POINTS = 475,
    TIME_STEP = 476,
    TIME_STEP_MAX = 477,
    TO = 478,
    TOP = 479,
    TRAIN_DURATION = 480,
    TRAIN_INTERVAL = 481,
    TRANSLATE = 482,
    TRANSPARENT = 483,
    TRIGGER = 484,
    TRUE = 485,
    UNLIMITED = 486,
    USELESS_VOLUME_ORIENTATION = 487,
    VACANCY_SEARCH_DISTANCE = 488,
    VAR = 489,
    VARYING_PROBABILITY_REPORT = 490,
    VERTEX_LIST = 491,
    VIZ_MESH_FORMAT = 492,
    VIZ_MOLECULE_FORMAT = 493,
    VIZ_OUTPUT = 494,
    VIZ_OUTPUT_REPORT = 495,
    VIZ_VALUE = 496,
    VOLUME_DATA_OUTPUT = 497,
    VOLUME_OUTPUT_REPORT = 498,
    VOLUME_DEPENDENT_RELEASE_NUMBER = 499,
    VOLUME_ONLY = 500,
    VOXEL_COUNT = 501,
    VOXEL_LIST = 502,
    VOXEL_SIZE = 503,
    WARNING = 504,
    WARNINGS = 505,
    WORLD = 506,
    YES = 507,
    UNARYMINUS = 508
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
#define CONCENTRATION 295
#define CORNERS 296
#define COS 297
#define COUNT 298
#define CUBIC 299
#define CUBIC_RELEASE_SITE 300
#define CUSTOM_SPACE_STEP 301
#define CUSTOM_TIME_STEP 302
#define DEFAULT 303
#define DEFINE_MOLECULE 304
#define DEFINE_MOLECULES 305
#define DEFINE_REACTIONS 306
#define DEFINE_RELEASE_PATTERN 307
#define DEFINE_SURFACE_CLASS 308
#define DEFINE_SURFACE_CLASSES 309
#define DEFINE_SURFACE_REGIONS 310
#define DEGENERATE_POLYGONS 311
#define DELAY 312
#define DENSITY 313
#define DIFFUSION_CONSTANT_2D 314
#define DIFFUSION_CONSTANT_3D 315
#define DIFFUSION_CONSTANT_REPORT 316
#define EFFECTOR_GRID_DENSITY 317
#define ELEMENT_CONNECTIONS 318
#define ELLIPTIC 319
#define ELLIPTIC_RELEASE_SITE 320
#define EQUAL 321
#define ERROR 322
#define ESTIMATE_CONCENTRATION 323
#define EXCLUDE_ELEMENTS 324
#define EXCLUDE_PATCH 325
#define EXCLUDE_REGION 326
#define EXIT 327
#define EXP 328
#define EXPRESSION 329
#define FALSE 330
#define FCLOSE 331
#define FILENAME 332
#define FILENAME_PREFIX 333
#define FILE_OUTPUT_REPORT 334
#define FINAL_SUMMARY 335
#define FLOOR 336
#define FOPEN 337
#define FORMAT 338
#define FPRINTF 339
#define FPRINT_TIME 340
#define FRONT 341
#define FRONT_CROSSINGS 342
#define FRONT_HITS 343
#define GAUSSIAN_RELEASE_NUMBER 344
#define GEOMETRY 345
#define HEADER 346
#define HIGH_PROBABILITY_THRESHOLD 347
#define HIGH_REACTION_PROBABILITY 348
#define IGNORED 349
#define INCLUDE_ELEMENTS 350
#define INCLUDE_FILE 351
#define INCLUDE_PATCH 352
#define INCLUDE_REGION 353
#define INPUT_FILE 354
#define INSTANTIATE 355
#define LLINTEGER 356
#define FULLY_RANDOM 357
#define INTERACTION_RADIUS 358
#define ITERATION_LIST 359
#define ITERATION_NUMBERS 360
#define ITERATION_REPORT 361
#define ITERATIONS 362
#define KEEP_CHECKPOINT_FILES 363
#define LEFT 364
#define LIFETIME_THRESHOLD 365
#define LIFETIME_TOO_SHORT 366
#define LIST 367
#define LOCATION 368
#define LOG 369
#define LOG10 370
#define MAX_TOK 371
#define MAXIMUM_STEP_LENGTH 372
#define MEAN_DIAMETER 373
#define MEAN_NUMBER 374
#define MEMORY_PARTITION_X 375
#define MEMORY_PARTITION_Y 376
#define MEMORY_PARTITION_Z 377
#define MEMORY_PARTITION_POOL 378
#define MESHES 379
#define MICROSCOPIC_REVERSIBILITY 380
#define MIN_TOK 381
#define MISSED_REACTIONS 382
#define MISSED_REACTION_THRESHOLD 383
#define MISSING_SURFACE_ORIENTATION 384
#define MOD 385
#define MODE 386
#define MODIFY_SURFACE_REGIONS 387
#define MOLECULE 388
#define MOLECULE_COLLISION_REPORT 389
#define MOLECULE_DENSITY 390
#define MOLECULE_NUMBER 391
#define MOLECULE_POSITIONS 392
#define MOLECULES 393
#define MOLECULE_PLACEMENT_FAILURE 394
#define NAME_LIST 395
#define NEGATIVE_DIFFUSION_CONSTANT 396
#define NEGATIVE_REACTION_RATE 397
#define NO 398
#define NOEXIT 399
#define NONE 400
#define NO_SPECIES 401
#define NOT_EQUAL 402
#define NOTIFICATIONS 403
#define NUMBER_OF_SUBUNITS 404
#define NUMBER_OF_TRAINS 405
#define NUMBER_TO_RELEASE 406
#define OBJECT 407
#define OFF 408
#define ON 409
#define ORIENTATIONS 410
#define OUTPUT_BUFFER_SIZE 411
#define INVALID_OUTPUT_STEP_TIME 412
#define OVERWRITTEN_OUTPUT_FILE 413
#define PARTITION_LOCATION_REPORT 414
#define PARTITION_X 415
#define PARTITION_Y 416
#define PARTITION_Z 417
#define PI_TOK 418
#define POLYGON_LIST 419
#define POSITIONS 420
#define PRINTF 421
#define PRINT_TIME 422
#define PROBABILITY_REPORT 423
#define PROBABILITY_REPORT_THRESHOLD 424
#define PROGRESS_REPORT 425
#define RADIAL_DIRECTIONS 426
#define RADIAL_SUBDIVISIONS 427
#define RAND_GAUSSIAN 428
#define RAND_UNIFORM 429
#define RATE_RULES 430
#define REACTION_DATA_OUTPUT 431
#define REACTION_OUTPUT_REPORT 432
#define REAL 433
#define RECTANGULAR_RELEASE_SITE 434
#define RECTANGULAR_TOKEN 435
#define REFLECTIVE 436
#define REGION_DATA 437
#define RELEASE_EVENT_REPORT 438
#define RELEASE_INTERVAL 439
#define RELEASE_PATTERN 440
#define RELEASE_PROBABILITY 441
#define RELEASE_SITE 442
#define REMOVE_ELEMENTS 443
#define RIGHT 444
#define ROTATE 445
#define ROUND_OFF 446
#define SCALE 447
#define SEED 448
#define SHAPE 449
#define SHOW_EXACT_TIME 450
#define SIN 451
#define SITE_DIAMETER 452
#define SITE_RADIUS 453
#define SPACE_STEP 454
#define SPHERICAL 455
#define SPHERICAL_RELEASE_SITE 456
#define SPHERICAL_SHELL 457
#define SPHERICAL_SHELL_SITE 458
#define SPRINTF 459
#define SQRT 460
#define STANDARD_DEVIATION 461
#define STEP 462
#define STRING_TO_NUM 463
#define STR_VALUE 464
#define SUBUNIT 465
#define SUBUNIT_RELATIONSHIPS 466
#define SUMMATION_OPERATOR 467
#define SURFACE_CLASS 468
#define SURFACE_ONLY 469
#define TAN 470
#define TARGET_ONLY 471
#define TET_ELEMENT_CONNECTIONS 472
#define THROUGHPUT_REPORT 473
#define TIME_LIST 474
#define TIME_POINTS 475
#define TIME_STEP 476
#define TIME_STEP_MAX 477
#define TO 478
#define TOP 479
#define TRAIN_DURATION 480
#define TRAIN_INTERVAL 481
#define TRANSLATE 482
#define TRANSPARENT 483
#define TRIGGER 484
#define TRUE 485
#define UNLIMITED 486
#define USELESS_VOLUME_ORIENTATION 487
#define VACANCY_SEARCH_DISTANCE 488
#define VAR 489
#define VARYING_PROBABILITY_REPORT 490
#define VERTEX_LIST 491
#define VIZ_MESH_FORMAT 492
#define VIZ_MOLECULE_FORMAT 493
#define VIZ_OUTPUT 494
#define VIZ_OUTPUT_REPORT 495
#define VIZ_VALUE 496
#define VOLUME_DATA_OUTPUT 497
#define VOLUME_OUTPUT_REPORT 498
#define VOLUME_DEPENDENT_RELEASE_NUMBER 499
#define VOLUME_ONLY 500
#define VOXEL_COUNT 501
#define VOXEL_LIST 502
#define VOXEL_SIZE 503
#define WARNING 504
#define WARNINGS 505
#define WORLD 506
#define YES 507
#define UNARYMINUS 508

/* Value type.  */
#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED
typedef union YYSTYPE YYSTYPE;
union YYSTYPE
{
#line 68 "../src/../src/mdlparse.y" /* yacc.c:355  */

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


#line 750 "mdlparse.c" /* yacc.c:355  */
};
# define YYSTYPE_IS_TRIVIAL 1
# define YYSTYPE_IS_DECLARED 1
#endif



int mdlparse (struct mdlparse_vars *parse_state, yyscan_t scanner);

#endif /* !YY_MDL_MDLPARSE_H_INCLUDED  */

/* Copy the second part of user declarations.  */

#line 764 "mdlparse.c" /* yacc.c:358  */

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
#define YYFINAL  142
/* YYLAST -- Last index in YYTABLE.  */
#define YYLAST   2588

/* YYNTOKENS -- Number of terminals.  */
#define YYNTOKENS  274
/* YYNNTS -- Number of nonterminals.  */
#define YYNNTS  279
/* YYNRULES -- Number of rules.  */
#define YYNRULES  596
/* YYNSTATES -- Number of states.  */
#define YYNSTATES  1172

/* YYTRANSLATE[YYX] -- Symbol number corresponding to YYX as returned
   by yylex, with out-of-bounds checking.  */
#define YYUNDEFTOK  2
#define YYMAXUTOK   508

#define YYTRANSLATE(YYX)                                                \
  ((unsigned int) (YYX) <= YYMAXUTOK ? yytranslate[YYX] : YYUNDEFTOK)

/* YYTRANSLATE[TOKEN-NUM] -- Symbol number corresponding to TOKEN-NUM
   as returned by yylex, without out-of-bounds checking.  */
static const yytype_uint16 yytranslate[] =
{
       0,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,   254,   265,
     269,   270,   258,   256,   266,   257,     2,   259,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,   255,   264,
     272,   253,   271,     2,   273,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,   262,     2,   263,   260,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,   267,     2,   268,     2,     2,     2,     2,
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
     245,   246,   247,   248,   249,   250,   251,   252,   261
};

#if YYDEBUG
  /* YYRLINE[YYN] -- Source line where rule number YYN was defined.  */
static const yytype_uint16 yyrline[] =
{
       0,   588,   588,   592,   593,   598,   599,   600,   601,   602,
     603,   604,   605,   606,   607,   608,   609,   610,   611,   612,
     613,   614,   615,   616,   617,   622,   625,   628,   631,   634,
     637,   640,   641,   644,   645,   646,   647,   648,   649,   652,
     653,   654,   655,   659,   660,   661,   671,   683,   686,   689,
     701,   702,   715,   716,   722,   745,   746,   747,   748,   751,
     754,   757,   758,   769,   772,   775,   776,   779,   780,   783,
     784,   787,   788,   791,   795,   796,   797,   798,   799,   800,
     801,   802,   803,   804,   805,   806,   807,   808,   809,   810,
     811,   812,   813,   814,   815,   816,   817,   818,   819,   820,
     821,   822,   823,   824,   828,   829,   833,   834,   835,   836,
     839,   845,   846,   847,   848,   849,   850,   851,   854,   858,
     861,   864,   867,   870,   873,   874,   884,   885,   886,   898,
     902,   908,   913,   917,   926,   930,   931,   935,   936,   937,
     938,   939,   940,   941,   942,   943,   944,   945,   946,   947,
     948,   949,   950,   951,   955,   956,   960,   964,   965,   972,
     976,   977,   981,   982,   983,   984,   985,   986,   987,   988,
     989,   990,   991,   992,   993,   994,   995,   996,  1000,  1001,
    1002,  1008,  1009,  1010,  1011,  1012,  1016,  1017,  1018,  1022,
    1023,  1024,  1025,  1033,  1034,  1035,  1036,  1037,  1038,  1039,
    1040,  1041,  1042,  1043,  1044,  1045,  1046,  1047,  1054,  1055,
    1056,  1057,  1061,  1065,  1066,  1067,  1074,  1075,  1078,  1082,
    1086,  1087,  1091,  1099,  1102,  1106,  1107,  1111,  1112,  1121,
    1132,  1133,  1137,  1138,  1148,  1152,  1156,  1169,  1170,  1175,
    1179,  1185,  1186,  1191,  1191,  1196,  1199,  1201,  1206,  1207,
    1211,  1214,  1220,  1225,  1226,  1227,  1230,  1231,  1234,  1238,
    1242,  1249,  1253,  1262,  1266,  1275,  1282,  1287,  1288,  1291,
    1294,  1295,  1298,  1299,  1300,  1303,  1308,  1314,  1315,  1316,
    1317,  1320,  1321,  1325,  1330,  1331,  1334,  1338,  1339,  1343,
    1346,  1347,  1350,  1351,  1355,  1356,  1359,  1370,  1387,  1388,
    1389,  1393,  1394,  1395,  1402,  1409,  1412,  1416,  1423,  1425,
    1427,  1429,  1431,  1435,  1436,  1443,  1443,  1454,  1457,  1458,
    1459,  1460,  1461,  1470,  1473,  1476,  1479,  1481,  1485,  1489,
    1490,  1491,  1496,  1509,  1510,  1513,  1514,  1519,  1518,  1525,
    1526,  1531,  1530,  1538,  1539,  1540,  1541,  1542,  1543,  1544,
    1545,  1552,  1553,  1554,  1555,  1556,  1561,  1560,  1567,  1568,
    1569,  1570,  1571,  1575,  1576,  1579,  1583,  1584,  1585,  1592,
    1593,  1594,  1595,  1596,  1598,  1603,  1604,  1608,  1609,  1610,
    1611,  1616,  1617,  1623,  1630,  1638,  1639,  1643,  1644,  1649,
    1652,  1660,  1657,  1675,  1678,  1681,  1682,  1686,  1691,  1692,
    1696,  1699,  1701,  1707,  1708,  1712,  1712,  1724,  1725,  1728,
    1729,  1730,  1731,  1732,  1733,  1734,  1738,  1739,  1744,  1745,
    1746,  1747,  1751,  1756,  1760,  1764,  1765,  1768,  1769,  1770,
    1773,  1776,  1777,  1780,  1783,  1784,  1788,  1794,  1795,  1800,
    1801,  1800,  1811,  1808,  1821,  1825,  1829,  1833,  1844,  1845,
    1841,  1853,  1854,  1868,  1874,  1875,  1880,  1881,  1883,  1880,
    1892,  1895,  1897,  1902,  1903,  1907,  1914,  1920,  1921,  1926,
    1926,  1936,  1935,  1946,  1947,  1958,  1959,  1960,  1963,  1967,
    1975,  1982,  1983,  1997,  1998,  1999,  2003,  2003,  2009,  2010,
    2011,  2015,  2019,  2023,  2024,  2033,  2037,  2038,  2039,  2040,
    2041,  2042,  2043,  2044,  2045,  2050,  2050,  2052,  2053,  2053,
    2057,  2058,  2059,  2060,  2061,  2064,  2067,  2071,  2081,  2082,
    2083,  2087,  2092,  2097,  2101,  2102,  2103,  2106,  2107,  2110,
    2111,  2112,  2113,  2114,  2115,  2116,  2117,  2120,  2121,  2128,
    2128,  2135,  2136,  2140,  2141,  2144,  2145,  2146,  2150,  2151,
    2161,  2164,  2168,  2174,  2175,  2190,  2191,  2192,  2196,  2202,
    2203,  2207,  2208,  2212,  2214,  2218,  2219,  2223,  2224,  2227,
    2233,  2234,  2251,  2256,  2257,  2261,  2267,  2268,  2285,  2289,
    2290,  2291,  2298,  2314,  2318,  2319,  2329,  2332,  2350,  2351,
    2360,  2364,  2368,  2389,  2390,  2391,  2392
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
  "CLAMP_CONCENTRATION", "CLOSE_PARTITION_SPACING", "CONCENTRATION",
  "CORNERS", "COS", "COUNT", "CUBIC", "CUBIC_RELEASE_SITE",
  "CUSTOM_SPACE_STEP", "CUSTOM_TIME_STEP", "DEFAULT", "DEFINE_MOLECULE",
  "DEFINE_MOLECULES", "DEFINE_REACTIONS", "DEFINE_RELEASE_PATTERN",
  "DEFINE_SURFACE_CLASS", "DEFINE_SURFACE_CLASSES",
  "DEFINE_SURFACE_REGIONS", "DEGENERATE_POLYGONS", "DELAY", "DENSITY",
  "DIFFUSION_CONSTANT_2D", "DIFFUSION_CONSTANT_3D",
  "DIFFUSION_CONSTANT_REPORT", "EFFECTOR_GRID_DENSITY",
  "ELEMENT_CONNECTIONS", "ELLIPTIC", "ELLIPTIC_RELEASE_SITE", "EQUAL",
  "ERROR", "ESTIMATE_CONCENTRATION", "EXCLUDE_ELEMENTS", "EXCLUDE_PATCH",
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
  "REACTION_DATA_OUTPUT", "REACTION_OUTPUT_REPORT", "REAL",
  "RECTANGULAR_RELEASE_SITE", "RECTANGULAR_TOKEN", "REFLECTIVE",
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
  "existing_region", "point", "point_or_num", "boolean",
  "orientation_class", "list_orient_marks", "head_mark", "tail_mark",
  "orient_class_number", "list_range_specs", "range_spec", "include_stmt",
  "assignment_stmt", "assign_var", "existing_var_only", "array_value",
  "array_expr_only", "existing_array", "num_expr", "num_value",
  "intOrReal", "num_expr_only", "existing_num_var", "arith_expr",
  "str_expr", "str_expr_only", "existing_str_var", "io_stmt", "fopen_stmt",
  "new_file_stream", "file_mode", "fclose_stmt", "existing_file_stream",
  "format_string", "list_args", "list_arg", "printf_stmt", "fprintf_stmt",
  "sprintf_stmt", "print_time_stmt", "fprint_time_stmt",
  "notification_def", "notification_list", "notification_item_def",
  "notify_bilevel", "notify_level", "warnings_def", "warning_list",
  "warning_item_def", "warning_level", "chkpt_stmt", "exit_or_no",
  "time_expr", "parameter_def", "memory_partition_def", "partition_def",
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
  "list_tet_arrays", "box_def", "$@11", "$@12", "opt_aspect_ratio_def",
  "existing_obj_define_surface_regions",
  "list_existing_obj_surface_region_defs",
  "existing_obj_surface_region_def", "$@13", "$@14", "$@15", "new_region",
  "list_opt_surface_region_stmts", "opt_surface_region_stmt",
  "set_surface_class_stmt", "mod_surface_regions",
  "list_existing_surface_region_refs", "existing_surface_region_ref",
  "$@16", "output_def", "$@17", "output_buffer_size_def",
  "output_timer_def", "step_time_def", "iteration_time_def",
  "real_time_def", "list_count_cmds", "count_cmd", "count_stmt", "$@18",
  "custom_header_value", "custom_header", "exact_time_toggle",
  "list_count_exprs", "single_count_expr", "count_expr", "count_value",
  "$@19", "$@20", "file_arrow", "outfile_syntax",
  "existing_rxpn_or_molecule", "existing_molecule_required_orient_braces",
  "count_syntax", "count_syntax_1", "count_syntax_2", "count_syntax_3",
  "count_location_specifier", "opt_hit_spec", "hit_spec",
  "opt_custom_header", "viz_output_def", "$@21", "list_viz_output_cmds",
  "viz_output_maybe_mode_cmd", "viz_mode_def", "viz_output_cmd",
  "viz_frames_def", "viz_filename_prefix_def", "viz_molecules_block_def",
  "list_viz_molecules_block_cmds", "viz_molecules_block_cmd",
  "viz_molecules_name_list_cmd", "optional_state",
  "viz_include_mols_cmd_list", "viz_include_mols_cmd",
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
     505,   506,   507,    61,    38,    58,    43,    45,    42,    47,
      94,   508,    91,    93,    59,    39,    44,   123,   125,    40,
      41,    62,    60,    64
};
# endif

#define YYPACT_NINF -1009

#define yypact_value_is_default(Yystate) \
  (!!((Yystate) == (-1009)))

#define YYTABLE_NINF -391

#define yytable_value_is_error(Yytable_value) \
  0

  /* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
     STATE-NUM.  */
static const yytype_int16 yypact[] =
{
    2338,  -108,   -47,    36,    92,   109,   130,   -20,    29,    41,
     -20,   -20,   132,   148,   181,   138,   158,   179,   210, -1009,
     232,   239,   272,   291,   302,   305,   307,   310,   191,   300,
   -1009, -1009, -1009,   306,   316,   330,   334,   337,   359,   340,
     372,   374,   378, -1009,   396,   399,   414,   499,  2338, -1009,
     -42, -1009, -1009,   391, -1009, -1009,   564, -1009, -1009, -1009,
   -1009, -1009, -1009, -1009, -1009, -1009, -1009, -1009, -1009,   395,
   -1009, -1009, -1009, -1009, -1009, -1009, -1009, -1009, -1009, -1009,
   -1009,   174, -1009, -1009, -1009, -1009,   495, -1009, -1009, -1009,
   -1009, -1009, -1009, -1009, -1009,    67,    67,   -40,  2171,   -40,
    2171, -1009, -1009,   422,   -20,   -20, -1009,   432, -1009,   433,
   -1009,   -20,   -20,  2171,   -20,   -20,   -20,   -40,   -20,  2171,
    2171,    67,  2171,  2171,  2171,  2171,   381,   -20,   339,   -40,
     -40,  1810,  2171,   549,  2171,   -20,  2171,  2171,  2171, -1009,
     624,   420, -1009, -1009,  1585,   437,  -162,   440, -1009, -1009,
     440, -1009,   440, -1009, -1009,   440,   440,   440, -1009, -1009,
   -1009, -1009, -1009, -1009, -1009, -1009,   442, -1009, -1009, -1009,
   -1009, -1009,   455, -1009, -1009,   444,   445,   446,   448,   449,
     452,   453,   454, -1009,   460,   464,   465,   473,   474, -1009,
   -1009, -1009, -1009,   475, -1009,   476,   479,   480,   481,  2171,
    2171,  2171, -1009,   128, -1009, -1009, -1009, -1009, -1009,   952,
      -7,   530,  -145, -1009, -1009,   356, -1009,    71, -1009, -1009,
      54, -1009, -1009, -1009,    75, -1009, -1009, -1009,   129, -1009,
     781, -1009,   470,   458,   486,   455, -1009,   606, -1009,   781,
     781, -1009,   781,   781,   781,   781, -1009, -1009, -1009,   502,
     501,   140, -1009,   518,   520,   524,   525,   527,   538,   544,
     551,   556,   557,   566,   567,   577,   583,   584,   593,   594,
     596,   234, -1009,   455, -1009,   519, -1009,   781,   781,   597,
   -1009,   781, -1009,   588,   781,   781,   781,   686,   598,   718,
     604,   605,   607,   608,   609,   613,   614,   616,   617,   620,
     621,   622,   626,   627,   633,   655,    47, -1009,  1860,   631,
   -1009, -1009,   781,   905, -1009,  1027,   455,   623,   -40, -1009,
   -1009, -1009, -1009, -1009,   870,   -20, -1009,   676, -1009,   676,
     -40,   -40,  2171,  2171,  2171,  2171,  2171,  2171,  2171,  2171,
    2171,  2171,  2171,  2171,  2171,  2171,  2171,  2171,   -40,  2171,
   -1009, -1009,   343, -1009, -1009,  2171,  2171,  2171,  2171,  2171,
   -1009,  2171, -1009,   660,   661,   550, -1009, -1009, -1009, -1009,
   -1009,  2171, -1009,   377, -1009, -1009, -1009, -1009, -1009,   -20,
     -20,  -177,     1, -1009, -1009, -1009,   658, -1009, -1009, -1009,
     -40,   -40,   -20, -1009, -1009, -1009,    67,    67,    67,    18,
      67,    67,  1498,    67,    67,    67,  2171,    67,    18,    67,
      67,    67,    18,    18, -1009, -1009,    85, -1009,  2171,    28,
     -40,   677,   164, -1009,   -40,   682,   -30, -1009,   -34,   -34,
     -34,  2171,   -34,  2171,   -34,   -34,  2171,   -34,   -34,   -34,
     -34,   -34,   -34,   -34, -1009, -1009,  2171,  -181, -1009,   781,
     670,   684,   774, -1009,   246,   -20, -1009, -1009,   750,   693,
     746,   443,   898, -1009, -1009,   379,   497,   528,   543,   585,
     669,   699,   715,   744,   766,   872,   883,   888,   900,   772,
     789,   112,   811, -1009,   261,   261,   704,   704, -1009,  1022,
    2171,  2171,   712,   717,   752,   397, -1009, -1009, -1009, -1009,
     356, -1009, -1009,   719,  -161, -1009,  -147, -1009, -1009, -1009,
     -12,   724,   727,   728,   729,   730, -1009,    32,   -20, -1009,
     720,   725, -1009, -1009, -1009, -1009, -1009, -1009, -1009, -1009,
   -1009, -1009,   781, -1009, -1009, -1009, -1009,   781, -1009, -1009,
   -1009, -1009, -1009, -1009, -1009,  1723, -1009,   781,   734,   736,
     745,   -29, -1009, -1009, -1009, -1009,    50,   754,   738,   -31,
   -1009, -1009, -1009, -1009,   455,   -20,   755, -1009,   749, -1009,
   -1009, -1009, -1009, -1009, -1009,   781, -1009,   781, -1009, -1009,
     781, -1009, -1009, -1009, -1009, -1009, -1009, -1009,   230, -1009,
    1860,   -40,  -162,  -103,  -121, -1009,   757,   443,  -162,   748,
   -1009,   760,   763,   753,   764,   780,   767,   790,   799,   800,
   -1009, -1009,   787,   443, -1009,   802, -1009, -1009, -1009, -1009,
   -1009,   791, -1009,   139, -1009, -1009, -1009, -1009, -1009, -1009,
   -1009, -1009, -1009, -1009,  2171,  2171,  2171,  2171, -1009, -1009,
   -1009, -1009,  2171,   781,   781,  2171,  2171, -1009,   902, -1009,
   -1009,   804, -1009, -1009,   719, -1009,   719, -1009, -1009,  -172,
   -1009,  2171,  1947,  2171,  2171,  2171, -1009,   -20,   795,   796,
   -1009, -1009, -1009, -1009, -1009,    83, -1009, -1009, -1009,   801,
     205, -1009, -1009,   -59, -1009, -1009,   623, -1009,  -162,  2171,
    -162,   819,   821, -1009,     7, -1009, -1009, -1009, -1009,   225,
   -1009, -1009, -1009,   -40,   -39, -1009, -1009, -1009, -1009,   826,
    -162,   830,   852,  2171, -1009,   455,   827,   840, -1009,   440,
     859,   860,   862, -1009, -1009, -1009, -1009,    -3,   443, -1009,
   -1009,   -16,  -162, -1009,  2171,  2171,   997,  -162,   -20,   -20,
    2171,   -20,  2171,   999,  -121, -1009,  2084,  -162, -1009, -1009,
     820,   829,   835,   844,  1112,   781,   781,   880,   884,    81,
   -1009, -1009,   -12,  1054,   896, -1009, -1009,   781, -1009,   781,
   -1009,   781,   781,   781,   914,   -20,   -20, -1009, -1009,    24,
   -1009, -1009,   923, -1009, -1009, -1009, -1009, -1009,   781, -1009,
     629,    67,    78, -1009, -1009, -1009,   455,   910,   912,   916,
     -38, -1009, -1009, -1009, -1009,   -20, -1009,  2084,   928,    58,
     222, -1009,  -162, -1009,  -162,  2084,  -162, -1009, -1009, -1009,
   -1009, -1009, -1009,    52,   502, -1009,   358,  -121, -1009, -1009,
   -1009, -1009,   176,  -121,   781,   781,   937, -1009, -1009,  -162,
     144, -1009,   781, -1009, -1009,   781,   938, -1009,  1033, -1009,
   -1009, -1009, -1009,   192, -1009,    -1, -1009, -1009, -1009, -1009,
    2171,  2171, -1009, -1009,  1723,  1723, -1009, -1009,   623,   169,
   -1009,   -20, -1009,  2171,   356,   940,   160, -1009,   162, -1009,
     356, -1009,   931,   -20, -1009, -1009,   455, -1009, -1009, -1009,
     920,   930, -1009,    78,    78, -1009,   204, -1009,  1118, -1009,
      43, -1009,    43, -1009, -1009, -1009,  1033, -1009, -1009, -1009,
    2084,   947,   948,   950,   936,  2171,  1193, -1009,   949, -1009,
   -1009,   349,    52,    52,    52, -1009, -1009, -1009, -1009,  2171,
   -1009, -1009, -1009,  2171, -1009, -1009,   951,   953,   154, -1009,
   -1009, -1009,   781,   781, -1009, -1009, -1009,  1054, -1009,   781,
   -1009,  2171, -1009, -1009, -1009, -1009, -1009,   471, -1009,   954,
    2171,    78,   957, -1009,   314,    78,   244,   -40,    78,    78,
      78,    78, -1009, -1009, -1009, -1009,     8, -1009,   941,    20,
      15, -1009,   956, -1009,  -162,  2171,  -162, -1009,   915,   973,
   -1009,  -121,  2171, -1009,   972,   972, -1009,   183,   208,   -20,
   -1009, -1009,   968,   781,   980, -1009, -1009,   981, -1009, -1009,
     471, -1009, -1009, -1009, -1009,   982, -1009,   983,    90,   929,
     536,    90, -1009, -1009,   967,   970,   971,   -40,   455,   392,
     392, -1009, -1009, -1009, -1009,    16,   986, -1009, -1009, -1009,
   -1009,   986, -1009, -1009,    23, -1009,   781, -1009, -1009,  2171,
   -1009, -1009,   781,   991, -1009,   993,   189, -1009,   984,  1234,
   -1009,   990,   992, -1009, -1009,   -20,  -162,   995,   998,  1000,
    1001,   985, -1009, -1009, -1009, -1009, -1009,   994, -1009, -1009,
     987, -1009, -1009, -1009, -1009, -1009,  2171, -1009, -1009, -1009,
   -1009, -1009,   781,    -1,  2171,  2171, -1009, -1009, -1009, -1009,
   -1009, -1009, -1009, -1009, -1009, -1009,   294,  1002, -1009,   471,
   -1009,  1005, -1009,  1336,  1336,    48, -1009,  1006,    69, -1009,
      69,    69, -1009, -1009, -1009,   781, -1009,   851,   113,   471,
    2171, -1009,  1336,   248,   280, -1009,  -162, -1009,   502, -1009,
    1008,  1008,  1008,  -121, -1009,  1015,   471,   781, -1009, -1009,
   -1009, -1009,    68, -1009, -1009, -1009, -1009,  2171, -1009, -1009,
   -1009, -1009, -1009, -1009, -1009, -1009, -1009, -1009,   866,    55,
   -1009, -1009
};

  /* YYDEFACT[STATE-NUM] -- Default reduction number in state STATE-NUM.
     Performed when YYTABLE does not specify something else to do.  Zero
     means the default is an error.  */
static const yytype_uint16 yydefact[] =
{
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,   315,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     213,   214,   215,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,    26,     0,     0,     0,     0,     2,     3,
     323,     5,     6,     0,     7,   111,     0,   112,   113,   114,
     115,   116,   117,     8,     9,    10,    11,    13,    12,     0,
      14,   216,   217,    15,   237,   238,    16,    17,    19,    18,
     317,     0,   318,   319,   340,   339,     0,   321,   322,   320,
      20,    21,    22,    23,    24,     0,     0,     0,     0,     0,
       0,   223,   218,     0,     0,     0,   305,     0,   224,     0,
     239,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,   473,     0,     0,     0,     0,     0,   539,
       0,     0,     1,     4,     0,     0,     0,     0,   359,   360,
       0,   361,     0,   358,   362,     0,     0,     0,    34,    36,
      38,    37,    33,    35,   198,   197,     0,   107,    25,   106,
     110,   181,    27,   104,   105,     0,     0,     0,     0,     0,
       0,     0,     0,    69,     0,     0,     0,     0,     0,    92,
      94,    93,    70,     0,    95,     0,     0,     0,     0,     0,
       0,     0,    73,   186,    65,    67,    68,    66,   182,   189,
     186,     0,     0,   220,   234,    39,   286,     0,   267,   269,
     287,   284,   307,   243,     0,   241,    28,   456,     0,   454,
     203,   122,     0,     0,     0,    54,   323,     0,   316,   204,
     196,   184,   208,   209,   210,   211,   206,   207,   205,     0,
       0,     0,   467,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,   135,   123,   124,     0,   201,   200,   202,     0,
     471,   194,    59,     0,   193,   195,   199,   543,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,   160,     0,    60,
      57,    58,     0,    71,    55,    72,     0,    56,     0,    64,
     212,    61,    62,   324,     0,     0,   341,     0,   356,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     103,   102,     0,   188,   187,     0,     0,     0,     0,     0,
     183,     0,   185,     0,     0,   227,   219,   221,    42,    47,
      48,     0,   236,    40,    43,    44,    41,   266,   268,     0,
       0,     0,     0,   246,   240,   242,     0,   453,   455,   121,
       0,     0,     0,   469,   466,   468,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,   134,   136,     0,   132,     0,     0,
       0,     0,     0,   544,     0,     0,     0,   584,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,   159,   161,     0,     0,    50,    52,
       0,     0,   323,   336,     0,   326,   333,   335,     0,     0,
       0,     0,     0,   124,   108,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,    74,    97,    98,    99,   100,   101,   190,
       0,     0,     0,     0,   230,     0,    45,    46,   285,   245,
      39,   288,   270,     0,     0,   277,     0,   279,   278,   280,
       0,     0,     0,     0,     0,     0,   304,     0,     0,   124,
       0,     0,   461,   156,   137,   144,   152,   158,   157,   139,
     146,   147,   154,   153,   155,   143,   140,   142,   138,   149,
     145,   148,   141,   151,   150,     0,   129,   474,     0,     0,
       0,     0,   475,   476,   477,   124,     0,     0,     0,     0,
     541,   549,   548,   550,   583,     0,     0,   585,     0,   180,
     178,   179,   162,   167,   168,   166,   165,   171,   170,   172,
     173,   174,   176,   163,   164,   177,   169,   175,     0,    63,
       0,     0,     0,     0,     0,   334,     0,     0,     0,     0,
     442,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     375,   376,     0,   326,   363,     0,   368,   377,   378,   379,
     380,     0,   391,     0,    90,    87,    86,    88,    82,    84,
      75,    81,    76,    77,     0,     0,     0,     0,    83,    89,
      96,    85,     0,   226,   225,     0,     0,   231,   232,    49,
     289,   273,   271,   272,     0,   274,     0,   292,   293,     0,
     290,     0,     0,     0,     0,     0,   255,     0,     0,     0,
     253,   254,   244,   247,   248,     0,   249,   258,   460,     0,
       0,   133,    29,     0,   128,   126,   127,   125,     0,     0,
       0,     0,     0,   486,     0,   481,   483,   484,   485,     0,
     546,   547,   545,     0,     0,   540,   542,   587,   588,   586,
       0,     0,     0,     0,    51,   120,     0,     0,    30,     0,
       0,     0,     0,   325,   332,   327,   328,     0,   326,   394,
     395,     0,     0,   326,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,   364,     0,     0,   401,   109,
       0,     0,     0,     0,   191,   229,   228,     0,     0,     0,
     275,   276,     0,     0,   281,   294,   295,   308,   314,   313,
     312,   309,   311,   310,     0,     0,     0,   257,   256,     0,
     457,   130,     0,   470,   464,   462,   463,   479,   478,   480,
       0,     0,     0,   472,   482,   131,   551,     0,     0,     0,
       0,   553,   555,   556,   557,     0,   590,     0,     0,   593,
       0,   118,     0,   337,     0,     0,     0,   346,   347,   350,
     348,   345,   349,     0,   344,   351,   343,     0,   393,   396,
     445,   446,     0,     0,   385,   386,     0,   366,   367,     0,
       0,   387,   381,   306,   373,   372,     0,   357,   365,   370,
     369,   371,   400,     0,   398,   326,    78,    79,    91,    80,
       0,     0,   222,   291,     0,     0,   303,   301,   302,     0,
     298,     0,   283,     0,    39,     0,     0,   261,     0,   263,
      39,   250,     0,     0,   488,   489,   490,   491,   492,   505,
       0,     0,   508,     0,     0,   496,     0,   493,   537,   497,
       0,   562,     0,   552,   554,   589,    64,    31,   591,    32,
       0,     0,     0,     0,     0,     0,   451,   326,     0,   330,
     329,     0,     0,     0,     0,   342,   444,   447,   443,     0,
     389,   374,   388,     0,   397,   399,     0,     0,     0,   402,
     403,   404,   192,   233,   299,   300,   296,     0,   282,   252,
     235,     0,   259,   262,   260,   264,   251,     0,   465,     0,
       0,     0,     0,   503,     0,     0,     0,     0,     0,     0,
       0,     0,   495,   579,   581,   580,     0,   576,     0,     0,
       0,   570,     0,   592,     0,     0,     0,   582,     0,     0,
     448,     0,     0,   352,   353,   354,   355,     0,     0,     0,
     405,   392,     0,   265,     0,   435,   432,     0,   434,   431,
     458,   416,   418,   419,   420,     0,   421,     0,     0,     0,
       0,     0,   498,   494,     0,     0,   510,     0,   538,   499,
     500,   501,   502,   575,   577,     0,   560,   558,   566,   565,
     561,   560,   569,   571,     0,   595,   594,   596,    53,     0,
     401,   338,   331,     0,   382,     0,     0,   437,     0,     0,
     297,     0,     0,   417,   461,     0,     0,     0,   516,     0,
       0,     0,   518,   519,   520,   507,   504,     0,   511,   514,
     512,   515,   487,   573,   574,   578,     0,   564,   563,   567,
     568,   572,   452,   449,     0,     0,   436,   438,   439,   415,
     412,   410,   411,   413,   414,   409,   427,     0,   429,   407,
     408,   424,   425,     0,     0,     0,   430,     0,     0,   517,
       0,     0,   506,   509,   513,   559,   326,     0,     0,     0,
       0,   406,     0,     0,     0,   459,     0,   524,   526,   525,
     527,   527,   527,     0,   383,     0,   440,   428,   426,   423,
     422,   433,     0,   523,   521,   522,   450,     0,   461,   534,
     536,   531,   533,   530,   535,   532,   529,   528,     0,     0,
     384,   441
};

  /* YYPGOTO[NTERM-NUM].  */
static const yytype_int16 yypgoto[] =
{
   -1009, -1009, -1009,  1227,  -720,     0,   -93,  -109,  -126,  -572,
    -762,   -71,  -479, -1009,   924,   925,   228, -1009,   709, -1009,
   -1009,  1166,  -132,  -102,  -136, -1009,   560,  -432,  -129,  -114,
   -1009,   -99,   -95,  -139, -1009, -1009, -1009, -1009, -1009, -1009,
     546,   -98,  -406, -1009, -1009, -1009, -1009, -1009, -1009, -1009,
   -1009,  1035,   424,   115, -1009, -1009,  1003,   511, -1009,  1098,
   -1009, -1009, -1009, -1009, -1009, -1009, -1009, -1009, -1009,   -52,
   -1009, -1009, -1009, -1009, -1009, -1009, -1009, -1009,  -475, -1009,
   -1009, -1009, -1009,   -17, -1009,   429, -1009, -1009, -1009, -1009,
   -1009, -1009,   797, -1009, -1009,  -633, -1009, -1009,  1096,  -307,
    -151, -1009, -1009, -1009, -1009, -1009, -1009, -1009, -1009,   939,
   -1009, -1009, -1009,   554, -1009, -1009, -1009,   370,  -195, -1009,
   -1009, -1009, -1009, -1009, -1009, -1009, -1009,  -262,   -91,  -138,
    -707,  -600, -1009, -1009,  1201, -1009,   867, -1009, -1009, -1009,
   -1009, -1009, -1009,  -590, -1009, -1009, -1009,   731, -1009,  -525,
   -1009, -1009, -1009, -1009, -1009, -1009, -1009,   484, -1009, -1009,
   -1009,  1004,   590, -1009, -1009, -1009,   477,   277, -1009, -1009,
   -1009, -1009, -1009,  -779,  -968, -1009, -1009, -1009,  -442,   202,
   -1009, -1009, -1009, -1009, -1009, -1009,   281, -1009, -1009, -1009,
   -1009, -1009,   504, -1009, -1009, -1009, -1009, -1009, -1009, -1009,
    1110, -1009, -1009, -1009,   822, -1008, -1009, -1009, -1009, -1009,
    1090, -1009, -1009, -1009, -1009, -1009, -1009, -1009, -1009, -1009,
     650, -1009, -1009, -1009, -1009, -1009, -1009,   382,  -612, -1009,
   -1009, -1009, -1009, -1009, -1009, -1009,   325, -1009, -1009, -1009,
    -446,  -465, -1009, -1009, -1009, -1009, -1009, -1009, -1009,   792,
   -1009, -1009, -1009, -1009,   552, -1009,   312, -1009, -1009, -1009,
   -1009, -1009, -1009,   376, -1009, -1009, -1009,   383,  -686, -1009,
   -1009, -1009,   935,   558, -1009, -1009, -1009, -1009, -1009
};

  /* YYDEFGOTO[NTERM-NUM].  */
static const yytype_int16 yydefgoto[] =
{
      -1,    47,    48,    49,   169,   202,   171,   249,   825,   907,
     908,   523,   372,   373,   374,   375,   376,   447,   448,    51,
      52,    53,   866,   718,   321,   322,   312,   204,   205,   867,
     206,   207,   273,   173,   174,    54,    55,    56,   716,    57,
     232,   274,   416,   687,    58,    59,    60,    61,    62,    63,
     271,   272,   524,   529,    64,   306,   307,   572,    65,   360,
     210,    66,    67,    68,    69,    70,    71,    72,   212,   102,
     103,   109,   365,   494,   648,   758,   215,   875,   216,    73,
      74,    75,   224,   110,   383,   500,   517,   673,   674,   675,
     779,   676,   784,   876,   878,   877,    76,   217,   218,   759,
     505,   506,   507,   508,   509,   510,   872,   219,   220,   221,
     381,   501,   659,   660,   764,   765,   766,   869,   870,    77,
     107,   844,   382,   770,    78,   118,    79,    80,    81,   325,
     724,   594,   725,   726,    82,   455,   456,   457,   917,    83,
      84,   458,   597,   826,    85,   461,   156,   613,   851,   614,
     615,   616,   617,   618,   619,   620,   840,   841,    86,    87,
     748,   460,   730,   731,   622,   853,   854,   855,   939,   940,
    1059,  1107,  1108,  1010,  1011,  1012,  1013,  1110,  1111,  1112,
    1014,  1015,  1016,  1017,   941,  1056,  1057,  1129,  1158,    88,
     733,   600,   831,   832,    89,  1050,  1126,   990,    90,   228,
     229,   386,   882,  1064,  1058,   683,   785,   786,    91,   251,
     252,   522,    92,   419,   280,   551,   552,   553,   554,   694,
     695,   696,   792,   887,   697,   698,   896,   897,   898,   899,
     959,   962,  1027,  1082,  1069,  1070,  1071,  1072,  1073,  1074,
    1140,  1153,  1167,   972,    93,   287,   559,   422,   423,   560,
     561,   562,   563,   800,   801,   802,  1087,   979,  1040,  1041,
    1091,   803,   980,   981,  1085,   804,   976,   977,   978,    94,
     289,   426,   427,   708,   709,   568,   712,   809,   914
};

  /* YYTABLE[YYPACT[STATE-NUM]] -- What to do in state STATE-NUM.  If
     positive, shift that token.  If negative, reduce the rule whose
     number is the opposite.  If YYTABLE_NINF, syntax error.  */
static const yytype_int16 yytable[] =
{
      50,   250,   172,   227,   172,   317,   208,   101,   311,   324,
     106,   108,   310,   744,   326,   313,   973,   327,   328,   329,
     717,   650,   235,   973,   164,   165,   729,   237,  1083,   654,
     314,   656,   275,   569,  1036,   658,   666,   847,   880,  1089,
    -119,   817,  1063,   166,   320,   315,   557,   527,    50,   316,
     241,   973,   213,   919,   936,   248,  1115,   623,   511,   167,
     570,   818,   691,   453,   290,   353,   797,   797,   700,   720,
     667,   721,    43,    43,   504,  1159,   668,   669,  1160,  1161,
     502,   701,   589,   566,   762,   590,   291,    43,   745,    43,
     763,  1162,  1163,   158,   225,   503,   652,   170,   691,   170,
     308,   798,   798,   292,   101,   214,   722,   558,   425,   819,
     653,   108,   226,   680,   231,   231,   231,   170,   236,   227,
     925,   889,  -390,   366,   655,   250,   928,   226,   827,   170,
     170,    43,   548,   833,   657,   282,  1164,   354,   806,   293,
     294,  1063,   158,   879,   309,    95,   319,   723,   983,   699,
    1169,   512,   890,  1145,   782,  1165,  1166,   295,   296,   729,
     367,   159,   911,   974,   323,   837,   692,   668,   669,   168,
     974,   160,   161,   975,   297,   298,   299,   820,  1063,   183,
     975,   799,   799,   668,   669,   513,   300,   937,   301,   302,
     668,   669,   774,   453,    43,   702,   651,   821,   974,   822,
     353,   147,   692,   745,   303,   304,    96,   385,   975,   783,
     159,   -59,   101,   670,    43,   571,   982,   214,    43,   148,
     160,   161,    43,   172,   108,   450,   514,   515,   226,   168,
     903,    43,   463,   921,   454,   549,   464,   705,   693,   149,
     916,   557,   918,   953,   920,   955,   308,   550,   162,   253,
      43,   226,   828,   481,    43,   938,   192,    43,    43,  1038,
     671,   782,   254,   838,   839,   912,   823,   930,   782,   516,
     163,   255,   354,   147,   693,   793,  1033,   913,   308,   305,
    1109,   963,   964,  1042,  1051,   308,    43,   658,  1037,    97,
     891,   148,   519,   520,   982,   256,   104,   162,  1067,   168,
     672,  1067,   558,    43,   881,    43,   777,   892,   105,    43,
     379,   149,    43,   257,   258,   444,  1135,   991,   170,   163,
    1137,   823,   555,  1171,    43,   452,   150,   380,   528,   564,
     170,   170,   994,   995,   996,   893,   778,   528,   652,   377,
     259,   528,   528,   384,   720,    98,   721,   894,   170,  1020,
    1146,   545,   653,   151,   253,   546,  1029,  1030,  1031,  1032,
     895,   152,    99,    43,   454,   839,   331,   254,   260,   355,
     356,   357,   358,   359,    43,   153,   255,   154,    43,   214,
     499,   722,   640,   100,   355,   356,   357,   358,   359,  1053,
     170,   170,   521,   261,    43,   950,    43,   387,   593,   111,
     256,   956,   262,   263,   264,   545,   686,   114,   394,   749,
      43,   265,   931,   684,  1055,   112,   313,   266,   257,   258,
     170,   155,  1001,    43,   170,   151,    43,   115,   952,   915,
     954,   685,   946,   152,   113,   947,  1156,   290,   308,   355,
     356,   357,   358,   359,   926,   259,   315,   153,   116,   154,
     316,  1054,   267,   713,   308,   452,   158,  1096,   127,   291,
     934,   895,   895,   117,   355,   356,   357,   358,   359,   268,
     965,   545,   966,   260,   269,   781,   292,   270,   355,   356,
     357,   358,   359,   601,   719,   119,   355,   356,   357,   358,
     359,   545,   120,   155,  1117,   795,   715,  1024,   261,   142,
    1025,   602,   414,   760,   214,   761,   214,   262,   263,   264,
     214,  1149,   293,   294,  1132,  1026,   265,  1130,   678,   357,
     358,   359,   266,   539,   159,   121,  1143,   543,   544,   895,
     295,   296,   603,   895,   160,   161,   895,   895,   895,   895,
    1004,  1005,  1006,  1150,   122,   309,  1132,   297,   298,   299,
     355,   356,   357,   358,   359,   123,   604,   267,   124,   300,
     125,   301,   302,   126,  1151,   707,  1007,   128,  1008,  1009,
     968,   969,   970,   971,   268,   129,   605,   303,   304,   269,
     606,   813,   270,   131,  1022,   130,   787,   132,   789,   363,
     364,   170,   319,   226,   607,   246,   492,   493,   319,   355,
     356,   357,   358,   359,   133,   922,   923,   924,   796,   135,
     849,   162,   134,   483,   922,   923,   924,   313,   824,   993,
     368,   369,   370,   371,   868,   136,   247,   137,   608,   609,
     830,   138,   850,   163,   313,   355,   356,   357,   358,   359,
     610,   611,   369,   370,   144,   852,   145,   315,   146,   624,
     970,   971,   305,   355,   356,   357,   358,   359,   203,   157,
     209,   233,   234,   139,   315,   649,   140,   214,   316,   944,
     945,  1133,  1134,   230,  1141,  1142,  1154,  1155,   313,   239,
     240,   141,   242,   243,   244,   245,   313,   612,   319,   211,
     319,   277,   278,   909,   281,   886,   284,   285,   286,   222,
     223,   909,   288,   170,   158,   279,   318,   323,   315,   331,
     319,   330,   166,   332,   333,   334,   315,   335,   336,   885,
     888,   337,   338,   339,   390,   868,   868,   226,   167,   340,
     830,   319,   319,   341,   342,   313,   313,   319,   214,   214,
     389,   843,   343,   344,   345,   346,   848,   319,   347,   348,
     349,   852,   391,   355,   356,   357,   358,   359,   150,   350,
     351,   352,   214,   309,   392,   315,   315,   625,   393,   316,
     316,   396,   159,   397,   884,   874,   874,   398,   399,   214,
     400,   313,   160,   161,   355,   356,   357,   358,   359,   417,
     170,   401,   968,   969,   970,   971,   909,   402,   626,   355,
     356,   357,   358,   359,   403,   707,  1076,   906,   868,   404,
     405,   315,   319,   627,   319,   906,   319,   421,   313,   406,
     407,   525,   526,   226,   530,   531,   533,   534,   535,   536,
     408,   538,   319,   540,   541,   542,   409,   410,   168,   319,
     214,   355,   356,   357,   358,   359,   411,   412,   315,   413,
     418,   424,   316,   319,   420,   628,   425,   428,   429,   162,
     430,   431,   432,    43,   309,   309,   433,   434,   449,   435,
     436,   948,  1028,   437,   438,   439,   874,  -104,   874,   440,
     441,   163,  1045,   499,  1047,  -110,   442,   -73,   -73,   -73,
     -73,   -73,   465,   466,   467,   468,   469,   470,   471,   472,
     473,   474,   475,   476,   477,   478,   479,   480,   443,   482,
     906,   451,   459,   490,   491,   484,   485,   486,   487,   488,
     518,   489,   226,   226,   226,   355,   356,   357,   358,   359,
     556,   495,   172,  1084,  1081,   565,   591,   592,  -390,   629,
     573,   574,  1090,   576,   596,   578,   579,   309,   581,   582,
     583,   584,   585,   586,   587,   355,   356,   357,   358,   359,
     598,   621,   532,   599,   359,   645,   537,   170,   647,   630,
     646,   355,   356,   357,   358,   359,   502,   661,   547,  1039,
     662,   663,   664,   665,   319,   631,   319,   688,   682,   689,
     681,   575,  1139,   577,  1139,  1139,   580,   711,   690,   678,
     355,   356,   357,   358,   359,   704,   588,   703,   710,  1138,
     727,  1138,  1138,   734,   632,   732,   735,   737,  1068,   757,
     736,  1068,   355,   356,   357,   358,   359,   170,   355,   356,
     357,   358,   359,   738,   739,   319,   633,   355,   356,   357,
     358,   359,   638,   740,   319,   355,   356,   357,   358,   359,
     643,   644,   741,   742,   743,   746,   678,   175,   747,   639,
     176,   652,   775,   776,   780,  1116,   319,   355,   356,   357,
     358,   359,   790,   177,   791,   178,   355,   356,   357,   358,
     359,   641,   805,   807,   179,   355,   356,   357,   358,   359,
     856,   355,   356,   357,   358,   359,   180,   811,   808,   857,
     355,   356,   357,   358,   359,   858,   812,   355,   356,   357,
     358,   359,   814,   815,   859,   816,   836,   846,   226,  1144,
     226,   226,   355,   356,   357,   358,   359,   181,   355,   356,
     357,   358,   359,   861,  1170,   182,   319,   166,   634,   355,
     356,   357,   358,   359,   355,   356,   357,   358,   359,   635,
     449,   871,   862,   167,   636,   183,   355,   356,   357,   358,
     359,   -67,   -67,   -67,   -67,   -67,   637,   873,   184,   185,
     186,   355,   356,   357,   358,   359,   883,   900,  1048,   901,
     187,   910,   960,   902,   188,   355,   356,   357,   358,   359,
     929,   933,  1075,   951,   750,   751,   752,   753,   957,   961,
     984,   985,   754,   986,   987,   755,   756,   361,   355,   356,
     357,   358,   359,   989,  1035,   992,  1018,   189,   999,  1021,
    1000,   767,   769,   771,   772,   773,  1049,   190,   191,  1044,
     924,  1060,   192,  1061,  1062,  1065,  1066,   175,  1078,  1086,
     176,  1079,  1080,  1099,  1094,   193,  1095,   194,  1122,   788,
     195,  1098,  1113,   177,  1114,   178,  1100,  1123,  1124,   196,
    1101,  1118,   197,   168,   179,   371,  1120,  1121,  1157,   198,
    1131,  1132,  1136,   810,  1152,   143,   180,   642,   355,   356,
     357,   358,   359,   -66,   -66,   -66,   -66,   -66,    43,   -73,
     -73,   -73,   -73,   -73,   834,   835,  1119,   496,   497,   714,
     842,   283,   845,  1004,  1005,  1006,   415,   181,   362,   445,
     199,   200,   958,   378,   677,   182,   863,  1002,   498,   238,
    1102,   829,   595,   201,   932,   864,   865,  1093,   728,  1007,
     935,  1008,  1009,   462,  1148,   183,   927,  1097,   388,   175,
     679,   395,   176,  1103,   794,  1099,  1077,  1023,   184,   185,
     186,   706,   904,  1088,     0,   177,  1043,   178,  1100,  1034,
     187,   567,  1101,   905,   188,     0,   179,   860,   355,   356,
     357,   358,   359,   967,   968,   969,   970,   971,   180,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,   189,     0,     0,
       0,     0,     0,     0,     0,     0,     0,   190,   191,   181,
       0,     0,   192,     0,     0,     0,     0,   182,     0,     0,
     942,   943,  1102,  1104,     0,   193,     0,   194,     0,     0,
     195,     0,     0,   949,     0,     0,     0,   183,     0,   196,
       0,     0,   197,     0,     0,  1103,     0,     0,     0,   198,
     184,   185,   186,     0,     0,     0,     0,     0,  1105,     0,
       0,     0,   187,     0,     0,     0,   188,     0,    43,     0,
       0,     0,     0,     0,     0,   988,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,   997,
     199,   200,     0,   998,     0,     0,     0,     0,     0,   189,
       0,   175,     0,   201,   176,     0,     0,     0,     0,   190,
     191,  1003,     0,     0,   192,     0,     0,   177,     0,   178,
    1019,     0,     0,     0,     0,  1104,     0,   193,   179,   194,
       0,     0,   195,     0,     0,     0,     0,     0,     0,     0,
     180,   196,     0,     0,   197,  1046,     0,     0,     0,     0,
       0,   198,  1052,     0,     0,     0,     0,     0,     0,     0,
    1105,     0,     0,     0,     0,     0,     0,     0,     0,     0,
      43,   181,     0,   158,     0,     0,     0,     0,     0,   182,
       0,     0,     0,     0,     0,     0,     0,     0,   175,     0,
       0,   176,   199,   200,     0,     0,     0,     0,     0,   183,
       0,     0,     0,     0,   177,   201,   178,     0,     0,  1092,
       0,     0,   184,   185,   186,   179,     0,     0,     0,  1106,
       0,     0,     0,     0,   187,     0,     0,   180,   188,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,   159,     0,     0,     0,     0,  1125,     0,     0,     0,
       0,   160,   161,     0,  1127,  1128,     0,     0,   181,     0,
       0,   189,     0,     0,     0,     0,   182,     0,   166,     0,
       0,   190,   191,  1106,  1106,     0,   192,     0,     0,     0,
       0,     0,     0,     0,   167,     0,   183,     0,     0,   193,
    1147,   194,  1106,     0,   195,     0,     0,     0,     0,   184,
     185,   186,     0,   196,     0,     0,   197,     0,     0,     0,
       0,   187,     0,   198,     0,   188,     0,  1168,     0,     0,
       0,     0,     0,     0,     0,     0,   175,     0,   162,   176,
       0,     0,    43,     0,     0,     0,     0,     0,     0,     0,
       0,     0,   177,     0,   178,     0,     0,     0,   189,     0,
     163,     0,     0,   179,   199,   200,     0,     0,   190,   191,
       0,     0,     0,   192,     0,   180,     0,   201,     0,     0,
       0,     0,     0,     0,     0,     0,   193,     0,   194,     0,
       0,   195,     0,     0,     0,     0,     0,     0,     0,     0,
     196,     0,     0,   197,   168,     0,   181,     0,     0,     0,
     198,     0,     0,     0,   182,     0,   166,     0,     0,     0,
       0,     0,     0,   175,     0,     0,   176,     0,     0,    43,
       0,     0,   167,     0,   183,     0,     0,     0,     0,   177,
       0,   178,     0,     0,     0,     0,     0,   184,   185,   186,
     179,   199,   200,     0,     0,     0,     0,   308,     0,   187,
       0,     0,   180,   188,   201,     0,     0,     0,     0,     0,
       0,     0,     0,   175,     0,     0,   176,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,   177,
       0,   178,     0,   181,     0,     0,   189,     0,     0,     0,
     179,   182,     0,     0,     0,     0,   190,   191,     0,     0,
       0,   192,   180,     0,     0,     0,     0,     0,     0,     0,
       0,   183,   276,     0,   193,     0,   194,     0,     0,   195,
       0,     0,     0,     0,   184,   185,   186,     0,   196,     0,
       0,   197,   168,   181,     0,     0,   187,     0,   198,     0,
     188,   182,     0,     0,     0,     0,     0,     0,     0,     0,
     175,     0,     0,   176,     0,     0,     0,    43,     0,     0,
       0,   183,     0,     0,     0,     0,   177,     0,   178,     0,
       0,     0,     0,   189,   184,   185,   186,   179,     0,   199,
     200,     0,     0,   190,   191,     0,   187,     0,   192,   180,
     188,     0,   201,     0,     0,     0,     0,     0,     0,     0,
       0,   193,     0,   194,     0,     0,   195,     0,     0,     0,
       0,     0,     0,     0,     0,   196,     0,     0,   197,     0,
     181,     0,     0,   189,     0,   198,     0,     0,   182,     0,
       0,     0,     0,   190,   191,     0,     0,     0,   192,     0,
       0,     0,     0,     0,    43,     0,     0,     0,   183,     0,
       0,   193,     0,   194,     0,     0,   195,     0,     0,     0,
       0,   184,   185,   186,     0,   196,   199,   200,   197,     0,
       0,     0,     0,   187,     0,   198,     0,   188,     0,   201,
       0,     0,     0,     0,     0,     0,     0,   175,     0,     0,
     176,     0,     0,     0,    43,     0,     0,     0,     0,     0,
       0,     0,     0,   177,     0,   178,     0,     0,     0,     0,
     189,     0,     0,     0,   179,     0,   199,   200,     0,     0,
     190,   191,   446,     0,     0,   192,   180,     0,     0,   201,
       0,     0,     0,     0,     0,     0,     0,     0,   193,     0,
     194,     0,     0,   195,     0,     0,     0,     0,     0,     0,
       0,     0,   196,     0,     0,   197,     0,   181,     0,     0,
       0,     0,   198,     0,     0,   182,     0,     0,     0,     0,
       0,     0,     0,     0,   175,     0,     0,   176,   768,     0,
       0,    43,     0,     0,     0,   183,     0,     0,     0,     0,
     177,     0,   178,     0,     0,     0,     0,     0,   184,   185,
     186,   179,     0,   199,   200,     0,     0,     0,     0,     0,
     187,     0,     0,   180,   188,     0,   201,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,   181,     0,     0,   189,     0,     0,
       0,     0,   182,     0,     0,     0,     0,   190,   191,     0,
       0,     0,   192,     0,     0,     0,     0,     0,     0,     0,
       0,     0,   183,     0,     0,   193,     0,   194,     0,     0,
     195,     0,     0,     0,     0,   184,   185,   186,     0,   196,
       0,     0,   197,     0,     0,     0,     0,   187,     0,   198,
       0,   188,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,    43,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,   189,     0,     0,     0,     0,     0,
     199,   200,     0,     1,   190,   191,   308,     0,     0,   192,
       0,     0,     0,   201,     0,     0,     0,     0,     0,     0,
       0,     0,   193,     0,   194,     0,     0,   195,     0,     0,
       2,     3,     4,     5,     6,     0,   196,     0,     0,   197,
       0,     0,     0,     0,     0,     0,   198,     7,     8,     9,
      10,    11,    12,    13,     0,     0,     0,     0,     0,     0,
      14,     0,     0,     0,     0,    43,     0,     0,     0,     0,
       0,     0,     0,     0,    15,     0,     0,     0,     0,     0,
       0,     0,    16,    17,     0,     0,     0,   199,   200,     0,
       0,     0,     0,     0,    18,     0,     0,     0,    19,     0,
     201,    20,     0,     0,     0,    21,    22,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,    23,    24,
      25,    26,     0,    27,     0,     0,     0,     0,     0,     0,
      28,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,    29,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,    30,    31,
      32,     0,     0,     0,    33,    34,     0,     0,     0,    35,
      36,     0,     0,     0,    37,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,    38,     0,     0,
       0,     0,    39,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,    40,
      41,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,    42,    43,     0,     0,     0,     0,    44,     0,     0,
      45,     0,     0,     0,     0,     0,     0,     0,    46
};

static const yytype_int16 yycheck[] =
{
       0,   127,    97,   112,    99,   144,    99,     7,   144,   147,
      10,    11,   144,   613,   152,   144,     8,   155,   156,   157,
     592,   500,   117,     8,    95,    96,   598,   118,    12,   504,
     144,   506,   130,    67,    14,   510,     4,   744,    14,    16,
      82,    44,  1010,    83,   146,   144,    77,    29,    48,   144,
     121,     8,   104,   815,    55,   126,  1064,   463,    57,    99,
      94,    64,    91,   325,    17,    72,   105,   105,    18,   190,
      38,   192,   234,   234,   381,     7,   135,   136,    10,    11,
     257,    31,   263,   113,   256,   266,    39,   234,   613,   234,
     262,    23,    24,    75,   111,   272,   257,    97,    91,    99,
     262,   140,   140,    56,   104,   105,   227,   138,   138,   112,
     271,   111,   112,   519,   114,   115,   116,   117,   118,   228,
     827,    43,   164,   268,   271,   251,   833,   127,   728,   129,
     130,   234,   104,   733,   146,   135,    68,   144,   710,    92,
      93,  1109,    75,   776,   144,   253,   146,   268,   910,   555,
    1158,   150,    74,    40,   213,    87,    88,   110,   111,   731,
     212,   143,   104,   155,   267,   737,   195,   135,   136,   209,
     155,   153,   154,   165,   127,   128,   129,   180,  1146,   101,
     165,   220,   220,   135,   136,   184,   139,   188,   141,   142,
     135,   136,   667,   455,   234,   145,   503,   200,   155,   202,
      72,    27,   195,   728,   157,   158,   253,   224,   165,   268,
     143,   253,   212,   181,   234,   249,   902,   217,   234,    45,
     153,   154,   234,   318,   224,   318,   225,   226,   228,   209,
     268,   234,   330,   823,   325,   207,   331,   268,   267,    65,
     812,    77,   814,   876,   816,   878,   262,   219,   230,    15,
     234,   251,   268,   348,   234,   855,   178,   234,   234,   979,
     228,   213,    28,   738,   739,   207,   269,   839,   213,   268,
     252,    37,   144,    27,   267,   268,   268,   219,   262,   232,
    1059,   893,   894,   268,   991,   262,   234,   762,   268,   253,
     212,    45,   390,   391,   980,    61,   267,   230,  1018,   209,
     268,  1021,   138,   234,   779,   234,   223,   229,   267,   234,
     256,    65,   234,    79,    80,   268,   268,   917,   318,   252,
     251,   269,   420,   268,   234,   325,   152,   273,   399,   424,
     330,   331,   922,   923,   924,   257,   253,   408,   257,   268,
     106,   412,   413,   268,   190,   253,   192,   269,   348,   961,
    1129,   266,   271,   179,    15,   270,   968,   969,   970,   971,
     792,   187,   253,   234,   455,   840,   254,    28,   134,   256,
     257,   258,   259,   260,   234,   201,    37,   203,   234,   379,
     380,   227,   270,   253,   256,   257,   258,   259,   260,   206,
     390,   391,   392,   159,   234,   874,   234,   268,   152,   267,
      61,   880,   168,   169,   170,   266,   545,   269,   268,   270,
     234,   177,   268,   545,   206,   267,   545,   183,    79,    80,
     420,   247,   268,   234,   424,   179,   234,   269,   268,   207,
     268,   545,   263,   187,   253,   266,  1143,    17,   262,   256,
     257,   258,   259,   260,   268,   106,   545,   201,   269,   203,
     545,   268,   218,   223,   262,   455,    75,   268,   267,    39,
     268,   893,   894,   253,   256,   257,   258,   259,   260,   235,
     266,   266,   268,   134,   240,   270,    56,   243,   256,   257,
     258,   259,   260,    40,   593,   253,   256,   257,   258,   259,
     260,   266,   253,   247,  1066,   270,   591,   253,   159,     0,
     256,    58,   268,   654,   504,   656,   506,   168,   169,   170,
     510,   263,    92,    93,   266,   271,   177,   223,   518,   258,
     259,   260,   183,   408,   143,   253,  1126,   412,   413,   961,
     110,   111,    89,   965,   153,   154,   968,   969,   970,   971,
      69,    70,    71,   263,   253,   545,   266,   127,   128,   129,
     256,   257,   258,   259,   260,   253,   113,   218,   253,   139,
     253,   141,   142,   253,  1136,   565,    95,   267,    97,    98,
     256,   257,   258,   259,   235,   269,   133,   157,   158,   240,
     137,   719,   243,   253,   270,   269,   688,   253,   690,    59,
      60,   591,   592,   593,   151,   214,    46,    47,   598,   256,
     257,   258,   259,   260,   267,   256,   257,   258,   703,   269,
     746,   230,   253,   270,   256,   257,   258,   746,   727,   270,
     264,   265,   266,   267,   763,   253,   245,   253,   185,   186,
     732,   253,   746,   252,   763,   256,   257,   258,   259,   260,
     197,   198,   265,   266,   253,   747,    82,   746,   253,   270,
     258,   259,   232,   256,   257,   258,   259,   260,    98,   164,
     100,   115,   116,   267,   763,   268,   267,   667,   763,   864,
     865,  1113,  1114,   113,  1120,  1121,  1141,  1142,   807,   119,
     120,   267,   122,   123,   124,   125,   815,   244,   688,   267,
     690,   131,   132,   807,   134,   790,   136,   137,   138,   267,
     267,   815,    78,   703,    75,   156,   269,   267,   807,   254,
     710,   269,    83,   269,   269,   269,   815,   269,   269,   790,
     791,   269,   269,   269,   266,   864,   865,   727,    99,   269,
     832,   731,   732,   269,   269,   864,   865,   737,   738,   739,
     270,   741,   269,   269,   269,   269,   746,   747,   269,   269,
     269,   853,   266,   256,   257,   258,   259,   260,   152,   199,
     200,   201,   762,   763,   262,   864,   865,   270,   267,   864,
     865,   253,   143,   253,   145,   775,   776,   253,   253,   779,
     253,   910,   153,   154,   256,   257,   258,   259,   260,   270,
     790,   253,   256,   257,   258,   259,   910,   253,   270,   256,
     257,   258,   259,   260,   253,   805,   270,   807,   947,   253,
     253,   910,   812,   270,   814,   815,   816,   131,   947,   253,
     253,   397,   398,   823,   400,   401,   402,   403,   404,   405,
     253,   407,   832,   409,   410,   411,   253,   253,   209,   839,
     840,   256,   257,   258,   259,   260,   253,   253,   947,   253,
     253,   253,   947,   853,   266,   270,   138,   253,   253,   230,
     253,   253,   253,   234,   864,   865,   253,   253,   308,   253,
     253,   871,   967,   253,   253,   253,   876,   254,   878,   253,
     253,   252,   984,   883,   986,   254,   253,   256,   257,   258,
     259,   260,   332,   333,   334,   335,   336,   337,   338,   339,
     340,   341,   342,   343,   344,   345,   346,   347,   253,   349,
     910,    41,   236,   253,   253,   355,   356,   357,   358,   359,
     262,   361,   922,   923,   924,   256,   257,   258,   259,   260,
     253,   371,  1027,  1035,  1027,   253,   266,   253,   164,   270,
     429,   430,  1044,   432,   194,   434,   435,   947,   437,   438,
     439,   440,   441,   442,   443,   256,   257,   258,   259,   260,
     267,    63,   402,   217,   260,   253,   406,   967,   216,   270,
     253,   256,   257,   258,   259,   260,   257,   253,   418,   979,
     253,   253,   253,   253,   984,   270,   986,   253,   263,   253,
     270,   431,  1118,   433,  1120,  1121,   436,   248,   253,   999,
     256,   257,   258,   259,   260,   267,   446,   253,   253,  1118,
     253,  1120,  1121,   253,   270,   267,   253,   253,  1018,   117,
     267,  1021,   256,   257,   258,   259,   260,  1027,   256,   257,
     258,   259,   260,   253,   267,  1035,   270,   256,   257,   258,
     259,   260,   270,   253,  1044,   256,   257,   258,   259,   260,
     490,   491,   253,   253,   267,   253,  1056,     3,   267,   270,
       6,   257,   267,   267,   263,  1065,  1066,   256,   257,   258,
     259,   260,   253,    19,   253,    21,   256,   257,   258,   259,
     260,   270,   256,   253,    30,   256,   257,   258,   259,   260,
     270,   256,   257,   258,   259,   260,    42,   270,   246,   270,
     256,   257,   258,   259,   260,   270,   266,   256,   257,   258,
     259,   260,   253,   253,   270,   253,   119,   118,  1118,   268,
    1120,  1121,   256,   257,   258,   259,   260,    73,   256,   257,
     258,   259,   260,   253,   268,    81,  1136,    83,   266,   256,
     257,   258,   259,   260,   256,   257,   258,   259,   260,   266,
     590,   255,   268,    99,   266,   101,   256,   257,   258,   259,
     260,   256,   257,   258,   259,   260,   266,   253,   114,   115,
     116,   256,   257,   258,   259,   260,   253,   267,   263,   267,
     126,   253,   262,   267,   130,   256,   257,   258,   259,   260,
     253,   253,   263,   253,   634,   635,   636,   637,   267,   269,
     253,   253,   642,   253,   268,   645,   646,   255,   256,   257,
     258,   259,   260,    20,   273,   266,   262,   163,   267,   262,
     267,   661,   662,   663,   664,   665,   253,   173,   174,   273,
     258,   263,   178,   253,   253,   253,   253,     3,   271,   253,
       6,   271,   271,     9,   253,   191,   253,   193,   263,   689,
     196,   267,   262,    19,   262,    21,    22,   263,   271,   205,
      26,   266,   208,   209,    30,   267,   266,   266,   253,   215,
     268,   266,   266,   713,   266,    48,    42,   255,   256,   257,
     258,   259,   260,   256,   257,   258,   259,   260,   234,   256,
     257,   258,   259,   260,   734,   735,  1068,   373,   373,   590,
     740,   135,   742,    69,    70,    71,   271,    73,   210,   306,
     256,   257,   883,   217,   517,    81,   762,   947,   379,   118,
      86,   731,   455,   269,   840,   271,   272,  1050,   597,    95,
     853,    97,    98,   329,  1132,   101,   832,  1056,   228,     3,
     518,   251,     6,   109,   694,     9,  1021,   965,   114,   115,
     116,   559,   800,  1041,    -1,    19,   980,    21,    22,   976,
     126,   426,    26,   805,   130,    -1,    30,   255,   256,   257,
     258,   259,   260,   255,   256,   257,   258,   259,    42,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,   163,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,   173,   174,    73,
      -1,    -1,   178,    -1,    -1,    -1,    -1,    81,    -1,    -1,
     860,   861,    86,   189,    -1,   191,    -1,   193,    -1,    -1,
     196,    -1,    -1,   873,    -1,    -1,    -1,   101,    -1,   205,
      -1,    -1,   208,    -1,    -1,   109,    -1,    -1,    -1,   215,
     114,   115,   116,    -1,    -1,    -1,    -1,    -1,   224,    -1,
      -1,    -1,   126,    -1,    -1,    -1,   130,    -1,   234,    -1,
      -1,    -1,    -1,    -1,    -1,   915,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   929,
     256,   257,    -1,   933,    -1,    -1,    -1,    -1,    -1,   163,
      -1,     3,    -1,   269,     6,    -1,    -1,    -1,    -1,   173,
     174,   951,    -1,    -1,   178,    -1,    -1,    19,    -1,    21,
     960,    -1,    -1,    -1,    -1,   189,    -1,   191,    30,   193,
      -1,    -1,   196,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      42,   205,    -1,    -1,   208,   985,    -1,    -1,    -1,    -1,
      -1,   215,   992,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
     224,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
     234,    73,    -1,    75,    -1,    -1,    -1,    -1,    -1,    81,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,     3,    -1,
      -1,     6,   256,   257,    -1,    -1,    -1,    -1,    -1,   101,
      -1,    -1,    -1,    -1,    19,   269,    21,    -1,    -1,  1049,
      -1,    -1,   114,   115,   116,    30,    -1,    -1,    -1,  1059,
      -1,    -1,    -1,    -1,   126,    -1,    -1,    42,   130,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,   143,    -1,    -1,    -1,    -1,  1086,    -1,    -1,    -1,
      -1,   153,   154,    -1,  1094,  1095,    -1,    -1,    73,    -1,
      -1,   163,    -1,    -1,    -1,    -1,    81,    -1,    83,    -1,
      -1,   173,   174,  1113,  1114,    -1,   178,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    99,    -1,   101,    -1,    -1,   191,
    1130,   193,  1132,    -1,   196,    -1,    -1,    -1,    -1,   114,
     115,   116,    -1,   205,    -1,    -1,   208,    -1,    -1,    -1,
      -1,   126,    -1,   215,    -1,   130,    -1,  1157,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,     3,    -1,   230,     6,
      -1,    -1,   234,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    19,    -1,    21,    -1,    -1,    -1,   163,    -1,
     252,    -1,    -1,    30,   256,   257,    -1,    -1,   173,   174,
      -1,    -1,    -1,   178,    -1,    42,    -1,   269,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,   191,    -1,   193,    -1,
      -1,   196,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
     205,    -1,    -1,   208,   209,    -1,    73,    -1,    -1,    -1,
     215,    -1,    -1,    -1,    81,    -1,    83,    -1,    -1,    -1,
      -1,    -1,    -1,     3,    -1,    -1,     6,    -1,    -1,   234,
      -1,    -1,    99,    -1,   101,    -1,    -1,    -1,    -1,    19,
      -1,    21,    -1,    -1,    -1,    -1,    -1,   114,   115,   116,
      30,   256,   257,    -1,    -1,    -1,    -1,   262,    -1,   126,
      -1,    -1,    42,   130,   269,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,     3,    -1,    -1,     6,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    19,
      -1,    21,    -1,    73,    -1,    -1,   163,    -1,    -1,    -1,
      30,    81,    -1,    -1,    -1,    -1,   173,   174,    -1,    -1,
      -1,   178,    42,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,   101,   102,    -1,   191,    -1,   193,    -1,    -1,   196,
      -1,    -1,    -1,    -1,   114,   115,   116,    -1,   205,    -1,
      -1,   208,   209,    73,    -1,    -1,   126,    -1,   215,    -1,
     130,    81,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
       3,    -1,    -1,     6,    -1,    -1,    -1,   234,    -1,    -1,
      -1,   101,    -1,    -1,    -1,    -1,    19,    -1,    21,    -1,
      -1,    -1,    -1,   163,   114,   115,   116,    30,    -1,   256,
     257,    -1,    -1,   173,   174,    -1,   126,    -1,   178,    42,
     130,    -1,   269,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,   191,    -1,   193,    -1,    -1,   196,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,   205,    -1,    -1,   208,    -1,
      73,    -1,    -1,   163,    -1,   215,    -1,    -1,    81,    -1,
      -1,    -1,    -1,   173,   174,    -1,    -1,    -1,   178,    -1,
      -1,    -1,    -1,    -1,   234,    -1,    -1,    -1,   101,    -1,
      -1,   191,    -1,   193,    -1,    -1,   196,    -1,    -1,    -1,
      -1,   114,   115,   116,    -1,   205,   256,   257,   208,    -1,
      -1,    -1,    -1,   126,    -1,   215,    -1,   130,    -1,   269,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,     3,    -1,    -1,
       6,    -1,    -1,    -1,   234,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    19,    -1,    21,    -1,    -1,    -1,    -1,
     163,    -1,    -1,    -1,    30,    -1,   256,   257,    -1,    -1,
     173,   174,   262,    -1,    -1,   178,    42,    -1,    -1,   269,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   191,    -1,
     193,    -1,    -1,   196,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,   205,    -1,    -1,   208,    -1,    73,    -1,    -1,
      -1,    -1,   215,    -1,    -1,    81,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,     3,    -1,    -1,     6,   231,    -1,
      -1,   234,    -1,    -1,    -1,   101,    -1,    -1,    -1,    -1,
      19,    -1,    21,    -1,    -1,    -1,    -1,    -1,   114,   115,
     116,    30,    -1,   256,   257,    -1,    -1,    -1,    -1,    -1,
     126,    -1,    -1,    42,   130,    -1,   269,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    73,    -1,    -1,   163,    -1,    -1,
      -1,    -1,    81,    -1,    -1,    -1,    -1,   173,   174,    -1,
      -1,    -1,   178,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,   101,    -1,    -1,   191,    -1,   193,    -1,    -1,
     196,    -1,    -1,    -1,    -1,   114,   115,   116,    -1,   205,
      -1,    -1,   208,    -1,    -1,    -1,    -1,   126,    -1,   215,
      -1,   130,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   234,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,   163,    -1,    -1,    -1,    -1,    -1,
     256,   257,    -1,     5,   173,   174,   262,    -1,    -1,   178,
      -1,    -1,    -1,   269,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,   191,    -1,   193,    -1,    -1,   196,    -1,    -1,
      32,    33,    34,    35,    36,    -1,   205,    -1,    -1,   208,
      -1,    -1,    -1,    -1,    -1,    -1,   215,    49,    50,    51,
      52,    53,    54,    55,    -1,    -1,    -1,    -1,    -1,    -1,
      62,    -1,    -1,    -1,    -1,   234,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    76,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    84,    85,    -1,    -1,    -1,   256,   257,    -1,
      -1,    -1,    -1,    -1,    96,    -1,    -1,    -1,   100,    -1,
     269,   103,    -1,    -1,    -1,   107,   108,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   120,   121,
     122,   123,    -1,   125,    -1,    -1,    -1,    -1,    -1,    -1,
     132,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,   148,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   160,   161,
     162,    -1,    -1,    -1,   166,   167,    -1,    -1,    -1,   171,
     172,    -1,    -1,    -1,   176,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,   199,    -1,    -1,
      -1,    -1,   204,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   221,
     222,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,   233,   234,    -1,    -1,    -1,    -1,   239,    -1,    -1,
     242,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   250
};

  /* YYSTOS[STATE-NUM] -- The (internal number of the) accessing
     symbol of state STATE-NUM.  */
static const yytype_uint16 yystos[] =
{
       0,     5,    32,    33,    34,    35,    36,    49,    50,    51,
      52,    53,    54,    55,    62,    76,    84,    85,    96,   100,
     103,   107,   108,   120,   121,   122,   123,   125,   132,   148,
     160,   161,   162,   166,   167,   171,   172,   176,   199,   204,
     221,   222,   233,   234,   239,   242,   250,   275,   276,   277,
     279,   293,   294,   295,   309,   310,   311,   313,   318,   319,
     320,   321,   322,   323,   328,   332,   335,   336,   337,   338,
     339,   340,   341,   353,   354,   355,   370,   393,   398,   400,
     401,   402,   408,   413,   414,   418,   432,   433,   463,   468,
     472,   482,   486,   518,   543,   253,   253,   253,   253,   253,
     253,   279,   343,   344,   267,   267,   279,   394,   279,   345,
     357,   267,   267,   253,   269,   269,   269,   253,   399,   253,
     253,   253,   253,   253,   253,   253,   253,   267,   267,   269,
     269,   253,   253,   267,   253,   269,   253,   253,   253,   267,
     267,   267,     0,   277,   253,    82,   253,    27,    45,    65,
     152,   179,   187,   201,   203,   247,   420,   164,    75,   143,
     153,   154,   230,   252,   285,   285,    83,    99,   209,   278,
     279,   280,   306,   307,   308,     3,     6,    19,    21,    30,
      42,    73,    81,   101,   114,   115,   116,   126,   130,   163,
     173,   174,   178,   191,   193,   196,   205,   208,   215,   256,
     257,   269,   279,   300,   301,   302,   304,   305,   280,   300,
     334,   267,   342,   343,   279,   350,   352,   371,   372,   381,
     382,   383,   267,   267,   356,   357,   279,   281,   473,   474,
     300,   279,   314,   314,   314,   306,   279,   402,   408,   300,
     300,   285,   300,   300,   300,   300,   214,   245,   285,   281,
     282,   483,   484,    15,    28,    37,    61,    79,    80,   106,
     134,   159,   168,   169,   170,   177,   183,   218,   235,   240,
     243,   324,   325,   306,   315,   315,   102,   300,   300,   156,
     488,   300,   279,   295,   300,   300,   300,   519,    78,   544,
      17,    39,    56,    92,    93,   110,   111,   127,   128,   129,
     139,   141,   142,   157,   158,   232,   329,   330,   262,   279,
     296,   298,   300,   302,   303,   305,   306,   307,   269,   279,
     297,   298,   299,   267,   403,   403,   403,   403,   403,   403,
     269,   254,   269,   269,   269,   269,   269,   269,   269,   269,
     269,   269,   269,   269,   269,   269,   269,   269,   269,   269,
     300,   300,   300,    72,   144,   256,   257,   258,   259,   260,
     333,   255,   333,    59,    60,   346,   268,   343,   264,   265,
     266,   267,   286,   287,   288,   289,   290,   268,   372,   256,
     273,   384,   396,   358,   268,   357,   475,   268,   474,   270,
     266,   266,   262,   267,   268,   484,   253,   253,   253,   253,
     253,   253,   253,   253,   253,   253,   253,   253,   253,   253,
     253,   253,   253,   253,   268,   325,   316,   270,   253,   487,
     266,   131,   521,   522,   253,   138,   545,   546,   253,   253,
     253,   253,   253,   253,   253,   253,   253,   253,   253,   253,
     253,   253,   253,   253,   268,   330,   262,   291,   292,   300,
     280,    41,   279,   401,   402,   409,   410,   411,   415,   236,
     435,   419,   435,   315,   306,   300,   300,   300,   300,   300,
     300,   300,   300,   300,   300,   300,   300,   300,   300,   300,
     300,   306,   300,   270,   300,   300,   300,   300,   300,   300,
     253,   253,    46,    47,   347,   300,   288,   289,   383,   279,
     359,   385,   257,   272,   373,   374,   375,   376,   377,   378,
     379,    57,   150,   184,   225,   226,   268,   360,   262,   315,
     315,   279,   485,   285,   326,   326,   326,    29,   285,   327,
     326,   326,   300,   326,   326,   326,   326,   300,   326,   327,
     326,   326,   326,   327,   327,   266,   270,   300,   104,   207,
     219,   489,   490,   491,   492,   315,   253,    77,   138,   520,
     523,   524,   525,   526,   306,   253,   113,   546,   549,    67,
      94,   249,   331,   331,   331,   300,   331,   300,   331,   331,
     300,   331,   331,   331,   331,   331,   331,   331,   300,   263,
     266,   266,   253,   152,   405,   410,   194,   416,   267,   217,
     465,    40,    58,    89,   113,   133,   137,   151,   185,   186,
     197,   198,   244,   421,   423,   424,   425,   426,   427,   428,
     429,    63,   438,   316,   270,   270,   270,   270,   270,   270,
     270,   270,   270,   270,   266,   266,   266,   266,   270,   270,
     270,   270,   255,   300,   300,   253,   253,   216,   348,   268,
     286,   373,   257,   271,   352,   271,   352,   146,   352,   386,
     387,   253,   253,   253,   253,   253,     4,    38,   135,   136,
     181,   228,   268,   361,   362,   363,   365,   366,   279,   478,
     316,   270,   263,   479,   296,   303,   307,   317,   253,   253,
     253,    91,   195,   267,   493,   494,   495,   498,   499,   316,
      18,    31,   145,   253,   267,   268,   523,   279,   547,   548,
     253,   248,   550,   223,   292,   306,   312,   283,   297,   281,
     190,   192,   227,   268,   404,   406,   407,   253,   421,   283,
     436,   437,   267,   464,   253,   253,   267,   253,   253,   267,
     253,   253,   253,   267,   405,   423,   253,   267,   434,   270,
     300,   300,   300,   300,   300,   300,   300,   117,   349,   373,
     374,   374,   256,   262,   388,   389,   390,   300,   231,   300,
     397,   300,   300,   300,   352,   267,   267,   223,   253,   364,
     263,   270,   213,   268,   366,   480,   481,   297,   300,   297,
     253,   253,   496,   268,   494,   270,   306,   105,   140,   220,
     527,   528,   529,   535,   539,   256,   283,   253,   246,   551,
     300,   270,   266,   403,   253,   253,   253,    44,    64,   112,
     180,   200,   202,   269,   281,   282,   417,   405,   268,   436,
     297,   466,   467,   405,   300,   300,   119,   283,   352,   352,
     430,   431,   300,   279,   395,   300,   118,   404,   279,   298,
     303,   422,   297,   439,   440,   441,   270,   270,   270,   270,
     255,   253,   268,   387,   271,   272,   296,   303,   307,   391,
     392,   255,   380,   253,   279,   351,   367,   369,   368,   369,
      14,   352,   476,   253,   145,   285,   306,   497,   285,    43,
      74,   212,   229,   257,   269,   301,   500,   501,   502,   503,
     267,   267,   267,   268,   528,   547,   279,   283,   284,   303,
     253,   104,   207,   219,   552,   207,   283,   412,   283,   284,
     283,   417,   256,   257,   258,   404,   268,   466,   404,   253,
     283,   268,   431,   253,   268,   440,    55,   188,   405,   442,
     443,   458,   300,   300,   392,   392,   263,   266,   279,   300,
     286,   253,   268,   369,   268,   369,   286,   267,   359,   504,
     262,   269,   505,   502,   502,   266,   268,   255,   256,   257,
     258,   259,   517,     8,   155,   165,   540,   541,   542,   531,
     536,   537,   542,   284,   253,   253,   253,   268,   300,    20,
     471,   405,   266,   270,   417,   417,   417,   300,   300,   267,
     267,   268,   391,   300,    69,    70,    71,    95,    97,    98,
     447,   448,   449,   450,   454,   455,   456,   457,   262,   300,
     502,   262,   270,   501,   253,   256,   271,   506,   306,   502,
     502,   502,   502,   268,   541,   273,    14,   268,   278,   279,
     532,   533,   268,   537,   273,   297,   300,   297,   263,   253,
     469,   404,   300,   206,   268,   206,   459,   460,   478,   444,
     263,   253,   253,   448,   477,   253,   253,   278,   279,   508,
     509,   510,   511,   512,   513,   263,   270,   510,   271,   271,
     271,   280,   507,    12,   297,   538,   253,   530,   530,    16,
     297,   534,   300,   441,   253,   253,   268,   460,   267,     9,
      22,    26,    86,   109,   189,   224,   300,   445,   446,   447,
     451,   452,   453,   262,   262,   479,   279,   283,   266,   290,
     266,   266,   263,   263,   271,   300,   470,   300,   300,   461,
     223,   268,   266,   452,   452,   268,   266,   251,   281,   282,
     514,   514,   514,   405,   268,    40,   447,   300,   453,   263,
     263,   283,   266,   515,   515,   515,   404,   253,   462,     7,
      10,    11,    23,    24,    68,    87,    88,   516,   300,   479,
     268,   268
};

  /* YYR1[YYN] -- Symbol number of symbol that rule YYN derives.  */
static const yytype_uint16 yyr1[] =
{
       0,   274,   275,   276,   276,   277,   277,   277,   277,   277,
     277,   277,   277,   277,   277,   277,   277,   277,   277,   277,
     277,   277,   277,   277,   277,   278,   279,   280,   281,   282,
     283,   284,   284,   285,   285,   285,   285,   285,   285,   286,
     286,   286,   286,   287,   287,   287,   287,   288,   289,   290,
     291,   291,   292,   292,   293,   294,   294,   294,   294,   295,
     296,   297,   297,   298,   299,   300,   300,   301,   301,   302,
     302,   303,   303,   304,   305,   305,   305,   305,   305,   305,
     305,   305,   305,   305,   305,   305,   305,   305,   305,   305,
     305,   305,   305,   305,   305,   305,   305,   305,   305,   305,
     305,   305,   305,   305,   306,   306,   307,   307,   307,   307,
     308,   309,   309,   309,   309,   309,   309,   309,   310,   311,
     312,   313,   314,   315,   316,   316,   317,   317,   317,   318,
     319,   320,   321,   322,   323,   324,   324,   325,   325,   325,
     325,   325,   325,   325,   325,   325,   325,   325,   325,   325,
     325,   325,   325,   325,   325,   325,   326,   327,   327,   328,
     329,   329,   330,   330,   330,   330,   330,   330,   330,   330,
     330,   330,   330,   330,   330,   330,   330,   330,   331,   331,
     331,   332,   332,   332,   332,   332,   333,   333,   333,   334,
     334,   334,   334,   335,   335,   335,   335,   335,   335,   335,
     335,   335,   335,   335,   335,   335,   335,   335,   336,   336,
     336,   336,   337,   338,   338,   338,   339,   339,   340,   341,
     342,   342,   343,   344,   345,   346,   346,   347,   347,   347,
     348,   348,   349,   349,   350,   351,   352,   353,   353,   354,
     355,   356,   356,   358,   357,   359,   360,   360,   361,   361,
     362,   362,   362,   363,   363,   363,   364,   364,   365,   366,
     366,   367,   367,   368,   368,   369,   370,   371,   371,   372,
     373,   373,   374,   375,   376,   377,   378,   379,   379,   379,
     379,   380,   380,   381,   382,   382,   383,   384,   384,   385,
     386,   386,   387,   387,   388,   388,   389,   390,   391,   391,
     391,   392,   392,   392,   393,   394,   395,   396,   396,   396,
     396,   396,   396,   397,   397,   399,   398,   400,   401,   401,
     401,   401,   401,   402,   403,   404,   405,   405,   406,   407,
     407,   407,   408,   409,   409,   410,   410,   412,   411,   413,
     413,   415,   414,   416,   416,   416,   416,   416,   416,   416,
     416,   417,   417,   417,   417,   417,   419,   418,   420,   420,
     420,   420,   420,   421,   421,   422,   423,   423,   423,   423,
     423,   423,   423,   423,   423,   424,   424,   425,   425,   425,
     425,   426,   426,   427,   428,   429,   429,   430,   430,   431,
     432,   434,   433,   435,   436,   437,   437,   438,   439,   439,
     440,   441,   441,   442,   442,   444,   443,   445,   445,   446,
     446,   446,   446,   446,   446,   446,   447,   447,   448,   448,
     448,   448,   449,   450,   451,   452,   452,   453,   453,   453,
     454,   455,   455,   456,   457,   457,   458,   459,   459,   461,
     462,   460,   464,   463,   465,   466,   467,   467,   469,   470,
     468,   471,   471,   472,   473,   473,   475,   476,   477,   474,
     478,   479,   479,   480,   480,   481,   482,   483,   483,   485,
     484,   487,   486,   488,   488,   489,   489,   489,   490,   491,
     492,   493,   493,   494,   494,   494,   496,   495,   497,   497,
     497,   498,   499,   500,   500,   501,   502,   502,   502,   502,
     502,   502,   502,   502,   502,   504,   503,   503,   505,   503,
     506,   506,   506,   506,   506,   507,   508,   509,   510,   510,
     510,   511,   512,   513,   514,   514,   514,   515,   515,   516,
     516,   516,   516,   516,   516,   516,   516,   517,   517,   519,
     518,   520,   520,   521,   521,   522,   522,   522,   523,   523,
     524,   525,   526,   527,   527,   528,   528,   528,   529,   530,
     530,   531,   531,   532,   532,   533,   533,   534,   534,   535,
     536,   536,   537,   538,   538,   539,   540,   540,   541,   542,
     542,   542,   543,   544,   545,   545,   546,   547,   548,   548,
     549,   550,   551,   552,   552,   552,   552
};

  /* YYR2[YYN] -- Number of symbols on the right hand side of rule YYN.  */
static const yytype_uint8 yyr2[] =
{
       0,     2,     1,     1,     2,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     4,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     0,
       1,     1,     1,     1,     1,     2,     2,     1,     1,     3,
       1,     3,     1,     7,     3,     3,     3,     3,     3,     1,
       1,     1,     1,     3,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     3,     4,     4,     4,     6,     6,
       6,     4,     4,     4,     4,     4,     4,     4,     4,     4,
       4,     6,     1,     1,     1,     1,     4,     3,     3,     3,
       3,     3,     2,     2,     1,     1,     1,     1,     3,     5,
       1,     1,     1,     1,     1,     1,     1,     1,     7,     1,
       1,     4,     1,     1,     0,     3,     1,     1,     1,     5,
       7,     7,     4,     6,     4,     1,     2,     3,     3,     3,
       3,     3,     3,     3,     3,     3,     3,     3,     3,     3,
       3,     3,     3,     3,     3,     3,     1,     1,     1,     4,
       1,     2,     3,     3,     3,     3,     3,     3,     3,     3,
       3,     3,     3,     3,     3,     3,     3,     3,     1,     1,
       1,     3,     3,     4,     3,     4,     0,     1,     1,     1,
       3,     5,     7,     3,     3,     3,     3,     3,     3,     3,
       3,     3,     3,     3,     3,     3,     3,     3,     3,     3,
       3,     3,     3,     1,     1,     1,     1,     1,     2,     4,
       1,     2,     7,     1,     1,     3,     3,     0,     3,     3,
       0,     1,     0,     3,     1,     2,     2,     1,     1,     2,
       4,     1,     2,     0,     5,     1,     0,     2,     1,     1,
       3,     4,     4,     1,     1,     1,     1,     1,     1,     4,
       4,     1,     2,     1,     2,     3,     4,     1,     2,     1,
       1,     2,     2,     2,     2,     3,     3,     1,     1,     1,
       1,     0,     2,     6,     1,     3,     1,     0,     2,     2,
       1,     3,     1,     1,     1,     1,     3,     5,     1,     2,
       2,     1,     1,     1,     5,     1,     1,     0,     4,     4,
       4,     4,     4,     1,     1,     0,     3,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     0,     2,     1,     3,
       3,     5,     6,     1,     2,     1,     1,     0,     7,     1,
       1,     0,     8,     3,     3,     3,     3,     3,     3,     3,
       3,     1,     3,     3,     3,     3,     0,     7,     1,     1,
       1,     1,     1,     1,     2,     1,     3,     3,     1,     3,
       3,     3,     3,     3,     4,     1,     1,     1,     1,     1,
       1,     3,     6,     9,    12,     3,     3,     1,     2,     2,
       1,     0,     9,     4,     1,     1,     2,     4,     1,     2,
       1,     0,     2,     1,     1,     0,     5,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     2,     1,     1,
       1,     1,     5,     5,     1,     1,     3,     1,     3,     1,
       3,     1,     1,     5,     1,     1,     4,     1,     2,     0,
       0,     7,     0,     8,     4,     1,     1,     2,     0,     0,
      14,     0,     3,     4,     1,     2,     0,     0,     0,    11,
       1,     0,     2,     1,     1,     3,     4,     1,     2,     0,
       5,     0,     7,     0,     3,     1,     1,     1,     3,     3,
       3,     1,     2,     1,     1,     1,     0,     6,     1,     1,
       1,     3,     3,     1,     3,     2,     1,     1,     3,     3,
       3,     3,     3,     2,     4,     0,     5,     4,     0,     5,
       1,     2,     2,     3,     2,     1,     1,     2,     1,     1,
       1,     4,     4,     4,     1,     1,     1,     0,     2,     1,
       1,     1,     1,     1,     1,     1,     1,     0,     2,     0,
       6,     1,     2,     0,     1,     3,     3,     3,     1,     1,
       1,     3,     4,     1,     2,     1,     1,     1,     4,     2,
       0,     2,     0,     2,     2,     1,     1,     1,     1,     4,
       1,     2,     3,     1,     1,     4,     1,     2,     3,     1,
       1,     1,     9,     3,     1,     2,     3,     1,     1,     3,
       3,     3,     3,     0,     3,     3,     3
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
#line 631 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.sym) = mdl_existing_object(parse_state, (yyvsp[0].str))); }
#line 3144 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 29:
#line 634 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.sym) = mdl_existing_region(parse_state, (yyvsp[-3].sym), (yyvsp[-1].str))); }
#line 3150 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 30:
#line 637 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.vec3) = mdl_point(parse_state, &(yyvsp[0].nlist))); }
#line 3156 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 32:
#line 641 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.vec3) = mdl_point_scalar((yyvsp[0].dbl))); }
#line 3162 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 33:
#line 644 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = 1; }
#line 3168 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 34:
#line 645 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = 0; }
#line 3174 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 35:
#line 646 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = 1; }
#line 3180 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 36:
#line 647 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = 0; }
#line 3186 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 37:
#line 648 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = 1; }
#line 3192 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 38:
#line 649 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = 0; }
#line 3198 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 39:
#line 652 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.mol_type).orient_set = 0; }
#line 3204 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 42:
#line 655 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.mol_type).orient_set = 1; (yyval.mol_type).orient = 0; }
#line 3210 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 43:
#line 659 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.mol_type).orient = 1; (yyval.mol_type).orient_set = 1; }
#line 3216 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 44:
#line 660 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.mol_type).orient = -1; (yyval.mol_type).orient_set = 1; }
#line 3222 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 45:
#line 661 "../src/../src/mdlparse.y" /* yacc.c:1646  */
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
#line 3237 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 46:
#line 671 "../src/../src/mdlparse.y" /* yacc.c:1646  */
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
#line 3252 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 49:
#line 689 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                          (yyval.mol_type).orient = (int) (yyvsp[-1].dbl);
                                                          (yyval.mol_type).orient_set = 1;
                                                          if ((yyval.mol_type).orient != (yyvsp[-1].dbl))
                                                          {
                                                            mdlerror(parse_state, "molecule orientation specified inside braces must be an integer between -32768 and 32767.");
                                                            return 1;
                                                          }
                                                      }
#line 3266 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 51:
#line 702 "../src/../src/mdlparse.y" /* yacc.c:1646  */
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
#line 3282 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 52:
#line 715 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mcell_generate_range_singleton(&(yyval.nlist), (yyvsp[0].dbl))); }
#line 3288 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 53:
#line 716 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_generate_range(parse_state, &(yyval.nlist), (yyvsp[-5].dbl), (yyvsp[-3].dbl), (yyvsp[-1].dbl))); }
#line 3294 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 54:
#line 722 "../src/../src/mdlparse.y" /* yacc.c:1646  */
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
#line 3316 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 55:
#line 745 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_assign_variable_double(parse_state, (yyvsp[-2].sym), (yyvsp[0].dbl))); }
#line 3322 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 56:
#line 746 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_assign_variable_string(parse_state, (yyvsp[-2].sym), (yyvsp[0].str))); }
#line 3328 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 57:
#line 747 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_assign_variable(parse_state, (yyvsp[-2].sym), (yyvsp[0].sym))); }
#line 3334 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 58:
#line 748 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_assign_variable_array(parse_state, (yyvsp[-2].sym), (yyvsp[0].nlist).value_head)); }
#line 3340 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 59:
#line 751 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.sym) = mdl_get_or_create_variable(parse_state, (yyvsp[0].str))); }
#line 3346 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 60:
#line 754 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.sym) = mdl_existing_variable(parse_state, (yyvsp[0].str))); }
#line 3352 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 62:
#line 758 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                          struct num_expr_list *elp;
                                                          (yyval.nlist).value_head = (struct num_expr_list *) (yyvsp[0].sym)->value;
                                                          (yyval.nlist).value_count = 1;
                                                          for (elp = (yyval.nlist).value_head; elp->next != NULL; elp = elp->next)
                                                            ++ (yyval.nlist).value_count;
                                                          (yyval.nlist).value_tail = elp;
                                                          (yyval.nlist).shared = 1;
                                                      }
#line 3366 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 63:
#line 769 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { mdl_debug_dump_array((yyvsp[-1].nlist).value_head); (yyval.nlist) = (yyvsp[-1].nlist); }
#line 3372 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 64:
#line 772 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.sym) = mdl_existing_array(parse_state, (yyvsp[0].str))); }
#line 3378 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 68:
#line 780 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.dbl) = *(double *) (yyvsp[0].sym)->value; }
#line 3384 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 69:
#line 783 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.dbl) = (yyvsp[0].llival); }
#line 3390 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 73:
#line 791 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.sym) = mdl_existing_double(parse_state, (yyvsp[0].str))); }
#line 3396 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 74:
#line 795 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.dbl) = (yyvsp[-1].dbl); }
#line 3402 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 75:
#line 796 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKF((yyval.dbl) = exp((yyvsp[-1].dbl))); }
#line 3408 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 76:
#line 797 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_expr_log(parse_state, (yyvsp[-1].dbl), &(yyval.dbl))); }
#line 3414 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 77:
#line 798 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_expr_log10(parse_state, (yyvsp[-1].dbl), &(yyval.dbl))); }
#line 3420 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 78:
#line 799 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.dbl) = max2d((yyvsp[-3].dbl), (yyvsp[-1].dbl)); }
#line 3426 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 79:
#line 800 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.dbl) = min2d((yyvsp[-3].dbl), (yyvsp[-1].dbl)); }
#line 3432 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 80:
#line 801 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.dbl) = mdl_expr_roundoff((yyvsp[-1].dbl), (int) (yyvsp[-3].dbl)); }
#line 3438 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 81:
#line 802 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.dbl) = floor((yyvsp[-1].dbl)); }
#line 3444 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 82:
#line 803 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.dbl) = ceil((yyvsp[-1].dbl)); }
#line 3450 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 83:
#line 804 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.dbl) = sin((yyvsp[-1].dbl)); }
#line 3456 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 84:
#line 805 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.dbl) = cos((yyvsp[-1].dbl)); }
#line 3462 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 85:
#line 806 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKF((yyval.dbl) = tan((yyvsp[-1].dbl))); }
#line 3468 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 86:
#line 807 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKF((yyval.dbl) = asin((yyvsp[-1].dbl))); }
#line 3474 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 87:
#line 808 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKF((yyval.dbl) = acos((yyvsp[-1].dbl))); }
#line 3480 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 88:
#line 809 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.dbl) = atan((yyvsp[-1].dbl)); }
#line 3486 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 89:
#line 810 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKF((yyval.dbl) = sqrt((yyvsp[-1].dbl))); }
#line 3492 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 90:
#line 811 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.dbl) = fabs((yyvsp[-1].dbl)); }
#line 3498 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 91:
#line 812 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_expr_mod(parse_state, (yyvsp[-3].dbl), (yyvsp[-1].dbl), &(yyval.dbl))); }
#line 3504 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 92:
#line 813 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.dbl) = MY_PI; }
#line 3510 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 93:
#line 814 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.dbl) = mdl_expr_rng_uniform(parse_state); }
#line 3516 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 94:
#line 815 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.dbl) = rng_gauss(parse_state->vol->rng); }
#line 3522 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 95:
#line 816 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.dbl) = parse_state->vol->seed_seq; }
#line 3528 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 96:
#line 817 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_expr_string_to_double(parse_state, (yyvsp[-1].str), &(yyval.dbl))); }
#line 3534 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 97:
#line 818 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKF((yyval.dbl) = (yyvsp[-2].dbl) + (yyvsp[0].dbl)); }
#line 3540 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 98:
#line 819 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKF((yyval.dbl) = (yyvsp[-2].dbl) - (yyvsp[0].dbl)); }
#line 3546 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 99:
#line 820 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKF((yyval.dbl) = (yyvsp[-2].dbl) * (yyvsp[0].dbl)); }
#line 3552 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 100:
#line 821 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_expr_div(parse_state, (yyvsp[-2].dbl), (yyvsp[0].dbl), &(yyval.dbl))); }
#line 3558 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 101:
#line 822 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_expr_pow(parse_state, (yyvsp[-2].dbl), (yyvsp[0].dbl), &(yyval.dbl))); }
#line 3564 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 102:
#line 823 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.dbl) = -(yyvsp[0].dbl); }
#line 3570 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 103:
#line 824 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.dbl) = (yyvsp[0].dbl); }
#line 3576 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 105:
#line 829 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.str) = mdl_strdup((char const *) (yyvsp[0].sym)->value)); }
#line 3582 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 106:
#line 833 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.str) = mdl_strip_quotes((yyvsp[0].str))); }
#line 3588 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 107:
#line 834 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.str) = mdl_strdup(parse_state->vol->mdl_infile_name)); }
#line 3594 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 108:
#line 835 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.str) = mdl_strcat((yyvsp[-2].str), (yyvsp[0].str))); }
#line 3600 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 109:
#line 836 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.str) = mdl_string_format(parse_state, (yyvsp[-2].str), (yyvsp[-1].printfargs).arg_head)); }
#line 3606 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 110:
#line 839 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.sym) = mdl_existing_string(parse_state, (yyvsp[0].str))); }
#line 3612 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 118:
#line 855 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_fopen(parse_state, (yyvsp[-6].sym), (yyvsp[-3].str), (yyvsp[-1].str))); }
#line 3618 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 119:
#line 858 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.sym) = mdl_new_filehandle(parse_state, (yyvsp[0].str))); }
#line 3624 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 120:
#line 861 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.str) = (yyvsp[0].str); CHECK(mdl_valid_file_mode(parse_state, (yyvsp[0].str))); }
#line 3630 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 121:
#line 864 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_fclose(parse_state, (yyvsp[-1].sym))); }
#line 3636 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 122:
#line 867 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.sym) = mdl_existing_file_stream(parse_state, (yyvsp[0].str))); }
#line 3642 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 123:
#line 870 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.str) = mdl_expand_string_escapes((yyvsp[0].str))); }
#line 3648 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 124:
#line 873 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.printfargs).arg_head = (yyval.printfargs).arg_tail = NULL; }
#line 3654 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 125:
#line 874 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                        (yyval.printfargs) = (yyvsp[-2].printfargs);
                                                        if ((yyval.printfargs).arg_tail)
                                                          (yyval.printfargs).arg_tail = (yyval.printfargs).arg_tail->next = (yyvsp[0].printfarg);
                                                        else
                                                          (yyval.printfargs).arg_tail = (yyval.printfargs).arg_head = (yyvsp[0].printfarg);
                                                        (yyvsp[0].printfarg)->next = NULL;
                                                      }
#line 3667 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 126:
#line 884 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.printfarg) = mdl_new_printf_arg_double((yyvsp[0].dbl))); }
#line 3673 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 127:
#line 885 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.printfarg) = mdl_new_printf_arg_string((yyvsp[0].str))); }
#line 3679 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 128:
#line 886 "../src/../src/mdlparse.y" /* yacc.c:1646  */
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
#line 3694 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 129:
#line 898 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_printf(parse_state, (yyvsp[-2].str), (yyvsp[-1].printfargs).arg_head)); }
#line 3700 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 130:
#line 904 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_fprintf(parse_state, (struct file_stream *) (yyvsp[-4].sym)->value, (yyvsp[-2].str), (yyvsp[-1].printfargs).arg_head)); }
#line 3706 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 131:
#line 910 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_sprintf(parse_state, (yyvsp[-4].sym), (yyvsp[-2].str), (yyvsp[-1].printfargs).arg_head)); }
#line 3712 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 132:
#line 913 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { mdl_print_time(parse_state, (yyvsp[-1].str)); }
#line 3718 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 133:
#line 919 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_fprint_time(parse_state, (yyvsp[-3].sym), (yyvsp[-1].str))); }
#line 3724 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 137:
#line 935 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { if (!parse_state->vol->quiet_flag) mdl_set_all_notifications(parse_state->vol, (yyvsp[0].tok)); }
#line 3730 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 138:
#line 936 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { if (!parse_state->vol->quiet_flag) parse_state->vol->notify->progress_report        = (yyvsp[0].tok); }
#line 3736 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 139:
#line 937 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { if (!parse_state->vol->quiet_flag) parse_state->vol->notify->diffusion_constants    = (yyvsp[0].tok); }
#line 3742 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 140:
#line 938 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { if (!parse_state->vol->quiet_flag) parse_state->vol->notify->reaction_probabilities = (yyvsp[0].tok); }
#line 3748 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 141:
#line 939 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { if (!parse_state->vol->quiet_flag) parse_state->vol->notify->time_varying_reactions = (yyvsp[0].tok); }
#line 3754 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 142:
#line 940 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { if (!parse_state->vol->quiet_flag) parse_state->vol->notify->reaction_prob_notify   = (yyvsp[0].dbl); }
#line 3760 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 143:
#line 941 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { if (!parse_state->vol->quiet_flag) parse_state->vol->notify->partition_location     = (yyvsp[0].tok); }
#line 3766 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 144:
#line 942 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { if (!parse_state->vol->quiet_flag) parse_state->vol->notify->box_triangulation      = (yyvsp[0].tok); }
#line 3772 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 145:
#line 943 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { if (!parse_state->vol->quiet_flag) parse_state->vol->notify->release_events         = (yyvsp[0].tok); }
#line 3778 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 146:
#line 944 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { if (!parse_state->vol->quiet_flag) parse_state->vol->notify->file_writes            = (yyvsp[0].tok); }
#line 3784 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 147:
#line 945 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { if (!parse_state->vol->quiet_flag) parse_state->vol->notify->final_summary          = (yyvsp[0].tok); }
#line 3790 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 148:
#line 946 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { if (!parse_state->vol->quiet_flag) parse_state->vol->notify->throughput_report      = (yyvsp[0].tok); }
#line 3796 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 149:
#line 947 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { if (!parse_state->vol->quiet_flag) parse_state->vol->notify->reaction_output_report = (yyvsp[0].tok); }
#line 3802 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 150:
#line 948 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { if (!parse_state->vol->quiet_flag) parse_state->vol->notify->volume_output_report   = (yyvsp[0].tok); }
#line 3808 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 151:
#line 949 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { if (!parse_state->vol->quiet_flag) parse_state->vol->notify->viz_output_report      = (yyvsp[0].tok); }
#line 3814 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 152:
#line 950 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { if (!parse_state->vol->quiet_flag) parse_state->vol->notify->checkpoint_report      = (yyvsp[0].tok); }
#line 3820 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 153:
#line 951 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                          if (!parse_state->vol->quiet_flag && parse_state->vol->log_freq == ULONG_MAX)
                                                            parse_state->vol->notify->iteration_report = (yyvsp[0].tok);
                                                      }
#line 3829 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 154:
#line 955 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { if (!parse_state->vol->quiet_flag) CHECK(mdl_set_iteration_report_freq(parse_state, (long long) (yyvsp[0].dbl))); }
#line 3835 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 155:
#line 956 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { if (!parse_state->vol->quiet_flag) parse_state->vol->notify->molecule_collision_report    = (yyvsp[0].tok); }
#line 3841 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 156:
#line 960 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = ((yyvsp[0].tok) ? NOTIFY_FULL : NOTIFY_NONE); }
#line 3847 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 157:
#line 964 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = ((yyvsp[0].tok) ? NOTIFY_FULL : NOTIFY_NONE); }
#line 3853 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 158:
#line 965 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = NOTIFY_BRIEF; }
#line 3859 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 162:
#line 981 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { mdl_set_all_warnings(parse_state->vol, (byte) (yyvsp[0].tok)); }
#line 3865 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 163:
#line 982 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->vol->notify->neg_diffusion = (byte)(yyvsp[0].tok); }
#line 3871 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 164:
#line 983 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->vol->notify->neg_reaction = (byte)(yyvsp[0].tok); }
#line 3877 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 165:
#line 984 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->vol->notify->high_reaction_prob = (byte)(yyvsp[0].tok); }
#line 3883 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 166:
#line 985 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->vol->notify->reaction_prob_warn = (yyvsp[0].dbl); }
#line 3889 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 167:
#line 986 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->vol->notify->close_partitions = (byte)(yyvsp[0].tok); }
#line 3895 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 168:
#line 987 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->vol->notify->degenerate_polys = (byte)(yyvsp[0].tok); }
#line 3901 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 169:
#line 988 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->vol->notify->overwritten_file = (byte)(yyvsp[0].tok); }
#line 3907 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 170:
#line 989 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->vol->notify->short_lifetime = (byte)(yyvsp[0].tok); }
#line 3913 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 171:
#line 990 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_lifetime_warning_threshold(parse_state, (long long) (yyvsp[0].dbl))); }
#line 3919 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 172:
#line 991 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->vol->notify->missed_reactions = (byte)(yyvsp[0].tok); }
#line 3925 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 173:
#line 992 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_missed_reaction_warning_threshold(parse_state, (yyvsp[0].dbl))); }
#line 3931 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 174:
#line 993 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->vol->notify->missed_surf_orient = (byte)(yyvsp[0].tok); }
#line 3937 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 175:
#line 994 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->vol->notify->useless_vol_orient = (byte)(yyvsp[0].tok); }
#line 3943 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 176:
#line 995 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->vol->notify->mol_placement_failure = (byte) (yyvsp[0].tok); }
#line 3949 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 177:
#line 996 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->vol->notify->invalid_output_step_time = (byte) (yyvsp[0].tok); }
#line 3955 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 178:
#line 1000 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = WARN_COPE;  }
#line 3961 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 179:
#line 1001 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = WARN_WARN;  }
#line 3967 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 180:
#line 1002 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = WARN_ERROR; }
#line 3973 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 181:
#line 1008 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_checkpoint_infile(parse_state, (yyvsp[0].str))); }
#line 3979 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 182:
#line 1009 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_checkpoint_outfile(parse_state, (yyvsp[0].str))); }
#line 3985 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 183:
#line 1010 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_checkpoint_interval(parse_state, (yyvsp[-1].dbl), (yyvsp[0].tok))); }
#line 3991 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 184:
#line 1011 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_keep_checkpoint_files(parse_state, (yyvsp[0].tok))); }
#line 3997 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 185:
#line 1013 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_realtime_checkpoint(parse_state, (long) (yyvsp[-1].dbl), (yyvsp[0].tok))); }
#line 4003 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 186:
#line 1016 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = 0; }
#line 4009 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 187:
#line 1017 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = 1; }
#line 4015 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 188:
#line 1018 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = 0; }
#line 4021 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 189:
#line 1022 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { /* seconds */     (yyval.dbl) = (yyvsp[0].dbl); }
#line 4027 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 190:
#line 1023 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { /* mm:ss */       (yyval.dbl) = (yyvsp[-2].dbl) * 60 + (yyvsp[0].dbl); }
#line 4033 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 191:
#line 1024 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { /* hh:mm:ss */    (yyval.dbl) = (yyvsp[-4].dbl) * 3600 + (yyvsp[-2].dbl) * 60 + (yyvsp[0].dbl); }
#line 4039 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 192:
#line 1026 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { /* dd:hh:mm:ss */ (yyval.dbl) = (yyvsp[-6].dbl) * 86400 + (yyvsp[-4].dbl) * 3600 + (yyvsp[-2].dbl) * 60 + (yyvsp[0].dbl); }
#line 4045 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 193:
#line 1033 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_time_step(parse_state, (yyvsp[0].dbl))); }
#line 4051 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 194:
#line 1034 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_space_step(parse_state, (yyvsp[0].dbl))); }
#line 4057 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 195:
#line 1035 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_max_time_step(parse_state, (yyvsp[0].dbl))); }
#line 4063 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 196:
#line 1036 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_num_iterations(parse_state, (long long) (yyvsp[0].dbl))); }
#line 4069 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 197:
#line 1037 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->vol->randomize_smol_pos = !((yyvsp[0].tok)); }
#line 4075 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 198:
#line 1038 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->vol->use_expanded_list = (yyvsp[0].tok); }
#line 4081 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 199:
#line 1039 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->vol->vacancy_search_dist2 = max2d((yyvsp[0].dbl), 0.0); }
#line 4087 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 200:
#line 1040 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_num_radial_directions(parse_state, (int) (yyvsp[0].dbl))); }
#line 4093 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 201:
#line 1041 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->vol->fully_random = 1; }
#line 4099 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 202:
#line 1042 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_num_radial_subdivisions(parse_state, (int) (yyvsp[0].dbl))); }
#line 4105 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 203:
#line 1043 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_grid_density(parse_state, (yyvsp[0].dbl))); }
#line 4111 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 204:
#line 1044 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_interaction_radius(parse_state, (yyvsp[0].dbl))); }
#line 4117 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 205:
#line 1045 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->vol->surface_reversibility=(yyvsp[0].tok); parse_state->vol->volume_reversibility=(yyvsp[0].tok); }
#line 4123 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 206:
#line 1046 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->vol->surface_reversibility=1;  parse_state->vol->volume_reversibility=0;  }
#line 4129 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 207:
#line 1047 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->vol->surface_reversibility=0;  parse_state->vol->volume_reversibility=1;  }
#line 4135 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 208:
#line 1054 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->vol->mem_part_x = (int) (yyvsp[0].dbl); }
#line 4141 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 209:
#line 1055 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->vol->mem_part_y = (int) (yyvsp[0].dbl); }
#line 4147 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 210:
#line 1056 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->vol->mem_part_z = (int) (yyvsp[0].dbl); }
#line 4153 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 211:
#line 1057 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->vol->mem_part_pool = (int) (yyvsp[0].dbl); }
#line 4159 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 212:
#line 1061 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mcell_set_partition(parse_state->vol, (yyvsp[-2].tok), & (yyvsp[0].nlist))); }
#line 4165 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 213:
#line 1065 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = X_PARTS; }
#line 4171 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 214:
#line 1066 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = Y_PARTS; }
#line 4177 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 215:
#line 1067 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = Z_PARTS; }
#line 4183 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 218:
#line 1078 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { mdl_print_species_summary(parse_state->vol, (yyvsp[0].mcell_mol_spec)); }
#line 4189 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 219:
#line 1082 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { mdl_print_species_summaries(parse_state->vol, (yyvsp[-1].mcell_species_lst).species_head); }
#line 4195 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 220:
#line 1086 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.mcell_species_lst).species_count = 0; CHECK(mdl_add_to_species_list(&(yyval.mcell_species_lst), (yyvsp[0].mcell_mol_spec))); }
#line 4201 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 221:
#line 1087 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.mcell_species_lst) = (yyvsp[-1].mcell_species_lst); CHECK(mdl_add_to_species_list(&(yyval.mcell_species_lst), (yyvsp[0].mcell_mol_spec))); }
#line 4207 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 222:
#line 1096 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.mcell_mol_spec) = mdl_create_species(parse_state, (yyvsp[-6].str), (yyvsp[-4].diff_const).D, (yyvsp[-4].diff_const).is_2d, (yyvsp[-3].dbl), (yyvsp[-2].ival), (yyvsp[-1].dbl) )); }
#line 4213 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 224:
#line 1102 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.sym) = mdl_new_mol_species(parse_state, (yyvsp[0].str))); }
#line 4219 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 225:
#line 1106 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.diff_const).is_2d = 0; (yyval.diff_const).D = (yyvsp[0].dbl); CHECK(mdl_check_diffusion_constant(parse_state, & (yyval.diff_const).D)); }
#line 4225 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 226:
#line 1107 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.diff_const).is_2d = 1; (yyval.diff_const).D = (yyvsp[0].dbl); CHECK(mdl_check_diffusion_constant(parse_state, & (yyval.diff_const).D)); }
#line 4231 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 227:
#line 1111 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.dbl) = 0.0; }
#line 4237 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 228:
#line 1112 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                          if ((yyvsp[0].dbl) <= 0)
                                                          {
                                                            mdlerror_fmt(parse_state, "Requested custom time step of %.15g; custom time step must be positive.", (yyvsp[0].dbl));
                                                            return 1;
                                                          }

                                                          (yyval.dbl) = (yyvsp[0].dbl);
                                                      }
#line 4251 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 229:
#line 1121 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                          if ((yyvsp[0].dbl) <= 0)
                                                          {
                                                            mdlerror_fmt(parse_state, "Requested custom space step of %.15g; custom space step must be positive.", (yyvsp[0].dbl));
                                                            return 1;
                                                          }

                                                          (yyval.dbl) = -(yyvsp[0].dbl);
                                                      }
#line 4265 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 230:
#line 1132 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.ival) = 0; }
#line 4271 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 231:
#line 1133 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.ival) = 1; }
#line 4277 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 232:
#line 1137 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.dbl) = 0; }
#line 4283 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 233:
#line 1138 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                        if ((yyvsp[0].dbl) <= 0)
                                                        {
                                                          mdlerror_fmt(parse_state, "Requested maximum step length of %.15g; maximum step length must be positive.", (yyvsp[0].dbl));
                                                          return 1;
                                                        }
                                                        (yyval.dbl) = (yyvsp[0].dbl);
                                                      }
#line 4296 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 234:
#line 1148 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.sym) = mdl_existing_molecule(parse_state, (yyvsp[0].str))); }
#line 4302 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 235:
#line 1152 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.mol_type) = (yyvsp[0].mol_type); CHECKN((yyval.mol_type).mol_type = mdl_existing_surface_molecule(parse_state, (yyvsp[-1].str))); }
#line 4308 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 236:
#line 1156 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                        (yyval.mol_type) = (yyvsp[0].mol_type);
                                                        if (! (yyval.mol_type).orient_set)
                                                          (yyval.mol_type).orient = 0;
                                                        (yyval.mol_type).mol_type = (yyvsp[-1].sym);
                                                      }
#line 4319 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 243:
#line 1191 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { mdl_start_surface_class(parse_state, (yyvsp[-1].sym)); }
#line 4325 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 244:
#line 1193 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { mdl_finish_surface_class(parse_state); }
#line 4331 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 245:
#line 1196 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.sym) = mdl_existing_surface_class(parse_state, (yyvsp[0].str))); }
#line 4337 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 250:
#line 1213 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN(mdl_assemble_surface_reaction(parse_state, (yyvsp[-2].tok), parse_state->current_surface_class, (yyvsp[0].mol_type).mol_type, (yyvsp[0].mol_type).orient)); }
#line 4343 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 251:
#line 1216 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
              struct sym_entry *mol_sym = retrieve_sym("ALL_MOLECULES", parse_state->vol->mol_sym_table);
              if(!(yyvsp[0].mol_type).orient_set) (yyvsp[0].mol_type).orient = 0;
              CHECKN(mdl_assemble_surface_reaction(parse_state, (yyvsp[-3].tok), parse_state->current_surface_class, mol_sym, (yyvsp[0].mol_type).orient));}
#line 4352 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 252:
#line 1222 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN(mdl_assemble_concentration_clamp_reaction(parse_state, parse_state->current_surface_class, (yyvsp[-2].mol_type).mol_type, (yyvsp[-2].mol_type).orient, (yyvsp[0].dbl))); }
#line 4358 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 253:
#line 1225 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = RFLCT; }
#line 4364 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 254:
#line 1226 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = TRANSP; }
#line 4370 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 255:
#line 1227 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = SINK; }
#line 4376 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 258:
#line 1234 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->current_surface_class->sm_dat_head = (yyvsp[0].surf_mol_dat_list).sm_head; }
#line 4382 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 259:
#line 1241 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.surf_mol_dat_list) = (yyvsp[-1].surf_mol_dat_list); }
#line 4388 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 260:
#line 1245 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.surf_mol_dat_list) = (yyvsp[-1].surf_mol_dat_list); }
#line 4394 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 261:
#line 1249 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                          (yyvsp[0].surf_mol_dat)->quantity_type = SURFMOLDENS;
                                                          (yyval.surf_mol_dat_list).sm_tail = (yyval.surf_mol_dat_list).sm_head = (yyvsp[0].surf_mol_dat);
                                                      }
#line 4403 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 262:
#line 1254 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                          (yyval.surf_mol_dat_list) = (yyvsp[-1].surf_mol_dat_list);
                                                          (yyvsp[0].surf_mol_dat)->quantity_type = SURFMOLDENS;
                                                          (yyval.surf_mol_dat_list).sm_tail = (yyval.surf_mol_dat_list).sm_tail->next = (yyvsp[0].surf_mol_dat);
                                                      }
#line 4413 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 263:
#line 1262 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                          (yyvsp[0].surf_mol_dat)->quantity_type = SURFMOLNUM;
                                                          (yyval.surf_mol_dat_list).sm_tail = (yyval.surf_mol_dat_list).sm_head = (yyvsp[0].surf_mol_dat);
                                                      }
#line 4422 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 264:
#line 1267 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                          (yyval.surf_mol_dat_list) = (yyvsp[-1].surf_mol_dat_list);
                                                          (yyvsp[0].surf_mol_dat)->quantity_type = SURFMOLNUM;
                                                          (yyval.surf_mol_dat_list).sm_tail = (yyval.surf_mol_dat_list).sm_tail->next = (yyvsp[0].surf_mol_dat);
                                                      }
#line 4432 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 265:
#line 1275 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.surf_mol_dat) = mdl_new_surf_mol_data(parse_state, &(yyvsp[-2].mol_type), (yyvsp[0].dbl))); }
#line 4438 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 275:
#line 1304 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.react_arrow).catalyst = (yyvsp[-1].mol_type); (yyval.react_arrow).flags = ARROW_CATALYTIC; }
#line 4444 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 276:
#line 1309 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.react_arrow).catalyst = (yyvsp[-1].mol_type); (yyval.react_arrow).flags = ARROW_CATALYTIC | ARROW_BIDIRECTIONAL; }
#line 4450 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 277:
#line 1314 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.react_arrow).catalyst.mol_type = NULL; (yyval.react_arrow).flags = 0; }
#line 4456 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 279:
#line 1316 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.react_arrow).catalyst.mol_type = NULL; (yyval.react_arrow).flags = ARROW_BIDIRECTIONAL; }
#line 4462 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 281:
#line 1320 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.sym) = NULL; }
#line 4468 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 282:
#line 1321 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.sym) = mdl_new_rxn_pathname(parse_state, (yyvsp[0].str))); }
#line 4474 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 283:
#line 1327 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN(mdl_assemble_reaction(parse_state, (yyvsp[-5].mol_type_list).mol_type_head, &(yyvsp[-4].mol_type), &(yyvsp[-3].react_arrow), (yyvsp[-2].mol_type_list).mol_type_head, &(yyvsp[-1].react_rates), (yyvsp[0].sym))); }
#line 4480 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 284:
#line 1330 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_reaction_player_singleton(parse_state, & (yyval.mol_type_list), & (yyvsp[0].mol_type))); }
#line 4486 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 285:
#line 1331 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.mol_type_list) = (yyvsp[-2].mol_type_list); CHECK(mdl_add_reaction_player(parse_state, & (yyval.mol_type_list), & (yyvsp[0].mol_type))); }
#line 4492 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 287:
#line 1338 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.mol_type).mol_type = NULL; }
#line 4498 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 288:
#line 1339 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.mol_type) = (yyvsp[0].mol_type); }
#line 4504 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 289:
#line 1343 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.mol_type) = (yyvsp[0].mol_type); (yyval.mol_type).mol_type = (yyvsp[-1].sym); }
#line 4510 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 290:
#line 1346 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_reaction_player_singleton(parse_state, & (yyval.mol_type_list), & (yyvsp[0].mol_type))); }
#line 4516 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 291:
#line 1347 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.mol_type_list) = (yyvsp[-2].mol_type_list); CHECK(mdl_add_reaction_player(parse_state, & (yyval.mol_type_list), & (yyvsp[0].mol_type))); }
#line 4522 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 292:
#line 1350 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.mol_type).mol_type = NULL; (yyval.mol_type).orient_set = 0; }
#line 4528 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 296:
#line 1359 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                        if ((yyvsp[-1].react_rates).forward_rate.rate_type == RATE_UNSET)
                                                        {
                                                          mdlerror(parse_state, "invalid reaction rate specification: must specify a forward rate.");
                                                          return 1;
                                                        }

                                                        (yyval.react_rates) = (yyvsp[-1].react_rates);
                                                      }
#line 4542 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 297:
#line 1370 "../src/../src/mdlparse.y" /* yacc.c:1646  */
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
#line 4561 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 298:
#line 1387 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.react_rates).forward_rate = (yyvsp[0].react_rate); (yyval.react_rates).backward_rate.rate_type = RATE_UNSET; CHECK(mdl_valid_rate(parse_state, &(yyvsp[0].react_rate))); }
#line 4567 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 299:
#line 1388 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.react_rates).forward_rate = (yyvsp[0].react_rate); (yyval.react_rates).backward_rate.rate_type = RATE_UNSET; CHECK(mdl_valid_rate(parse_state, &(yyvsp[0].react_rate))); }
#line 4573 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 300:
#line 1389 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.react_rates).backward_rate = (yyvsp[0].react_rate); (yyval.react_rates).forward_rate.rate_type = RATE_UNSET; CHECK(mdl_valid_rate(parse_state, &(yyvsp[0].react_rate))); }
#line 4579 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 301:
#line 1393 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.react_rate).rate_type = RATE_CONSTANT; (yyval.react_rate).v.rate_constant = (yyvsp[0].dbl); }
#line 4585 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 302:
#line 1394 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.react_rate).rate_type = RATE_FILE; (yyval.react_rate).v.rate_file = (yyvsp[0].str); }
#line 4591 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 303:
#line 1395 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_reaction_rate_from_var(parse_state, & (yyval.react_rate), (yyvsp[0].sym))); }
#line 4597 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 304:
#line 1406 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_release_pattern(parse_state, (yyvsp[-3].sym), &(yyvsp[-1].rpat))); }
#line 4603 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 305:
#line 1409 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.sym) = mdl_new_release_pattern(parse_state, (yyvsp[0].str))); }
#line 4609 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 306:
#line 1412 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.sym) = mdl_existing_release_pattern_or_rxn_pathname(parse_state, (yyvsp[0].str))); }
#line 4615 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 307:
#line 1416 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                        (yyval.rpat).delay = 0;
                                                        (yyval.rpat).release_interval = FOREVER;
                                                        (yyval.rpat).train_interval = FOREVER;
                                                        (yyval.rpat).train_duration = FOREVER;
                                                        (yyval.rpat).number_of_trains = 1;
                                                      }
#line 4627 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 308:
#line 1424 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.rpat) = (yyvsp[-3].rpat); (yyval.rpat).delay = (yyvsp[0].dbl) / parse_state->vol->time_unit; }
#line 4633 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 309:
#line 1426 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.rpat) = (yyvsp[-3].rpat); (yyval.rpat).release_interval = (yyvsp[0].dbl) / parse_state->vol->time_unit; }
#line 4639 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 310:
#line 1428 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.rpat) = (yyvsp[-3].rpat); (yyval.rpat).train_interval = (yyvsp[0].dbl) / parse_state->vol->time_unit; }
#line 4645 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 311:
#line 1430 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.rpat) = (yyvsp[-3].rpat); (yyval.rpat).train_duration = (yyvsp[0].dbl) / parse_state->vol->time_unit; }
#line 4651 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 312:
#line 1432 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.rpat) = (yyvsp[-3].rpat); (yyval.rpat).number_of_trains = (yyvsp[0].ival); }
#line 4657 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 313:
#line 1435 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.ival) = (int) (yyvsp[0].dbl); }
#line 4663 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 314:
#line 1436 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.ival) = INT_MAX; }
#line 4669 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 315:
#line 1443 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->current_object = parse_state->vol->root_instance; }
#line 4675 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 316:
#line 1444 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                        check_regions(parse_state->vol->root_instance, (yyvsp[0].obj));
                                                        add_child_objects(parse_state->vol->root_instance, (yyvsp[0].obj), (yyvsp[0].obj));
                                                        parse_state->current_object = parse_state->vol->root_object;
                                                      }
#line 4685 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 317:
#line 1454 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { add_child_objects(parse_state->vol->root_object, (yyvsp[0].obj), (yyvsp[0].obj)); }
#line 4691 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 323:
#line 1470 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.sym) = mdl_start_object(parse_state, (yyvsp[0].str))); }
#line 4697 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 325:
#line 1476 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { mdl_finish_object(parse_state); }
#line 4703 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 329:
#line 1489 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { transform_translate(parse_state->vol, parse_state->current_object->t_matrix, (yyvsp[0].vec3)); }
#line 4709 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 330:
#line 1490 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { transform_scale(parse_state->current_object->t_matrix, (yyvsp[0].vec3)); }
#line 4715 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 331:
#line 1491 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_transform_rotate(parse_state, parse_state->current_object->t_matrix, (yyvsp[-2].vec3), (yyvsp[0].dbl))); }
#line 4721 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 332:
#line 1500 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                          struct object *the_object = (struct object *) (yyvsp[-5].sym)->value;
                                                          the_object->object_type = META_OBJ;
                                                          add_child_objects(the_object, (yyvsp[-2].obj_list).obj_head, (yyvsp[-2].obj_list).obj_tail);
                                                          (yyval.obj) = the_object;
                                                      }
#line 4732 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 333:
#line 1509 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { mdl_object_list_singleton(& (yyval.obj_list), (yyvsp[0].obj)); }
#line 4738 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 334:
#line 1510 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.obj_list) = (yyvsp[-1].obj_list); mdl_add_object_to_list(& (yyval.obj_list), (yyvsp[0].obj)); }
#line 4744 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 337:
#line 1519 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_deep_copy_object(parse_state, (struct object *) (yyvsp[-3].sym)->value, (struct object *) (yyvsp[-1].sym)->value)); }
#line 4750 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 338:
#line 1521 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.obj) = (struct object *) (yyvsp[-6].sym)->value; }
#line 4756 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 341:
#line 1531 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_start_release_site(parse_state, (yyvsp[-2].sym), SHAPE_UNDEFINED)); }
#line 4762 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 342:
#line 1535 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.obj) = mdl_finish_release_site(parse_state, (yyvsp[-7].sym))); }
#line 4768 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 343:
#line 1538 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_release_site_geometry_region(parse_state, parse_state->current_release_site, parse_state->current_object, (yyvsp[0].rev))); }
#line 4774 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 344:
#line 1539 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_release_site_geometry_object(parse_state, parse_state->current_release_site, (struct object *) (yyvsp[0].sym)->value)); }
#line 4780 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 345:
#line 1540 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->current_release_site->release_shape = SHAPE_SPHERICAL; }
#line 4786 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 346:
#line 1541 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->current_release_site->release_shape = SHAPE_CUBIC; }
#line 4792 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 347:
#line 1542 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->current_release_site->release_shape = SHAPE_ELLIPTIC; }
#line 4798 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 348:
#line 1543 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->current_release_site->release_shape = SHAPE_RECTANGULAR; }
#line 4804 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 349:
#line 1544 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->current_release_site->release_shape = SHAPE_SPHERICAL_SHELL; }
#line 4810 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 350:
#line 1545 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                          parse_state->current_release_site->release_shape = SHAPE_LIST;
                                                          parse_state->current_release_site->release_number_method = CONSTNUM;
                                                      }
#line 4819 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 351:
#line 1552 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.rev) = new_release_region_expr_term((yyvsp[0].sym))); }
#line 4825 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 352:
#line 1553 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.rev) = (yyvsp[-1].rev); }
#line 4831 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 353:
#line 1554 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.rev) = new_release_region_expr_binary((yyvsp[-2].rev), (yyvsp[0].rev), REXP_UNION)); }
#line 4837 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 354:
#line 1555 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.rev) = new_release_region_expr_binary((yyvsp[-2].rev), (yyvsp[0].rev), REXP_SUBTRACTION)); }
#line 4843 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 355:
#line 1556 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.rev) = new_release_region_expr_binary((yyvsp[-2].rev), (yyvsp[0].rev), REXP_INTERSECTION)); }
#line 4849 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 356:
#line 1561 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_start_release_site(parse_state, (yyvsp[-2].sym), (yyvsp[-1].tok))); }
#line 4855 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 357:
#line 1564 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.obj) = mdl_finish_release_site(parse_state, (yyvsp[-6].sym))); }
#line 4861 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 358:
#line 1567 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = SHAPE_SPHERICAL; }
#line 4867 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 359:
#line 1568 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = SHAPE_CUBIC; }
#line 4873 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 360:
#line 1569 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = SHAPE_ELLIPTIC; }
#line 4879 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 361:
#line 1570 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = SHAPE_RECTANGULAR; }
#line 4885 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 362:
#line 1571 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = SHAPE_SPHERICAL_SHELL; }
#line 4891 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 365:
#line 1579 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.sym) = mdl_existing_num_or_array(parse_state, (yyvsp[0].str))); }
#line 4897 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 366:
#line 1583 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { set_release_site_location(parse_state->vol, parse_state->current_release_site, (yyvsp[0].vec3)); }
#line 4903 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 367:
#line 1584 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_release_site_molecule(parse_state, parse_state->current_release_site, & (yyvsp[0].mol_type))); }
#line 4909 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 368:
#line 1585 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                        if (parse_state->current_release_site->release_shape == SHAPE_LIST)
                                                        {
                                                          mdlerror(parse_state, "molecules are already specified in a list--cannot set number or density.");
                                                          return 1;
                                                        }
                                                      }
#line 4921 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 369:
#line 1592 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_release_site_diameter(parse_state, parse_state->current_release_site, (yyvsp[0].dbl) * (((yyvsp[-2].tok) == SITE_RADIUS) ? 2.0 : 1.0))); }
#line 4927 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 370:
#line 1593 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_release_site_diameter_array(parse_state, parse_state->current_release_site, (yyvsp[0].nlist).value_count, (yyvsp[0].nlist).value_head, ((yyvsp[-2].tok) == SITE_RADIUS) ? 2.0 : 1.0)); }
#line 4933 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 371:
#line 1594 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_release_site_diameter_var(parse_state, parse_state->current_release_site, ((yyvsp[-2].tok) == SITE_RADIUS) ? 2.0 : 1.0, (yyvsp[0].sym))); }
#line 4939 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 372:
#line 1595 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_release_site_probability(parse_state, parse_state->current_release_site, (yyvsp[0].dbl))); }
#line 4945 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 373:
#line 1597 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_release_site_pattern(parse_state, parse_state->current_release_site, (yyvsp[0].sym))); }
#line 4951 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 374:
#line 1599 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_release_site_molecule_positions(parse_state, parse_state->current_release_site, & (yyvsp[-1].rsm_list))); }
#line 4957 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 375:
#line 1603 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = SITE_DIAMETER; }
#line 4963 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 376:
#line 1604 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = SITE_RADIUS; }
#line 4969 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 381:
#line 1616 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { set_release_site_constant_number(parse_state->current_release_site, (yyvsp[0].dbl)); }
#line 4975 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 382:
#line 1619 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { set_release_site_constant_number(parse_state->current_release_site, (yyvsp[-1].dbl)); }
#line 4981 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 383:
#line 1626 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { set_release_site_gaussian_number(parse_state->current_release_site, (yyvsp[-4].dbl), (yyvsp[-1].dbl)); }
#line 4987 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 384:
#line 1634 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { set_release_site_volume_dependent_number(parse_state->current_release_site, (yyvsp[-7].dbl), (yyvsp[-4].dbl), (yyvsp[-1].dbl)); }
#line 4993 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 385:
#line 1638 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_release_site_concentration(parse_state, parse_state->current_release_site, (yyvsp[0].dbl))); }
#line 4999 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 386:
#line 1639 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(set_release_site_density(parse_state->current_release_site, (yyvsp[0].dbl))); }
#line 5005 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 387:
#line 1643 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { release_single_molecule_singleton(& (yyval.rsm_list), (yyvsp[0].rsm)); }
#line 5011 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 388:
#line 1645 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.rsm_list) = (yyvsp[-1].rsm_list); add_release_single_molecule_to_list(& (yyval.rsm_list), (yyvsp[0].rsm)); }
#line 5017 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 389:
#line 1649 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.rsm) = mdl_new_release_single_molecule(parse_state, &(yyvsp[-1].mol_type), (yyvsp[0].vec3))); }
#line 5023 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 391:
#line 1660 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                        CHECKN((yyval.obj) = mdl_new_polygon_list(
                                                          parse_state, (yyvsp[-4].str), (yyvsp[-1].vertlist).vertex_count, (yyvsp[-1].vertlist).vertex_head,
                                                          (yyvsp[0].ecl).connection_count, (yyvsp[0].ecl).connection_head));
                                                      }
#line 5033 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 392:
#line 1669 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                          (yyval.obj) = (struct object *) (yyvsp[-3].obj);
                                                          CHECK(mdl_finish_polygon_list(parse_state, (yyval.obj)));
                                                      }
#line 5042 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 393:
#line 1675 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.vertlist) = (yyvsp[-1].vertlist); }
#line 5048 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 394:
#line 1678 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.vertlistitem) = mdl_new_vertex_list_item((yyvsp[0].vec3))); }
#line 5054 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 395:
#line 1681 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { mdl_vertex_list_singleton(& (yyval.vertlist), (yyvsp[0].vertlistitem)); }
#line 5060 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 396:
#line 1682 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.vertlist) = (yyvsp[-1].vertlist); mdl_add_vertex_to_list(& (yyval.vertlist), (yyvsp[0].vertlistitem)); }
#line 5066 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 397:
#line 1687 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.ecl) = (yyvsp[-1].ecl); }
#line 5072 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 398:
#line 1691 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { mdl_element_connection_list_singleton(& (yyval.ecl), (yyvsp[0].elem_conn)); }
#line 5078 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 399:
#line 1693 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.ecl) = (yyvsp[-1].ecl); mdl_add_element_connection_to_list(& (yyval.ecl), (yyvsp[0].elem_conn)); }
#line 5084 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 400:
#line 1696 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.elem_conn) = mdl_new_element_connection(parse_state, & (yyvsp[0].nlist))); }
#line 5090 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 405:
#line 1712 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN(parse_state->current_region = mdl_get_region(parse_state, parse_state->current_object, "REMOVED")); }
#line 5096 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 406:
#line 1714 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                          parse_state->current_region->element_list_head = (yyvsp[-1].elem_list).elml_head;
                                                          if (parse_state->current_object->object_type == POLY_OBJ)
                                                          {
                                                            CHECK(mdl_normalize_elements(parse_state, parse_state->current_region,0));
                                                          }
                                                      }
#line 5108 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 409:
#line 1728 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = Z_POS; }
#line 5114 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 410:
#line 1729 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = Z_NEG; }
#line 5120 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 411:
#line 1730 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = Y_NEG; }
#line 5126 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 412:
#line 1731 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = Y_POS; }
#line 5132 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 413:
#line 1732 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = X_NEG; }
#line 5138 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 414:
#line 1733 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = X_POS; }
#line 5144 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 415:
#line 1734 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = ALL_SIDES; }
#line 5150 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 417:
#line 1740 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.elem_list) = (yyvsp[-1].elem_list); mdl_add_elements_to_list(& (yyval.elem_list), (yyvsp[0].elem_list).elml_head, (yyvsp[0].elem_list).elml_tail); }
#line 5156 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 420:
#line 1746 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.elem_list).elml_tail = (yyval.elem_list).elml_head = (yyvsp[0].elem_list_item); }
#line 5162 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 421:
#line 1747 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.elem_list).elml_tail = (yyval.elem_list).elml_head = (yyvsp[0].elem_list_item); }
#line 5168 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 422:
#line 1752 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.elem_list) = (yyvsp[-1].elem_list); }
#line 5174 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 423:
#line 1757 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.elem_list) = (yyvsp[-1].elem_list); mdl_set_elements_to_exclude((yyval.elem_list).elml_head); }
#line 5180 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 425:
#line 1764 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.elem_list).elml_tail = (yyval.elem_list).elml_head = (yyvsp[0].elem_list_item); }
#line 5186 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 426:
#line 1765 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.elem_list) = (yyvsp[-2].elem_list); mdl_add_elements_to_list(& (yyval.elem_list), (yyvsp[0].elem_list_item), (yyvsp[0].elem_list_item)); }
#line 5192 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 427:
#line 1768 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.elem_list_item) = new_element_list((unsigned int) (yyvsp[0].dbl), (unsigned int) (yyvsp[0].dbl))); }
#line 5198 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 428:
#line 1769 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.elem_list_item) = new_element_list((unsigned int) (yyvsp[-2].dbl), (unsigned int) (yyvsp[0].dbl))); }
#line 5204 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 429:
#line 1770 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.elem_list_item) = mdl_new_element_side(parse_state, (yyvsp[0].tok))); }
#line 5210 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 430:
#line 1773 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.elem_list_item) = mdl_new_element_previous_region(parse_state, parse_state->current_object, parse_state->current_region, (yyvsp[0].str), (yyvsp[-2].tok))); }
#line 5216 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 431:
#line 1776 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = 0; }
#line 5222 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 432:
#line 1777 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = 1; }
#line 5228 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 433:
#line 1780 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.elem_list_item) = mdl_new_element_patch(parse_state, parse_state->current_polygon, (yyvsp[-2].vec3), (yyvsp[0].vec3), (yyvsp[-4].tok))); }
#line 5234 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 434:
#line 1783 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = 0; }
#line 5240 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 435:
#line 1784 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = 1; }
#line 5246 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 439:
#line 1800 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->current_region = (yyvsp[-1].reg); }
#line 5252 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 440:
#line 1801 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_region_elements(parse_state, (yyvsp[-3].reg), (yyvsp[0].elem_list).elml_head, (yyvsp[-3].reg)->parent->object_type == POLY_OBJ)); }
#line 5258 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 441:
#line 1803 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->current_region = NULL; }
#line 5264 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 442:
#line 1811 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                        CHECKN(mdl_new_voxel_list(parse_state, (yyvsp[-4].sym),
                                                                                  (yyvsp[-1].vertlist).vertex_count, (yyvsp[-1].vertlist).vertex_head,
                                                                                  (yyvsp[0].ecl).connection_count, (yyvsp[0].ecl).connection_head));
                                                      }
#line 5274 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 443:
#line 1817 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.obj) = (struct object *) (yyvsp[-7].sym)->value; }
#line 5280 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 444:
#line 1822 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.ecl) = (yyvsp[-1].ecl); }
#line 5286 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 445:
#line 1825 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.elem_conn) = mdl_new_tet_element_connection(parse_state, & (yyvsp[0].nlist))); }
#line 5292 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 446:
#line 1829 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                          (yyval.ecl).connection_head = (yyval.ecl).connection_tail = (yyvsp[0].elem_conn);
                                                          (yyval.ecl).connection_count = 1;
                                                      }
#line 5301 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 447:
#line 1833 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                          (yyval.ecl) = (yyvsp[-1].ecl);
                                                          (yyval.ecl).connection_tail = (yyval.ecl).connection_tail->next = (yyvsp[0].elem_conn);
                                                          ++ (yyval.ecl).connection_count;
                                                      }
#line 5311 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 448:
#line 1844 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN(mdl_new_box_object(parse_state, (yyvsp[-8].sym), (yyvsp[-3].vec3), (yyvsp[-1].vec3))); }
#line 5317 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 449:
#line 1845 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_triangulate_box_object(parse_state, (yyvsp[-10].sym), parse_state->current_polygon, (yyvsp[-2].dbl))); }
#line 5323 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 450:
#line 1847 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                          CHECK(mdl_finish_box_object(parse_state, (yyvsp[-13].sym)));
                                                          (yyval.obj) = (struct object *) (yyvsp[-13].sym)->value;
                                                      }
#line 5332 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 451:
#line 1853 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.dbl) = 0.0; }
#line 5338 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 452:
#line 1854 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                        (yyval.dbl) = (yyvsp[0].dbl);
                                                        if ((yyval.dbl) < 2.0)
                                                        {
                                                          mdlerror(parse_state, "invalid aspect ratio requested (must be greater than or equal to 2.0)");
                                                          return 1;
                                                        }
                                                      }
#line 5351 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 456:
#line 1880 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_start_existing_obj_region_def(parse_state, (yyvsp[0].sym))); }
#line 5357 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 457:
#line 1881 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->current_region = (yyvsp[-1].reg); }
#line 5363 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 458:
#line 1883 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { mdl_set_region_elements(parse_state, (yyvsp[-4].reg), (yyvsp[0].elem_list).elml_head, 1); }
#line 5369 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 459:
#line 1885 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                          parse_state->current_region = NULL;
                                                          parse_state->current_polygon = NULL;
                                                          parse_state->current_object = parse_state->vol->root_object;
                                                      }
#line 5379 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 460:
#line 1892 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.reg) = mdl_create_region(parse_state, parse_state->current_object, (yyvsp[0].str))); }
#line 5385 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 464:
#line 1903 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { mdl_add_surf_mol_to_region(parse_state->current_region, & (yyvsp[0].surf_mol_dat_list)); }
#line 5391 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 465:
#line 1907 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { mdl_set_region_surface_class(parse_state, parse_state->current_region, (yyvsp[0].sym)); }
#line 5397 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 469:
#line 1926 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->current_region = (struct region *) (yyvsp[-1].sym)->value; }
#line 5403 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 470:
#line 1928 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->current_region = NULL; }
#line 5409 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 471:
#line 1936 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                          parse_state->header_comment = NULL;  /* No header by default */
                                                          parse_state->exact_time_flag = 1;    /* Print exact_time column in TRIGGER output by default */
                                                      }
#line 5418 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 472:
#line 1942 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_add_reaction_output_block_to_world(parse_state, (int) (yyvsp[-4].dbl), & (yyvsp[-2].ro_otimes), & (yyvsp[-1].ro_sets))); }
#line 5424 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 473:
#line 1946 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.dbl) = COUNTBUFFERSIZE; }
#line 5430 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 474:
#line 1947 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                          double temp_value = (yyvsp[0].dbl);
                                                          if (!(temp_value >= 1.0 && temp_value < UINT_MAX))
                                                          {
                                                            mdlerror_fmt(parse_state, "Requested buffer size of %.15g lines is invalid.  Suggested range is 100-1000000.", temp_value);
                                                            return 1;
                                                          }
                                                          (yyval.dbl) = (yyvsp[0].dbl);
                                                      }
#line 5444 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 478:
#line 1963 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.ro_otimes).type = OUTPUT_BY_STEP; (yyval.ro_otimes).step = (yyvsp[0].dbl); }
#line 5450 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 479:
#line 1967 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                        (yyval.ro_otimes).type = OUTPUT_BY_ITERATION_LIST;
                                                        (yyval.ro_otimes).values = (yyvsp[0].nlist);
                                                      }
#line 5459 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 480:
#line 1975 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                        (yyval.ro_otimes).type = OUTPUT_BY_TIME_LIST;
                                                        (yyval.ro_otimes).values = (yyvsp[0].nlist);
                                                      }
#line 5468 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 481:
#line 1982 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.ro_sets).set_head = (yyval.ro_sets).set_tail = (yyvsp[0].ro_set); }
#line 5474 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 482:
#line 1984 "../src/../src/mdlparse.y" /* yacc.c:1646  */
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
#line 5489 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 484:
#line 1998 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.ro_set) = NULL; }
#line 5495 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 485:
#line 1999 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.ro_set) = NULL; }
#line 5501 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 486:
#line 2003 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {  parse_state->count_flags = 0; }
#line 5507 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 487:
#line 2005 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.ro_set) = mdl_populate_output_set(parse_state, parse_state->header_comment, parse_state->exact_time_flag, (yyvsp[-3].ro_cols).column_head, (yyvsp[-1].tok), (yyvsp[0].str))); }
#line 5513 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 488:
#line 2009 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.str) = NULL; }
#line 5519 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 489:
#line 2010 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.str) = ((yyvsp[0].tok) ? "" : NULL); }
#line 5525 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 490:
#line 2011 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.str) = (yyvsp[0].str); }
#line 5531 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 491:
#line 2015 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->header_comment = (yyvsp[0].str); }
#line 5537 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 492:
#line 2019 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->exact_time_flag = (yyvsp[0].tok); }
#line 5543 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 494:
#line 2025 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                          (yyval.ro_cols) = (yyvsp[-2].ro_cols);
                                                          (yyval.ro_cols).column_tail->next = (yyvsp[0].ro_cols).column_head;
                                                          (yyval.ro_cols).column_tail = (yyvsp[0].ro_cols).column_tail;
                                                      }
#line 5553 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 495:
#line 2033 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_single_count_expr(parse_state, & (yyval.ro_cols), (yyvsp[-1].cnt), (yyvsp[0].str))); }
#line 5559 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 496:
#line 2037 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.cnt) = mdl_new_oexpr_constant(parse_state, (yyvsp[0].dbl))); }
#line 5565 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 498:
#line 2039 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.cnt) = mdl_join_oexpr_tree(parse_state, (yyvsp[-1].cnt), NULL, '(')); }
#line 5571 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 499:
#line 2040 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.cnt) = mdl_join_oexpr_tree(parse_state, (yyvsp[-2].cnt),   (yyvsp[0].cnt), '+')); }
#line 5577 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 500:
#line 2041 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.cnt) = mdl_join_oexpr_tree(parse_state, (yyvsp[-2].cnt),   (yyvsp[0].cnt), '-')); }
#line 5583 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 501:
#line 2042 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.cnt) = mdl_join_oexpr_tree(parse_state, (yyvsp[-2].cnt),   (yyvsp[0].cnt), '*')); }
#line 5589 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 502:
#line 2043 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.cnt) = mdl_join_oexpr_tree(parse_state, (yyvsp[-2].cnt),   (yyvsp[0].cnt), '/')); }
#line 5595 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 503:
#line 2044 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.cnt) = mdl_join_oexpr_tree(parse_state, (yyvsp[0].cnt), NULL, '_')); }
#line 5601 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 504:
#line 2045 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.cnt) = mdl_sum_oexpr((yyvsp[-1].cnt))); }
#line 5607 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 505:
#line 2050 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->count_flags |= COUNT_PRESENT; }
#line 5613 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 506:
#line 2051 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.cnt) = (yyvsp[-1].cnt); }
#line 5619 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 507:
#line 2052 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.cnt) = mdl_new_oexpr_constant(parse_state, (yyvsp[-1].dbl))); }
#line 5625 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 508:
#line 2053 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->count_flags |= TRIGGER_PRESENT; }
#line 5631 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 509:
#line 2054 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.cnt) = (yyvsp[-1].cnt); }
#line 5637 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 510:
#line 2057 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = FILE_OVERWRITE; }
#line 5643 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 511:
#line 2058 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = FILE_SUBSTITUTE; }
#line 5649 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 512:
#line 2059 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = FILE_APPEND; }
#line 5655 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 513:
#line 2060 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = FILE_APPEND_HEADER; }
#line 5661 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 514:
#line 2061 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = FILE_CREATE; }
#line 5667 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 516:
#line 2067 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.sym) = mdl_existing_rxn_pathname_or_molecule(parse_state, (yyvsp[0].str))); }
#line 5673 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 517:
#line 2071 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                        (yyval.mol_type) = (yyvsp[0].mol_type);
                                                        if ((yyval.mol_type).orient > 0)
                                                          (yyval.mol_type).orient = 1;
                                                        else if ((yyval.mol_type).orient < 0)
                                                          (yyval.mol_type).orient = -1;
                                                        CHECKN((yyval.mol_type).mol_type = mdl_existing_molecule(parse_state, (yyvsp[-1].str)));
                                                      }
#line 5686 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 521:
#line 2088 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.cnt) = mdl_count_syntax_1(parse_state, (yyvsp[-3].sym), (yyvsp[-1].sym), (yyvsp[0].tok), parse_state->count_flags)); }
#line 5692 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 522:
#line 2093 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.cnt) = mdl_count_syntax_2(parse_state, (yyvsp[-3].mol_type).mol_type, (yyvsp[-3].mol_type).orient, (yyvsp[-1].sym), (yyvsp[0].tok), parse_state->count_flags)); }
#line 5698 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 523:
#line 2098 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.cnt) = mdl_count_syntax_3(parse_state, (yyvsp[-3].str), (yyvsp[-1].sym), (yyvsp[0].tok), parse_state->count_flags)); }
#line 5704 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 524:
#line 2101 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.sym) = NULL; }
#line 5710 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 525:
#line 2102 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.sym) = (yyvsp[0].sym); }
#line 5716 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 526:
#line 2103 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.sym) = (yyvsp[0].sym); }
#line 5722 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 527:
#line 2106 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = REPORT_NOTHING; }
#line 5728 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 528:
#line 2107 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = (yyvsp[0].tok); }
#line 5734 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 529:
#line 2110 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = REPORT_FRONT_HITS; }
#line 5740 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 530:
#line 2111 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = REPORT_BACK_HITS; }
#line 5746 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 531:
#line 2112 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = REPORT_ALL_HITS; }
#line 5752 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 532:
#line 2113 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = REPORT_FRONT_CROSSINGS; }
#line 5758 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 533:
#line 2114 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = REPORT_BACK_CROSSINGS; }
#line 5764 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 534:
#line 2115 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = REPORT_ALL_CROSSINGS; }
#line 5770 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 535:
#line 2116 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = REPORT_CONCENTRATION; }
#line 5776 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 536:
#line 2117 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = REPORT_ENCLOSED; }
#line 5782 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 537:
#line 2120 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.str) = NULL; }
#line 5788 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 538:
#line 2121 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.str) = (yyvsp[0].str); }
#line 5794 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 539:
#line 2128 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_new_viz_output_block(parse_state)); }
#line 5800 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 540:
#line 2131 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { }
#line 5806 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 543:
#line 2140 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_viz_mode(parse_state->vol->viz_blocks, CELLBLENDER_MODE)); }
#line 5812 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 544:
#line 2141 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_viz_mode(parse_state->vol->viz_blocks, (yyvsp[0].ival))); }
#line 5818 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 545:
#line 2144 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.ival) = NO_VIZ_MODE; }
#line 5824 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 546:
#line 2145 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.ival) = ASCII_MODE; }
#line 5830 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 547:
#line 2146 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.ival) = CELLBLENDER_MODE; }
#line 5836 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 549:
#line 2151 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                        if ((yyvsp[0].frame_list).frame_head)
                                                        {
                                                          (yyvsp[0].frame_list).frame_tail->next = parse_state->vol->viz_blocks->frame_data_head;
                                                          parse_state->vol->viz_blocks->frame_data_head = (yyvsp[0].frame_list).frame_head;
                                                        }
                                                      }
#line 5848 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 551:
#line 2164 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_viz_filename_prefix(parse_state, parse_state->vol->viz_blocks, (yyvsp[0].str))); }
#line 5854 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 552:
#line 2170 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.frame_list) = (yyvsp[-1].frame_list); }
#line 5860 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 554:
#line 2176 "../src/../src/mdlparse.y" /* yacc.c:1646  */
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
#line 5876 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 555:
#line 2190 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.frame_list).frame_head = (yyval.frame_list).frame_tail = NULL; }
#line 5882 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 559:
#line 2202 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_viz_state(parse_state, & (yyval.ival), (yyvsp[0].dbl))); }
#line 5888 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 560:
#line 2203 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.ival) = INCLUDE_OBJ; }
#line 5894 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 563:
#line 2213 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_viz_include_molecules(parse_state, parse_state->vol->viz_blocks, (yyvsp[-1].symlist), (yyvsp[0].ival))); }
#line 5900 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 564:
#line 2214 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_viz_include_all_molecules(parse_state->vol->viz_blocks, (yyvsp[0].ival))); }
#line 5906 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 565:
#line 2218 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.symlist) = mdl_existing_molecule_list(parse_state, (yyvsp[0].str))); }
#line 5912 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 566:
#line 2219 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.symlist) = mdl_existing_molecules_wildcard(parse_state, (yyvsp[0].str))); }
#line 5918 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 567:
#line 2223 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_new_viz_all_times(parse_state, & (yyval.nlist))); }
#line 5924 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 569:
#line 2229 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.frame_list) = (yyvsp[-1].frame_list); }
#line 5930 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 571:
#line 2235 "../src/../src/mdlparse.y" /* yacc.c:1646  */
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
#line 5948 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 572:
#line 2252 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_new_viz_mol_frames(parse_state, parse_state->vol->viz_blocks, & (yyval.frame_list), OUTPUT_BY_TIME_LIST, (yyvsp[-2].tok), & (yyvsp[0].nlist))); }
#line 5954 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 573:
#line 2256 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_new_viz_all_iterations(parse_state, & (yyval.nlist))); }
#line 5960 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 575:
#line 2263 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.frame_list) = (yyvsp[-1].frame_list); }
#line 5966 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 577:
#line 2269 "../src/../src/mdlparse.y" /* yacc.c:1646  */
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
#line 5984 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 578:
#line 2286 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_new_viz_mol_frames(parse_state, parse_state->vol->viz_blocks, & (yyval.frame_list), OUTPUT_BY_ITERATION_LIST, (yyvsp[-2].tok), & (yyvsp[0].nlist))); }
#line 5990 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 579:
#line 2289 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = ALL_MOL_DATA; }
#line 5996 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 580:
#line 2290 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = MOL_POS; }
#line 6002 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 581:
#line 2291 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = MOL_ORIENT; }
#line 6008 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 582:
#line 2305 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                          struct volume_output_item *vo;
                                                          CHECKN(vo = mdl_new_volume_output_item(parse_state, (yyvsp[-6].str), & (yyvsp[-5].species_lst), (yyvsp[-4].vec3), (yyvsp[-3].vec3), (yyvsp[-2].vec3), (yyvsp[-1].otimes)));
                                                          vo->next = parse_state->vol->volume_output_head;
                                                          parse_state->vol->volume_output_head = vo;
                                                      }
#line 6019 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 583:
#line 2314 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.str) = (yyvsp[0].str); }
#line 6025 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 585:
#line 2320 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                          (yyval.species_lst) = (yyvsp[-1].species_lst);
                                                          (yyval.species_lst).species_count += (yyvsp[0].species_lst).species_count;
                                                          (yyval.species_lst).species_tail->next = (yyvsp[0].species_lst).species_head;
                                                          (yyval.species_lst).species_tail = (yyvsp[0].species_lst).species_tail;
                                                      }
#line 6036 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 586:
#line 2329 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.species_lst) = (yyvsp[0].species_lst); }
#line 6042 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 587:
#line 2332 "../src/../src/mdlparse.y" /* yacc.c:1646  */
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
#line 6062 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 588:
#line 2350 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.species_lst).species_tail = (yyval.species_lst).species_head = (yyvsp[0].species_lst_item); (yyval.species_lst).species_count = 1; }
#line 6068 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 589:
#line 2352 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                        (yyval.species_lst) = (yyvsp[-2].species_lst);
                                                        (yyval.species_lst).species_tail = (yyval.species_lst).species_tail->next = (yyvsp[0].species_lst_item);
                                                        ++ (yyval.species_lst).species_count;
                                                      }
#line 6078 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 590:
#line 2360 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.vec3) = (yyvsp[0].vec3); }
#line 6084 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 591:
#line 2364 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.vec3) = (yyvsp[0].vec3); }
#line 6090 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 592:
#line 2368 "../src/../src/mdlparse.y" /* yacc.c:1646  */
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
#line 6113 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 593:
#line 2389 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.otimes) = mdl_new_output_times_default(parse_state)); }
#line 6119 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 594:
#line 2390 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.otimes) = mdl_new_output_times_step(parse_state, (yyvsp[0].dbl))); }
#line 6125 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 595:
#line 2391 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.otimes) = mdl_new_output_times_iterations(parse_state, & (yyvsp[0].nlist))); }
#line 6131 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 596:
#line 2392 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.otimes) = mdl_new_output_times_time(parse_state, & (yyvsp[0].nlist))); }
#line 6137 "mdlparse.c" /* yacc.c:1646  */
    break;


#line 6141 "mdlparse.c" /* yacc.c:1646  */
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
#line 2395 "../src/../src/mdlparse.y" /* yacc.c:1906  */






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

  if ((mpv.header_comment != 0) || (mpv.header_comment != 0)) {
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
