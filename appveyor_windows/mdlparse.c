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
    EFFECTOR_GRID_DENSITY = 314,
    ELEMENT_CONNECTIONS = 315,
    ELLIPTIC = 316,
    ELLIPTIC_RELEASE_SITE = 317,
    EQUAL = 318,
    ERROR = 319,
    ESTIMATE_CONCENTRATION = 320,
    EXCLUDE_ELEMENTS = 321,
    EXCLUDE_PATCH = 322,
    EXCLUDE_REGION = 323,
    EXIT = 324,
    EXP = 325,
    EXPRESSION = 326,
    FALSE = 327,
    FCLOSE = 328,
    FILENAME = 329,
    FILENAME_PREFIX = 330,
    FILE_OUTPUT_REPORT = 331,
    FINAL_SUMMARY = 332,
    FLOOR = 333,
    FOPEN = 334,
    FORMAT = 335,
    FPRINTF = 336,
    FPRINT_TIME = 337,
    FRONT = 338,
    FRONT_CROSSINGS = 339,
    FRONT_HITS = 340,
    GAUSSIAN_RELEASE_NUMBER = 341,
    HEADER = 342,
    HIGH_PROBABILITY_THRESHOLD = 343,
    HIGH_REACTION_PROBABILITY = 344,
    IGNORED = 345,
    INCLUDE_ELEMENTS = 346,
    INCLUDE_FILE = 347,
    INCLUDE_PATCH = 348,
    INCLUDE_REGION = 349,
    INPUT_FILE = 350,
    INSTANTIATE = 351,
    LLINTEGER = 352,
    FULLY_RANDOM = 353,
    INTERACTION_RADIUS = 354,
    ITERATION_LIST = 355,
    ITERATION_NUMBERS = 356,
    ITERATION_REPORT = 357,
    ITERATIONS = 358,
    KEEP_CHECKPOINT_FILES = 359,
    LEFT = 360,
    LIFETIME_THRESHOLD = 361,
    LIFETIME_TOO_SHORT = 362,
    LIST = 363,
    LOCATION = 364,
    LOG = 365,
    LOG10 = 366,
    MAX_TOK = 367,
    MAXIMUM_STEP_LENGTH = 368,
    MEAN_DIAMETER = 369,
    MEAN_NUMBER = 370,
    MEMORY_PARTITION_X = 371,
    MEMORY_PARTITION_Y = 372,
    MEMORY_PARTITION_Z = 373,
    MEMORY_PARTITION_POOL = 374,
    MICROSCOPIC_REVERSIBILITY = 375,
    MIN_TOK = 376,
    MISSED_REACTIONS = 377,
    MISSED_REACTION_THRESHOLD = 378,
    MISSING_SURFACE_ORIENTATION = 379,
    MOD = 380,
    MODE = 381,
    MODIFY_SURFACE_REGIONS = 382,
    MOLECULE = 383,
    MOLECULE_COLLISION_REPORT = 384,
    MOLECULE_DENSITY = 385,
    MOLECULE_NUMBER = 386,
    MOLECULE_POSITIONS = 387,
    MOLECULES = 388,
    MOLECULE_PLACEMENT_FAILURE = 389,
    NAME_LIST = 390,
    NEGATIVE_DIFFUSION_CONSTANT = 391,
    NEGATIVE_REACTION_RATE = 392,
    NO = 393,
    NOEXIT = 394,
    NONE = 395,
    NO_SPECIES = 396,
    NOT_EQUAL = 397,
    NOTIFICATIONS = 398,
    NUMBER_OF_SUBUNITS = 399,
    NUMBER_OF_TRAINS = 400,
    NUMBER_TO_RELEASE = 401,
    OBJECT = 402,
    OFF = 403,
    ON = 404,
    ORIENTATIONS = 405,
    OUTPUT_BUFFER_SIZE = 406,
    INVALID_OUTPUT_STEP_TIME = 407,
    OVERWRITTEN_OUTPUT_FILE = 408,
    PARTITION_LOCATION_REPORT = 409,
    PARTITION_X = 410,
    PARTITION_Y = 411,
    PARTITION_Z = 412,
    PERIODIC_BOX = 413,
    PERIODIC_X = 414,
    PERIODIC_Y = 415,
    PERIODIC_Z = 416,
    PERIODIC_TRADITIONAL = 417,
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
    REACTION_DATA_OUTPUT = 430,
    REACTION_OUTPUT_REPORT = 431,
    REAL = 432,
    RECTANGULAR_RELEASE_SITE = 433,
    RECTANGULAR_TOKEN = 434,
    REFLECTIVE = 435,
    RELEASE_EVENT_REPORT = 436,
    RELEASE_INTERVAL = 437,
    RELEASE_PATTERN = 438,
    RELEASE_PROBABILITY = 439,
    RELEASE_SITE = 440,
    REMOVE_ELEMENTS = 441,
    RIGHT = 442,
    ROTATE = 443,
    ROUND_OFF = 444,
    SCALE = 445,
    SEED = 446,
    SHAPE = 447,
    SHOW_EXACT_TIME = 448,
    SIN = 449,
    SITE_DIAMETER = 450,
    SITE_RADIUS = 451,
    SPACE_STEP = 452,
    SPHERICAL = 453,
    SPHERICAL_RELEASE_SITE = 454,
    SPHERICAL_SHELL = 455,
    SPHERICAL_SHELL_SITE = 456,
    SPRINTF = 457,
    SQRT = 458,
    STANDARD_DEVIATION = 459,
    PERIODIC_BOX_INITIAL = 460,
    STEP = 461,
    STRING_TO_NUM = 462,
    STR_VALUE = 463,
    SUBUNIT = 464,
    SUBUNIT_RELATIONSHIPS = 465,
    SUMMATION_OPERATOR = 466,
    SURFACE_CLASS = 467,
    SURFACE_ONLY = 468,
    TAN = 469,
    TARGET_ONLY = 470,
    TET_ELEMENT_CONNECTIONS = 471,
    THROUGHPUT_REPORT = 472,
    TIME_LIST = 473,
    TIME_POINTS = 474,
    TIME_STEP = 475,
    TIME_STEP_MAX = 476,
    TO = 477,
    TOP = 478,
    TRAIN_DURATION = 479,
    TRAIN_INTERVAL = 480,
    TRANSLATE = 481,
    TRANSPARENT = 482,
    TRIGGER = 483,
    TRUE = 484,
    UNLIMITED = 485,
    USELESS_VOLUME_ORIENTATION = 486,
    VACANCY_SEARCH_DISTANCE = 487,
    VAR = 488,
    VARYING_PROBABILITY_REPORT = 489,
    VERTEX_LIST = 490,
    VIZ_OUTPUT = 491,
    VIZ_OUTPUT_REPORT = 492,
    VIZ_VALUE = 493,
    VOLUME_DATA_OUTPUT = 494,
    VOLUME_OUTPUT_REPORT = 495,
    VOLUME_DEPENDENT_RELEASE_NUMBER = 496,
    VOLUME_ONLY = 497,
    VOXEL_COUNT = 498,
    VOXEL_LIST = 499,
    VOXEL_SIZE = 500,
    WARNING = 501,
    WARNINGS = 502,
    WORLD = 503,
    YES = 504,
    UNARYMINUS = 505
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
#define EFFECTOR_GRID_DENSITY 314
#define ELEMENT_CONNECTIONS 315
#define ELLIPTIC 316
#define ELLIPTIC_RELEASE_SITE 317
#define EQUAL 318
#define ERROR 319
#define ESTIMATE_CONCENTRATION 320
#define EXCLUDE_ELEMENTS 321
#define EXCLUDE_PATCH 322
#define EXCLUDE_REGION 323
#define EXIT 324
#define EXP 325
#define EXPRESSION 326
#define FALSE 327
#define FCLOSE 328
#define FILENAME 329
#define FILENAME_PREFIX 330
#define FILE_OUTPUT_REPORT 331
#define FINAL_SUMMARY 332
#define FLOOR 333
#define FOPEN 334
#define FORMAT 335
#define FPRINTF 336
#define FPRINT_TIME 337
#define FRONT 338
#define FRONT_CROSSINGS 339
#define FRONT_HITS 340
#define GAUSSIAN_RELEASE_NUMBER 341
#define HEADER 342
#define HIGH_PROBABILITY_THRESHOLD 343
#define HIGH_REACTION_PROBABILITY 344
#define IGNORED 345
#define INCLUDE_ELEMENTS 346
#define INCLUDE_FILE 347
#define INCLUDE_PATCH 348
#define INCLUDE_REGION 349
#define INPUT_FILE 350
#define INSTANTIATE 351
#define LLINTEGER 352
#define FULLY_RANDOM 353
#define INTERACTION_RADIUS 354
#define ITERATION_LIST 355
#define ITERATION_NUMBERS 356
#define ITERATION_REPORT 357
#define ITERATIONS 358
#define KEEP_CHECKPOINT_FILES 359
#define LEFT 360
#define LIFETIME_THRESHOLD 361
#define LIFETIME_TOO_SHORT 362
#define LIST 363
#define LOCATION 364
#define LOG 365
#define LOG10 366
#define MAX_TOK 367
#define MAXIMUM_STEP_LENGTH 368
#define MEAN_DIAMETER 369
#define MEAN_NUMBER 370
#define MEMORY_PARTITION_X 371
#define MEMORY_PARTITION_Y 372
#define MEMORY_PARTITION_Z 373
#define MEMORY_PARTITION_POOL 374
#define MICROSCOPIC_REVERSIBILITY 375
#define MIN_TOK 376
#define MISSED_REACTIONS 377
#define MISSED_REACTION_THRESHOLD 378
#define MISSING_SURFACE_ORIENTATION 379
#define MOD 380
#define MODE 381
#define MODIFY_SURFACE_REGIONS 382
#define MOLECULE 383
#define MOLECULE_COLLISION_REPORT 384
#define MOLECULE_DENSITY 385
#define MOLECULE_NUMBER 386
#define MOLECULE_POSITIONS 387
#define MOLECULES 388
#define MOLECULE_PLACEMENT_FAILURE 389
#define NAME_LIST 390
#define NEGATIVE_DIFFUSION_CONSTANT 391
#define NEGATIVE_REACTION_RATE 392
#define NO 393
#define NOEXIT 394
#define NONE 395
#define NO_SPECIES 396
#define NOT_EQUAL 397
#define NOTIFICATIONS 398
#define NUMBER_OF_SUBUNITS 399
#define NUMBER_OF_TRAINS 400
#define NUMBER_TO_RELEASE 401
#define OBJECT 402
#define OFF 403
#define ON 404
#define ORIENTATIONS 405
#define OUTPUT_BUFFER_SIZE 406
#define INVALID_OUTPUT_STEP_TIME 407
#define OVERWRITTEN_OUTPUT_FILE 408
#define PARTITION_LOCATION_REPORT 409
#define PARTITION_X 410
#define PARTITION_Y 411
#define PARTITION_Z 412
#define PERIODIC_BOX 413
#define PERIODIC_X 414
#define PERIODIC_Y 415
#define PERIODIC_Z 416
#define PERIODIC_TRADITIONAL 417
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
#define REACTION_DATA_OUTPUT 430
#define REACTION_OUTPUT_REPORT 431
#define REAL 432
#define RECTANGULAR_RELEASE_SITE 433
#define RECTANGULAR_TOKEN 434
#define REFLECTIVE 435
#define RELEASE_EVENT_REPORT 436
#define RELEASE_INTERVAL 437
#define RELEASE_PATTERN 438
#define RELEASE_PROBABILITY 439
#define RELEASE_SITE 440
#define REMOVE_ELEMENTS 441
#define RIGHT 442
#define ROTATE 443
#define ROUND_OFF 444
#define SCALE 445
#define SEED 446
#define SHAPE 447
#define SHOW_EXACT_TIME 448
#define SIN 449
#define SITE_DIAMETER 450
#define SITE_RADIUS 451
#define SPACE_STEP 452
#define SPHERICAL 453
#define SPHERICAL_RELEASE_SITE 454
#define SPHERICAL_SHELL 455
#define SPHERICAL_SHELL_SITE 456
#define SPRINTF 457
#define SQRT 458
#define STANDARD_DEVIATION 459
#define PERIODIC_BOX_INITIAL 460
#define STEP 461
#define STRING_TO_NUM 462
#define STR_VALUE 463
#define SUBUNIT 464
#define SUBUNIT_RELATIONSHIPS 465
#define SUMMATION_OPERATOR 466
#define SURFACE_CLASS 467
#define SURFACE_ONLY 468
#define TAN 469
#define TARGET_ONLY 470
#define TET_ELEMENT_CONNECTIONS 471
#define THROUGHPUT_REPORT 472
#define TIME_LIST 473
#define TIME_POINTS 474
#define TIME_STEP 475
#define TIME_STEP_MAX 476
#define TO 477
#define TOP 478
#define TRAIN_DURATION 479
#define TRAIN_INTERVAL 480
#define TRANSLATE 481
#define TRANSPARENT 482
#define TRIGGER 483
#define TRUE 484
#define UNLIMITED 485
#define USELESS_VOLUME_ORIENTATION 486
#define VACANCY_SEARCH_DISTANCE 487
#define VAR 488
#define VARYING_PROBABILITY_REPORT 489
#define VERTEX_LIST 490
#define VIZ_OUTPUT 491
#define VIZ_OUTPUT_REPORT 492
#define VIZ_VALUE 493
#define VOLUME_DATA_OUTPUT 494
#define VOLUME_OUTPUT_REPORT 495
#define VOLUME_DEPENDENT_RELEASE_NUMBER 496
#define VOLUME_ONLY 497
#define VOXEL_COUNT 498
#define VOXEL_LIST 499
#define VOXEL_SIZE 500
#define WARNING 501
#define WARNINGS 502
#define WORLD 503
#define YES 504
#define UNARYMINUS 505

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


#line 744 "mdlparse.c" /* yacc.c:355  */
};
# define YYSTYPE_IS_TRIVIAL 1
# define YYSTYPE_IS_DECLARED 1
#endif



int mdlparse (struct mdlparse_vars *parse_state, yyscan_t scanner);

#endif /* !YY_MDL_MDLPARSE_H_INCLUDED  */

/* Copy the second part of user declarations.  */

#line 758 "mdlparse.c" /* yacc.c:358  */

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
#define YYFINAL  146
/* YYLAST -- Last index in YYTABLE.  */
#define YYLAST   2485

/* YYNTOKENS -- Number of terminals.  */
#define YYNTOKENS  271
/* YYNNTS -- Number of nonterminals.  */
#define YYNNTS  289
/* YYNRULES -- Number of rules.  */
#define YYNRULES  615
/* YYNSTATES -- Number of states.  */
#define YYNSTATES  1214

/* YYTRANSLATE[YYX] -- Symbol number corresponding to YYX as returned
   by yylex, with out-of-bounds checking.  */
#define YYUNDEFTOK  2
#define YYMAXUTOK   505

#define YYTRANSLATE(YYX)                                                \
  ((unsigned int) (YYX) <= YYMAXUTOK ? yytranslate[YYX] : YYUNDEFTOK)

/* YYTRANSLATE[TOKEN-NUM] -- Symbol number corresponding to TOKEN-NUM
   as returned by yylex, without out-of-bounds checking.  */
static const yytype_uint16 yytranslate[] =
{
       0,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,   251,   262,
     266,   267,   255,   253,   263,   254,     2,   256,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,   252,   261,
     269,   250,   268,     2,   270,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,   259,     2,   260,   258,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,   264,     2,   265,     2,     2,     2,     2,
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
     245,   246,   247,   248,   249,   257
};

#if YYDEBUG
  /* YYRLINE[YYN] -- Source line where rule number YYN was defined.  */
static const yytype_uint16 yyrline[] =
{
       0,   593,   593,   597,   598,   603,   604,   605,   606,   607,
     608,   609,   610,   611,   612,   613,   614,   615,   616,   617,
     618,   619,   620,   621,   622,   623,   628,   631,   634,   637,
     640,   643,   646,   647,   650,   651,   652,   653,   654,   655,
     658,   659,   660,   661,   665,   666,   667,   677,   689,   692,
     695,   707,   708,   721,   722,   728,   751,   752,   753,   754,
     757,   760,   763,   764,   775,   778,   781,   782,   785,   786,
     789,   790,   793,   794,   797,   801,   802,   803,   804,   805,
     806,   807,   808,   809,   810,   811,   812,   813,   814,   815,
     816,   817,   818,   819,   820,   821,   822,   823,   824,   825,
     826,   827,   828,   829,   830,   834,   835,   839,   840,   841,
     842,   845,   851,   852,   853,   854,   855,   856,   857,   860,
     864,   867,   870,   873,   876,   879,   880,   890,   891,   892,
     904,   908,   914,   919,   923,   932,   936,   937,   941,   942,
     943,   944,   945,   946,   947,   948,   949,   950,   951,   952,
     953,   954,   955,   956,   957,   961,   962,   966,   970,   971,
     978,   982,   983,   987,   988,   989,   990,   991,   992,   993,
     994,   995,   996,   997,   998,   999,  1000,  1001,  1002,  1006,
    1007,  1008,  1014,  1015,  1016,  1017,  1018,  1022,  1023,  1024,
    1028,  1029,  1030,  1031,  1039,  1040,  1041,  1042,  1043,  1044,
    1045,  1046,  1047,  1048,  1049,  1050,  1051,  1052,  1053,  1060,
    1061,  1062,  1063,  1067,  1071,  1072,  1073,  1080,  1081,  1084,
    1088,  1092,  1093,  1097,  1105,  1108,  1112,  1113,  1117,  1118,
    1127,  1138,  1139,  1143,  1144,  1154,  1158,  1162,  1175,  1176,
    1181,  1185,  1191,  1192,  1197,  1197,  1202,  1205,  1207,  1212,
    1213,  1217,  1220,  1226,  1231,  1232,  1233,  1236,  1237,  1240,
    1244,  1248,  1255,  1259,  1268,  1272,  1281,  1288,  1293,  1294,
    1297,  1300,  1301,  1304,  1305,  1306,  1309,  1314,  1320,  1321,
    1322,  1323,  1326,  1327,  1331,  1336,  1337,  1340,  1344,  1345,
    1349,  1352,  1353,  1356,  1357,  1361,  1362,  1365,  1376,  1393,
    1394,  1395,  1399,  1400,  1401,  1408,  1415,  1418,  1422,  1429,
    1431,  1433,  1435,  1437,  1441,  1442,  1449,  1449,  1460,  1463,
    1464,  1465,  1466,  1467,  1476,  1479,  1482,  1485,  1487,  1491,
    1495,  1496,  1497,  1502,  1515,  1516,  1519,  1520,  1525,  1524,
    1531,  1532,  1537,  1536,  1544,  1545,  1546,  1547,  1548,  1549,
    1550,  1551,  1558,  1559,  1560,  1561,  1562,  1567,  1566,  1573,
    1574,  1575,  1576,  1577,  1581,  1582,  1585,  1589,  1590,  1591,
    1598,  1599,  1600,  1601,  1602,  1603,  1605,  1610,  1611,  1615,
    1616,  1617,  1618,  1623,  1624,  1630,  1637,  1645,  1646,  1650,
    1651,  1656,  1659,  1667,  1664,  1682,  1685,  1688,  1689,  1693,
    1698,  1699,  1703,  1706,  1708,  1714,  1715,  1719,  1719,  1731,
    1732,  1735,  1736,  1737,  1738,  1739,  1740,  1741,  1745,  1746,
    1751,  1752,  1753,  1754,  1758,  1763,  1767,  1771,  1772,  1775,
    1776,  1777,  1780,  1783,  1784,  1787,  1790,  1791,  1795,  1801,
    1802,  1807,  1808,  1807,  1818,  1815,  1828,  1832,  1836,  1840,
    1850,  1853,  1847,  1860,  1861,  1857,  1870,  1871,  1875,  1876,
    1880,  1881,  1885,  1886,  1889,  1890,  1904,  1910,  1911,  1916,
    1917,  1919,  1916,  1928,  1931,  1933,  1938,  1939,  1943,  1950,
    1956,  1957,  1962,  1962,  1972,  1971,  1982,  1983,  1994,  1995,
    1996,  1999,  2003,  2011,  2018,  2019,  2033,  2034,  2035,  2039,
    2039,  2045,  2046,  2047,  2051,  2055,  2059,  2060,  2069,  2073,
    2074,  2075,  2076,  2077,  2078,  2079,  2080,  2081,  2086,  2086,
    2088,  2089,  2089,  2093,  2094,  2095,  2096,  2097,  2100,  2103,
    2107,  2117,  2118,  2119,  2120,  2121,  2122,  2126,  2131,  2136,
    2141,  2146,  2151,  2155,  2156,  2157,  2160,  2161,  2164,  2165,
    2166,  2167,  2168,  2169,  2170,  2171,  2174,  2175,  2182,  2182,
    2189,  2190,  2194,  2195,  2198,  2199,  2200,  2204,  2205,  2215,
    2218,  2222,  2228,  2229,  2244,  2245,  2246,  2250,  2256,  2257,
    2261,  2262,  2266,  2268,  2272,  2273,  2277,  2278,  2281,  2287,
    2288,  2305,  2310,  2311,  2315,  2321,  2322,  2339,  2343,  2344,
    2345,  2352,  2368,  2372,  2373,  2383,  2386,  2404,  2405,  2414,
    2418,  2422,  2443,  2444,  2445,  2446
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
  "DIFFUSION_CONSTANT_3D", "DIFFUSION_CONSTANT_REPORT",
  "EFFECTOR_GRID_DENSITY", "ELEMENT_CONNECTIONS", "ELLIPTIC",
  "ELLIPTIC_RELEASE_SITE", "EQUAL", "ERROR", "ESTIMATE_CONCENTRATION",
  "EXCLUDE_ELEMENTS", "EXCLUDE_PATCH", "EXCLUDE_REGION", "EXIT", "EXP",
  "EXPRESSION", "FALSE", "FCLOSE", "FILENAME", "FILENAME_PREFIX",
  "FILE_OUTPUT_REPORT", "FINAL_SUMMARY", "FLOOR", "FOPEN", "FORMAT",
  "FPRINTF", "FPRINT_TIME", "FRONT", "FRONT_CROSSINGS", "FRONT_HITS",
  "GAUSSIAN_RELEASE_NUMBER", "HEADER", "HIGH_PROBABILITY_THRESHOLD",
  "HIGH_REACTION_PROBABILITY", "IGNORED", "INCLUDE_ELEMENTS",
  "INCLUDE_FILE", "INCLUDE_PATCH", "INCLUDE_REGION", "INPUT_FILE",
  "INSTANTIATE", "LLINTEGER", "FULLY_RANDOM", "INTERACTION_RADIUS",
  "ITERATION_LIST", "ITERATION_NUMBERS", "ITERATION_REPORT", "ITERATIONS",
  "KEEP_CHECKPOINT_FILES", "LEFT", "LIFETIME_THRESHOLD",
  "LIFETIME_TOO_SHORT", "LIST", "LOCATION", "LOG", "LOG10", "MAX_TOK",
  "MAXIMUM_STEP_LENGTH", "MEAN_DIAMETER", "MEAN_NUMBER",
  "MEMORY_PARTITION_X", "MEMORY_PARTITION_Y", "MEMORY_PARTITION_Z",
  "MEMORY_PARTITION_POOL", "MICROSCOPIC_REVERSIBILITY", "MIN_TOK",
  "MISSED_REACTIONS", "MISSED_REACTION_THRESHOLD",
  "MISSING_SURFACE_ORIENTATION", "MOD", "MODE", "MODIFY_SURFACE_REGIONS",
  "MOLECULE", "MOLECULE_COLLISION_REPORT", "MOLECULE_DENSITY",
  "MOLECULE_NUMBER", "MOLECULE_POSITIONS", "MOLECULES",
  "MOLECULE_PLACEMENT_FAILURE", "NAME_LIST", "NEGATIVE_DIFFUSION_CONSTANT",
  "NEGATIVE_REACTION_RATE", "NO", "NOEXIT", "NONE", "NO_SPECIES",
  "NOT_EQUAL", "NOTIFICATIONS", "NUMBER_OF_SUBUNITS", "NUMBER_OF_TRAINS",
  "NUMBER_TO_RELEASE", "OBJECT", "OFF", "ON", "ORIENTATIONS",
  "OUTPUT_BUFFER_SIZE", "INVALID_OUTPUT_STEP_TIME",
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
      61,    38,    58,    43,    45,    42,    47,   505,    94,    91,
      93,    59,    39,    44,   123,   125,    40,    41,    62,    60,
      64
};
# endif

#define YYPACT_NINF -1041

#define yypact_value_is_default(Yystate) \
  (!!((Yystate) == (-1041)))

#define YYTABLE_NINF -393

#define yytable_value_is_error(Yytable_value) \
  0

  /* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
     STATE-NUM.  */
static const yytype_int16 yypact[] =
{
    2238,  -157,  -104,   -56,   -31,    -9,     4,   -44,  -179,  -153,
     -44,   -44,     1,    68,    96,   108,   116,   158,   181, -1041,
     191,   202,   232,   256,   260,   266,   270,   283,   304,   323,
   -1041, -1041, -1041,   343,   195,   277,   364,   369,   353,   374,
     360,   386,   390,   394, -1041,   384,   385,   387,   652,  2238,
   -1041,   -52, -1041, -1041,   405, -1041, -1041,   578, -1041, -1041,
   -1041, -1041, -1041, -1041, -1041, -1041, -1041, -1041, -1041, -1041,
     409, -1041, -1041, -1041, -1041, -1041, -1041, -1041, -1041, -1041,
   -1041, -1041,    28, -1041, -1041, -1041, -1041,   498, -1041, -1041,
   -1041, -1041, -1041, -1041, -1041, -1041, -1041,   279,   279,   -43,
    2060,   -43,  2060, -1041, -1041,   399,   -44,   -44, -1041,   401,
   -1041,   412, -1041,   -44,   -44,  2060,   -44,   -44,   -44,   -43,
     -44,  2060,  2060,   279,  2060,  2060,  2060,  2060,   425,   -44,
     751, -1041,   639,   -43,   -43,  1762,  2060,   528,  2060,   -44,
    2060,  2060,  2060, -1041,   605,  1054, -1041, -1041,  1498,   420,
      -3,   343, -1041, -1041,   343, -1041,   343, -1041, -1041,   343,
     343,   343, -1041, -1041, -1041, -1041, -1041, -1041, -1041, -1041,
     422, -1041, -1041, -1041, -1041, -1041,   432, -1041, -1041,   424,
     427,   434,   441,   442,   444,   447,   449, -1041,   450,   454,
     458,   465,   471, -1041, -1041, -1041, -1041,   472, -1041,   475,
     478,   482,   484,  2060,  2060,  2060, -1041,   111, -1041, -1041,
   -1041, -1041, -1041,  1168,   -14,    66,  -177, -1041, -1041,   333,
   -1041,  -160, -1041, -1041,    63, -1041, -1041, -1041,   -66, -1041,
   -1041, -1041,   -16, -1041,  1231, -1041,   476,   428,   462,   432,
   -1041,   610, -1041,  1231,  1231, -1041,  1231,  1231,  1231,  1231,
   -1041, -1041, -1041,   493,   494,   135, -1041,   511,   520,   523,
     524,   525,   529,   531,   532,   537,   539,   541,   542,   543,
     544,   547,   548,   549,   552,   477, -1041,   562,   432, -1041,
     551, -1041,  1231,  1231,   564, -1041,  1231, -1041,   515,  1231,
    1231,  1231,   693,   570,   688,   585,   586,   590,   595,   596,
     597,   602,   608,   612,   623,   628,   629,   634,   636,   638,
     650,   433, -1041,  1768,   983, -1041, -1041,  1231,  1237, -1041,
    1298,   432,   600,   -43, -1041, -1041, -1041, -1041,   821,   -44,
   -1041,   667, -1041,   667,   -43,   -43,  2060,  2060,  2060,  2060,
    2060,  2060,  2060,  2060,  2060,  2060,  2060,  2060,  2060,  2060,
    2060,  2060,   -43,  2060,   654,   654,   576, -1041, -1041,  2060,
    2060,  2060,  2060,  2060, -1041,  2060, -1041,   663,   666,   249,
   -1041, -1041, -1041, -1041, -1041,  2060, -1041,   263, -1041, -1041,
   -1041, -1041, -1041,   -44,   -44,    -8,    38, -1041, -1041, -1041,
     658, -1041, -1041, -1041,   -43,   -43,   -44, -1041, -1041, -1041,
     279,   279,   279,    99,   279,   279,  1424,   279,   279,   279,
    2060,   279,    99,   279,   279,   279,    99,    99, -1041, -1041,
      -3,   175, -1041,  2060,   -18,   -43,   672,    -5, -1041,   -43,
     676,   -39, -1041,   -42,   -42,   -42,  2060,   -42,  2060,   -42,
     -42,  2060,   -42,   -42,   -42,   -42,   -42,   -42,   -42, -1041,
   -1041,  2060,   244, -1041,  1231,   662,   684,   774, -1041,   162,
     -44, -1041, -1041,   748,   677,   727,   918,   884, -1041, -1041,
     601,   616,   637,   643,   675,   699,   709,   725,   742,   757,
     550,   832,  1043,  1095,   775,   850,   -33,   874, -1041,   217,
     217,   654,   654, -1041,  1213,  2060,  2060,   695,   697,   734,
     653, -1041, -1041, -1041, -1041,   333, -1041, -1041,   704,   183,
   -1041,  -165, -1041, -1041, -1041,   -41,   701,   719,   720,   721,
     724, -1041,    42,   -44, -1041,   692,   715, -1041, -1041, -1041,
   -1041, -1041, -1041, -1041, -1041, -1041, -1041,  1231, -1041, -1041,
   -1041, -1041,  1231, -1041, -1041, -1041, -1041, -1041, -1041, -1041,
     698, -1041,  1648, -1041,  1231,   736,   737,   739,   -57, -1041,
   -1041, -1041, -1041,    69,   743,   713,   -50, -1041, -1041, -1041,
   -1041,   432,   -44,   752, -1041,   758, -1041, -1041, -1041, -1041,
   -1041, -1041,  1231, -1041,  1231, -1041, -1041,  1231, -1041, -1041,
   -1041, -1041, -1041, -1041, -1041,   291, -1041,  1768,   -43,    -3,
    -123,   182, -1041,   755,   918,    -3,   744, -1041,   767,   769,
     759,   772,   776,   761,   782,   784,   787, -1041, -1041,   793,
     788,   918, -1041,   801, -1041, -1041, -1041, -1041, -1041,   790,
   -1041,   196, -1041, -1041, -1041, -1041, -1041, -1041, -1041, -1041,
   -1041, -1041,  2060,  2060,  2060,  2060, -1041, -1041, -1041, -1041,
    2060,  1231,  1231,  2060,  2060, -1041,   871, -1041, -1041,   753,
   -1041, -1041,   704, -1041,   704, -1041, -1041,   145, -1041,  2060,
    1888,  2060,  2060,  2060, -1041,   -44,   791,   792, -1041, -1041,
   -1041, -1041, -1041,  -118, -1041, -1041, -1041,   800,   216, -1041,
   -1041,   -72,    -3, -1041, -1041,   600, -1041,    -3,  2060,    -3,
     813,   815, -1041,   -62, -1041, -1041, -1041, -1041,   218, -1041,
   -1041, -1041,   -43,    14, -1041, -1041, -1041, -1041,   814,    -3,
     816,   825,  2060, -1041,   432,   805,   810,   343,   828,   829,
     830, -1041, -1041, -1041, -1041,   -21,   918, -1041, -1041,    94,
      -3, -1041,  2060,  2060,   960,    -3,   -44,   -44,  2060,   -44,
    2060,    -3,   968,   182, -1041,  1986,    -3, -1041, -1041,   914,
     941,   966,  1068,  1225,  1231,  1231,   839,   827,    51, -1041,
   -1041,   -41,  1041,   841, -1041, -1041,  1231, -1041,  1231, -1041,
    1231,  1231,  1231,   844,   -44,   -44, -1041, -1041,    10, -1041,
   -1041,   846, -1041, -1041, -1041, -1041,   936, -1041,  1231, -1041,
    1044,   279,   347, -1041, -1041, -1041,   432,   836,   845,   848,
     -51, -1041, -1041, -1041, -1041,   -44, -1041,  1986,   860,   125,
     240, -1041,    -3, -1041,    -3,  1986,    -3, -1041, -1041, -1041,
   -1041, -1041, -1041,  -171,   493, -1041,   375,   182, -1041, -1041,
   -1041, -1041,   127,   182,  1231,  1231,   865, -1041, -1041,    -3,
     140, -1041,  1231, -1041, -1041,  1231, -1041,   868, -1041,  1326,
   -1041, -1041, -1041, -1041,   212, -1041,    24, -1041, -1041, -1041,
   -1041,  2060,  2060, -1041, -1041,  1648,  1648, -1041, -1041,   600,
     267, -1041,   -44, -1041,  2060,   333,   870,   148, -1041,   156,
   -1041,   333, -1041,   858,   -44,   875, -1041, -1041, -1041,   432,
   -1041, -1041, -1041,   867,   869, -1041,   347,   347, -1041,    93,
   -1041,   902, -1041,    36, -1041,    36, -1041, -1041, -1041,  1326,
   -1041, -1041, -1041,  1986,   883,   887,   890,   879,  2060,  1115,
   -1041,   882, -1041, -1041,   311,  -171,  -171,  -171, -1041, -1041,
   -1041, -1041,  2060, -1041, -1041, -1041,  2060, -1041, -1041,   886,
     899,   197, -1041, -1041, -1041,  1231,  1231, -1041, -1041, -1041,
    1041, -1041,  1231, -1041,  2060, -1041, -1041, -1041, -1041, -1041,
     575, -1041,   279,   988,   906,  2060,   347,   912, -1041,   417,
     347,   200,   -43,   347,   347,   347,   347, -1041, -1041, -1041,
   -1041,     9, -1041,   903,     3,    26, -1041,   904, -1041,    -3,
    2060,    -3, -1041,  1133,   925, -1041,   182,  2060, -1041,   924,
     924, -1041,   367,   211,   -44, -1041, -1041,   920,  1231,   933,
   -1041, -1041,   948, -1041, -1041,   575, -1041, -1041, -1041, -1041,
     959, -1041,   961, -1041,   963,  1050,    34,  1147,   513,    34,
   -1041, -1041,   949,   955,   958,   -43,   432,   295,   295, -1041,
   -1041, -1041, -1041,     7,   977, -1041, -1041, -1041, -1041,   977,
   -1041, -1041,     5, -1041,  1231, -1041, -1041,  2060, -1041, -1041,
    1231,   979, -1041,   995,   189, -1041,   922,  1222, -1041,   957,
     992, -1041, -1041,   -44,    -3,   279,   997,  1103,  1005,  1006,
    1008,  1009,  1015, -1041, -1041, -1041, -1041, -1041, -1041, -1041,
   -1041,  1016, -1041, -1041,  1001, -1041, -1041, -1041, -1041, -1041,
    2060, -1041, -1041, -1041, -1041, -1041,  1231,    24,  2060,  2060,
   -1041, -1041, -1041, -1041, -1041, -1041, -1041, -1041, -1041, -1041,
     474,  1013, -1041,   575, -1041,  1020, -1041,  1336,  1336,   -34,
   -1041,  1023, -1041,   279,  1037, -1041,   269, -1041,   269,   269,
   -1041, -1041, -1041,  1231, -1041,   947,    83,   575,  2060, -1041,
    1336,   274,   278, -1041,    -3, -1041,   279,  1038, -1041,   493,
   -1041,  1039,  1045,  1048,   182, -1041,  1062,   575,  1231, -1041,
   -1041, -1041, -1041, -1041, -1041,   174, -1041,   174, -1041,   174,
   -1041, -1041,  2060, -1041, -1041, -1041, -1041, -1041, -1041, -1041,
   -1041, -1041,  1051, -1041,  1051,  1051,  1026,    32,   711, -1041,
   -1041, -1041, -1041, -1041
};

  /* YYDEFACT[STATE-NUM] -- Default reduction number in state STATE-NUM.
     Performed when YYTABLE does not specify something else to do.  Zero
     means the default is an error.  */
static const yytype_uint16 yydefact[] =
{
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,   316,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     214,   215,   216,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,    27,     0,     0,     0,     0,     2,
       3,   324,     5,     6,     0,     7,   112,     0,   113,   114,
     115,   116,   117,   118,     8,     9,    10,    11,    14,    12,
       0,    15,   217,   218,    16,   238,   239,    17,    18,    20,
      19,   318,     0,   319,   320,   341,   340,     0,   322,   323,
      13,   321,    21,    22,    23,    24,    25,     0,     0,     0,
       0,     0,     0,   224,   219,     0,     0,     0,   306,     0,
     225,     0,   240,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,   325,     0,     0,     0,     0,     0,   486,     0,     0,
       0,     0,     0,   558,     0,     0,     1,     4,     0,     0,
       0,     0,   360,   361,     0,   362,     0,   359,   363,     0,
       0,     0,    35,    37,    39,    38,    34,    36,   199,   198,
       0,   108,    26,   107,   111,   182,    28,   105,   106,     0,
       0,     0,     0,     0,     0,     0,     0,    70,     0,     0,
       0,     0,     0,    93,    95,    94,    71,     0,    96,     0,
       0,     0,     0,     0,     0,     0,    74,   187,    66,    68,
      69,    67,   183,   190,   187,     0,     0,   221,   235,    40,
     287,     0,   268,   270,   288,   285,   308,   244,     0,   242,
      29,   469,     0,   467,   204,   123,     0,     0,     0,    55,
     324,     0,   317,   205,   197,   185,   209,   210,   211,   212,
     207,   208,   206,     0,     0,     0,   480,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,   136,     0,   124,   125,
       0,   202,   201,   203,     0,   484,   195,    60,     0,   194,
     196,   200,   562,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,   161,     0,    61,    58,    59,     0,    72,    56,
      73,     0,    57,     0,    65,   213,    62,    63,     0,     0,
     342,     0,   357,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,   104,   103,     0,   189,   188,     0,
       0,     0,     0,     0,   184,     0,   186,     0,     0,   228,
     220,   222,    43,    48,    49,     0,   237,    41,    44,    45,
      42,   267,   269,     0,     0,     0,     0,   247,   241,   243,
       0,   466,   468,   122,     0,     0,     0,   482,   479,   481,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,   135,   137,
       0,     0,   133,     0,     0,     0,     0,     0,   563,     0,
       0,     0,   603,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,   160,
     162,     0,     0,    51,    53,     0,     0,   324,   337,     0,
     327,   334,   336,     0,     0,     0,     0,     0,   125,   109,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,    75,    98,
      99,   100,   101,   102,   191,     0,     0,     0,     0,   231,
       0,    46,    47,   286,   246,    40,   289,   271,     0,     0,
     278,     0,   280,   279,   281,     0,     0,     0,     0,     0,
       0,   305,     0,     0,   125,     0,     0,   474,   157,   138,
     145,   153,   159,   158,   140,   147,   148,   155,   154,   156,
     144,   141,   143,   139,   150,   146,   149,   142,   152,   151,
       0,    31,     0,   130,   487,     0,     0,     0,     0,   488,
     489,   490,   125,     0,     0,     0,     0,   560,   568,   567,
     569,   602,     0,     0,   604,     0,   181,   179,   180,   163,
     168,   169,   167,   166,   172,   171,   173,   174,   175,   177,
     164,   165,   178,   170,   176,     0,    64,     0,     0,     0,
       0,     0,   335,     0,     0,     0,     0,   444,     0,     0,
       0,     0,     0,     0,     0,     0,     0,   377,   378,     0,
       0,   327,   364,     0,   369,   379,   380,   381,   382,     0,
     393,     0,    91,    88,    87,    89,    83,    85,    76,    82,
      77,    78,     0,     0,     0,     0,    84,    90,    97,    86,
       0,   227,   226,     0,     0,   232,   233,    50,   290,   274,
     272,   273,     0,   275,     0,   293,   294,     0,   291,     0,
       0,     0,     0,     0,   256,     0,     0,     0,   254,   255,
     245,   248,   249,     0,   250,   259,   473,     0,     0,   134,
      30,     0,     0,   129,   127,   128,   126,     0,     0,     0,
       0,     0,   499,     0,   494,   496,   497,   498,     0,   565,
     566,   564,     0,     0,   559,   561,   606,   607,   605,     0,
       0,     0,     0,    52,   121,     0,     0,     0,     0,     0,
       0,   326,   333,   328,   329,     0,   327,   396,   397,     0,
       0,   327,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,   365,     0,     0,   403,   110,     0,
       0,     0,     0,   192,   230,   229,     0,     0,     0,   276,
     277,     0,     0,   282,   295,   296,   309,   315,   314,   313,
     310,   312,   311,     0,     0,     0,   258,   257,     0,   470,
     131,     0,   483,   477,   475,   476,   462,   492,   491,   493,
       0,     0,     0,   485,   495,   132,   570,     0,     0,     0,
       0,   572,   574,   575,   576,     0,   609,     0,     0,   612,
       0,   119,     0,   338,     0,     0,     0,   347,   348,   351,
     349,   346,   350,     0,   345,   352,   344,     0,   395,   398,
     447,   448,     0,     0,   387,   388,     0,   367,   368,     0,
       0,   389,   383,   307,   375,   374,   373,     0,   358,   366,
     371,   370,   372,   402,     0,   400,   327,    79,    80,    92,
      81,     0,     0,   223,   292,     0,     0,   304,   302,   303,
       0,   299,     0,   284,     0,    40,     0,     0,   262,     0,
     264,    40,   251,     0,     0,     0,   450,   501,   502,   503,
     504,   505,   518,     0,     0,   521,     0,     0,   509,     0,
     506,   556,   510,     0,   581,     0,   571,   573,   608,    65,
      32,   610,    33,     0,     0,     0,     0,     0,     0,   464,
     327,     0,   331,   330,     0,     0,     0,     0,   343,   446,
     449,   445,     0,   391,   376,   390,     0,   399,   401,     0,
       0,     0,   404,   405,   406,   193,   234,   300,   301,   297,
       0,   283,   253,   236,     0,   260,   263,   261,   265,   252,
       0,   478,     0,   456,     0,     0,     0,     0,   516,     0,
       0,     0,     0,     0,     0,     0,     0,   508,   598,   600,
     599,     0,   595,     0,     0,     0,   589,     0,   611,     0,
       0,     0,   601,     0,     0,   453,     0,     0,   353,   354,
     355,   356,     0,     0,     0,   407,   394,     0,   266,     0,
     437,   434,     0,   436,   433,   471,   418,   420,   421,   422,
       0,   423,     0,   463,     0,   458,     0,     0,     0,     0,
     511,   507,     0,     0,   523,     0,   557,   512,   513,   514,
     515,   594,   596,     0,   579,   577,   585,   584,   580,   579,
     588,   590,     0,   614,   613,   615,    54,     0,   403,   339,
     332,     0,   384,     0,     0,   439,     0,     0,   298,     0,
       0,   419,   474,     0,     0,     0,     0,   460,     0,   529,
       0,     0,     0,   531,   532,   533,   534,   535,   536,   520,
     517,     0,   524,   527,   525,   528,   500,   592,   593,   597,
       0,   583,   582,   586,   587,   591,   465,   454,     0,     0,
     438,   440,   441,   417,   414,   412,   413,   415,   416,   411,
     429,     0,   431,   409,   410,   426,   427,     0,     0,     0,
     432,     0,   457,     0,     0,   451,     0,   530,     0,     0,
     519,   522,   526,   578,   327,     0,     0,     0,     0,   408,
       0,     0,     0,   472,     0,   459,     0,     0,   543,   545,
     544,   546,   546,   546,     0,   385,     0,   442,   430,   428,
     425,   424,   435,   461,   452,     0,   539,     0,   537,     0,
     538,   455,     0,   474,   553,   555,   550,   552,   549,   554,
     551,   548,   546,   547,   546,   546,     0,     0,     0,   542,
     540,   541,   386,   443
};

  /* YYPGOTO[NTERM-NUM].  */
static const yytype_int16 yypgoto[] =
{
   -1041, -1041, -1041,  1255,  -721,     0,   -97,  -108,  -128,  -348,
    -778,   -95,  -492, -1041,   940,   943,   229, -1041,   728, -1041,
   -1041,  1189,  -143,  -138,  -139, -1041,   -61,  -496,  -122,  -140,
   -1041,  -116,   -10,  -112, -1041, -1041, -1041, -1041, -1041, -1041,
     459,   -99,  -148, -1041, -1041, -1041, -1041, -1041, -1041, -1041,
   -1041,  1055,   852,   102, -1041, -1041,  1018,   935, -1041,  1117,
   -1041, -1041, -1041, -1041, -1041, -1041, -1041, -1041, -1041,   -55,
   -1041, -1041, -1041, -1041, -1041, -1041, -1041, -1041,  -345, -1041,
   -1041, -1041, -1041,   -68, -1041,   443, -1041, -1041, -1041, -1041,
   -1041, -1041,   818, -1041, -1041,  -532, -1041, -1041,  1120,  -354,
    -102, -1041, -1041, -1041, -1041, -1041, -1041, -1041, -1041,   953,
   -1041, -1041, -1041,   567, -1041, -1041, -1041,   392,  -285, -1041,
   -1041, -1041, -1041, -1041, -1041, -1041, -1041,  -291,  -106,    -4,
    -735,  -606, -1041, -1041,  1224, -1041,   895, -1041, -1041, -1041,
   -1041, -1041, -1041,  -625, -1041, -1041, -1041,   762, -1041,  -554,
   -1041, -1041, -1041, -1041, -1041, -1041, -1041,   496, -1041, -1041,
   -1041,  1028,   620, -1041, -1041, -1041,   499,   294, -1041, -1041,
   -1041, -1041, -1041, -1020,  -982, -1041, -1041, -1041,  -534,   205,
   -1041, -1041, -1041, -1041, -1041, -1041,   293, -1041, -1041, -1041,
   -1041, -1041,   526, -1041, -1041, -1041, -1041, -1041, -1041, -1041,
   -1041, -1041, -1041, -1041, -1041, -1041, -1041,  1139, -1041, -1041,
   -1041,   861, -1040, -1041, -1041, -1041, -1041,  1118, -1041, -1041,
   -1041, -1041, -1041, -1041, -1041, -1041, -1041,   687, -1041, -1041,
   -1041, -1041, -1041, -1041,   414,  -528, -1041, -1041, -1041, -1041,
   -1041, -1041, -1041,   358, -1041, -1041, -1041, -1041, -1041, -1041,
    -514,  -673, -1041, -1041, -1041, -1041, -1041, -1041, -1041,   826,
   -1041, -1041, -1041, -1041,   588, -1041,   345, -1041, -1041, -1041,
   -1041, -1041, -1041,   413, -1041, -1041, -1041,   419,  -645, -1041,
   -1041, -1041,   981,   603, -1041, -1041, -1041, -1041, -1041
};

  /* YYDEFGOTO[NTERM-NUM].  */
static const yytype_int16 yydefgoto[] =
{
      -1,    48,    49,    50,   173,   206,   175,   253,   835,   920,
     921,   528,   376,   377,   378,   379,   380,   452,   453,    52,
      53,    54,   877,   551,   326,   327,   317,   208,   209,   878,
     210,   211,   278,   177,   178,    55,    56,    57,   725,    58,
     236,   279,   421,   696,    59,    60,    61,    62,    63,    64,
     275,   276,   529,   534,    65,   311,   312,   579,    66,   364,
     214,    67,    68,    69,    70,    71,    72,    73,   216,   104,
     105,   111,   369,   499,   656,   767,   219,   886,   220,    74,
      75,    76,   228,   112,   387,   505,   522,   681,   682,   683,
     788,   684,   793,   887,   889,   888,    77,   221,   222,   768,
     510,   511,   512,   513,   514,   515,   883,   223,   224,   225,
     385,   506,   667,   668,   773,   774,   775,   880,   881,    78,
     109,   854,   386,   779,    79,   120,    80,    81,    82,   329,
     732,   601,   733,   734,    83,   460,   461,   462,   930,    84,
      85,   463,   604,   836,    86,   466,   160,   621,   862,   622,
     623,   624,   625,   626,   627,   628,   850,   851,    87,    88,
     757,   465,   738,   739,   630,   864,   865,   866,   952,   953,
    1077,  1131,  1132,  1025,  1026,  1027,  1028,  1134,  1135,  1136,
    1029,  1030,  1031,  1032,   954,  1074,  1075,  1157,  1193,    89,
     741,   607,   841,   842,    90,   973,  1167,    91,  1068,  1154,
    1035,  1087,  1145,   896,  1005,    92,   232,   233,   390,   893,
    1082,  1076,   691,   794,   795,    93,   255,   256,   527,    94,
     424,   285,   558,   559,   560,   561,   703,   704,   705,   802,
     900,   706,   707,   909,   910,   911,   912,   974,   977,  1045,
    1106,  1090,  1091,  1092,  1093,  1094,  1095,  1096,  1097,  1098,
    1171,  1186,  1203,   987,    95,   292,   566,   427,   428,   567,
     568,   569,   570,   810,   811,   812,  1111,   994,  1058,  1059,
    1115,   813,   995,   996,  1109,   814,   991,   992,   993,    96,
     294,   431,   432,   717,   718,   575,   721,   819,   927
};

  /* YYTABLE[YYPACT[STATE-NUM]] -- What to do in state STATE-NUM.  If
     positive, shift that token.  If negative, reduce the rule whose
     number is the opposite.  If YYTABLE_NINF, syntax error.  */
static const yytype_int16 yytable[] =
{
      51,   254,   168,   169,   212,   315,   231,   103,   319,   316,
     108,   110,   325,   658,   241,   753,  1054,   988,   858,  1107,
    1113,   827,   576,   891,   564,   700,   318,  -120,   245,   132,
     700,   509,   320,   252,   988,   280,   322,   170,   458,   207,
     828,   213,  1139,  1081,   988,   229,   674,   932,   577,    51,
     807,   217,   171,   151,   234,   357,    44,  1133,   676,   677,
     243,   244,    44,   246,   247,   248,   249,   754,    44,   564,
     573,   152,   550,    44,   282,   283,   949,   286,   675,   289,
     290,   291,   555,   565,   808,   106,   709,   829,   370,   176,
     153,   176,   516,    97,   430,   833,   676,   677,   710,   174,
     665,   174,   938,   663,   786,   381,   103,   218,   941,   239,
      44,   107,  -392,   110,   230,   807,   235,   235,   235,   174,
     240,  1176,   367,   368,   231,   358,   532,   254,   565,   230,
     837,   701,   787,   174,   174,   843,   701,  1177,   321,   287,
     791,   131,   354,   355,   356,   998,    98,   328,   314,   808,
     324,  1081,   330,  1207,   659,   331,   332,   333,   830,   989,
     389,   371,   676,   677,   662,   172,   664,    44,   809,   458,
     666,   162,   676,   677,   990,   154,   989,   831,   791,   832,
     357,  1194,   754,   517,  1195,  1196,   989,   151,   556,    44,
      44,   990,    44,   792,    99,  1081,  1197,  1198,   -60,   388,
     557,   990,   702,   803,   578,   152,   155,   702,   934,   711,
     950,   172,    44,   156,   916,   714,   103,    44,   335,   100,
     518,   218,   678,   459,   153,   924,   455,   157,   110,   158,
      44,  1163,   230,   809,   648,   468,    44,   163,    44,  1199,
      44,   101,   172,    44,   791,   833,   507,   164,   165,   391,
     358,   726,   454,   890,   102,   230,   313,   737,  1200,  1201,
     951,   508,   519,   520,   313,   113,   313,    44,  1055,   679,
     997,  1069,   159,  1056,  1051,   470,   471,   472,   473,   474,
     475,   476,   477,   478,   479,   480,   481,   482,   483,   484,
     485,  1060,   487,   497,   498,   524,   525,  1213,   489,   490,
     491,   492,   493,   521,   494,   660,   908,   680,   533,   600,
    1009,  1010,  1011,   176,   500,  1088,   383,   533,  1088,   661,
     631,   533,   533,   174,  1006,   469,   562,    44,   166,   457,
     783,   925,   114,   384,   174,   174,   359,   360,   361,   362,
     155,   363,   486,   926,   796,   537,   115,   156,   167,   542,
     997,   162,   174,   313,   459,   966,   980,   968,   981,   838,
      44,   157,   554,   158,   359,   360,   361,   362,    44,   363,
     728,   816,   729,    44,   116,   582,   688,   584,   978,   979,
     587,    44,   117,   218,   504,   728,   313,   729,   902,    44,
     595,   737,   939,   963,   174,   174,   526,   847,   771,   969,
     398,   848,   849,   856,   772,   944,   159,    44,   730,   693,
     908,   908,   694,   965,   708,  1073,    44,   163,   903,   571,
     324,   967,    44,   730,   118,   174,   666,   164,   165,   174,
     318,   119,  1184,   313,   651,   652,   320,   660,   552,  1191,
     695,   121,   553,   892,   187,    44,   928,   731,  1038,   295,
    1042,   661,   122,  1043,  1120,  1047,  1048,  1049,  1050,   552,
     457,   133,  1016,   758,   359,   360,   361,   362,  1044,   363,
     296,   313,   361,   362,   929,   363,   931,   947,   933,   552,
     908,   552,   123,   790,   908,   805,   297,   908,   908,   908,
     908,   257,   727,   359,   360,   361,   362,   162,   363,  1188,
    1190,   943,    44,   258,   596,   849,   124,   597,   166,   218,
     125,   218,   259,   722,   544,   218,   126,  1168,   548,   549,
     127,   298,   299,   686,   196,   373,   374,   959,   167,  1209,
     960,  1210,  1211,   128,  1180,   260,   454,  1160,  1181,   300,
     301,  1160,   321,   134,   359,   360,   361,   362,  1174,   363,
     985,   986,   314,   261,   262,   302,   303,   304,   904,   797,
     769,   799,   770,   163,   935,   936,   937,   305,   129,   306,
     307,  1071,   716,   164,   165,   905,   237,   238,  1008,   263,
      44,   759,   760,   761,   762,   308,   309,   130,   724,   763,
     957,   958,   764,   765,   372,   373,   374,   375,   174,   324,
     230,   906,   840,  1161,  1162,   324,   264,   131,   776,   778,
     780,   781,   782,   907,   135,   861,   860,   137,   863,   136,
     359,   360,   361,   362,   138,   363,   139,   834,   935,   936,
     937,   265,  1072,   318,  1172,  1173,   140,   798,   250,   320,
     141,  1019,  1020,  1021,   142,   266,   267,   268,   143,   144,
     318,   145,   146,   269,   166,   148,   320,   149,   270,   150,
     879,   820,   161,   215,   310,   226,  1022,   251,  1023,  1024,
     983,   984,   985,   986,   167,   218,   227,   922,   277,   284,
     293,   844,   845,   335,  1040,   922,   323,   852,   334,   855,
     336,   394,   324,   337,   271,   318,  1158,   324,   449,   324,
     338,   320,   806,   318,   840,   898,   901,   339,   340,   320,
     341,   272,   174,   342,   273,   343,   344,   274,  1194,   324,
     345,  1195,  1196,   823,   346,   395,   863,   359,   360,   361,
     362,   347,   363,  1197,  1198,   230,  1141,   348,   349,   324,
     324,   350,   418,   393,   351,   324,   218,   218,   352,   853,
     353,   324,   396,   318,   318,   859,   324,   154,   397,   320,
     320,   400,   321,   879,   879,   257,   983,   984,   985,   986,
     401,   218,   314,   402,   403,   404,  1199,   258,   425,   405,
    1100,   406,   407,   922,   885,   885,   259,   408,   218,   409,
     899,   410,   411,   412,   413,  1200,  1201,   414,   415,   416,
     174,   318,   417,   359,   360,   361,   362,   320,   363,   260,
     955,   956,   420,   642,   423,   716,  1182,   919,   422,   426,
     429,   430,   324,   962,   324,   919,   324,   261,   262,   359,
     360,   361,   362,   230,   363,   433,   434,  1202,   318,  1204,
     435,  1205,   324,   488,   320,   436,   437,   438,   879,   324,
     218,  -105,   439,   263,   359,   360,   361,   362,   440,   363,
     456,  1063,   441,  1065,   324,   321,   321,  1003,   632,   359,
     360,   361,   362,   442,   363,   314,   314,  1033,   443,   444,
     264,  1012,   961,   633,   445,  1013,   446,   885,   447,   885,
     359,   360,   361,   362,   504,   363,   359,   360,   361,   362,
     448,   363,   464,  1018,   634,   265,   359,   360,   361,   362,
     635,   363,   363,   495,  1037,  1108,   496,   523,   657,   266,
     267,   268,   563,   919,  1114,   598,   572,   269,   359,   360,
     361,   362,   270,   363,   599,   230,   230,   230,  -392,  1064,
     603,   605,   636,   606,   629,   653,  1070,   654,  1105,   655,
     321,   669,   359,   360,   361,   362,   608,   363,   507,   689,
     314,   692,   359,   360,   361,   362,   637,   363,   271,   670,
     671,   672,  1046,   609,   673,   690,   638,   713,   359,   360,
     361,   362,   174,   363,   766,   272,   697,   698,   273,   699,
    1142,   274,   639,   712,  1057,   359,   360,   361,   362,   324,
     363,   324,   719,   720,   610,   735,  1116,   660,   740,   640,
     359,   360,   361,   362,   686,   363,  1130,   742,  1170,   743,
    1170,  1170,   745,   744,   641,   747,   746,   611,   359,   360,
     361,   362,   748,   363,   749,   176,  1089,   750,  1169,  1089,
    1169,  1169,   646,   751,   179,   174,   612,   180,  1165,  1153,
     613,   755,   752,   324,   756,   784,   785,  1155,  1156,   181,
     789,   182,   324,   800,   614,   801,   817,   815,   818,   183,
     295,  1183,   821,   822,   686,   846,  1130,  1130,   824,   825,
     826,   184,   857,  1140,   324,   359,   360,   361,   362,   872,
     363,   296,   873,   882,   884,   643,   894,  1178,   895,  1130,
     913,   615,   616,   359,   360,   361,   362,   297,   363,   914,
     923,   185,   915,   617,   618,   942,   162,   647,   946,   186,
     964,   170,   970,   619,   170,   972,   975,   359,   360,   361,
     362,  1206,   363,   999,  1004,   976,   171,  1000,   187,   171,
    1001,   649,   298,   299,  1002,  1007,   230,  1034,   230,   230,
    1014,   188,   189,   190,   982,   983,   984,   985,   986,   620,
     300,   301,   191,  1015,   324,  1036,   192,   359,   360,   361,
     362,  1039,   363,  1053,  1062,  1067,   302,   303,   304,   937,
    1078,   867,   163,  1079,   897,   324,  1122,   324,   305,   324,
     306,   307,   164,   165,   359,   360,   361,   362,  1080,   363,
     359,   360,   361,   362,   193,   363,   308,   309,   868,  1083,
    1086,  1084,  1175,  1085,   194,   195,  1137,  1102,   196,   359,
     360,   361,   362,  1103,   363,   179,  1104,  1110,   180,  1118,
     197,  1123,   198,   869,  -111,   199,   -74,   -74,   -74,   -74,
     181,   -74,   182,  1124,   200,  1119,  1125,  1143,   201,   172,
     183,  1138,   172,   530,   531,   202,   535,   536,   538,   539,
     540,   541,   184,   543,  1144,   545,   546,   547,  1146,  1152,
     375,  1148,  1149,   166,    44,  1150,  1151,    44,  1159,   359,
     360,   361,   362,  1160,   363,   310,  1164,  1166,  1019,  1020,
    1021,  1212,   185,   167,   203,   204,   359,   360,   361,   362,
     186,   363,  1185,   731,   147,  1126,   644,   205,  1187,   875,
     876,  1189,  1192,  1022,  1208,  1023,  1024,   501,  1147,   187,
     502,   359,   360,   361,   362,   723,   363,  1127,   288,   450,
     419,   366,   188,   189,   190,   870,   503,   971,   874,   179,
     685,   382,   180,   191,   242,  1123,   945,   192,   359,   360,
     361,   362,  1017,   363,   181,   602,   182,  1124,   645,   839,
    1125,   467,  1117,   948,   183,  1179,   736,  1121,   940,   580,
     581,   392,   583,   399,   585,   586,   184,   588,   589,   590,
     591,   592,   593,   594,   687,   193,   359,   360,   361,   362,
     804,   363,   715,  1066,  1041,   194,   195,  1101,   917,   196,
     359,   360,   361,   362,  1112,   363,   185,  1099,  1061,  1128,
    1052,   197,   574,   198,   186,     0,   199,     0,   918,  1126,
     365,   359,   360,   361,   362,   200,   363,   179,     0,   201,
     180,     0,     0,   187,     0,     0,   202,     0,     0,     0,
       0,  1127,   181,     0,   182,  1129,   188,   189,   190,     0,
       0,     0,   183,     0,     0,    44,     0,   191,     0,     0,
       0,   192,     0,     0,   184,   650,   359,   360,   361,   362,
       0,   363,     0,     0,     0,   203,   204,   871,   359,   360,
     361,   362,     0,   363,   359,   360,   361,   362,   205,   363,
     -68,   -68,   -68,   -68,   185,   -68,   162,     0,     0,   193,
       0,   179,   186,     0,   180,     0,     0,     0,     0,   194,
     195,     0,     0,   196,     0,     0,   181,     0,   182,     0,
       0,   187,     0,  1128,     0,   197,   183,   198,     0,     0,
     199,     0,     0,     0,   188,   189,   190,     0,   184,   200,
       0,     0,     0,   201,     0,   191,     0,     0,     0,   192,
     202,   -67,   -67,   -67,   -67,     0,   -67,     0,     0,  1129,
       0,     0,   163,     0,     0,     0,     0,     0,   185,    44,
       0,     0,   164,   165,     0,     0,   186,     0,   170,   -74,
     -74,   -74,   -74,     0,   -74,     0,     0,   193,     0,   203,
     204,     0,     0,   171,     0,   187,     0,   194,   195,     0,
       0,   196,   205,     0,     0,     0,     0,     0,   188,   189,
     190,     0,     0,   197,     0,   198,     0,     0,   199,   191,
       0,     0,     0,   192,     0,     0,     0,   200,     0,     0,
       0,   201,     0,     0,     0,     0,     0,     0,   202,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,   179,     0,   166,   180,     0,     0,    44,     0,     0,
       0,   193,     0,     0,     0,     0,   181,     0,   182,     0,
       0,   194,   195,   167,     0,   196,   183,   203,   204,     0,
       0,     0,     0,     0,     0,     0,     0,   197,   184,   198,
     205,     0,   199,     0,     0,     0,     0,     0,     0,     0,
       0,   200,     0,     0,     0,   201,   172,     0,     0,     0,
       0,     0,   202,     0,     0,     0,     0,     0,   185,     0,
       0,     0,     0,     0,     0,     0,   186,     0,   170,     0,
       0,    44,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,   171,     0,   187,     0,     0,     0,     0,
       0,   203,   204,     0,     0,     0,     0,   313,   188,   189,
     190,     0,     0,     0,   205,   179,     0,     0,   180,   191,
       0,   179,     0,   192,   180,     0,     0,     0,     0,     0,
     181,     0,   182,     0,     0,     0,   181,     0,   182,     0,
     183,     0,     0,     0,     0,     0,   183,     0,     0,     0,
       0,     0,   184,     0,     0,     0,     0,     0,   184,     0,
       0,   193,     0,     0,     0,     0,     0,     0,     0,     0,
       0,   194,   195,     0,     0,   196,     0,     0,     0,     0,
       0,     0,   185,     0,     0,     0,     0,   197,   185,   198,
     186,     0,   199,     0,     0,     0,   186,     0,     0,     0,
       0,   200,     0,     0,     0,   201,   172,     0,     0,   187,
     281,     0,   202,     0,     0,   187,     0,     0,     0,     0,
       0,     0,   188,   189,   190,     0,     0,     0,   188,   189,
     190,    44,     0,   191,     0,     0,     0,   192,     0,   191,
       0,   179,     0,   192,   180,     0,     0,     0,     0,     0,
       0,   203,   204,     0,     0,     0,   181,     0,   182,     0,
       0,     0,     0,     0,   205,     0,   183,     0,     0,     0,
       0,     0,     0,     0,     0,   193,     0,     0,   184,     0,
       0,   193,     0,     0,     0,   194,   195,     0,     0,   196,
       0,   194,   195,     0,     0,   196,     0,     0,     0,     0,
       0,   197,     0,   198,     0,     0,   199,   197,   185,   198,
       0,     0,   199,     0,     0,   200,   186,     0,     0,   201,
       0,   200,     0,     0,     0,   201,   202,     0,     0,     0,
       0,     0,   202,     0,     0,   187,     0,     0,     0,   179,
       0,     0,   180,     0,     0,    44,     0,     0,   188,   189,
     190,    44,     0,     0,   181,     0,   182,     0,     0,   191,
       0,     0,     0,   192,   183,   203,   204,     0,     0,     0,
       0,   203,   204,     0,     0,     0,   184,   451,   205,     0,
       0,     0,     0,     0,   205,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,   193,     0,     0,     0,     0,   185,     0,     0,     0,
       0,   194,   195,   179,   186,   196,   180,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,   197,   181,   198,
     182,     0,   199,   187,     0,     0,     0,     0,   183,     0,
       0,   200,     0,     0,     0,   201,   188,   189,   190,     0,
     184,     0,   202,     0,     0,     0,     0,   191,     0,     0,
       0,   192,     0,     0,     0,     0,     0,     0,   777,     0,
       0,    44,     0,     0,     0,     0,     0,     0,     0,     0,
     185,     0,     0,     0,     0,     0,     0,     0,   186,     0,
       0,   203,   204,     0,     0,     0,     0,     0,     0,   193,
       0,     0,     0,     0,   205,     0,     0,   187,     0,   194,
     195,     0,     0,   196,     0,     0,     0,     0,     0,     0,
     188,   189,   190,     0,     0,   197,     0,   198,     0,     0,
     199,   191,     0,     0,     0,   192,     0,     0,     0,   200,
       0,     0,     0,   201,     0,     0,     0,     0,     0,     0,
     202,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,    44,
       0,     0,     0,   193,     0,     0,     0,     0,     0,     0,
       0,     0,     0,   194,   195,     0,     0,   196,     0,   203,
     204,     0,     0,     1,     0,   313,     0,     0,     0,   197,
       0,   198,   205,     0,   199,     0,     0,     0,     0,     0,
       0,     0,     0,   200,     0,     0,     0,   201,     2,     3,
       4,     5,     6,     0,   202,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     7,     8,     9,    10,    11,    12,
      13,     0,     0,    44,     0,     0,     0,    14,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,    15,     0,   203,   204,     0,     0,     0,     0,    16,
      17,     0,     0,     0,     0,     0,   205,     0,     0,     0,
      18,     0,     0,     0,    19,     0,     0,    20,     0,     0,
       0,    21,    22,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,    23,    24,    25,    26,    27,     0,
       0,     0,     0,     0,     0,    28,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,    29,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,    30,    31,    32,    33,     0,     0,     0,
       0,     0,     0,     0,    34,    35,     0,     0,     0,    36,
      37,     0,     0,    38,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,    39,     0,     0,     0,     0,
      40,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,    41,    42,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
      43,    44,     0,     0,    45,     0,     0,    46,     0,     0,
       0,     0,     0,     0,     0,    47
};

static const yytype_int16 yycheck[] =
{
       0,   129,    97,    98,   101,   148,   114,     7,   148,   148,
      10,    11,   150,   505,   120,   621,    13,     8,   753,    12,
      15,    42,    64,    13,    74,    87,   148,    79,   123,    33,
      87,   385,   148,   128,     8,   134,   148,    80,   329,   100,
      61,   102,  1082,  1025,     8,   113,     4,   825,    90,    49,
     101,   106,    95,    25,   115,    69,   233,  1077,   130,   131,
     121,   122,   233,   124,   125,   126,   127,   621,   233,    74,
     109,    43,   420,   233,   135,   136,    52,   138,    36,   140,
     141,   142,   100,   133,   135,   264,    17,   108,   265,    99,
      62,   101,    54,   250,   133,   266,   130,   131,    29,    99,
     141,   101,   837,   268,   222,   265,   106,   107,   843,   119,
     233,   264,   164,   113,   114,   101,   116,   117,   118,   119,
     120,    38,    56,    57,   232,   139,    27,   255,   133,   129,
     736,   193,   250,   133,   134,   741,   193,  1157,   148,   139,
     212,   264,   203,   204,   205,   923,   250,   151,   148,   135,
     150,  1133,   156,  1193,   508,   159,   160,   161,   179,   150,
     228,   216,   130,   131,   509,   208,   511,   233,   219,   460,
     515,    72,   130,   131,   165,   147,   150,   198,   212,   200,
      69,     7,   736,   145,    10,    11,   150,    25,   206,   233,
     233,   165,   233,   265,   250,  1177,    22,    23,   250,   265,
     218,   165,   264,   265,   246,    43,   178,   264,   833,   140,
     186,   208,   233,   185,   265,   265,   216,   233,   251,   250,
     182,   221,   180,   329,    62,   100,   323,   199,   228,   201,
     233,   265,   232,   219,   267,   334,   233,   138,   233,    65,
     233,   250,   208,   233,   212,   266,   254,   148,   149,   265,
     139,   599,   313,   785,   250,   255,   259,   605,    84,    85,
     866,   269,   224,   225,   259,   264,   259,   233,   265,   227,
     915,  1006,   244,   994,   265,   336,   337,   338,   339,   340,
     341,   342,   343,   344,   345,   346,   347,   348,   349,   350,
     351,   265,   353,    44,    45,   394,   395,   265,   359,   360,
     361,   362,   363,   265,   365,   254,   802,   265,   403,   147,
     935,   936,   937,   323,   375,  1036,   253,   412,  1039,   268,
     468,   416,   417,   323,   930,   335,   425,   233,   229,   329,
     675,   206,   264,   270,   334,   335,   253,   254,   255,   256,
     178,   258,   352,   218,   692,   406,   250,   185,   249,   410,
     995,    72,   352,   259,   460,   887,   263,   889,   265,   265,
     233,   199,   423,   201,   253,   254,   255,   256,   233,   258,
     188,   719,   190,   233,   266,   436,   524,   438,   906,   907,
     441,   233,   266,   383,   384,   188,   259,   190,    41,   233,
     451,   739,   265,   885,   394,   395,   396,   745,   253,   891,
     265,   746,   747,   751,   259,   265,   244,   233,   226,   552,
     906,   907,   552,   265,   562,   204,   233,   138,    71,   429,
     420,   265,   233,   226,   266,   425,   771,   148,   149,   429,
     552,   250,  1167,   259,   495,   496,   552,   254,   263,  1174,
     552,   250,   267,   788,    97,   233,   206,   265,   976,    16,
     250,   268,   250,   253,   265,   983,   984,   985,   986,   263,
     460,   266,   265,   267,   253,   254,   255,   256,   268,   258,
      37,   259,   255,   256,   822,   258,   824,   265,   826,   263,
     976,   263,   250,   267,   980,   267,    53,   983,   984,   985,
     986,    14,   600,   253,   254,   255,   256,    72,   258,  1172,
    1173,   849,   233,    26,   260,   850,   250,   263,   229,   509,
     250,   511,    35,   222,   412,   515,   250,   248,   416,   417,
     250,    88,    89,   523,   177,   262,   263,   260,   249,  1202,
     263,  1204,  1205,   250,   260,    58,   597,   263,   260,   106,
     107,   263,   552,   266,   253,   254,   255,   256,  1154,   258,
     255,   256,   552,    76,    77,   122,   123,   124,   211,   697,
     662,   699,   664,   138,   253,   254,   255,   134,   264,   136,
     137,   204,   572,   148,   149,   228,   117,   118,   267,   102,
     233,   642,   643,   644,   645,   152,   153,   264,   598,   650,
     875,   876,   653,   654,   261,   262,   263,   264,   598,   599,
     600,   254,   740,  1137,  1138,   605,   129,   264,   669,   670,
     671,   672,   673,   266,   250,   755,   755,   264,   756,   250,
     253,   254,   255,   256,   250,   258,   266,   735,   253,   254,
     255,   154,   265,   755,  1148,  1149,   250,   698,   213,   755,
     250,    66,    67,    68,   250,   168,   169,   170,   264,   264,
     772,   264,     0,   176,   229,   250,   772,    79,   181,   250,
     772,   722,   164,   264,   231,   264,    91,   242,    93,    94,
     253,   254,   255,   256,   249,   675,   264,   817,    39,   151,
      75,   742,   743,   251,   267,   825,   266,   748,   266,   750,
     266,   263,   692,   266,   217,   817,   222,   697,   265,   699,
     266,   817,   712,   825,   842,   800,   801,   266,   266,   825,
     266,   234,   712,   266,   237,   266,   266,   240,     7,   719,
     266,    10,    11,   727,   266,   263,   864,   253,   254,   255,
     256,   266,   258,    22,    23,   735,  1084,   266,   266,   739,
     740,   266,   265,   267,   266,   745,   746,   747,   266,   749,
     266,   751,   259,   875,   876,   755,   756,   147,   264,   875,
     876,   250,   772,   875,   876,    14,   253,   254,   255,   256,
     250,   771,   772,   250,   250,   250,    65,    26,   263,   250,
     267,   250,   250,   923,   784,   785,    35,   250,   788,   250,
     800,   250,   250,   250,   250,    84,    85,   250,   250,   250,
     800,   923,   250,   253,   254,   255,   256,   923,   258,    58,
     871,   872,   250,   263,   250,   815,  1164,   817,   267,   126,
     250,   133,   822,   884,   824,   825,   826,    76,    77,   253,
     254,   255,   256,   833,   258,   250,   250,  1185,   960,  1187,
     250,  1189,   842,   267,   960,   250,   250,   250,   960,   849,
     850,   251,   250,   102,   253,   254,   255,   256,   250,   258,
      39,   999,   250,  1001,   864,   875,   876,   928,   267,   253,
     254,   255,   256,   250,   258,   875,   876,   972,   250,   250,
     129,   942,   882,   267,   250,   946,   250,   887,   250,   889,
     253,   254,   255,   256,   894,   258,   253,   254,   255,   256,
     250,   258,   235,   964,   267,   154,   253,   254,   255,   256,
     267,   258,   258,   250,   975,  1053,   250,   259,   265,   168,
     169,   170,   250,   923,  1062,   263,   250,   176,   253,   254,
     255,   256,   181,   258,   250,   935,   936,   937,   164,  1000,
     192,   264,   267,   216,    60,   250,  1007,   250,  1045,   215,
     960,   250,   253,   254,   255,   256,    38,   258,   254,   267,
     960,   263,   253,   254,   255,   256,   267,   258,   217,   250,
     250,   250,   982,    55,   250,   260,   267,   264,   253,   254,
     255,   256,   982,   258,   113,   234,   250,   250,   237,   250,
    1085,   240,   267,   250,   994,   253,   254,   255,   256,   999,
     258,  1001,   250,   245,    86,   250,  1067,   254,   264,   267,
     253,   254,   255,   256,  1014,   258,  1077,   250,  1146,   250,
    1148,  1149,   250,   264,   267,   264,   250,   109,   253,   254,
     255,   256,   250,   258,   250,  1045,  1036,   250,  1146,  1039,
    1148,  1149,   267,   250,     3,  1045,   128,     6,  1143,  1110,
     132,   250,   264,  1053,   264,   264,   264,  1118,  1119,    18,
     260,    20,  1062,   250,   146,   250,   250,   253,   243,    28,
      16,  1166,   267,   263,  1074,   115,  1137,  1138,   250,   250,
     250,    40,   114,  1083,  1084,   253,   254,   255,   256,   250,
     258,    37,   265,   252,   250,   263,   250,  1158,   162,  1160,
     264,   183,   184,   253,   254,   255,   256,    53,   258,   264,
     250,    70,   264,   195,   196,   250,    72,   267,   250,    78,
     250,    80,   264,   205,    80,   250,   259,   253,   254,   255,
     256,  1192,   258,   250,    19,   266,    95,   250,    97,    95,
     250,   267,    88,    89,   265,   263,  1146,   159,  1148,  1149,
     264,   110,   111,   112,   252,   253,   254,   255,   256,   241,
     106,   107,   121,   264,  1164,   259,   125,   253,   254,   255,
     256,   259,   258,   270,   270,   250,   122,   123,   124,   255,
     260,   267,   138,   250,   140,  1185,   264,  1187,   134,  1189,
     136,   137,   148,   149,   253,   254,   255,   256,   250,   258,
     253,   254,   255,   256,   163,   258,   152,   153,   267,   250,
     160,   250,   265,   250,   173,   174,   259,   268,   177,   253,
     254,   255,   256,   268,   258,     3,   268,   250,     6,   250,
     189,     9,   191,   267,   251,   194,   253,   254,   255,   256,
      18,   258,    20,    21,   203,   250,    24,   250,   207,   208,
      28,   259,   208,   401,   402,   214,   404,   405,   406,   407,
     408,   409,    40,   411,   161,   413,   414,   415,   263,   268,
     264,   263,   263,   229,   233,   260,   260,   233,   265,   253,
     254,   255,   256,   263,   258,   231,   263,   250,    66,    67,
      68,   265,    70,   249,   253,   254,   253,   254,   255,   256,
      78,   258,   263,   265,    49,    83,   263,   266,   263,   268,
     269,   263,   250,    91,   263,    93,    94,   377,  1089,    97,
     377,   253,   254,   255,   256,   597,   258,   105,   139,   311,
     275,   214,   110,   111,   112,   267,   383,   894,   771,     3,
     522,   221,     6,   121,   120,     9,   850,   125,   253,   254,
     255,   256,   960,   258,    18,   460,    20,    21,   263,   739,
      24,   333,  1068,   864,    28,  1160,   604,  1074,   842,   434,
     435,   232,   437,   255,   439,   440,    40,   442,   443,   444,
     445,   446,   447,   448,   523,   163,   253,   254,   255,   256,
     703,   258,   566,   260,   980,   173,   174,  1039,   810,   177,
     253,   254,   255,   256,  1059,   258,    70,   260,   995,   187,
     991,   189,   431,   191,    78,    -1,   194,    -1,   815,    83,
     252,   253,   254,   255,   256,   203,   258,     3,    -1,   207,
       6,    -1,    -1,    97,    -1,    -1,   214,    -1,    -1,    -1,
      -1,   105,    18,    -1,    20,   223,   110,   111,   112,    -1,
      -1,    -1,    28,    -1,    -1,   233,    -1,   121,    -1,    -1,
      -1,   125,    -1,    -1,    40,   252,   253,   254,   255,   256,
      -1,   258,    -1,    -1,    -1,   253,   254,   252,   253,   254,
     255,   256,    -1,   258,   253,   254,   255,   256,   266,   258,
     253,   254,   255,   256,    70,   258,    72,    -1,    -1,   163,
      -1,     3,    78,    -1,     6,    -1,    -1,    -1,    -1,   173,
     174,    -1,    -1,   177,    -1,    -1,    18,    -1,    20,    -1,
      -1,    97,    -1,   187,    -1,   189,    28,   191,    -1,    -1,
     194,    -1,    -1,    -1,   110,   111,   112,    -1,    40,   203,
      -1,    -1,    -1,   207,    -1,   121,    -1,    -1,    -1,   125,
     214,   253,   254,   255,   256,    -1,   258,    -1,    -1,   223,
      -1,    -1,   138,    -1,    -1,    -1,    -1,    -1,    70,   233,
      -1,    -1,   148,   149,    -1,    -1,    78,    -1,    80,   253,
     254,   255,   256,    -1,   258,    -1,    -1,   163,    -1,   253,
     254,    -1,    -1,    95,    -1,    97,    -1,   173,   174,    -1,
      -1,   177,   266,    -1,    -1,    -1,    -1,    -1,   110,   111,
     112,    -1,    -1,   189,    -1,   191,    -1,    -1,   194,   121,
      -1,    -1,    -1,   125,    -1,    -1,    -1,   203,    -1,    -1,
      -1,   207,    -1,    -1,    -1,    -1,    -1,    -1,   214,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,     3,    -1,   229,     6,    -1,    -1,   233,    -1,    -1,
      -1,   163,    -1,    -1,    -1,    -1,    18,    -1,    20,    -1,
      -1,   173,   174,   249,    -1,   177,    28,   253,   254,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,   189,    40,   191,
     266,    -1,   194,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,   203,    -1,    -1,    -1,   207,   208,    -1,    -1,    -1,
      -1,    -1,   214,    -1,    -1,    -1,    -1,    -1,    70,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    78,    -1,    80,    -1,
      -1,   233,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    95,    -1,    97,    -1,    -1,    -1,    -1,
      -1,   253,   254,    -1,    -1,    -1,    -1,   259,   110,   111,
     112,    -1,    -1,    -1,   266,     3,    -1,    -1,     6,   121,
      -1,     3,    -1,   125,     6,    -1,    -1,    -1,    -1,    -1,
      18,    -1,    20,    -1,    -1,    -1,    18,    -1,    20,    -1,
      28,    -1,    -1,    -1,    -1,    -1,    28,    -1,    -1,    -1,
      -1,    -1,    40,    -1,    -1,    -1,    -1,    -1,    40,    -1,
      -1,   163,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,   173,   174,    -1,    -1,   177,    -1,    -1,    -1,    -1,
      -1,    -1,    70,    -1,    -1,    -1,    -1,   189,    70,   191,
      78,    -1,   194,    -1,    -1,    -1,    78,    -1,    -1,    -1,
      -1,   203,    -1,    -1,    -1,   207,   208,    -1,    -1,    97,
      98,    -1,   214,    -1,    -1,    97,    -1,    -1,    -1,    -1,
      -1,    -1,   110,   111,   112,    -1,    -1,    -1,   110,   111,
     112,   233,    -1,   121,    -1,    -1,    -1,   125,    -1,   121,
      -1,     3,    -1,   125,     6,    -1,    -1,    -1,    -1,    -1,
      -1,   253,   254,    -1,    -1,    -1,    18,    -1,    20,    -1,
      -1,    -1,    -1,    -1,   266,    -1,    28,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,   163,    -1,    -1,    40,    -1,
      -1,   163,    -1,    -1,    -1,   173,   174,    -1,    -1,   177,
      -1,   173,   174,    -1,    -1,   177,    -1,    -1,    -1,    -1,
      -1,   189,    -1,   191,    -1,    -1,   194,   189,    70,   191,
      -1,    -1,   194,    -1,    -1,   203,    78,    -1,    -1,   207,
      -1,   203,    -1,    -1,    -1,   207,   214,    -1,    -1,    -1,
      -1,    -1,   214,    -1,    -1,    97,    -1,    -1,    -1,     3,
      -1,    -1,     6,    -1,    -1,   233,    -1,    -1,   110,   111,
     112,   233,    -1,    -1,    18,    -1,    20,    -1,    -1,   121,
      -1,    -1,    -1,   125,    28,   253,   254,    -1,    -1,    -1,
      -1,   253,   254,    -1,    -1,    -1,    40,   259,   266,    -1,
      -1,    -1,    -1,    -1,   266,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,   163,    -1,    -1,    -1,    -1,    70,    -1,    -1,    -1,
      -1,   173,   174,     3,    78,   177,     6,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,   189,    18,   191,
      20,    -1,   194,    97,    -1,    -1,    -1,    -1,    28,    -1,
      -1,   203,    -1,    -1,    -1,   207,   110,   111,   112,    -1,
      40,    -1,   214,    -1,    -1,    -1,    -1,   121,    -1,    -1,
      -1,   125,    -1,    -1,    -1,    -1,    -1,    -1,   230,    -1,
      -1,   233,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      70,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    78,    -1,
      -1,   253,   254,    -1,    -1,    -1,    -1,    -1,    -1,   163,
      -1,    -1,    -1,    -1,   266,    -1,    -1,    97,    -1,   173,
     174,    -1,    -1,   177,    -1,    -1,    -1,    -1,    -1,    -1,
     110,   111,   112,    -1,    -1,   189,    -1,   191,    -1,    -1,
     194,   121,    -1,    -1,    -1,   125,    -1,    -1,    -1,   203,
      -1,    -1,    -1,   207,    -1,    -1,    -1,    -1,    -1,    -1,
     214,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   233,
      -1,    -1,    -1,   163,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,   173,   174,    -1,    -1,   177,    -1,   253,
     254,    -1,    -1,     5,    -1,   259,    -1,    -1,    -1,   189,
      -1,   191,   266,    -1,   194,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,   203,    -1,    -1,    -1,   207,    30,    31,
      32,    33,    34,    -1,   214,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    46,    47,    48,    49,    50,    51,
      52,    -1,    -1,   233,    -1,    -1,    -1,    59,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    73,    -1,   253,   254,    -1,    -1,    -1,    -1,    81,
      82,    -1,    -1,    -1,    -1,    -1,   266,    -1,    -1,    -1,
      92,    -1,    -1,    -1,    96,    -1,    -1,    99,    -1,    -1,
      -1,   103,   104,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,   116,   117,   118,   119,   120,    -1,
      -1,    -1,    -1,    -1,    -1,   127,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,   143,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,   155,   156,   157,   158,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,   166,   167,    -1,    -1,    -1,   171,
     172,    -1,    -1,   175,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,   197,    -1,    -1,    -1,    -1,
     202,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   220,   221,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
     232,   233,    -1,    -1,   236,    -1,    -1,   239,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,   247
};

  /* YYSTOS[STATE-NUM] -- The (internal number of the) accessing
     symbol of state STATE-NUM.  */
static const yytype_uint16 yystos[] =
{
       0,     5,    30,    31,    32,    33,    34,    46,    47,    48,
      49,    50,    51,    52,    59,    73,    81,    82,    92,    96,
      99,   103,   104,   116,   117,   118,   119,   120,   127,   143,
     155,   156,   157,   158,   166,   167,   171,   172,   175,   197,
     202,   220,   221,   232,   233,   236,   239,   247,   272,   273,
     274,   276,   290,   291,   292,   306,   307,   308,   310,   315,
     316,   317,   318,   319,   320,   325,   329,   332,   333,   334,
     335,   336,   337,   338,   350,   351,   352,   367,   390,   395,
     397,   398,   399,   405,   410,   411,   415,   429,   430,   460,
     465,   468,   476,   486,   490,   525,   550,   250,   250,   250,
     250,   250,   250,   276,   340,   341,   264,   264,   276,   391,
     276,   342,   354,   264,   264,   250,   266,   266,   266,   250,
     396,   250,   250,   250,   250,   250,   250,   250,   250,   264,
     264,   264,   400,   266,   266,   250,   250,   264,   250,   266,
     250,   250,   250,   264,   264,   264,     0,   274,   250,    79,
     250,    25,    43,    62,   147,   178,   185,   199,   201,   244,
     417,   164,    72,   138,   148,   149,   229,   249,   282,   282,
      80,    95,   208,   275,   276,   277,   303,   304,   305,     3,
       6,    18,    20,    28,    40,    70,    78,    97,   110,   111,
     112,   121,   125,   163,   173,   174,   177,   189,   191,   194,
     203,   207,   214,   253,   254,   266,   276,   297,   298,   299,
     301,   302,   277,   297,   331,   264,   339,   340,   276,   347,
     349,   368,   369,   378,   379,   380,   264,   264,   353,   354,
     276,   278,   477,   478,   297,   276,   311,   311,   311,   303,
     276,   399,   405,   297,   297,   282,   297,   297,   297,   297,
     213,   242,   282,   278,   279,   487,   488,    14,    26,    35,
      58,    76,    77,   102,   129,   154,   168,   169,   170,   176,
     181,   217,   234,   237,   240,   321,   322,    39,   303,   312,
     312,    98,   297,   297,   151,   492,   297,   276,   292,   297,
     297,   297,   526,    75,   551,    16,    37,    53,    88,    89,
     106,   107,   122,   123,   124,   134,   136,   137,   152,   153,
     231,   326,   327,   259,   276,   293,   295,   297,   299,   300,
     302,   303,   304,   266,   276,   294,   295,   296,   400,   400,
     400,   400,   400,   400,   266,   251,   266,   266,   266,   266,
     266,   266,   266,   266,   266,   266,   266,   266,   266,   266,
     266,   266,   266,   266,   297,   297,   297,    69,   139,   253,
     254,   255,   256,   258,   330,   252,   330,    56,    57,   343,
     265,   340,   261,   262,   263,   264,   283,   284,   285,   286,
     287,   265,   369,   253,   270,   381,   393,   355,   265,   354,
     479,   265,   478,   267,   263,   263,   259,   264,   265,   488,
     250,   250,   250,   250,   250,   250,   250,   250,   250,   250,
     250,   250,   250,   250,   250,   250,   250,   250,   265,   322,
     250,   313,   267,   250,   491,   263,   126,   528,   529,   250,
     133,   552,   553,   250,   250,   250,   250,   250,   250,   250,
     250,   250,   250,   250,   250,   250,   250,   250,   250,   265,
     327,   259,   288,   289,   297,   277,    39,   276,   398,   399,
     406,   407,   408,   412,   235,   432,   416,   432,   312,   303,
     297,   297,   297,   297,   297,   297,   297,   297,   297,   297,
     297,   297,   297,   297,   297,   297,   303,   297,   267,   297,
     297,   297,   297,   297,   297,   250,   250,    44,    45,   344,
     297,   285,   286,   380,   276,   356,   382,   254,   269,   370,
     371,   372,   373,   374,   375,   376,    54,   145,   182,   224,
     225,   265,   357,   259,   312,   312,   276,   489,   282,   323,
     323,   323,    27,   282,   324,   323,   323,   297,   323,   323,
     323,   323,   297,   323,   324,   323,   323,   323,   324,   324,
     280,   294,   263,   267,   297,   100,   206,   218,   493,   494,
     495,   496,   312,   250,    74,   133,   527,   530,   531,   532,
     533,   303,   250,   109,   553,   556,    64,    90,   246,   328,
     328,   328,   297,   328,   297,   328,   328,   297,   328,   328,
     328,   328,   328,   328,   328,   297,   260,   263,   263,   250,
     147,   402,   407,   192,   413,   264,   216,   462,    38,    55,
      86,   109,   128,   132,   146,   183,   184,   195,   196,   205,
     241,   418,   420,   421,   422,   423,   424,   425,   426,    60,
     435,   313,   267,   267,   267,   267,   267,   267,   267,   267,
     267,   267,   263,   263,   263,   263,   267,   267,   267,   267,
     252,   297,   297,   250,   250,   215,   345,   265,   283,   370,
     254,   268,   349,   268,   349,   141,   349,   383,   384,   250,
     250,   250,   250,   250,     4,    36,   130,   131,   180,   227,
     265,   358,   359,   360,   362,   363,   276,   482,   313,   267,
     260,   483,   263,   293,   300,   304,   314,   250,   250,   250,
      87,   193,   264,   497,   498,   499,   502,   503,   313,    17,
      29,   140,   250,   264,   265,   530,   276,   554,   555,   250,
     245,   557,   222,   289,   303,   309,   280,   278,   188,   190,
     226,   265,   401,   403,   404,   250,   418,   280,   433,   434,
     264,   461,   250,   250,   264,   250,   250,   264,   250,   250,
     250,   250,   264,   402,   420,   250,   264,   431,   267,   297,
     297,   297,   297,   297,   297,   297,   113,   346,   370,   371,
     371,   253,   259,   385,   386,   387,   297,   230,   297,   394,
     297,   297,   297,   349,   264,   264,   222,   250,   361,   260,
     267,   212,   265,   363,   484,   485,   280,   294,   297,   294,
     250,   250,   500,   265,   498,   267,   303,   101,   135,   219,
     534,   535,   536,   542,   546,   253,   280,   250,   243,   558,
     297,   267,   263,   400,   250,   250,   250,    42,    61,   108,
     179,   198,   200,   266,   278,   279,   414,   402,   265,   433,
     294,   463,   464,   402,   297,   297,   115,   280,   349,   349,
     427,   428,   297,   276,   392,   297,   280,   114,   401,   276,
     295,   300,   419,   294,   436,   437,   438,   267,   267,   267,
     267,   252,   250,   265,   384,   268,   269,   293,   300,   304,
     388,   389,   252,   377,   250,   276,   348,   364,   366,   365,
     366,    13,   349,   480,   250,   162,   474,   140,   282,   303,
     501,   282,    41,    71,   211,   228,   254,   266,   298,   504,
     505,   506,   507,   264,   264,   264,   265,   535,   554,   276,
     280,   281,   300,   250,   100,   206,   218,   559,   206,   280,
     409,   280,   281,   280,   414,   253,   254,   255,   401,   265,
     463,   401,   250,   280,   265,   428,   250,   265,   437,    52,
     186,   402,   439,   440,   455,   297,   297,   389,   389,   260,
     263,   276,   297,   283,   250,   265,   366,   265,   366,   283,
     264,   356,   250,   466,   508,   259,   266,   509,   506,   506,
     263,   265,   252,   253,   254,   255,   256,   524,     8,   150,
     165,   547,   548,   549,   538,   543,   544,   549,   281,   250,
     250,   250,   265,   297,    19,   475,   402,   263,   267,   414,
     414,   414,   297,   297,   264,   264,   265,   388,   297,    66,
      67,    68,    91,    93,    94,   444,   445,   446,   447,   451,
     452,   453,   454,   282,   159,   471,   259,   297,   506,   259,
     267,   505,   250,   253,   268,   510,   303,   506,   506,   506,
     506,   265,   548,   270,    13,   265,   275,   276,   539,   540,
     265,   544,   270,   294,   297,   294,   260,   250,   469,   401,
     297,   204,   265,   204,   456,   457,   482,   441,   260,   250,
     250,   445,   481,   250,   250,   250,   160,   472,   275,   276,
     512,   513,   514,   515,   516,   517,   518,   519,   520,   260,
     267,   514,   268,   268,   268,   277,   511,    12,   294,   545,
     250,   537,   537,    15,   294,   541,   297,   438,   250,   250,
     265,   457,   264,     9,    21,    24,    83,   105,   187,   223,
     297,   442,   443,   444,   448,   449,   450,   259,   259,   483,
     276,   280,   282,   250,   161,   473,   263,   287,   263,   263,
     260,   260,   268,   297,   470,   297,   297,   458,   222,   265,
     263,   449,   449,   265,   263,   282,   250,   467,   248,   278,
     279,   521,   521,   521,   402,   265,    38,   444,   297,   450,
     260,   260,   280,   282,   401,   263,   522,   263,   522,   263,
     522,   401,   250,   459,     7,    10,    11,    22,    23,    65,
      84,    85,   280,   523,   280,   280,   297,   483,   263,   522,
     522,   522,   265,   265
};

  /* YYR1[YYN] -- Symbol number of symbol that rule YYN derives.  */
static const yytype_uint16 yyr1[] =
{
       0,   271,   272,   273,   273,   274,   274,   274,   274,   274,
     274,   274,   274,   274,   274,   274,   274,   274,   274,   274,
     274,   274,   274,   274,   274,   274,   275,   276,   277,   278,
     279,   280,   281,   281,   282,   282,   282,   282,   282,   282,
     283,   283,   283,   283,   284,   284,   284,   284,   285,   286,
     287,   288,   288,   289,   289,   290,   291,   291,   291,   291,
     292,   293,   294,   294,   295,   296,   297,   297,   298,   298,
     299,   299,   300,   300,   301,   302,   302,   302,   302,   302,
     302,   302,   302,   302,   302,   302,   302,   302,   302,   302,
     302,   302,   302,   302,   302,   302,   302,   302,   302,   302,
     302,   302,   302,   302,   302,   303,   303,   304,   304,   304,
     304,   305,   306,   306,   306,   306,   306,   306,   306,   307,
     308,   309,   310,   311,   312,   313,   313,   314,   314,   314,
     315,   316,   317,   318,   319,   320,   321,   321,   322,   322,
     322,   322,   322,   322,   322,   322,   322,   322,   322,   322,
     322,   322,   322,   322,   322,   322,   322,   323,   324,   324,
     325,   326,   326,   327,   327,   327,   327,   327,   327,   327,
     327,   327,   327,   327,   327,   327,   327,   327,   327,   328,
     328,   328,   329,   329,   329,   329,   329,   330,   330,   330,
     331,   331,   331,   331,   332,   332,   332,   332,   332,   332,
     332,   332,   332,   332,   332,   332,   332,   332,   332,   333,
     333,   333,   333,   334,   335,   335,   335,   336,   336,   337,
     338,   339,   339,   340,   341,   342,   343,   343,   344,   344,
     344,   345,   345,   346,   346,   347,   348,   349,   350,   350,
     351,   352,   353,   353,   355,   354,   356,   357,   357,   358,
     358,   359,   359,   359,   360,   360,   360,   361,   361,   362,
     363,   363,   364,   364,   365,   365,   366,   367,   368,   368,
     369,   370,   370,   371,   372,   373,   374,   375,   376,   376,
     376,   376,   377,   377,   378,   379,   379,   380,   381,   381,
     382,   383,   383,   384,   384,   385,   385,   386,   387,   388,
     388,   388,   389,   389,   389,   390,   391,   392,   393,   393,
     393,   393,   393,   393,   394,   394,   396,   395,   397,   398,
     398,   398,   398,   398,   399,   400,   401,   402,   402,   403,
     404,   404,   404,   405,   406,   406,   407,   407,   409,   408,
     410,   410,   412,   411,   413,   413,   413,   413,   413,   413,
     413,   413,   414,   414,   414,   414,   414,   416,   415,   417,
     417,   417,   417,   417,   418,   418,   419,   420,   420,   420,
     420,   420,   420,   420,   420,   420,   420,   421,   421,   422,
     422,   422,   422,   423,   423,   424,   425,   426,   426,   427,
     427,   428,   429,   431,   430,   432,   433,   434,   434,   435,
     436,   436,   437,   438,   438,   439,   439,   441,   440,   442,
     442,   443,   443,   443,   443,   443,   443,   443,   444,   444,
     445,   445,   445,   445,   446,   447,   448,   449,   449,   450,
     450,   450,   451,   452,   452,   453,   454,   454,   455,   456,
     456,   458,   459,   457,   461,   460,   462,   463,   464,   464,
     466,   467,   465,   469,   470,   468,   471,   471,   472,   472,
     473,   473,   474,   474,   475,   475,   476,   477,   477,   479,
     480,   481,   478,   482,   483,   483,   484,   484,   485,   486,
     487,   487,   489,   488,   491,   490,   492,   492,   493,   493,
     493,   494,   495,   496,   497,   497,   498,   498,   498,   500,
     499,   501,   501,   501,   502,   503,   504,   504,   505,   506,
     506,   506,   506,   506,   506,   506,   506,   506,   508,   507,
     507,   509,   507,   510,   510,   510,   510,   510,   511,   512,
     513,   514,   514,   514,   514,   514,   514,   515,   516,   517,
     518,   519,   520,   521,   521,   521,   522,   522,   523,   523,
     523,   523,   523,   523,   523,   523,   524,   524,   526,   525,
     527,   527,   528,   528,   529,   529,   529,   530,   530,   531,
     532,   533,   534,   534,   535,   535,   535,   536,   537,   537,
     538,   538,   539,   539,   540,   540,   541,   541,   542,   543,
     543,   544,   545,   545,   546,   547,   547,   548,   549,   549,
     549,   550,   551,   552,   552,   553,   554,   555,   555,   556,
     557,   558,   559,   559,   559,   559
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
       3,     3,     3,     3,     3,     3,     3,     3,     3,     1,
       1,     1,     3,     3,     4,     3,     4,     0,     1,     1,
       1,     3,     5,     7,     3,     3,     3,     3,     3,     3,
       3,     3,     3,     3,     3,     3,     3,     3,     3,     3,
       3,     3,     3,     3,     1,     1,     1,     1,     1,     2,
       4,     1,     2,     7,     1,     1,     3,     3,     0,     3,
       3,     0,     1,     0,     3,     1,     2,     2,     1,     1,
       2,     4,     1,     2,     0,     5,     1,     0,     2,     1,
       1,     3,     4,     4,     1,     1,     1,     1,     1,     1,
       4,     4,     1,     2,     1,     2,     3,     4,     1,     2,
       1,     1,     2,     2,     2,     2,     3,     3,     1,     1,
       1,     1,     0,     2,     6,     1,     3,     1,     0,     2,
       2,     1,     3,     1,     1,     1,     1,     3,     5,     1,
       2,     2,     1,     1,     1,     5,     1,     1,     0,     4,
       4,     4,     4,     4,     1,     1,     0,     3,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     0,     2,     1,
       3,     3,     5,     6,     1,     2,     1,     1,     0,     7,
       1,     1,     0,     8,     3,     3,     3,     3,     3,     3,
       3,     3,     1,     3,     3,     3,     3,     0,     7,     1,
       1,     1,     1,     1,     1,     2,     1,     3,     3,     1,
       3,     3,     3,     3,     3,     3,     4,     1,     1,     1,
       1,     1,     1,     3,     6,     9,    12,     3,     3,     1,
       2,     2,     1,     0,     9,     4,     1,     1,     2,     4,
       1,     2,     1,     0,     2,     1,     1,     0,     5,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     2,
       1,     1,     1,     1,     5,     5,     1,     1,     3,     1,
       3,     1,     3,     1,     1,     5,     1,     1,     4,     1,
       2,     0,     0,     7,     0,     8,     4,     1,     1,     2,
       0,     0,    14,     0,     0,    14,     0,     3,     0,     3,
       0,     3,     0,     3,     0,     3,     4,     1,     2,     0,
       0,     0,    11,     1,     0,     2,     1,     1,     3,     4,
       1,     2,     0,     5,     0,     7,     0,     3,     1,     1,
       1,     3,     3,     3,     1,     2,     1,     1,     1,     0,
       6,     1,     1,     1,     3,     3,     1,     3,     2,     1,
       1,     3,     3,     3,     3,     3,     2,     4,     0,     5,
       4,     0,     5,     1,     2,     2,     3,     2,     1,     1,
       2,     1,     1,     1,     1,     1,     1,     4,     4,     4,
       6,     6,     6,     1,     1,     1,     0,     2,     1,     1,
       1,     1,     1,     1,     1,     1,     0,     2,     0,     6,
       1,     2,     0,     1,     3,     3,     3,     1,     1,     1,
       3,     4,     1,     2,     1,     1,     1,     4,     2,     0,
       2,     0,     2,     2,     1,     1,     1,     1,     4,     1,
       2,     3,     1,     1,     4,     1,     2,     3,     1,     1,
       1,     9,     3,     1,     2,     3,     1,     1,     3,     3,
       3,     3,     0,     3,     3,     3
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
#line 637 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.sym) = mdl_existing_object(parse_state, (yyvsp[0].str))); }
#line 3141 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 30:
#line 640 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.sym) = mdl_existing_region(parse_state, (yyvsp[-3].sym), (yyvsp[-1].str))); }
#line 3147 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 31:
#line 643 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.vec3) = mdl_point(parse_state, &(yyvsp[0].nlist))); }
#line 3153 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 33:
#line 647 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.vec3) = mdl_point_scalar((yyvsp[0].dbl))); }
#line 3159 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 34:
#line 650 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = 1; }
#line 3165 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 35:
#line 651 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = 0; }
#line 3171 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 36:
#line 652 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = 1; }
#line 3177 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 37:
#line 653 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = 0; }
#line 3183 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 38:
#line 654 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = 1; }
#line 3189 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 39:
#line 655 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = 0; }
#line 3195 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 40:
#line 658 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.mol_type).orient_set = 0; }
#line 3201 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 43:
#line 661 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.mol_type).orient_set = 1; (yyval.mol_type).orient = 0; }
#line 3207 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 44:
#line 665 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.mol_type).orient = 1; (yyval.mol_type).orient_set = 1; }
#line 3213 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 45:
#line 666 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.mol_type).orient = -1; (yyval.mol_type).orient_set = 1; }
#line 3219 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 46:
#line 667 "../src/../src/mdlparse.y" /* yacc.c:1646  */
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
#line 3234 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 47:
#line 677 "../src/../src/mdlparse.y" /* yacc.c:1646  */
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
#line 3249 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 50:
#line 695 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                          (yyval.mol_type).orient = (int) (yyvsp[-1].dbl);
                                                          (yyval.mol_type).orient_set = 1;
                                                          if ((yyval.mol_type).orient != (yyvsp[-1].dbl))
                                                          {
                                                            mdlerror(parse_state, "molecule orientation specified inside braces must be an integer between -32768 and 32767.");
                                                            return 1;
                                                          }
                                                      }
#line 3263 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 52:
#line 708 "../src/../src/mdlparse.y" /* yacc.c:1646  */
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
#line 3279 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 53:
#line 721 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mcell_generate_range_singleton(&(yyval.nlist), (yyvsp[0].dbl))); }
#line 3285 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 54:
#line 722 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_generate_range(parse_state, &(yyval.nlist), (yyvsp[-5].dbl), (yyvsp[-3].dbl), (yyvsp[-1].dbl))); }
#line 3291 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 55:
#line 728 "../src/../src/mdlparse.y" /* yacc.c:1646  */
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
#line 3313 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 56:
#line 751 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_assign_variable_double(parse_state, (yyvsp[-2].sym), (yyvsp[0].dbl))); }
#line 3319 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 57:
#line 752 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_assign_variable_string(parse_state, (yyvsp[-2].sym), (yyvsp[0].str))); }
#line 3325 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 58:
#line 753 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_assign_variable(parse_state, (yyvsp[-2].sym), (yyvsp[0].sym))); }
#line 3331 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 59:
#line 754 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_assign_variable_array(parse_state, (yyvsp[-2].sym), (yyvsp[0].nlist).value_head)); }
#line 3337 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 60:
#line 757 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.sym) = mdl_get_or_create_variable(parse_state, (yyvsp[0].str))); }
#line 3343 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 61:
#line 760 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.sym) = mdl_existing_variable(parse_state, (yyvsp[0].str))); }
#line 3349 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 63:
#line 764 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                          struct num_expr_list *elp;
                                                          (yyval.nlist).value_head = (struct num_expr_list *) (yyvsp[0].sym)->value;
                                                          (yyval.nlist).value_count = 1;
                                                          for (elp = (yyval.nlist).value_head; elp->next != NULL; elp = elp->next)
                                                            ++ (yyval.nlist).value_count;
                                                          (yyval.nlist).value_tail = elp;
                                                          (yyval.nlist).shared = 1;
                                                      }
#line 3363 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 64:
#line 775 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { mdl_debug_dump_array((yyvsp[-1].nlist).value_head); (yyval.nlist) = (yyvsp[-1].nlist); }
#line 3369 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 65:
#line 778 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.sym) = mdl_existing_array(parse_state, (yyvsp[0].str))); }
#line 3375 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 69:
#line 786 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.dbl) = *(double *) (yyvsp[0].sym)->value; }
#line 3381 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 70:
#line 789 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.dbl) = (yyvsp[0].llival); }
#line 3387 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 74:
#line 797 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.sym) = mdl_existing_double(parse_state, (yyvsp[0].str))); }
#line 3393 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 75:
#line 801 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.dbl) = (yyvsp[-1].dbl); }
#line 3399 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 76:
#line 802 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKF((yyval.dbl) = exp((yyvsp[-1].dbl))); }
#line 3405 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 77:
#line 803 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_expr_log(parse_state, (yyvsp[-1].dbl), &(yyval.dbl))); }
#line 3411 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 78:
#line 804 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_expr_log10(parse_state, (yyvsp[-1].dbl), &(yyval.dbl))); }
#line 3417 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 79:
#line 805 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.dbl) = max2d((yyvsp[-3].dbl), (yyvsp[-1].dbl)); }
#line 3423 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 80:
#line 806 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.dbl) = min2d((yyvsp[-3].dbl), (yyvsp[-1].dbl)); }
#line 3429 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 81:
#line 807 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.dbl) = mdl_expr_roundoff((yyvsp[-1].dbl), (int) (yyvsp[-3].dbl)); }
#line 3435 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 82:
#line 808 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.dbl) = floor((yyvsp[-1].dbl)); }
#line 3441 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 83:
#line 809 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.dbl) = ceil((yyvsp[-1].dbl)); }
#line 3447 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 84:
#line 810 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.dbl) = sin((yyvsp[-1].dbl)); }
#line 3453 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 85:
#line 811 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.dbl) = cos((yyvsp[-1].dbl)); }
#line 3459 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 86:
#line 812 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKF((yyval.dbl) = tan((yyvsp[-1].dbl))); }
#line 3465 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 87:
#line 813 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKF((yyval.dbl) = asin((yyvsp[-1].dbl))); }
#line 3471 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 88:
#line 814 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKF((yyval.dbl) = acos((yyvsp[-1].dbl))); }
#line 3477 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 89:
#line 815 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.dbl) = atan((yyvsp[-1].dbl)); }
#line 3483 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 90:
#line 816 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKF((yyval.dbl) = sqrt((yyvsp[-1].dbl))); }
#line 3489 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 91:
#line 817 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.dbl) = fabs((yyvsp[-1].dbl)); }
#line 3495 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 92:
#line 818 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_expr_mod(parse_state, (yyvsp[-3].dbl), (yyvsp[-1].dbl), &(yyval.dbl))); }
#line 3501 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 93:
#line 819 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.dbl) = MY_PI; }
#line 3507 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 94:
#line 820 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.dbl) = mdl_expr_rng_uniform(parse_state); }
#line 3513 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 95:
#line 821 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.dbl) = rng_gauss(parse_state->vol->rng); }
#line 3519 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 96:
#line 822 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.dbl) = parse_state->vol->seed_seq; }
#line 3525 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 97:
#line 823 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_expr_string_to_double(parse_state, (yyvsp[-1].str), &(yyval.dbl))); }
#line 3531 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 98:
#line 824 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKF((yyval.dbl) = (yyvsp[-2].dbl) + (yyvsp[0].dbl)); }
#line 3537 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 99:
#line 825 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKF((yyval.dbl) = (yyvsp[-2].dbl) - (yyvsp[0].dbl)); }
#line 3543 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 100:
#line 826 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKF((yyval.dbl) = (yyvsp[-2].dbl) * (yyvsp[0].dbl)); }
#line 3549 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 101:
#line 827 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_expr_div(parse_state, (yyvsp[-2].dbl), (yyvsp[0].dbl), &(yyval.dbl))); }
#line 3555 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 102:
#line 828 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_expr_pow(parse_state, (yyvsp[-2].dbl), (yyvsp[0].dbl), &(yyval.dbl))); }
#line 3561 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 103:
#line 829 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.dbl) = -(yyvsp[0].dbl); }
#line 3567 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 104:
#line 830 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.dbl) = (yyvsp[0].dbl); }
#line 3573 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 106:
#line 835 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.str) = mdl_strdup((char const *) (yyvsp[0].sym)->value)); }
#line 3579 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 107:
#line 839 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.str) = mdl_strip_quotes((yyvsp[0].str))); }
#line 3585 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 108:
#line 840 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.str) = mdl_strdup(parse_state->vol->mdl_infile_name)); }
#line 3591 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 109:
#line 841 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.str) = mdl_strcat((yyvsp[-2].str), (yyvsp[0].str))); }
#line 3597 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 110:
#line 842 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.str) = mdl_string_format(parse_state, (yyvsp[-2].str), (yyvsp[-1].printfargs).arg_head)); }
#line 3603 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 111:
#line 845 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.sym) = mdl_existing_string(parse_state, (yyvsp[0].str))); }
#line 3609 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 119:
#line 861 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_fopen(parse_state, (yyvsp[-6].sym), (yyvsp[-3].str), (yyvsp[-1].str))); }
#line 3615 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 120:
#line 864 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.sym) = mdl_new_filehandle(parse_state, (yyvsp[0].str))); }
#line 3621 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 121:
#line 867 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.str) = (yyvsp[0].str); CHECK(mdl_valid_file_mode(parse_state, (yyvsp[0].str))); }
#line 3627 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 122:
#line 870 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_fclose(parse_state, (yyvsp[-1].sym))); }
#line 3633 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 123:
#line 873 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.sym) = mdl_existing_file_stream(parse_state, (yyvsp[0].str))); }
#line 3639 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 124:
#line 876 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.str) = mdl_expand_string_escapes((yyvsp[0].str))); }
#line 3645 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 125:
#line 879 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.printfargs).arg_head = (yyval.printfargs).arg_tail = NULL; }
#line 3651 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 126:
#line 880 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                        (yyval.printfargs) = (yyvsp[-2].printfargs);
                                                        if ((yyval.printfargs).arg_tail)
                                                          (yyval.printfargs).arg_tail = (yyval.printfargs).arg_tail->next = (yyvsp[0].printfarg);
                                                        else
                                                          (yyval.printfargs).arg_tail = (yyval.printfargs).arg_head = (yyvsp[0].printfarg);
                                                        (yyvsp[0].printfarg)->next = NULL;
                                                      }
#line 3664 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 127:
#line 890 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.printfarg) = mdl_new_printf_arg_double((yyvsp[0].dbl))); }
#line 3670 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 128:
#line 891 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.printfarg) = mdl_new_printf_arg_string((yyvsp[0].str))); }
#line 3676 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 129:
#line 892 "../src/../src/mdlparse.y" /* yacc.c:1646  */
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
#line 3691 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 130:
#line 904 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_printf(parse_state, (yyvsp[-2].str), (yyvsp[-1].printfargs).arg_head)); }
#line 3697 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 131:
#line 910 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_fprintf(parse_state, (struct file_stream *) (yyvsp[-4].sym)->value, (yyvsp[-2].str), (yyvsp[-1].printfargs).arg_head)); }
#line 3703 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 132:
#line 916 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_sprintf(parse_state, (yyvsp[-4].sym), (yyvsp[-2].str), (yyvsp[-1].printfargs).arg_head)); }
#line 3709 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 133:
#line 919 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { mdl_print_time(parse_state, (yyvsp[-1].str)); }
#line 3715 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 134:
#line 925 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_fprint_time(parse_state, (yyvsp[-3].sym), (yyvsp[-1].str))); }
#line 3721 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 138:
#line 941 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { if (!parse_state->vol->quiet_flag) mdl_set_all_notifications(parse_state->vol, (yyvsp[0].tok)); }
#line 3727 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 139:
#line 942 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { if (!parse_state->vol->quiet_flag) parse_state->vol->notify->progress_report        = (yyvsp[0].tok); }
#line 3733 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 140:
#line 943 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { if (!parse_state->vol->quiet_flag) parse_state->vol->notify->diffusion_constants    = (yyvsp[0].tok); }
#line 3739 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 141:
#line 944 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { if (!parse_state->vol->quiet_flag) parse_state->vol->notify->reaction_probabilities = (yyvsp[0].tok); }
#line 3745 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 142:
#line 945 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { if (!parse_state->vol->quiet_flag) parse_state->vol->notify->time_varying_reactions = (yyvsp[0].tok); }
#line 3751 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 143:
#line 946 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { if (!parse_state->vol->quiet_flag) parse_state->vol->notify->reaction_prob_notify   = (yyvsp[0].dbl); }
#line 3757 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 144:
#line 947 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { if (!parse_state->vol->quiet_flag) parse_state->vol->notify->partition_location     = (yyvsp[0].tok); }
#line 3763 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 145:
#line 948 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { if (!parse_state->vol->quiet_flag) parse_state->vol->notify->box_triangulation      = (yyvsp[0].tok); }
#line 3769 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 146:
#line 949 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { if (!parse_state->vol->quiet_flag) parse_state->vol->notify->release_events         = (yyvsp[0].tok); }
#line 3775 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 147:
#line 950 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { if (!parse_state->vol->quiet_flag) parse_state->vol->notify->file_writes            = (yyvsp[0].tok); }
#line 3781 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 148:
#line 951 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { if (!parse_state->vol->quiet_flag) parse_state->vol->notify->final_summary          = (yyvsp[0].tok); }
#line 3787 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 149:
#line 952 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { if (!parse_state->vol->quiet_flag) parse_state->vol->notify->throughput_report      = (yyvsp[0].tok); }
#line 3793 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 150:
#line 953 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { if (!parse_state->vol->quiet_flag) parse_state->vol->notify->reaction_output_report = (yyvsp[0].tok); }
#line 3799 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 151:
#line 954 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { if (!parse_state->vol->quiet_flag) parse_state->vol->notify->volume_output_report   = (yyvsp[0].tok); }
#line 3805 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 152:
#line 955 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { if (!parse_state->vol->quiet_flag) parse_state->vol->notify->viz_output_report      = (yyvsp[0].tok); }
#line 3811 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 153:
#line 956 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { if (!parse_state->vol->quiet_flag) parse_state->vol->notify->checkpoint_report      = (yyvsp[0].tok); }
#line 3817 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 154:
#line 957 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                          if (!parse_state->vol->quiet_flag && parse_state->vol->log_freq == ULONG_MAX)
                                                            parse_state->vol->notify->iteration_report = (yyvsp[0].tok);
                                                      }
#line 3826 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 155:
#line 961 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { if (!parse_state->vol->quiet_flag) CHECK(mdl_set_iteration_report_freq(parse_state, (long long) (yyvsp[0].dbl))); }
#line 3832 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 156:
#line 962 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { if (!parse_state->vol->quiet_flag) parse_state->vol->notify->molecule_collision_report    = (yyvsp[0].tok); }
#line 3838 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 157:
#line 966 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = ((yyvsp[0].tok) ? NOTIFY_FULL : NOTIFY_NONE); }
#line 3844 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 158:
#line 970 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = ((yyvsp[0].tok) ? NOTIFY_FULL : NOTIFY_NONE); }
#line 3850 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 159:
#line 971 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = NOTIFY_BRIEF; }
#line 3856 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 163:
#line 987 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { mdl_set_all_warnings(parse_state->vol, (byte) (yyvsp[0].tok)); }
#line 3862 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 164:
#line 988 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->vol->notify->neg_diffusion = (byte)(yyvsp[0].tok); }
#line 3868 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 165:
#line 989 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->vol->notify->neg_reaction = (byte)(yyvsp[0].tok); }
#line 3874 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 166:
#line 990 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->vol->notify->high_reaction_prob = (byte)(yyvsp[0].tok); }
#line 3880 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 167:
#line 991 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->vol->notify->reaction_prob_warn = (yyvsp[0].dbl); }
#line 3886 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 168:
#line 992 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->vol->notify->close_partitions = (byte)(yyvsp[0].tok); }
#line 3892 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 169:
#line 993 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->vol->notify->degenerate_polys = (byte)(yyvsp[0].tok); }
#line 3898 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 170:
#line 994 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->vol->notify->overwritten_file = (byte)(yyvsp[0].tok); }
#line 3904 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 171:
#line 995 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->vol->notify->short_lifetime = (byte)(yyvsp[0].tok); }
#line 3910 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 172:
#line 996 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_lifetime_warning_threshold(parse_state, (long long) (yyvsp[0].dbl))); }
#line 3916 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 173:
#line 997 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->vol->notify->missed_reactions = (byte)(yyvsp[0].tok); }
#line 3922 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 174:
#line 998 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_missed_reaction_warning_threshold(parse_state, (yyvsp[0].dbl))); }
#line 3928 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 175:
#line 999 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->vol->notify->missed_surf_orient = (byte)(yyvsp[0].tok); }
#line 3934 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 176:
#line 1000 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->vol->notify->useless_vol_orient = (byte)(yyvsp[0].tok); }
#line 3940 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 177:
#line 1001 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->vol->notify->mol_placement_failure = (byte) (yyvsp[0].tok); }
#line 3946 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 178:
#line 1002 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->vol->notify->invalid_output_step_time = (byte) (yyvsp[0].tok); }
#line 3952 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 179:
#line 1006 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = WARN_COPE;  }
#line 3958 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 180:
#line 1007 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = WARN_WARN;  }
#line 3964 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 181:
#line 1008 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = WARN_ERROR; }
#line 3970 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 182:
#line 1014 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_checkpoint_infile(parse_state, (yyvsp[0].str))); }
#line 3976 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 183:
#line 1015 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_checkpoint_outfile(parse_state, (yyvsp[0].str))); }
#line 3982 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 184:
#line 1016 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_checkpoint_interval(parse_state, (yyvsp[-1].dbl), (yyvsp[0].tok))); }
#line 3988 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 185:
#line 1017 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_keep_checkpoint_files(parse_state, (yyvsp[0].tok))); }
#line 3994 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 186:
#line 1019 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_realtime_checkpoint(parse_state, (long) (yyvsp[-1].dbl), (yyvsp[0].tok))); }
#line 4000 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 187:
#line 1022 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = 0; }
#line 4006 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 188:
#line 1023 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = 1; }
#line 4012 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 189:
#line 1024 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = 0; }
#line 4018 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 190:
#line 1028 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { /* seconds */     (yyval.dbl) = (yyvsp[0].dbl); }
#line 4024 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 191:
#line 1029 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { /* mm:ss */       (yyval.dbl) = (yyvsp[-2].dbl) * 60 + (yyvsp[0].dbl); }
#line 4030 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 192:
#line 1030 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { /* hh:mm:ss */    (yyval.dbl) = (yyvsp[-4].dbl) * 3600 + (yyvsp[-2].dbl) * 60 + (yyvsp[0].dbl); }
#line 4036 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 193:
#line 1032 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { /* dd:hh:mm:ss */ (yyval.dbl) = (yyvsp[-6].dbl) * 86400 + (yyvsp[-4].dbl) * 3600 + (yyvsp[-2].dbl) * 60 + (yyvsp[0].dbl); }
#line 4042 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 194:
#line 1039 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_time_step(parse_state, (yyvsp[0].dbl))); }
#line 4048 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 195:
#line 1040 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_space_step(parse_state, (yyvsp[0].dbl))); }
#line 4054 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 196:
#line 1041 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_max_time_step(parse_state, (yyvsp[0].dbl))); }
#line 4060 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 197:
#line 1042 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_num_iterations(parse_state, (long long) (yyvsp[0].dbl))); }
#line 4066 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 198:
#line 1043 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->vol->randomize_smol_pos = !((yyvsp[0].tok)); }
#line 4072 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 199:
#line 1044 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->vol->use_expanded_list = (yyvsp[0].tok); }
#line 4078 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 200:
#line 1045 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->vol->vacancy_search_dist2 = max2d((yyvsp[0].dbl), 0.0); }
#line 4084 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 201:
#line 1046 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_num_radial_directions(parse_state, (int) (yyvsp[0].dbl))); }
#line 4090 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 202:
#line 1047 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->vol->fully_random = 1; }
#line 4096 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 203:
#line 1048 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_num_radial_subdivisions(parse_state, (int) (yyvsp[0].dbl))); }
#line 4102 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 204:
#line 1049 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_grid_density(parse_state, (yyvsp[0].dbl))); }
#line 4108 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 205:
#line 1050 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_interaction_radius(parse_state, (yyvsp[0].dbl))); }
#line 4114 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 206:
#line 1051 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->vol->surface_reversibility=(yyvsp[0].tok); parse_state->vol->volume_reversibility=(yyvsp[0].tok); }
#line 4120 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 207:
#line 1052 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->vol->surface_reversibility=1;  parse_state->vol->volume_reversibility=0;  }
#line 4126 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 208:
#line 1053 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->vol->surface_reversibility=0;  parse_state->vol->volume_reversibility=1;  }
#line 4132 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 209:
#line 1060 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->vol->mem_part_x = (int) (yyvsp[0].dbl); }
#line 4138 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 210:
#line 1061 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->vol->mem_part_y = (int) (yyvsp[0].dbl); }
#line 4144 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 211:
#line 1062 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->vol->mem_part_z = (int) (yyvsp[0].dbl); }
#line 4150 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 212:
#line 1063 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->vol->mem_part_pool = (int) (yyvsp[0].dbl); }
#line 4156 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 213:
#line 1067 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mcell_set_partition(parse_state->vol, (yyvsp[-2].tok), & (yyvsp[0].nlist))); }
#line 4162 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 214:
#line 1071 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = X_PARTS; }
#line 4168 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 215:
#line 1072 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = Y_PARTS; }
#line 4174 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 216:
#line 1073 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = Z_PARTS; }
#line 4180 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 219:
#line 1084 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { mdl_print_species_summary(parse_state->vol, (yyvsp[0].mcell_mol_spec)); }
#line 4186 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 220:
#line 1088 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { mdl_print_species_summaries(parse_state->vol, (yyvsp[-1].mcell_species_lst).species_head); }
#line 4192 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 221:
#line 1092 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.mcell_species_lst).species_count = 0; CHECK(mdl_add_to_species_list(&(yyval.mcell_species_lst), (yyvsp[0].mcell_mol_spec))); }
#line 4198 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 222:
#line 1093 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.mcell_species_lst) = (yyvsp[-1].mcell_species_lst); CHECK(mdl_add_to_species_list(&(yyval.mcell_species_lst), (yyvsp[0].mcell_mol_spec))); }
#line 4204 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 223:
#line 1102 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.mcell_mol_spec) = mdl_create_species(parse_state, (yyvsp[-6].str), (yyvsp[-4].diff_const).D, (yyvsp[-4].diff_const).is_2d, (yyvsp[-3].dbl), (yyvsp[-2].ival), (yyvsp[-1].dbl) )); }
#line 4210 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 225:
#line 1108 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.sym) = mdl_new_mol_species(parse_state, (yyvsp[0].str))); }
#line 4216 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 226:
#line 1112 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.diff_const).is_2d = 0; (yyval.diff_const).D = (yyvsp[0].dbl); CHECK(mdl_check_diffusion_constant(parse_state, & (yyval.diff_const).D)); }
#line 4222 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 227:
#line 1113 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.diff_const).is_2d = 1; (yyval.diff_const).D = (yyvsp[0].dbl); CHECK(mdl_check_diffusion_constant(parse_state, & (yyval.diff_const).D)); }
#line 4228 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 228:
#line 1117 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.dbl) = 0.0; }
#line 4234 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 229:
#line 1118 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                          if ((yyvsp[0].dbl) <= 0)
                                                          {
                                                            mdlerror_fmt(parse_state, "Requested custom time step of %.15g; custom time step must be positive.", (yyvsp[0].dbl));
                                                            return 1;
                                                          }

                                                          (yyval.dbl) = (yyvsp[0].dbl);
                                                      }
#line 4248 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 230:
#line 1127 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                          if ((yyvsp[0].dbl) <= 0)
                                                          {
                                                            mdlerror_fmt(parse_state, "Requested custom space step of %.15g; custom space step must be positive.", (yyvsp[0].dbl));
                                                            return 1;
                                                          }

                                                          (yyval.dbl) = -(yyvsp[0].dbl);
                                                      }
#line 4262 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 231:
#line 1138 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.ival) = 0; }
#line 4268 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 232:
#line 1139 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.ival) = 1; }
#line 4274 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 233:
#line 1143 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.dbl) = 0; }
#line 4280 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 234:
#line 1144 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                        if ((yyvsp[0].dbl) <= 0)
                                                        {
                                                          mdlerror_fmt(parse_state, "Requested maximum step length of %.15g; maximum step length must be positive.", (yyvsp[0].dbl));
                                                          return 1;
                                                        }
                                                        (yyval.dbl) = (yyvsp[0].dbl);
                                                      }
#line 4293 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 235:
#line 1154 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.sym) = mdl_existing_molecule(parse_state, (yyvsp[0].str))); }
#line 4299 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 236:
#line 1158 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.mol_type) = (yyvsp[0].mol_type); CHECKN((yyval.mol_type).mol_type = mdl_existing_surface_molecule(parse_state, (yyvsp[-1].str))); }
#line 4305 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 237:
#line 1162 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                        (yyval.mol_type) = (yyvsp[0].mol_type);
                                                        if (! (yyval.mol_type).orient_set)
                                                          (yyval.mol_type).orient = 0;
                                                        (yyval.mol_type).mol_type = (yyvsp[-1].sym);
                                                      }
#line 4316 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 244:
#line 1197 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { mdl_start_surface_class(parse_state, (yyvsp[-1].sym)); }
#line 4322 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 245:
#line 1199 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { mdl_finish_surface_class(parse_state); }
#line 4328 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 246:
#line 1202 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.sym) = mdl_existing_surface_class(parse_state, (yyvsp[0].str))); }
#line 4334 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 251:
#line 1219 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN(mdl_assemble_surface_reaction(parse_state, (yyvsp[-2].tok), parse_state->current_surface_class, (yyvsp[0].mol_type).mol_type, (yyvsp[0].mol_type).orient)); }
#line 4340 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 252:
#line 1222 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
              struct sym_entry *mol_sym = retrieve_sym("ALL_MOLECULES", parse_state->vol->mol_sym_table);
              if(!(yyvsp[0].mol_type).orient_set) (yyvsp[0].mol_type).orient = 0;
              CHECKN(mdl_assemble_surface_reaction(parse_state, (yyvsp[-3].tok), parse_state->current_surface_class, mol_sym, (yyvsp[0].mol_type).orient));}
#line 4349 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 253:
#line 1228 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN(mdl_assemble_concentration_clamp_reaction(parse_state, parse_state->current_surface_class, (yyvsp[-2].mol_type).mol_type, (yyvsp[-2].mol_type).orient, (yyvsp[0].dbl))); }
#line 4355 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 254:
#line 1231 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = RFLCT; }
#line 4361 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 255:
#line 1232 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = TRANSP; }
#line 4367 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 256:
#line 1233 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = SINK; }
#line 4373 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 259:
#line 1240 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->current_surface_class->sm_dat_head = (yyvsp[0].surf_mol_dat_list).sm_head; }
#line 4379 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 260:
#line 1247 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.surf_mol_dat_list) = (yyvsp[-1].surf_mol_dat_list); }
#line 4385 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 261:
#line 1251 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.surf_mol_dat_list) = (yyvsp[-1].surf_mol_dat_list); }
#line 4391 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 262:
#line 1255 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                          (yyvsp[0].surf_mol_dat)->quantity_type = SURFMOLDENS;
                                                          (yyval.surf_mol_dat_list).sm_tail = (yyval.surf_mol_dat_list).sm_head = (yyvsp[0].surf_mol_dat);
                                                      }
#line 4400 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 263:
#line 1260 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                          (yyval.surf_mol_dat_list) = (yyvsp[-1].surf_mol_dat_list);
                                                          (yyvsp[0].surf_mol_dat)->quantity_type = SURFMOLDENS;
                                                          (yyval.surf_mol_dat_list).sm_tail = (yyval.surf_mol_dat_list).sm_tail->next = (yyvsp[0].surf_mol_dat);
                                                      }
#line 4410 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 264:
#line 1268 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                          (yyvsp[0].surf_mol_dat)->quantity_type = SURFMOLNUM;
                                                          (yyval.surf_mol_dat_list).sm_tail = (yyval.surf_mol_dat_list).sm_head = (yyvsp[0].surf_mol_dat);
                                                      }
#line 4419 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 265:
#line 1273 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                          (yyval.surf_mol_dat_list) = (yyvsp[-1].surf_mol_dat_list);
                                                          (yyvsp[0].surf_mol_dat)->quantity_type = SURFMOLNUM;
                                                          (yyval.surf_mol_dat_list).sm_tail = (yyval.surf_mol_dat_list).sm_tail->next = (yyvsp[0].surf_mol_dat);
                                                      }
#line 4429 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 266:
#line 1281 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.surf_mol_dat) = mdl_new_surf_mol_data(parse_state, &(yyvsp[-2].mol_type), (yyvsp[0].dbl))); }
#line 4435 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 276:
#line 1310 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.react_arrow).catalyst = (yyvsp[-1].mol_type); (yyval.react_arrow).flags = ARROW_CATALYTIC; }
#line 4441 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 277:
#line 1315 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.react_arrow).catalyst = (yyvsp[-1].mol_type); (yyval.react_arrow).flags = ARROW_CATALYTIC | ARROW_BIDIRECTIONAL; }
#line 4447 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 278:
#line 1320 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.react_arrow).catalyst.mol_type = NULL; (yyval.react_arrow).flags = 0; }
#line 4453 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 280:
#line 1322 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.react_arrow).catalyst.mol_type = NULL; (yyval.react_arrow).flags = ARROW_BIDIRECTIONAL; }
#line 4459 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 282:
#line 1326 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.sym) = NULL; }
#line 4465 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 283:
#line 1327 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.sym) = mdl_new_rxn_pathname(parse_state, (yyvsp[0].str))); }
#line 4471 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 284:
#line 1333 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN(mdl_assemble_reaction(parse_state, (yyvsp[-5].mol_type_list).mol_type_head, &(yyvsp[-4].mol_type), &(yyvsp[-3].react_arrow), (yyvsp[-2].mol_type_list).mol_type_head, &(yyvsp[-1].react_rates), (yyvsp[0].sym))); }
#line 4477 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 285:
#line 1336 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_reaction_player_singleton(parse_state, & (yyval.mol_type_list), & (yyvsp[0].mol_type))); }
#line 4483 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 286:
#line 1337 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.mol_type_list) = (yyvsp[-2].mol_type_list); CHECK(mdl_add_reaction_player(parse_state, & (yyval.mol_type_list), & (yyvsp[0].mol_type))); }
#line 4489 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 288:
#line 1344 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.mol_type).mol_type = NULL; }
#line 4495 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 289:
#line 1345 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.mol_type) = (yyvsp[0].mol_type); }
#line 4501 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 290:
#line 1349 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.mol_type) = (yyvsp[0].mol_type); (yyval.mol_type).mol_type = (yyvsp[-1].sym); }
#line 4507 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 291:
#line 1352 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_reaction_player_singleton(parse_state, & (yyval.mol_type_list), & (yyvsp[0].mol_type))); }
#line 4513 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 292:
#line 1353 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.mol_type_list) = (yyvsp[-2].mol_type_list); CHECK(mdl_add_reaction_player(parse_state, & (yyval.mol_type_list), & (yyvsp[0].mol_type))); }
#line 4519 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 293:
#line 1356 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.mol_type).mol_type = NULL; (yyval.mol_type).orient_set = 0; }
#line 4525 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 297:
#line 1365 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                        if ((yyvsp[-1].react_rates).forward_rate.rate_type == RATE_UNSET)
                                                        {
                                                          mdlerror(parse_state, "invalid reaction rate specification: must specify a forward rate.");
                                                          return 1;
                                                        }

                                                        (yyval.react_rates) = (yyvsp[-1].react_rates);
                                                      }
#line 4539 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 298:
#line 1376 "../src/../src/mdlparse.y" /* yacc.c:1646  */
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
#line 4558 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 299:
#line 1393 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.react_rates).forward_rate = (yyvsp[0].react_rate); (yyval.react_rates).backward_rate.rate_type = RATE_UNSET; CHECK(mdl_valid_rate(parse_state, &(yyvsp[0].react_rate))); }
#line 4564 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 300:
#line 1394 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.react_rates).forward_rate = (yyvsp[0].react_rate); (yyval.react_rates).backward_rate.rate_type = RATE_UNSET; CHECK(mdl_valid_rate(parse_state, &(yyvsp[0].react_rate))); }
#line 4570 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 301:
#line 1395 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.react_rates).backward_rate = (yyvsp[0].react_rate); (yyval.react_rates).forward_rate.rate_type = RATE_UNSET; CHECK(mdl_valid_rate(parse_state, &(yyvsp[0].react_rate))); }
#line 4576 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 302:
#line 1399 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.react_rate).rate_type = RATE_CONSTANT; (yyval.react_rate).v.rate_constant = (yyvsp[0].dbl); }
#line 4582 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 303:
#line 1400 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.react_rate).rate_type = RATE_FILE; (yyval.react_rate).v.rate_file = (yyvsp[0].str); }
#line 4588 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 304:
#line 1401 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_reaction_rate_from_var(parse_state, & (yyval.react_rate), (yyvsp[0].sym))); }
#line 4594 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 305:
#line 1412 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_release_pattern(parse_state, (yyvsp[-3].sym), &(yyvsp[-1].rpat))); }
#line 4600 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 306:
#line 1415 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.sym) = mdl_new_release_pattern(parse_state, (yyvsp[0].str))); }
#line 4606 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 307:
#line 1418 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.sym) = mdl_existing_release_pattern_or_rxn_pathname(parse_state, (yyvsp[0].str))); }
#line 4612 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 308:
#line 1422 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                        (yyval.rpat).delay = 0;
                                                        (yyval.rpat).release_interval = FOREVER;
                                                        (yyval.rpat).train_interval = FOREVER;
                                                        (yyval.rpat).train_duration = FOREVER;
                                                        (yyval.rpat).number_of_trains = 1;
                                                      }
#line 4624 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 309:
#line 1430 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.rpat) = (yyvsp[-3].rpat); (yyval.rpat).delay = (yyvsp[0].dbl) / parse_state->vol->time_unit; }
#line 4630 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 310:
#line 1432 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.rpat) = (yyvsp[-3].rpat); (yyval.rpat).release_interval = (yyvsp[0].dbl) / parse_state->vol->time_unit; }
#line 4636 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 311:
#line 1434 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.rpat) = (yyvsp[-3].rpat); (yyval.rpat).train_interval = (yyvsp[0].dbl) / parse_state->vol->time_unit; }
#line 4642 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 312:
#line 1436 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.rpat) = (yyvsp[-3].rpat); (yyval.rpat).train_duration = (yyvsp[0].dbl) / parse_state->vol->time_unit; }
#line 4648 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 313:
#line 1438 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.rpat) = (yyvsp[-3].rpat); (yyval.rpat).number_of_trains = (yyvsp[0].ival); }
#line 4654 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 314:
#line 1441 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.ival) = (int) (yyvsp[0].dbl); }
#line 4660 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 315:
#line 1442 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.ival) = INT_MAX; }
#line 4666 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 316:
#line 1449 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->current_object = parse_state->vol->root_instance; }
#line 4672 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 317:
#line 1450 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                        check_regions(parse_state->vol->root_instance, (yyvsp[0].obj));
                                                        add_child_objects(parse_state->vol->root_instance, (yyvsp[0].obj), (yyvsp[0].obj));
                                                        parse_state->current_object = parse_state->vol->root_object;
                                                      }
#line 4682 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 318:
#line 1460 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { add_child_objects(parse_state->vol->root_object, (yyvsp[0].obj), (yyvsp[0].obj)); }
#line 4688 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 324:
#line 1476 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.sym) = mdl_start_object(parse_state, (yyvsp[0].str))); }
#line 4694 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 326:
#line 1482 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { mdl_finish_object(parse_state); }
#line 4700 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 330:
#line 1495 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { transform_translate(parse_state->vol, parse_state->current_object->t_matrix, (yyvsp[0].vec3)); }
#line 4706 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 331:
#line 1496 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { transform_scale(parse_state->current_object->t_matrix, (yyvsp[0].vec3)); }
#line 4712 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 332:
#line 1497 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_transform_rotate(parse_state, parse_state->current_object->t_matrix, (yyvsp[-2].vec3), (yyvsp[0].dbl))); }
#line 4718 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 333:
#line 1506 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                          struct object *the_object = (struct object *) (yyvsp[-5].sym)->value;
                                                          the_object->object_type = META_OBJ;
                                                          add_child_objects(the_object, (yyvsp[-2].obj_list).obj_head, (yyvsp[-2].obj_list).obj_tail);
                                                          (yyval.obj) = the_object;
                                                      }
#line 4729 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 334:
#line 1515 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { mdl_object_list_singleton(& (yyval.obj_list), (yyvsp[0].obj)); }
#line 4735 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 335:
#line 1516 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.obj_list) = (yyvsp[-1].obj_list); mdl_add_object_to_list(& (yyval.obj_list), (yyvsp[0].obj)); }
#line 4741 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 338:
#line 1525 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_deep_copy_object(parse_state, (struct object *) (yyvsp[-3].sym)->value, (struct object *) (yyvsp[-1].sym)->value)); }
#line 4747 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 339:
#line 1527 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.obj) = (struct object *) (yyvsp[-6].sym)->value; }
#line 4753 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 342:
#line 1537 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_start_release_site(parse_state, (yyvsp[-2].sym), SHAPE_UNDEFINED)); }
#line 4759 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 343:
#line 1541 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.obj) = mdl_finish_release_site(parse_state, (yyvsp[-7].sym))); }
#line 4765 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 344:
#line 1544 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_release_site_geometry_region(parse_state, parse_state->current_release_site, parse_state->current_object, (yyvsp[0].rev))); }
#line 4771 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 345:
#line 1545 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_release_site_geometry_object(parse_state, parse_state->current_release_site, (struct object *) (yyvsp[0].sym)->value)); }
#line 4777 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 346:
#line 1546 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->current_release_site->release_shape = SHAPE_SPHERICAL; }
#line 4783 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 347:
#line 1547 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->current_release_site->release_shape = SHAPE_CUBIC; }
#line 4789 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 348:
#line 1548 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->current_release_site->release_shape = SHAPE_ELLIPTIC; }
#line 4795 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 349:
#line 1549 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->current_release_site->release_shape = SHAPE_RECTANGULAR; }
#line 4801 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 350:
#line 1550 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->current_release_site->release_shape = SHAPE_SPHERICAL_SHELL; }
#line 4807 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 351:
#line 1551 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                          parse_state->current_release_site->release_shape = SHAPE_LIST;
                                                          parse_state->current_release_site->release_number_method = CONSTNUM;
                                                      }
#line 4816 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 352:
#line 1558 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.rev) = new_release_region_expr_term((yyvsp[0].sym))); }
#line 4822 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 353:
#line 1559 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.rev) = (yyvsp[-1].rev); }
#line 4828 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 354:
#line 1560 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.rev) = new_release_region_expr_binary((yyvsp[-2].rev), (yyvsp[0].rev), REXP_UNION)); }
#line 4834 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 355:
#line 1561 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.rev) = new_release_region_expr_binary((yyvsp[-2].rev), (yyvsp[0].rev), REXP_SUBTRACTION)); }
#line 4840 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 356:
#line 1562 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.rev) = new_release_region_expr_binary((yyvsp[-2].rev), (yyvsp[0].rev), REXP_INTERSECTION)); }
#line 4846 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 357:
#line 1567 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_start_release_site(parse_state, (yyvsp[-2].sym), (yyvsp[-1].tok))); }
#line 4852 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 358:
#line 1570 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.obj) = mdl_finish_release_site(parse_state, (yyvsp[-6].sym))); }
#line 4858 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 359:
#line 1573 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = SHAPE_SPHERICAL; }
#line 4864 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 360:
#line 1574 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = SHAPE_CUBIC; }
#line 4870 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 361:
#line 1575 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = SHAPE_ELLIPTIC; }
#line 4876 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 362:
#line 1576 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = SHAPE_RECTANGULAR; }
#line 4882 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 363:
#line 1577 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = SHAPE_SPHERICAL_SHELL; }
#line 4888 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 366:
#line 1585 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.sym) = mdl_existing_num_or_array(parse_state, (yyvsp[0].str))); }
#line 4894 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 367:
#line 1589 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { set_release_site_location(parse_state->vol, parse_state->current_release_site, (yyvsp[0].vec3)); }
#line 4900 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 368:
#line 1590 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_release_site_molecule(parse_state, parse_state->current_release_site, & (yyvsp[0].mol_type))); }
#line 4906 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 369:
#line 1591 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                        if (parse_state->current_release_site->release_shape == SHAPE_LIST)
                                                        {
                                                          mdlerror(parse_state, "molecules are already specified in a list--cannot set number or density.");
                                                          return 1;
                                                        }
                                                      }
#line 4918 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 370:
#line 1598 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_release_site_diameter(parse_state, parse_state->current_release_site, (yyvsp[0].dbl) * (((yyvsp[-2].tok) == SITE_RADIUS) ? 2.0 : 1.0))); }
#line 4924 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 371:
#line 1599 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_release_site_diameter_array(parse_state, parse_state->current_release_site, (yyvsp[0].nlist).value_count, (yyvsp[0].nlist).value_head, ((yyvsp[-2].tok) == SITE_RADIUS) ? 2.0 : 1.0)); }
#line 4930 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 372:
#line 1600 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_release_site_diameter_var(parse_state, parse_state->current_release_site, ((yyvsp[-2].tok) == SITE_RADIUS) ? 2.0 : 1.0, (yyvsp[0].sym))); }
#line 4936 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 373:
#line 1601 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_release_site_periodic_box(parse_state, parse_state->current_release_site, (yyvsp[0].vec3))); }
#line 4942 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 374:
#line 1602 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_release_site_probability(parse_state, parse_state->current_release_site, (yyvsp[0].dbl))); }
#line 4948 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 375:
#line 1604 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_release_site_pattern(parse_state, parse_state->current_release_site, (yyvsp[0].sym))); }
#line 4954 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 376:
#line 1606 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_release_site_molecule_positions(parse_state, parse_state->current_release_site, & (yyvsp[-1].rsm_list))); }
#line 4960 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 377:
#line 1610 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = SITE_DIAMETER; }
#line 4966 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 378:
#line 1611 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = SITE_RADIUS; }
#line 4972 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 383:
#line 1623 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { set_release_site_constant_number(parse_state->current_release_site, (yyvsp[0].dbl)); }
#line 4978 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 384:
#line 1626 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { set_release_site_constant_number(parse_state->current_release_site, (yyvsp[-1].dbl)); }
#line 4984 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 385:
#line 1633 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { set_release_site_gaussian_number(parse_state->current_release_site, (yyvsp[-4].dbl), (yyvsp[-1].dbl)); }
#line 4990 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 386:
#line 1641 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { set_release_site_volume_dependent_number(parse_state->current_release_site, (yyvsp[-7].dbl), (yyvsp[-4].dbl), (yyvsp[-1].dbl)); }
#line 4996 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 387:
#line 1645 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_release_site_concentration(parse_state, parse_state->current_release_site, (yyvsp[0].dbl))); }
#line 5002 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 388:
#line 1646 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(set_release_site_density(parse_state->current_release_site, (yyvsp[0].dbl))); }
#line 5008 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 389:
#line 1650 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { release_single_molecule_singleton(& (yyval.rsm_list), (yyvsp[0].rsm)); }
#line 5014 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 390:
#line 1652 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.rsm_list) = (yyvsp[-1].rsm_list); add_release_single_molecule_to_list(& (yyval.rsm_list), (yyvsp[0].rsm)); }
#line 5020 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 391:
#line 1656 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.rsm) = mdl_new_release_single_molecule(parse_state, &(yyvsp[-1].mol_type), (yyvsp[0].vec3))); }
#line 5026 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 393:
#line 1667 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                        CHECKN((yyval.obj) = mdl_new_polygon_list(
                                                          parse_state, (yyvsp[-4].str), (yyvsp[-1].vertlist).vertex_count, (yyvsp[-1].vertlist).vertex_head,
                                                          (yyvsp[0].ecl).connection_count, (yyvsp[0].ecl).connection_head));
                                                      }
#line 5036 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 394:
#line 1676 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                          (yyval.obj) = (struct object *) (yyvsp[-3].obj);
                                                          CHECK(mdl_finish_polygon_list(parse_state, (yyval.obj)));
                                                      }
#line 5045 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 395:
#line 1682 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.vertlist) = (yyvsp[-1].vertlist); }
#line 5051 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 396:
#line 1685 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.vertlistitem) = mdl_new_vertex_list_item((yyvsp[0].vec3))); }
#line 5057 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 397:
#line 1688 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { mdl_vertex_list_singleton(& (yyval.vertlist), (yyvsp[0].vertlistitem)); }
#line 5063 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 398:
#line 1689 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.vertlist) = (yyvsp[-1].vertlist); mdl_add_vertex_to_list(& (yyval.vertlist), (yyvsp[0].vertlistitem)); }
#line 5069 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 399:
#line 1694 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.ecl) = (yyvsp[-1].ecl); }
#line 5075 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 400:
#line 1698 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { mdl_element_connection_list_singleton(& (yyval.ecl), (yyvsp[0].elem_conn)); }
#line 5081 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 401:
#line 1700 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.ecl) = (yyvsp[-1].ecl); mdl_add_element_connection_to_list(& (yyval.ecl), (yyvsp[0].elem_conn)); }
#line 5087 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 402:
#line 1703 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.elem_conn) = mdl_new_element_connection(parse_state, & (yyvsp[0].nlist))); }
#line 5093 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 407:
#line 1719 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN(parse_state->current_region = mdl_get_region(parse_state, parse_state->current_object, "REMOVED")); }
#line 5099 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 408:
#line 1721 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                          parse_state->current_region->element_list_head = (yyvsp[-1].elem_list).elml_head;
                                                          if (parse_state->current_object->object_type == POLY_OBJ)
                                                          {
                                                            CHECK(mdl_normalize_elements(parse_state, parse_state->current_region,0));
                                                          }
                                                      }
#line 5111 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 411:
#line 1735 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = Z_POS; }
#line 5117 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 412:
#line 1736 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = Z_NEG; }
#line 5123 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 413:
#line 1737 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = Y_NEG; }
#line 5129 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 414:
#line 1738 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = Y_POS; }
#line 5135 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 415:
#line 1739 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = X_NEG; }
#line 5141 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 416:
#line 1740 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = X_POS; }
#line 5147 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 417:
#line 1741 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = ALL_SIDES; }
#line 5153 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 419:
#line 1747 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.elem_list) = (yyvsp[-1].elem_list); mdl_add_elements_to_list(& (yyval.elem_list), (yyvsp[0].elem_list).elml_head, (yyvsp[0].elem_list).elml_tail); }
#line 5159 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 422:
#line 1753 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.elem_list).elml_tail = (yyval.elem_list).elml_head = (yyvsp[0].elem_list_item); }
#line 5165 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 423:
#line 1754 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.elem_list).elml_tail = (yyval.elem_list).elml_head = (yyvsp[0].elem_list_item); }
#line 5171 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 424:
#line 1759 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.elem_list) = (yyvsp[-1].elem_list); }
#line 5177 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 425:
#line 1764 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.elem_list) = (yyvsp[-1].elem_list); mdl_set_elements_to_exclude((yyval.elem_list).elml_head); }
#line 5183 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 427:
#line 1771 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.elem_list).elml_tail = (yyval.elem_list).elml_head = (yyvsp[0].elem_list_item); }
#line 5189 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 428:
#line 1772 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.elem_list) = (yyvsp[-2].elem_list); mdl_add_elements_to_list(& (yyval.elem_list), (yyvsp[0].elem_list_item), (yyvsp[0].elem_list_item)); }
#line 5195 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 429:
#line 1775 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.elem_list_item) = new_element_list((unsigned int) (yyvsp[0].dbl), (unsigned int) (yyvsp[0].dbl))); }
#line 5201 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 430:
#line 1776 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.elem_list_item) = new_element_list((unsigned int) (yyvsp[-2].dbl), (unsigned int) (yyvsp[0].dbl))); }
#line 5207 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 431:
#line 1777 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.elem_list_item) = mdl_new_element_side(parse_state, (yyvsp[0].tok))); }
#line 5213 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 432:
#line 1780 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.elem_list_item) = mdl_new_element_previous_region(parse_state, parse_state->current_object, parse_state->current_region, (yyvsp[0].str), (yyvsp[-2].tok))); }
#line 5219 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 433:
#line 1783 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = 0; }
#line 5225 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 434:
#line 1784 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = 1; }
#line 5231 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 435:
#line 1787 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.elem_list_item) = mdl_new_element_patch(parse_state, parse_state->current_polygon, (yyvsp[-2].vec3), (yyvsp[0].vec3), (yyvsp[-4].tok))); }
#line 5237 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 436:
#line 1790 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = 0; }
#line 5243 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 437:
#line 1791 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = 1; }
#line 5249 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 441:
#line 1807 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->current_region = (yyvsp[-1].reg); }
#line 5255 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 442:
#line 1808 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_region_elements(parse_state, (yyvsp[-3].reg), (yyvsp[0].elem_list).elml_head, (yyvsp[-3].reg)->parent->object_type == POLY_OBJ)); }
#line 5261 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 443:
#line 1810 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->current_region = NULL; }
#line 5267 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 444:
#line 1818 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                        CHECKN(mdl_new_voxel_list(parse_state, (yyvsp[-4].sym),
                                                                                  (yyvsp[-1].vertlist).vertex_count, (yyvsp[-1].vertlist).vertex_head,
                                                                                  (yyvsp[0].ecl).connection_count, (yyvsp[0].ecl).connection_head));
                                                      }
#line 5277 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 445:
#line 1824 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.obj) = (struct object *) (yyvsp[-7].sym)->value; }
#line 5283 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 446:
#line 1829 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.ecl) = (yyvsp[-1].ecl); }
#line 5289 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 447:
#line 1832 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.elem_conn) = mdl_new_tet_element_connection(parse_state, & (yyvsp[0].nlist))); }
#line 5295 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 448:
#line 1836 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                          (yyval.ecl).connection_head = (yyval.ecl).connection_tail = (yyvsp[0].elem_conn);
                                                          (yyval.ecl).connection_count = 1;
                                                      }
#line 5304 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 449:
#line 1840 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                          (yyval.ecl) = (yyvsp[-1].ecl);
                                                          (yyval.ecl).connection_tail = (yyval.ecl).connection_tail->next = (yyvsp[0].elem_conn);
                                                          ++ (yyval.ecl).connection_count;
                                                      }
#line 5314 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 450:
#line 1850 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->vol->periodic_traditional = (yyvsp[0].tok); }
#line 5320 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 451:
#line 1853 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN(mdl_create_periodic_box(parse_state, (yyvsp[-7].vec3), (yyvsp[-5].vec3), (yyvsp[-2].tok), (yyvsp[-1].tok), (yyvsp[0].tok))); }
#line 5326 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 452:
#line 1854 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_finish_periodic_box(parse_state)); }
#line 5332 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 453:
#line 1860 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN(mdl_new_box_object(parse_state, (yyvsp[-8].sym), (yyvsp[-3].vec3), (yyvsp[-1].vec3))); }
#line 5338 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 454:
#line 1861 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_triangulate_box_object(parse_state, (yyvsp[-10].sym), parse_state->current_polygon, (yyvsp[-2].dbl))); }
#line 5344 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 455:
#line 1863 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                          CHECK(mdl_finish_box_object(parse_state, (yyvsp[-13].sym)));
                                                          (yyval.obj) = (struct object *) (yyvsp[-13].sym)->value;
                                                      }
#line 5353 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 456:
#line 1870 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = 0; }
#line 5359 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 457:
#line 1871 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = (yyvsp[0].tok); }
#line 5365 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 458:
#line 1875 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = 0; }
#line 5371 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 459:
#line 1876 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = (yyvsp[0].tok); }
#line 5377 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 460:
#line 1880 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = 0; }
#line 5383 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 461:
#line 1881 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = (yyvsp[0].tok); }
#line 5389 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 462:
#line 1885 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = 0; }
#line 5395 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 463:
#line 1886 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = (yyvsp[0].tok); }
#line 5401 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 464:
#line 1889 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.dbl) = 0.0; }
#line 5407 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 465:
#line 1890 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                        (yyval.dbl) = (yyvsp[0].dbl);
                                                        if ((yyval.dbl) < 2.0)
                                                        {
                                                          mdlerror(parse_state, "invalid aspect ratio requested (must be greater than or equal to 2.0)");
                                                          return 1;
                                                        }
                                                      }
#line 5420 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 469:
#line 1916 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_start_existing_obj_region_def(parse_state, (yyvsp[0].sym))); }
#line 5426 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 470:
#line 1917 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->current_region = (yyvsp[-1].reg); }
#line 5432 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 471:
#line 1919 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { mdl_set_region_elements(parse_state, (yyvsp[-4].reg), (yyvsp[0].elem_list).elml_head, 1); }
#line 5438 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 472:
#line 1921 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                          parse_state->current_region = NULL;
                                                          parse_state->current_polygon = NULL;
                                                          parse_state->current_object = parse_state->vol->root_object;
                                                      }
#line 5448 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 473:
#line 1928 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.reg) = mdl_create_region(parse_state, parse_state->current_object, (yyvsp[0].str))); }
#line 5454 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 477:
#line 1939 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { mdl_add_surf_mol_to_region(parse_state->current_region, & (yyvsp[0].surf_mol_dat_list)); }
#line 5460 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 478:
#line 1943 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { mdl_set_region_surface_class(parse_state, parse_state->current_region, (yyvsp[0].sym)); }
#line 5466 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 482:
#line 1962 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->current_region = (struct region *) (yyvsp[-1].sym)->value; }
#line 5472 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 483:
#line 1964 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->current_region = NULL; }
#line 5478 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 484:
#line 1972 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                          parse_state->header_comment = NULL;  /* No header by default */
                                                          parse_state->exact_time_flag = 1;    /* Print exact_time column in TRIGGER output by default */
                                                      }
#line 5487 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 485:
#line 1978 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_add_reaction_output_block_to_world(parse_state, (int) (yyvsp[-4].dbl), & (yyvsp[-2].ro_otimes), & (yyvsp[-1].ro_sets))); }
#line 5493 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 486:
#line 1982 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.dbl) = COUNTBUFFERSIZE; }
#line 5499 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 487:
#line 1983 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                          double temp_value = (yyvsp[0].dbl);
                                                          if (!(temp_value >= 1.0 && temp_value < UINT_MAX))
                                                          {
                                                            mdlerror_fmt(parse_state, "Requested buffer size of %.15g lines is invalid.  Suggested range is 100-1000000.", temp_value);
                                                            return 1;
                                                          }
                                                          (yyval.dbl) = (yyvsp[0].dbl);
                                                      }
#line 5513 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 491:
#line 1999 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.ro_otimes).type = OUTPUT_BY_STEP; (yyval.ro_otimes).step = (yyvsp[0].dbl); }
#line 5519 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 492:
#line 2003 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                        (yyval.ro_otimes).type = OUTPUT_BY_ITERATION_LIST;
                                                        (yyval.ro_otimes).values = (yyvsp[0].nlist);
                                                      }
#line 5528 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 493:
#line 2011 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                        (yyval.ro_otimes).type = OUTPUT_BY_TIME_LIST;
                                                        (yyval.ro_otimes).values = (yyvsp[0].nlist);
                                                      }
#line 5537 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 494:
#line 2018 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.ro_sets).set_head = (yyval.ro_sets).set_tail = (yyvsp[0].ro_set); }
#line 5543 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 495:
#line 2020 "../src/../src/mdlparse.y" /* yacc.c:1646  */
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
#line 5558 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 497:
#line 2034 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.ro_set) = NULL; }
#line 5564 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 498:
#line 2035 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.ro_set) = NULL; }
#line 5570 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 499:
#line 2039 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {  parse_state->count_flags = 0; }
#line 5576 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 500:
#line 2041 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.ro_set) = mdl_populate_output_set(parse_state, parse_state->header_comment, parse_state->exact_time_flag, (yyvsp[-3].ro_cols).column_head, (yyvsp[-1].tok), (yyvsp[0].str))); }
#line 5582 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 501:
#line 2045 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.str) = NULL; }
#line 5588 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 502:
#line 2046 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.str) = ((yyvsp[0].tok) ? "" : NULL); }
#line 5594 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 503:
#line 2047 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.str) = (yyvsp[0].str); }
#line 5600 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 504:
#line 2051 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->header_comment = (yyvsp[0].str); }
#line 5606 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 505:
#line 2055 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->exact_time_flag = (yyvsp[0].tok); }
#line 5612 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 507:
#line 2061 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                          (yyval.ro_cols) = (yyvsp[-2].ro_cols);
                                                          (yyval.ro_cols).column_tail->next = (yyvsp[0].ro_cols).column_head;
                                                          (yyval.ro_cols).column_tail = (yyvsp[0].ro_cols).column_tail;
                                                      }
#line 5622 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 508:
#line 2069 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_single_count_expr(parse_state, & (yyval.ro_cols), (yyvsp[-1].cnt), (yyvsp[0].str))); }
#line 5628 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 509:
#line 2073 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.cnt) = mdl_new_oexpr_constant(parse_state, (yyvsp[0].dbl))); }
#line 5634 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 511:
#line 2075 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.cnt) = mdl_join_oexpr_tree(parse_state, (yyvsp[-1].cnt), NULL, '(')); }
#line 5640 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 512:
#line 2076 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.cnt) = mdl_join_oexpr_tree(parse_state, (yyvsp[-2].cnt),   (yyvsp[0].cnt), '+')); }
#line 5646 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 513:
#line 2077 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.cnt) = mdl_join_oexpr_tree(parse_state, (yyvsp[-2].cnt),   (yyvsp[0].cnt), '-')); }
#line 5652 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 514:
#line 2078 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.cnt) = mdl_join_oexpr_tree(parse_state, (yyvsp[-2].cnt),   (yyvsp[0].cnt), '*')); }
#line 5658 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 515:
#line 2079 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.cnt) = mdl_join_oexpr_tree(parse_state, (yyvsp[-2].cnt),   (yyvsp[0].cnt), '/')); }
#line 5664 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 516:
#line 2080 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.cnt) = mdl_join_oexpr_tree(parse_state, (yyvsp[0].cnt), NULL, '_')); }
#line 5670 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 517:
#line 2081 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.cnt) = mdl_sum_oexpr((yyvsp[-1].cnt))); }
#line 5676 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 518:
#line 2086 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->count_flags |= COUNT_PRESENT; }
#line 5682 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 519:
#line 2087 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.cnt) = (yyvsp[-1].cnt); }
#line 5688 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 520:
#line 2088 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.cnt) = mdl_new_oexpr_constant(parse_state, (yyvsp[-1].dbl))); }
#line 5694 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 521:
#line 2089 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->count_flags |= TRIGGER_PRESENT; }
#line 5700 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 522:
#line 2090 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.cnt) = (yyvsp[-1].cnt); }
#line 5706 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 523:
#line 2093 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = FILE_OVERWRITE; }
#line 5712 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 524:
#line 2094 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = FILE_SUBSTITUTE; }
#line 5718 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 525:
#line 2095 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = FILE_APPEND; }
#line 5724 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 526:
#line 2096 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = FILE_APPEND_HEADER; }
#line 5730 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 527:
#line 2097 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = FILE_CREATE; }
#line 5736 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 529:
#line 2103 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.sym) = mdl_existing_rxn_pathname_or_molecule(parse_state, (yyvsp[0].str))); }
#line 5742 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 530:
#line 2107 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                        (yyval.mol_type) = (yyvsp[0].mol_type);
                                                        if ((yyval.mol_type).orient > 0)
                                                          (yyval.mol_type).orient = 1;
                                                        else if ((yyval.mol_type).orient < 0)
                                                          (yyval.mol_type).orient = -1;
                                                        CHECKN((yyval.mol_type).mol_type = mdl_existing_molecule(parse_state, (yyvsp[-1].str)));
                                                      }
#line 5755 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 537:
#line 2127 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.cnt) = mdl_count_syntax_1(parse_state, (yyvsp[-3].sym), (yyvsp[-1].sym), (yyvsp[0].tok), parse_state->count_flags)); }
#line 5761 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 538:
#line 2132 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.cnt) = mdl_count_syntax_2(parse_state, (yyvsp[-3].mol_type).mol_type, (yyvsp[-3].mol_type).orient, (yyvsp[-1].sym), (yyvsp[0].tok), parse_state->count_flags)); }
#line 5767 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 539:
#line 2137 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.cnt) = mdl_count_syntax_3(parse_state, (yyvsp[-3].str), (yyvsp[-1].sym), (yyvsp[0].tok), parse_state->count_flags)); }
#line 5773 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 540:
#line 2143 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.cnt) = mdl_count_syntax_periodic_1(parse_state, (yyvsp[-5].sym), (yyvsp[-3].sym), (yyvsp[-1].vec3), (yyvsp[0].tok), parse_state->count_flags)); }
#line 5779 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 541:
#line 2147 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.cnt) = mdl_count_syntax_periodic_2(parse_state, (yyvsp[-5].mol_type).mol_type, (yyvsp[-5].mol_type).orient, (yyvsp[-3].sym), (yyvsp[-1].vec3), (yyvsp[0].tok), parse_state->count_flags)); }
#line 5785 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 542:
#line 2152 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.cnt) = mdl_count_syntax_periodic_3(parse_state, (yyvsp[-5].str), (yyvsp[-3].sym), (yyvsp[-1].vec3), (yyvsp[0].tok), parse_state->count_flags)); }
#line 5791 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 543:
#line 2155 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.sym) = NULL; }
#line 5797 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 544:
#line 2156 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.sym) = (yyvsp[0].sym); }
#line 5803 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 545:
#line 2157 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.sym) = (yyvsp[0].sym); }
#line 5809 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 546:
#line 2160 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = REPORT_NOTHING; }
#line 5815 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 547:
#line 2161 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = (yyvsp[0].tok); }
#line 5821 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 548:
#line 2164 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = REPORT_FRONT_HITS; }
#line 5827 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 549:
#line 2165 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = REPORT_BACK_HITS; }
#line 5833 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 550:
#line 2166 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = REPORT_ALL_HITS; }
#line 5839 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 551:
#line 2167 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = REPORT_FRONT_CROSSINGS; }
#line 5845 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 552:
#line 2168 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = REPORT_BACK_CROSSINGS; }
#line 5851 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 553:
#line 2169 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = REPORT_ALL_CROSSINGS; }
#line 5857 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 554:
#line 2170 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = REPORT_CONCENTRATION; }
#line 5863 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 555:
#line 2171 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = REPORT_ENCLOSED; }
#line 5869 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 556:
#line 2174 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.str) = NULL; }
#line 5875 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 557:
#line 2175 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.str) = (yyvsp[0].str); }
#line 5881 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 558:
#line 2182 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_new_viz_output_block(parse_state)); }
#line 5887 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 559:
#line 2185 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { }
#line 5893 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 562:
#line 2194 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_viz_mode(parse_state->vol->viz_blocks, CELLBLENDER_MODE)); }
#line 5899 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 563:
#line 2195 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_viz_mode(parse_state->vol->viz_blocks, (yyvsp[0].ival))); }
#line 5905 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 564:
#line 2198 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.ival) = NO_VIZ_MODE; }
#line 5911 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 565:
#line 2199 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.ival) = ASCII_MODE; }
#line 5917 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 566:
#line 2200 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.ival) = CELLBLENDER_MODE; }
#line 5923 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 568:
#line 2205 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                        if ((yyvsp[0].frame_list).frame_head)
                                                        {
                                                          (yyvsp[0].frame_list).frame_tail->next = parse_state->vol->viz_blocks->frame_data_head;
                                                          parse_state->vol->viz_blocks->frame_data_head = (yyvsp[0].frame_list).frame_head;
                                                        }
                                                      }
#line 5935 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 570:
#line 2218 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_viz_filename_prefix(parse_state, parse_state->vol->viz_blocks, (yyvsp[0].str))); }
#line 5941 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 571:
#line 2224 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.frame_list) = (yyvsp[-1].frame_list); }
#line 5947 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 573:
#line 2230 "../src/../src/mdlparse.y" /* yacc.c:1646  */
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
#line 5963 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 574:
#line 2244 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.frame_list).frame_head = (yyval.frame_list).frame_tail = NULL; }
#line 5969 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 578:
#line 2256 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_viz_state(parse_state, & (yyval.ival), (yyvsp[0].dbl))); }
#line 5975 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 579:
#line 2257 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.ival) = INCLUDE_OBJ; }
#line 5981 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 582:
#line 2267 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_viz_include_molecules(parse_state, parse_state->vol->viz_blocks, (yyvsp[-1].symlist), (yyvsp[0].ival))); }
#line 5987 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 583:
#line 2268 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_viz_include_all_molecules(parse_state->vol->viz_blocks, (yyvsp[0].ival))); }
#line 5993 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 584:
#line 2272 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.symlist) = mdl_existing_molecule_list(parse_state, (yyvsp[0].str))); }
#line 5999 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 585:
#line 2273 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.symlist) = mdl_existing_molecules_wildcard(parse_state, (yyvsp[0].str))); }
#line 6005 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 586:
#line 2277 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_new_viz_all_times(parse_state, & (yyval.nlist))); }
#line 6011 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 588:
#line 2283 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.frame_list) = (yyvsp[-1].frame_list); }
#line 6017 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 590:
#line 2289 "../src/../src/mdlparse.y" /* yacc.c:1646  */
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
#line 6035 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 591:
#line 2306 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_new_viz_mol_frames(parse_state, parse_state->vol->viz_blocks, & (yyval.frame_list), OUTPUT_BY_TIME_LIST, (yyvsp[-2].tok), & (yyvsp[0].nlist))); }
#line 6041 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 592:
#line 2310 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_new_viz_all_iterations(parse_state, & (yyval.nlist))); }
#line 6047 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 594:
#line 2317 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.frame_list) = (yyvsp[-1].frame_list); }
#line 6053 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 596:
#line 2323 "../src/../src/mdlparse.y" /* yacc.c:1646  */
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
#line 6071 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 597:
#line 2340 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_new_viz_mol_frames(parse_state, parse_state->vol->viz_blocks, & (yyval.frame_list), OUTPUT_BY_ITERATION_LIST, (yyvsp[-2].tok), & (yyvsp[0].nlist))); }
#line 6077 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 598:
#line 2343 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = ALL_MOL_DATA; }
#line 6083 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 599:
#line 2344 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = MOL_POS; }
#line 6089 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 600:
#line 2345 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = MOL_ORIENT; }
#line 6095 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 601:
#line 2359 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                          struct volume_output_item *vo;
                                                          CHECKN(vo = mdl_new_volume_output_item(parse_state, (yyvsp[-6].str), & (yyvsp[-5].species_lst), (yyvsp[-4].vec3), (yyvsp[-3].vec3), (yyvsp[-2].vec3), (yyvsp[-1].otimes)));
                                                          vo->next = parse_state->vol->volume_output_head;
                                                          parse_state->vol->volume_output_head = vo;
                                                      }
#line 6106 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 602:
#line 2368 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.str) = (yyvsp[0].str); }
#line 6112 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 604:
#line 2374 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                          (yyval.species_lst) = (yyvsp[-1].species_lst);
                                                          (yyval.species_lst).species_count += (yyvsp[0].species_lst).species_count;
                                                          (yyval.species_lst).species_tail->next = (yyvsp[0].species_lst).species_head;
                                                          (yyval.species_lst).species_tail = (yyvsp[0].species_lst).species_tail;
                                                      }
#line 6123 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 605:
#line 2383 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.species_lst) = (yyvsp[0].species_lst); }
#line 6129 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 606:
#line 2386 "../src/../src/mdlparse.y" /* yacc.c:1646  */
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
#line 6149 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 607:
#line 2404 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.species_lst).species_tail = (yyval.species_lst).species_head = (yyvsp[0].species_lst_item); (yyval.species_lst).species_count = 1; }
#line 6155 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 608:
#line 2406 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    {
                                                        (yyval.species_lst) = (yyvsp[-2].species_lst);
                                                        (yyval.species_lst).species_tail = (yyval.species_lst).species_tail->next = (yyvsp[0].species_lst_item);
                                                        ++ (yyval.species_lst).species_count;
                                                      }
#line 6165 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 609:
#line 2414 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.vec3) = (yyvsp[0].vec3); }
#line 6171 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 610:
#line 2418 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.vec3) = (yyvsp[0].vec3); }
#line 6177 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 611:
#line 2422 "../src/../src/mdlparse.y" /* yacc.c:1646  */
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
#line 6200 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 612:
#line 2443 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.otimes) = mdl_new_output_times_default(parse_state)); }
#line 6206 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 613:
#line 2444 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.otimes) = mdl_new_output_times_step(parse_state, (yyvsp[0].dbl))); }
#line 6212 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 614:
#line 2445 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.otimes) = mdl_new_output_times_iterations(parse_state, & (yyvsp[0].nlist))); }
#line 6218 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 615:
#line 2446 "../src/../src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.otimes) = mdl_new_output_times_time(parse_state, & (yyvsp[0].nlist))); }
#line 6224 "mdlparse.c" /* yacc.c:1646  */
    break;


#line 6228 "mdlparse.c" /* yacc.c:1646  */
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
#line 2449 "../src/../src/mdlparse.y" /* yacc.c:1906  */






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
