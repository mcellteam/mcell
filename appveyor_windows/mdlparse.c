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
#line 1 "./mdlparse.y" /* yacc.c:339  */

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
    PI_TOK = 413,
    POLYGON_LIST = 414,
    POSITIONS = 415,
    PRINTF = 416,
    PRINT_TIME = 417,
    PROBABILITY_REPORT = 418,
    PROBABILITY_REPORT_THRESHOLD = 419,
    PROGRESS_REPORT = 420,
    RADIAL_DIRECTIONS = 421,
    RADIAL_SUBDIVISIONS = 422,
    RAND_GAUSSIAN = 423,
    RAND_UNIFORM = 424,
    REACTION_DATA_OUTPUT = 425,
    REACTION_OUTPUT_REPORT = 426,
    REAL = 427,
    RECTANGULAR_RELEASE_SITE = 428,
    RECTANGULAR_TOKEN = 429,
    REFLECTIVE = 430,
    RELEASE_EVENT_REPORT = 431,
    RELEASE_INTERVAL = 432,
    RELEASE_PATTERN = 433,
    RELEASE_PROBABILITY = 434,
    RELEASE_SITE = 435,
    REMOVE_ELEMENTS = 436,
    RIGHT = 437,
    ROTATE = 438,
    ROUND_OFF = 439,
    SCALE = 440,
    SEED = 441,
    SHAPE = 442,
    SHOW_EXACT_TIME = 443,
    SIN = 444,
    SITE_DIAMETER = 445,
    SITE_RADIUS = 446,
    SPACE_STEP = 447,
    SPHERICAL = 448,
    SPHERICAL_RELEASE_SITE = 449,
    SPHERICAL_SHELL = 450,
    SPHERICAL_SHELL_SITE = 451,
    SPRINTF = 452,
    SQRT = 453,
    STANDARD_DEVIATION = 454,
    STEP = 455,
    STRING_TO_NUM = 456,
    STR_VALUE = 457,
    SUBUNIT = 458,
    SUBUNIT_RELATIONSHIPS = 459,
    SUMMATION_OPERATOR = 460,
    SURFACE_CLASS = 461,
    SURFACE_ONLY = 462,
    TAN = 463,
    TARGET_ONLY = 464,
    TET_ELEMENT_CONNECTIONS = 465,
    THROUGHPUT_REPORT = 466,
    TIME_LIST = 467,
    TIME_POINTS = 468,
    TIME_STEP = 469,
    TIME_STEP_MAX = 470,
    TO = 471,
    TOP = 472,
    TRAIN_DURATION = 473,
    TRAIN_INTERVAL = 474,
    TRANSLATE = 475,
    TRANSPARENT = 476,
    TRIGGER = 477,
    TRUE = 478,
    UNLIMITED = 479,
    USELESS_VOLUME_ORIENTATION = 480,
    VACANCY_SEARCH_DISTANCE = 481,
    VAR = 482,
    VARYING_PROBABILITY_REPORT = 483,
    VERTEX_LIST = 484,
    VIZ_OUTPUT = 485,
    VIZ_OUTPUT_REPORT = 486,
    VIZ_VALUE = 487,
    VOLUME_DATA_OUTPUT = 488,
    VOLUME_OUTPUT_REPORT = 489,
    VOLUME_DEPENDENT_RELEASE_NUMBER = 490,
    VOLUME_ONLY = 491,
    VOXEL_COUNT = 492,
    VOXEL_LIST = 493,
    VOXEL_SIZE = 494,
    WARNING = 495,
    WARNINGS = 496,
    WORLD = 497,
    YES = 498,
    UNARYMINUS = 499
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
#define PI_TOK 413
#define POLYGON_LIST 414
#define POSITIONS 415
#define PRINTF 416
#define PRINT_TIME 417
#define PROBABILITY_REPORT 418
#define PROBABILITY_REPORT_THRESHOLD 419
#define PROGRESS_REPORT 420
#define RADIAL_DIRECTIONS 421
#define RADIAL_SUBDIVISIONS 422
#define RAND_GAUSSIAN 423
#define RAND_UNIFORM 424
#define REACTION_DATA_OUTPUT 425
#define REACTION_OUTPUT_REPORT 426
#define REAL 427
#define RECTANGULAR_RELEASE_SITE 428
#define RECTANGULAR_TOKEN 429
#define REFLECTIVE 430
#define RELEASE_EVENT_REPORT 431
#define RELEASE_INTERVAL 432
#define RELEASE_PATTERN 433
#define RELEASE_PROBABILITY 434
#define RELEASE_SITE 435
#define REMOVE_ELEMENTS 436
#define RIGHT 437
#define ROTATE 438
#define ROUND_OFF 439
#define SCALE 440
#define SEED 441
#define SHAPE 442
#define SHOW_EXACT_TIME 443
#define SIN 444
#define SITE_DIAMETER 445
#define SITE_RADIUS 446
#define SPACE_STEP 447
#define SPHERICAL 448
#define SPHERICAL_RELEASE_SITE 449
#define SPHERICAL_SHELL 450
#define SPHERICAL_SHELL_SITE 451
#define SPRINTF 452
#define SQRT 453
#define STANDARD_DEVIATION 454
#define STEP 455
#define STRING_TO_NUM 456
#define STR_VALUE 457
#define SUBUNIT 458
#define SUBUNIT_RELATIONSHIPS 459
#define SUMMATION_OPERATOR 460
#define SURFACE_CLASS 461
#define SURFACE_ONLY 462
#define TAN 463
#define TARGET_ONLY 464
#define TET_ELEMENT_CONNECTIONS 465
#define THROUGHPUT_REPORT 466
#define TIME_LIST 467
#define TIME_POINTS 468
#define TIME_STEP 469
#define TIME_STEP_MAX 470
#define TO 471
#define TOP 472
#define TRAIN_DURATION 473
#define TRAIN_INTERVAL 474
#define TRANSLATE 475
#define TRANSPARENT 476
#define TRIGGER 477
#define TRUE 478
#define UNLIMITED 479
#define USELESS_VOLUME_ORIENTATION 480
#define VACANCY_SEARCH_DISTANCE 481
#define VAR 482
#define VARYING_PROBABILITY_REPORT 483
#define VERTEX_LIST 484
#define VIZ_OUTPUT 485
#define VIZ_OUTPUT_REPORT 486
#define VIZ_VALUE 487
#define VOLUME_DATA_OUTPUT 488
#define VOLUME_OUTPUT_REPORT 489
#define VOLUME_DEPENDENT_RELEASE_NUMBER 490
#define VOLUME_ONLY 491
#define VOXEL_COUNT 492
#define VOXEL_LIST 493
#define VOXEL_SIZE 494
#define WARNING 495
#define WARNINGS 496
#define WORLD 497
#define YES 498
#define UNARYMINUS 499

/* Value type.  */
#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED
typedef union YYSTYPE YYSTYPE;
union YYSTYPE
{
#line 68 "./mdlparse.y" /* yacc.c:355  */

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


#line 732 "mdlparse.c" /* yacc.c:355  */
};
# define YYSTYPE_IS_TRIVIAL 1
# define YYSTYPE_IS_DECLARED 1
#endif



int mdlparse (struct mdlparse_vars *parse_state, yyscan_t scanner);

#endif /* !YY_MDL_MDLPARSE_H_INCLUDED  */

/* Copy the second part of user declarations.  */

#line 746 "mdlparse.c" /* yacc.c:358  */

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
#define YYLAST   2521

/* YYNTOKENS -- Number of terminals.  */
#define YYNTOKENS  265
/* YYNNTS -- Number of nonterminals.  */
#define YYNNTS  279
/* YYNRULES -- Number of rules.  */
#define YYNRULES  596
/* YYNSTATES -- Number of states.  */
#define YYNSTATES  1172

/* YYTRANSLATE[YYX] -- Symbol number corresponding to YYX as returned
   by yylex, with out-of-bounds checking.  */
#define YYUNDEFTOK  2
#define YYMAXUTOK   499

#define YYTRANSLATE(YYX)                                                \
  ((unsigned int) (YYX) <= YYMAXUTOK ? yytranslate[YYX] : YYUNDEFTOK)

/* YYTRANSLATE[TOKEN-NUM] -- Symbol number corresponding to TOKEN-NUM
   as returned by yylex, without out-of-bounds checking.  */
static const yytype_uint16 yytranslate[] =
{
       0,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,   245,   256,
     260,   261,   249,   247,   257,   248,     2,   250,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,   246,   255,
     263,   244,   262,     2,   264,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,   253,     2,   254,   252,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,   258,     2,   259,     2,     2,     2,     2,
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
     235,   236,   237,   238,   239,   240,   241,   242,   243,   251
};

#if YYDEBUG
  /* YYRLINE[YYN] -- Source line where rule number YYN was defined.  */
static const yytype_uint16 yyrline[] =
{
       0,   579,   579,   583,   584,   589,   590,   591,   592,   593,
     594,   595,   596,   597,   598,   599,   600,   601,   602,   603,
     604,   605,   606,   607,   608,   613,   616,   619,   622,   625,
     628,   631,   632,   635,   636,   637,   638,   639,   640,   643,
     644,   645,   646,   650,   651,   652,   662,   674,   677,   680,
     692,   693,   706,   707,   713,   736,   737,   738,   739,   742,
     745,   748,   749,   760,   763,   766,   767,   770,   771,   774,
     775,   778,   779,   782,   786,   787,   788,   789,   790,   791,
     792,   793,   794,   795,   796,   797,   798,   799,   800,   801,
     802,   803,   804,   805,   806,   807,   808,   809,   810,   811,
     812,   813,   814,   815,   819,   820,   824,   825,   826,   827,
     830,   836,   837,   838,   839,   840,   841,   842,   845,   849,
     852,   855,   858,   861,   864,   865,   875,   876,   877,   889,
     893,   899,   904,   908,   917,   921,   922,   926,   927,   928,
     929,   930,   931,   932,   933,   934,   935,   936,   937,   938,
     939,   940,   941,   942,   946,   947,   951,   955,   956,   963,
     967,   968,   972,   973,   974,   975,   976,   977,   978,   979,
     980,   981,   982,   983,   984,   985,   986,   987,   991,   992,
     993,   999,  1000,  1001,  1002,  1003,  1007,  1008,  1009,  1013,
    1014,  1015,  1016,  1024,  1025,  1026,  1027,  1028,  1029,  1030,
    1031,  1032,  1033,  1034,  1035,  1036,  1037,  1038,  1045,  1046,
    1047,  1048,  1052,  1056,  1057,  1058,  1065,  1066,  1069,  1073,
    1077,  1078,  1082,  1090,  1093,  1097,  1098,  1102,  1103,  1112,
    1123,  1124,  1128,  1129,  1139,  1143,  1147,  1160,  1161,  1166,
    1170,  1176,  1177,  1182,  1182,  1187,  1190,  1192,  1197,  1198,
    1202,  1205,  1211,  1216,  1217,  1218,  1221,  1222,  1225,  1229,
    1233,  1240,  1244,  1253,  1257,  1266,  1273,  1278,  1279,  1282,
    1285,  1286,  1289,  1290,  1291,  1294,  1299,  1305,  1306,  1307,
    1308,  1311,  1312,  1316,  1321,  1322,  1325,  1329,  1330,  1334,
    1337,  1338,  1341,  1342,  1346,  1347,  1350,  1361,  1378,  1379,
    1380,  1384,  1385,  1386,  1393,  1400,  1403,  1407,  1414,  1416,
    1418,  1420,  1422,  1426,  1427,  1434,  1434,  1445,  1448,  1449,
    1450,  1451,  1452,  1461,  1464,  1467,  1470,  1472,  1476,  1480,
    1481,  1482,  1487,  1500,  1501,  1504,  1505,  1510,  1509,  1516,
    1517,  1522,  1521,  1529,  1530,  1531,  1532,  1533,  1534,  1535,
    1536,  1543,  1544,  1545,  1546,  1547,  1552,  1551,  1558,  1559,
    1560,  1561,  1562,  1566,  1567,  1570,  1574,  1575,  1576,  1583,
    1584,  1585,  1586,  1587,  1589,  1594,  1595,  1599,  1600,  1601,
    1602,  1607,  1608,  1614,  1621,  1629,  1630,  1634,  1635,  1640,
    1643,  1651,  1648,  1666,  1669,  1672,  1673,  1677,  1682,  1683,
    1687,  1690,  1692,  1698,  1699,  1703,  1703,  1715,  1716,  1719,
    1720,  1721,  1722,  1723,  1724,  1725,  1729,  1730,  1735,  1736,
    1737,  1738,  1742,  1747,  1751,  1755,  1756,  1759,  1760,  1761,
    1764,  1767,  1768,  1771,  1774,  1775,  1779,  1785,  1786,  1791,
    1792,  1791,  1802,  1799,  1812,  1816,  1820,  1824,  1835,  1836,
    1832,  1844,  1845,  1859,  1865,  1866,  1871,  1872,  1874,  1871,
    1883,  1886,  1888,  1893,  1894,  1898,  1905,  1911,  1912,  1917,
    1917,  1927,  1926,  1937,  1938,  1949,  1950,  1951,  1954,  1958,
    1966,  1973,  1974,  1988,  1989,  1990,  1994,  1994,  2000,  2001,
    2002,  2006,  2010,  2014,  2015,  2024,  2028,  2029,  2030,  2031,
    2032,  2033,  2034,  2035,  2036,  2041,  2041,  2043,  2044,  2044,
    2048,  2049,  2050,  2051,  2052,  2055,  2058,  2062,  2072,  2073,
    2074,  2078,  2083,  2088,  2092,  2093,  2094,  2097,  2098,  2101,
    2102,  2103,  2104,  2105,  2106,  2107,  2108,  2111,  2112,  2119,
    2119,  2126,  2127,  2131,  2132,  2135,  2136,  2137,  2141,  2142,
    2152,  2155,  2159,  2165,  2166,  2181,  2182,  2183,  2187,  2193,
    2194,  2198,  2199,  2203,  2205,  2209,  2210,  2214,  2215,  2218,
    2224,  2225,  2242,  2247,  2248,  2252,  2258,  2259,  2276,  2280,
    2281,  2282,  2289,  2305,  2309,  2310,  2320,  2323,  2341,  2342,
    2351,  2355,  2359,  2380,  2381,  2382,  2383
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
  "PARTITION_Y", "PARTITION_Z", "PI_TOK", "POLYGON_LIST", "POSITIONS",
  "PRINTF", "PRINT_TIME", "PROBABILITY_REPORT",
  "PROBABILITY_REPORT_THRESHOLD", "PROGRESS_REPORT", "RADIAL_DIRECTIONS",
  "RADIAL_SUBDIVISIONS", "RAND_GAUSSIAN", "RAND_UNIFORM",
  "REACTION_DATA_OUTPUT", "REACTION_OUTPUT_REPORT", "REAL",
  "RECTANGULAR_RELEASE_SITE", "RECTANGULAR_TOKEN", "REFLECTIVE",
  "RELEASE_EVENT_REPORT", "RELEASE_INTERVAL", "RELEASE_PATTERN",
  "RELEASE_PROBABILITY", "RELEASE_SITE", "REMOVE_ELEMENTS", "RIGHT",
  "ROTATE", "ROUND_OFF", "SCALE", "SEED", "SHAPE", "SHOW_EXACT_TIME",
  "SIN", "SITE_DIAMETER", "SITE_RADIUS", "SPACE_STEP", "SPHERICAL",
  "SPHERICAL_RELEASE_SITE", "SPHERICAL_SHELL", "SPHERICAL_SHELL_SITE",
  "SPRINTF", "SQRT", "STANDARD_DEVIATION", "STEP", "STRING_TO_NUM",
  "STR_VALUE", "SUBUNIT", "SUBUNIT_RELATIONSHIPS", "SUMMATION_OPERATOR",
  "SURFACE_CLASS", "SURFACE_ONLY", "TAN", "TARGET_ONLY",
  "TET_ELEMENT_CONNECTIONS", "THROUGHPUT_REPORT", "TIME_LIST",
  "TIME_POINTS", "TIME_STEP", "TIME_STEP_MAX", "TO", "TOP",
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
     495,   496,   497,   498,    61,    38,    58,    43,    45,    42,
      47,   499,    94,    91,    93,    59,    39,    44,   123,   125,
      40,    41,    62,    60,    64
};
# endif

#define YYPACT_NINF -975

#define yypact_value_is_default(Yystate) \
  (!!((Yystate) == (-975)))

#define YYTABLE_NINF -391

#define yytable_value_is_error(Yytable_value) \
  0

  /* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
     STATE-NUM.  */
static const yytype_int16 yypact[] =
{
    2280,  -151,  -141,   -13,    28,    32,    72,   -78,  -134,    69,
     -78,   -78,    81,   116,   137,   129,   182,   188,   194,  -975,
     209,   242,   265,   269,   277,   279,   281,   285,   219,   273,
    -975,  -975,  -975,   275,   288,   293,   308,   297,   322,   323,
     324,   346,   361,  -975,   349,   351,   375,   311,  2280,  -975,
     -38,  -975,  -975,   395,  -975,  -975,   562,  -975,  -975,  -975,
    -975,  -975,  -975,  -975,  -975,  -975,  -975,  -975,  -975,   398,
    -975,  -975,  -975,  -975,  -975,  -975,  -975,  -975,  -975,  -975,
    -975,    22,  -975,  -975,  -975,  -975,   484,  -975,  -975,  -975,
    -975,  -975,  -975,  -975,  -975,   214,   214,   -42,  2108,   -42,
    2108,  -975,  -975,   386,   -78,   -78,  -975,   390,  -975,   391,
    -975,   -78,   -78,  2108,   -78,   -78,   -78,   -42,   -78,  2108,
    2108,   214,  2108,  2108,  2108,  2108,   333,   -78,   811,   -42,
     -42,  1208,  2108,   503,  2108,   -78,  2108,  2108,  2108,  -975,
     582,  1010,  -975,  -975,  1802,   399,   161,   404,  -975,  -975,
     404,  -975,   404,  -975,  -975,   404,   404,   404,  -975,  -975,
    -975,  -975,  -975,  -975,  -975,  -975,   405,  -975,  -975,  -975,
    -975,  -975,   418,  -975,  -975,   406,   411,   412,   414,   415,
     416,   417,   427,  -975,   429,   439,   440,   442,   445,  -975,
    -975,  -975,  -975,   447,  -975,   449,   451,   452,   453,  2108,
    2108,  2108,  -975,   220,  -975,  -975,  -975,  -975,  -975,  1115,
      11,   276,    76,  -975,  -975,   358,  -975,    92,  -975,  -975,
      -8,  -975,  -975,  -975,    99,  -975,  -975,  -975,   139,  -975,
    1230,  -975,   409,   424,   460,   418,  -975,   571,  -975,  1230,
    1230,  -975,  1230,  1230,  1230,  1230,  -975,  -975,  -975,   468,
     464,   143,  -975,   479,   480,   485,   489,   490,   496,   498,
     499,   500,   501,   504,   505,   508,   513,   520,   527,   528,
     529,   280,  -975,   418,  -975,   467,  -975,  1230,  1230,   530,
    -975,  1230,  -975,   493,  1230,  1230,  1230,   651,   534,   647,
     538,   543,   547,   550,   551,   554,   559,   573,   587,   590,
     591,   592,   594,   605,   606,   607,   410,  -975,    61,  1200,
    -975,  -975,  1230,  1236,  -975,  1247,   418,   561,   -42,  -975,
    -975,  -975,  -975,  -975,   771,   -78,  -975,   563,  -975,   563,
     -42,   -42,  2108,  2108,  2108,  2108,  2108,  2108,  2108,  2108,
    2108,  2108,  2108,  2108,  2108,  2108,  2108,  2108,   -42,  2108,
     574,   574,   379,  -975,  -975,  2108,  2108,  2108,  2108,  2108,
    -975,  2108,  -975,   610,   611,   374,  -975,  -975,  -975,  -975,
    -975,  2108,  -975,   205,  -975,  -975,  -975,  -975,  -975,   -78,
     -78,   -48,   125,  -975,  -975,  -975,   604,  -975,  -975,  -975,
     -42,   -42,   -78,  -975,  -975,  -975,   214,   214,   214,    15,
     214,   214,  1642,   214,   214,   214,  2108,   214,    15,   214,
     214,   214,    15,    15,  -975,  -975,  -163,  -975,  2108,   -57,
     -42,   614,    14,  -975,   -42,   619,    -1,  -975,   -30,   -30,
     -30,  2108,   -30,  2108,   -30,   -30,  2108,   -30,   -30,   -30,
     -30,   -30,   -30,   -30,  -975,  -975,  2108,  -103,  -975,  1230,
     609,   623,   711,  -975,   342,   -78,  -975,  -975,   687,   617,
     667,   421,   819,  -975,  -975,   506,   536,   552,   572,   595,
     678,   753,   782,   810,   825,  1065,  1097,  1121,  1134,   858,
     876,    51,   905,  -975,   183,   183,   574,   574,  -975,  1212,
    2108,  2108,   636,   637,   676,   721,  -975,  -975,  -975,  -975,
     358,  -975,  -975,   638,   -87,  -975,  -120,  -975,  -975,  -975,
     -50,   645,   646,   664,   668,   670,  -975,    36,   -78,  -975,
     630,   657,  -975,  -975,  -975,  -975,  -975,  -975,  -975,  -975,
    -975,  -975,  1230,  -975,  -975,  -975,  -975,  1230,  -975,  -975,
    -975,  -975,  -975,  -975,  -975,  1855,  -975,  1230,   685,   691,
     692,   -31,  -975,  -975,  -975,  -975,    49,   693,   662,   -58,
    -975,  -975,  -975,  -975,   418,   -78,   694,  -975,   702,  -975,
    -975,  -975,  -975,  -975,  -975,  1230,  -975,  1230,  -975,  -975,
    1230,  -975,  -975,  -975,  -975,  -975,  -975,  -975,   373,  -975,
      61,   -42,   161,  -145,   151,  -975,   699,   421,   161,   686,
    -975,   710,   712,   703,   716,   719,   706,   728,   733,   737,
    -975,  -975,   725,   421,  -975,   741,  -975,  -975,  -975,  -975,
    -975,   730,  -975,    93,  -975,  -975,  -975,  -975,  -975,  -975,
    -975,  -975,  -975,  -975,  2108,  2108,  2108,  2108,  -975,  -975,
    -975,  -975,  2108,  1230,  1230,  2108,  2108,  -975,   877,  -975,
    -975,   749,  -975,  -975,   638,  -975,   638,  -975,  -975,    -4,
    -975,  2108,  1955,  2108,  2108,  2108,  -975,   -78,   731,   740,
    -975,  -975,  -975,  -975,  -975,  -121,  -975,  -975,  -975,   750,
     112,  -975,  -975,   -72,  -975,  -975,   561,  -975,   161,  2108,
     161,   763,   764,  -975,   -51,  -975,  -975,  -975,  -975,   226,
    -975,  -975,  -975,   -42,     1,  -975,  -975,  -975,  -975,   766,
     161,   772,   773,  2108,  -975,   418,   754,   760,  -975,   404,
     775,   776,   779,  -975,  -975,  -975,  -975,    80,   421,  -975,
    -975,  -153,   161,  -975,  2108,  2108,   909,   161,   -78,   -78,
    2108,   -78,  2108,   911,   151,  -975,  2008,   161,  -975,  -975,
     927,   951,   990,  1005,  1224,  1230,  1230,   784,   774,    29,
    -975,  -975,   -50,  1061,   790,  -975,  -975,  1230,  -975,  1230,
    -975,  1230,  1230,  1230,   793,   -78,   -78,  -975,  -975,    38,
    -975,  -975,   794,  -975,  -975,  -975,  -975,  -975,  1230,  -975,
     566,   214,   181,  -975,  -975,  -975,   418,   783,   788,   791,
     -39,  -975,  -975,  -975,  -975,   -78,  -975,  2008,   796,     9,
     354,  -975,   161,  -975,   161,  2008,   161,  -975,  -975,  -975,
    -975,  -975,  -975,   -23,   468,  -975,   294,   151,  -975,  -975,
    -975,  -975,    96,   151,  1230,  1230,   804,  -975,  -975,   161,
     156,  -975,  1230,  -975,  -975,  1230,   808,  -975,  1258,  -975,
    -975,  -975,  -975,   134,  -975,     2,  -975,  -975,  -975,  -975,
    2108,  2108,  -975,  -975,  1855,  1855,  -975,  -975,   561,   140,
    -975,   -78,  -975,  2108,   358,   809,   180,  -975,   190,  -975,
     358,  -975,   797,   -78,  -975,  -975,   418,  -975,  -975,  -975,
     801,   816,  -975,   181,   181,  -975,   233,  -975,   865,  -975,
      53,  -975,    53,  -975,  -975,  -975,  1258,  -975,  -975,  -975,
    2008,   817,   824,   826,   821,  2108,  1050,  -975,   827,  -975,
    -975,   217,   -23,   -23,   -23,  -975,  -975,  -975,  -975,  2108,
    -975,  -975,  -975,  2108,  -975,  -975,   820,   839,   238,  -975,
    -975,  -975,  1230,  1230,  -975,  -975,  -975,  1061,  -975,  1230,
    -975,  2108,  -975,  -975,  -975,  -975,  -975,   491,  -975,   847,
    2108,   181,   849,  -975,   580,   181,   178,   -42,   181,   181,
     181,   181,  -975,  -975,  -975,  -975,    20,  -975,   840,    39,
      31,  -975,   845,  -975,   161,  2108,   161,  -975,  1049,   859,
    -975,   151,  2108,  -975,   873,   873,  -975,   253,   325,   -78,
    -975,  -975,   875,  1230,   883,  -975,  -975,   886,  -975,  -975,
     491,  -975,  -975,  -975,  -975,   891,  -975,   894,  -117,  1078,
     612,  -117,  -975,  -975,   878,   880,   881,   -42,   418,   278,
     278,  -975,  -975,  -975,  -975,    21,   901,  -975,  -975,  -975,
    -975,   901,  -975,  -975,     8,  -975,  1230,  -975,  -975,  2108,
    -975,  -975,  1230,   904,  -975,   907,   201,  -975,   902,  1332,
    -975,   896,   906,  -975,  -975,   -78,   161,   908,   903,   910,
     912,   914,  -975,  -975,  -975,  -975,  -975,   916,  -975,  -975,
     918,  -975,  -975,  -975,  -975,  -975,  2108,  -975,  -975,  -975,
    -975,  -975,  1230,     2,  2108,  2108,  -975,  -975,  -975,  -975,
    -975,  -975,  -975,  -975,  -975,  -975,   403,   919,  -975,   491,
    -975,   924,  -975,  1547,  1547,   -61,  -975,   926,    78,  -975,
      78,    78,  -975,  -975,  -975,  1230,  -975,   968,    33,   491,
    2108,  -975,  1547,   173,   231,  -975,   161,  -975,   468,  -975,
     928,   928,   928,   151,  -975,   920,   491,  1230,  -975,  -975,
    -975,  -975,   389,  -975,  -975,  -975,  -975,  2108,  -975,  -975,
    -975,  -975,  -975,  -975,  -975,  -975,  -975,  -975,  1032,   -54,
    -975,  -975
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
    -975,  -975,  -975,  1136,  -708,     0,   -93,  -109,  -126,  -572,
    -742,   -71,  -479,  -975,   814,   818,   121,  -975,   600,  -975,
    -975,  1057,  -132,  -102,  -136,  -975,   560,  -593,  -129,  -114,
    -975,   -99,   -95,  -139,  -975,  -975,  -975,  -975,  -975,  -975,
     435,   -98,  -195,  -975,  -975,  -975,  -975,  -975,  -975,  -975,
    -975,   922,   548,    67,  -975,  -975,   898,   653,  -975,   997,
    -975,  -975,  -975,  -975,  -975,  -975,  -975,  -975,  -975,   -18,
    -975,  -975,  -975,  -975,  -975,  -975,  -975,  -975,  -475,  -975,
    -975,  -975,  -975,   -28,  -975,   326,  -975,  -975,  -975,  -975,
    -975,  -975,   696,  -975,  -975,  -579,  -975,  -975,   991,  -313,
    -161,  -975,  -975,  -975,  -975,  -975,  -975,  -975,  -975,   831,
    -975,  -975,  -975,   469,  -975,  -975,  -975,   287,  -304,  -975,
    -975,  -975,  -975,  -975,  -975,  -975,  -975,  -279,   -91,  -138,
    -707,  -600,  -975,  -975,  1114,  -975,   786,  -975,  -975,  -975,
    -975,  -975,  -975,  -577,  -975,  -975,  -975,   649,  -975,  -550,
    -975,  -975,  -975,  -975,  -975,  -975,  -975,   420,  -975,  -975,
    -975,   915,   512,  -975,  -975,  -975,   408,   206,  -975,  -975,
    -975,  -975,  -975,  -967,  -953,  -975,  -975,  -975,  -543,   126,
    -975,  -975,  -975,  -975,  -975,  -975,   208,  -975,  -975,  -975,
    -975,  -975,   433,  -975,  -975,  -975,  -975,  -975,  -975,  -975,
    1039,  -975,  -975,  -975,   752,  -974,  -975,  -975,  -975,  -975,
    1017,  -975,  -975,  -975,  -975,  -975,  -975,  -975,  -975,  -975,
     577,  -975,  -975,  -975,  -975,  -975,  -975,   307,  -374,  -975,
    -975,  -975,  -975,  -975,  -975,  -975,   254,  -975,  -975,  -975,
    -542,  -505,  -975,  -975,  -975,  -975,  -975,  -975,  -975,   715,
    -975,  -975,  -975,  -975,   476,  -975,   236,  -975,  -975,  -975,
    -975,  -975,  -975,   303,  -975,  -975,  -975,   309,  -754,  -975,
    -975,  -975,   861,   487,  -975,  -975,  -975,  -975,  -975
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
     106,   108,   310,   744,   326,   313,   557,   327,   328,   329,
     717,   650,   235,  1089,   164,   165,   729,   237,   973,   654,
     314,   656,   275,  1083,   569,   658,   691,   847,   166,   973,
     666,  -119,   527,   548,   320,   315,   453,   147,    50,   316,
     241,   880,  1036,   167,   936,   248,   691,  1063,   668,   669,
     570,   973,   797,   745,   175,   148,   700,   176,   504,   668,
     669,  1145,   667,   919,    43,   558,   668,   669,   701,   177,
     353,   178,    43,   225,   149,   168,   213,   158,   557,   179,
    1115,   657,  1109,    95,   545,   777,   798,   170,   546,   170,
     308,   180,   797,    96,   101,   214,   828,    43,   566,   911,
      43,   108,   226,   323,   231,   231,   231,   170,   236,   227,
     925,  -390,   817,   778,   104,   250,   928,   226,   827,   170,
     170,   181,   425,   833,   782,   282,   798,   692,   806,   182,
      43,   818,   655,   549,   309,   782,   319,   558,   982,    43,
     354,   589,   782,   159,   590,   550,  1063,   692,   183,   729,
     168,   652,  1146,   160,   161,   837,   668,   669,   983,   150,
     974,   184,   185,   186,   799,   653,   453,    43,   745,   511,
     975,   974,   187,   937,  1169,    43,   188,   783,   819,   702,
     651,   975,   774,  1063,   367,   151,   385,   879,  1135,   895,
     502,   705,   152,   974,    43,  1171,   -59,   693,   793,   912,
     571,   670,   101,   975,   799,   503,   153,   214,   154,   189,
     903,   913,   889,   172,   108,   450,   982,   693,   226,   190,
     191,    97,   463,   192,   454,    43,   464,   823,   162,   379,
     916,   168,   918,   762,   920,   193,   921,   194,    43,   763,
     195,   226,   890,   481,   820,   938,   380,   671,   163,   196,
     155,   308,   197,   838,   839,    43,    43,   930,   623,   198,
     512,  1038,    98,   821,   308,   822,    99,   652,   183,  1033,
     355,   356,   357,   358,  1051,   359,   158,   658,    43,   353,
    1042,   653,   519,   520,   253,   672,   331,   953,  1037,   955,
     895,   895,   513,    43,   881,    43,   254,    43,   199,   200,
    1067,   142,   640,  1067,   446,   255,   100,   991,   170,    43,
    1137,   201,   555,    43,   680,   452,    43,   105,   528,   564,
     170,   170,   363,   364,   720,   366,   721,   528,   256,   111,
     823,   528,   528,   514,   515,   994,   995,   996,   170,   308,
     545,   377,   159,   192,   749,   926,   257,   258,   384,   354,
     699,    43,   160,   161,   454,   839,    43,   147,   895,   545,
      43,   722,   895,   781,   112,   895,   895,   895,   895,   214,
     499,   113,   259,    43,   516,   148,   891,   308,    43,   114,
     170,   170,   521,   934,   946,   950,  1159,   947,   387,  1160,
    1161,   956,   394,   892,   149,   158,   686,    43,    43,   260,
     723,  1162,  1163,   684,   308,   931,   313,    43,   492,   493,
     170,   720,  1024,   721,   170,  1025,   290,  1149,    43,   893,
    1132,   685,   357,   358,   261,   359,  1156,   162,   117,   952,
    1026,   894,   115,   262,   263,   264,   315,   291,   116,   954,
     316,   265,  1053,   119,  1164,   452,   266,   163,   722,   601,
    1096,   369,   370,   292,   922,   923,   924,   355,   356,   357,
     358,   159,   359,  1165,  1166,   539,   602,   127,   993,   543,
     544,   160,   161,   545,   719,  1150,   120,   795,  1132,   593,
     965,   267,   966,   760,  1117,   761,   715,  1001,   293,   294,
     355,   356,   357,   358,   214,   359,   214,   603,   268,   121,
     214,   269,  1054,   122,   270,   151,   295,   296,   678,   963,
     964,   123,   152,   124,  1055,   125,  1143,   970,   971,   126,
     604,   128,   297,   298,   299,   129,   153,   131,   154,   414,
     246,   922,   923,   924,   300,   309,   301,   302,   130,   605,
     233,   234,   132,   606,   915,   133,   162,  1004,  1005,  1006,
     944,   945,   303,   304,  1151,   707,   134,   607,   136,   247,
    1133,  1134,   355,   356,   357,   358,   163,   359,  1141,  1142,
     155,   813,  1007,   135,  1008,  1009,   787,  1020,   789,   713,
     137,   170,   319,   226,  1029,  1030,  1031,  1032,   319,   608,
     609,   355,   356,   357,   358,   138,   359,   139,   796,   140,
     849,   610,   611,   368,   369,   370,   371,   313,   824,  1130,
     355,   356,   357,   358,   868,   359,   355,   356,   357,   358,
     830,   359,   850,   141,   313,   305,  1154,  1155,   158,   144,
     483,   145,   146,   157,   211,   852,   166,   315,   222,   223,
     355,   356,   357,   358,   279,   359,   612,   288,   203,   318,
     209,   167,   323,   331,   315,   330,   332,   214,   316,   444,
     389,   333,   334,   230,   335,   336,   337,   338,   313,   239,
     240,   390,   242,   243,   244,   245,   313,   339,   319,   340,
     319,   277,   278,   909,   281,   886,   284,   285,   286,   341,
     342,   909,   343,   170,   159,   344,   884,   345,   315,   346,
     319,   347,   348,   349,   160,   161,   315,   391,   150,   885,
     888,   392,   393,   396,   397,   868,   868,   226,   417,   398,
     830,   319,   319,   399,   400,   313,   313,   319,   214,   214,
     401,   843,   402,   403,   404,   405,   848,   319,   406,   407,
     420,   852,   408,   355,   356,   357,   358,   409,   359,   350,
     351,   352,   214,   309,   410,   315,   315,   624,   168,   316,
     316,   411,   412,   413,   418,   874,   874,   421,   424,   214,
     425,   313,   428,   355,   356,   357,   358,   429,   359,   162,
     170,   430,   459,    43,   431,   432,   909,   625,   433,   355,
     356,   357,   358,   434,   359,   707,  -104,   906,   868,   163,
     451,   315,   319,   626,   319,   906,   319,   435,   313,   355,
     356,   357,   358,   226,   359,   253,   359,   968,   969,   970,
     971,   436,   319,   627,   437,   438,   439,   254,   440,   319,
     214,  1022,   355,   356,   357,   358,   255,   359,   315,   441,
     442,   443,   316,   319,   490,   491,   628,   518,   556,   968,
     969,   970,   971,   565,   309,   309,   591,   592,   449,   256,
    -390,   948,  1028,  1076,   596,   598,   874,   599,   874,   621,
     645,   646,  1045,   499,  1047,   647,   502,   257,   258,   661,
     662,   681,   465,   466,   467,   468,   469,   470,   471,   472,
     473,   474,   475,   476,   477,   478,   479,   480,   663,   482,
     906,   682,   664,   259,   665,   484,   485,   486,   487,   488,
     704,   489,   226,   226,   226,   355,   356,   357,   358,   688,
     359,   495,   172,  1084,  1081,   689,   690,   703,   710,   629,
     260,   711,  1090,   727,   732,   525,   526,   309,   530,   531,
     533,   534,   535,   536,   734,   538,   735,   540,   541,   542,
     737,   736,   532,   738,   739,   261,   537,   170,   355,   356,
     357,   358,   740,   359,   262,   263,   264,   741,   547,  1039,
     649,   742,   265,   743,   319,   746,   319,   266,   747,   775,
     757,   575,  1139,   577,  1139,  1139,   580,   652,   776,   678,
     355,   356,   357,   358,   780,   359,   588,   790,   791,  1138,
     808,  1138,  1138,   805,   630,   811,   807,   812,  1068,   814,
     815,  1068,   267,   816,   836,   846,   290,   170,   861,   355,
     356,   357,   358,   862,   359,   319,   871,   873,   883,   268,
     910,   900,   269,   631,   319,   270,   901,   291,   929,   902,
     643,   644,   933,   951,   960,   957,   678,   355,   356,   357,
     358,   984,   359,   292,   175,  1116,   319,   176,   985,   989,
     986,   632,   355,   356,   357,   358,   961,   359,   999,   177,
     987,   178,   573,   574,   992,   576,   633,   578,   579,   179,
     581,   582,   583,   584,   585,   586,   587,  1000,   293,   294,
    1018,   180,  1021,  1049,  1035,   355,   356,   357,   358,  1044,
     359,   967,   968,   969,   970,   971,   295,   296,   226,   638,
     226,   226,   924,   355,   356,   357,   358,  1061,   359,  1060,
    1062,   181,   297,   298,   299,  1065,   319,   639,  1066,   182,
    1078,   166,  1079,  1080,   300,  1086,   301,   302,  1094,  1113,
     449,  1095,   355,   356,   357,   358,   167,   359,   183,  1114,
    1098,   371,   303,   304,  1157,  1118,   641,  1120,  1122,  1121,
    1123,   184,   185,   186,   355,   356,   357,   358,  1131,   359,
    1124,  1132,   187,  1136,   143,  1152,   188,   496,   856,  1119,
     714,   497,   283,   415,   750,   751,   752,   753,   355,   356,
     357,   358,   754,   359,   445,   755,   756,   362,   378,   958,
     498,   175,   857,   677,   176,   355,   356,   357,   358,   189,
     359,   767,   769,   771,   772,   773,   177,  1144,   178,   190,
     191,   863,   238,   192,  1002,   305,   179,   355,   356,   357,
     358,   595,   359,   829,   462,   193,   728,   194,   180,   788,
     195,   858,   355,   356,   357,   358,  1093,   359,  1148,   196,
     932,   935,   197,   168,  1097,   927,   859,   388,   395,   198,
     679,   794,  1023,   810,   706,  1077,   904,  1088,   181,   355,
     356,   357,   358,  1043,   359,  1034,   182,   567,    43,     0,
       0,  1170,   905,     0,   834,   835,   355,   356,   357,   358,
     842,   359,   845,  1048,     0,   183,   276,     0,   199,   200,
       0,     0,   355,   356,   357,   358,     0,   359,   184,   185,
     186,   201,   634,   864,   865,   355,   356,   357,   358,   187,
     359,     0,  1075,   188,     0,   175,     0,     0,   176,     0,
       0,  1099,     0,     0,   355,   356,   357,   358,     0,   359,
     177,     0,   178,  1100,   635,     0,  1101,     0,     0,     0,
     179,   361,   355,   356,   357,   358,   189,   359,   355,   356,
     357,   358,   180,   359,     0,     0,   190,   191,   636,     0,
     192,   355,   356,   357,   358,     0,   359,     0,     0,     0,
       0,   637,   193,     0,   194,     0,     0,   195,  1004,  1005,
    1006,     0,   181,     0,     0,     0,   196,     0,     0,   197,
     182,     0,     0,     0,     0,  1102,   198,     0,     0,     0,
     942,   943,     0,  1007,     0,  1008,  1009,     0,     0,   183,
       0,     0,     0,   949,     0,    43,     0,  1103,     0,     0,
       0,     0,   184,   185,   186,  -110,     0,   -73,   -73,   -73,
     -73,     0,   -73,   187,     0,   199,   200,   188,   642,   355,
     356,   357,   358,     0,   359,     0,     0,     0,   201,     0,
     860,   355,   356,   357,   358,   988,   359,   355,   356,   357,
     358,     0,   359,   -67,   -67,   -67,   -67,     0,   -67,   997,
     189,     0,     0,   998,   -66,   -66,   -66,   -66,     0,   -66,
     190,   191,     0,     0,   192,   -73,   -73,   -73,   -73,     0,
     -73,  1003,     0,     0,  1104,     0,   193,     0,   194,     0,
    1019,   195,     0,     0,     0,     0,     0,     0,     0,     0,
     196,     0,     0,   197,     0,     0,     0,     0,     0,     0,
     198,     0,     0,     0,     0,  1046,     0,     0,     0,  1105,
     175,     0,  1052,   176,     0,     0,  1099,     0,     0,    43,
       0,     0,     0,     0,     0,   177,     0,   178,  1100,     0,
       0,  1101,     0,     0,     0,   179,     0,     0,     0,   199,
     200,     0,     0,     0,     0,     0,     0,   180,     0,     0,
       0,     0,   201,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,  1092,
       0,     0,     0,     0,     0,     0,     0,   181,     0,  1106,
       0,     0,     0,     0,     0,   182,     0,     0,     0,     0,
    1102,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,   183,   175,  1125,     0,   176,     0,
       0,     0,  1103,     0,  1127,  1128,     0,   184,   185,   186,
     177,     0,   178,     0,     0,     0,     0,     0,   187,     0,
     179,     0,   188,  1106,  1106,     0,     0,     0,     0,     0,
       0,     0,   180,     0,     0,     0,     0,     0,     0,     0,
    1147,     0,  1106,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,   189,     0,     0,     0,     0,
       0,     0,   181,     0,   158,   190,   191,  1168,     0,   192,
     182,     0,     0,     0,     0,     0,     0,     0,     0,  1104,
       0,   193,     0,   194,     0,     0,   195,     0,     0,   183,
       0,     0,     0,     0,     0,   196,     0,     0,   197,     0,
       0,     0,   184,   185,   186,   198,     0,     0,     0,     0,
       0,     0,     0,   187,  1105,     0,     0,   188,     0,     0,
       0,     0,     0,     0,    43,     0,     0,     0,     0,     0,
     159,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     160,   161,     0,     0,   199,   200,     0,     0,     0,     0,
     189,     0,     0,     0,     0,   175,     0,   201,   176,     0,
     190,   191,     0,     0,   192,     0,     0,     0,     0,     0,
     177,     0,   178,     0,     0,     0,   193,     0,   194,     0,
     179,   195,     0,     0,     0,     0,     0,     0,     0,     0,
     196,     0,   180,   197,     0,     0,     0,     0,     0,     0,
     198,     0,     0,     0,     0,     0,     0,     0,   175,     0,
       0,   176,     0,     0,     0,   162,     0,     0,     0,    43,
       0,     0,   181,   177,     0,   178,     0,     0,     0,     0,
     182,     0,   166,   179,     0,   163,     0,     0,     0,   199,
     200,     0,     0,     0,     0,   180,     0,   167,     0,   183,
       0,     0,   201,     0,     0,     0,     0,     0,     0,     0,
       0,     0,   184,   185,   186,     0,     0,     0,     0,     0,
       0,     0,     0,   187,     0,   181,     0,   188,     0,     0,
       0,     0,     0,   182,     0,   166,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     167,     0,   183,     0,     0,     0,     0,     0,   175,     0,
     189,   176,     0,     0,     0,   184,   185,   186,     0,     0,
     190,   191,     0,   177,   192,   178,   187,     0,     0,     0,
     188,     0,     0,   179,     0,     0,   193,     0,   194,     0,
       0,   195,     0,     0,     0,   180,     0,     0,     0,     0,
     196,     0,     0,   197,   168,     0,     0,     0,     0,     0,
     198,   175,     0,   189,   176,     0,     0,     0,     0,     0,
       0,     0,     0,   190,   191,   181,   177,   192,   178,    43,
       0,     0,     0,   182,     0,     0,   179,     0,     0,   193,
       0,   194,     0,     0,   195,     0,     0,     0,   180,   199,
     200,     0,   183,   196,     0,   308,   197,   168,     0,     0,
       0,     0,   201,   198,     0,   184,   185,   186,     0,     0,
       0,     0,     0,     0,     0,     0,   187,     0,   181,     0,
     188,     0,    43,     0,     0,     0,   182,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,   199,   200,     0,   183,     0,     0,     0,     0,
       0,   175,     0,   189,   176,   201,     0,     0,   184,   185,
     186,     0,     0,   190,   191,     0,   177,   192,   178,   187,
       0,     0,     0,   188,     0,     0,   179,     0,     0,   193,
       0,   194,     0,     0,   195,     0,     0,     0,   180,     0,
       0,     0,     0,   196,     0,     0,   197,     0,     0,     0,
       0,     0,     0,   198,     0,     0,   189,     0,     0,     0,
       0,     0,     0,     0,     0,     0,   190,   191,   181,   768,
     192,     0,    43,     0,     0,     0,   182,     0,     0,     0,
       0,     0,   193,     0,   194,     0,     0,   195,     0,     0,
       0,     0,   199,   200,     0,   183,   196,     0,     0,   197,
       0,     0,     0,     0,     0,   201,   198,     0,   184,   185,
     186,     0,     0,     0,     0,     0,     0,     0,     0,   187,
       0,     0,     0,   188,     0,    43,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,   199,   200,     0,     0,     0,
       0,   308,     0,     0,     0,     0,   189,     0,   201,     0,
       0,     0,     0,     0,     0,     0,   190,   191,     0,     0,
     192,     0,     0,     0,     0,     1,     0,     0,     0,     0,
       0,     0,   193,     0,   194,     0,     0,   195,     0,     0,
       0,     0,     0,     0,     0,     0,   196,     0,     0,   197,
       2,     3,     4,     5,     6,     0,   198,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     7,     8,     9,    10,
      11,    12,    13,     0,     0,    43,     0,     0,     0,    14,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,    15,     0,   199,   200,     0,     0,     0,
       0,    16,    17,     0,     0,     0,     0,     0,   201,     0,
       0,     0,    18,     0,     0,     0,    19,     0,     0,    20,
       0,     0,     0,    21,    22,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,    23,    24,    25,    26,
      27,     0,     0,     0,     0,     0,     0,    28,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,    29,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,    30,    31,    32,     0,     0,
       0,    33,    34,     0,     0,     0,    35,    36,     0,     0,
      37,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,    38,     0,     0,     0,     0,    39,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,    40,    41,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,    42,    43,     0,     0,
      44,     0,     0,    45,     0,     0,     0,     0,     0,     0,
       0,    46
};

static const yytype_int16 yycheck[] =
{
       0,   127,    97,   112,    99,   144,    99,     7,   144,   147,
      10,    11,   144,   613,   152,   144,    74,   155,   156,   157,
     592,   500,   117,    15,    95,    96,   598,   118,     8,   504,
     144,   506,   130,    12,    64,   510,    87,   744,    80,     8,
       4,    79,    27,   100,   146,   144,   325,    25,    48,   144,
     121,    13,    13,    95,    52,   126,    87,  1010,   130,   131,
      90,     8,   101,   613,     3,    43,    17,     6,   381,   130,
     131,    38,    36,   815,   227,   133,   130,   131,    29,    18,
      69,    20,   227,   111,    62,   202,   104,    72,    74,    28,
    1064,   141,  1059,   244,   257,   216,   135,    97,   261,    99,
     253,    40,   101,   244,   104,   105,   259,   227,   109,   100,
     227,   111,   112,   258,   114,   115,   116,   117,   118,   228,
     827,   159,    42,   244,   258,   251,   833,   127,   728,   129,
     130,    70,   133,   733,   206,   135,   135,   188,   710,    78,
     227,    61,   262,   200,   144,   206,   146,   133,   902,   227,
     139,   254,   206,   138,   257,   212,  1109,   188,    97,   731,
     202,   248,  1129,   148,   149,   737,   130,   131,   910,   147,
     150,   110,   111,   112,   213,   262,   455,   227,   728,    54,
     160,   150,   121,   181,  1158,   227,   125,   259,   108,   140,
     503,   160,   667,  1146,   212,   173,   224,   776,   259,   792,
     248,   259,   180,   150,   227,   259,   244,   258,   259,   200,
     240,   175,   212,   160,   213,   263,   194,   217,   196,   158,
     259,   212,    41,   318,   224,   318,   980,   258,   228,   168,
     169,   244,   330,   172,   325,   227,   331,   260,   223,   247,
     812,   202,   814,   247,   816,   184,   823,   186,   227,   253,
     189,   251,    71,   348,   174,   855,   264,   221,   243,   198,
     238,   253,   201,   738,   739,   227,   227,   839,   463,   208,
     145,   979,   244,   193,   253,   195,   244,   248,    97,   259,
     247,   248,   249,   250,   991,   252,    72,   762,   227,    69,
     259,   262,   390,   391,    14,   259,   245,   876,   259,   878,
     893,   894,   177,   227,   779,   227,    26,   227,   247,   248,
    1018,     0,   261,  1021,   253,    35,   244,   917,   318,   227,
     242,   260,   420,   227,   519,   325,   227,   258,   399,   424,
     330,   331,    56,    57,   183,   259,   185,   408,    58,   258,
     260,   412,   413,   218,   219,   922,   923,   924,   348,   253,
     257,   259,   138,   172,   261,   259,    76,    77,   259,   139,
     555,   227,   148,   149,   455,   840,   227,    25,   961,   257,
     227,   220,   965,   261,   258,   968,   969,   970,   971,   379,
     380,   244,   102,   227,   259,    43,   205,   253,   227,   260,
     390,   391,   392,   259,   254,   874,     7,   257,   259,    10,
      11,   880,   259,   222,    62,    72,   545,   227,   227,   129,
     259,    22,    23,   545,   253,   259,   545,   227,    44,    45,
     420,   183,   244,   185,   424,   247,    16,   254,   227,   248,
     257,   545,   249,   250,   154,   252,  1143,   223,   244,   259,
     262,   260,   260,   163,   164,   165,   545,    37,   260,   259,
     545,   171,   199,   244,    65,   455,   176,   243,   220,    38,
     259,   256,   257,    53,   247,   248,   249,   247,   248,   249,
     250,   138,   252,    84,    85,   408,    55,   258,   261,   412,
     413,   148,   149,   257,   593,   254,   244,   261,   257,   147,
     257,   211,   259,   654,  1066,   656,   591,   259,    88,    89,
     247,   248,   249,   250,   504,   252,   506,    86,   228,   244,
     510,   231,   259,   244,   234,   173,   106,   107,   518,   893,
     894,   244,   180,   244,   199,   244,  1126,   249,   250,   244,
     109,   258,   122,   123,   124,   260,   194,   244,   196,   259,
     207,   247,   248,   249,   134,   545,   136,   137,   260,   128,
     115,   116,   244,   132,   200,   258,   223,    66,    67,    68,
     864,   865,   152,   153,  1136,   565,   244,   146,   244,   236,
    1113,  1114,   247,   248,   249,   250,   243,   252,  1120,  1121,
     238,   719,    91,   260,    93,    94,   688,   961,   690,   216,
     244,   591,   592,   593,   968,   969,   970,   971,   598,   178,
     179,   247,   248,   249,   250,   244,   252,   258,   703,   258,
     746,   190,   191,   255,   256,   257,   258,   746,   727,   216,
     247,   248,   249,   250,   763,   252,   247,   248,   249,   250,
     732,   252,   746,   258,   763,   225,  1141,  1142,    72,   244,
     261,    79,   244,   159,   258,   747,    80,   746,   258,   258,
     247,   248,   249,   250,   151,   252,   235,    75,    98,   260,
     100,    95,   258,   245,   763,   260,   260,   667,   763,   259,
     261,   260,   260,   113,   260,   260,   260,   260,   807,   119,
     120,   257,   122,   123,   124,   125,   815,   260,   688,   260,
     690,   131,   132,   807,   134,   790,   136,   137,   138,   260,
     260,   815,   260,   703,   138,   260,   140,   260,   807,   260,
     710,   260,   260,   260,   148,   149,   815,   257,   147,   790,
     791,   253,   258,   244,   244,   864,   865,   727,   261,   244,
     832,   731,   732,   244,   244,   864,   865,   737,   738,   739,
     244,   741,   244,   244,   244,   244,   746,   747,   244,   244,
     257,   853,   244,   247,   248,   249,   250,   244,   252,   199,
     200,   201,   762,   763,   244,   864,   865,   261,   202,   864,
     865,   244,   244,   244,   244,   775,   776,   126,   244,   779,
     133,   910,   244,   247,   248,   249,   250,   244,   252,   223,
     790,   244,   229,   227,   244,   244,   910,   261,   244,   247,
     248,   249,   250,   244,   252,   805,   245,   807,   947,   243,
      39,   910,   812,   261,   814,   815,   816,   244,   947,   247,
     248,   249,   250,   823,   252,    14,   252,   247,   248,   249,
     250,   244,   832,   261,   244,   244,   244,    26,   244,   839,
     840,   261,   247,   248,   249,   250,    35,   252,   947,   244,
     244,   244,   947,   853,   244,   244,   261,   253,   244,   247,
     248,   249,   250,   244,   864,   865,   257,   244,   308,    58,
     159,   871,   967,   261,   187,   258,   876,   210,   878,    60,
     244,   244,   984,   883,   986,   209,   248,    76,    77,   244,
     244,   261,   332,   333,   334,   335,   336,   337,   338,   339,
     340,   341,   342,   343,   344,   345,   346,   347,   244,   349,
     910,   254,   244,   102,   244,   355,   356,   357,   358,   359,
     258,   361,   922,   923,   924,   247,   248,   249,   250,   244,
     252,   371,  1027,  1035,  1027,   244,   244,   244,   244,   261,
     129,   239,  1044,   244,   258,   397,   398,   947,   400,   401,
     402,   403,   404,   405,   244,   407,   244,   409,   410,   411,
     244,   258,   402,   244,   258,   154,   406,   967,   247,   248,
     249,   250,   244,   252,   163,   164,   165,   244,   418,   979,
     259,   244,   171,   258,   984,   244,   986,   176,   258,   258,
     113,   431,  1118,   433,  1120,  1121,   436,   248,   258,   999,
     247,   248,   249,   250,   254,   252,   446,   244,   244,  1118,
     237,  1120,  1121,   247,   261,   261,   244,   257,  1018,   244,
     244,  1021,   211,   244,   115,   114,    16,  1027,   244,   247,
     248,   249,   250,   259,   252,  1035,   246,   244,   244,   228,
     244,   258,   231,   261,  1044,   234,   258,    37,   244,   258,
     490,   491,   244,   244,   253,   258,  1056,   247,   248,   249,
     250,   244,   252,    53,     3,  1065,  1066,     6,   244,    19,
     244,   261,   247,   248,   249,   250,   260,   252,   258,    18,
     259,    20,   429,   430,   257,   432,   261,   434,   435,    28,
     437,   438,   439,   440,   441,   442,   443,   258,    88,    89,
     253,    40,   253,   244,   264,   247,   248,   249,   250,   264,
     252,   246,   247,   248,   249,   250,   106,   107,  1118,   261,
    1120,  1121,   249,   247,   248,   249,   250,   244,   252,   254,
     244,    70,   122,   123,   124,   244,  1136,   261,   244,    78,
     262,    80,   262,   262,   134,   244,   136,   137,   244,   253,
     590,   244,   247,   248,   249,   250,    95,   252,    97,   253,
     258,   258,   152,   153,   244,   257,   261,   257,   254,   257,
     254,   110,   111,   112,   247,   248,   249,   250,   259,   252,
     262,   257,   121,   257,    48,   257,   125,   373,   261,  1068,
     590,   373,   135,   271,   634,   635,   636,   637,   247,   248,
     249,   250,   642,   252,   306,   645,   646,   210,   217,   883,
     379,     3,   261,   517,     6,   247,   248,   249,   250,   158,
     252,   661,   662,   663,   664,   665,    18,   259,    20,   168,
     169,   762,   118,   172,   947,   225,    28,   247,   248,   249,
     250,   455,   252,   731,   329,   184,   597,   186,    40,   689,
     189,   261,   247,   248,   249,   250,  1050,   252,  1132,   198,
     840,   853,   201,   202,  1056,   832,   261,   228,   251,   208,
     518,   694,   965,   713,   559,  1021,   800,  1041,    70,   247,
     248,   249,   250,   980,   252,   976,    78,   426,   227,    -1,
      -1,   259,   805,    -1,   734,   735,   247,   248,   249,   250,
     740,   252,   742,   254,    -1,    97,    98,    -1,   247,   248,
      -1,    -1,   247,   248,   249,   250,    -1,   252,   110,   111,
     112,   260,   257,   262,   263,   247,   248,   249,   250,   121,
     252,    -1,   254,   125,    -1,     3,    -1,    -1,     6,    -1,
      -1,     9,    -1,    -1,   247,   248,   249,   250,    -1,   252,
      18,    -1,    20,    21,   257,    -1,    24,    -1,    -1,    -1,
      28,   246,   247,   248,   249,   250,   158,   252,   247,   248,
     249,   250,    40,   252,    -1,    -1,   168,   169,   257,    -1,
     172,   247,   248,   249,   250,    -1,   252,    -1,    -1,    -1,
      -1,   257,   184,    -1,   186,    -1,    -1,   189,    66,    67,
      68,    -1,    70,    -1,    -1,    -1,   198,    -1,    -1,   201,
      78,    -1,    -1,    -1,    -1,    83,   208,    -1,    -1,    -1,
     860,   861,    -1,    91,    -1,    93,    94,    -1,    -1,    97,
      -1,    -1,    -1,   873,    -1,   227,    -1,   105,    -1,    -1,
      -1,    -1,   110,   111,   112,   245,    -1,   247,   248,   249,
     250,    -1,   252,   121,    -1,   247,   248,   125,   246,   247,
     248,   249,   250,    -1,   252,    -1,    -1,    -1,   260,    -1,
     246,   247,   248,   249,   250,   915,   252,   247,   248,   249,
     250,    -1,   252,   247,   248,   249,   250,    -1,   252,   929,
     158,    -1,    -1,   933,   247,   248,   249,   250,    -1,   252,
     168,   169,    -1,    -1,   172,   247,   248,   249,   250,    -1,
     252,   951,    -1,    -1,   182,    -1,   184,    -1,   186,    -1,
     960,   189,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
     198,    -1,    -1,   201,    -1,    -1,    -1,    -1,    -1,    -1,
     208,    -1,    -1,    -1,    -1,   985,    -1,    -1,    -1,   217,
       3,    -1,   992,     6,    -1,    -1,     9,    -1,    -1,   227,
      -1,    -1,    -1,    -1,    -1,    18,    -1,    20,    21,    -1,
      -1,    24,    -1,    -1,    -1,    28,    -1,    -1,    -1,   247,
     248,    -1,    -1,    -1,    -1,    -1,    -1,    40,    -1,    -1,
      -1,    -1,   260,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,  1049,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    70,    -1,  1059,
      -1,    -1,    -1,    -1,    -1,    78,    -1,    -1,    -1,    -1,
      83,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    97,     3,  1086,    -1,     6,    -1,
      -1,    -1,   105,    -1,  1094,  1095,    -1,   110,   111,   112,
      18,    -1,    20,    -1,    -1,    -1,    -1,    -1,   121,    -1,
      28,    -1,   125,  1113,  1114,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    40,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
    1130,    -1,  1132,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,   158,    -1,    -1,    -1,    -1,
      -1,    -1,    70,    -1,    72,   168,   169,  1157,    -1,   172,
      78,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   182,
      -1,   184,    -1,   186,    -1,    -1,   189,    -1,    -1,    97,
      -1,    -1,    -1,    -1,    -1,   198,    -1,    -1,   201,    -1,
      -1,    -1,   110,   111,   112,   208,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,   121,   217,    -1,    -1,   125,    -1,    -1,
      -1,    -1,    -1,    -1,   227,    -1,    -1,    -1,    -1,    -1,
     138,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
     148,   149,    -1,    -1,   247,   248,    -1,    -1,    -1,    -1,
     158,    -1,    -1,    -1,    -1,     3,    -1,   260,     6,    -1,
     168,   169,    -1,    -1,   172,    -1,    -1,    -1,    -1,    -1,
      18,    -1,    20,    -1,    -1,    -1,   184,    -1,   186,    -1,
      28,   189,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
     198,    -1,    40,   201,    -1,    -1,    -1,    -1,    -1,    -1,
     208,    -1,    -1,    -1,    -1,    -1,    -1,    -1,     3,    -1,
      -1,     6,    -1,    -1,    -1,   223,    -1,    -1,    -1,   227,
      -1,    -1,    70,    18,    -1,    20,    -1,    -1,    -1,    -1,
      78,    -1,    80,    28,    -1,   243,    -1,    -1,    -1,   247,
     248,    -1,    -1,    -1,    -1,    40,    -1,    95,    -1,    97,
      -1,    -1,   260,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,   110,   111,   112,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,   121,    -1,    70,    -1,   125,    -1,    -1,
      -1,    -1,    -1,    78,    -1,    80,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      95,    -1,    97,    -1,    -1,    -1,    -1,    -1,     3,    -1,
     158,     6,    -1,    -1,    -1,   110,   111,   112,    -1,    -1,
     168,   169,    -1,    18,   172,    20,   121,    -1,    -1,    -1,
     125,    -1,    -1,    28,    -1,    -1,   184,    -1,   186,    -1,
      -1,   189,    -1,    -1,    -1,    40,    -1,    -1,    -1,    -1,
     198,    -1,    -1,   201,   202,    -1,    -1,    -1,    -1,    -1,
     208,     3,    -1,   158,     6,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,   168,   169,    70,    18,   172,    20,   227,
      -1,    -1,    -1,    78,    -1,    -1,    28,    -1,    -1,   184,
      -1,   186,    -1,    -1,   189,    -1,    -1,    -1,    40,   247,
     248,    -1,    97,   198,    -1,   253,   201,   202,    -1,    -1,
      -1,    -1,   260,   208,    -1,   110,   111,   112,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,   121,    -1,    70,    -1,
     125,    -1,   227,    -1,    -1,    -1,    78,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,   247,   248,    -1,    97,    -1,    -1,    -1,    -1,
      -1,     3,    -1,   158,     6,   260,    -1,    -1,   110,   111,
     112,    -1,    -1,   168,   169,    -1,    18,   172,    20,   121,
      -1,    -1,    -1,   125,    -1,    -1,    28,    -1,    -1,   184,
      -1,   186,    -1,    -1,   189,    -1,    -1,    -1,    40,    -1,
      -1,    -1,    -1,   198,    -1,    -1,   201,    -1,    -1,    -1,
      -1,    -1,    -1,   208,    -1,    -1,   158,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,   168,   169,    70,   224,
     172,    -1,   227,    -1,    -1,    -1,    78,    -1,    -1,    -1,
      -1,    -1,   184,    -1,   186,    -1,    -1,   189,    -1,    -1,
      -1,    -1,   247,   248,    -1,    97,   198,    -1,    -1,   201,
      -1,    -1,    -1,    -1,    -1,   260,   208,    -1,   110,   111,
     112,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   121,
      -1,    -1,    -1,   125,    -1,   227,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,   247,   248,    -1,    -1,    -1,
      -1,   253,    -1,    -1,    -1,    -1,   158,    -1,   260,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,   168,   169,    -1,    -1,
     172,    -1,    -1,    -1,    -1,     5,    -1,    -1,    -1,    -1,
      -1,    -1,   184,    -1,   186,    -1,    -1,   189,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,   198,    -1,    -1,   201,
      30,    31,    32,    33,    34,    -1,   208,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    46,    47,    48,    49,
      50,    51,    52,    -1,    -1,   227,    -1,    -1,    -1,    59,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    73,    -1,   247,   248,    -1,    -1,    -1,
      -1,    81,    82,    -1,    -1,    -1,    -1,    -1,   260,    -1,
      -1,    -1,    92,    -1,    -1,    -1,    96,    -1,    -1,    99,
      -1,    -1,    -1,   103,   104,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,   116,   117,   118,   119,
     120,    -1,    -1,    -1,    -1,    -1,    -1,   127,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,   143,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,   155,   156,   157,    -1,    -1,
      -1,   161,   162,    -1,    -1,    -1,   166,   167,    -1,    -1,
     170,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,   192,    -1,    -1,    -1,    -1,   197,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,   214,   215,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,   226,   227,    -1,    -1,
     230,    -1,    -1,   233,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,   241
};

  /* YYSTOS[STATE-NUM] -- The (internal number of the) accessing
     symbol of state STATE-NUM.  */
static const yytype_uint16 yystos[] =
{
       0,     5,    30,    31,    32,    33,    34,    46,    47,    48,
      49,    50,    51,    52,    59,    73,    81,    82,    92,    96,
      99,   103,   104,   116,   117,   118,   119,   120,   127,   143,
     155,   156,   157,   161,   162,   166,   167,   170,   192,   197,
     214,   215,   226,   227,   230,   233,   241,   266,   267,   268,
     270,   284,   285,   286,   300,   301,   302,   304,   309,   310,
     311,   312,   313,   314,   319,   323,   326,   327,   328,   329,
     330,   331,   332,   344,   345,   346,   361,   384,   389,   391,
     392,   393,   399,   404,   405,   409,   423,   424,   454,   459,
     463,   473,   477,   509,   534,   244,   244,   244,   244,   244,
     244,   270,   334,   335,   258,   258,   270,   385,   270,   336,
     348,   258,   258,   244,   260,   260,   260,   244,   390,   244,
     244,   244,   244,   244,   244,   244,   244,   258,   258,   260,
     260,   244,   244,   258,   244,   260,   244,   244,   244,   258,
     258,   258,     0,   268,   244,    79,   244,    25,    43,    62,
     147,   173,   180,   194,   196,   238,   411,   159,    72,   138,
     148,   149,   223,   243,   276,   276,    80,    95,   202,   269,
     270,   271,   297,   298,   299,     3,     6,    18,    20,    28,
      40,    70,    78,    97,   110,   111,   112,   121,   125,   158,
     168,   169,   172,   184,   186,   189,   198,   201,   208,   247,
     248,   260,   270,   291,   292,   293,   295,   296,   271,   291,
     325,   258,   333,   334,   270,   341,   343,   362,   363,   372,
     373,   374,   258,   258,   347,   348,   270,   272,   464,   465,
     291,   270,   305,   305,   305,   297,   270,   393,   399,   291,
     291,   276,   291,   291,   291,   291,   207,   236,   276,   272,
     273,   474,   475,    14,    26,    35,    58,    76,    77,   102,
     129,   154,   163,   164,   165,   171,   176,   211,   228,   231,
     234,   315,   316,   297,   306,   306,    98,   291,   291,   151,
     479,   291,   270,   286,   291,   291,   291,   510,    75,   535,
      16,    37,    53,    88,    89,   106,   107,   122,   123,   124,
     134,   136,   137,   152,   153,   225,   320,   321,   253,   270,
     287,   289,   291,   293,   294,   296,   297,   298,   260,   270,
     288,   289,   290,   258,   394,   394,   394,   394,   394,   394,
     260,   245,   260,   260,   260,   260,   260,   260,   260,   260,
     260,   260,   260,   260,   260,   260,   260,   260,   260,   260,
     291,   291,   291,    69,   139,   247,   248,   249,   250,   252,
     324,   246,   324,    56,    57,   337,   259,   334,   255,   256,
     257,   258,   277,   278,   279,   280,   281,   259,   363,   247,
     264,   375,   387,   349,   259,   348,   466,   259,   465,   261,
     257,   257,   253,   258,   259,   475,   244,   244,   244,   244,
     244,   244,   244,   244,   244,   244,   244,   244,   244,   244,
     244,   244,   244,   244,   259,   316,   307,   261,   244,   478,
     257,   126,   512,   513,   244,   133,   536,   537,   244,   244,
     244,   244,   244,   244,   244,   244,   244,   244,   244,   244,
     244,   244,   244,   244,   259,   321,   253,   282,   283,   291,
     271,    39,   270,   392,   393,   400,   401,   402,   406,   229,
     426,   410,   426,   306,   297,   291,   291,   291,   291,   291,
     291,   291,   291,   291,   291,   291,   291,   291,   291,   291,
     291,   297,   291,   261,   291,   291,   291,   291,   291,   291,
     244,   244,    44,    45,   338,   291,   279,   280,   374,   270,
     350,   376,   248,   263,   364,   365,   366,   367,   368,   369,
     370,    54,   145,   177,   218,   219,   259,   351,   253,   306,
     306,   270,   476,   276,   317,   317,   317,    27,   276,   318,
     317,   317,   291,   317,   317,   317,   317,   291,   317,   318,
     317,   317,   317,   318,   318,   257,   261,   291,   100,   200,
     212,   480,   481,   482,   483,   306,   244,    74,   133,   511,
     514,   515,   516,   517,   297,   244,   109,   537,   540,    64,
      90,   240,   322,   322,   322,   291,   322,   291,   322,   322,
     291,   322,   322,   322,   322,   322,   322,   322,   291,   254,
     257,   257,   244,   147,   396,   401,   187,   407,   258,   210,
     456,    38,    55,    86,   109,   128,   132,   146,   178,   179,
     190,   191,   235,   412,   414,   415,   416,   417,   418,   419,
     420,    60,   429,   307,   261,   261,   261,   261,   261,   261,
     261,   261,   261,   261,   257,   257,   257,   257,   261,   261,
     261,   261,   246,   291,   291,   244,   244,   209,   339,   259,
     277,   364,   248,   262,   343,   262,   343,   141,   343,   377,
     378,   244,   244,   244,   244,   244,     4,    36,   130,   131,
     175,   221,   259,   352,   353,   354,   356,   357,   270,   469,
     307,   261,   254,   470,   287,   294,   298,   308,   244,   244,
     244,    87,   188,   258,   484,   485,   486,   489,   490,   307,
      17,    29,   140,   244,   258,   259,   514,   270,   538,   539,
     244,   239,   541,   216,   283,   297,   303,   274,   288,   272,
     183,   185,   220,   259,   395,   397,   398,   244,   412,   274,
     427,   428,   258,   455,   244,   244,   258,   244,   244,   258,
     244,   244,   244,   258,   396,   414,   244,   258,   425,   261,
     291,   291,   291,   291,   291,   291,   291,   113,   340,   364,
     365,   365,   247,   253,   379,   380,   381,   291,   224,   291,
     388,   291,   291,   291,   343,   258,   258,   216,   244,   355,
     254,   261,   206,   259,   357,   471,   472,   288,   291,   288,
     244,   244,   487,   259,   485,   261,   297,   101,   135,   213,
     518,   519,   520,   526,   530,   247,   274,   244,   237,   542,
     291,   261,   257,   394,   244,   244,   244,    42,    61,   108,
     174,   193,   195,   260,   272,   273,   408,   396,   259,   427,
     288,   457,   458,   396,   291,   291,   115,   274,   343,   343,
     421,   422,   291,   270,   386,   291,   114,   395,   270,   289,
     294,   413,   288,   430,   431,   432,   261,   261,   261,   261,
     246,   244,   259,   378,   262,   263,   287,   294,   298,   382,
     383,   246,   371,   244,   270,   342,   358,   360,   359,   360,
      13,   343,   467,   244,   140,   276,   297,   488,   276,    41,
      71,   205,   222,   248,   260,   292,   491,   492,   493,   494,
     258,   258,   258,   259,   519,   538,   270,   274,   275,   294,
     244,   100,   200,   212,   543,   200,   274,   403,   274,   275,
     274,   408,   247,   248,   249,   395,   259,   457,   395,   244,
     274,   259,   422,   244,   259,   431,    52,   181,   396,   433,
     434,   449,   291,   291,   383,   383,   254,   257,   270,   291,
     277,   244,   259,   360,   259,   360,   277,   258,   350,   495,
     253,   260,   496,   493,   493,   257,   259,   246,   247,   248,
     249,   250,   508,     8,   150,   160,   531,   532,   533,   522,
     527,   528,   533,   275,   244,   244,   244,   259,   291,    19,
     462,   396,   257,   261,   408,   408,   408,   291,   291,   258,
     258,   259,   382,   291,    66,    67,    68,    91,    93,    94,
     438,   439,   440,   441,   445,   446,   447,   448,   253,   291,
     493,   253,   261,   492,   244,   247,   262,   497,   297,   493,
     493,   493,   493,   259,   532,   264,    13,   259,   269,   270,
     523,   524,   259,   528,   264,   288,   291,   288,   254,   244,
     460,   395,   291,   199,   259,   199,   450,   451,   469,   435,
     254,   244,   244,   439,   468,   244,   244,   269,   270,   499,
     500,   501,   502,   503,   504,   254,   261,   501,   262,   262,
     262,   271,   498,    12,   288,   529,   244,   521,   521,    15,
     288,   525,   291,   432,   244,   244,   259,   451,   258,     9,
      21,    24,    83,   105,   182,   217,   291,   436,   437,   438,
     442,   443,   444,   253,   253,   470,   270,   274,   257,   281,
     257,   257,   254,   254,   262,   291,   461,   291,   291,   452,
     216,   259,   257,   443,   443,   259,   257,   242,   272,   273,
     505,   505,   505,   396,   259,    38,   438,   291,   444,   254,
     254,   274,   257,   506,   506,   506,   395,   244,   453,     7,
      10,    11,    22,    23,    65,    84,    85,   507,   291,   470,
     259,   259
};

  /* YYR1[YYN] -- Symbol number of symbol that rule YYN derives.  */
static const yytype_uint16 yyr1[] =
{
       0,   265,   266,   267,   267,   268,   268,   268,   268,   268,
     268,   268,   268,   268,   268,   268,   268,   268,   268,   268,
     268,   268,   268,   268,   268,   269,   270,   271,   272,   273,
     274,   275,   275,   276,   276,   276,   276,   276,   276,   277,
     277,   277,   277,   278,   278,   278,   278,   279,   280,   281,
     282,   282,   283,   283,   284,   285,   285,   285,   285,   286,
     287,   288,   288,   289,   290,   291,   291,   292,   292,   293,
     293,   294,   294,   295,   296,   296,   296,   296,   296,   296,
     296,   296,   296,   296,   296,   296,   296,   296,   296,   296,
     296,   296,   296,   296,   296,   296,   296,   296,   296,   296,
     296,   296,   296,   296,   297,   297,   298,   298,   298,   298,
     299,   300,   300,   300,   300,   300,   300,   300,   301,   302,
     303,   304,   305,   306,   307,   307,   308,   308,   308,   309,
     310,   311,   312,   313,   314,   315,   315,   316,   316,   316,
     316,   316,   316,   316,   316,   316,   316,   316,   316,   316,
     316,   316,   316,   316,   316,   316,   317,   318,   318,   319,
     320,   320,   321,   321,   321,   321,   321,   321,   321,   321,
     321,   321,   321,   321,   321,   321,   321,   321,   322,   322,
     322,   323,   323,   323,   323,   323,   324,   324,   324,   325,
     325,   325,   325,   326,   326,   326,   326,   326,   326,   326,
     326,   326,   326,   326,   326,   326,   326,   326,   327,   327,
     327,   327,   328,   329,   329,   329,   330,   330,   331,   332,
     333,   333,   334,   335,   336,   337,   337,   338,   338,   338,
     339,   339,   340,   340,   341,   342,   343,   344,   344,   345,
     346,   347,   347,   349,   348,   350,   351,   351,   352,   352,
     353,   353,   353,   354,   354,   354,   355,   355,   356,   357,
     357,   358,   358,   359,   359,   360,   361,   362,   362,   363,
     364,   364,   365,   366,   367,   368,   369,   370,   370,   370,
     370,   371,   371,   372,   373,   373,   374,   375,   375,   376,
     377,   377,   378,   378,   379,   379,   380,   381,   382,   382,
     382,   383,   383,   383,   384,   385,   386,   387,   387,   387,
     387,   387,   387,   388,   388,   390,   389,   391,   392,   392,
     392,   392,   392,   393,   394,   395,   396,   396,   397,   398,
     398,   398,   399,   400,   400,   401,   401,   403,   402,   404,
     404,   406,   405,   407,   407,   407,   407,   407,   407,   407,
     407,   408,   408,   408,   408,   408,   410,   409,   411,   411,
     411,   411,   411,   412,   412,   413,   414,   414,   414,   414,
     414,   414,   414,   414,   414,   415,   415,   416,   416,   416,
     416,   417,   417,   418,   419,   420,   420,   421,   421,   422,
     423,   425,   424,   426,   427,   428,   428,   429,   430,   430,
     431,   432,   432,   433,   433,   435,   434,   436,   436,   437,
     437,   437,   437,   437,   437,   437,   438,   438,   439,   439,
     439,   439,   440,   441,   442,   443,   443,   444,   444,   444,
     445,   446,   446,   447,   448,   448,   449,   450,   450,   452,
     453,   451,   455,   454,   456,   457,   458,   458,   460,   461,
     459,   462,   462,   463,   464,   464,   466,   467,   468,   465,
     469,   470,   470,   471,   471,   472,   473,   474,   474,   476,
     475,   478,   477,   479,   479,   480,   480,   480,   481,   482,
     483,   484,   484,   485,   485,   485,   487,   486,   488,   488,
     488,   489,   490,   491,   491,   492,   493,   493,   493,   493,
     493,   493,   493,   493,   493,   495,   494,   494,   496,   494,
     497,   497,   497,   497,   497,   498,   499,   500,   501,   501,
     501,   502,   503,   504,   505,   505,   505,   506,   506,   507,
     507,   507,   507,   507,   507,   507,   507,   508,   508,   510,
     509,   511,   511,   512,   512,   513,   513,   513,   514,   514,
     515,   516,   517,   518,   518,   519,   519,   519,   520,   521,
     521,   522,   522,   523,   523,   524,   524,   525,   525,   526,
     527,   527,   528,   529,   529,   530,   531,   531,   532,   533,
     533,   533,   534,   535,   536,   536,   537,   538,   539,   539,
     540,   541,   542,   543,   543,   543,   543
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
#line 622 "./mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.sym) = mdl_existing_object(parse_state, (yyvsp[0].str))); }
#line 3111 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 29:
#line 625 "./mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.sym) = mdl_existing_region(parse_state, (yyvsp[-3].sym), (yyvsp[-1].str))); }
#line 3117 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 30:
#line 628 "./mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.vec3) = mdl_point(parse_state, &(yyvsp[0].nlist))); }
#line 3123 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 32:
#line 632 "./mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.vec3) = mdl_point_scalar((yyvsp[0].dbl))); }
#line 3129 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 33:
#line 635 "./mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = 1; }
#line 3135 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 34:
#line 636 "./mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = 0; }
#line 3141 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 35:
#line 637 "./mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = 1; }
#line 3147 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 36:
#line 638 "./mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = 0; }
#line 3153 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 37:
#line 639 "./mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = 1; }
#line 3159 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 38:
#line 640 "./mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = 0; }
#line 3165 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 39:
#line 643 "./mdlparse.y" /* yacc.c:1646  */
    { (yyval.mol_type).orient_set = 0; }
#line 3171 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 42:
#line 646 "./mdlparse.y" /* yacc.c:1646  */
    { (yyval.mol_type).orient_set = 1; (yyval.mol_type).orient = 0; }
#line 3177 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 43:
#line 650 "./mdlparse.y" /* yacc.c:1646  */
    { (yyval.mol_type).orient = 1; (yyval.mol_type).orient_set = 1; }
#line 3183 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 44:
#line 651 "./mdlparse.y" /* yacc.c:1646  */
    { (yyval.mol_type).orient = -1; (yyval.mol_type).orient_set = 1; }
#line 3189 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 45:
#line 652 "./mdlparse.y" /* yacc.c:1646  */
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
#line 3204 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 46:
#line 662 "./mdlparse.y" /* yacc.c:1646  */
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
#line 3219 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 49:
#line 680 "./mdlparse.y" /* yacc.c:1646  */
    {
                                                          (yyval.mol_type).orient = (int) (yyvsp[-1].dbl);
                                                          (yyval.mol_type).orient_set = 1;
                                                          if ((yyval.mol_type).orient != (yyvsp[-1].dbl))
                                                          {
                                                            mdlerror(parse_state, "molecule orientation specified inside braces must be an integer between -32768 and 32767.");
                                                            return 1;
                                                          }
                                                      }
#line 3233 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 51:
#line 693 "./mdlparse.y" /* yacc.c:1646  */
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
#line 3249 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 52:
#line 706 "./mdlparse.y" /* yacc.c:1646  */
    { CHECK(mcell_generate_range_singleton(&(yyval.nlist), (yyvsp[0].dbl))); }
#line 3255 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 53:
#line 707 "./mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_generate_range(parse_state, &(yyval.nlist), (yyvsp[-5].dbl), (yyvsp[-3].dbl), (yyvsp[-1].dbl))); }
#line 3261 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 54:
#line 713 "./mdlparse.y" /* yacc.c:1646  */
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
#line 3283 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 55:
#line 736 "./mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_assign_variable_double(parse_state, (yyvsp[-2].sym), (yyvsp[0].dbl))); }
#line 3289 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 56:
#line 737 "./mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_assign_variable_string(parse_state, (yyvsp[-2].sym), (yyvsp[0].str))); }
#line 3295 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 57:
#line 738 "./mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_assign_variable(parse_state, (yyvsp[-2].sym), (yyvsp[0].sym))); }
#line 3301 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 58:
#line 739 "./mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_assign_variable_array(parse_state, (yyvsp[-2].sym), (yyvsp[0].nlist).value_head)); }
#line 3307 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 59:
#line 742 "./mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.sym) = mdl_get_or_create_variable(parse_state, (yyvsp[0].str))); }
#line 3313 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 60:
#line 745 "./mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.sym) = mdl_existing_variable(parse_state, (yyvsp[0].str))); }
#line 3319 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 62:
#line 749 "./mdlparse.y" /* yacc.c:1646  */
    {
                                                          struct num_expr_list *elp;
                                                          (yyval.nlist).value_head = (struct num_expr_list *) (yyvsp[0].sym)->value;
                                                          (yyval.nlist).value_count = 1;
                                                          for (elp = (yyval.nlist).value_head; elp->next != NULL; elp = elp->next)
                                                            ++ (yyval.nlist).value_count;
                                                          (yyval.nlist).value_tail = elp;
                                                          (yyval.nlist).shared = 1;
                                                      }
#line 3333 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 63:
#line 760 "./mdlparse.y" /* yacc.c:1646  */
    { mdl_debug_dump_array((yyvsp[-1].nlist).value_head); (yyval.nlist) = (yyvsp[-1].nlist); }
#line 3339 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 64:
#line 763 "./mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.sym) = mdl_existing_array(parse_state, (yyvsp[0].str))); }
#line 3345 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 68:
#line 771 "./mdlparse.y" /* yacc.c:1646  */
    { (yyval.dbl) = *(double *) (yyvsp[0].sym)->value; }
#line 3351 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 69:
#line 774 "./mdlparse.y" /* yacc.c:1646  */
    { (yyval.dbl) = (yyvsp[0].llival); }
#line 3357 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 73:
#line 782 "./mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.sym) = mdl_existing_double(parse_state, (yyvsp[0].str))); }
#line 3363 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 74:
#line 786 "./mdlparse.y" /* yacc.c:1646  */
    { (yyval.dbl) = (yyvsp[-1].dbl); }
#line 3369 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 75:
#line 787 "./mdlparse.y" /* yacc.c:1646  */
    { CHECKF((yyval.dbl) = exp((yyvsp[-1].dbl))); }
#line 3375 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 76:
#line 788 "./mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_expr_log(parse_state, (yyvsp[-1].dbl), &(yyval.dbl))); }
#line 3381 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 77:
#line 789 "./mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_expr_log10(parse_state, (yyvsp[-1].dbl), &(yyval.dbl))); }
#line 3387 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 78:
#line 790 "./mdlparse.y" /* yacc.c:1646  */
    { (yyval.dbl) = max2d((yyvsp[-3].dbl), (yyvsp[-1].dbl)); }
#line 3393 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 79:
#line 791 "./mdlparse.y" /* yacc.c:1646  */
    { (yyval.dbl) = min2d((yyvsp[-3].dbl), (yyvsp[-1].dbl)); }
#line 3399 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 80:
#line 792 "./mdlparse.y" /* yacc.c:1646  */
    { (yyval.dbl) = mdl_expr_roundoff((yyvsp[-1].dbl), (int) (yyvsp[-3].dbl)); }
#line 3405 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 81:
#line 793 "./mdlparse.y" /* yacc.c:1646  */
    { (yyval.dbl) = floor((yyvsp[-1].dbl)); }
#line 3411 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 82:
#line 794 "./mdlparse.y" /* yacc.c:1646  */
    { (yyval.dbl) = ceil((yyvsp[-1].dbl)); }
#line 3417 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 83:
#line 795 "./mdlparse.y" /* yacc.c:1646  */
    { (yyval.dbl) = sin((yyvsp[-1].dbl)); }
#line 3423 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 84:
#line 796 "./mdlparse.y" /* yacc.c:1646  */
    { (yyval.dbl) = cos((yyvsp[-1].dbl)); }
#line 3429 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 85:
#line 797 "./mdlparse.y" /* yacc.c:1646  */
    { CHECKF((yyval.dbl) = tan((yyvsp[-1].dbl))); }
#line 3435 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 86:
#line 798 "./mdlparse.y" /* yacc.c:1646  */
    { CHECKF((yyval.dbl) = asin((yyvsp[-1].dbl))); }
#line 3441 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 87:
#line 799 "./mdlparse.y" /* yacc.c:1646  */
    { CHECKF((yyval.dbl) = acos((yyvsp[-1].dbl))); }
#line 3447 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 88:
#line 800 "./mdlparse.y" /* yacc.c:1646  */
    { (yyval.dbl) = atan((yyvsp[-1].dbl)); }
#line 3453 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 89:
#line 801 "./mdlparse.y" /* yacc.c:1646  */
    { CHECKF((yyval.dbl) = sqrt((yyvsp[-1].dbl))); }
#line 3459 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 90:
#line 802 "./mdlparse.y" /* yacc.c:1646  */
    { (yyval.dbl) = fabs((yyvsp[-1].dbl)); }
#line 3465 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 91:
#line 803 "./mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_expr_mod(parse_state, (yyvsp[-3].dbl), (yyvsp[-1].dbl), &(yyval.dbl))); }
#line 3471 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 92:
#line 804 "./mdlparse.y" /* yacc.c:1646  */
    { (yyval.dbl) = MY_PI; }
#line 3477 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 93:
#line 805 "./mdlparse.y" /* yacc.c:1646  */
    { (yyval.dbl) = mdl_expr_rng_uniform(parse_state); }
#line 3483 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 94:
#line 806 "./mdlparse.y" /* yacc.c:1646  */
    { (yyval.dbl) = rng_gauss(parse_state->vol->rng); }
#line 3489 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 95:
#line 807 "./mdlparse.y" /* yacc.c:1646  */
    { (yyval.dbl) = parse_state->vol->seed_seq; }
#line 3495 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 96:
#line 808 "./mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_expr_string_to_double(parse_state, (yyvsp[-1].str), &(yyval.dbl))); }
#line 3501 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 97:
#line 809 "./mdlparse.y" /* yacc.c:1646  */
    { CHECKF((yyval.dbl) = (yyvsp[-2].dbl) + (yyvsp[0].dbl)); }
#line 3507 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 98:
#line 810 "./mdlparse.y" /* yacc.c:1646  */
    { CHECKF((yyval.dbl) = (yyvsp[-2].dbl) - (yyvsp[0].dbl)); }
#line 3513 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 99:
#line 811 "./mdlparse.y" /* yacc.c:1646  */
    { CHECKF((yyval.dbl) = (yyvsp[-2].dbl) * (yyvsp[0].dbl)); }
#line 3519 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 100:
#line 812 "./mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_expr_div(parse_state, (yyvsp[-2].dbl), (yyvsp[0].dbl), &(yyval.dbl))); }
#line 3525 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 101:
#line 813 "./mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_expr_pow(parse_state, (yyvsp[-2].dbl), (yyvsp[0].dbl), &(yyval.dbl))); }
#line 3531 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 102:
#line 814 "./mdlparse.y" /* yacc.c:1646  */
    { (yyval.dbl) = -(yyvsp[0].dbl); }
#line 3537 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 103:
#line 815 "./mdlparse.y" /* yacc.c:1646  */
    { (yyval.dbl) = (yyvsp[0].dbl); }
#line 3543 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 105:
#line 820 "./mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.str) = mdl_strdup((char const *) (yyvsp[0].sym)->value)); }
#line 3549 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 106:
#line 824 "./mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.str) = mdl_strip_quotes((yyvsp[0].str))); }
#line 3555 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 107:
#line 825 "./mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.str) = mdl_strdup(parse_state->vol->mdl_infile_name)); }
#line 3561 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 108:
#line 826 "./mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.str) = mdl_strcat((yyvsp[-2].str), (yyvsp[0].str))); }
#line 3567 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 109:
#line 827 "./mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.str) = mdl_string_format(parse_state, (yyvsp[-2].str), (yyvsp[-1].printfargs).arg_head)); }
#line 3573 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 110:
#line 830 "./mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.sym) = mdl_existing_string(parse_state, (yyvsp[0].str))); }
#line 3579 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 118:
#line 846 "./mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_fopen(parse_state, (yyvsp[-6].sym), (yyvsp[-3].str), (yyvsp[-1].str))); }
#line 3585 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 119:
#line 849 "./mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.sym) = mdl_new_filehandle(parse_state, (yyvsp[0].str))); }
#line 3591 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 120:
#line 852 "./mdlparse.y" /* yacc.c:1646  */
    { (yyval.str) = (yyvsp[0].str); CHECK(mdl_valid_file_mode(parse_state, (yyvsp[0].str))); }
#line 3597 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 121:
#line 855 "./mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_fclose(parse_state, (yyvsp[-1].sym))); }
#line 3603 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 122:
#line 858 "./mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.sym) = mdl_existing_file_stream(parse_state, (yyvsp[0].str))); }
#line 3609 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 123:
#line 861 "./mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.str) = mdl_expand_string_escapes((yyvsp[0].str))); }
#line 3615 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 124:
#line 864 "./mdlparse.y" /* yacc.c:1646  */
    { (yyval.printfargs).arg_head = (yyval.printfargs).arg_tail = NULL; }
#line 3621 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 125:
#line 865 "./mdlparse.y" /* yacc.c:1646  */
    {
                                                        (yyval.printfargs) = (yyvsp[-2].printfargs);
                                                        if ((yyval.printfargs).arg_tail)
                                                          (yyval.printfargs).arg_tail = (yyval.printfargs).arg_tail->next = (yyvsp[0].printfarg);
                                                        else
                                                          (yyval.printfargs).arg_tail = (yyval.printfargs).arg_head = (yyvsp[0].printfarg);
                                                        (yyvsp[0].printfarg)->next = NULL;
                                                      }
#line 3634 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 126:
#line 875 "./mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.printfarg) = mdl_new_printf_arg_double((yyvsp[0].dbl))); }
#line 3640 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 127:
#line 876 "./mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.printfarg) = mdl_new_printf_arg_string((yyvsp[0].str))); }
#line 3646 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 128:
#line 877 "./mdlparse.y" /* yacc.c:1646  */
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
#line 3661 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 129:
#line 889 "./mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_printf(parse_state, (yyvsp[-2].str), (yyvsp[-1].printfargs).arg_head)); }
#line 3667 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 130:
#line 895 "./mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_fprintf(parse_state, (struct file_stream *) (yyvsp[-4].sym)->value, (yyvsp[-2].str), (yyvsp[-1].printfargs).arg_head)); }
#line 3673 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 131:
#line 901 "./mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_sprintf(parse_state, (yyvsp[-4].sym), (yyvsp[-2].str), (yyvsp[-1].printfargs).arg_head)); }
#line 3679 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 132:
#line 904 "./mdlparse.y" /* yacc.c:1646  */
    { mdl_print_time(parse_state, (yyvsp[-1].str)); }
#line 3685 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 133:
#line 910 "./mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_fprint_time(parse_state, (yyvsp[-3].sym), (yyvsp[-1].str))); }
#line 3691 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 137:
#line 926 "./mdlparse.y" /* yacc.c:1646  */
    { if (!parse_state->vol->quiet_flag) mdl_set_all_notifications(parse_state->vol, (yyvsp[0].tok)); }
#line 3697 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 138:
#line 927 "./mdlparse.y" /* yacc.c:1646  */
    { if (!parse_state->vol->quiet_flag) parse_state->vol->notify->progress_report        = (yyvsp[0].tok); }
#line 3703 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 139:
#line 928 "./mdlparse.y" /* yacc.c:1646  */
    { if (!parse_state->vol->quiet_flag) parse_state->vol->notify->diffusion_constants    = (yyvsp[0].tok); }
#line 3709 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 140:
#line 929 "./mdlparse.y" /* yacc.c:1646  */
    { if (!parse_state->vol->quiet_flag) parse_state->vol->notify->reaction_probabilities = (yyvsp[0].tok); }
#line 3715 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 141:
#line 930 "./mdlparse.y" /* yacc.c:1646  */
    { if (!parse_state->vol->quiet_flag) parse_state->vol->notify->time_varying_reactions = (yyvsp[0].tok); }
#line 3721 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 142:
#line 931 "./mdlparse.y" /* yacc.c:1646  */
    { if (!parse_state->vol->quiet_flag) parse_state->vol->notify->reaction_prob_notify   = (yyvsp[0].dbl); }
#line 3727 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 143:
#line 932 "./mdlparse.y" /* yacc.c:1646  */
    { if (!parse_state->vol->quiet_flag) parse_state->vol->notify->partition_location     = (yyvsp[0].tok); }
#line 3733 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 144:
#line 933 "./mdlparse.y" /* yacc.c:1646  */
    { if (!parse_state->vol->quiet_flag) parse_state->vol->notify->box_triangulation      = (yyvsp[0].tok); }
#line 3739 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 145:
#line 934 "./mdlparse.y" /* yacc.c:1646  */
    { if (!parse_state->vol->quiet_flag) parse_state->vol->notify->release_events         = (yyvsp[0].tok); }
#line 3745 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 146:
#line 935 "./mdlparse.y" /* yacc.c:1646  */
    { if (!parse_state->vol->quiet_flag) parse_state->vol->notify->file_writes            = (yyvsp[0].tok); }
#line 3751 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 147:
#line 936 "./mdlparse.y" /* yacc.c:1646  */
    { if (!parse_state->vol->quiet_flag) parse_state->vol->notify->final_summary          = (yyvsp[0].tok); }
#line 3757 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 148:
#line 937 "./mdlparse.y" /* yacc.c:1646  */
    { if (!parse_state->vol->quiet_flag) parse_state->vol->notify->throughput_report      = (yyvsp[0].tok); }
#line 3763 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 149:
#line 938 "./mdlparse.y" /* yacc.c:1646  */
    { if (!parse_state->vol->quiet_flag) parse_state->vol->notify->reaction_output_report = (yyvsp[0].tok); }
#line 3769 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 150:
#line 939 "./mdlparse.y" /* yacc.c:1646  */
    { if (!parse_state->vol->quiet_flag) parse_state->vol->notify->volume_output_report   = (yyvsp[0].tok); }
#line 3775 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 151:
#line 940 "./mdlparse.y" /* yacc.c:1646  */
    { if (!parse_state->vol->quiet_flag) parse_state->vol->notify->viz_output_report      = (yyvsp[0].tok); }
#line 3781 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 152:
#line 941 "./mdlparse.y" /* yacc.c:1646  */
    { if (!parse_state->vol->quiet_flag) parse_state->vol->notify->checkpoint_report      = (yyvsp[0].tok); }
#line 3787 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 153:
#line 942 "./mdlparse.y" /* yacc.c:1646  */
    {
                                                          if (!parse_state->vol->quiet_flag && parse_state->vol->log_freq == ULONG_MAX)
                                                            parse_state->vol->notify->iteration_report = (yyvsp[0].tok);
                                                      }
#line 3796 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 154:
#line 946 "./mdlparse.y" /* yacc.c:1646  */
    { if (!parse_state->vol->quiet_flag) CHECK(mdl_set_iteration_report_freq(parse_state, (long long) (yyvsp[0].dbl))); }
#line 3802 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 155:
#line 947 "./mdlparse.y" /* yacc.c:1646  */
    { if (!parse_state->vol->quiet_flag) parse_state->vol->notify->molecule_collision_report    = (yyvsp[0].tok); }
#line 3808 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 156:
#line 951 "./mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = ((yyvsp[0].tok) ? NOTIFY_FULL : NOTIFY_NONE); }
#line 3814 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 157:
#line 955 "./mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = ((yyvsp[0].tok) ? NOTIFY_FULL : NOTIFY_NONE); }
#line 3820 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 158:
#line 956 "./mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = NOTIFY_BRIEF; }
#line 3826 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 162:
#line 972 "./mdlparse.y" /* yacc.c:1646  */
    { mdl_set_all_warnings(parse_state->vol, (byte) (yyvsp[0].tok)); }
#line 3832 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 163:
#line 973 "./mdlparse.y" /* yacc.c:1646  */
    { parse_state->vol->notify->neg_diffusion = (byte)(yyvsp[0].tok); }
#line 3838 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 164:
#line 974 "./mdlparse.y" /* yacc.c:1646  */
    { parse_state->vol->notify->neg_reaction = (byte)(yyvsp[0].tok); }
#line 3844 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 165:
#line 975 "./mdlparse.y" /* yacc.c:1646  */
    { parse_state->vol->notify->high_reaction_prob = (byte)(yyvsp[0].tok); }
#line 3850 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 166:
#line 976 "./mdlparse.y" /* yacc.c:1646  */
    { parse_state->vol->notify->reaction_prob_warn = (yyvsp[0].dbl); }
#line 3856 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 167:
#line 977 "./mdlparse.y" /* yacc.c:1646  */
    { parse_state->vol->notify->close_partitions = (byte)(yyvsp[0].tok); }
#line 3862 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 168:
#line 978 "./mdlparse.y" /* yacc.c:1646  */
    { parse_state->vol->notify->degenerate_polys = (byte)(yyvsp[0].tok); }
#line 3868 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 169:
#line 979 "./mdlparse.y" /* yacc.c:1646  */
    { parse_state->vol->notify->overwritten_file = (byte)(yyvsp[0].tok); }
#line 3874 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 170:
#line 980 "./mdlparse.y" /* yacc.c:1646  */
    { parse_state->vol->notify->short_lifetime = (byte)(yyvsp[0].tok); }
#line 3880 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 171:
#line 981 "./mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_lifetime_warning_threshold(parse_state, (long long) (yyvsp[0].dbl))); }
#line 3886 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 172:
#line 982 "./mdlparse.y" /* yacc.c:1646  */
    { parse_state->vol->notify->missed_reactions = (byte)(yyvsp[0].tok); }
#line 3892 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 173:
#line 983 "./mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_missed_reaction_warning_threshold(parse_state, (yyvsp[0].dbl))); }
#line 3898 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 174:
#line 984 "./mdlparse.y" /* yacc.c:1646  */
    { parse_state->vol->notify->missed_surf_orient = (byte)(yyvsp[0].tok); }
#line 3904 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 175:
#line 985 "./mdlparse.y" /* yacc.c:1646  */
    { parse_state->vol->notify->useless_vol_orient = (byte)(yyvsp[0].tok); }
#line 3910 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 176:
#line 986 "./mdlparse.y" /* yacc.c:1646  */
    { parse_state->vol->notify->mol_placement_failure = (byte) (yyvsp[0].tok); }
#line 3916 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 177:
#line 987 "./mdlparse.y" /* yacc.c:1646  */
    { parse_state->vol->notify->invalid_output_step_time = (byte) (yyvsp[0].tok); }
#line 3922 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 178:
#line 991 "./mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = WARN_COPE;  }
#line 3928 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 179:
#line 992 "./mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = WARN_WARN;  }
#line 3934 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 180:
#line 993 "./mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = WARN_ERROR; }
#line 3940 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 181:
#line 999 "./mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_checkpoint_infile(parse_state, (yyvsp[0].str))); }
#line 3946 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 182:
#line 1000 "./mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_checkpoint_outfile(parse_state, (yyvsp[0].str))); }
#line 3952 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 183:
#line 1001 "./mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_checkpoint_interval(parse_state, (yyvsp[-1].dbl), (yyvsp[0].tok))); }
#line 3958 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 184:
#line 1002 "./mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_keep_checkpoint_files(parse_state, (yyvsp[0].tok))); }
#line 3964 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 185:
#line 1004 "./mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_realtime_checkpoint(parse_state, (long) (yyvsp[-1].dbl), (yyvsp[0].tok))); }
#line 3970 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 186:
#line 1007 "./mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = 0; }
#line 3976 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 187:
#line 1008 "./mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = 1; }
#line 3982 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 188:
#line 1009 "./mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = 0; }
#line 3988 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 189:
#line 1013 "./mdlparse.y" /* yacc.c:1646  */
    { /* seconds */     (yyval.dbl) = (yyvsp[0].dbl); }
#line 3994 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 190:
#line 1014 "./mdlparse.y" /* yacc.c:1646  */
    { /* mm:ss */       (yyval.dbl) = (yyvsp[-2].dbl) * 60 + (yyvsp[0].dbl); }
#line 4000 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 191:
#line 1015 "./mdlparse.y" /* yacc.c:1646  */
    { /* hh:mm:ss */    (yyval.dbl) = (yyvsp[-4].dbl) * 3600 + (yyvsp[-2].dbl) * 60 + (yyvsp[0].dbl); }
#line 4006 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 192:
#line 1017 "./mdlparse.y" /* yacc.c:1646  */
    { /* dd:hh:mm:ss */ (yyval.dbl) = (yyvsp[-6].dbl) * 86400 + (yyvsp[-4].dbl) * 3600 + (yyvsp[-2].dbl) * 60 + (yyvsp[0].dbl); }
#line 4012 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 193:
#line 1024 "./mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_time_step(parse_state, (yyvsp[0].dbl))); }
#line 4018 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 194:
#line 1025 "./mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_space_step(parse_state, (yyvsp[0].dbl))); }
#line 4024 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 195:
#line 1026 "./mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_max_time_step(parse_state, (yyvsp[0].dbl))); }
#line 4030 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 196:
#line 1027 "./mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_num_iterations(parse_state, (long long) (yyvsp[0].dbl))); }
#line 4036 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 197:
#line 1028 "./mdlparse.y" /* yacc.c:1646  */
    { parse_state->vol->randomize_smol_pos = !((yyvsp[0].tok)); }
#line 4042 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 198:
#line 1029 "./mdlparse.y" /* yacc.c:1646  */
    { parse_state->vol->use_expanded_list = (yyvsp[0].tok); }
#line 4048 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 199:
#line 1030 "./mdlparse.y" /* yacc.c:1646  */
    { parse_state->vol->vacancy_search_dist2 = max2d((yyvsp[0].dbl), 0.0); }
#line 4054 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 200:
#line 1031 "./mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_num_radial_directions(parse_state, (int) (yyvsp[0].dbl))); }
#line 4060 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 201:
#line 1032 "./mdlparse.y" /* yacc.c:1646  */
    { parse_state->vol->fully_random = 1; }
#line 4066 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 202:
#line 1033 "./mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_num_radial_subdivisions(parse_state, (int) (yyvsp[0].dbl))); }
#line 4072 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 203:
#line 1034 "./mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_grid_density(parse_state, (yyvsp[0].dbl))); }
#line 4078 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 204:
#line 1035 "./mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_interaction_radius(parse_state, (yyvsp[0].dbl))); }
#line 4084 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 205:
#line 1036 "./mdlparse.y" /* yacc.c:1646  */
    { parse_state->vol->surface_reversibility=(yyvsp[0].tok); parse_state->vol->volume_reversibility=(yyvsp[0].tok); }
#line 4090 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 206:
#line 1037 "./mdlparse.y" /* yacc.c:1646  */
    { parse_state->vol->surface_reversibility=1;  parse_state->vol->volume_reversibility=0;  }
#line 4096 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 207:
#line 1038 "./mdlparse.y" /* yacc.c:1646  */
    { parse_state->vol->surface_reversibility=0;  parse_state->vol->volume_reversibility=1;  }
#line 4102 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 208:
#line 1045 "./mdlparse.y" /* yacc.c:1646  */
    { parse_state->vol->mem_part_x = (int) (yyvsp[0].dbl); }
#line 4108 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 209:
#line 1046 "./mdlparse.y" /* yacc.c:1646  */
    { parse_state->vol->mem_part_y = (int) (yyvsp[0].dbl); }
#line 4114 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 210:
#line 1047 "./mdlparse.y" /* yacc.c:1646  */
    { parse_state->vol->mem_part_z = (int) (yyvsp[0].dbl); }
#line 4120 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 211:
#line 1048 "./mdlparse.y" /* yacc.c:1646  */
    { parse_state->vol->mem_part_pool = (int) (yyvsp[0].dbl); }
#line 4126 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 212:
#line 1052 "./mdlparse.y" /* yacc.c:1646  */
    { CHECK(mcell_set_partition(parse_state->vol, (yyvsp[-2].tok), & (yyvsp[0].nlist))); }
#line 4132 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 213:
#line 1056 "./mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = X_PARTS; }
#line 4138 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 214:
#line 1057 "./mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = Y_PARTS; }
#line 4144 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 215:
#line 1058 "./mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = Z_PARTS; }
#line 4150 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 218:
#line 1069 "./mdlparse.y" /* yacc.c:1646  */
    { mdl_print_species_summary(parse_state->vol, (yyvsp[0].mcell_mol_spec)); }
#line 4156 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 219:
#line 1073 "./mdlparse.y" /* yacc.c:1646  */
    { mdl_print_species_summaries(parse_state->vol, (yyvsp[-1].mcell_species_lst).species_head); }
#line 4162 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 220:
#line 1077 "./mdlparse.y" /* yacc.c:1646  */
    { (yyval.mcell_species_lst).species_count = 0; CHECK(mdl_add_to_species_list(&(yyval.mcell_species_lst), (yyvsp[0].mcell_mol_spec))); }
#line 4168 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 221:
#line 1078 "./mdlparse.y" /* yacc.c:1646  */
    { (yyval.mcell_species_lst) = (yyvsp[-1].mcell_species_lst); CHECK(mdl_add_to_species_list(&(yyval.mcell_species_lst), (yyvsp[0].mcell_mol_spec))); }
#line 4174 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 222:
#line 1087 "./mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.mcell_mol_spec) = mdl_create_species(parse_state, (yyvsp[-6].str), (yyvsp[-4].diff_const).D, (yyvsp[-4].diff_const).is_2d, (yyvsp[-3].dbl), (yyvsp[-2].ival), (yyvsp[-1].dbl) )); }
#line 4180 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 224:
#line 1093 "./mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.sym) = mdl_new_mol_species(parse_state, (yyvsp[0].str))); }
#line 4186 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 225:
#line 1097 "./mdlparse.y" /* yacc.c:1646  */
    { (yyval.diff_const).is_2d = 0; (yyval.diff_const).D = (yyvsp[0].dbl); CHECK(mdl_check_diffusion_constant(parse_state, & (yyval.diff_const).D)); }
#line 4192 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 226:
#line 1098 "./mdlparse.y" /* yacc.c:1646  */
    { (yyval.diff_const).is_2d = 1; (yyval.diff_const).D = (yyvsp[0].dbl); CHECK(mdl_check_diffusion_constant(parse_state, & (yyval.diff_const).D)); }
#line 4198 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 227:
#line 1102 "./mdlparse.y" /* yacc.c:1646  */
    { (yyval.dbl) = 0.0; }
#line 4204 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 228:
#line 1103 "./mdlparse.y" /* yacc.c:1646  */
    {
                                                          if ((yyvsp[0].dbl) <= 0)
                                                          {
                                                            mdlerror_fmt(parse_state, "Requested custom time step of %.15g; custom time step must be positive.", (yyvsp[0].dbl));
                                                            return 1;
                                                          }

                                                          (yyval.dbl) = (yyvsp[0].dbl);
                                                      }
#line 4218 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 229:
#line 1112 "./mdlparse.y" /* yacc.c:1646  */
    {
                                                          if ((yyvsp[0].dbl) <= 0)
                                                          {
                                                            mdlerror_fmt(parse_state, "Requested custom space step of %.15g; custom space step must be positive.", (yyvsp[0].dbl));
                                                            return 1;
                                                          }

                                                          (yyval.dbl) = -(yyvsp[0].dbl);
                                                      }
#line 4232 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 230:
#line 1123 "./mdlparse.y" /* yacc.c:1646  */
    { (yyval.ival) = 0; }
#line 4238 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 231:
#line 1124 "./mdlparse.y" /* yacc.c:1646  */
    { (yyval.ival) = 1; }
#line 4244 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 232:
#line 1128 "./mdlparse.y" /* yacc.c:1646  */
    { (yyval.dbl) = 0; }
#line 4250 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 233:
#line 1129 "./mdlparse.y" /* yacc.c:1646  */
    {
                                                        if ((yyvsp[0].dbl) <= 0)
                                                        {
                                                          mdlerror_fmt(parse_state, "Requested maximum step length of %.15g; maximum step length must be positive.", (yyvsp[0].dbl));
                                                          return 1;
                                                        }
                                                        (yyval.dbl) = (yyvsp[0].dbl);
                                                      }
#line 4263 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 234:
#line 1139 "./mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.sym) = mdl_existing_molecule(parse_state, (yyvsp[0].str))); }
#line 4269 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 235:
#line 1143 "./mdlparse.y" /* yacc.c:1646  */
    { (yyval.mol_type) = (yyvsp[0].mol_type); CHECKN((yyval.mol_type).mol_type = mdl_existing_surface_molecule(parse_state, (yyvsp[-1].str))); }
#line 4275 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 236:
#line 1147 "./mdlparse.y" /* yacc.c:1646  */
    {
                                                        (yyval.mol_type) = (yyvsp[0].mol_type);
                                                        if (! (yyval.mol_type).orient_set)
                                                          (yyval.mol_type).orient = 0;
                                                        (yyval.mol_type).mol_type = (yyvsp[-1].sym);
                                                      }
#line 4286 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 243:
#line 1182 "./mdlparse.y" /* yacc.c:1646  */
    { mdl_start_surface_class(parse_state, (yyvsp[-1].sym)); }
#line 4292 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 244:
#line 1184 "./mdlparse.y" /* yacc.c:1646  */
    { mdl_finish_surface_class(parse_state); }
#line 4298 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 245:
#line 1187 "./mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.sym) = mdl_existing_surface_class(parse_state, (yyvsp[0].str))); }
#line 4304 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 250:
#line 1204 "./mdlparse.y" /* yacc.c:1646  */
    { CHECKN(mdl_assemble_surface_reaction(parse_state, (yyvsp[-2].tok), parse_state->current_surface_class, (yyvsp[0].mol_type).mol_type, (yyvsp[0].mol_type).orient)); }
#line 4310 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 251:
#line 1207 "./mdlparse.y" /* yacc.c:1646  */
    {
              struct sym_entry *mol_sym = retrieve_sym("ALL_MOLECULES", parse_state->vol->mol_sym_table);
              if(!(yyvsp[0].mol_type).orient_set) (yyvsp[0].mol_type).orient = 0;
              CHECKN(mdl_assemble_surface_reaction(parse_state, (yyvsp[-3].tok), parse_state->current_surface_class, mol_sym, (yyvsp[0].mol_type).orient));}
#line 4319 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 252:
#line 1213 "./mdlparse.y" /* yacc.c:1646  */
    { CHECKN(mdl_assemble_concentration_clamp_reaction(parse_state, parse_state->current_surface_class, (yyvsp[-2].mol_type).mol_type, (yyvsp[-2].mol_type).orient, (yyvsp[0].dbl))); }
#line 4325 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 253:
#line 1216 "./mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = RFLCT; }
#line 4331 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 254:
#line 1217 "./mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = TRANSP; }
#line 4337 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 255:
#line 1218 "./mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = SINK; }
#line 4343 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 258:
#line 1225 "./mdlparse.y" /* yacc.c:1646  */
    { parse_state->current_surface_class->sm_dat_head = (yyvsp[0].surf_mol_dat_list).sm_head; }
#line 4349 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 259:
#line 1232 "./mdlparse.y" /* yacc.c:1646  */
    { (yyval.surf_mol_dat_list) = (yyvsp[-1].surf_mol_dat_list); }
#line 4355 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 260:
#line 1236 "./mdlparse.y" /* yacc.c:1646  */
    { (yyval.surf_mol_dat_list) = (yyvsp[-1].surf_mol_dat_list); }
#line 4361 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 261:
#line 1240 "./mdlparse.y" /* yacc.c:1646  */
    {
                                                          (yyvsp[0].surf_mol_dat)->quantity_type = SURFMOLDENS;
                                                          (yyval.surf_mol_dat_list).sm_tail = (yyval.surf_mol_dat_list).sm_head = (yyvsp[0].surf_mol_dat);
                                                      }
#line 4370 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 262:
#line 1245 "./mdlparse.y" /* yacc.c:1646  */
    {
                                                          (yyval.surf_mol_dat_list) = (yyvsp[-1].surf_mol_dat_list);
                                                          (yyvsp[0].surf_mol_dat)->quantity_type = SURFMOLDENS;
                                                          (yyval.surf_mol_dat_list).sm_tail = (yyval.surf_mol_dat_list).sm_tail->next = (yyvsp[0].surf_mol_dat);
                                                      }
#line 4380 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 263:
#line 1253 "./mdlparse.y" /* yacc.c:1646  */
    {
                                                          (yyvsp[0].surf_mol_dat)->quantity_type = SURFMOLNUM;
                                                          (yyval.surf_mol_dat_list).sm_tail = (yyval.surf_mol_dat_list).sm_head = (yyvsp[0].surf_mol_dat);
                                                      }
#line 4389 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 264:
#line 1258 "./mdlparse.y" /* yacc.c:1646  */
    {
                                                          (yyval.surf_mol_dat_list) = (yyvsp[-1].surf_mol_dat_list);
                                                          (yyvsp[0].surf_mol_dat)->quantity_type = SURFMOLNUM;
                                                          (yyval.surf_mol_dat_list).sm_tail = (yyval.surf_mol_dat_list).sm_tail->next = (yyvsp[0].surf_mol_dat);
                                                      }
#line 4399 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 265:
#line 1266 "./mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.surf_mol_dat) = mdl_new_surf_mol_data(parse_state, &(yyvsp[-2].mol_type), (yyvsp[0].dbl))); }
#line 4405 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 275:
#line 1295 "./mdlparse.y" /* yacc.c:1646  */
    { (yyval.react_arrow).catalyst = (yyvsp[-1].mol_type); (yyval.react_arrow).flags = ARROW_CATALYTIC; }
#line 4411 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 276:
#line 1300 "./mdlparse.y" /* yacc.c:1646  */
    { (yyval.react_arrow).catalyst = (yyvsp[-1].mol_type); (yyval.react_arrow).flags = ARROW_CATALYTIC | ARROW_BIDIRECTIONAL; }
#line 4417 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 277:
#line 1305 "./mdlparse.y" /* yacc.c:1646  */
    { (yyval.react_arrow).catalyst.mol_type = NULL; (yyval.react_arrow).flags = 0; }
#line 4423 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 279:
#line 1307 "./mdlparse.y" /* yacc.c:1646  */
    { (yyval.react_arrow).catalyst.mol_type = NULL; (yyval.react_arrow).flags = ARROW_BIDIRECTIONAL; }
#line 4429 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 281:
#line 1311 "./mdlparse.y" /* yacc.c:1646  */
    { (yyval.sym) = NULL; }
#line 4435 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 282:
#line 1312 "./mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.sym) = mdl_new_rxn_pathname(parse_state, (yyvsp[0].str))); }
#line 4441 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 283:
#line 1318 "./mdlparse.y" /* yacc.c:1646  */
    { CHECKN(mdl_assemble_reaction(parse_state, (yyvsp[-5].mol_type_list).mol_type_head, &(yyvsp[-4].mol_type), &(yyvsp[-3].react_arrow), (yyvsp[-2].mol_type_list).mol_type_head, &(yyvsp[-1].react_rates), (yyvsp[0].sym))); }
#line 4447 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 284:
#line 1321 "./mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_reaction_player_singleton(parse_state, & (yyval.mol_type_list), & (yyvsp[0].mol_type))); }
#line 4453 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 285:
#line 1322 "./mdlparse.y" /* yacc.c:1646  */
    { (yyval.mol_type_list) = (yyvsp[-2].mol_type_list); CHECK(mdl_add_reaction_player(parse_state, & (yyval.mol_type_list), & (yyvsp[0].mol_type))); }
#line 4459 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 287:
#line 1329 "./mdlparse.y" /* yacc.c:1646  */
    { (yyval.mol_type).mol_type = NULL; }
#line 4465 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 288:
#line 1330 "./mdlparse.y" /* yacc.c:1646  */
    { (yyval.mol_type) = (yyvsp[0].mol_type); }
#line 4471 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 289:
#line 1334 "./mdlparse.y" /* yacc.c:1646  */
    { (yyval.mol_type) = (yyvsp[0].mol_type); (yyval.mol_type).mol_type = (yyvsp[-1].sym); }
#line 4477 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 290:
#line 1337 "./mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_reaction_player_singleton(parse_state, & (yyval.mol_type_list), & (yyvsp[0].mol_type))); }
#line 4483 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 291:
#line 1338 "./mdlparse.y" /* yacc.c:1646  */
    { (yyval.mol_type_list) = (yyvsp[-2].mol_type_list); CHECK(mdl_add_reaction_player(parse_state, & (yyval.mol_type_list), & (yyvsp[0].mol_type))); }
#line 4489 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 292:
#line 1341 "./mdlparse.y" /* yacc.c:1646  */
    { (yyval.mol_type).mol_type = NULL; (yyval.mol_type).orient_set = 0; }
#line 4495 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 296:
#line 1350 "./mdlparse.y" /* yacc.c:1646  */
    {
                                                        if ((yyvsp[-1].react_rates).forward_rate.rate_type == RATE_UNSET)
                                                        {
                                                          mdlerror(parse_state, "invalid reaction rate specification: must specify a forward rate.");
                                                          return 1;
                                                        }

                                                        (yyval.react_rates) = (yyvsp[-1].react_rates);
                                                      }
#line 4509 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 297:
#line 1361 "./mdlparse.y" /* yacc.c:1646  */
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
#line 4528 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 298:
#line 1378 "./mdlparse.y" /* yacc.c:1646  */
    { (yyval.react_rates).forward_rate = (yyvsp[0].react_rate); (yyval.react_rates).backward_rate.rate_type = RATE_UNSET; CHECK(mdl_valid_rate(parse_state, &(yyvsp[0].react_rate))); }
#line 4534 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 299:
#line 1379 "./mdlparse.y" /* yacc.c:1646  */
    { (yyval.react_rates).forward_rate = (yyvsp[0].react_rate); (yyval.react_rates).backward_rate.rate_type = RATE_UNSET; CHECK(mdl_valid_rate(parse_state, &(yyvsp[0].react_rate))); }
#line 4540 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 300:
#line 1380 "./mdlparse.y" /* yacc.c:1646  */
    { (yyval.react_rates).backward_rate = (yyvsp[0].react_rate); (yyval.react_rates).forward_rate.rate_type = RATE_UNSET; CHECK(mdl_valid_rate(parse_state, &(yyvsp[0].react_rate))); }
#line 4546 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 301:
#line 1384 "./mdlparse.y" /* yacc.c:1646  */
    { (yyval.react_rate).rate_type = RATE_CONSTANT; (yyval.react_rate).v.rate_constant = (yyvsp[0].dbl); }
#line 4552 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 302:
#line 1385 "./mdlparse.y" /* yacc.c:1646  */
    { (yyval.react_rate).rate_type = RATE_FILE; (yyval.react_rate).v.rate_file = (yyvsp[0].str); }
#line 4558 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 303:
#line 1386 "./mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_reaction_rate_from_var(parse_state, & (yyval.react_rate), (yyvsp[0].sym))); }
#line 4564 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 304:
#line 1397 "./mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_release_pattern(parse_state, (yyvsp[-3].sym), &(yyvsp[-1].rpat))); }
#line 4570 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 305:
#line 1400 "./mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.sym) = mdl_new_release_pattern(parse_state, (yyvsp[0].str))); }
#line 4576 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 306:
#line 1403 "./mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.sym) = mdl_existing_release_pattern_or_rxn_pathname(parse_state, (yyvsp[0].str))); }
#line 4582 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 307:
#line 1407 "./mdlparse.y" /* yacc.c:1646  */
    {
                                                        (yyval.rpat).delay = 0;
                                                        (yyval.rpat).release_interval = FOREVER;
                                                        (yyval.rpat).train_interval = FOREVER;
                                                        (yyval.rpat).train_duration = FOREVER;
                                                        (yyval.rpat).number_of_trains = 1;
                                                      }
#line 4594 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 308:
#line 1415 "./mdlparse.y" /* yacc.c:1646  */
    { (yyval.rpat) = (yyvsp[-3].rpat); (yyval.rpat).delay = (yyvsp[0].dbl) / parse_state->vol->time_unit; }
#line 4600 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 309:
#line 1417 "./mdlparse.y" /* yacc.c:1646  */
    { (yyval.rpat) = (yyvsp[-3].rpat); (yyval.rpat).release_interval = (yyvsp[0].dbl) / parse_state->vol->time_unit; }
#line 4606 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 310:
#line 1419 "./mdlparse.y" /* yacc.c:1646  */
    { (yyval.rpat) = (yyvsp[-3].rpat); (yyval.rpat).train_interval = (yyvsp[0].dbl) / parse_state->vol->time_unit; }
#line 4612 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 311:
#line 1421 "./mdlparse.y" /* yacc.c:1646  */
    { (yyval.rpat) = (yyvsp[-3].rpat); (yyval.rpat).train_duration = (yyvsp[0].dbl) / parse_state->vol->time_unit; }
#line 4618 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 312:
#line 1423 "./mdlparse.y" /* yacc.c:1646  */
    { (yyval.rpat) = (yyvsp[-3].rpat); (yyval.rpat).number_of_trains = (yyvsp[0].ival); }
#line 4624 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 313:
#line 1426 "./mdlparse.y" /* yacc.c:1646  */
    { (yyval.ival) = (int) (yyvsp[0].dbl); }
#line 4630 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 314:
#line 1427 "./mdlparse.y" /* yacc.c:1646  */
    { (yyval.ival) = INT_MAX; }
#line 4636 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 315:
#line 1434 "./mdlparse.y" /* yacc.c:1646  */
    { parse_state->current_object = parse_state->vol->root_instance; }
#line 4642 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 316:
#line 1435 "./mdlparse.y" /* yacc.c:1646  */
    {
                                                        check_regions(parse_state->vol->root_instance, (yyvsp[0].obj));
                                                        add_child_objects(parse_state->vol->root_instance, (yyvsp[0].obj), (yyvsp[0].obj));
                                                        parse_state->current_object = parse_state->vol->root_object;
                                                      }
#line 4652 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 317:
#line 1445 "./mdlparse.y" /* yacc.c:1646  */
    { add_child_objects(parse_state->vol->root_object, (yyvsp[0].obj), (yyvsp[0].obj)); }
#line 4658 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 323:
#line 1461 "./mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.sym) = mdl_start_object(parse_state, (yyvsp[0].str))); }
#line 4664 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 325:
#line 1467 "./mdlparse.y" /* yacc.c:1646  */
    { mdl_finish_object(parse_state); }
#line 4670 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 329:
#line 1480 "./mdlparse.y" /* yacc.c:1646  */
    { transform_translate(parse_state->vol, parse_state->current_object->t_matrix, (yyvsp[0].vec3)); }
#line 4676 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 330:
#line 1481 "./mdlparse.y" /* yacc.c:1646  */
    { transform_scale(parse_state->current_object->t_matrix, (yyvsp[0].vec3)); }
#line 4682 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 331:
#line 1482 "./mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_transform_rotate(parse_state, parse_state->current_object->t_matrix, (yyvsp[-2].vec3), (yyvsp[0].dbl))); }
#line 4688 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 332:
#line 1491 "./mdlparse.y" /* yacc.c:1646  */
    {
                                                          struct object *the_object = (struct object *) (yyvsp[-5].sym)->value;
                                                          the_object->object_type = META_OBJ;
                                                          add_child_objects(the_object, (yyvsp[-2].obj_list).obj_head, (yyvsp[-2].obj_list).obj_tail);
                                                          (yyval.obj) = the_object;
                                                      }
#line 4699 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 333:
#line 1500 "./mdlparse.y" /* yacc.c:1646  */
    { mdl_object_list_singleton(& (yyval.obj_list), (yyvsp[0].obj)); }
#line 4705 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 334:
#line 1501 "./mdlparse.y" /* yacc.c:1646  */
    { (yyval.obj_list) = (yyvsp[-1].obj_list); mdl_add_object_to_list(& (yyval.obj_list), (yyvsp[0].obj)); }
#line 4711 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 337:
#line 1510 "./mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_deep_copy_object(parse_state, (struct object *) (yyvsp[-3].sym)->value, (struct object *) (yyvsp[-1].sym)->value)); }
#line 4717 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 338:
#line 1512 "./mdlparse.y" /* yacc.c:1646  */
    { (yyval.obj) = (struct object *) (yyvsp[-6].sym)->value; }
#line 4723 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 341:
#line 1522 "./mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_start_release_site(parse_state, (yyvsp[-2].sym), SHAPE_UNDEFINED)); }
#line 4729 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 342:
#line 1526 "./mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.obj) = mdl_finish_release_site(parse_state, (yyvsp[-7].sym))); }
#line 4735 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 343:
#line 1529 "./mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_release_site_geometry_region(parse_state, parse_state->current_release_site, parse_state->current_object, (yyvsp[0].rev))); }
#line 4741 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 344:
#line 1530 "./mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_release_site_geometry_object(parse_state, parse_state->current_release_site, (struct object *) (yyvsp[0].sym)->value)); }
#line 4747 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 345:
#line 1531 "./mdlparse.y" /* yacc.c:1646  */
    { parse_state->current_release_site->release_shape = SHAPE_SPHERICAL; }
#line 4753 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 346:
#line 1532 "./mdlparse.y" /* yacc.c:1646  */
    { parse_state->current_release_site->release_shape = SHAPE_CUBIC; }
#line 4759 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 347:
#line 1533 "./mdlparse.y" /* yacc.c:1646  */
    { parse_state->current_release_site->release_shape = SHAPE_ELLIPTIC; }
#line 4765 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 348:
#line 1534 "./mdlparse.y" /* yacc.c:1646  */
    { parse_state->current_release_site->release_shape = SHAPE_RECTANGULAR; }
#line 4771 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 349:
#line 1535 "./mdlparse.y" /* yacc.c:1646  */
    { parse_state->current_release_site->release_shape = SHAPE_SPHERICAL_SHELL; }
#line 4777 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 350:
#line 1536 "./mdlparse.y" /* yacc.c:1646  */
    {
                                                          parse_state->current_release_site->release_shape = SHAPE_LIST;
                                                          parse_state->current_release_site->release_number_method = CONSTNUM;
                                                      }
#line 4786 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 351:
#line 1543 "./mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.rev) = new_release_region_expr_term((yyvsp[0].sym))); }
#line 4792 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 352:
#line 1544 "./mdlparse.y" /* yacc.c:1646  */
    { (yyval.rev) = (yyvsp[-1].rev); }
#line 4798 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 353:
#line 1545 "./mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.rev) = new_release_region_expr_binary((yyvsp[-2].rev), (yyvsp[0].rev), REXP_UNION)); }
#line 4804 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 354:
#line 1546 "./mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.rev) = new_release_region_expr_binary((yyvsp[-2].rev), (yyvsp[0].rev), REXP_SUBTRACTION)); }
#line 4810 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 355:
#line 1547 "./mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.rev) = new_release_region_expr_binary((yyvsp[-2].rev), (yyvsp[0].rev), REXP_INTERSECTION)); }
#line 4816 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 356:
#line 1552 "./mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_start_release_site(parse_state, (yyvsp[-2].sym), (yyvsp[-1].tok))); }
#line 4822 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 357:
#line 1555 "./mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.obj) = mdl_finish_release_site(parse_state, (yyvsp[-6].sym))); }
#line 4828 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 358:
#line 1558 "./mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = SHAPE_SPHERICAL; }
#line 4834 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 359:
#line 1559 "./mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = SHAPE_CUBIC; }
#line 4840 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 360:
#line 1560 "./mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = SHAPE_ELLIPTIC; }
#line 4846 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 361:
#line 1561 "./mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = SHAPE_RECTANGULAR; }
#line 4852 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 362:
#line 1562 "./mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = SHAPE_SPHERICAL_SHELL; }
#line 4858 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 365:
#line 1570 "./mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.sym) = mdl_existing_num_or_array(parse_state, (yyvsp[0].str))); }
#line 4864 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 366:
#line 1574 "./mdlparse.y" /* yacc.c:1646  */
    { set_release_site_location(parse_state->vol, parse_state->current_release_site, (yyvsp[0].vec3)); }
#line 4870 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 367:
#line 1575 "./mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_release_site_molecule(parse_state, parse_state->current_release_site, & (yyvsp[0].mol_type))); }
#line 4876 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 368:
#line 1576 "./mdlparse.y" /* yacc.c:1646  */
    {
                                                        if (parse_state->current_release_site->release_shape == SHAPE_LIST)
                                                        {
                                                          mdlerror(parse_state, "molecules are already specified in a list--cannot set number or density.");
                                                          return 1;
                                                        }
                                                      }
#line 4888 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 369:
#line 1583 "./mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_release_site_diameter(parse_state, parse_state->current_release_site, (yyvsp[0].dbl) * (((yyvsp[-2].tok) == SITE_RADIUS) ? 2.0 : 1.0))); }
#line 4894 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 370:
#line 1584 "./mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_release_site_diameter_array(parse_state, parse_state->current_release_site, (yyvsp[0].nlist).value_count, (yyvsp[0].nlist).value_head, ((yyvsp[-2].tok) == SITE_RADIUS) ? 2.0 : 1.0)); }
#line 4900 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 371:
#line 1585 "./mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_release_site_diameter_var(parse_state, parse_state->current_release_site, ((yyvsp[-2].tok) == SITE_RADIUS) ? 2.0 : 1.0, (yyvsp[0].sym))); }
#line 4906 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 372:
#line 1586 "./mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_release_site_probability(parse_state, parse_state->current_release_site, (yyvsp[0].dbl))); }
#line 4912 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 373:
#line 1588 "./mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_release_site_pattern(parse_state, parse_state->current_release_site, (yyvsp[0].sym))); }
#line 4918 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 374:
#line 1590 "./mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_release_site_molecule_positions(parse_state, parse_state->current_release_site, & (yyvsp[-1].rsm_list))); }
#line 4924 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 375:
#line 1594 "./mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = SITE_DIAMETER; }
#line 4930 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 376:
#line 1595 "./mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = SITE_RADIUS; }
#line 4936 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 381:
#line 1607 "./mdlparse.y" /* yacc.c:1646  */
    { set_release_site_constant_number(parse_state->current_release_site, (yyvsp[0].dbl)); }
#line 4942 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 382:
#line 1610 "./mdlparse.y" /* yacc.c:1646  */
    { set_release_site_constant_number(parse_state->current_release_site, (yyvsp[-1].dbl)); }
#line 4948 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 383:
#line 1617 "./mdlparse.y" /* yacc.c:1646  */
    { set_release_site_gaussian_number(parse_state->current_release_site, (yyvsp[-4].dbl), (yyvsp[-1].dbl)); }
#line 4954 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 384:
#line 1625 "./mdlparse.y" /* yacc.c:1646  */
    { set_release_site_volume_dependent_number(parse_state->current_release_site, (yyvsp[-7].dbl), (yyvsp[-4].dbl), (yyvsp[-1].dbl)); }
#line 4960 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 385:
#line 1629 "./mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_release_site_concentration(parse_state, parse_state->current_release_site, (yyvsp[0].dbl))); }
#line 4966 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 386:
#line 1630 "./mdlparse.y" /* yacc.c:1646  */
    { CHECK(set_release_site_density(parse_state->current_release_site, (yyvsp[0].dbl))); }
#line 4972 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 387:
#line 1634 "./mdlparse.y" /* yacc.c:1646  */
    { release_single_molecule_singleton(& (yyval.rsm_list), (yyvsp[0].rsm)); }
#line 4978 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 388:
#line 1636 "./mdlparse.y" /* yacc.c:1646  */
    { (yyval.rsm_list) = (yyvsp[-1].rsm_list); add_release_single_molecule_to_list(& (yyval.rsm_list), (yyvsp[0].rsm)); }
#line 4984 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 389:
#line 1640 "./mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.rsm) = mdl_new_release_single_molecule(parse_state, &(yyvsp[-1].mol_type), (yyvsp[0].vec3))); }
#line 4990 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 391:
#line 1651 "./mdlparse.y" /* yacc.c:1646  */
    {
                                                        CHECKN((yyval.obj) = mdl_new_polygon_list(
                                                          parse_state, (yyvsp[-4].str), (yyvsp[-1].vertlist).vertex_count, (yyvsp[-1].vertlist).vertex_head,
                                                          (yyvsp[0].ecl).connection_count, (yyvsp[0].ecl).connection_head));
                                                      }
#line 5000 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 392:
#line 1660 "./mdlparse.y" /* yacc.c:1646  */
    {
                                                          (yyval.obj) = (struct object *) (yyvsp[-3].obj);
                                                          CHECK(mdl_finish_polygon_list(parse_state, (yyval.obj)));
                                                      }
#line 5009 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 393:
#line 1666 "./mdlparse.y" /* yacc.c:1646  */
    { (yyval.vertlist) = (yyvsp[-1].vertlist); }
#line 5015 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 394:
#line 1669 "./mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.vertlistitem) = mdl_new_vertex_list_item((yyvsp[0].vec3))); }
#line 5021 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 395:
#line 1672 "./mdlparse.y" /* yacc.c:1646  */
    { mdl_vertex_list_singleton(& (yyval.vertlist), (yyvsp[0].vertlistitem)); }
#line 5027 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 396:
#line 1673 "./mdlparse.y" /* yacc.c:1646  */
    { (yyval.vertlist) = (yyvsp[-1].vertlist); mdl_add_vertex_to_list(& (yyval.vertlist), (yyvsp[0].vertlistitem)); }
#line 5033 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 397:
#line 1678 "./mdlparse.y" /* yacc.c:1646  */
    { (yyval.ecl) = (yyvsp[-1].ecl); }
#line 5039 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 398:
#line 1682 "./mdlparse.y" /* yacc.c:1646  */
    { mdl_element_connection_list_singleton(& (yyval.ecl), (yyvsp[0].elem_conn)); }
#line 5045 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 399:
#line 1684 "./mdlparse.y" /* yacc.c:1646  */
    { (yyval.ecl) = (yyvsp[-1].ecl); mdl_add_element_connection_to_list(& (yyval.ecl), (yyvsp[0].elem_conn)); }
#line 5051 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 400:
#line 1687 "./mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.elem_conn) = mdl_new_element_connection(parse_state, & (yyvsp[0].nlist))); }
#line 5057 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 405:
#line 1703 "./mdlparse.y" /* yacc.c:1646  */
    { CHECKN(parse_state->current_region = mdl_get_region(parse_state, parse_state->current_object, "REMOVED")); }
#line 5063 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 406:
#line 1705 "./mdlparse.y" /* yacc.c:1646  */
    {
                                                          parse_state->current_region->element_list_head = (yyvsp[-1].elem_list).elml_head;
                                                          if (parse_state->current_object->object_type == POLY_OBJ)
                                                          {
                                                            CHECK(mdl_normalize_elements(parse_state, parse_state->current_region,0));
                                                          }
                                                      }
#line 5075 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 409:
#line 1719 "./mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = Z_POS; }
#line 5081 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 410:
#line 1720 "./mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = Z_NEG; }
#line 5087 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 411:
#line 1721 "./mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = Y_NEG; }
#line 5093 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 412:
#line 1722 "./mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = Y_POS; }
#line 5099 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 413:
#line 1723 "./mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = X_NEG; }
#line 5105 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 414:
#line 1724 "./mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = X_POS; }
#line 5111 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 415:
#line 1725 "./mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = ALL_SIDES; }
#line 5117 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 417:
#line 1731 "./mdlparse.y" /* yacc.c:1646  */
    { (yyval.elem_list) = (yyvsp[-1].elem_list); mdl_add_elements_to_list(& (yyval.elem_list), (yyvsp[0].elem_list).elml_head, (yyvsp[0].elem_list).elml_tail); }
#line 5123 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 420:
#line 1737 "./mdlparse.y" /* yacc.c:1646  */
    { (yyval.elem_list).elml_tail = (yyval.elem_list).elml_head = (yyvsp[0].elem_list_item); }
#line 5129 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 421:
#line 1738 "./mdlparse.y" /* yacc.c:1646  */
    { (yyval.elem_list).elml_tail = (yyval.elem_list).elml_head = (yyvsp[0].elem_list_item); }
#line 5135 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 422:
#line 1743 "./mdlparse.y" /* yacc.c:1646  */
    { (yyval.elem_list) = (yyvsp[-1].elem_list); }
#line 5141 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 423:
#line 1748 "./mdlparse.y" /* yacc.c:1646  */
    { (yyval.elem_list) = (yyvsp[-1].elem_list); mdl_set_elements_to_exclude((yyval.elem_list).elml_head); }
#line 5147 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 425:
#line 1755 "./mdlparse.y" /* yacc.c:1646  */
    { (yyval.elem_list).elml_tail = (yyval.elem_list).elml_head = (yyvsp[0].elem_list_item); }
#line 5153 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 426:
#line 1756 "./mdlparse.y" /* yacc.c:1646  */
    { (yyval.elem_list) = (yyvsp[-2].elem_list); mdl_add_elements_to_list(& (yyval.elem_list), (yyvsp[0].elem_list_item), (yyvsp[0].elem_list_item)); }
#line 5159 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 427:
#line 1759 "./mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.elem_list_item) = new_element_list((unsigned int) (yyvsp[0].dbl), (unsigned int) (yyvsp[0].dbl))); }
#line 5165 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 428:
#line 1760 "./mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.elem_list_item) = new_element_list((unsigned int) (yyvsp[-2].dbl), (unsigned int) (yyvsp[0].dbl))); }
#line 5171 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 429:
#line 1761 "./mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.elem_list_item) = mdl_new_element_side(parse_state, (yyvsp[0].tok))); }
#line 5177 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 430:
#line 1764 "./mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.elem_list_item) = mdl_new_element_previous_region(parse_state, parse_state->current_object, parse_state->current_region, (yyvsp[0].str), (yyvsp[-2].tok))); }
#line 5183 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 431:
#line 1767 "./mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = 0; }
#line 5189 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 432:
#line 1768 "./mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = 1; }
#line 5195 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 433:
#line 1771 "./mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.elem_list_item) = mdl_new_element_patch(parse_state, parse_state->current_polygon, (yyvsp[-2].vec3), (yyvsp[0].vec3), (yyvsp[-4].tok))); }
#line 5201 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 434:
#line 1774 "./mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = 0; }
#line 5207 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 435:
#line 1775 "./mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = 1; }
#line 5213 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 439:
#line 1791 "./mdlparse.y" /* yacc.c:1646  */
    { parse_state->current_region = (yyvsp[-1].reg); }
#line 5219 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 440:
#line 1792 "./mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_region_elements(parse_state, (yyvsp[-3].reg), (yyvsp[0].elem_list).elml_head, (yyvsp[-3].reg)->parent->object_type == POLY_OBJ)); }
#line 5225 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 441:
#line 1794 "./mdlparse.y" /* yacc.c:1646  */
    { parse_state->current_region = NULL; }
#line 5231 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 442:
#line 1802 "./mdlparse.y" /* yacc.c:1646  */
    {
                                                        CHECKN(mdl_new_voxel_list(parse_state, (yyvsp[-4].sym),
                                                                                  (yyvsp[-1].vertlist).vertex_count, (yyvsp[-1].vertlist).vertex_head,
                                                                                  (yyvsp[0].ecl).connection_count, (yyvsp[0].ecl).connection_head));
                                                      }
#line 5241 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 443:
#line 1808 "./mdlparse.y" /* yacc.c:1646  */
    { (yyval.obj) = (struct object *) (yyvsp[-7].sym)->value; }
#line 5247 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 444:
#line 1813 "./mdlparse.y" /* yacc.c:1646  */
    { (yyval.ecl) = (yyvsp[-1].ecl); }
#line 5253 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 445:
#line 1816 "./mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.elem_conn) = mdl_new_tet_element_connection(parse_state, & (yyvsp[0].nlist))); }
#line 5259 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 446:
#line 1820 "./mdlparse.y" /* yacc.c:1646  */
    {
                                                          (yyval.ecl).connection_head = (yyval.ecl).connection_tail = (yyvsp[0].elem_conn);
                                                          (yyval.ecl).connection_count = 1;
                                                      }
#line 5268 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 447:
#line 1824 "./mdlparse.y" /* yacc.c:1646  */
    {
                                                          (yyval.ecl) = (yyvsp[-1].ecl);
                                                          (yyval.ecl).connection_tail = (yyval.ecl).connection_tail->next = (yyvsp[0].elem_conn);
                                                          ++ (yyval.ecl).connection_count;
                                                      }
#line 5278 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 448:
#line 1835 "./mdlparse.y" /* yacc.c:1646  */
    { CHECKN(mdl_new_box_object(parse_state, (yyvsp[-8].sym), (yyvsp[-3].vec3), (yyvsp[-1].vec3))); }
#line 5284 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 449:
#line 1836 "./mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_triangulate_box_object(parse_state, (yyvsp[-10].sym), parse_state->current_polygon, (yyvsp[-2].dbl))); }
#line 5290 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 450:
#line 1838 "./mdlparse.y" /* yacc.c:1646  */
    {
                                                          CHECK(mdl_finish_box_object(parse_state, (yyvsp[-13].sym)));
                                                          (yyval.obj) = (struct object *) (yyvsp[-13].sym)->value;
                                                      }
#line 5299 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 451:
#line 1844 "./mdlparse.y" /* yacc.c:1646  */
    { (yyval.dbl) = 0.0; }
#line 5305 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 452:
#line 1845 "./mdlparse.y" /* yacc.c:1646  */
    {
                                                        (yyval.dbl) = (yyvsp[0].dbl);
                                                        if ((yyval.dbl) < 2.0)
                                                        {
                                                          mdlerror(parse_state, "invalid aspect ratio requested (must be greater than or equal to 2.0)");
                                                          return 1;
                                                        }
                                                      }
#line 5318 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 456:
#line 1871 "./mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_start_existing_obj_region_def(parse_state, (yyvsp[0].sym))); }
#line 5324 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 457:
#line 1872 "./mdlparse.y" /* yacc.c:1646  */
    { parse_state->current_region = (yyvsp[-1].reg); }
#line 5330 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 458:
#line 1874 "./mdlparse.y" /* yacc.c:1646  */
    { mdl_set_region_elements(parse_state, (yyvsp[-4].reg), (yyvsp[0].elem_list).elml_head, 1); }
#line 5336 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 459:
#line 1876 "./mdlparse.y" /* yacc.c:1646  */
    {
                                                          parse_state->current_region = NULL;
                                                          parse_state->current_polygon = NULL;
                                                          parse_state->current_object = parse_state->vol->root_object;
                                                      }
#line 5346 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 460:
#line 1883 "./mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.reg) = mdl_create_region(parse_state, parse_state->current_object, (yyvsp[0].str))); }
#line 5352 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 464:
#line 1894 "./mdlparse.y" /* yacc.c:1646  */
    { mdl_add_surf_mol_to_region(parse_state->current_region, & (yyvsp[0].surf_mol_dat_list)); }
#line 5358 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 465:
#line 1898 "./mdlparse.y" /* yacc.c:1646  */
    { mdl_set_region_surface_class(parse_state, parse_state->current_region, (yyvsp[0].sym)); }
#line 5364 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 469:
#line 1917 "./mdlparse.y" /* yacc.c:1646  */
    { parse_state->current_region = (struct region *) (yyvsp[-1].sym)->value; }
#line 5370 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 470:
#line 1919 "./mdlparse.y" /* yacc.c:1646  */
    { parse_state->current_region = NULL; }
#line 5376 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 471:
#line 1927 "./mdlparse.y" /* yacc.c:1646  */
    {
                                                          parse_state->header_comment = NULL;  /* No header by default */
                                                          parse_state->exact_time_flag = 1;    /* Print exact_time column in TRIGGER output by default */
                                                      }
#line 5385 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 472:
#line 1933 "./mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_add_reaction_output_block_to_world(parse_state, (int) (yyvsp[-4].dbl), & (yyvsp[-2].ro_otimes), & (yyvsp[-1].ro_sets))); }
#line 5391 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 473:
#line 1937 "./mdlparse.y" /* yacc.c:1646  */
    { (yyval.dbl) = COUNTBUFFERSIZE; }
#line 5397 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 474:
#line 1938 "./mdlparse.y" /* yacc.c:1646  */
    {
                                                          double temp_value = (yyvsp[0].dbl);
                                                          if (!(temp_value >= 1.0 && temp_value < UINT_MAX))
                                                          {
                                                            mdlerror_fmt(parse_state, "Requested buffer size of %.15g lines is invalid.  Suggested range is 100-1000000.", temp_value);
                                                            return 1;
                                                          }
                                                          (yyval.dbl) = (yyvsp[0].dbl);
                                                      }
#line 5411 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 478:
#line 1954 "./mdlparse.y" /* yacc.c:1646  */
    { (yyval.ro_otimes).type = OUTPUT_BY_STEP; (yyval.ro_otimes).step = (yyvsp[0].dbl); }
#line 5417 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 479:
#line 1958 "./mdlparse.y" /* yacc.c:1646  */
    {
                                                        (yyval.ro_otimes).type = OUTPUT_BY_ITERATION_LIST;
                                                        (yyval.ro_otimes).values = (yyvsp[0].nlist);
                                                      }
#line 5426 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 480:
#line 1966 "./mdlparse.y" /* yacc.c:1646  */
    {
                                                        (yyval.ro_otimes).type = OUTPUT_BY_TIME_LIST;
                                                        (yyval.ro_otimes).values = (yyvsp[0].nlist);
                                                      }
#line 5435 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 481:
#line 1973 "./mdlparse.y" /* yacc.c:1646  */
    { (yyval.ro_sets).set_head = (yyval.ro_sets).set_tail = (yyvsp[0].ro_set); }
#line 5441 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 482:
#line 1975 "./mdlparse.y" /* yacc.c:1646  */
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
#line 5456 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 484:
#line 1989 "./mdlparse.y" /* yacc.c:1646  */
    { (yyval.ro_set) = NULL; }
#line 5462 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 485:
#line 1990 "./mdlparse.y" /* yacc.c:1646  */
    { (yyval.ro_set) = NULL; }
#line 5468 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 486:
#line 1994 "./mdlparse.y" /* yacc.c:1646  */
    {  parse_state->count_flags = 0; }
#line 5474 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 487:
#line 1996 "./mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.ro_set) = mdl_populate_output_set(parse_state, parse_state->header_comment, parse_state->exact_time_flag, (yyvsp[-3].ro_cols).column_head, (yyvsp[-1].tok), (yyvsp[0].str))); }
#line 5480 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 488:
#line 2000 "./mdlparse.y" /* yacc.c:1646  */
    { (yyval.str) = NULL; }
#line 5486 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 489:
#line 2001 "./mdlparse.y" /* yacc.c:1646  */
    { (yyval.str) = ((yyvsp[0].tok) ? "" : NULL); }
#line 5492 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 490:
#line 2002 "./mdlparse.y" /* yacc.c:1646  */
    { (yyval.str) = (yyvsp[0].str); }
#line 5498 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 491:
#line 2006 "./mdlparse.y" /* yacc.c:1646  */
    { parse_state->header_comment = (yyvsp[0].str); }
#line 5504 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 492:
#line 2010 "./mdlparse.y" /* yacc.c:1646  */
    { parse_state->exact_time_flag = (yyvsp[0].tok); }
#line 5510 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 494:
#line 2016 "./mdlparse.y" /* yacc.c:1646  */
    {
                                                          (yyval.ro_cols) = (yyvsp[-2].ro_cols);
                                                          (yyval.ro_cols).column_tail->next = (yyvsp[0].ro_cols).column_head;
                                                          (yyval.ro_cols).column_tail = (yyvsp[0].ro_cols).column_tail;
                                                      }
#line 5520 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 495:
#line 2024 "./mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_single_count_expr(parse_state, & (yyval.ro_cols), (yyvsp[-1].cnt), (yyvsp[0].str))); }
#line 5526 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 496:
#line 2028 "./mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.cnt) = mdl_new_oexpr_constant(parse_state, (yyvsp[0].dbl))); }
#line 5532 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 498:
#line 2030 "./mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.cnt) = mdl_join_oexpr_tree(parse_state, (yyvsp[-1].cnt), NULL, '(')); }
#line 5538 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 499:
#line 2031 "./mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.cnt) = mdl_join_oexpr_tree(parse_state, (yyvsp[-2].cnt),   (yyvsp[0].cnt), '+')); }
#line 5544 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 500:
#line 2032 "./mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.cnt) = mdl_join_oexpr_tree(parse_state, (yyvsp[-2].cnt),   (yyvsp[0].cnt), '-')); }
#line 5550 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 501:
#line 2033 "./mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.cnt) = mdl_join_oexpr_tree(parse_state, (yyvsp[-2].cnt),   (yyvsp[0].cnt), '*')); }
#line 5556 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 502:
#line 2034 "./mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.cnt) = mdl_join_oexpr_tree(parse_state, (yyvsp[-2].cnt),   (yyvsp[0].cnt), '/')); }
#line 5562 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 503:
#line 2035 "./mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.cnt) = mdl_join_oexpr_tree(parse_state, (yyvsp[0].cnt), NULL, '_')); }
#line 5568 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 504:
#line 2036 "./mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.cnt) = mdl_sum_oexpr((yyvsp[-1].cnt))); }
#line 5574 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 505:
#line 2041 "./mdlparse.y" /* yacc.c:1646  */
    { parse_state->count_flags |= COUNT_PRESENT; }
#line 5580 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 506:
#line 2042 "./mdlparse.y" /* yacc.c:1646  */
    { (yyval.cnt) = (yyvsp[-1].cnt); }
#line 5586 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 507:
#line 2043 "./mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.cnt) = mdl_new_oexpr_constant(parse_state, (yyvsp[-1].dbl))); }
#line 5592 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 508:
#line 2044 "./mdlparse.y" /* yacc.c:1646  */
    { parse_state->count_flags |= TRIGGER_PRESENT; }
#line 5598 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 509:
#line 2045 "./mdlparse.y" /* yacc.c:1646  */
    { (yyval.cnt) = (yyvsp[-1].cnt); }
#line 5604 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 510:
#line 2048 "./mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = FILE_OVERWRITE; }
#line 5610 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 511:
#line 2049 "./mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = FILE_SUBSTITUTE; }
#line 5616 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 512:
#line 2050 "./mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = FILE_APPEND; }
#line 5622 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 513:
#line 2051 "./mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = FILE_APPEND_HEADER; }
#line 5628 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 514:
#line 2052 "./mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = FILE_CREATE; }
#line 5634 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 516:
#line 2058 "./mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.sym) = mdl_existing_rxn_pathname_or_molecule(parse_state, (yyvsp[0].str))); }
#line 5640 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 517:
#line 2062 "./mdlparse.y" /* yacc.c:1646  */
    {
                                                        (yyval.mol_type) = (yyvsp[0].mol_type);
                                                        if ((yyval.mol_type).orient > 0)
                                                          (yyval.mol_type).orient = 1;
                                                        else if ((yyval.mol_type).orient < 0)
                                                          (yyval.mol_type).orient = -1;
                                                        CHECKN((yyval.mol_type).mol_type = mdl_existing_molecule(parse_state, (yyvsp[-1].str)));
                                                      }
#line 5653 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 521:
#line 2079 "./mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.cnt) = mdl_count_syntax_1(parse_state, (yyvsp[-3].sym), (yyvsp[-1].sym), (yyvsp[0].tok), parse_state->count_flags)); }
#line 5659 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 522:
#line 2084 "./mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.cnt) = mdl_count_syntax_2(parse_state, (yyvsp[-3].mol_type).mol_type, (yyvsp[-3].mol_type).orient, (yyvsp[-1].sym), (yyvsp[0].tok), parse_state->count_flags)); }
#line 5665 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 523:
#line 2089 "./mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.cnt) = mdl_count_syntax_3(parse_state, (yyvsp[-3].str), (yyvsp[-1].sym), (yyvsp[0].tok), parse_state->count_flags)); }
#line 5671 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 524:
#line 2092 "./mdlparse.y" /* yacc.c:1646  */
    { (yyval.sym) = NULL; }
#line 5677 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 525:
#line 2093 "./mdlparse.y" /* yacc.c:1646  */
    { (yyval.sym) = (yyvsp[0].sym); }
#line 5683 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 526:
#line 2094 "./mdlparse.y" /* yacc.c:1646  */
    { (yyval.sym) = (yyvsp[0].sym); }
#line 5689 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 527:
#line 2097 "./mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = REPORT_NOTHING; }
#line 5695 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 528:
#line 2098 "./mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = (yyvsp[0].tok); }
#line 5701 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 529:
#line 2101 "./mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = REPORT_FRONT_HITS; }
#line 5707 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 530:
#line 2102 "./mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = REPORT_BACK_HITS; }
#line 5713 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 531:
#line 2103 "./mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = REPORT_ALL_HITS; }
#line 5719 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 532:
#line 2104 "./mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = REPORT_FRONT_CROSSINGS; }
#line 5725 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 533:
#line 2105 "./mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = REPORT_BACK_CROSSINGS; }
#line 5731 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 534:
#line 2106 "./mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = REPORT_ALL_CROSSINGS; }
#line 5737 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 535:
#line 2107 "./mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = REPORT_CONCENTRATION; }
#line 5743 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 536:
#line 2108 "./mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = REPORT_ENCLOSED; }
#line 5749 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 537:
#line 2111 "./mdlparse.y" /* yacc.c:1646  */
    { (yyval.str) = NULL; }
#line 5755 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 538:
#line 2112 "./mdlparse.y" /* yacc.c:1646  */
    { (yyval.str) = (yyvsp[0].str); }
#line 5761 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 539:
#line 2119 "./mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_new_viz_output_block(parse_state)); }
#line 5767 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 540:
#line 2122 "./mdlparse.y" /* yacc.c:1646  */
    { }
#line 5773 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 543:
#line 2131 "./mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_viz_mode(parse_state->vol->viz_blocks, CELLBLENDER_MODE)); }
#line 5779 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 544:
#line 2132 "./mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_viz_mode(parse_state->vol->viz_blocks, (yyvsp[0].ival))); }
#line 5785 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 545:
#line 2135 "./mdlparse.y" /* yacc.c:1646  */
    { (yyval.ival) = NO_VIZ_MODE; }
#line 5791 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 546:
#line 2136 "./mdlparse.y" /* yacc.c:1646  */
    { (yyval.ival) = ASCII_MODE; }
#line 5797 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 547:
#line 2137 "./mdlparse.y" /* yacc.c:1646  */
    { (yyval.ival) = CELLBLENDER_MODE; }
#line 5803 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 549:
#line 2142 "./mdlparse.y" /* yacc.c:1646  */
    {
                                                        if ((yyvsp[0].frame_list).frame_head)
                                                        {
                                                          (yyvsp[0].frame_list).frame_tail->next = parse_state->vol->viz_blocks->frame_data_head;
                                                          parse_state->vol->viz_blocks->frame_data_head = (yyvsp[0].frame_list).frame_head;
                                                        }
                                                      }
#line 5815 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 551:
#line 2155 "./mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_viz_filename_prefix(parse_state, parse_state->vol->viz_blocks, (yyvsp[0].str))); }
#line 5821 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 552:
#line 2161 "./mdlparse.y" /* yacc.c:1646  */
    { (yyval.frame_list) = (yyvsp[-1].frame_list); }
#line 5827 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 554:
#line 2167 "./mdlparse.y" /* yacc.c:1646  */
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
#line 5843 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 555:
#line 2181 "./mdlparse.y" /* yacc.c:1646  */
    { (yyval.frame_list).frame_head = (yyval.frame_list).frame_tail = NULL; }
#line 5849 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 559:
#line 2193 "./mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_viz_state(parse_state, & (yyval.ival), (yyvsp[0].dbl))); }
#line 5855 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 560:
#line 2194 "./mdlparse.y" /* yacc.c:1646  */
    { (yyval.ival) = INCLUDE_OBJ; }
#line 5861 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 563:
#line 2204 "./mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_viz_include_molecules(parse_state, parse_state->vol->viz_blocks, (yyvsp[-1].symlist), (yyvsp[0].ival))); }
#line 5867 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 564:
#line 2205 "./mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_viz_include_all_molecules(parse_state->vol->viz_blocks, (yyvsp[0].ival))); }
#line 5873 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 565:
#line 2209 "./mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.symlist) = mdl_existing_molecule_list(parse_state, (yyvsp[0].str))); }
#line 5879 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 566:
#line 2210 "./mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.symlist) = mdl_existing_molecules_wildcard(parse_state, (yyvsp[0].str))); }
#line 5885 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 567:
#line 2214 "./mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_new_viz_all_times(parse_state, & (yyval.nlist))); }
#line 5891 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 569:
#line 2220 "./mdlparse.y" /* yacc.c:1646  */
    { (yyval.frame_list) = (yyvsp[-1].frame_list); }
#line 5897 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 571:
#line 2226 "./mdlparse.y" /* yacc.c:1646  */
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
#line 5915 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 572:
#line 2243 "./mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_new_viz_mol_frames(parse_state, parse_state->vol->viz_blocks, & (yyval.frame_list), OUTPUT_BY_TIME_LIST, (yyvsp[-2].tok), & (yyvsp[0].nlist))); }
#line 5921 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 573:
#line 2247 "./mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_new_viz_all_iterations(parse_state, & (yyval.nlist))); }
#line 5927 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 575:
#line 2254 "./mdlparse.y" /* yacc.c:1646  */
    { (yyval.frame_list) = (yyvsp[-1].frame_list); }
#line 5933 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 577:
#line 2260 "./mdlparse.y" /* yacc.c:1646  */
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
#line 5951 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 578:
#line 2277 "./mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_new_viz_mol_frames(parse_state, parse_state->vol->viz_blocks, & (yyval.frame_list), OUTPUT_BY_ITERATION_LIST, (yyvsp[-2].tok), & (yyvsp[0].nlist))); }
#line 5957 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 579:
#line 2280 "./mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = ALL_MOL_DATA; }
#line 5963 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 580:
#line 2281 "./mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = MOL_POS; }
#line 5969 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 581:
#line 2282 "./mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = MOL_ORIENT; }
#line 5975 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 582:
#line 2296 "./mdlparse.y" /* yacc.c:1646  */
    {
                                                          struct volume_output_item *vo;
                                                          CHECKN(vo = mdl_new_volume_output_item(parse_state, (yyvsp[-6].str), & (yyvsp[-5].species_lst), (yyvsp[-4].vec3), (yyvsp[-3].vec3), (yyvsp[-2].vec3), (yyvsp[-1].otimes)));
                                                          vo->next = parse_state->vol->volume_output_head;
                                                          parse_state->vol->volume_output_head = vo;
                                                      }
#line 5986 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 583:
#line 2305 "./mdlparse.y" /* yacc.c:1646  */
    { (yyval.str) = (yyvsp[0].str); }
#line 5992 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 585:
#line 2311 "./mdlparse.y" /* yacc.c:1646  */
    {
                                                          (yyval.species_lst) = (yyvsp[-1].species_lst);
                                                          (yyval.species_lst).species_count += (yyvsp[0].species_lst).species_count;
                                                          (yyval.species_lst).species_tail->next = (yyvsp[0].species_lst).species_head;
                                                          (yyval.species_lst).species_tail = (yyvsp[0].species_lst).species_tail;
                                                      }
#line 6003 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 586:
#line 2320 "./mdlparse.y" /* yacc.c:1646  */
    { (yyval.species_lst) = (yyvsp[0].species_lst); }
#line 6009 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 587:
#line 2323 "./mdlparse.y" /* yacc.c:1646  */
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
#line 6029 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 588:
#line 2341 "./mdlparse.y" /* yacc.c:1646  */
    { (yyval.species_lst).species_tail = (yyval.species_lst).species_head = (yyvsp[0].species_lst_item); (yyval.species_lst).species_count = 1; }
#line 6035 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 589:
#line 2343 "./mdlparse.y" /* yacc.c:1646  */
    {
                                                        (yyval.species_lst) = (yyvsp[-2].species_lst);
                                                        (yyval.species_lst).species_tail = (yyval.species_lst).species_tail->next = (yyvsp[0].species_lst_item);
                                                        ++ (yyval.species_lst).species_count;
                                                      }
#line 6045 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 590:
#line 2351 "./mdlparse.y" /* yacc.c:1646  */
    { (yyval.vec3) = (yyvsp[0].vec3); }
#line 6051 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 591:
#line 2355 "./mdlparse.y" /* yacc.c:1646  */
    { (yyval.vec3) = (yyvsp[0].vec3); }
#line 6057 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 592:
#line 2359 "./mdlparse.y" /* yacc.c:1646  */
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
#line 6080 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 593:
#line 2380 "./mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.otimes) = mdl_new_output_times_default(parse_state)); }
#line 6086 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 594:
#line 2381 "./mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.otimes) = mdl_new_output_times_step(parse_state, (yyvsp[0].dbl))); }
#line 6092 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 595:
#line 2382 "./mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.otimes) = mdl_new_output_times_iterations(parse_state, & (yyvsp[0].nlist))); }
#line 6098 "mdlparse.c" /* yacc.c:1646  */
    break;

  case 596:
#line 2383 "./mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.otimes) = mdl_new_output_times_time(parse_state, & (yyvsp[0].nlist))); }
#line 6104 "mdlparse.c" /* yacc.c:1646  */
    break;


#line 6108 "mdlparse.c" /* yacc.c:1646  */
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
#line 2386 "./mdlparse.y" /* yacc.c:1906  */






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
