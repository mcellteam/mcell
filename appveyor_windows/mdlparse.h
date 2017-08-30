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
#line 67 "../src/mdlparse.y" /* yacc.c:1909  */

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


#line 632 "mdlparse.h" /* yacc.c:1909  */
};

typedef union YYSTYPE YYSTYPE;
# define YYSTYPE_IS_TRIVIAL 1
# define YYSTYPE_IS_DECLARED 1
#endif



int mdlparse (struct mdlparse_vars *parse_state, yyscan_t scanner);

#endif /* !YY_MDL_MDLPARSE_H_INCLUDED  */
