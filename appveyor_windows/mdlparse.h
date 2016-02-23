/* A Bison parser, made by GNU Bison 3.0.2.  */

/* Bison interface for Yacc-like parsers in C

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
#line 68 "../src/mdlparse.y" /* yacc.c:1909  */

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


#line 626 "mdlparse.h" /* yacc.c:1909  */
};
# define YYSTYPE_IS_TRIVIAL 1
# define YYSTYPE_IS_DECLARED 1
#endif



int mdlparse (struct mdlparse_vars *parse_state, yyscan_t scanner);

#endif /* !YY_MDL_MDLPARSE_H_INCLUDED  */
