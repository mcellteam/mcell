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
#line 68 "mdlparse.y" /* yacc.c:1909  */

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


#line 608 "mdlparse.h" /* yacc.c:1909  */
};
# define YYSTYPE_IS_TRIVIAL 1
# define YYSTYPE_IS_DECLARED 1
#endif



int mdlparse (struct mdlparse_vars *parse_state, yyscan_t scanner);

#endif /* !YY_MDL_MDLPARSE_H_INCLUDED  */
