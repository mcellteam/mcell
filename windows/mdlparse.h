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
typedef union YYSTYPE YYSTYPE;
union YYSTYPE
{
#line 69 "../src/mdlparse.y" /* yacc.c:1909  */

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

#line 655 "mdlparse.h" /* yacc.c:1909  */
};
# define YYSTYPE_IS_TRIVIAL 1
# define YYSTYPE_IS_DECLARED 1
#endif



int mdlparse (struct mdlparse_vars *parse_state, yyscan_t scanner);

#endif /* !YY_MDL_MDLPARSE_H_INCLUDED  */
