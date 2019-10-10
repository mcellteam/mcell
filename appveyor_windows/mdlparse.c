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
#line 1 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:339  */

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

#line 137 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:339  */

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
   by #include "mdlparse.h".  */
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
#line 67 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:355  */

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
struct geom_object *obj;
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


#line 503 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:355  */
};

typedef union YYSTYPE YYSTYPE;
# define YYSTYPE_IS_TRIVIAL 1
# define YYSTYPE_IS_DECLARED 1
#endif



int mdlparse (struct mdlparse_vars *parse_state, yyscan_t scanner);

#endif /* !YY_MDL_HOME_JCZECH_MCELL_BUILD_DEPS_MDLPARSE_H_INCLUDED  */

/* Copy the second part of user declarations.  */

#line 519 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:358  */

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
#define YYLAST   2815

/* YYNTOKENS -- Number of terminals.  */
#define YYNTOKENS  280
/* YYNNTS -- Number of nonterminals.  */
#define YYNNTS  290
/* YYNRULES -- Number of rules.  */
#define YYNRULES  624
/* YYNSTATES -- Number of states.  */
#define YYNSTATES  1232

/* YYTRANSLATE[YYX] -- Symbol number corresponding to YYX as returned
   by yylex, with out-of-bounds checking.  */
#define YYUNDEFTOK  2
#define YYMAXUTOK   514

#define YYTRANSLATE(YYX)                                                \
  ((unsigned int) (YYX) <= YYMAXUTOK ? yytranslate[YYX] : YYUNDEFTOK)

/* YYTRANSLATE[TOKEN-NUM] -- Symbol number corresponding to TOKEN-NUM
   as returned by yylex, without out-of-bounds checking.  */
static const yytype_uint16 yytranslate[] =
{
       0,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,   260,   271,
     275,   276,   264,   262,   272,   263,     2,   265,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,   261,   270,
     278,   259,   277,     2,   279,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,   268,     2,   269,   267,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,   273,     2,   274,     2,     2,     2,     2,
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
     255,   256,   257,   258,   266
};

#if YYDEBUG
  /* YYRLINE[YYN] -- Source line where rule number YYN was defined.  */
static const yytype_uint16 yyrline[] =
{
       0,   602,   602,   606,   607,   612,   613,   614,   615,   616,
     617,   618,   619,   620,   621,   622,   623,   624,   625,   626,
     627,   628,   629,   630,   631,   632,   637,   640,   643,   646,
     649,   652,   655,   656,   659,   660,   661,   662,   663,   664,
     667,   668,   669,   670,   674,   675,   676,   686,   698,   701,
     704,   716,   717,   730,   731,   737,   760,   761,   762,   763,
     766,   769,   772,   773,   784,   787,   790,   791,   794,   795,
     798,   799,   802,   803,   806,   810,   811,   812,   813,   814,
     815,   816,   817,   818,   819,   820,   821,   822,   823,   824,
     825,   826,   827,   828,   829,   830,   831,   832,   833,   834,
     835,   836,   837,   838,   839,   843,   844,   848,   849,   850,
     851,   854,   860,   861,   862,   863,   864,   865,   866,   869,
     873,   876,   879,   882,   885,   888,   889,   899,   900,   901,
     913,   917,   923,   928,   932,   941,   945,   946,   950,   951,
     952,   953,   954,   955,   956,   957,   958,   959,   960,   961,
     962,   963,   964,   965,   966,   970,   971,   975,   979,   980,
     987,   991,   992,   996,   997,   998,   999,  1000,  1001,  1002,
    1003,  1004,  1005,  1006,  1007,  1008,  1009,  1010,  1011,  1012,
    1013,  1017,  1018,  1019,  1025,  1026,  1027,  1028,  1029,  1033,
    1034,  1035,  1039,  1040,  1041,  1042,  1050,  1051,  1052,  1053,
    1054,  1055,  1056,  1057,  1058,  1059,  1060,  1061,  1062,  1063,
    1064,  1065,  1066,  1067,  1074,  1075,  1076,  1077,  1081,  1085,
    1086,  1087,  1094,  1095,  1098,  1102,  1106,  1107,  1111,  1120,
    1123,  1127,  1128,  1132,  1133,  1142,  1153,  1154,  1158,  1159,
    1169,  1170,  1173,  1177,  1181,  1194,  1195,  1200,  1204,  1210,
    1211,  1216,  1216,  1221,  1224,  1226,  1231,  1232,  1236,  1239,
    1245,  1250,  1251,  1252,  1255,  1256,  1259,  1263,  1267,  1274,
    1278,  1287,  1291,  1300,  1307,  1312,  1313,  1316,  1319,  1320,
    1323,  1324,  1325,  1328,  1333,  1339,  1340,  1341,  1342,  1345,
    1346,  1350,  1355,  1356,  1359,  1363,  1364,  1368,  1371,  1372,
    1375,  1376,  1380,  1381,  1384,  1395,  1412,  1413,  1414,  1418,
    1419,  1420,  1427,  1434,  1437,  1441,  1448,  1450,  1452,  1454,
    1456,  1460,  1461,  1468,  1468,  1479,  1482,  1483,  1484,  1485,
    1486,  1495,  1498,  1501,  1504,  1506,  1510,  1514,  1515,  1516,
    1521,  1534,  1535,  1538,  1539,  1544,  1543,  1550,  1551,  1556,
    1555,  1563,  1564,  1565,  1566,  1567,  1568,  1569,  1570,  1577,
    1578,  1579,  1580,  1581,  1586,  1585,  1592,  1593,  1594,  1595,
    1596,  1600,  1601,  1604,  1608,  1609,  1610,  1617,  1618,  1619,
    1620,  1621,  1622,  1624,  1626,  1630,  1631,  1635,  1636,  1637,
    1638,  1643,  1644,  1650,  1657,  1665,  1666,  1670,  1671,  1676,
    1679,  1687,  1684,  1702,  1705,  1708,  1709,  1713,  1718,  1719,
    1723,  1726,  1728,  1734,  1735,  1739,  1739,  1751,  1752,  1755,
    1756,  1757,  1758,  1759,  1760,  1761,  1765,  1766,  1771,  1772,
    1773,  1774,  1778,  1783,  1787,  1791,  1792,  1795,  1796,  1797,
    1800,  1803,  1804,  1807,  1810,  1811,  1815,  1821,  1822,  1827,
    1828,  1827,  1838,  1835,  1848,  1852,  1856,  1860,  1870,  1873,
    1867,  1880,  1881,  1877,  1890,  1891,  1895,  1896,  1900,  1901,
    1905,  1906,  1909,  1910,  1924,  1930,  1931,  1936,  1937,  1939,
    1936,  1948,  1951,  1953,  1958,  1959,  1963,  1970,  1976,  1977,
    1982,  1982,  1992,  1991,  2002,  2003,  2014,  2015,  2016,  2019,
    2023,  2031,  2038,  2039,  2040,  2054,  2055,  2056,  2060,  2060,
    2066,  2067,  2068,  2072,  2076,  2080,  2081,  2090,  2094,  2095,
    2096,  2097,  2098,  2099,  2100,  2101,  2102,  2107,  2107,  2109,
    2110,  2110,  2114,  2115,  2116,  2117,  2118,  2121,  2124,  2128,
    2138,  2139,  2140,  2141,  2142,  2143,  2147,  2152,  2157,  2162,
    2167,  2172,  2176,  2177,  2178,  2181,  2182,  2185,  2186,  2187,
    2188,  2189,  2190,  2191,  2192,  2195,  2196,  2203,  2203,  2210,
    2211,  2215,  2216,  2219,  2220,  2221,  2225,  2226,  2236,  2239,
    2243,  2249,  2250,  2265,  2266,  2267,  2271,  2277,  2278,  2282,
    2283,  2287,  2289,  2293,  2294,  2298,  2299,  2302,  2308,  2309,
    2326,  2331,  2332,  2336,  2342,  2343,  2360,  2364,  2365,  2366,
    2373,  2389,  2393,  2394,  2404,  2407,  2425,  2426,  2435,  2439,
    2443,  2464,  2465,  2466,  2467
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
  "EXCLUDE_REGION", "EXIT", "EXP", "EXPRESSION", "EXTERN", "FALSE",
  "FCLOSE", "FILENAME", "FILENAME_PREFIX", "FILE_OUTPUT_REPORT",
  "FINAL_SUMMARY", "FLOOR", "FOPEN", "FORMAT", "FPRINTF", "FPRINT_TIME",
  "FRONT", "FRONT_CROSSINGS", "FRONT_HITS", "GAUSSIAN_RELEASE_NUMBER",
  "GEOMETRY", "GRAPH_PATTERN", "HEADER", "HIGH_PROBABILITY_THRESHOLD",
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
  "MOLECULE_PLACEMENT_FAILURE", "NAME_LIST", "NEAREST_POINT",
  "NEAREST_TRIANGLE", "NEGATIVE_DIFFUSION_CONSTANT",
  "NEGATIVE_REACTION_RATE", "NO", "NOEXIT", "NONE", "NO_SPECIES",
  "NOT_EQUAL", "NOTIFICATIONS", "NUMBER_OF_SUBUNITS", "NUMBER_OF_TRAINS",
  "NUMBER_TO_RELEASE", "OBJECT", "OFF", "ON", "ORIENTATIONS",
  "OUTPUT_BUFFER_SIZE", "INVALID_OUTPUT_STEP_TIME",
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
  "target_def", "maximum_step_length_def", "extern_def",
  "existing_molecule", "existing_surface_molecule",
  "existing_molecule_opt_orient", "surface_classes_def",
  "define_one_surface_class", "define_multiple_surface_classes",
  "list_surface_class_stmts", "surface_class_stmt", "$@1",
  "existing_surface_class", "list_surface_prop_stmts", "surface_prop_stmt",
  "surface_rxn_stmt", "surface_rxn_type", "equals_or_to",
  "surface_class_mol_stmt", "surface_mol_stmt", "list_surface_mol_density",
  "list_surface_mol_num", "surface_mol_quant", "rx_net_def",
  "list_rx_stmts", "rx_stmt", "list_dashes", "right_arrow", "left_arrow",
  "double_arrow", "right_cat_arrow", "double_cat_arrow", "reaction_arrow",
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
     505,   506,   507,   508,   509,   510,   511,   512,   513,    61,
      38,    58,    43,    45,    42,    47,   514,    94,    91,    93,
      59,    39,    44,   123,   125,    40,    41,    62,    60,    64
};
# endif

#define YYPACT_NINF -1031

#define yypact_value_is_default(Yystate) \
  (!!((Yystate) == (-1031)))

#define YYTABLE_NINF -401

#define yytable_value_is_error(Yytable_value) \
  0

  /* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
     STATE-NUM.  */
static const yytype_int16 yypact[] =
{
    2559,  -175,  -162,   -93,   -41,   -17,    -6,  -114,  -138,    82,
    -114,  -114,   115,   135,    15,   178,   200,   177,   224,   228,
     246, -1031,   251,   254,   262,   274,   281,   288,   291,   300,
     282,   287, -1031, -1031, -1031,   292,   244,   279,   308,   310,
     301,   323,   309,   326,   327,   329, -1031,   316,   319,   322,
     596,  2559, -1031,   -42, -1031, -1031,   339, -1031, -1031,   517,
   -1031, -1031, -1031, -1031, -1031, -1031, -1031, -1031, -1031, -1031,
   -1031, -1031,   342, -1031, -1031, -1031, -1031, -1031, -1031, -1031,
   -1031, -1031, -1031, -1031,   165, -1031, -1031, -1031, -1031,   431,
   -1031, -1031, -1031, -1031, -1031, -1031, -1031, -1031, -1031,   144,
     144,   -29,  2430,   -29,  2430, -1031, -1031,   341,  -114,  -114,
   -1031,   355, -1031,   358, -1031,  -114,  -114,   -29,   -51,  2430,
    -114,  -114,  -114,   -29,  -114,  2430,  2430,   144,  2430,  2430,
    2430,  2430,   397,  -114,   980, -1031,   576,   -29,   -29,  2029,
    2430,   458,  2430,  -114,  2430,  2430,  2430, -1031,   555,   685,
   -1031, -1031,  1873,   361,   104,   292, -1031, -1031,   292, -1031,
     292, -1031, -1031,   292,   292,   292, -1031, -1031, -1031, -1031,
   -1031, -1031, -1031, -1031,   362, -1031, -1031, -1031, -1031, -1031,
     383, -1031, -1031,   369,   370,   374,   378,   379,   385,   387,
     389, -1031,   390,   391,   398,   399,   401, -1031, -1031, -1031,
   -1031,   402, -1031,   410,   412,   413,   415,  2430,  2430,  2430,
   -1031,   -19, -1031, -1031, -1031, -1031, -1031,  1006,     4,   142,
    -179, -1031, -1031,   350, -1031,  -172, -1031, -1031,  -116, -1031,
   -1031, -1031,   -71, -1031, -1031, -1031,    22, -1031,   383,   409,
   -1031, -1031,   528, -1031,   396,   420,   423,   383, -1031,   543,
   -1031,   528,   528, -1031,   528,   528,   528,   528, -1031, -1031,
   -1031,   432,   426,    46, -1031,   443,   444,   448,   450,   451,
     455,   456,   459,   462,   465,   468,   469,   470,   473,   475,
     477,   478,   481,   351, -1031,   482,   383, -1031,   441, -1031,
     528,   528,   485, -1031,   528, -1031,   479,   528,   528,   528,
     624,   504,   627,   508,   521,   522,   523,   539,   542,   546,
     548,   550,   552,   559,   560,   567,   568,   572,   584,   590,
     599,   881, -1031,  2144,   638, -1031, -1031,   528,   764, -1031,
     841,   409,   -29, -1031, -1031, -1031, -1031,   719,  -114, -1031,
     532, -1031,   532,   -29,   -29,  2430,  2430,  2430,  2430,  2430,
    2430,  2430,  2430,  2430,  2430,  2430,  2430,  2430,  2430,  2430,
    2430,   -29,  2430,   589,   589,   394, -1031, -1031,  2430,  2430,
    2430,  2430,  2430, -1031,  2430, -1031,   604,   608,   289, -1031,
   -1031, -1031, -1031, -1031,  2430, -1031,   108, -1031, -1031, -1031,
   -1031, -1031,  -114,  -114,  -137,     1, -1031, -1031, -1031,   549,
   -1031, -1031, -1031,   -29,   -29,  -114, -1031, -1031, -1031,   144,
     144,   144,    39,   144,   144,  1680,   144,   144,   144,  2430,
     144,    39,   144,   144,   144,    39,    39, -1031, -1031,   104,
     -25, -1031,  2430,   -13,   -29,   609,    17, -1031,   -29,   611,
     -32, -1031,   -31,   -31,   -31,  2430,   -31,  2430,   -31,   -31,
    2430,   -31,   -31,   -31,   -31,   -31,   -31,   -31,   -31,   -31,
   -1031, -1031,  2430,   102, -1031,   528,   603,   612,   713, -1031,
     296,  -114, -1031, -1031,   687,   616,   666,   538,   833, -1031,
   -1031,   417,   588,   614,   620,   665,   684,   805,   831,   848,
     872,  1102,  1108,  1114,  1135,   977,   984,  -183,   999, -1031,
     260,   260,   589,   589, -1031,  1153,  2430,  2430,   648,   649,
     696,   811, -1031, -1031, -1031, -1031,   350, -1031, -1031,   668,
    -159, -1031,  -164, -1031, -1031, -1031,   -81,   677,   678,   679,
     680,   686, -1031,    37,  -114, -1031,   657,   681, -1031, -1031,
   -1031, -1031, -1031, -1031, -1031, -1031, -1031, -1031,   528, -1031,
   -1031, -1031, -1031,   528, -1031, -1031, -1031, -1031, -1031, -1031,
   -1031,   671, -1031,  1911, -1031,   528,   697,   698,   699,   -24,
   -1031, -1031, -1031, -1031,    28,   703,   682,   -33, -1031, -1031,
   -1031, -1031,   383,  -114,   704, -1031,   710, -1031, -1031, -1031,
   -1031, -1031, -1031,   528, -1031,   528, -1031, -1031,   528, -1031,
   -1031, -1031, -1031, -1031, -1031, -1031, -1031, -1031,   558, -1031,
    2144,   -29,   104,    18,   152, -1031,   708,   538,   104,   700,
   -1031,   717,   718,   715,   722,   724,   730,   729,   731,   732,
     736, -1031, -1031,   745,   738,   538, -1031,   746, -1031, -1031,
   -1031, -1031, -1031,   740, -1031,   204, -1031, -1031, -1031, -1031,
   -1031, -1031, -1031, -1031, -1031, -1031,  2430,  2430,  2430,  2430,
   -1031, -1031, -1031, -1031,  2430,   528,   528,  2430,  2430, -1031,
     854, -1031, -1031,   751, -1031, -1031,   668, -1031,   668, -1031,
   -1031,  -181, -1031,  2430,  2182,  2430,  2430,  2430, -1031,  -114,
     743,   748, -1031, -1031, -1031, -1031, -1031,    86, -1031, -1031,
   -1031,   749,   215, -1031, -1031,   -47,   104, -1031, -1031,   409,
   -1031,   104,  2430,   104,   763,   771, -1031,    67, -1031, -1031,
   -1031, -1031,   220, -1031, -1031, -1031,   -29,   -21, -1031, -1031,
   -1031, -1031,   761,   104,   775,   785,  2430, -1031,   383,   768,
     767,   292,   786,   787,   791, -1031, -1031, -1031, -1031,    -3,
     538, -1031, -1031,   -63,   104, -1031,  2430,  2430,   931,   -29,
     104,  -114,  -114,  2430,  -114,  2430,   104,   933,   152, -1031,
    2300,   104, -1031, -1031,  1052,  1067,  1086,  1092,  1167,   528,
     528,   799,   987,   -23, -1031, -1031,   -81,   295,   801, -1031,
   -1031,   528, -1031,   528, -1031,   528,   528,   528,   807,  -114,
    -114, -1031, -1031,    16, -1031, -1031,   818, -1031, -1031, -1031,
   -1031,   893, -1031,   528, -1031,   306,   144,    -4, -1031, -1031,
   -1031,   383,   806,   809,   810,   -64, -1031, -1031, -1031, -1031,
    -114, -1031,  2300,   825,    55,   483, -1031,   104, -1031,   104,
    2300,   104, -1031, -1031, -1031, -1031, -1031, -1031,  -168,   432,
   -1031,   202,   152, -1031, -1031, -1031, -1031,   -46,   152,   528,
     528,   827,   383, -1031, -1031,   104,    68, -1031,   528, -1031,
   -1031,   528, -1031,   830, -1031,  1173, -1031, -1031, -1031, -1031,
      50, -1031,   -14, -1031, -1031, -1031, -1031,  2430,  2430, -1031,
     816, -1031,  1911,  1911, -1031, -1031,   409,   163, -1031,  -114,
   -1031,  2430,   350,   832,    80, -1031,   154, -1031,   350, -1031,
     824,  -114,   850, -1031, -1031, -1031,   383, -1031, -1031, -1031,
     851,   863, -1031,    -4,    -4, -1031,    21, -1031,   889, -1031,
      26, -1031,    26, -1031, -1031, -1031,  1173, -1031, -1031, -1031,
    2300,   857,   883,   901,   826,  2430,  1080, -1031,   868, -1031,
   -1031,   193,  -168,  -168,  -168, -1031, -1031, -1031, -1031,  2430,
   -1031, -1031, -1031,  2430, -1031, -1031,   888,   895,   201, -1031,
   -1031, -1031,   528,   528, -1031, -1031, -1031, -1031,   295, -1031,
     528, -1031,  2430, -1031, -1031, -1031, -1031, -1031,   447, -1031,
     144,   973,   894,  2430,    -4,   907, -1031,   184,    -4,    91,
     -29,    -4,    -4,    -4,    -4, -1031, -1031, -1031, -1031,    12,
   -1031,   890,    20,    13, -1031,   898, -1031,   104,  2430,   104,
   -1031,   734,   919, -1031,   152,  2430, -1031,   915,   915, -1031,
     343,   510,  -114, -1031, -1031,   911,   528,   922, -1031, -1031,
     925, -1031, -1031,   447, -1031, -1031, -1031, -1031,   942, -1031,
     945, -1031,   956,  1048,    95,  1030,   597,    95, -1031, -1031,
     941,   943,   944,   -29,   383,   151,   151, -1031, -1031, -1031,
   -1031,    10,   963, -1031, -1031, -1031, -1031,   963, -1031, -1031,
       8, -1031,   528, -1031, -1031,  2430, -1031, -1031,   528,   966,
   -1031,   968,   171, -1031,   959,  1210, -1031,   965,   967, -1031,
   -1031,  -114,   104,   144,   978,  1066,   971,   972,   982,   985,
     983, -1031, -1031, -1031, -1031, -1031, -1031, -1031, -1031,   989,
   -1031, -1031,   979, -1031, -1031, -1031, -1031, -1031,  2430, -1031,
   -1031, -1031, -1031, -1031,   528,   -14,  2430,  2430, -1031, -1031,
   -1031, -1031, -1031, -1031, -1031, -1031, -1031, -1031,   571,   991,
   -1031,   447, -1031,  1002, -1031,  1462,  1462,   145, -1031,  1004,
   -1031,   144,  1018, -1031,   157, -1031,   157,   157, -1031, -1031,
   -1031,   528, -1031,   882,    41,   447,  2430, -1031,  1462,   210,
     225, -1031,   104, -1031,   144,  1007, -1031,   432, -1031,  1011,
    1012,  1014,   152, -1031,  1029,   447,   528, -1031, -1031, -1031,
   -1031, -1031, -1031,    89, -1031,    89, -1031,    89, -1031, -1031,
    2430, -1031, -1031, -1031, -1031, -1031, -1031, -1031, -1031, -1031,
    1017, -1031,  1017,  1017,   909,   212,   580, -1031, -1031, -1031,
   -1031, -1031
};

  /* YYDEFACT[STATE-NUM] -- Default reduction number in state STATE-NUM.
     Performed when YYTABLE does not specify something else to do.  Zero
     means the default is an error.  */
static const yytype_uint16 yydefact[] =
{
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,   323,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,   219,   220,   221,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,    27,     0,     0,     0,
       0,     2,     3,   331,     5,     6,     0,     7,   112,     0,
     113,   114,   115,   116,   117,   118,     8,     9,    10,    11,
      14,    12,     0,    15,   222,   223,    16,   245,   246,    17,
      18,    20,    19,   325,     0,   326,   327,   348,   347,     0,
     329,   330,    13,   328,    21,    22,    23,    24,    25,     0,
       0,     0,     0,     0,     0,   229,   224,     0,     0,     0,
     313,     0,   230,     0,   247,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,   332,     0,     0,     0,     0,
       0,   494,     0,     0,     0,     0,     0,   567,     0,     0,
       1,     4,     0,     0,     0,     0,   367,   368,     0,   369,
       0,   366,   370,     0,     0,     0,    35,    37,    39,    38,
      34,    36,   201,   200,     0,   108,    26,   107,   111,   184,
      28,   105,   106,     0,     0,     0,     0,     0,     0,     0,
       0,    70,     0,     0,     0,     0,     0,    93,    95,    94,
      71,     0,    96,     0,     0,     0,     0,     0,     0,     0,
      74,   189,    66,    68,    69,    67,   185,   192,   189,     0,
       0,   226,   242,    40,   294,     0,   275,   277,   295,   292,
     315,   251,     0,   249,    29,   477,     0,   475,     0,   211,
     212,   213,   206,   123,     0,     0,     0,    55,   331,     0,
     324,   207,   199,   187,   214,   215,   216,   217,   209,   210,
     208,     0,     0,     0,   488,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,   136,     0,   124,   125,     0,   204,
     203,   205,     0,   492,   197,    60,     0,   196,   198,   202,
     571,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,   161,     0,    61,    58,    59,     0,    72,    56,
      73,    57,     0,    65,   218,    62,    63,     0,     0,   349,
       0,   364,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,   104,   103,     0,   191,   190,     0,     0,
       0,     0,     0,   186,     0,   188,     0,     0,   233,   225,
     227,    43,    48,    49,     0,   244,    41,    44,    45,    42,
     274,   276,     0,     0,     0,     0,   254,   248,   250,     0,
     474,   476,   122,     0,     0,     0,   490,   487,   489,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,   135,   137,     0,
       0,   133,     0,     0,     0,     0,     0,   572,     0,     0,
       0,   612,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     160,   162,     0,     0,    51,    53,     0,     0,   331,   344,
       0,   334,   341,   343,     0,     0,     0,     0,     0,   125,
     109,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,    75,
      98,    99,   100,   101,   102,   193,     0,     0,     0,     0,
     236,     0,    46,    47,   293,   253,    40,   296,   278,     0,
       0,   285,     0,   287,   286,   288,     0,     0,     0,     0,
       0,     0,   312,     0,     0,   125,     0,     0,   482,   157,
     138,   145,   153,   159,   158,   140,   147,   148,   155,   154,
     156,   144,   141,   143,   139,   150,   146,   149,   142,   152,
     151,     0,    31,     0,   130,   495,     0,     0,     0,   502,
     496,   497,   498,   125,     0,     0,     0,     0,   569,   577,
     576,   578,   611,     0,     0,   613,     0,   183,   181,   182,
     163,   168,   169,   167,   166,   172,   171,   173,   174,   175,
     177,   164,   165,   178,   179,   180,   170,   176,     0,    64,
       0,     0,     0,     0,     0,   342,     0,     0,     0,     0,
     452,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,   385,   386,     0,     0,   334,   371,     0,   376,   387,
     388,   389,   390,     0,   401,     0,    91,    88,    87,    89,
      83,    85,    76,    82,    77,    78,     0,     0,     0,     0,
      84,    90,    97,    86,     0,   232,   231,     0,     0,   237,
     238,    50,   297,   281,   279,   280,     0,   282,     0,   300,
     301,     0,   298,     0,     0,     0,     0,     0,   263,     0,
       0,     0,   261,   262,   252,   255,   256,     0,   257,   266,
     481,     0,     0,   134,    30,     0,     0,   129,   127,   128,
     126,     0,     0,     0,     0,     0,   508,     0,   503,   505,
     506,   507,     0,   574,   575,   573,     0,     0,   568,   570,
     615,   616,   614,     0,     0,     0,     0,    52,   121,     0,
       0,     0,     0,     0,     0,   333,   340,   335,   336,     0,
     334,   404,   405,     0,     0,   334,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,   372,
       0,     0,   411,   110,     0,     0,     0,     0,   194,   235,
     234,     0,   240,     0,   283,   284,     0,     0,   289,   302,
     303,   316,   322,   321,   320,   317,   319,   318,     0,     0,
       0,   265,   264,     0,   478,   131,     0,   491,   485,   483,
     484,   470,   500,   499,   501,     0,     0,     0,   493,   504,
     132,   579,     0,     0,     0,     0,   581,   583,   584,   585,
       0,   618,     0,     0,   621,     0,   119,     0,   345,     0,
       0,     0,   354,   355,   358,   356,   353,   357,     0,   352,
     359,   351,     0,   403,   406,   455,   456,     0,     0,   395,
     396,     0,   384,   374,   375,     0,     0,   397,   391,   314,
     382,   381,   380,     0,   365,   373,   378,   377,   379,   410,
       0,   408,   334,    79,    80,    92,    81,     0,     0,   241,
       0,   299,     0,     0,   311,   309,   310,     0,   306,     0,
     291,     0,    40,     0,     0,   269,     0,   271,    40,   258,
       0,     0,     0,   458,   510,   511,   512,   513,   514,   527,
       0,     0,   530,     0,     0,   518,     0,   515,   565,   519,
       0,   590,     0,   580,   582,   617,    65,    32,   619,    33,
       0,     0,     0,     0,     0,     0,   472,   334,     0,   338,
     337,     0,     0,     0,     0,   350,   454,   457,   453,     0,
     399,   383,   398,     0,   407,   409,     0,     0,     0,   412,
     413,   414,   195,   239,   228,   307,   308,   304,     0,   290,
     260,   243,     0,   267,   270,   268,   272,   259,     0,   486,
       0,   464,     0,     0,     0,     0,   525,     0,     0,     0,
       0,     0,     0,     0,     0,   517,   607,   609,   608,     0,
     604,     0,     0,     0,   598,     0,   620,     0,     0,     0,
     610,     0,     0,   461,     0,     0,   360,   361,   362,   363,
       0,     0,     0,   415,   402,     0,   273,     0,   445,   442,
       0,   444,   441,   479,   426,   428,   429,   430,     0,   431,
       0,   471,     0,   466,     0,     0,     0,     0,   520,   516,
       0,     0,   532,     0,   566,   521,   522,   523,   524,   603,
     605,     0,   588,   586,   594,   593,   589,   588,   597,   599,
       0,   623,   622,   624,    54,     0,   411,   346,   339,     0,
     392,     0,     0,   447,     0,     0,   305,     0,     0,   427,
     482,     0,     0,     0,     0,   468,     0,   538,     0,     0,
       0,   540,   541,   542,   543,   544,   545,   529,   526,     0,
     533,   536,   534,   537,   509,   601,   602,   606,     0,   592,
     591,   595,   596,   600,   473,   462,     0,     0,   446,   448,
     449,   425,   422,   420,   421,   423,   424,   419,   437,     0,
     439,   417,   418,   434,   435,     0,     0,     0,   440,     0,
     465,     0,     0,   459,     0,   539,     0,     0,   528,   531,
     535,   587,   334,     0,     0,     0,     0,   416,     0,     0,
       0,   480,     0,   467,     0,     0,   552,   554,   553,   555,
     555,   555,     0,   393,     0,   450,   438,   436,   433,   432,
     443,   469,   460,     0,   548,     0,   546,     0,   547,   463,
       0,   482,   562,   564,   559,   561,   558,   563,   560,   557,
     555,   556,   555,   555,     0,     0,     0,   551,   549,   550,
     394,   451
};

  /* YYPGOTO[NTERM-NUM].  */
static const yytype_int16 yypgoto[] =
{
   -1031, -1031, -1031,  1247,  -845,     0,  -102,  -111,  -131,  -397,
    -782,   -96,  -501, -1031,   914,   917,   194, -1031,   694, -1031,
   -1031,  1162,  -139,  -145,  -140, -1031,   840,  -423,  -136,  -146,
   -1031,  -124,   -76,  -109, -1031, -1031, -1031, -1031, -1031, -1031,
     352,  -120,  -393, -1031, -1031, -1031, -1031, -1031, -1031, -1031,
   -1031,  1024,   499,    63, -1031, -1031,   990,   674, -1031,  1095,
   -1031, -1031, -1031, -1031, -1031, -1031, -1031, -1031, -1031,   -58,
   -1031, -1031, -1031, -1031, -1031, -1031, -1031, -1031, -1031,  -496,
   -1031, -1031, -1031, -1031,   -67, -1031,   407, -1031, -1031, -1031,
   -1031, -1031, -1031,   777, -1031, -1031,  -649, -1031, -1031,  1096,
    -322,  -235, -1031, -1031, -1031, -1031, -1031, -1031, -1031, -1031,
     930, -1031, -1031, -1031,   537, -1031, -1031, -1031,   346,  -384,
   -1031, -1031, -1031, -1031, -1031, -1031, -1031, -1031,  -282,  -107,
     -16,  -722,  -621, -1031, -1031,  1209, -1031,   864, -1031, -1031,
   -1031, -1031, -1031, -1031,  -569, -1031, -1031, -1031,   720, -1031,
    -582, -1031, -1031, -1031, -1031, -1031, -1031, -1031,   472, -1031,
   -1031, -1031,   997,   591, -1031, -1031, -1031,   461,   256, -1031,
   -1031, -1031, -1031, -1031, -1030,  -994, -1031, -1031, -1031,  -624,
     167, -1031, -1031, -1031, -1031, -1031, -1031,   266, -1031, -1031,
   -1031, -1031, -1031,   495, -1031, -1031, -1031, -1031, -1031, -1031,
   -1031, -1031, -1031, -1031, -1031, -1031, -1031, -1031,  1124, -1031,
   -1031, -1031,   829, -1020, -1031, -1031, -1031, -1031,  1098, -1031,
   -1031, -1031, -1031, -1031, -1031, -1031, -1031, -1031,   667, -1031,
   -1031, -1031, -1031, -1031, -1031,   392,  -362, -1031, -1031, -1031,
   -1031, -1031, -1031, -1031,   328, -1031, -1031, -1031, -1031, -1031,
   -1031,  -628,  -800, -1031, -1031, -1031, -1031, -1031, -1031, -1031,
     812, -1031, -1031, -1031, -1031,   562, -1031,   311, -1031, -1031,
   -1031, -1031, -1031, -1031,   381, -1031, -1031, -1031,   382,  -873,
   -1031, -1031, -1031,   955,   553, -1031, -1031, -1031, -1031, -1031
};

  /* YYDEFGOTO[NTERM-NUM].  */
static const yytype_int16 yydefgoto[] =
{
      -1,    50,    51,    52,   177,   210,   179,   261,   850,   937,
     938,   539,   385,   386,   387,   388,   389,   463,   464,    54,
      55,    56,   894,   562,   335,   336,   327,   212,   213,   895,
     214,   215,   238,   181,   182,    57,    58,    59,   739,    60,
     244,   287,   430,   710,    61,    62,    63,    64,    65,    66,
     283,   284,   540,   545,    67,   321,   322,   590,    68,   373,
     218,    69,    70,    71,    72,    73,    74,    75,   220,   106,
     107,   113,   378,   510,   670,   782,   890,   223,   903,   224,
      76,    77,    78,   232,   114,   396,   516,   533,   695,   696,
     697,   803,   698,   808,   904,   906,   905,    79,   225,   226,
     783,   521,   522,   523,   524,   525,   526,   900,   227,   228,
     229,   394,   517,   681,   682,   788,   789,   790,   897,   898,
      80,   111,   870,   395,   794,    81,   124,    82,    83,    84,
     338,   746,   614,   747,   748,    85,   471,   472,   473,   947,
      86,    87,   474,   617,   851,    88,   477,   164,   635,   878,
     636,   637,   638,   639,   640,   641,   642,   866,   867,    89,
      90,   772,   476,   752,   753,   644,   880,   881,   882,   969,
     970,  1095,  1149,  1150,  1043,  1044,  1045,  1046,  1152,  1153,
    1154,  1047,  1048,  1049,  1050,   971,  1092,  1093,  1175,  1211,
      91,   755,   620,   856,   857,    92,   991,  1185,    93,  1086,
    1172,  1053,  1105,  1163,   913,  1023,    94,   236,   237,   399,
     910,  1100,  1094,   705,   809,   810,    95,   263,   264,   538,
      96,   433,   293,   569,   570,   571,   572,   717,   718,   719,
     817,   917,   720,   721,   926,   927,   928,   929,   992,   995,
    1063,  1124,  1108,  1109,  1110,  1111,  1112,  1113,  1114,  1115,
    1116,  1189,  1204,  1221,  1005,    97,   300,   577,   436,   437,
     578,   579,   580,   581,   825,   826,   827,  1129,  1012,  1076,
    1077,  1133,   828,  1013,  1014,  1127,   829,  1009,  1010,  1011,
      98,   302,   440,   441,   731,   732,   586,   735,   834,   944
};

  /* YYTABLE[YYPACT[STATE-NUM]] -- What to do in state STATE-NUM.  If
     positive, shift that token.  If negative, reduce the rule whose
     number is the opposite.  If YYTABLE_NINF, syntax error.  */
static const yytype_int16 yytable[] =
{
      53,   216,   262,   172,   173,   235,   329,   105,   239,   334,
     110,   112,   326,   325,   768,   672,   328,   249,   288,   136,
    1006,  1006,  1125,  1131,   676,   180,   678,   180,   330,   908,
     680,   253,   561,  1072,  1006,   587,   260,   919,   966,   842,
    -120,   688,   822,   331,   575,   723,   874,   247,   233,  1099,
     221,    53,   366,   769,   174,   527,   469,   724,   949,  1015,
     843,   286,   286,    46,   588,  1151,   543,   679,   714,   920,
      46,   175,   520,   689,    46,   366,   823,   344,    46,  1194,
    1157,   786,   584,    46,    99,   822,   645,   787,   690,   691,
     240,   241,   566,   662,   575,   379,  1212,   100,   191,  1213,
    1214,   178,   390,   178,   674,   576,   439,   848,   105,   222,
     844,  1215,  1216,   677,   166,   112,   234,   178,   675,   823,
     243,   243,   243,   178,   248,   235,   518,   367,    46,   852,
     955,  -400,   262,   234,   858,   108,   958,   178,   178,   337,
    1015,   519,   702,   295,   339,  1195,   392,   340,   341,   342,
     367,   907,   324,   528,   333,   576,  1217,  1099,  1016,   714,
     941,    46,   380,   393,   824,   398,   101,  1074,   769,  1007,
    1007,    46,   690,   691,   806,   725,  1218,  1219,   715,    46,
     722,   967,   200,  1007,   167,   845,  1008,  1008,   176,   469,
     155,  1225,   529,   798,   168,   169,    46,   673,   376,   377,
    1008,  1099,   567,   397,   846,   323,   847,   824,   156,  1106,
     933,   853,  1106,    46,   568,   740,   921,   -60,   102,   166,
     105,   751,   323,   479,   589,   222,   692,   807,   956,   157,
     466,   470,   112,   922,   530,   531,   234,   176,    46,    46,
     674,   728,   103,   368,   369,   370,   371,   563,   372,   716,
      46,   564,    46,   104,   675,   984,   180,   986,    46,   923,
      46,   968,    46,   234,    46,   864,   865,   286,   480,   715,
     942,   924,   848,   693,   117,   532,   323,   170,   323,   951,
     690,   691,   943,   535,   536,   497,  1069,  1078,    46,   167,
     680,   135,    46,   998,  1073,   999,   400,   171,   183,   168,
     169,   184,  1087,   368,   369,   370,   371,   909,   372,   811,
      46,   694,   176,   185,   573,   186,   544,   801,   323,   158,
     407,   155,    46,   187,   964,   544,  1024,   286,   286,   544,
     544,    46,   178,   508,   509,   188,   831,    46,   468,   156,
     716,   818,   961,   178,   178,   802,    46,   690,   691,   742,
    1060,   743,   159,  1061,   983,   109,   751,   323,   286,   160,
     157,   178,   582,   863,   470,   265,   806,   189,  1062,   872,
     865,   609,   323,   161,   610,   162,   190,   266,   174,   382,
     383,   166,   170,  1027,  1028,  1029,   267,   744,   115,   174,
    1206,  1208,   222,   515,   925,   175,    46,   191,   742,    46,
     743,   981,   171,   178,   178,   537,   175,   987,   116,   268,
     192,   193,   194,    46,  1186,  1003,  1004,   708,   163,  1181,
    1227,   195,  1228,  1229,   707,   196,   745,   328,   985,   333,
     269,   270,   977,   806,   178,   978,   744,   118,   178,   330,
     946,   784,   948,   785,   950,  1138,  1001,  1002,  1003,  1004,
     613,   167,   120,   914,   709,   952,   953,   954,   271,   119,
    1058,   168,   169,  1202,   952,   953,   954,   197,   960,  1026,
    1209,   468,   166,   245,   246,  1034,   563,   198,   199,  1198,
     773,   200,  1178,   159,   555,   272,  1231,   563,   559,   560,
     160,   805,   563,   201,  1199,   202,   820,  1178,   203,   121,
     925,   925,   741,   122,   161,   123,   162,   204,   975,   976,
     125,   205,   176,   126,   273,  1037,  1038,  1039,   206,   137,
     222,   127,   222,   176,   370,   371,   222,   372,   274,   275,
     276,  1179,  1180,   128,   700,   738,   277,    46,  1190,  1191,
     129,   278,   167,  1040,   170,  1041,  1042,   130,    46,   163,
     131,  1192,   168,   169,   138,   133,  1089,   207,   208,   132,
     134,   996,   997,   324,   171,   135,   812,   139,   814,   140,
     209,   925,   892,   893,   141,   925,   621,   279,   925,   925,
     925,   925,   142,   730,   143,   144,   145,  1212,   146,   147,
    1213,  1214,   148,   622,   280,   149,   150,   281,   152,   153,
     282,   154,  1215,  1216,   165,   368,   369,   370,   371,   855,
     372,   178,   333,   234,   219,   285,   292,  1090,   333,   258,
     381,   382,   383,   384,   877,   427,   879,   623,   230,   624,
     876,   231,  1056,   301,   328,   170,   332,   343,   849,  1065,
    1066,  1067,  1068,   344,   345,   346,   330,  1217,   259,   347,
     821,   328,   625,   348,   349,   171,   368,   369,   370,   371,
     350,   372,   351,   330,   352,   353,   354,  1218,  1219,  -105,
     499,   626,   402,   355,   356,   627,   357,   358,   896,   368,
     369,   370,   371,   862,   372,   359,   939,   360,   361,   222,
     362,   628,   403,   646,   939,   404,   328,   158,   945,   406,
     405,   303,   409,   410,   328,  1159,   333,   411,   330,   412,
     413,   333,   855,   333,   414,   415,   330,   431,   416,   915,
     918,   417,   304,  1091,   418,   838,   178,   419,   420,   421,
     629,   630,   422,   333,   423,   879,   424,   425,   305,   916,
     426,   429,   631,   632,   432,   368,   369,   370,   371,   234,
     372,   434,   633,   333,   333,   435,   328,   328,   467,   178,
     333,   222,   222,   438,   869,   439,   333,   442,   330,   330,
     875,   333,   368,   369,   370,   371,   475,   372,   306,   307,
     443,   444,   445,   896,   896,  1200,   222,   324,   634,   736,
     368,   369,   370,   371,   939,   372,   308,   309,   446,   902,
     902,   447,  1176,   222,   328,   448,  1220,   449,  1222,   450,
    1223,   451,   310,   311,   312,   178,   330,   534,   452,   453,
     368,   369,   370,   371,   313,   372,   454,   455,   314,   315,
     730,   456,   936,   368,   369,   370,   371,   333,   372,   333,
     936,   333,   328,   457,   316,   317,   318,   319,   234,   458,
     368,   369,   370,   371,   330,   372,   372,   333,   459,  1001,
    1002,  1003,  1004,   506,   647,   333,   222,   507,   574,   896,
     583,   612,  1081,  1118,  1083,   611,   368,   369,   370,   371,
     333,   372,   368,   369,   370,   371,  -400,   372,   616,   618,
     648,   619,   324,   324,  1051,   643,   649,   303,  -111,   979,
     -74,   -74,   -74,   -74,   902,   -74,   902,   667,   668,   541,
     542,   515,   546,   547,   549,   550,   551,   552,   304,   554,
     669,   556,   557,   558,  1064,   320,  1126,   368,   369,   370,
     371,   518,   372,   703,   305,  1132,   683,   684,   685,   686,
     936,   650,   211,   706,   217,   687,   368,   369,   370,   371,
     704,   372,   234,   234,   234,   727,   711,   712,   713,   242,
     651,  1123,   726,   733,   734,   251,   252,   749,   254,   255,
     256,   257,   781,   754,   306,   307,   756,   757,   324,   290,
     291,   759,   294,   760,   297,   298,   299,   180,   758,   761,
     763,   764,   308,   309,   265,   765,   368,   369,   370,   371,
     178,   372,   762,  1084,   766,   770,   266,  1160,   310,   311,
     312,   767,  1075,   771,   674,   267,   799,   333,   804,   333,
     313,   800,   815,   830,   314,   315,   -68,   -68,   -68,   -68,
     816,   -68,   700,  1188,   832,  1188,  1188,   833,   268,   837,
     316,   317,   318,   319,   836,   839,   840,   363,   364,   365,
     841,   861,   873,  1187,  1107,  1187,  1187,  1107,   888,   269,
     270,   889,   899,   178,   912,  1183,   901,   368,   369,   370,
     371,   333,   372,   368,   369,   370,   371,   911,   372,   930,
     333,   652,   931,   932,   940,   671,   959,   271,  1201,   963,
     974,   982,   700,   368,   369,   370,   371,   988,   372,  1022,
    1020,  1158,   333,   -67,   -67,   -67,   -67,   653,   -67,   990,
     368,   369,   370,   371,   272,   372,  1017,   591,   592,   993,
     594,   320,   596,   597,   654,   599,   600,   601,   602,   603,
     604,   605,   606,   607,   368,   369,   370,   371,   994,   372,
    1025,  1052,  1018,   273,   368,   369,   370,   371,   655,   372,
    1000,  1001,  1002,  1003,  1004,   460,  1193,   274,   275,   276,
    1019,  1032,  1054,   465,   234,   277,   234,   234,  1033,  1071,
     278,   368,   369,   370,   371,  1057,   372,  1080,  1085,   954,
    1096,  1097,   333,  1230,  1098,   481,   482,   483,   484,   485,
     486,   487,   488,   489,   490,   491,   492,   493,   494,   495,
     496,  1101,   498,   333,  1102,   333,   279,   333,   500,   501,
     502,   503,   504,   183,   505,  1103,   184,  1104,  1120,  1141,
    1121,  1122,  1128,   280,   511,  1136,   281,  1137,   185,   282,
     186,  1142,  1140,  1155,  1143,  1156,  1162,  1161,   187,   368,
     369,   370,   371,  1164,   372,   384,   368,   369,   370,   371,
     188,   372,  1168,   660,  1166,   548,  1170,  1167,  1169,   553,
     661,   368,   369,   370,   371,  1177,   372,   374,   368,   369,
     370,   371,   565,   372,  1178,   663,  1182,  1184,  1037,  1038,
    1039,   745,   189,  1203,  1205,   593,  1207,   595,  1210,  1226,
     598,   190,   368,   369,   370,   371,  1144,   372,   151,  1117,
     512,  1165,   608,   513,   737,   296,  1040,   428,  1041,  1042,
     699,   461,   191,   375,   368,   369,   370,   371,   989,   372,
    1145,   391,   514,   891,  1035,   192,   193,   194,   883,   368,
     369,   370,   371,   250,   372,   615,   195,   750,   962,   478,
     196,   965,  1135,   884,   854,  1197,   665,   666,   368,   369,
     370,   371,   957,   372,   368,   369,   370,   371,  1139,   372,
     401,   408,   885,   701,   368,   369,   370,   371,   886,   372,
     368,   369,   370,   371,   656,   372,   368,   369,   370,   371,
     657,   372,   197,   935,   819,  1119,   658,   934,  1130,   729,
    1059,  1070,   198,   199,  1079,   585,   200,   368,   369,   370,
     371,     0,   372,     0,     0,     0,  1146,   659,   201,     0,
     202,     0,     0,   203,   664,   368,   369,   370,   371,     0,
     372,     0,   204,     0,     0,     0,   205,     0,   887,   368,
     369,   370,   371,   206,   372,   -74,   -74,   -74,   -74,     0,
     -74,     0,  1147,     0,     0,     0,     0,     0,     0,     0,
     465,     0,    46,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,   183,     0,     0,   184,     0,
       0,  1141,   207,   208,     0,     0,     0,     0,     0,     0,
     185,     0,   186,  1142,     0,   209,  1143,     0,     0,     0,
     187,     0,     0,     0,     0,     0,   774,   775,   776,   777,
       0,     0,   188,     0,   778,     0,     0,   779,   780,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,   791,   793,   795,   796,   797,     0,     0,
       0,     0,     0,     0,   189,     0,     0,     0,     0,     0,
       0,     0,     0,   190,     0,     0,     0,     0,  1144,     0,
       0,     0,   813,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,   191,     0,     0,     0,     0,     0,
       0,     0,  1145,     0,     0,     0,   835,   192,   193,   194,
       0,     0,     0,     0,     0,     0,     0,     0,   195,     0,
       0,     0,   196,     0,     0,     0,   859,   860,     0,     0,
       0,     0,     0,   868,     0,   871,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,   197,     0,     0,     0,     0,     0,
       0,     0,     0,     0,   198,   199,     0,     0,   200,     0,
       0,     0,     0,     0,     0,     0,     0,     0,  1146,     0,
     201,     0,   202,     0,     0,   203,     0,     0,     0,     0,
       0,     0,     0,     0,   204,     0,     0,     0,   205,     0,
       0,     0,     0,   183,     0,   206,   184,     0,     0,     0,
       0,     0,     0,     0,  1147,     0,     0,     0,   185,     0,
     186,     0,     0,     0,    46,     0,     0,     0,   187,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     188,     0,     0,     0,   207,   208,     0,   972,   973,     0,
       0,     0,     0,     0,     0,     0,     0,   209,     0,     0,
       0,   980,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,   189,     0,     0,   166,     0,     0,     0,     0,
       0,   190,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,   191,     0,     0,  1021,     0,     0,     0,     0,
       0,     0,     0,     0,     0,   192,   193,   194,     0,  1030,
       0,     0,     0,  1031,     0,     0,   195,     0,     0,     0,
     196,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,  1036,     0,     0,   167,     0,     0,     0,     0,
       0,     0,     0,  1055,     0,   168,   169,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,   197,     0,     0,     0,     0,     0,  1082,     0,
       0,     0,   198,   199,     0,  1088,   200,     0,     0,     0,
       0,     0,     0,     0,     0,     0,   183,     0,   201,   184,
     202,     0,     0,   203,     0,     0,     0,     0,     0,     0,
       0,   185,   204,   186,     0,     0,   205,     0,     0,     0,
       0,   187,     0,   206,     0,     0,     0,     0,     0,     0,
       0,     0,     0,   188,   183,     0,     0,   184,   170,     0,
       0,     0,    46,     0,     0,  1134,     0,     0,     0,   185,
       0,   186,     0,     0,     0,  1148,     0,     0,   171,   187,
       0,     0,   207,   208,     0,   189,     0,     0,     0,     0,
       0,   188,     0,     0,   190,   209,   174,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,  1171,     0,
       0,     0,     0,   175,     0,   191,  1173,  1174,     0,     0,
       0,     0,     0,   189,     0,     0,     0,     0,   192,   193,
     194,     0,   190,     0,   174,  1148,  1148,     0,     0,   195,
       0,     0,     0,   196,     0,     0,     0,     0,     0,     0,
       0,   175,     0,   191,     0,     0,  1196,     0,  1148,     0,
       0,     0,     0,     0,     0,     0,   192,   193,   194,     0,
       0,     0,   183,     0,     0,   184,     0,   195,     0,     0,
       0,   196,     0,     0,     0,   197,     0,   185,     0,   186,
    1224,     0,     0,     0,     0,   198,   199,   187,     0,   200,
       0,     0,     0,     0,     0,     0,     0,     0,     0,   188,
       0,   201,     0,   202,     0,     0,   203,     0,     0,     0,
       0,     0,     0,   197,     0,   204,     0,     0,     0,   205,
     176,     0,     0,   198,   199,     0,   206,   200,     0,     0,
       0,   189,     0,     0,     0,     0,     0,     0,     0,   201,
     190,   202,     0,     0,   203,    46,     0,     0,     0,     0,
       0,     0,     0,   204,     0,     0,     0,   205,   176,     0,
       0,   191,   289,     0,   206,   207,   208,     0,     0,     0,
       0,   323,     0,     0,   192,   193,   194,   183,   209,     0,
     184,     0,     0,    46,     0,   195,     0,     0,     0,   196,
       0,     0,   185,     0,   186,     0,     0,     0,     0,     0,
       0,     0,   187,   207,   208,     0,     0,     0,     0,     0,
       0,     0,     0,     0,   188,   183,   209,     0,   184,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     185,   197,   186,     0,     0,     0,     0,     0,     0,     0,
     187,   198,   199,     0,     0,   200,   189,     0,     0,     0,
       0,     0,   188,     0,     0,   190,     0,   201,     0,   202,
       0,     0,   203,     0,     0,     0,     0,     0,     0,     0,
       0,   204,     0,     0,     0,   205,   191,     0,     0,     0,
       0,     0,   206,     0,   189,     0,     0,     0,     0,   192,
     193,   194,     0,   190,     0,     0,     0,     0,     0,     0,
     195,    46,     0,     0,   196,     0,     0,     0,     0,     0,
       0,     0,     0,     0,   191,     0,     0,     0,     0,     0,
       0,   207,   208,     0,     0,     0,     0,   192,   193,   194,
       0,     0,     0,   183,   209,     0,   184,     0,   195,     0,
       0,     0,   196,     0,     0,     0,   197,     0,   185,     0,
     186,     0,     0,     0,     0,     0,   198,   199,   187,     0,
     200,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     188,     0,   201,     0,   202,     0,     0,   203,     0,     0,
       0,     0,     0,     0,   197,     0,   204,     0,     0,     0,
     205,     0,     0,     0,   198,   199,     0,   206,   200,     0,
       0,     0,   189,     0,     0,     0,     0,     0,     0,     0,
     201,   190,   202,     0,     0,   203,    46,     0,     0,     0,
       0,     0,     0,     0,   204,     0,     0,     0,   205,     0,
       0,     0,   191,     0,     0,   206,   207,   208,     0,     0,
       0,     0,   462,     0,     0,   192,   193,   194,     0,   209,
       0,   792,     0,     0,    46,     0,   195,     0,     0,     0,
     196,     0,     0,   183,     0,     0,   184,     0,     0,     0,
       0,     0,     0,     0,   207,   208,     0,     0,   185,     0,
     186,     0,     0,     0,     0,     0,     0,   209,   187,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     188,     0,   197,     0,     0,     0,     0,     0,     0,     0,
       0,     0,   198,   199,     0,     0,   200,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,   201,     0,
     202,     0,   189,   203,     0,     0,     0,     0,     0,     0,
       0,   190,   204,     0,     0,     0,   205,     0,     0,     0,
       0,     0,     0,   206,     0,     0,     0,     0,     0,     0,
       0,     0,   191,     0,     0,     0,     0,     0,     0,     0,
       0,     0,    46,     0,     0,   192,   193,   194,     0,     0,
       0,     0,     0,     0,     0,     0,   195,     0,     0,     0,
     196,     0,   207,   208,     1,     0,     0,     0,   323,     0,
       0,     0,     0,     0,     0,   209,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     2,
       3,     4,     5,     6,     0,     0,     0,     0,     0,     0,
       0,     0,   197,     0,     0,     7,     8,     9,    10,    11,
      12,    13,   198,   199,     0,     0,   200,     0,    14,    15,
      16,     0,     0,     0,     0,     0,     0,     0,   201,     0,
     202,     0,     0,   203,     0,    17,     0,     0,     0,     0,
       0,     0,   204,    18,    19,     0,   205,     0,     0,     0,
       0,     0,     0,   206,     0,     0,    20,     0,     0,     0,
      21,     0,     0,    22,     0,     0,     0,    23,    24,     0,
       0,     0,    46,     0,     0,     0,     0,     0,     0,     0,
      25,    26,    27,    28,    29,     0,     0,     0,     0,     0,
       0,    30,   207,   208,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,   209,     0,     0,     0,    31,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,    32,    33,    34,    35,     0,     0,     0,
       0,     0,     0,     0,    36,    37,     0,     0,     0,    38,
      39,     0,     0,    40,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,    41,     0,     0,     0,     0,
      42,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,    43,    44,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
      45,    46,     0,     0,    47,     0,     0,    48,     0,     0,
       0,     0,     0,     0,     0,    49
};

static const yytype_int16 yycheck[] =
{
       0,   103,   133,    99,   100,   116,   152,     7,   117,   154,
      10,    11,   152,   152,   635,   516,   152,   124,   138,    35,
       8,     8,    12,    15,   520,   101,   522,   103,   152,    13,
     526,   127,   429,    13,     8,    66,   132,    41,    52,    42,
      82,     4,   106,   152,    77,    17,   768,   123,   115,  1043,
     108,    51,    71,   635,    83,    54,   338,    29,   840,   932,
      63,   137,   138,   242,    95,  1095,    27,   148,    92,    73,
     242,   100,   394,    36,   242,    71,   140,   260,   242,    38,
    1100,   262,   114,   242,   259,   106,   479,   268,   135,   136,
     141,   142,   105,   276,    77,   274,     7,   259,   102,    10,
      11,   101,   274,   103,   263,   138,   138,   275,   108,   109,
     113,    22,    23,   277,    75,   115,   116,   117,   277,   140,
     120,   121,   122,   123,   124,   236,   263,   146,   242,   750,
     852,   173,   263,   133,   755,   273,   858,   137,   138,   155,
    1013,   278,   535,   143,   160,  1175,   262,   163,   164,   165,
     146,   800,   152,   152,   154,   138,    67,  1151,   940,    92,
     105,   242,   220,   279,   228,   232,   259,  1012,   750,   157,
     157,   242,   135,   136,   221,   147,    87,    88,   202,   242,
     573,   195,   186,   157,   145,   188,   174,   174,   217,   471,
      25,  1211,   191,   689,   155,   156,   242,   519,    56,    57,
     174,  1195,   215,   274,   207,   268,   209,   228,    43,  1054,
     274,   274,  1057,   242,   227,   612,   220,   259,   259,    75,
     220,   618,   268,   343,   255,   225,   189,   274,   274,    64,
     332,   338,   232,   237,   233,   234,   236,   217,   242,   242,
     263,   274,   259,   262,   263,   264,   265,   272,   267,   273,
     242,   276,   242,   259,   277,   904,   332,   906,   242,   263,
     242,   882,   242,   263,   242,   761,   762,   343,   344,   202,
     215,   275,   275,   236,   259,   274,   268,   238,   268,   848,
     135,   136,   227,   403,   404,   361,   274,   274,   242,   145,
     786,   273,   242,   272,   274,   274,   274,   258,     3,   155,
     156,     6,  1024,   262,   263,   264,   265,   803,   267,   706,
     242,   274,   217,    18,   434,    20,   412,   231,   268,   154,
     274,    25,   242,    28,   274,   421,   947,   403,   404,   425,
     426,   242,   332,    44,    45,    40,   733,   242,   338,    43,
     273,   274,   274,   343,   344,   259,   242,   135,   136,   197,
     259,   199,   187,   262,   274,   273,   753,   268,   434,   194,
      64,   361,   438,   760,   471,    14,   221,    72,   277,   766,
     866,   269,   268,   208,   272,   210,    81,    26,    83,   271,
     272,    75,   238,   952,   953,   954,    35,   235,   273,    83,
    1190,  1191,   392,   393,   817,   100,   242,   102,   197,   242,
     199,   902,   258,   403,   404,   405,   100,   908,   273,    58,
     115,   116,   117,   242,   257,   264,   265,   563,   253,   274,
    1220,   126,  1222,  1223,   563,   130,   274,   563,   274,   429,
      79,    80,   269,   221,   434,   272,   235,   259,   438,   563,
     837,   676,   839,   678,   841,   274,   262,   263,   264,   265,
     154,   145,   275,   147,   563,   262,   263,   264,   107,   259,
     276,   155,   156,  1185,   262,   263,   264,   172,   865,   276,
    1192,   471,    75,   121,   122,   274,   272,   182,   183,   269,
     276,   186,   272,   187,   421,   134,   274,   272,   425,   426,
     194,   276,   272,   198,   269,   200,   276,   272,   203,   275,
     923,   924,   613,   275,   208,   259,   210,   212,   892,   893,
     259,   216,   217,   259,   163,    68,    69,    70,   223,   275,
     520,   259,   522,   217,   264,   265,   526,   267,   177,   178,
     179,  1155,  1156,   259,   534,   611,   185,   242,  1166,  1167,
     259,   190,   145,    96,   238,    98,    99,   259,   242,   253,
     259,  1172,   155,   156,   275,   273,   213,   262,   263,   259,
     273,   923,   924,   563,   258,   273,   711,   259,   713,   259,
     275,   994,   277,   278,   273,   998,    38,   226,  1001,  1002,
    1003,  1004,   259,   583,   275,   259,   259,     7,   259,   273,
      10,    11,   273,    55,   243,   273,     0,   246,   259,    82,
     249,   259,    22,    23,   173,   262,   263,   264,   265,   754,
     267,   611,   612,   613,   273,    39,   158,   274,   618,   222,
     270,   271,   272,   273,   770,   274,   771,    89,   273,    91,
     770,   273,   994,    78,   770,   238,   275,   275,   749,  1001,
    1002,  1003,  1004,   260,   275,   275,   770,    67,   251,   275,
     726,   787,   114,   275,   275,   258,   262,   263,   264,   265,
     275,   267,   275,   787,   275,   275,   275,    87,    88,   260,
     276,   133,   276,   275,   275,   137,   275,   275,   787,   262,
     263,   264,   265,   759,   267,   275,   832,   275,   275,   689,
     275,   153,   272,   276,   840,   272,   832,   154,   215,   273,
     268,    16,   259,   259,   840,  1102,   706,   259,   832,   259,
     259,   711,   857,   713,   259,   259,   840,   276,   259,   815,
     816,   259,    37,   213,   259,   741,   726,   259,   259,   259,
     192,   193,   259,   733,   259,   880,   259,   259,    53,   815,
     259,   259,   204,   205,   259,   262,   263,   264,   265,   749,
     267,   272,   214,   753,   754,   131,   892,   893,    39,   759,
     760,   761,   762,   259,   764,   138,   766,   259,   892,   893,
     770,   771,   262,   263,   264,   265,   244,   267,    93,    94,
     259,   259,   259,   892,   893,  1182,   786,   787,   250,   231,
     262,   263,   264,   265,   940,   267,   111,   112,   259,   799,
     800,   259,   231,   803,   940,   259,  1203,   259,  1205,   259,
    1207,   259,   127,   128,   129,   815,   940,   268,   259,   259,
     262,   263,   264,   265,   139,   267,   259,   259,   143,   144,
     830,   259,   832,   262,   263,   264,   265,   837,   267,   839,
     840,   841,   978,   259,   159,   160,   161,   162,   848,   259,
     262,   263,   264,   265,   978,   267,   267,   857,   259,   262,
     263,   264,   265,   259,   276,   865,   866,   259,   259,   978,
     259,   259,  1017,   276,  1019,   272,   262,   263,   264,   265,
     880,   267,   262,   263,   264,   265,   173,   267,   201,   273,
     276,   225,   892,   893,   990,    62,   276,    16,   260,   899,
     262,   263,   264,   265,   904,   267,   906,   259,   259,   410,
     411,   911,   413,   414,   415,   416,   417,   418,    37,   420,
     224,   422,   423,   424,  1000,   240,  1071,   262,   263,   264,
     265,   263,   267,   276,    53,  1080,   259,   259,   259,   259,
     940,   276,   102,   272,   104,   259,   262,   263,   264,   265,
     269,   267,   952,   953,   954,   273,   259,   259,   259,   119,
     276,  1063,   259,   259,   254,   125,   126,   259,   128,   129,
     130,   131,   118,   273,    93,    94,   259,   259,   978,   139,
     140,   259,   142,   259,   144,   145,   146,  1063,   273,   259,
     259,   259,   111,   112,    14,   259,   262,   263,   264,   265,
    1000,   267,   273,   269,   259,   259,    26,  1103,   127,   128,
     129,   273,  1012,   273,   263,    35,   273,  1017,   269,  1019,
     139,   273,   259,   262,   143,   144,   262,   263,   264,   265,
     259,   267,  1032,  1164,   259,  1166,  1167,   252,    58,   272,
     159,   160,   161,   162,   276,   259,   259,   207,   208,   209,
     259,   120,   119,  1164,  1054,  1166,  1167,  1057,   259,    79,
      80,    74,   261,  1063,   171,  1161,   259,   262,   263,   264,
     265,  1071,   267,   262,   263,   264,   265,   259,   267,   273,
    1080,   276,   273,   273,   259,   274,   259,   107,  1184,   259,
     274,   259,  1092,   262,   263,   264,   265,   273,   267,    19,
     274,  1101,  1102,   262,   263,   264,   265,   276,   267,   259,
     262,   263,   264,   265,   134,   267,   259,   443,   444,   268,
     446,   240,   448,   449,   276,   451,   452,   453,   454,   455,
     456,   457,   458,   459,   262,   263,   264,   265,   275,   267,
     272,   168,   259,   163,   262,   263,   264,   265,   276,   267,
     261,   262,   263,   264,   265,   274,   274,   177,   178,   179,
     259,   273,   268,   323,  1164,   185,  1166,  1167,   273,   279,
     190,   262,   263,   264,   265,   268,   267,   279,   259,   264,
     269,   259,  1182,   274,   259,   345,   346,   347,   348,   349,
     350,   351,   352,   353,   354,   355,   356,   357,   358,   359,
     360,   259,   362,  1203,   259,  1205,   226,  1207,   368,   369,
     370,   371,   372,     3,   374,   259,     6,   169,   277,     9,
     277,   277,   259,   243,   384,   259,   246,   259,    18,   249,
      20,    21,   273,   268,    24,   268,   170,   259,    28,   262,
     263,   264,   265,   272,   267,   273,   262,   263,   264,   265,
      40,   267,   269,   276,   272,   415,   277,   272,   269,   419,
     276,   262,   263,   264,   265,   274,   267,   261,   262,   263,
     264,   265,   432,   267,   272,   276,   272,   259,    68,    69,
      70,   274,    72,   272,   272,   445,   272,   447,   259,   272,
     450,    81,   262,   263,   264,   265,    86,   267,    51,   269,
     386,  1107,   462,   386,   610,   143,    96,   283,    98,    99,
     533,   321,   102,   218,   262,   263,   264,   265,   911,   267,
     110,   225,   392,   786,   978,   115,   116,   117,   276,   262,
     263,   264,   265,   124,   267,   471,   126,   617,   866,   342,
     130,   880,  1086,   276,   753,  1178,   506,   507,   262,   263,
     264,   265,   857,   267,   262,   263,   264,   265,  1092,   267,
     236,   263,   276,   534,   262,   263,   264,   265,   276,   267,
     262,   263,   264,   265,   272,   267,   262,   263,   264,   265,
     272,   267,   172,   830,   717,  1057,   272,   825,  1077,   577,
     998,  1009,   182,   183,  1013,   440,   186,   262,   263,   264,
     265,    -1,   267,    -1,    -1,    -1,   196,   272,   198,    -1,
     200,    -1,    -1,   203,   261,   262,   263,   264,   265,    -1,
     267,    -1,   212,    -1,    -1,    -1,   216,    -1,   261,   262,
     263,   264,   265,   223,   267,   262,   263,   264,   265,    -1,
     267,    -1,   232,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
     610,    -1,   242,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,     3,    -1,    -1,     6,    -1,
      -1,     9,   262,   263,    -1,    -1,    -1,    -1,    -1,    -1,
      18,    -1,    20,    21,    -1,   275,    24,    -1,    -1,    -1,
      28,    -1,    -1,    -1,    -1,    -1,   656,   657,   658,   659,
      -1,    -1,    40,    -1,   664,    -1,    -1,   667,   668,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,   683,   684,   685,   686,   687,    -1,    -1,
      -1,    -1,    -1,    -1,    72,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    81,    -1,    -1,    -1,    -1,    86,    -1,
      -1,    -1,   712,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,   102,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,   110,    -1,    -1,    -1,   736,   115,   116,   117,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   126,    -1,
      -1,    -1,   130,    -1,    -1,    -1,   756,   757,    -1,    -1,
      -1,    -1,    -1,   763,    -1,   765,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,   172,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,   182,   183,    -1,    -1,   186,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   196,    -1,
     198,    -1,   200,    -1,    -1,   203,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,   212,    -1,    -1,    -1,   216,    -1,
      -1,    -1,    -1,     3,    -1,   223,     6,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,   232,    -1,    -1,    -1,    18,    -1,
      20,    -1,    -1,    -1,   242,    -1,    -1,    -1,    28,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      40,    -1,    -1,    -1,   262,   263,    -1,   887,   888,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,   275,    -1,    -1,
      -1,   901,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    72,    -1,    -1,    75,    -1,    -1,    -1,    -1,
      -1,    81,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,   102,    -1,    -1,   945,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,   115,   116,   117,    -1,   959,
      -1,    -1,    -1,   963,    -1,    -1,   126,    -1,    -1,    -1,
     130,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,   982,    -1,    -1,   145,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,   993,    -1,   155,   156,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,   172,    -1,    -1,    -1,    -1,    -1,  1018,    -1,
      -1,    -1,   182,   183,    -1,  1025,   186,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,     3,    -1,   198,     6,
     200,    -1,    -1,   203,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    18,   212,    20,    -1,    -1,   216,    -1,    -1,    -1,
      -1,    28,    -1,   223,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    40,     3,    -1,    -1,     6,   238,    -1,
      -1,    -1,   242,    -1,    -1,  1085,    -1,    -1,    -1,    18,
      -1,    20,    -1,    -1,    -1,  1095,    -1,    -1,   258,    28,
      -1,    -1,   262,   263,    -1,    72,    -1,    -1,    -1,    -1,
      -1,    40,    -1,    -1,    81,   275,    83,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,  1128,    -1,
      -1,    -1,    -1,   100,    -1,   102,  1136,  1137,    -1,    -1,
      -1,    -1,    -1,    72,    -1,    -1,    -1,    -1,   115,   116,
     117,    -1,    81,    -1,    83,  1155,  1156,    -1,    -1,   126,
      -1,    -1,    -1,   130,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,   100,    -1,   102,    -1,    -1,  1176,    -1,  1178,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,   115,   116,   117,    -1,
      -1,    -1,     3,    -1,    -1,     6,    -1,   126,    -1,    -1,
      -1,   130,    -1,    -1,    -1,   172,    -1,    18,    -1,    20,
    1210,    -1,    -1,    -1,    -1,   182,   183,    28,    -1,   186,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    40,
      -1,   198,    -1,   200,    -1,    -1,   203,    -1,    -1,    -1,
      -1,    -1,    -1,   172,    -1,   212,    -1,    -1,    -1,   216,
     217,    -1,    -1,   182,   183,    -1,   223,   186,    -1,    -1,
      -1,    72,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   198,
      81,   200,    -1,    -1,   203,   242,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,   212,    -1,    -1,    -1,   216,   217,    -1,
      -1,   102,   103,    -1,   223,   262,   263,    -1,    -1,    -1,
      -1,   268,    -1,    -1,   115,   116,   117,     3,   275,    -1,
       6,    -1,    -1,   242,    -1,   126,    -1,    -1,    -1,   130,
      -1,    -1,    18,    -1,    20,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    28,   262,   263,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    40,     3,   275,    -1,     6,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      18,   172,    20,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      28,   182,   183,    -1,    -1,   186,    72,    -1,    -1,    -1,
      -1,    -1,    40,    -1,    -1,    81,    -1,   198,    -1,   200,
      -1,    -1,   203,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,   212,    -1,    -1,    -1,   216,   102,    -1,    -1,    -1,
      -1,    -1,   223,    -1,    72,    -1,    -1,    -1,    -1,   115,
     116,   117,    -1,    81,    -1,    -1,    -1,    -1,    -1,    -1,
     126,   242,    -1,    -1,   130,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,   102,    -1,    -1,    -1,    -1,    -1,
      -1,   262,   263,    -1,    -1,    -1,    -1,   115,   116,   117,
      -1,    -1,    -1,     3,   275,    -1,     6,    -1,   126,    -1,
      -1,    -1,   130,    -1,    -1,    -1,   172,    -1,    18,    -1,
      20,    -1,    -1,    -1,    -1,    -1,   182,   183,    28,    -1,
     186,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      40,    -1,   198,    -1,   200,    -1,    -1,   203,    -1,    -1,
      -1,    -1,    -1,    -1,   172,    -1,   212,    -1,    -1,    -1,
     216,    -1,    -1,    -1,   182,   183,    -1,   223,   186,    -1,
      -1,    -1,    72,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
     198,    81,   200,    -1,    -1,   203,   242,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,   212,    -1,    -1,    -1,   216,    -1,
      -1,    -1,   102,    -1,    -1,   223,   262,   263,    -1,    -1,
      -1,    -1,   268,    -1,    -1,   115,   116,   117,    -1,   275,
      -1,   239,    -1,    -1,   242,    -1,   126,    -1,    -1,    -1,
     130,    -1,    -1,     3,    -1,    -1,     6,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,   262,   263,    -1,    -1,    18,    -1,
      20,    -1,    -1,    -1,    -1,    -1,    -1,   275,    28,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      40,    -1,   172,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,   182,   183,    -1,    -1,   186,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   198,    -1,
     200,    -1,    72,   203,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    81,   212,    -1,    -1,    -1,   216,    -1,    -1,    -1,
      -1,    -1,    -1,   223,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,   102,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,   242,    -1,    -1,   115,   116,   117,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,   126,    -1,    -1,    -1,
     130,    -1,   262,   263,     5,    -1,    -1,    -1,   268,    -1,
      -1,    -1,    -1,    -1,    -1,   275,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    30,
      31,    32,    33,    34,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,   172,    -1,    -1,    46,    47,    48,    49,    50,
      51,    52,   182,   183,    -1,    -1,   186,    -1,    59,    60,
      61,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   198,    -1,
     200,    -1,    -1,   203,    -1,    76,    -1,    -1,    -1,    -1,
      -1,    -1,   212,    84,    85,    -1,   216,    -1,    -1,    -1,
      -1,    -1,    -1,   223,    -1,    -1,    97,    -1,    -1,    -1,
     101,    -1,    -1,   104,    -1,    -1,    -1,   108,   109,    -1,
      -1,    -1,   242,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
     121,   122,   123,   124,   125,    -1,    -1,    -1,    -1,    -1,
      -1,   132,   262,   263,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,   275,    -1,    -1,    -1,   150,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,   164,   165,   166,   167,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,   175,   176,    -1,    -1,    -1,   180,
     181,    -1,    -1,   184,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,   206,    -1,    -1,    -1,    -1,
     211,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   229,   230,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
     241,   242,    -1,    -1,   245,    -1,    -1,   248,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,   256
};

  /* YYSTOS[STATE-NUM] -- The (internal number of the) accessing
     symbol of state STATE-NUM.  */
static const yytype_uint16 yystos[] =
{
       0,     5,    30,    31,    32,    33,    34,    46,    47,    48,
      49,    50,    51,    52,    59,    60,    61,    76,    84,    85,
      97,   101,   104,   108,   109,   121,   122,   123,   124,   125,
     132,   150,   164,   165,   166,   167,   175,   176,   180,   181,
     184,   206,   211,   229,   230,   241,   242,   245,   248,   256,
     281,   282,   283,   285,   299,   300,   301,   315,   316,   317,
     319,   324,   325,   326,   327,   328,   329,   334,   338,   341,
     342,   343,   344,   345,   346,   347,   360,   361,   362,   377,
     400,   405,   407,   408,   409,   415,   420,   421,   425,   439,
     440,   470,   475,   478,   486,   496,   500,   535,   560,   259,
     259,   259,   259,   259,   259,   285,   349,   350,   273,   273,
     285,   401,   285,   351,   364,   273,   273,   259,   259,   259,
     275,   275,   275,   259,   406,   259,   259,   259,   259,   259,
     259,   259,   259,   273,   273,   273,   410,   275,   275,   259,
     259,   273,   259,   275,   259,   259,   259,   273,   273,   273,
       0,   283,   259,    82,   259,    25,    43,    64,   154,   187,
     194,   208,   210,   253,   427,   173,    75,   145,   155,   156,
     238,   258,   291,   291,    83,   100,   217,   284,   285,   286,
     312,   313,   314,     3,     6,    18,    20,    28,    40,    72,
      81,   102,   115,   116,   117,   126,   130,   172,   182,   183,
     186,   198,   200,   203,   212,   216,   223,   262,   263,   275,
     285,   306,   307,   308,   310,   311,   286,   306,   340,   273,
     348,   349,   285,   357,   359,   378,   379,   388,   389,   390,
     273,   273,   363,   364,   285,   287,   487,   488,   312,   313,
     141,   142,   306,   285,   320,   320,   320,   312,   285,   409,
     415,   306,   306,   291,   306,   306,   306,   306,   222,   251,
     291,   287,   288,   497,   498,    14,    26,    35,    58,    79,
      80,   107,   134,   163,   177,   178,   179,   185,   190,   226,
     243,   246,   249,   330,   331,    39,   312,   321,   321,   103,
     306,   306,   158,   502,   306,   285,   301,   306,   306,   306,
     536,    78,   561,    16,    37,    53,    93,    94,   111,   112,
     127,   128,   129,   139,   143,   144,   159,   160,   161,   162,
     240,   335,   336,   268,   285,   302,   304,   306,   308,   309,
     311,   313,   275,   285,   303,   304,   305,   410,   410,   410,
     410,   410,   410,   275,   260,   275,   275,   275,   275,   275,
     275,   275,   275,   275,   275,   275,   275,   275,   275,   275,
     275,   275,   275,   306,   306,   306,    71,   146,   262,   263,
     264,   265,   267,   339,   261,   339,    56,    57,   352,   274,
     349,   270,   271,   272,   273,   292,   293,   294,   295,   296,
     274,   379,   262,   279,   391,   403,   365,   274,   364,   489,
     274,   488,   276,   272,   272,   268,   273,   274,   498,   259,
     259,   259,   259,   259,   259,   259,   259,   259,   259,   259,
     259,   259,   259,   259,   259,   259,   259,   274,   331,   259,
     322,   276,   259,   501,   272,   131,   538,   539,   259,   138,
     562,   563,   259,   259,   259,   259,   259,   259,   259,   259,
     259,   259,   259,   259,   259,   259,   259,   259,   259,   259,
     274,   336,   268,   297,   298,   306,   286,    39,   285,   408,
     409,   416,   417,   418,   422,   244,   442,   426,   442,   321,
     312,   306,   306,   306,   306,   306,   306,   306,   306,   306,
     306,   306,   306,   306,   306,   306,   306,   312,   306,   276,
     306,   306,   306,   306,   306,   306,   259,   259,    44,    45,
     353,   306,   294,   295,   390,   285,   366,   392,   263,   278,
     380,   381,   382,   383,   384,   385,   386,    54,   152,   191,
     233,   234,   274,   367,   268,   321,   321,   285,   499,   291,
     332,   332,   332,    27,   291,   333,   332,   332,   306,   332,
     332,   332,   332,   306,   332,   333,   332,   332,   332,   333,
     333,   289,   303,   272,   276,   306,   105,   215,   227,   503,
     504,   505,   506,   321,   259,    77,   138,   537,   540,   541,
     542,   543,   312,   259,   114,   563,   566,    66,    95,   255,
     337,   337,   337,   306,   337,   306,   337,   337,   306,   337,
     337,   337,   337,   337,   337,   337,   337,   337,   306,   269,
     272,   272,   259,   154,   412,   417,   201,   423,   273,   225,
     472,    38,    55,    89,    91,   114,   133,   137,   153,   192,
     193,   204,   205,   214,   250,   428,   430,   431,   432,   433,
     434,   435,   436,    62,   445,   322,   276,   276,   276,   276,
     276,   276,   276,   276,   276,   276,   272,   272,   272,   272,
     276,   276,   276,   276,   261,   306,   306,   259,   259,   224,
     354,   274,   292,   380,   263,   277,   359,   277,   359,   148,
     359,   393,   394,   259,   259,   259,   259,   259,     4,    36,
     135,   136,   189,   236,   274,   368,   369,   370,   372,   373,
     285,   492,   322,   276,   269,   493,   272,   302,   309,   313,
     323,   259,   259,   259,    92,   202,   273,   507,   508,   509,
     512,   513,   322,    17,    29,   147,   259,   273,   274,   540,
     285,   564,   565,   259,   254,   567,   231,   298,   312,   318,
     289,   287,   197,   199,   235,   274,   411,   413,   414,   259,
     428,   289,   443,   444,   273,   471,   259,   259,   273,   259,
     259,   259,   273,   259,   259,   259,   259,   273,   412,   430,
     259,   273,   441,   276,   306,   306,   306,   306,   306,   306,
     306,   118,   355,   380,   381,   381,   262,   268,   395,   396,
     397,   306,   239,   306,   404,   306,   306,   306,   359,   273,
     273,   231,   259,   371,   269,   276,   221,   274,   373,   494,
     495,   289,   303,   306,   303,   259,   259,   510,   274,   508,
     276,   312,   106,   140,   228,   544,   545,   546,   552,   556,
     262,   289,   259,   252,   568,   306,   276,   272,   410,   259,
     259,   259,    42,    63,   113,   188,   207,   209,   275,   287,
     288,   424,   412,   274,   443,   303,   473,   474,   412,   306,
     306,   120,   312,   289,   359,   359,   437,   438,   306,   285,
     402,   306,   289,   119,   411,   285,   304,   309,   429,   303,
     446,   447,   448,   276,   276,   276,   276,   261,   259,    74,
     356,   394,   277,   278,   302,   309,   313,   398,   399,   261,
     387,   259,   285,   358,   374,   376,   375,   376,    13,   359,
     490,   259,   171,   484,   147,   291,   312,   511,   291,    41,
      73,   220,   237,   263,   275,   307,   514,   515,   516,   517,
     273,   273,   273,   274,   545,   564,   285,   289,   290,   309,
     259,   105,   215,   227,   569,   215,   289,   419,   289,   290,
     289,   424,   262,   263,   264,   411,   274,   473,   411,   259,
     289,   274,   438,   259,   274,   447,    52,   195,   412,   449,
     450,   465,   306,   306,   274,   399,   399,   269,   272,   285,
     306,   292,   259,   274,   376,   274,   376,   292,   273,   366,
     259,   476,   518,   268,   275,   519,   516,   516,   272,   274,
     261,   262,   263,   264,   265,   534,     8,   157,   174,   557,
     558,   559,   548,   553,   554,   559,   290,   259,   259,   259,
     274,   306,    19,   485,   412,   272,   276,   424,   424,   424,
     306,   306,   273,   273,   274,   398,   306,    68,    69,    70,
      96,    98,    99,   454,   455,   456,   457,   461,   462,   463,
     464,   291,   168,   481,   268,   306,   516,   268,   276,   515,
     259,   262,   277,   520,   312,   516,   516,   516,   516,   274,
     558,   279,    13,   274,   284,   285,   549,   550,   274,   554,
     279,   303,   306,   303,   269,   259,   479,   411,   306,   213,
     274,   213,   466,   467,   492,   451,   269,   259,   259,   455,
     491,   259,   259,   259,   169,   482,   284,   285,   522,   523,
     524,   525,   526,   527,   528,   529,   530,   269,   276,   524,
     277,   277,   277,   286,   521,    12,   303,   555,   259,   547,
     547,    15,   303,   551,   306,   448,   259,   259,   274,   467,
     273,     9,    21,    24,    86,   110,   196,   232,   306,   452,
     453,   454,   458,   459,   460,   268,   268,   493,   285,   289,
     291,   259,   170,   483,   272,   296,   272,   272,   269,   269,
     277,   306,   480,   306,   306,   468,   231,   274,   272,   459,
     459,   274,   272,   291,   259,   477,   257,   287,   288,   531,
     531,   531,   412,   274,    38,   454,   306,   460,   269,   269,
     289,   291,   411,   272,   532,   272,   532,   272,   532,   411,
     259,   469,     7,    10,    11,    22,    23,    67,    87,    88,
     289,   533,   289,   289,   306,   493,   272,   532,   532,   532,
     274,   274
};

  /* YYR1[YYN] -- Symbol number of symbol that rule YYN derives.  */
static const yytype_uint16 yyr1[] =
{
       0,   280,   281,   282,   282,   283,   283,   283,   283,   283,
     283,   283,   283,   283,   283,   283,   283,   283,   283,   283,
     283,   283,   283,   283,   283,   283,   284,   285,   286,   287,
     288,   289,   290,   290,   291,   291,   291,   291,   291,   291,
     292,   292,   292,   292,   293,   293,   293,   293,   294,   295,
     296,   297,   297,   298,   298,   299,   300,   300,   300,   300,
     301,   302,   303,   303,   304,   305,   306,   306,   307,   307,
     308,   308,   309,   309,   310,   311,   311,   311,   311,   311,
     311,   311,   311,   311,   311,   311,   311,   311,   311,   311,
     311,   311,   311,   311,   311,   311,   311,   311,   311,   311,
     311,   311,   311,   311,   311,   312,   312,   313,   313,   313,
     313,   314,   315,   315,   315,   315,   315,   315,   315,   316,
     317,   318,   319,   320,   321,   322,   322,   323,   323,   323,
     324,   325,   326,   327,   328,   329,   330,   330,   331,   331,
     331,   331,   331,   331,   331,   331,   331,   331,   331,   331,
     331,   331,   331,   331,   331,   331,   331,   332,   333,   333,
     334,   335,   335,   336,   336,   336,   336,   336,   336,   336,
     336,   336,   336,   336,   336,   336,   336,   336,   336,   336,
     336,   337,   337,   337,   338,   338,   338,   338,   338,   339,
     339,   339,   340,   340,   340,   340,   341,   341,   341,   341,
     341,   341,   341,   341,   341,   341,   341,   341,   341,   341,
     341,   341,   341,   341,   342,   342,   342,   342,   343,   344,
     344,   344,   345,   345,   346,   347,   348,   348,   349,   350,
     351,   352,   352,   353,   353,   353,   354,   354,   355,   355,
     356,   356,   357,   358,   359,   360,   360,   361,   362,   363,
     363,   365,   364,   366,   367,   367,   368,   368,   369,   369,
     369,   370,   370,   370,   371,   371,   372,   373,   373,   374,
     374,   375,   375,   376,   377,   378,   378,   379,   380,   380,
     381,   382,   383,   384,   385,   386,   386,   386,   386,   387,
     387,   388,   389,   389,   390,   391,   391,   392,   393,   393,
     394,   394,   395,   395,   396,   397,   398,   398,   398,   399,
     399,   399,   400,   401,   402,   403,   403,   403,   403,   403,
     403,   404,   404,   406,   405,   407,   408,   408,   408,   408,
     408,   409,   410,   411,   412,   412,   413,   414,   414,   414,
     415,   416,   416,   417,   417,   419,   418,   420,   420,   422,
     421,   423,   423,   423,   423,   423,   423,   423,   423,   424,
     424,   424,   424,   424,   426,   425,   427,   427,   427,   427,
     427,   428,   428,   429,   430,   430,   430,   430,   430,   430,
     430,   430,   430,   430,   430,   431,   431,   432,   432,   432,
     432,   433,   433,   434,   435,   436,   436,   437,   437,   438,
     439,   441,   440,   442,   443,   444,   444,   445,   446,   446,
     447,   448,   448,   449,   449,   451,   450,   452,   452,   453,
     453,   453,   453,   453,   453,   453,   454,   454,   455,   455,
     455,   455,   456,   457,   458,   459,   459,   460,   460,   460,
     461,   462,   462,   463,   464,   464,   465,   466,   466,   468,
     469,   467,   471,   470,   472,   473,   474,   474,   476,   477,
     475,   479,   480,   478,   481,   481,   482,   482,   483,   483,
     484,   484,   485,   485,   486,   487,   487,   489,   490,   491,
     488,   492,   493,   493,   494,   494,   495,   496,   497,   497,
     499,   498,   501,   500,   502,   502,   503,   503,   503,   504,
     505,   506,   507,   507,   507,   508,   508,   508,   510,   509,
     511,   511,   511,   512,   513,   514,   514,   515,   516,   516,
     516,   516,   516,   516,   516,   516,   516,   518,   517,   517,
     519,   517,   520,   520,   520,   520,   520,   521,   522,   523,
     524,   524,   524,   524,   524,   524,   525,   526,   527,   528,
     529,   530,   531,   531,   531,   532,   532,   533,   533,   533,
     533,   533,   533,   533,   533,   534,   534,   536,   535,   537,
     537,   538,   538,   539,   539,   539,   540,   540,   541,   542,
     543,   544,   544,   545,   545,   545,   546,   547,   547,   548,
     548,   549,   549,   550,   550,   551,   551,   552,   553,   553,
     554,   555,   555,   556,   557,   557,   558,   559,   559,   559,
     560,   561,   562,   562,   563,   564,   565,   565,   566,   567,
     568,   569,   569,   569,   569
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
       1,     1,     1,     1,     2,     4,     1,     2,     8,     1,
       1,     3,     3,     0,     3,     3,     0,     1,     0,     3,
       0,     1,     1,     2,     2,     1,     1,     2,     4,     1,
       2,     0,     5,     1,     0,     2,     1,     1,     3,     4,
       4,     1,     1,     1,     1,     1,     1,     4,     4,     1,
       2,     1,     2,     3,     4,     1,     2,     1,     1,     2,
       2,     2,     2,     3,     3,     1,     1,     1,     1,     0,
       2,     6,     1,     3,     1,     0,     2,     2,     1,     3,
       1,     1,     1,     1,     3,     5,     1,     2,     2,     1,
       1,     1,     5,     1,     1,     0,     4,     4,     4,     4,
       4,     1,     1,     0,     3,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     0,     2,     1,     3,     3,     5,
       6,     1,     2,     1,     1,     0,     7,     1,     1,     0,
       8,     3,     3,     3,     3,     3,     3,     3,     3,     1,
       3,     3,     3,     3,     0,     7,     1,     1,     1,     1,
       1,     1,     2,     1,     3,     3,     1,     3,     3,     3,
       3,     3,     3,     4,     3,     1,     1,     1,     1,     1,
       1,     3,     6,     9,    12,     3,     3,     1,     2,     2,
       1,     0,     9,     4,     1,     1,     2,     4,     1,     2,
       1,     0,     2,     1,     1,     0,     5,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     2,     1,     1,
       1,     1,     5,     5,     1,     1,     3,     1,     3,     1,
       3,     1,     1,     5,     1,     1,     4,     1,     2,     0,
       0,     7,     0,     8,     4,     1,     1,     2,     0,     0,
      14,     0,     0,    14,     0,     3,     0,     3,     0,     3,
       0,     3,     0,     3,     4,     1,     2,     0,     0,     0,
      11,     1,     0,     2,     1,     1,     3,     4,     1,     2,
       0,     5,     0,     7,     0,     3,     1,     1,     1,     3,
       3,     3,     0,     1,     2,     1,     1,     1,     0,     6,
       1,     1,     1,     3,     3,     1,     3,     2,     1,     1,
       3,     3,     3,     3,     3,     2,     4,     0,     5,     4,
       0,     5,     1,     2,     2,     3,     2,     1,     1,     2,
       1,     1,     1,     1,     1,     1,     4,     4,     4,     6,
       6,     6,     1,     1,     1,     0,     2,     1,     1,     1,
       1,     1,     1,     1,     1,     0,     2,     0,     6,     1,
       2,     0,     1,     3,     3,     3,     1,     1,     1,     3,
       4,     1,     2,     1,     1,     1,     4,     2,     0,     2,
       0,     2,     2,     1,     1,     1,     1,     4,     1,     2,
       3,     1,     1,     4,     1,     2,     3,     1,     1,     1,
       9,     3,     1,     2,     3,     1,     1,     3,     3,     3,
       3,     0,     3,     3,     3
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
#line 646 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.sym) = mdl_existing_object(parse_state, (yyvsp[0].str))); }
#line 2981 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 30:
#line 649 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.sym) = mdl_existing_region(parse_state, (yyvsp[-3].sym), (yyvsp[-1].str))); }
#line 2987 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 31:
#line 652 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.vec3) = mdl_point(parse_state, &(yyvsp[0].nlist))); }
#line 2993 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 33:
#line 656 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.vec3) = mdl_point_scalar((yyvsp[0].dbl))); }
#line 2999 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 34:
#line 659 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = 1; }
#line 3005 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 35:
#line 660 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = 0; }
#line 3011 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 36:
#line 661 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = 1; }
#line 3017 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 37:
#line 662 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = 0; }
#line 3023 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 38:
#line 663 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = 1; }
#line 3029 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 39:
#line 664 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = 0; }
#line 3035 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 40:
#line 667 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.mol_type).orient_set = 0; }
#line 3041 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 43:
#line 670 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.mol_type).orient_set = 1; (yyval.mol_type).orient = 0; }
#line 3047 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 44:
#line 674 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.mol_type).orient = 1; (yyval.mol_type).orient_set = 1; }
#line 3053 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 45:
#line 675 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.mol_type).orient = -1; (yyval.mol_type).orient_set = 1; }
#line 3059 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 46:
#line 676 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
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
#line 3074 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 47:
#line 686 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
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
#line 3089 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 50:
#line 704 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    {
                                                          (yyval.mol_type).orient = (int) (yyvsp[-1].dbl);
                                                          (yyval.mol_type).orient_set = 1;
                                                          if ((yyval.mol_type).orient != (yyvsp[-1].dbl))
                                                          {
                                                            mdlerror(parse_state, "molecule orientation specified inside braces must be an integer between -32768 and 32767.");
                                                            return 1;
                                                          }
                                                      }
#line 3103 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 52:
#line 717 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
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
#line 3119 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 53:
#line 730 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mcell_generate_range_singleton(&(yyval.nlist), (yyvsp[0].dbl))); }
#line 3125 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 54:
#line 731 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_generate_range(parse_state, &(yyval.nlist), (yyvsp[-5].dbl), (yyvsp[-3].dbl), (yyvsp[-1].dbl))); }
#line 3131 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 55:
#line 737 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
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
#line 3153 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 56:
#line 760 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_assign_variable_double(parse_state, (yyvsp[-2].sym), (yyvsp[0].dbl))); }
#line 3159 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 57:
#line 761 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_assign_variable_string(parse_state, (yyvsp[-2].sym), (yyvsp[0].str))); }
#line 3165 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 58:
#line 762 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_assign_variable(parse_state, (yyvsp[-2].sym), (yyvsp[0].sym))); }
#line 3171 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 59:
#line 763 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_assign_variable_array(parse_state, (yyvsp[-2].sym), (yyvsp[0].nlist).value_head)); }
#line 3177 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 60:
#line 766 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.sym) = mdl_get_or_create_variable(parse_state, (yyvsp[0].str))); }
#line 3183 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 61:
#line 769 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.sym) = mdl_existing_variable(parse_state, (yyvsp[0].str))); }
#line 3189 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 63:
#line 773 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    {
                                                          struct num_expr_list *elp;
                                                          (yyval.nlist).value_head = (struct num_expr_list *) (yyvsp[0].sym)->value;
                                                          (yyval.nlist).value_count = 1;
                                                          for (elp = (yyval.nlist).value_head; elp->next != NULL; elp = elp->next)
                                                            ++ (yyval.nlist).value_count;
                                                          (yyval.nlist).value_tail = elp;
                                                          (yyval.nlist).shared = 1;
                                                      }
#line 3203 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 64:
#line 784 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { mdl_debug_dump_array((yyvsp[-1].nlist).value_head); (yyval.nlist) = (yyvsp[-1].nlist); }
#line 3209 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 65:
#line 787 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.sym) = mdl_existing_array(parse_state, (yyvsp[0].str))); }
#line 3215 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 69:
#line 795 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.dbl) = *(double *) (yyvsp[0].sym)->value; }
#line 3221 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 70:
#line 798 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.dbl) = (yyvsp[0].llival); }
#line 3227 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 74:
#line 806 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.sym) = mdl_existing_double(parse_state, (yyvsp[0].str))); }
#line 3233 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 75:
#line 810 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.dbl) = (yyvsp[-1].dbl); }
#line 3239 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 76:
#line 811 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { CHECKF((yyval.dbl) = exp((yyvsp[-1].dbl))); }
#line 3245 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 77:
#line 812 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_expr_log(parse_state, (yyvsp[-1].dbl), &(yyval.dbl))); }
#line 3251 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 78:
#line 813 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_expr_log10(parse_state, (yyvsp[-1].dbl), &(yyval.dbl))); }
#line 3257 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 79:
#line 814 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.dbl) = max2d((yyvsp[-3].dbl), (yyvsp[-1].dbl)); }
#line 3263 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 80:
#line 815 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.dbl) = min2d((yyvsp[-3].dbl), (yyvsp[-1].dbl)); }
#line 3269 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 81:
#line 816 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.dbl) = mdl_expr_roundoff((yyvsp[-1].dbl), (int) (yyvsp[-3].dbl)); }
#line 3275 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 82:
#line 817 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.dbl) = floor((yyvsp[-1].dbl)); }
#line 3281 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 83:
#line 818 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.dbl) = ceil((yyvsp[-1].dbl)); }
#line 3287 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 84:
#line 819 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.dbl) = sin((yyvsp[-1].dbl)); }
#line 3293 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 85:
#line 820 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.dbl) = cos((yyvsp[-1].dbl)); }
#line 3299 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 86:
#line 821 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { CHECKF((yyval.dbl) = tan((yyvsp[-1].dbl))); }
#line 3305 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 87:
#line 822 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { CHECKF((yyval.dbl) = asin((yyvsp[-1].dbl))); }
#line 3311 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 88:
#line 823 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { CHECKF((yyval.dbl) = acos((yyvsp[-1].dbl))); }
#line 3317 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 89:
#line 824 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.dbl) = atan((yyvsp[-1].dbl)); }
#line 3323 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 90:
#line 825 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { CHECKF((yyval.dbl) = sqrt((yyvsp[-1].dbl))); }
#line 3329 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 91:
#line 826 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.dbl) = fabs((yyvsp[-1].dbl)); }
#line 3335 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 92:
#line 827 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_expr_mod(parse_state, (yyvsp[-3].dbl), (yyvsp[-1].dbl), &(yyval.dbl))); }
#line 3341 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 93:
#line 828 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.dbl) = MY_PI; }
#line 3347 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 94:
#line 829 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.dbl) = mdl_expr_rng_uniform(parse_state); }
#line 3353 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 95:
#line 830 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.dbl) = rng_gauss(parse_state->vol->rng); }
#line 3359 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 96:
#line 831 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.dbl) = parse_state->vol->seed_seq; }
#line 3365 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 97:
#line 832 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_expr_string_to_double(parse_state, (yyvsp[-1].str), &(yyval.dbl))); }
#line 3371 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 98:
#line 833 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { CHECKF((yyval.dbl) = (yyvsp[-2].dbl) + (yyvsp[0].dbl)); }
#line 3377 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 99:
#line 834 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { CHECKF((yyval.dbl) = (yyvsp[-2].dbl) - (yyvsp[0].dbl)); }
#line 3383 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 100:
#line 835 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { CHECKF((yyval.dbl) = (yyvsp[-2].dbl) * (yyvsp[0].dbl)); }
#line 3389 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 101:
#line 836 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_expr_div(parse_state, (yyvsp[-2].dbl), (yyvsp[0].dbl), &(yyval.dbl))); }
#line 3395 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 102:
#line 837 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_expr_pow(parse_state, (yyvsp[-2].dbl), (yyvsp[0].dbl), &(yyval.dbl))); }
#line 3401 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 103:
#line 838 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.dbl) = -(yyvsp[0].dbl); }
#line 3407 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 104:
#line 839 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.dbl) = (yyvsp[0].dbl); }
#line 3413 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 106:
#line 844 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.str) = mdl_strdup((char const *) (yyvsp[0].sym)->value)); }
#line 3419 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 107:
#line 848 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.str) = mdl_strip_quotes((yyvsp[0].str))); }
#line 3425 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 108:
#line 849 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.str) = mdl_strdup(parse_state->vol->mdl_infile_name)); }
#line 3431 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 109:
#line 850 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.str) = mdl_strcat((yyvsp[-2].str), (yyvsp[0].str))); }
#line 3437 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 110:
#line 851 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.str) = mdl_string_format(parse_state, (yyvsp[-2].str), (yyvsp[-1].printfargs).arg_head)); }
#line 3443 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 111:
#line 854 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.sym) = mdl_existing_string(parse_state, (yyvsp[0].str))); }
#line 3449 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 119:
#line 870 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_fopen(parse_state, (yyvsp[-6].sym), (yyvsp[-3].str), (yyvsp[-1].str))); }
#line 3455 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 120:
#line 873 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.sym) = mdl_new_filehandle(parse_state, (yyvsp[0].str))); }
#line 3461 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 121:
#line 876 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.str) = (yyvsp[0].str); CHECK(mdl_valid_file_mode(parse_state, (yyvsp[0].str))); }
#line 3467 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 122:
#line 879 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_fclose(parse_state, (yyvsp[-1].sym))); }
#line 3473 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 123:
#line 882 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.sym) = mdl_existing_file_stream(parse_state, (yyvsp[0].str))); }
#line 3479 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 124:
#line 885 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.str) = mdl_expand_string_escapes((yyvsp[0].str))); }
#line 3485 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 125:
#line 888 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.printfargs).arg_head = (yyval.printfargs).arg_tail = NULL; }
#line 3491 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 126:
#line 889 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    {
                                                        (yyval.printfargs) = (yyvsp[-2].printfargs);
                                                        if ((yyval.printfargs).arg_tail)
                                                          (yyval.printfargs).arg_tail = (yyval.printfargs).arg_tail->next = (yyvsp[0].printfarg);
                                                        else
                                                          (yyval.printfargs).arg_tail = (yyval.printfargs).arg_head = (yyvsp[0].printfarg);
                                                        (yyvsp[0].printfarg)->next = NULL;
                                                      }
#line 3504 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 127:
#line 899 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.printfarg) = mdl_new_printf_arg_double((yyvsp[0].dbl))); }
#line 3510 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 128:
#line 900 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.printfarg) = mdl_new_printf_arg_string((yyvsp[0].str))); }
#line 3516 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 129:
#line 901 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
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
#line 3531 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 130:
#line 913 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_printf(parse_state, (yyvsp[-2].str), (yyvsp[-1].printfargs).arg_head)); }
#line 3537 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 131:
#line 919 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_fprintf(parse_state, (struct file_stream *) (yyvsp[-4].sym)->value, (yyvsp[-2].str), (yyvsp[-1].printfargs).arg_head)); }
#line 3543 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 132:
#line 925 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_sprintf(parse_state, (yyvsp[-4].sym), (yyvsp[-2].str), (yyvsp[-1].printfargs).arg_head)); }
#line 3549 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 133:
#line 928 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { mdl_print_time(parse_state, (yyvsp[-1].str)); }
#line 3555 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 134:
#line 934 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_fprint_time(parse_state, (yyvsp[-3].sym), (yyvsp[-1].str))); }
#line 3561 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 138:
#line 950 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { if (!parse_state->vol->quiet_flag) mdl_set_all_notifications(parse_state->vol, (yyvsp[0].tok)); }
#line 3567 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 139:
#line 951 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { if (!parse_state->vol->quiet_flag) parse_state->vol->notify->progress_report        = (yyvsp[0].tok); }
#line 3573 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 140:
#line 952 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { if (!parse_state->vol->quiet_flag) parse_state->vol->notify->diffusion_constants    = (yyvsp[0].tok); }
#line 3579 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 141:
#line 953 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { if (!parse_state->vol->quiet_flag) parse_state->vol->notify->reaction_probabilities = (yyvsp[0].tok); }
#line 3585 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 142:
#line 954 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { if (!parse_state->vol->quiet_flag) parse_state->vol->notify->time_varying_reactions = (yyvsp[0].tok); }
#line 3591 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 143:
#line 955 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { if (!parse_state->vol->quiet_flag) parse_state->vol->notify->reaction_prob_notify   = (yyvsp[0].dbl); }
#line 3597 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 144:
#line 956 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { if (!parse_state->vol->quiet_flag) parse_state->vol->notify->partition_location     = (yyvsp[0].tok); }
#line 3603 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 145:
#line 957 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { if (!parse_state->vol->quiet_flag) parse_state->vol->notify->box_triangulation      = (yyvsp[0].tok); }
#line 3609 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 146:
#line 958 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { if (!parse_state->vol->quiet_flag) parse_state->vol->notify->release_events         = (yyvsp[0].tok); }
#line 3615 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 147:
#line 959 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { if (!parse_state->vol->quiet_flag) parse_state->vol->notify->file_writes            = (yyvsp[0].tok); }
#line 3621 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 148:
#line 960 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { if (!parse_state->vol->quiet_flag) parse_state->vol->notify->final_summary          = (yyvsp[0].tok); }
#line 3627 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 149:
#line 961 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { if (!parse_state->vol->quiet_flag) parse_state->vol->notify->throughput_report      = (yyvsp[0].tok); }
#line 3633 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 150:
#line 962 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { if (!parse_state->vol->quiet_flag) parse_state->vol->notify->reaction_output_report = (yyvsp[0].tok); }
#line 3639 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 151:
#line 963 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { if (!parse_state->vol->quiet_flag) parse_state->vol->notify->volume_output_report   = (yyvsp[0].tok); }
#line 3645 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 152:
#line 964 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { if (!parse_state->vol->quiet_flag) parse_state->vol->notify->viz_output_report      = (yyvsp[0].tok); }
#line 3651 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 153:
#line 965 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { if (!parse_state->vol->quiet_flag) parse_state->vol->notify->checkpoint_report      = (yyvsp[0].tok); }
#line 3657 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 154:
#line 966 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    {
                                                          if (!parse_state->vol->quiet_flag && parse_state->vol->log_freq == ULONG_MAX)
                                                            parse_state->vol->notify->iteration_report = (yyvsp[0].tok);
                                                      }
#line 3666 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 155:
#line 970 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { if (!parse_state->vol->quiet_flag) CHECK(mdl_set_iteration_report_freq(parse_state, (long long) (yyvsp[0].dbl))); }
#line 3672 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 156:
#line 971 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { if (!parse_state->vol->quiet_flag) parse_state->vol->notify->molecule_collision_report    = (yyvsp[0].tok); }
#line 3678 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 157:
#line 975 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = ((yyvsp[0].tok) ? NOTIFY_FULL : NOTIFY_NONE); }
#line 3684 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 158:
#line 979 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = ((yyvsp[0].tok) ? NOTIFY_FULL : NOTIFY_NONE); }
#line 3690 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 159:
#line 980 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = NOTIFY_BRIEF; }
#line 3696 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 163:
#line 996 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { mdl_set_all_warnings(parse_state->vol, (byte) (yyvsp[0].tok)); }
#line 3702 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 164:
#line 997 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->vol->notify->neg_diffusion = (byte)(yyvsp[0].tok); }
#line 3708 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 165:
#line 998 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->vol->notify->neg_reaction = (byte)(yyvsp[0].tok); }
#line 3714 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 166:
#line 999 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->vol->notify->high_reaction_prob = (byte)(yyvsp[0].tok); }
#line 3720 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 167:
#line 1000 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->vol->notify->reaction_prob_warn = (yyvsp[0].dbl); }
#line 3726 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 168:
#line 1001 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->vol->notify->close_partitions = (byte)(yyvsp[0].tok); }
#line 3732 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 169:
#line 1002 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->vol->notify->degenerate_polys = (byte)(yyvsp[0].tok); }
#line 3738 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 170:
#line 1003 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->vol->notify->overwritten_file = (byte)(yyvsp[0].tok); }
#line 3744 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 171:
#line 1004 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->vol->notify->short_lifetime = (byte)(yyvsp[0].tok); }
#line 3750 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 172:
#line 1005 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_lifetime_warning_threshold(parse_state, (long long) (yyvsp[0].dbl))); }
#line 3756 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 173:
#line 1006 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->vol->notify->missed_reactions = (byte)(yyvsp[0].tok); }
#line 3762 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 174:
#line 1007 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_missed_reaction_warning_threshold(parse_state, (yyvsp[0].dbl))); }
#line 3768 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 175:
#line 1008 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->vol->notify->missed_surf_orient = (byte)(yyvsp[0].tok); }
#line 3774 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 176:
#line 1009 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->vol->notify->useless_vol_orient = (byte)(yyvsp[0].tok); }
#line 3780 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 177:
#line 1010 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->vol->notify->mol_placement_failure = (byte) (yyvsp[0].tok); }
#line 3786 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 178:
#line 1011 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->vol->notify->invalid_output_step_time = (byte) (yyvsp[0].tok); }
#line 3792 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 179:
#line 1012 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->vol->notify->large_molecular_displacement = (byte) (yyvsp[0].tok); }
#line 3798 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 180:
#line 1013 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->vol->notify->add_remove_mesh_warning = (byte) (yyvsp[0].tok); }
#line 3804 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 181:
#line 1017 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = WARN_COPE;  }
#line 3810 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 182:
#line 1018 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = WARN_WARN;  }
#line 3816 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 183:
#line 1019 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = WARN_ERROR; }
#line 3822 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 184:
#line 1025 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_checkpoint_infile(parse_state, (yyvsp[0].str))); }
#line 3828 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 185:
#line 1026 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_checkpoint_outfile(parse_state, (yyvsp[0].str))); }
#line 3834 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 186:
#line 1027 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_checkpoint_interval(parse_state, (yyvsp[-1].dbl), (yyvsp[0].tok))); }
#line 3840 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 187:
#line 1028 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_keep_checkpoint_files(parse_state, (yyvsp[0].tok))); }
#line 3846 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 188:
#line 1030 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_realtime_checkpoint(parse_state, (long) (yyvsp[-1].dbl), (yyvsp[0].tok))); }
#line 3852 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 189:
#line 1033 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = 0; }
#line 3858 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 190:
#line 1034 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = 1; }
#line 3864 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 191:
#line 1035 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = 0; }
#line 3870 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 192:
#line 1039 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { /* seconds */     (yyval.dbl) = (yyvsp[0].dbl); }
#line 3876 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 193:
#line 1040 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { /* mm:ss */       (yyval.dbl) = (yyvsp[-2].dbl) * 60 + (yyvsp[0].dbl); }
#line 3882 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 194:
#line 1041 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { /* hh:mm:ss */    (yyval.dbl) = (yyvsp[-4].dbl) * 3600 + (yyvsp[-2].dbl) * 60 + (yyvsp[0].dbl); }
#line 3888 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 195:
#line 1043 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { /* dd:hh:mm:ss */ (yyval.dbl) = (yyvsp[-6].dbl) * 86400 + (yyvsp[-4].dbl) * 3600 + (yyvsp[-2].dbl) * 60 + (yyvsp[0].dbl); }
#line 3894 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 196:
#line 1050 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_time_step(parse_state, (yyvsp[0].dbl))); }
#line 3900 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 197:
#line 1051 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_space_step(parse_state, (yyvsp[0].dbl))); }
#line 3906 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 198:
#line 1052 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_max_time_step(parse_state, (yyvsp[0].dbl))); }
#line 3912 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 199:
#line 1053 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_num_iterations(parse_state, (long long) (yyvsp[0].dbl))); }
#line 3918 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 200:
#line 1054 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->vol->randomize_smol_pos = !((yyvsp[0].tok)); }
#line 3924 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 201:
#line 1055 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->vol->use_expanded_list = (yyvsp[0].tok); }
#line 3930 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 202:
#line 1056 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->vol->vacancy_search_dist2 = max2d((yyvsp[0].dbl), 0.0); }
#line 3936 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 203:
#line 1057 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_num_radial_directions(parse_state, (int) (yyvsp[0].dbl))); }
#line 3942 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 204:
#line 1058 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->vol->fully_random = 1; }
#line 3948 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 205:
#line 1059 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_num_radial_subdivisions(parse_state, (int) (yyvsp[0].dbl))); }
#line 3954 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 206:
#line 1060 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_grid_density(parse_state, (yyvsp[0].dbl))); }
#line 3960 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 207:
#line 1061 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_interaction_radius(parse_state, (yyvsp[0].dbl))); }
#line 3966 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 208:
#line 1062 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->vol->surface_reversibility=(yyvsp[0].tok); parse_state->vol->volume_reversibility=(yyvsp[0].tok); }
#line 3972 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 209:
#line 1063 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->vol->surface_reversibility=1;  parse_state->vol->volume_reversibility=0;  }
#line 3978 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 210:
#line 1064 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->vol->surface_reversibility=0;  parse_state->vol->volume_reversibility=1;  }
#line 3984 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 211:
#line 1065 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mcell_add_dynamic_geometry_file((yyvsp[0].str), parse_state)); }
#line 3990 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 212:
#line 1066 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->vol->dynamic_geometry_molecule_placement = 0; }
#line 3996 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 213:
#line 1067 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->vol->dynamic_geometry_molecule_placement = 1; }
#line 4002 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 214:
#line 1074 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->vol->mem_part_x = (int) (yyvsp[0].dbl); }
#line 4008 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 215:
#line 1075 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->vol->mem_part_y = (int) (yyvsp[0].dbl); }
#line 4014 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 216:
#line 1076 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->vol->mem_part_z = (int) (yyvsp[0].dbl); }
#line 4020 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 217:
#line 1077 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->vol->mem_part_pool = (int) (yyvsp[0].dbl); }
#line 4026 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 218:
#line 1081 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mcell_set_partition(parse_state->vol, (yyvsp[-2].tok), & (yyvsp[0].nlist))); }
#line 4032 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 219:
#line 1085 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = X_PARTS; }
#line 4038 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 220:
#line 1086 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = Y_PARTS; }
#line 4044 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 221:
#line 1087 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = Z_PARTS; }
#line 4050 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 224:
#line 1098 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { mdl_print_species_summary(parse_state->vol, (yyvsp[0].mcell_mol_spec)); }
#line 4056 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 225:
#line 1102 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { mdl_print_species_summaries(parse_state->vol, (yyvsp[-1].mcell_species_lst).species_head); }
#line 4062 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 226:
#line 1106 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.mcell_species_lst).species_count = 0; CHECK(mdl_add_to_species_list(&(yyval.mcell_species_lst), (yyvsp[0].mcell_mol_spec))); }
#line 4068 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 227:
#line 1107 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.mcell_species_lst) = (yyvsp[-1].mcell_species_lst); CHECK(mdl_add_to_species_list(&(yyval.mcell_species_lst), (yyvsp[0].mcell_mol_spec))); }
#line 4074 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 228:
#line 1117 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.mcell_mol_spec) = mdl_create_species(parse_state, (yyvsp[-7].str), (yyvsp[-5].diff_const).D, (yyvsp[-5].diff_const).is_2d, (yyvsp[-4].dbl), (yyvsp[-3].ival), (yyvsp[-2].dbl), (yyvsp[-1].ival) )); }
#line 4080 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 230:
#line 1123 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.sym) = mdl_new_mol_species(parse_state, (yyvsp[0].str))); }
#line 4086 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 231:
#line 1127 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.diff_const).is_2d = 0; (yyval.diff_const).D = (yyvsp[0].dbl); CHECK(mdl_check_diffusion_constant(parse_state, & (yyval.diff_const).D)); }
#line 4092 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 232:
#line 1128 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.diff_const).is_2d = 1; (yyval.diff_const).D = (yyvsp[0].dbl); CHECK(mdl_check_diffusion_constant(parse_state, & (yyval.diff_const).D)); }
#line 4098 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 233:
#line 1132 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.dbl) = 0.0; }
#line 4104 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 234:
#line 1133 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    {
                                                          if ((yyvsp[0].dbl) <= 0)
                                                          {
                                                            mdlerror_fmt(parse_state, "Requested custom time step of %.15g; custom time step must be positive.", (yyvsp[0].dbl));
                                                            return 1;
                                                          }

                                                          (yyval.dbl) = (yyvsp[0].dbl);
                                                      }
#line 4118 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 235:
#line 1142 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    {
                                                          if ((yyvsp[0].dbl) <= 0)
                                                          {
                                                            mdlerror_fmt(parse_state, "Requested custom space step of %.15g; custom space step must be positive.", (yyvsp[0].dbl));
                                                            return 1;
                                                          }

                                                          (yyval.dbl) = -(yyvsp[0].dbl);
                                                      }
#line 4132 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 236:
#line 1153 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.ival) = 0; }
#line 4138 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 237:
#line 1154 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.ival) = 1; }
#line 4144 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 238:
#line 1158 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.dbl) = 0; }
#line 4150 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 239:
#line 1159 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    {
                                                        if ((yyvsp[0].dbl) <= 0)
                                                        {
                                                          mdlerror_fmt(parse_state, "Requested maximum step length of %.15g; maximum step length must be positive.", (yyvsp[0].dbl));
                                                          return 1;
                                                        }
                                                        (yyval.dbl) = (yyvsp[0].dbl);
                                                      }
#line 4163 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 240:
#line 1169 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    {(yyval.ival) = 0;}
#line 4169 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 241:
#line 1170 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    {(yyval.ival) = 1;}
#line 4175 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 242:
#line 1173 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.sym) = mdl_existing_molecule(parse_state, (yyvsp[0].str))); }
#line 4181 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 243:
#line 1177 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.mol_type) = (yyvsp[0].mol_type); CHECKN((yyval.mol_type).mol_type = mdl_existing_surface_molecule(parse_state, (yyvsp[-1].str))); }
#line 4187 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 244:
#line 1181 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    {
                                                        (yyval.mol_type) = (yyvsp[0].mol_type);
                                                        if (! (yyval.mol_type).orient_set)
                                                          (yyval.mol_type).orient = 0;
                                                        (yyval.mol_type).mol_type = (yyvsp[-1].sym);
                                                      }
#line 4198 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 251:
#line 1216 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { mdl_start_surface_class(parse_state, (yyvsp[-1].sym)); }
#line 4204 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 252:
#line 1218 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { mdl_finish_surface_class(parse_state); }
#line 4210 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 253:
#line 1221 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.sym) = mdl_existing_surface_class(parse_state, (yyvsp[0].str))); }
#line 4216 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 258:
#line 1238 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN(mdl_assemble_surface_reaction(parse_state, (yyvsp[-2].tok), parse_state->current_surface_class, (yyvsp[0].mol_type).mol_type, (yyvsp[0].mol_type).orient)); }
#line 4222 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 259:
#line 1241 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    {
              struct sym_entry *mol_sym = retrieve_sym("ALL_MOLECULES", parse_state->vol->mol_sym_table);
              if(!(yyvsp[0].mol_type).orient_set) (yyvsp[0].mol_type).orient = 0;
              CHECKN(mdl_assemble_surface_reaction(parse_state, (yyvsp[-3].tok), parse_state->current_surface_class, mol_sym, (yyvsp[0].mol_type).orient));}
#line 4231 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 260:
#line 1247 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN(mdl_assemble_concentration_clamp_reaction(parse_state, parse_state->current_surface_class, (yyvsp[-2].mol_type).mol_type, (yyvsp[-2].mol_type).orient, (yyvsp[0].dbl))); }
#line 4237 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 261:
#line 1250 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = RFLCT; }
#line 4243 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 262:
#line 1251 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = TRANSP; }
#line 4249 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 263:
#line 1252 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = SINK; }
#line 4255 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 266:
#line 1259 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->current_surface_class->sm_dat_head = (yyvsp[0].surf_mol_dat_list).sm_head; }
#line 4261 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 267:
#line 1266 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.surf_mol_dat_list) = (yyvsp[-1].surf_mol_dat_list); }
#line 4267 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 268:
#line 1270 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.surf_mol_dat_list) = (yyvsp[-1].surf_mol_dat_list); }
#line 4273 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 269:
#line 1274 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    {
                                                          (yyvsp[0].surf_mol_dat)->quantity_type = SURFMOLDENS;
                                                          (yyval.surf_mol_dat_list).sm_tail = (yyval.surf_mol_dat_list).sm_head = (yyvsp[0].surf_mol_dat);
                                                      }
#line 4282 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 270:
#line 1279 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    {
                                                          (yyval.surf_mol_dat_list) = (yyvsp[-1].surf_mol_dat_list);
                                                          (yyvsp[0].surf_mol_dat)->quantity_type = SURFMOLDENS;
                                                          (yyval.surf_mol_dat_list).sm_tail = (yyval.surf_mol_dat_list).sm_tail->next = (yyvsp[0].surf_mol_dat);
                                                      }
#line 4292 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 271:
#line 1287 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    {
                                                          (yyvsp[0].surf_mol_dat)->quantity_type = SURFMOLNUM;
                                                          (yyval.surf_mol_dat_list).sm_tail = (yyval.surf_mol_dat_list).sm_head = (yyvsp[0].surf_mol_dat);
                                                      }
#line 4301 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 272:
#line 1292 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    {
                                                          (yyval.surf_mol_dat_list) = (yyvsp[-1].surf_mol_dat_list);
                                                          (yyvsp[0].surf_mol_dat)->quantity_type = SURFMOLNUM;
                                                          (yyval.surf_mol_dat_list).sm_tail = (yyval.surf_mol_dat_list).sm_tail->next = (yyvsp[0].surf_mol_dat);
                                                      }
#line 4311 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 273:
#line 1300 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.surf_mol_dat) = mdl_new_surf_mol_data(parse_state, &(yyvsp[-2].mol_type), (yyvsp[0].dbl))); }
#line 4317 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 283:
#line 1329 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.react_arrow).catalyst = (yyvsp[-1].mol_type); (yyval.react_arrow).flags = ARROW_CATALYTIC; }
#line 4323 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 284:
#line 1334 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.react_arrow).catalyst = (yyvsp[-1].mol_type); (yyval.react_arrow).flags = ARROW_CATALYTIC | ARROW_BIDIRECTIONAL; }
#line 4329 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 285:
#line 1339 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.react_arrow).catalyst.mol_type = NULL; (yyval.react_arrow).flags = 0; }
#line 4335 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 287:
#line 1341 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.react_arrow).catalyst.mol_type = NULL; (yyval.react_arrow).flags = ARROW_BIDIRECTIONAL; }
#line 4341 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 289:
#line 1345 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.sym) = NULL; }
#line 4347 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 290:
#line 1346 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.sym) = mdl_new_rxn_pathname(parse_state, (yyvsp[0].str))); }
#line 4353 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 291:
#line 1352 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN(mdl_assemble_reaction(parse_state, (yyvsp[-5].mol_type_list).mol_type_head, &(yyvsp[-4].mol_type), &(yyvsp[-3].react_arrow), (yyvsp[-2].mol_type_list).mol_type_head, &(yyvsp[-1].react_rates), (yyvsp[0].sym))); }
#line 4359 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 292:
#line 1355 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_reaction_player_singleton(parse_state, & (yyval.mol_type_list), & (yyvsp[0].mol_type))); }
#line 4365 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 293:
#line 1356 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.mol_type_list) = (yyvsp[-2].mol_type_list); CHECK(mdl_add_reaction_player(parse_state, & (yyval.mol_type_list), & (yyvsp[0].mol_type))); }
#line 4371 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 295:
#line 1363 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.mol_type).mol_type = NULL; }
#line 4377 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 296:
#line 1364 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.mol_type) = (yyvsp[0].mol_type); }
#line 4383 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 297:
#line 1368 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.mol_type) = (yyvsp[0].mol_type); (yyval.mol_type).mol_type = (yyvsp[-1].sym); }
#line 4389 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 298:
#line 1371 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_reaction_player_singleton(parse_state, & (yyval.mol_type_list), & (yyvsp[0].mol_type))); }
#line 4395 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 299:
#line 1372 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.mol_type_list) = (yyvsp[-2].mol_type_list); CHECK(mdl_add_reaction_player(parse_state, & (yyval.mol_type_list), & (yyvsp[0].mol_type))); }
#line 4401 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 300:
#line 1375 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.mol_type).mol_type = NULL; (yyval.mol_type).orient_set = 0; }
#line 4407 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 304:
#line 1384 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    {
                                                        if ((yyvsp[-1].react_rates).forward_rate.rate_type == RATE_UNSET)
                                                        {
                                                          mdlerror(parse_state, "invalid reaction rate specification: must specify a forward rate.");
                                                          return 1;
                                                        }

                                                        (yyval.react_rates) = (yyvsp[-1].react_rates);
                                                      }
#line 4421 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 305:
#line 1395 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
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
#line 4440 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 306:
#line 1412 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.react_rates).forward_rate = (yyvsp[0].react_rate); (yyval.react_rates).backward_rate.rate_type = RATE_UNSET; CHECK(mdl_valid_rate(parse_state, &(yyvsp[0].react_rate))); }
#line 4446 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 307:
#line 1413 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.react_rates).forward_rate = (yyvsp[0].react_rate); (yyval.react_rates).backward_rate.rate_type = RATE_UNSET; CHECK(mdl_valid_rate(parse_state, &(yyvsp[0].react_rate))); }
#line 4452 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 308:
#line 1414 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.react_rates).backward_rate = (yyvsp[0].react_rate); (yyval.react_rates).forward_rate.rate_type = RATE_UNSET; CHECK(mdl_valid_rate(parse_state, &(yyvsp[0].react_rate))); }
#line 4458 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 309:
#line 1418 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.react_rate).rate_type = RATE_CONSTANT; (yyval.react_rate).v.rate_constant = (yyvsp[0].dbl); }
#line 4464 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 310:
#line 1419 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.react_rate).rate_type = RATE_FILE; (yyval.react_rate).v.rate_file = (yyvsp[0].str); }
#line 4470 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 311:
#line 1420 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_reaction_rate_from_var(parse_state, & (yyval.react_rate), (yyvsp[0].sym))); }
#line 4476 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 312:
#line 1431 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_release_pattern(parse_state, (yyvsp[-3].sym), &(yyvsp[-1].rpat))); }
#line 4482 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 313:
#line 1434 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.sym) = mdl_new_release_pattern(parse_state, (yyvsp[0].str))); }
#line 4488 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 314:
#line 1437 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.sym) = mdl_existing_release_pattern_or_rxn_pathname(parse_state, (yyvsp[0].str))); }
#line 4494 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 315:
#line 1441 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    {
                                                        (yyval.rpat).delay = 0;
                                                        (yyval.rpat).release_interval = FOREVER;
                                                        (yyval.rpat).train_interval = FOREVER;
                                                        (yyval.rpat).train_duration = FOREVER;
                                                        (yyval.rpat).number_of_trains = 1;
                                                      }
#line 4506 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 316:
#line 1449 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.rpat) = (yyvsp[-3].rpat); (yyval.rpat).delay = (yyvsp[0].dbl) / parse_state->vol->time_unit; }
#line 4512 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 317:
#line 1451 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.rpat) = (yyvsp[-3].rpat); (yyval.rpat).release_interval = (yyvsp[0].dbl) / parse_state->vol->time_unit; }
#line 4518 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 318:
#line 1453 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.rpat) = (yyvsp[-3].rpat); (yyval.rpat).train_interval = (yyvsp[0].dbl) / parse_state->vol->time_unit; }
#line 4524 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 319:
#line 1455 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.rpat) = (yyvsp[-3].rpat); (yyval.rpat).train_duration = (yyvsp[0].dbl) / parse_state->vol->time_unit; }
#line 4530 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 320:
#line 1457 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.rpat) = (yyvsp[-3].rpat); (yyval.rpat).number_of_trains = (yyvsp[0].ival); }
#line 4536 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 321:
#line 1460 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.ival) = (int) (yyvsp[0].dbl); }
#line 4542 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 322:
#line 1461 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.ival) = INT_MAX; }
#line 4548 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 323:
#line 1468 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->current_object = parse_state->vol->root_instance; }
#line 4554 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 324:
#line 1469 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    {
                                                        check_regions(parse_state->vol->root_instance, (yyvsp[0].obj));
                                                        add_child_objects(parse_state->vol->root_instance, (yyvsp[0].obj), (yyvsp[0].obj));
                                                        parse_state->current_object = parse_state->vol->root_object;
                                                      }
#line 4564 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 325:
#line 1479 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { add_child_objects(parse_state->vol->root_object, (yyvsp[0].obj), (yyvsp[0].obj)); }
#line 4570 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 331:
#line 1495 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.sym) = mdl_start_object(parse_state, (yyvsp[0].str))); }
#line 4576 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 333:
#line 1501 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { mdl_finish_object(parse_state); }
#line 4582 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 337:
#line 1514 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { transform_translate(parse_state->vol, parse_state->current_object->t_matrix, (yyvsp[0].vec3)); }
#line 4588 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 338:
#line 1515 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { transform_scale(parse_state->current_object->t_matrix, (yyvsp[0].vec3)); }
#line 4594 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 339:
#line 1516 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_transform_rotate(parse_state, parse_state->current_object->t_matrix, (yyvsp[-2].vec3), (yyvsp[0].dbl))); }
#line 4600 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 340:
#line 1525 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    {
                                                          struct geom_object *the_object = (struct geom_object *) (yyvsp[-5].sym)->value;
                                                          the_object->object_type = META_OBJ;
                                                          add_child_objects(the_object, (yyvsp[-2].obj_list).obj_head, (yyvsp[-2].obj_list).obj_tail);
                                                          (yyval.obj) = the_object;
                                                      }
#line 4611 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 341:
#line 1534 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { mdl_object_list_singleton(& (yyval.obj_list), (yyvsp[0].obj)); }
#line 4617 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 342:
#line 1535 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.obj_list) = (yyvsp[-1].obj_list); mdl_add_object_to_list(& (yyval.obj_list), (yyvsp[0].obj)); }
#line 4623 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 345:
#line 1544 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_deep_copy_object(parse_state, (struct geom_object *) (yyvsp[-3].sym)->value, (struct geom_object *) (yyvsp[-1].sym)->value)); }
#line 4629 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 346:
#line 1546 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.obj) = (struct geom_object *) (yyvsp[-6].sym)->value; }
#line 4635 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 349:
#line 1556 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_start_release_site(parse_state, (yyvsp[-2].sym), SHAPE_UNDEFINED)); }
#line 4641 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 350:
#line 1560 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.obj) = mdl_finish_release_site(parse_state, (yyvsp[-7].sym))); }
#line 4647 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 351:
#line 1563 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_release_site_geometry_region(parse_state, parse_state->current_release_site, parse_state->current_object, (yyvsp[0].rev))); }
#line 4653 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 352:
#line 1564 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_release_site_geometry_object(parse_state, parse_state->current_release_site, (struct geom_object *) (yyvsp[0].sym)->value)); }
#line 4659 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 353:
#line 1565 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->current_release_site->release_shape = SHAPE_SPHERICAL; }
#line 4665 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 354:
#line 1566 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->current_release_site->release_shape = SHAPE_CUBIC; }
#line 4671 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 355:
#line 1567 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->current_release_site->release_shape = SHAPE_ELLIPTIC; }
#line 4677 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 356:
#line 1568 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->current_release_site->release_shape = SHAPE_RECTANGULAR; }
#line 4683 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 357:
#line 1569 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->current_release_site->release_shape = SHAPE_SPHERICAL_SHELL; }
#line 4689 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 358:
#line 1570 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    {
                                                          parse_state->current_release_site->release_shape = SHAPE_LIST;
                                                          parse_state->current_release_site->release_number_method = CONSTNUM;
                                                      }
#line 4698 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 359:
#line 1577 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.rev) = new_release_region_expr_term((yyvsp[0].sym))); }
#line 4704 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 360:
#line 1578 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.rev) = (yyvsp[-1].rev); }
#line 4710 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 361:
#line 1579 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.rev) = new_release_region_expr_binary((yyvsp[-2].rev), (yyvsp[0].rev), REXP_UNION)); }
#line 4716 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 362:
#line 1580 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.rev) = new_release_region_expr_binary((yyvsp[-2].rev), (yyvsp[0].rev), REXP_SUBTRACTION)); }
#line 4722 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 363:
#line 1581 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.rev) = new_release_region_expr_binary((yyvsp[-2].rev), (yyvsp[0].rev), REXP_INTERSECTION)); }
#line 4728 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 364:
#line 1586 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_start_release_site(parse_state, (yyvsp[-2].sym), (yyvsp[-1].tok))); }
#line 4734 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 365:
#line 1589 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.obj) = mdl_finish_release_site(parse_state, (yyvsp[-6].sym))); }
#line 4740 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 366:
#line 1592 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = SHAPE_SPHERICAL; }
#line 4746 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 367:
#line 1593 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = SHAPE_CUBIC; }
#line 4752 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 368:
#line 1594 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = SHAPE_ELLIPTIC; }
#line 4758 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 369:
#line 1595 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = SHAPE_RECTANGULAR; }
#line 4764 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 370:
#line 1596 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = SHAPE_SPHERICAL_SHELL; }
#line 4770 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 373:
#line 1604 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.sym) = mdl_existing_num_or_array(parse_state, (yyvsp[0].str))); }
#line 4776 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 374:
#line 1608 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { set_release_site_location(parse_state->vol, parse_state->current_release_site, (yyvsp[0].vec3)); }
#line 4782 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 375:
#line 1609 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_release_site_molecule(parse_state, parse_state->current_release_site, & (yyvsp[0].mol_type))); }
#line 4788 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 376:
#line 1610 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    {
                                                        if (parse_state->current_release_site->release_shape == SHAPE_LIST)
                                                        {
                                                          mdlerror(parse_state, "molecules are already specified in a list--cannot set number or density.");
                                                          return 1;
                                                        }
                                                      }
#line 4800 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 377:
#line 1617 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_release_site_diameter(parse_state, parse_state->current_release_site, (yyvsp[0].dbl) * (((yyvsp[-2].tok) == SITE_RADIUS) ? 2.0 : 1.0))); }
#line 4806 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 378:
#line 1618 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_release_site_diameter_array(parse_state, parse_state->current_release_site, (yyvsp[0].nlist).value_count, (yyvsp[0].nlist).value_head, ((yyvsp[-2].tok) == SITE_RADIUS) ? 2.0 : 1.0)); }
#line 4812 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 379:
#line 1619 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_release_site_diameter_var(parse_state, parse_state->current_release_site, ((yyvsp[-2].tok) == SITE_RADIUS) ? 2.0 : 1.0, (yyvsp[0].sym))); }
#line 4818 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 380:
#line 1620 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_release_site_periodic_box(parse_state, parse_state->current_release_site, (yyvsp[0].vec3))); }
#line 4824 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 381:
#line 1621 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_release_site_probability(parse_state, parse_state->current_release_site, (yyvsp[0].dbl))); }
#line 4830 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 382:
#line 1623 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_release_site_pattern(parse_state, parse_state->current_release_site, (yyvsp[0].sym))); }
#line 4836 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 383:
#line 1625 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_release_site_molecule_positions(parse_state, parse_state->current_release_site, & (yyvsp[-1].rsm_list))); }
#line 4842 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 384:
#line 1626 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    {CHECK(mdl_set_release_site_graph_pattern(parse_state, parse_state->current_release_site,  (yyvsp[0].str))); }
#line 4848 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 385:
#line 1630 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = SITE_DIAMETER; }
#line 4854 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 386:
#line 1631 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = SITE_RADIUS; }
#line 4860 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 391:
#line 1643 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { set_release_site_constant_number(parse_state->current_release_site, (yyvsp[0].dbl)); }
#line 4866 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 392:
#line 1646 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { set_release_site_constant_number(parse_state->current_release_site, (yyvsp[-1].dbl)); }
#line 4872 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 393:
#line 1653 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { set_release_site_gaussian_number(parse_state->current_release_site, (yyvsp[-4].dbl), (yyvsp[-1].dbl)); }
#line 4878 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 394:
#line 1661 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { set_release_site_volume_dependent_number(parse_state->current_release_site, (yyvsp[-7].dbl), (yyvsp[-4].dbl), (yyvsp[-1].dbl)); }
#line 4884 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 395:
#line 1665 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_release_site_concentration(parse_state, parse_state->current_release_site, (yyvsp[0].dbl))); }
#line 4890 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 396:
#line 1666 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(set_release_site_density(parse_state->current_release_site, (yyvsp[0].dbl))); }
#line 4896 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 397:
#line 1670 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { release_single_molecule_singleton(& (yyval.rsm_list), (yyvsp[0].rsm)); }
#line 4902 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 398:
#line 1672 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.rsm_list) = (yyvsp[-1].rsm_list); add_release_single_molecule_to_list(& (yyval.rsm_list), (yyvsp[0].rsm)); }
#line 4908 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 399:
#line 1676 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.rsm) = mdl_new_release_single_molecule(parse_state, &(yyvsp[-1].mol_type), (yyvsp[0].vec3))); }
#line 4914 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 401:
#line 1687 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    {
                                                        CHECKN((yyval.obj) = mdl_new_polygon_list(
                                                          parse_state, (yyvsp[-4].str), (yyvsp[-1].vertlist).vertex_count, (yyvsp[-1].vertlist).vertex_head,
                                                          (yyvsp[0].ecl).connection_count, (yyvsp[0].ecl).connection_head));
                                                      }
#line 4924 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 402:
#line 1696 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    {
                                                          (yyval.obj) = (struct geom_object *) (yyvsp[-3].obj);
                                                          CHECK(mdl_finish_polygon_list(parse_state, (yyval.obj)));
                                                      }
#line 4933 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 403:
#line 1702 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.vertlist) = (yyvsp[-1].vertlist); }
#line 4939 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 404:
#line 1705 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.vertlistitem) = mdl_new_vertex_list_item((yyvsp[0].vec3))); }
#line 4945 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 405:
#line 1708 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { mdl_vertex_list_singleton(& (yyval.vertlist), (yyvsp[0].vertlistitem)); }
#line 4951 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 406:
#line 1709 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.vertlist) = (yyvsp[-1].vertlist); mdl_add_vertex_to_list(& (yyval.vertlist), (yyvsp[0].vertlistitem)); }
#line 4957 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 407:
#line 1714 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.ecl) = (yyvsp[-1].ecl); }
#line 4963 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 408:
#line 1718 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { mdl_element_connection_list_singleton(& (yyval.ecl), (yyvsp[0].elem_conn)); }
#line 4969 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 409:
#line 1720 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.ecl) = (yyvsp[-1].ecl); mdl_add_element_connection_to_list(& (yyval.ecl), (yyvsp[0].elem_conn)); }
#line 4975 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 410:
#line 1723 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.elem_conn) = mdl_new_element_connection(parse_state, & (yyvsp[0].nlist))); }
#line 4981 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 415:
#line 1739 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN(parse_state->current_region = mdl_get_region(parse_state, parse_state->current_object, "REMOVED")); }
#line 4987 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 416:
#line 1741 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    {
                                                          parse_state->current_region->element_list_head = (yyvsp[-1].elem_list).elml_head;
                                                          if (parse_state->current_object->object_type == POLY_OBJ)
                                                          {
                                                            CHECK(mdl_normalize_elements(parse_state, parse_state->current_region,0));
                                                          }
                                                      }
#line 4999 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 419:
#line 1755 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = Z_POS; }
#line 5005 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 420:
#line 1756 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = Z_NEG; }
#line 5011 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 421:
#line 1757 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = Y_NEG; }
#line 5017 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 422:
#line 1758 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = Y_POS; }
#line 5023 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 423:
#line 1759 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = X_NEG; }
#line 5029 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 424:
#line 1760 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = X_POS; }
#line 5035 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 425:
#line 1761 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = ALL_SIDES; }
#line 5041 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 427:
#line 1767 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.elem_list) = (yyvsp[-1].elem_list); mdl_add_elements_to_list(& (yyval.elem_list), (yyvsp[0].elem_list).elml_head, (yyvsp[0].elem_list).elml_tail); }
#line 5047 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 430:
#line 1773 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.elem_list).elml_tail = (yyval.elem_list).elml_head = (yyvsp[0].elem_list_item); }
#line 5053 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 431:
#line 1774 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.elem_list).elml_tail = (yyval.elem_list).elml_head = (yyvsp[0].elem_list_item); }
#line 5059 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 432:
#line 1779 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.elem_list) = (yyvsp[-1].elem_list); }
#line 5065 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 433:
#line 1784 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.elem_list) = (yyvsp[-1].elem_list); mdl_set_elements_to_exclude((yyval.elem_list).elml_head); }
#line 5071 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 435:
#line 1791 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.elem_list).elml_tail = (yyval.elem_list).elml_head = (yyvsp[0].elem_list_item); }
#line 5077 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 436:
#line 1792 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.elem_list) = (yyvsp[-2].elem_list); mdl_add_elements_to_list(& (yyval.elem_list), (yyvsp[0].elem_list_item), (yyvsp[0].elem_list_item)); }
#line 5083 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 437:
#line 1795 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.elem_list_item) = new_element_list((unsigned int) (yyvsp[0].dbl), (unsigned int) (yyvsp[0].dbl))); }
#line 5089 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 438:
#line 1796 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.elem_list_item) = new_element_list((unsigned int) (yyvsp[-2].dbl), (unsigned int) (yyvsp[0].dbl))); }
#line 5095 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 439:
#line 1797 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.elem_list_item) = mdl_new_element_side(parse_state, (yyvsp[0].tok))); }
#line 5101 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 440:
#line 1800 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.elem_list_item) = mdl_new_element_previous_region(parse_state, parse_state->current_object, parse_state->current_region, (yyvsp[0].str), (yyvsp[-2].tok))); }
#line 5107 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 441:
#line 1803 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = 0; }
#line 5113 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 442:
#line 1804 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = 1; }
#line 5119 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 443:
#line 1807 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.elem_list_item) = mdl_new_element_patch(parse_state, parse_state->current_polygon, (yyvsp[-2].vec3), (yyvsp[0].vec3), (yyvsp[-4].tok))); }
#line 5125 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 444:
#line 1810 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = 0; }
#line 5131 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 445:
#line 1811 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = 1; }
#line 5137 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 449:
#line 1827 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->current_region = (yyvsp[-1].reg); }
#line 5143 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 450:
#line 1828 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_region_elements(parse_state, (yyvsp[-3].reg), (yyvsp[0].elem_list).elml_head, (yyvsp[-3].reg)->parent->object_type == POLY_OBJ)); }
#line 5149 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 451:
#line 1830 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->current_region = NULL; }
#line 5155 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 452:
#line 1838 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    {
                                                        CHECKN(mdl_new_voxel_list(parse_state, (yyvsp[-4].sym),
                                                                                  (yyvsp[-1].vertlist).vertex_count, (yyvsp[-1].vertlist).vertex_head,
                                                                                  (yyvsp[0].ecl).connection_count, (yyvsp[0].ecl).connection_head));
                                                      }
#line 5165 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 453:
#line 1844 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.obj) = (struct geom_object *) (yyvsp[-7].sym)->value; }
#line 5171 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 454:
#line 1849 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.ecl) = (yyvsp[-1].ecl); }
#line 5177 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 455:
#line 1852 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.elem_conn) = mdl_new_tet_element_connection(parse_state, & (yyvsp[0].nlist))); }
#line 5183 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 456:
#line 1856 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    {
                                                          (yyval.ecl).connection_head = (yyval.ecl).connection_tail = (yyvsp[0].elem_conn);
                                                          (yyval.ecl).connection_count = 1;
                                                      }
#line 5192 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 457:
#line 1860 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    {
                                                          (yyval.ecl) = (yyvsp[-1].ecl);
                                                          (yyval.ecl).connection_tail = (yyval.ecl).connection_tail->next = (yyvsp[0].elem_conn);
                                                          ++ (yyval.ecl).connection_count;
                                                      }
#line 5202 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 458:
#line 1870 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->vol->periodic_traditional = (yyvsp[0].tok); }
#line 5208 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 459:
#line 1873 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN(mdl_create_periodic_box(parse_state, (yyvsp[-7].vec3), (yyvsp[-5].vec3), (yyvsp[-2].tok), (yyvsp[-1].tok), (yyvsp[0].tok))); }
#line 5214 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 460:
#line 1874 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_finish_periodic_box(parse_state)); }
#line 5220 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 461:
#line 1880 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN(mdl_new_box_object(parse_state, (yyvsp[-8].sym), (yyvsp[-3].vec3), (yyvsp[-1].vec3))); }
#line 5226 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 462:
#line 1881 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_triangulate_box_object(parse_state, (yyvsp[-10].sym), parse_state->current_polygon, (yyvsp[-2].dbl))); }
#line 5232 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 463:
#line 1883 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    {
                                                          CHECK(mdl_finish_box_object(parse_state, (yyvsp[-13].sym)));
                                                          (yyval.obj) = (struct geom_object *) (yyvsp[-13].sym)->value;
                                                      }
#line 5241 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 464:
#line 1890 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = 0; }
#line 5247 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 465:
#line 1891 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = (yyvsp[0].tok); }
#line 5253 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 466:
#line 1895 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = 0; }
#line 5259 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 467:
#line 1896 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = (yyvsp[0].tok); }
#line 5265 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 468:
#line 1900 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = 0; }
#line 5271 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 469:
#line 1901 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = (yyvsp[0].tok); }
#line 5277 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 470:
#line 1905 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = 0; }
#line 5283 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 471:
#line 1906 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = (yyvsp[0].tok); }
#line 5289 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 472:
#line 1909 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.dbl) = 0.0; }
#line 5295 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 473:
#line 1910 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    {
                                                        (yyval.dbl) = (yyvsp[0].dbl);
                                                        if ((yyval.dbl) < 2.0)
                                                        {
                                                          mdlerror(parse_state, "invalid aspect ratio requested (must be greater than or equal to 2.0)");
                                                          return 1;
                                                        }
                                                      }
#line 5308 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 477:
#line 1936 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_start_existing_obj_region_def(parse_state, (yyvsp[0].sym))); }
#line 5314 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 478:
#line 1937 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->current_region = (yyvsp[-1].reg); }
#line 5320 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 479:
#line 1939 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { mdl_set_region_elements(parse_state, (yyvsp[-4].reg), (yyvsp[0].elem_list).elml_head, 1); }
#line 5326 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 480:
#line 1941 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    {
                                                          parse_state->current_region = NULL;
                                                          parse_state->current_polygon = NULL;
                                                          parse_state->current_object = parse_state->vol->root_object;
                                                      }
#line 5336 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 481:
#line 1948 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.reg) = mdl_create_region(parse_state, parse_state->current_object, (yyvsp[0].str))); }
#line 5342 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 485:
#line 1959 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { mdl_add_surf_mol_to_region(parse_state->current_region, & (yyvsp[0].surf_mol_dat_list)); }
#line 5348 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 486:
#line 1963 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { mdl_set_region_surface_class(parse_state, parse_state->current_region, (yyvsp[0].sym)); }
#line 5354 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 490:
#line 1982 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->current_region = (struct region *) (yyvsp[-1].sym)->value; }
#line 5360 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 491:
#line 1984 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->current_region = NULL; }
#line 5366 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 492:
#line 1992 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    {
                                                          parse_state->header_comment = NULL;  /* No header by default */
                                                          parse_state->exact_time_flag = 1;    /* Print exact_time column in TRIGGER output by default */
                                                      }
#line 5375 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 493:
#line 1998 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_add_reaction_output_block_to_world(parse_state, (int) (yyvsp[-4].dbl), & (yyvsp[-2].ro_otimes), & (yyvsp[-1].ro_sets))); }
#line 5381 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 494:
#line 2002 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.dbl) = COUNTBUFFERSIZE; }
#line 5387 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 495:
#line 2003 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    {
                                                          double temp_value = (yyvsp[0].dbl);
                                                          if (!(temp_value >= 1.0 && temp_value < UINT_MAX))
                                                          {
                                                            mdlerror_fmt(parse_state, "Requested buffer size of %.15g lines is invalid.  Suggested range is 100-1000000.", temp_value);
                                                            return 1;
                                                          }
                                                          (yyval.dbl) = (yyvsp[0].dbl);
                                                      }
#line 5401 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 499:
#line 2019 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.ro_otimes).type = OUTPUT_BY_STEP; (yyval.ro_otimes).step = (yyvsp[0].dbl); }
#line 5407 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 500:
#line 2023 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    {
                                                        (yyval.ro_otimes).type = OUTPUT_BY_ITERATION_LIST;
                                                        (yyval.ro_otimes).values = (yyvsp[0].nlist);
                                                      }
#line 5416 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 501:
#line 2031 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    {
                                                        (yyval.ro_otimes).type = OUTPUT_BY_TIME_LIST;
                                                        (yyval.ro_otimes).values = (yyvsp[0].nlist);
                                                      }
#line 5425 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 502:
#line 2038 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.ro_sets).set_head = (yyval.ro_sets).set_tail = NULL; }
#line 5431 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 503:
#line 2039 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.ro_sets).set_head = (yyval.ro_sets).set_tail = (yyvsp[0].ro_set); }
#line 5437 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 504:
#line 2041 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
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
#line 5452 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 506:
#line 2055 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.ro_set) = NULL; }
#line 5458 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 507:
#line 2056 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.ro_set) = NULL; }
#line 5464 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 508:
#line 2060 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    {  parse_state->count_flags = 0; }
#line 5470 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 509:
#line 2062 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.ro_set) = mdl_populate_output_set(parse_state, parse_state->header_comment, parse_state->exact_time_flag, (yyvsp[-3].ro_cols).column_head, (yyvsp[-1].tok), (yyvsp[0].str))); }
#line 5476 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 510:
#line 2066 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.str) = NULL; }
#line 5482 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 511:
#line 2067 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.str) = ((yyvsp[0].tok) ? "" : NULL); }
#line 5488 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 512:
#line 2068 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.str) = (yyvsp[0].str); }
#line 5494 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 513:
#line 2072 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->header_comment = (yyvsp[0].str); }
#line 5500 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 514:
#line 2076 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->exact_time_flag = (yyvsp[0].tok); }
#line 5506 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 516:
#line 2082 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    {
                                                          (yyval.ro_cols) = (yyvsp[-2].ro_cols);
                                                          (yyval.ro_cols).column_tail->next = (yyvsp[0].ro_cols).column_head;
                                                          (yyval.ro_cols).column_tail = (yyvsp[0].ro_cols).column_tail;
                                                      }
#line 5516 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 517:
#line 2090 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_single_count_expr(parse_state, & (yyval.ro_cols), (yyvsp[-1].cnt), (yyvsp[0].str))); }
#line 5522 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 518:
#line 2094 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.cnt) = mdl_new_oexpr_constant(parse_state, (yyvsp[0].dbl))); }
#line 5528 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 520:
#line 2096 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.cnt) = mdl_join_oexpr_tree(parse_state, (yyvsp[-1].cnt), NULL, '(')); }
#line 5534 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 521:
#line 2097 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.cnt) = mdl_join_oexpr_tree(parse_state, (yyvsp[-2].cnt),   (yyvsp[0].cnt), '+')); }
#line 5540 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 522:
#line 2098 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.cnt) = mdl_join_oexpr_tree(parse_state, (yyvsp[-2].cnt),   (yyvsp[0].cnt), '-')); }
#line 5546 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 523:
#line 2099 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.cnt) = mdl_join_oexpr_tree(parse_state, (yyvsp[-2].cnt),   (yyvsp[0].cnt), '*')); }
#line 5552 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 524:
#line 2100 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.cnt) = mdl_join_oexpr_tree(parse_state, (yyvsp[-2].cnt),   (yyvsp[0].cnt), '/')); }
#line 5558 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 525:
#line 2101 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.cnt) = mdl_join_oexpr_tree(parse_state, (yyvsp[0].cnt), NULL, '_')); }
#line 5564 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 526:
#line 2102 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.cnt) = mdl_sum_oexpr((yyvsp[-1].cnt))); }
#line 5570 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 527:
#line 2107 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->count_flags |= COUNT_PRESENT; }
#line 5576 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 528:
#line 2108 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.cnt) = (yyvsp[-1].cnt); }
#line 5582 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 529:
#line 2109 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.cnt) = mdl_new_oexpr_constant(parse_state, (yyvsp[-1].dbl))); }
#line 5588 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 530:
#line 2110 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { parse_state->count_flags |= TRIGGER_PRESENT; }
#line 5594 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 531:
#line 2111 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.cnt) = (yyvsp[-1].cnt); }
#line 5600 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 532:
#line 2114 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = FILE_OVERWRITE; }
#line 5606 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 533:
#line 2115 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = FILE_SUBSTITUTE; }
#line 5612 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 534:
#line 2116 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = FILE_APPEND; }
#line 5618 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 535:
#line 2117 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = FILE_APPEND_HEADER; }
#line 5624 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 536:
#line 2118 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = FILE_CREATE; }
#line 5630 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 538:
#line 2124 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.sym) = mdl_existing_rxn_pathname_or_molecule(parse_state, (yyvsp[0].str))); }
#line 5636 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 539:
#line 2128 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    {
                                                        (yyval.mol_type) = (yyvsp[0].mol_type);
                                                        if ((yyval.mol_type).orient > 0)
                                                          (yyval.mol_type).orient = 1;
                                                        else if ((yyval.mol_type).orient < 0)
                                                          (yyval.mol_type).orient = -1;
                                                        CHECKN((yyval.mol_type).mol_type = mdl_existing_molecule(parse_state, (yyvsp[-1].str)));
                                                      }
#line 5649 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 546:
#line 2148 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.cnt) = mdl_count_syntax_1(parse_state, (yyvsp[-3].sym), (yyvsp[-1].sym), (yyvsp[0].tok), parse_state->count_flags)); }
#line 5655 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 547:
#line 2153 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.cnt) = mdl_count_syntax_2(parse_state, (yyvsp[-3].mol_type).mol_type, (yyvsp[-3].mol_type).orient, (yyvsp[-1].sym), (yyvsp[0].tok), parse_state->count_flags)); }
#line 5661 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 548:
#line 2158 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.cnt) = mdl_count_syntax_3(parse_state, (yyvsp[-3].str), (yyvsp[-1].sym), (yyvsp[0].tok), parse_state->count_flags)); }
#line 5667 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 549:
#line 2164 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.cnt) = mdl_count_syntax_periodic_1(parse_state, (yyvsp[-5].sym), (yyvsp[-3].sym), (yyvsp[-1].vec3), (yyvsp[0].tok), parse_state->count_flags)); }
#line 5673 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 550:
#line 2168 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.cnt) = mdl_count_syntax_periodic_2(parse_state, (yyvsp[-5].mol_type).mol_type, (yyvsp[-5].mol_type).orient, (yyvsp[-3].sym), (yyvsp[-1].vec3), (yyvsp[0].tok), parse_state->count_flags)); }
#line 5679 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 551:
#line 2173 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.cnt) = mdl_count_syntax_periodic_3(parse_state, (yyvsp[-5].str), (yyvsp[-3].sym), (yyvsp[-1].vec3), (yyvsp[0].tok), parse_state->count_flags)); }
#line 5685 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 552:
#line 2176 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.sym) = NULL; }
#line 5691 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 553:
#line 2177 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.sym) = (yyvsp[0].sym); }
#line 5697 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 554:
#line 2178 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.sym) = (yyvsp[0].sym); }
#line 5703 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 555:
#line 2181 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = REPORT_NOTHING; }
#line 5709 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 556:
#line 2182 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = (yyvsp[0].tok); }
#line 5715 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 557:
#line 2185 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = REPORT_FRONT_HITS; }
#line 5721 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 558:
#line 2186 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = REPORT_BACK_HITS; }
#line 5727 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 559:
#line 2187 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = REPORT_ALL_HITS; }
#line 5733 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 560:
#line 2188 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = REPORT_FRONT_CROSSINGS; }
#line 5739 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 561:
#line 2189 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = REPORT_BACK_CROSSINGS; }
#line 5745 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 562:
#line 2190 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = REPORT_ALL_CROSSINGS; }
#line 5751 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 563:
#line 2191 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = REPORT_CONCENTRATION; }
#line 5757 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 564:
#line 2192 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = REPORT_ENCLOSED; }
#line 5763 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 565:
#line 2195 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.str) = NULL; }
#line 5769 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 566:
#line 2196 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.str) = (yyvsp[0].str); }
#line 5775 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 567:
#line 2203 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_new_viz_output_block(parse_state)); }
#line 5781 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 568:
#line 2206 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { }
#line 5787 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 571:
#line 2215 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_viz_mode(parse_state->vol->viz_blocks, CELLBLENDER_MODE)); }
#line 5793 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 572:
#line 2216 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_viz_mode(parse_state->vol->viz_blocks, (yyvsp[0].ival))); }
#line 5799 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 573:
#line 2219 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.ival) = NO_VIZ_MODE; }
#line 5805 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 574:
#line 2220 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.ival) = ASCII_MODE; }
#line 5811 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 575:
#line 2221 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.ival) = CELLBLENDER_MODE; }
#line 5817 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 577:
#line 2226 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    {
                                                        if ((yyvsp[0].frame_list).frame_head)
                                                        {
                                                          (yyvsp[0].frame_list).frame_tail->next = parse_state->vol->viz_blocks->frame_data_head;
                                                          parse_state->vol->viz_blocks->frame_data_head = (yyvsp[0].frame_list).frame_head;
                                                        }
                                                      }
#line 5829 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 579:
#line 2239 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_viz_filename_prefix(parse_state, parse_state->vol->viz_blocks, (yyvsp[0].str))); }
#line 5835 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 580:
#line 2245 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.frame_list) = (yyvsp[-1].frame_list); }
#line 5841 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 582:
#line 2251 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
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
#line 5857 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 583:
#line 2265 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.frame_list).frame_head = (yyval.frame_list).frame_tail = NULL; }
#line 5863 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 587:
#line 2277 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_viz_state(parse_state, & (yyval.ival), (yyvsp[0].dbl))); }
#line 5869 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 588:
#line 2278 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.ival) = INCLUDE_OBJ; }
#line 5875 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 591:
#line 2288 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_viz_include_molecules(parse_state, parse_state->vol->viz_blocks, (yyvsp[-1].symlist), (yyvsp[0].ival))); }
#line 5881 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 592:
#line 2289 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_set_viz_include_all_molecules(parse_state->vol->viz_blocks, (yyvsp[0].ival))); }
#line 5887 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 593:
#line 2293 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.symlist) = mdl_existing_molecule_list(parse_state, (yyvsp[0].str))); }
#line 5893 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 594:
#line 2294 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.symlist) = mdl_existing_molecules_wildcard(parse_state, (yyvsp[0].str))); }
#line 5899 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 595:
#line 2298 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_new_viz_all_times(parse_state, & (yyval.nlist))); }
#line 5905 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 597:
#line 2304 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.frame_list) = (yyvsp[-1].frame_list); }
#line 5911 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 599:
#line 2310 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
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
#line 5929 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 600:
#line 2327 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_new_viz_mol_frames(parse_state, parse_state->vol->viz_blocks, & (yyval.frame_list), OUTPUT_BY_TIME_LIST, (yyvsp[-2].tok), & (yyvsp[0].nlist))); }
#line 5935 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 601:
#line 2331 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_new_viz_all_iterations(parse_state, & (yyval.nlist))); }
#line 5941 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 603:
#line 2338 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.frame_list) = (yyvsp[-1].frame_list); }
#line 5947 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 605:
#line 2344 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
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
#line 5965 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 606:
#line 2361 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { CHECK(mdl_new_viz_mol_frames(parse_state, parse_state->vol->viz_blocks, & (yyval.frame_list), OUTPUT_BY_ITERATION_LIST, (yyvsp[-2].tok), & (yyvsp[0].nlist))); }
#line 5971 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 607:
#line 2364 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = ALL_MOL_DATA; }
#line 5977 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 608:
#line 2365 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = MOL_POS; }
#line 5983 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 609:
#line 2366 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.tok) = MOL_ORIENT; }
#line 5989 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 610:
#line 2380 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    {
                                                          struct volume_output_item *vo;
                                                          CHECKN(vo = mdl_new_volume_output_item(parse_state, (yyvsp[-6].str), & (yyvsp[-5].species_lst), (yyvsp[-4].vec3), (yyvsp[-3].vec3), (yyvsp[-2].vec3), (yyvsp[-1].otimes)));
                                                          vo->next = parse_state->vol->volume_output_head;
                                                          parse_state->vol->volume_output_head = vo;
                                                      }
#line 6000 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 611:
#line 2389 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.str) = (yyvsp[0].str); }
#line 6006 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 613:
#line 2395 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    {
                                                          (yyval.species_lst) = (yyvsp[-1].species_lst);
                                                          (yyval.species_lst).species_count += (yyvsp[0].species_lst).species_count;
                                                          (yyval.species_lst).species_tail->next = (yyvsp[0].species_lst).species_head;
                                                          (yyval.species_lst).species_tail = (yyvsp[0].species_lst).species_tail;
                                                      }
#line 6017 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 614:
#line 2404 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.species_lst) = (yyvsp[0].species_lst); }
#line 6023 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 615:
#line 2407 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
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
#line 6043 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 616:
#line 2425 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.species_lst).species_tail = (yyval.species_lst).species_head = (yyvsp[0].species_lst_item); (yyval.species_lst).species_count = 1; }
#line 6049 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 617:
#line 2427 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    {
                                                        (yyval.species_lst) = (yyvsp[-2].species_lst);
                                                        (yyval.species_lst).species_tail = (yyval.species_lst).species_tail->next = (yyvsp[0].species_lst_item);
                                                        ++ (yyval.species_lst).species_count;
                                                      }
#line 6059 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 618:
#line 2435 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.vec3) = (yyvsp[0].vec3); }
#line 6065 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 619:
#line 2439 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { (yyval.vec3) = (yyvsp[0].vec3); }
#line 6071 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 620:
#line 2443 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
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
#line 6094 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 621:
#line 2464 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.otimes) = mdl_new_output_times_default(parse_state)); }
#line 6100 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 622:
#line 2465 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.otimes) = mdl_new_output_times_step(parse_state, (yyvsp[0].dbl))); }
#line 6106 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 623:
#line 2466 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.otimes) = mdl_new_output_times_iterations(parse_state, & (yyvsp[0].nlist))); }
#line 6112 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;

  case 624:
#line 2467 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1646  */
    { CHECKN((yyval.otimes) = mdl_new_output_times_time(parse_state, & (yyvsp[0].nlist))); }
#line 6118 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
    break;


#line 6122 "/home/jczech/mcell/build/deps/mdlparse.c" /* yacc.c:1646  */
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
#line 2470 "/home/jczech/mcell/src/mdlparse.y" /* yacc.c:1906  */






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
