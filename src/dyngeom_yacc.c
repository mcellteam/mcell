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
#define YYPURE 0

/* Push parsers.  */
#define YYPUSH 0

/* Pull parsers.  */
#define YYPULL 1




/* Copy the first part of user declarations.  */
#line 1 "dyngeom_parse.y" /* yacc.c:339  */

  #include <stdlib.h>
  #include <string.h>

  #include "mcell_objects.h"
  #include "logging.h"
  #include "dyngeom_parse_extras.h"

  #define EPS_C 1e-12

  typedef void *yyscan_t;
  void yyerror(char *);
  int yylex(void);
  extern FILE *yyin;
  int yyparse(void);

  void object_list_singleton(struct object_list *head, struct object *objp) {
    objp->next = NULL;
    head->obj_tail = head->obj_head = objp;
  }

  void add_object_to_list(struct object_list *head, struct object *objp) {
    objp->next = NULL;
    head->obj_tail = head->obj_tail->next = objp;
  }

  int create_dg_parse() {
    dg_parse = (struct dyngeom_parse_vars *)malloc(sizeof(struct dyngeom_parse_vars));
    memset(dg_parse, 0, sizeof(struct dyngeom_parse_vars));
    init_top_level_objs(dg_parse);
    return 0;
  }

  int parse_dg(char *dynamic_geometry_filename) {
    FILE *fp=fopen(dynamic_geometry_filename,"r");
    if(!fp)
    {
      no_printf("Couldn't open file for reading\n");
      return 1;
    }
    yyin=fp;
    setup_root_obj_inst(dg_parse);
    yyparse();
    fclose(fp);
    return 0;
  }

#define CHECK(a)  do { if ((a) != 0) mcell_error_nodie("Parser fail: %s:%d\n", __FILE__, __LINE__); return 0; } while (0)


#line 117 "dyngeom_yacc.c" /* yacc.c:339  */

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
   by #include "dyngeom_yacc.h".  */
#ifndef YY_YY_DYNGEOM_YACC_H_INCLUDED
# define YY_YY_DYNGEOM_YACC_H_INCLUDED
/* Debug traces.  */
#ifndef YYDEBUG
# define YYDEBUG 0
#endif
#if YYDEBUG
extern int yydebug;
#endif

/* Token type.  */
#ifndef YYTOKENTYPE
# define YYTOKENTYPE
  enum yytokentype
  {
    OBJECT = 258,
    POLYGON_LIST = 259,
    VERTEX_LIST = 260,
    ELEMENT_CONNECTIONS = 261,
    INSTANTIATE = 262,
    PARTITION_X = 263,
    PARTITION_Y = 264,
    PARTITION_Z = 265,
    VAR = 266,
    REAL = 267,
    LLINTEGER = 268,
    TO = 269,
    STEP = 270,
    TRANSLATE = 271,
    SCALE = 272,
    ROTATE = 273,
    INCLUDE_ELEMENTS = 274,
    DEFINE_SURFACE_REGIONS = 275,
    UNARYMINUS = 276
  };
#endif
/* Tokens.  */
#define OBJECT 258
#define POLYGON_LIST 259
#define VERTEX_LIST 260
#define ELEMENT_CONNECTIONS 261
#define INSTANTIATE 262
#define PARTITION_X 263
#define PARTITION_Y 264
#define PARTITION_Z 265
#define VAR 266
#define REAL 267
#define LLINTEGER 268
#define TO 269
#define STEP 270
#define TRANSLATE 271
#define SCALE 272
#define ROTATE 273
#define INCLUDE_ELEMENTS 274
#define DEFINE_SURFACE_REGIONS 275
#define UNARYMINUS 276

/* Value type.  */
#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED
typedef union YYSTYPE YYSTYPE;
union YYSTYPE
{
#line 52 "dyngeom_parse.y" /* yacc.c:355  */

  int tok;
  double dbl;
  char *str;
  long long llival;
  struct sym_table *sym;
  struct vector3 *vec3;
  struct num_expr_list_head nlist;
  struct object *obj;
  struct object_list obj_list;
  struct region *reg;

#line 212 "dyngeom_yacc.c" /* yacc.c:355  */
};
# define YYSTYPE_IS_TRIVIAL 1
# define YYSTYPE_IS_DECLARED 1
#endif


extern YYSTYPE yylval;

int yyparse (void);

#endif /* !YY_YY_DYNGEOM_YACC_H_INCLUDED  */

/* Copy the second part of user declarations.  */

#line 227 "dyngeom_yacc.c" /* yacc.c:358  */

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
#define YYFINAL  20
/* YYLAST -- Last index in YYTABLE.  */
#define YYLAST   120

/* YYNTOKENS -- Number of terminals.  */
#define YYNTOKENS  37
/* YYNNTS -- Number of nonterminals.  */
#define YYNNTS  59
/* YYNRULES -- Number of rules.  */
#define YYNRULES  85
/* YYNSTATES -- Number of states.  */
#define YYNSTATES  143

/* YYTRANSLATE[YYX] -- Symbol number corresponding to YYX as returned
   by yylex, with out-of-bounds checking.  */
#define YYUNDEFTOK  2
#define YYMAXUTOK   276

#define YYTRANSLATE(YYX)                                                \
  ((unsigned int) (YYX) <= YYMAXUTOK ? yytranslate[YYX] : YYUNDEFTOK)

/* YYTRANSLATE[TOKEN-NUM] -- Symbol number corresponding to TOKEN-NUM
   as returned by yylex, without out-of-bounds checking.  */
static const yytype_uint8 yytranslate[] =
{
       0,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,    22,     2,
      35,    36,    26,    24,    30,    25,     2,    27,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,    23,     2,
       2,    21,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,    31,     2,    32,    28,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,    33,     2,    34,     2,     2,     2,     2,
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
      15,    16,    17,    18,    19,    20,    29
};

#if YYDEBUG
  /* YYRLINE[YYN] -- Source line where rule number YYN was defined.  */
static const yytype_uint16 yyrline[] =
{
       0,   121,   121,   125,   126,   130,   131,   132,   137,   140,
     144,   147,   148,   152,   153,   156,   157,   162,   176,   177,
     180,   181,   186,   185,   191,   198,   198,   209,   212,   213,
     216,   219,   222,   225,   226,   229,   233,   234,   235,   240,
     243,   240,   254,   254,   258,   261,   262,   266,   266,   271,
     272,   276,   279,   280,   286,   290,   291,   296,   300,   305,
     306,   309,   310,   315,   321,   322,   327,   327,   333,   340,
     344,   345,   346,   352,   355,   358,   359,   362,   363,   366,
     367,   370,   371,   375,   379,   380
};
#endif

#if YYDEBUG || YYERROR_VERBOSE || 0
/* YYTNAME[SYMBOL-NUM] -- String name of the symbol SYMBOL-NUM.
   First, the terminals, then, starting at YYNTOKENS, nonterminals.  */
static const char *const yytname[] =
{
  "$end", "error", "$undefined", "OBJECT", "POLYGON_LIST", "VERTEX_LIST",
  "ELEMENT_CONNECTIONS", "INSTANTIATE", "PARTITION_X", "PARTITION_Y",
  "PARTITION_Z", "VAR", "REAL", "LLINTEGER", "TO", "STEP", "TRANSLATE",
  "SCALE", "ROTATE", "INCLUDE_ELEMENTS", "DEFINE_SURFACE_REGIONS", "'='",
  "'&'", "':'", "'+'", "'-'", "'*'", "'/'", "'^'", "UNARYMINUS", "','",
  "'['", "']'", "'{'", "'}'", "'('", "')'", "$accept", "mdl_format",
  "mdl_stmt_list", "mdl_stmt", "var", "existing_object", "point",
  "point_or_num", "list_range_specs", "range_spec", "meta_object_def",
  "list_objects", "object_ref", "existing_object_ref", "$@1",
  "new_object_name", "instance_def", "$@2", "physical_object_def",
  "object_def", "new_object", "start_object", "end_object",
  "list_opt_object_cmds", "opt_object_cmd", "transformation",
  "polygon_list_def", "$@3", "@4", "vertex_list_cmd", "$@5",
  "single_vertex", "list_points", "element_connection_cmd", "$@6",
  "list_element_connections", "element_connection",
  "list_opt_polygon_object_cmds", "opt_polygon_object_cmd",
  "element_specifier_list", "element_specifier", "incl_element_list_stmt",
  "list_element_specs", "element_spec", "in_obj_define_surface_regions",
  "list_in_obj_surface_region_defs", "in_obj_surface_region_def", "$@7",
  "new_region", "partition_def", "partition_dimension", "array_value",
  "array_expr_only", "num_expr", "num_value", "intOrReal", "num_expr_only",
  "existing_num_var", "arith_expr", YY_NULLPTR
};
#endif

# ifdef YYPRINT
/* YYTOKNUM[NUM] -- (External) token number corresponding to the
   (internal) symbol number NUM (which must be that of a token).  */
static const yytype_uint16 yytoknum[] =
{
       0,   256,   257,   258,   259,   260,   261,   262,   263,   264,
     265,   266,   267,   268,   269,   270,   271,   272,   273,   274,
     275,    61,    38,    58,    43,    45,    42,    47,    94,   276,
      44,    91,    93,   123,   125,    40,    41
};
# endif

#define YYPACT_NINF -75

#define yypact_value_is_default(Yystate) \
  (!!((Yystate) == (-75)))

#define YYTABLE_NINF -31

#define yytable_value_is_error(Yytable_value) \
  0

  /* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
     STATE-NUM.  */
static const yytype_int8 yypact[] =
{
      64,   -75,   -75,   -75,   -75,   -75,     6,    64,   -75,    14,
     -75,    19,   -75,   -75,   -75,    24,   -75,   -75,   -10,    20,
     -75,   -75,   -75,     9,    10,   -75,   -75,     9,   -75,    20,
      -3,   -75,   -75,    48,    20,   -75,   -75,   -75,    56,   -75,
     -75,     1,     1,     1,   -75,    46,   -75,   -75,   -75,   -75,
     -75,   -75,   -75,    57,   -75,    28,     7,   -75,    63,    45,
      -3,   -75,    51,   -75,   -75,    43,    62,    65,   -75,   -75,
     -75,   -75,   -75,     9,     1,   -75,   -75,    10,    52,   -75,
      10,     8,    10,   -75,    72,   -75,   -75,    23,   -75,    10,
      68,   -75,   -75,   -75,   -75,   -75,   -75,    59,   -75,     1,
     -75,   -75,    27,   -75,   -75,    58,    31,   -75,   -75,     1,
      28,    60,   -75,   -75,    20,   -75,   -75,   -75,   -75,   -75,
       4,   -75,    61,   -75,   -75,   -75,    71,    75,    18,   -75,
     -75,    66,   -75,   -75,     1,    50,   -75,    79,     1,   -75,
       1,   -75,   -75
};

  /* YYDEFACT[STATE-NUM] -- Default reduction number in state STATE-NUM.
     Performed when YYTABLE does not specify something else to do.  Zero
     means the default is an error.  */
static const yytype_uint8 yydefact[] =
{
       0,    25,    70,    71,    72,     8,     0,     2,     3,    24,
      28,     0,     7,     6,    27,     0,    29,     5,     0,     0,
       1,     4,    39,     0,     0,    30,    26,     0,    31,     0,
       0,    69,    73,     0,    33,    18,    20,    21,     0,    80,
      79,     0,     0,     0,    83,     0,    13,    15,    75,    77,
      78,    76,    42,     0,    19,     0,     0,    85,     0,     0,
       0,    74,     0,    47,    40,     0,     0,     0,    32,    17,
      34,    35,     9,     0,     0,    84,    14,     0,     0,    52,
       0,     0,     0,    22,     0,    44,    45,     0,    10,     0,
      33,    36,    11,    37,    81,    12,    82,     0,    33,     0,
      43,    46,     0,    49,    51,     0,     0,    53,    54,     0,
       0,     0,    48,    50,     0,    41,    38,    23,    16,    68,
       0,    64,     0,    63,    65,    66,     0,     0,     0,    55,
      57,     0,    67,    56,     0,     0,    59,    61,     0,    58,
       0,    60,    62
};

  /* YYPGOTO[NTERM-NUM].  */
static const yytype_int8 yypgoto[] =
{
     -75,   -75,   -75,    91,     0,   -75,   -13,   -75,   -75,    40,
      83,   -75,    69,   -75,   -75,   -75,   -75,   -75,   -75,    -4,
      21,   -22,    -6,   -74,   -75,   -75,   -75,   -75,   -75,   -75,
     -75,    22,   -75,   -75,   -75,   -75,     3,   -75,   -75,   -75,
     -21,   -75,   -75,   -32,   -75,   -75,   -12,   -75,   -75,   -75,
     -75,   -23,   -75,   -39,   -75,    29,   -75,   -75,    30
};

  /* YYDEFGOTO[NTERM-NUM].  */
static const yytype_int16 yydefgoto[] =
{
      -1,     6,     7,     8,    44,    73,    85,    93,    45,    46,
      10,    34,    35,    36,    98,    11,    12,    19,    13,    14,
      15,    29,    69,    55,    70,    71,    16,    27,    79,    53,
      62,    86,    87,    64,    78,   102,   103,    90,   107,   128,
     129,   130,   135,   136,   108,   120,   121,   126,   122,    17,
      18,    88,    32,    47,    48,    49,    95,    50,    51
};

  /* YYTABLE[YYPACT[STATE-NUM]] -- What to do in state STATE-NUM.  If
     positive, shift that token.  If negative, reduce the rule whose
     number is the opposite.  If YYTABLE_NINF, syntax error.  */
static const yytype_int16 yytable[] =
{
       9,    31,    57,    58,    59,    33,    20,     9,     5,    39,
      40,    24,     5,    39,    40,     5,   106,   -30,     5,    25,
      39,    40,    41,    22,   110,    37,    41,    23,    42,     9,
      37,     5,    43,    41,     9,    84,    43,   127,   123,    30,
      28,    30,    28,    43,    65,    66,    67,    65,    66,    67,
      38,    83,   132,    52,    30,    38,    72,   100,    30,    56,
     111,   112,    68,    63,    80,   115,   104,    91,    92,    97,
     116,     1,     2,     3,     4,     5,    60,    74,    61,   104,
     138,    75,   139,    81,    77,    89,    82,    99,   105,   109,
     127,   114,   118,   140,   125,   137,   131,   134,    21,   137,
      76,   142,    26,    54,   117,   113,   141,   133,   124,   101,
      94,    96,     0,     0,   119,     0,     0,     0,     0,     0,
     119
};

static const yytype_int16 yycheck[] =
{
       0,    24,    41,    42,    43,    27,     0,     7,    11,    12,
      13,    21,    11,    12,    13,    11,    90,     3,    11,    19,
      12,    13,    25,     4,    98,    29,    25,     3,    31,    29,
      34,    11,    35,    25,    34,    74,    35,    19,    34,    31,
      33,    31,    33,    35,    16,    17,    18,    16,    17,    18,
      29,    73,    34,     5,    31,    34,    56,    34,    31,     3,
      99,    34,    34,     6,    21,    34,    89,    80,    81,    82,
     109,     7,     8,     9,    10,    11,    30,    14,    32,   102,
      30,    36,    32,    21,    33,    33,    21,    15,    20,    30,
      19,    33,    32,    14,    33,   134,    21,    31,     7,   138,
      60,   140,    19,    34,   110,   102,   138,   128,   120,    87,
      81,    81,    -1,    -1,   114,    -1,    -1,    -1,    -1,    -1,
     120
};

  /* YYSTOS[STATE-NUM] -- The (internal number of the) accessing
     symbol of state STATE-NUM.  */
static const yytype_uint8 yystos[] =
{
       0,     7,     8,     9,    10,    11,    38,    39,    40,    41,
      47,    52,    53,    55,    56,    57,    63,    86,    87,    54,
       0,    40,     4,     3,    21,    41,    47,    64,    33,    58,
      31,    88,    89,    58,    48,    49,    50,    56,    57,    12,
      13,    25,    31,    35,    41,    45,    46,    90,    91,    92,
      94,    95,     5,    66,    49,    60,     3,    90,    90,    90,
      30,    32,    67,     6,    70,    16,    17,    18,    34,    59,
      61,    62,    41,    42,    14,    36,    46,    33,    71,    65,
      21,    21,    21,    58,    90,    43,    68,    69,    88,    33,
      74,    43,    43,    44,    92,    93,    95,    43,    51,    15,
      34,    68,    72,    73,    88,    20,    60,    75,    81,    30,
      60,    90,    34,    73,    33,    34,    90,    59,    32,    41,
      82,    83,    85,    34,    83,    33,    84,    19,    76,    77,
      78,    21,    34,    77,    31,    79,    80,    90,    30,    32,
      14,    80,    90
};

  /* YYR1[YYN] -- Symbol number of symbol that rule YYN derives.  */
static const yytype_uint8 yyr1[] =
{
       0,    37,    38,    39,    39,    40,    40,    40,    41,    42,
      43,    44,    44,    45,    45,    46,    46,    47,    48,    48,
      49,    49,    51,    50,    52,    54,    53,    55,    56,    56,
      57,    58,    59,    60,    60,    61,    62,    62,    62,    64,
      65,    63,    67,    66,    68,    69,    69,    71,    70,    72,
      72,    73,    74,    74,    75,    76,    76,    77,    78,    79,
      79,    80,    80,    81,    82,    82,    84,    83,    85,    86,
      87,    87,    87,    88,    89,    90,    90,    91,    91,    92,
      92,    93,    93,    94,    95,    95
};

  /* YYR2[YYN] -- Number of symbols on the right hand side of rule YYN.  */
static const yytype_uint8 yyr2[] =
{
       0,     2,     1,     1,     2,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     3,     1,     7,     6,     1,     2,
       1,     1,     0,     7,     1,     0,     3,     1,     1,     1,
       1,     1,     1,     0,     2,     1,     3,     3,     5,     0,
       0,    10,     0,     5,     1,     1,     2,     0,     5,     1,
       2,     1,     0,     2,     1,     1,     2,     1,     5,     1,
       3,     1,     3,     4,     1,     2,     0,     5,     1,     3,
       1,     1,     1,     1,     3,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     3,     2
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
      yyerror (YY_("syntax error: cannot back up")); \
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
                  Type, Value); \
      YYFPRINTF (stderr, "\n");                                           \
    }                                                                     \
} while (0)


/*----------------------------------------.
| Print this symbol's value on YYOUTPUT.  |
`----------------------------------------*/

static void
yy_symbol_value_print (FILE *yyoutput, int yytype, YYSTYPE const * const yyvaluep)
{
  FILE *yyo = yyoutput;
  YYUSE (yyo);
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
yy_symbol_print (FILE *yyoutput, int yytype, YYSTYPE const * const yyvaluep)
{
  YYFPRINTF (yyoutput, "%s %s (",
             yytype < YYNTOKENS ? "token" : "nterm", yytname[yytype]);

  yy_symbol_value_print (yyoutput, yytype, yyvaluep);
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
yy_reduce_print (yytype_int16 *yyssp, YYSTYPE *yyvsp, int yyrule)
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
                                              );
      YYFPRINTF (stderr, "\n");
    }
}

# define YY_REDUCE_PRINT(Rule)          \
do {                                    \
  if (yydebug)                          \
    yy_reduce_print (yyssp, yyvsp, Rule); \
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
yydestruct (const char *yymsg, int yytype, YYSTYPE *yyvaluep)
{
  YYUSE (yyvaluep);
  if (!yymsg)
    yymsg = "Deleting";
  YY_SYMBOL_PRINT (yymsg, yytype, yyvaluep, yylocationp);

  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  YYUSE (yytype);
  YY_IGNORE_MAYBE_UNINITIALIZED_END
}




/* The lookahead symbol.  */
int yychar;

/* The semantic value of the lookahead symbol.  */
YYSTYPE yylval;
/* Number of syntax errors so far.  */
int yynerrs;


/*----------.
| yyparse.  |
`----------*/

int
yyparse (void)
{
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
      yychar = yylex ();
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
        case 9:
#line 140 "dyngeom_parse.y" /* yacc.c:1646  */
    { (yyval.sym) = dg_existing_object((yyvsp[0].str)); }
#line 1417 "dyngeom_yacc.c" /* yacc.c:1646  */
    break;

  case 10:
#line 144 "dyngeom_parse.y" /* yacc.c:1646  */
    {no_printf("point\n");}
#line 1423 "dyngeom_yacc.c" /* yacc.c:1646  */
    break;

  case 11:
#line 147 "dyngeom_parse.y" /* yacc.c:1646  */
    {no_printf("point_or_num\n");}
#line 1429 "dyngeom_yacc.c" /* yacc.c:1646  */
    break;

  case 12:
#line 148 "dyngeom_parse.y" /* yacc.c:1646  */
    { }
#line 1435 "dyngeom_yacc.c" /* yacc.c:1646  */
    break;

  case 14:
#line 153 "dyngeom_parse.y" /* yacc.c:1646  */
    { }
#line 1441 "dyngeom_yacc.c" /* yacc.c:1646  */
    break;

  case 15:
#line 156 "dyngeom_parse.y" /* yacc.c:1646  */
    { }
#line 1447 "dyngeom_yacc.c" /* yacc.c:1646  */
    break;

  case 16:
#line 157 "dyngeom_parse.y" /* yacc.c:1646  */
    { generate_range(&(yyval.nlist), (yyvsp[-5].dbl), (yyvsp[-3].dbl), (yyvsp[-1].dbl)); }
#line 1453 "dyngeom_yacc.c" /* yacc.c:1646  */
    break;

  case 17:
#line 166 "dyngeom_parse.y" /* yacc.c:1646  */
    {
                                                         struct object *the_object = (struct object *) (yyvsp[-5].sym)->value;
                                                         the_object->object_type = META_OBJ;
                                                         add_child_objects(the_object, (yyvsp[-2].obj_list).obj_head, (yyvsp[-2].obj_list).obj_tail);
                                                         (yyval.obj) = the_object;
                                                     }
#line 1464 "dyngeom_yacc.c" /* yacc.c:1646  */
    break;

  case 18:
#line 176 "dyngeom_parse.y" /* yacc.c:1646  */
    { object_list_singleton(& (yyval.obj_list), (yyvsp[0].obj)); }
#line 1470 "dyngeom_yacc.c" /* yacc.c:1646  */
    break;

  case 19:
#line 177 "dyngeom_parse.y" /* yacc.c:1646  */
    { (yyval.obj_list) = (yyvsp[-1].obj_list); add_object_to_list(& (yyval.obj_list), (yyvsp[0].obj)); }
#line 1476 "dyngeom_yacc.c" /* yacc.c:1646  */
    break;

  case 22:
#line 186 "dyngeom_parse.y" /* yacc.c:1646  */
    { dg_deep_copy_object((struct object *) (yyvsp[-3].sym)->value, (struct object *) (yyvsp[-1].sym)->value); }
#line 1482 "dyngeom_yacc.c" /* yacc.c:1646  */
    break;

  case 23:
#line 188 "dyngeom_parse.y" /* yacc.c:1646  */
    { (yyval.obj) = (struct object *) (yyvsp[-6].sym)->value; }
#line 1488 "dyngeom_yacc.c" /* yacc.c:1646  */
    break;

  case 24:
#line 191 "dyngeom_parse.y" /* yacc.c:1646  */
    {no_printf("new_object_name\n");}
#line 1494 "dyngeom_yacc.c" /* yacc.c:1646  */
    break;

  case 25:
#line 198 "dyngeom_parse.y" /* yacc.c:1646  */
    { no_printf("INSTANTIATE\n"); dg_parse->current_object = dg_parse->root_instance; }
#line 1500 "dyngeom_yacc.c" /* yacc.c:1646  */
    break;

  case 26:
#line 199 "dyngeom_parse.y" /* yacc.c:1646  */
    {
                                                        no_printf("meta_object_def\n");
                                                        add_child_objects(dg_parse->root_instance, (yyvsp[0].obj), (yyvsp[0].obj));
                                                        dg_parse->current_object = dg_parse->root_object;
                                                     }
#line 1510 "dyngeom_yacc.c" /* yacc.c:1646  */
    break;

  case 27:
#line 209 "dyngeom_parse.y" /* yacc.c:1646  */
    {add_child_objects(dg_parse->root_object, (yyvsp[0].obj), (yyvsp[0].obj)); no_printf("physical_object_def\n");}
#line 1516 "dyngeom_yacc.c" /* yacc.c:1646  */
    break;

  case 30:
#line 216 "dyngeom_parse.y" /* yacc.c:1646  */
    {no_printf("new_object\n"); (yyval.sym) = dg_start_object(dg_parse, (yyvsp[0].str));}
#line 1522 "dyngeom_yacc.c" /* yacc.c:1646  */
    break;

  case 31:
#line 219 "dyngeom_parse.y" /* yacc.c:1646  */
    {no_printf("start_object\n");}
#line 1528 "dyngeom_yacc.c" /* yacc.c:1646  */
    break;

  case 32:
#line 222 "dyngeom_parse.y" /* yacc.c:1646  */
    {no_printf("end_object\n"); dg_finish_object(dg_parse);}
#line 1534 "dyngeom_yacc.c" /* yacc.c:1646  */
    break;

  case 36:
#line 233 "dyngeom_parse.y" /* yacc.c:1646  */
    {no_printf("TRANSLATE\n");}
#line 1540 "dyngeom_yacc.c" /* yacc.c:1646  */
    break;

  case 37:
#line 234 "dyngeom_parse.y" /* yacc.c:1646  */
    {no_printf("SCALE\n");}
#line 1546 "dyngeom_yacc.c" /* yacc.c:1646  */
    break;

  case 38:
#line 235 "dyngeom_parse.y" /* yacc.c:1646  */
    {no_printf("ROTATE\n");}
#line 1552 "dyngeom_yacc.c" /* yacc.c:1646  */
    break;

  case 39:
#line 240 "dyngeom_parse.y" /* yacc.c:1646  */
    {no_printf("POLYGON_LIST\n");}
#line 1558 "dyngeom_yacc.c" /* yacc.c:1646  */
    break;

  case 40:
#line 243 "dyngeom_parse.y" /* yacc.c:1646  */
    { (yyval.obj) = dg_new_polygon_list(dg_parse, (yyvsp[-5].str)); }
#line 1564 "dyngeom_yacc.c" /* yacc.c:1646  */
    break;

  case 41:
#line 247 "dyngeom_parse.y" /* yacc.c:1646  */
    { 
                                                       (yyval.obj) = (struct object *) (yyvsp[-3].obj);
                                                       dg_finish_polygon_list(dg_parse, (yyval.obj));
                                                     }
#line 1573 "dyngeom_yacc.c" /* yacc.c:1646  */
    break;

  case 42:
#line 254 "dyngeom_parse.y" /* yacc.c:1646  */
    {no_printf("vertex_list_command\n");}
#line 1579 "dyngeom_yacc.c" /* yacc.c:1646  */
    break;

  case 44:
#line 258 "dyngeom_parse.y" /* yacc.c:1646  */
    {no_printf("single_vertex\n");}
#line 1585 "dyngeom_yacc.c" /* yacc.c:1646  */
    break;

  case 45:
#line 261 "dyngeom_parse.y" /* yacc.c:1646  */
    {no_printf("list_points\n");}
#line 1591 "dyngeom_yacc.c" /* yacc.c:1646  */
    break;

  case 47:
#line 266 "dyngeom_parse.y" /* yacc.c:1646  */
    {no_printf("element_connection_cmd\n");}
#line 1597 "dyngeom_yacc.c" /* yacc.c:1646  */
    break;

  case 49:
#line 271 "dyngeom_parse.y" /* yacc.c:1646  */
    {}
#line 1603 "dyngeom_yacc.c" /* yacc.c:1646  */
    break;

  case 50:
#line 273 "dyngeom_parse.y" /* yacc.c:1646  */
    {}
#line 1609 "dyngeom_yacc.c" /* yacc.c:1646  */
    break;

  case 51:
#line 276 "dyngeom_parse.y" /* yacc.c:1646  */
    {}
#line 1615 "dyngeom_yacc.c" /* yacc.c:1646  */
    break;

  case 56:
#line 292 "dyngeom_parse.y" /* yacc.c:1646  */
    {no_printf("element_specifier\n");}
#line 1621 "dyngeom_yacc.c" /* yacc.c:1646  */
    break;

  case 58:
#line 301 "dyngeom_parse.y" /* yacc.c:1646  */
    {no_printf("incl_element_list_stmt\n");}
#line 1627 "dyngeom_yacc.c" /* yacc.c:1646  */
    break;

  case 61:
#line 309 "dyngeom_parse.y" /* yacc.c:1646  */
    {no_printf("element_spec\n");}
#line 1633 "dyngeom_yacc.c" /* yacc.c:1646  */
    break;

  case 66:
#line 327 "dyngeom_parse.y" /* yacc.c:1646  */
    { dg_parse->current_region = (yyvsp[-1].reg); }
#line 1639 "dyngeom_yacc.c" /* yacc.c:1646  */
    break;

  case 67:
#line 329 "dyngeom_parse.y" /* yacc.c:1646  */
    { dg_parse->current_region = NULL; }
#line 1645 "dyngeom_yacc.c" /* yacc.c:1646  */
    break;

  case 68:
#line 333 "dyngeom_parse.y" /* yacc.c:1646  */
    {dg_create_region(dg_parse->reg_sym_table, dg_parse->current_object, (yyvsp[0].str)); no_printf("new_region\n");}
#line 1651 "dyngeom_yacc.c" /* yacc.c:1646  */
    break;

  case 70:
#line 344 "dyngeom_parse.y" /* yacc.c:1646  */
    { (yyval.tok) = X_PARTS; }
#line 1657 "dyngeom_yacc.c" /* yacc.c:1646  */
    break;

  case 71:
#line 345 "dyngeom_parse.y" /* yacc.c:1646  */
    { (yyval.tok) = Y_PARTS; }
#line 1663 "dyngeom_yacc.c" /* yacc.c:1646  */
    break;

  case 72:
#line 346 "dyngeom_parse.y" /* yacc.c:1646  */
    { (yyval.tok) = Z_PARTS; }
#line 1669 "dyngeom_yacc.c" /* yacc.c:1646  */
    break;

  case 73:
#line 352 "dyngeom_parse.y" /* yacc.c:1646  */
    { }
#line 1675 "dyngeom_yacc.c" /* yacc.c:1646  */
    break;

  case 74:
#line 355 "dyngeom_parse.y" /* yacc.c:1646  */
    { (yyval.nlist) = (yyvsp[-1].nlist); }
#line 1681 "dyngeom_yacc.c" /* yacc.c:1646  */
    break;

  case 75:
#line 358 "dyngeom_parse.y" /* yacc.c:1646  */
    { }
#line 1687 "dyngeom_yacc.c" /* yacc.c:1646  */
    break;

  case 77:
#line 362 "dyngeom_parse.y" /* yacc.c:1646  */
    { }
#line 1693 "dyngeom_yacc.c" /* yacc.c:1646  */
    break;

  case 78:
#line 363 "dyngeom_parse.y" /* yacc.c:1646  */
    { (yyval.dbl) = *(double *) (yyvsp[0].sym)->value; }
#line 1699 "dyngeom_yacc.c" /* yacc.c:1646  */
    break;

  case 79:
#line 366 "dyngeom_parse.y" /* yacc.c:1646  */
    { (yyval.dbl) = (yyvsp[0].llival); }
#line 1705 "dyngeom_yacc.c" /* yacc.c:1646  */
    break;

  case 80:
#line 367 "dyngeom_parse.y" /* yacc.c:1646  */
    { }
#line 1711 "dyngeom_yacc.c" /* yacc.c:1646  */
    break;

  case 83:
#line 375 "dyngeom_parse.y" /* yacc.c:1646  */
    { no_printf("existing_num_var\n"); }
#line 1717 "dyngeom_yacc.c" /* yacc.c:1646  */
    break;

  case 84:
#line 379 "dyngeom_parse.y" /* yacc.c:1646  */
    { (yyval.dbl) = (yyvsp[-1].dbl); }
#line 1723 "dyngeom_yacc.c" /* yacc.c:1646  */
    break;

  case 85:
#line 380 "dyngeom_parse.y" /* yacc.c:1646  */
    { (yyval.dbl) = -(yyvsp[0].dbl); }
#line 1729 "dyngeom_yacc.c" /* yacc.c:1646  */
    break;


#line 1733 "dyngeom_yacc.c" /* yacc.c:1646  */
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
      yyerror (YY_("syntax error"));
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
        yyerror (yymsgp);
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
                      yytoken, &yylval);
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
                  yystos[yystate], yyvsp);
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
  yyerror (YY_("memory exhausted"));
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
                  yytoken, &yylval);
    }
  /* Do not reclaim the symbols of the rule whose action triggered
     this YYABORT or YYACCEPT.  */
  YYPOPSTACK (yylen);
  YY_STACK_PRINT (yyss, yyssp);
  while (yyssp != yyss)
    {
      yydestruct ("Cleanup: popping",
                  yystos[*yyssp], yyvsp);
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
#line 384 "dyngeom_parse.y" /* yacc.c:1906  */


void yyerror(char *s) {
  mcell_error("%s\n", s);
}

/*int main(int argc, char *argv[])*/
/*{*/
/*  create_dg_parse();*/
/*  parse_dg(argv[1]);*/
/*}*/
