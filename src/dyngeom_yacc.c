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
#define yyparse         dgparse
#define yylex           dglex
#define yyerror         dgerror
#define yydebug         dgdebug
#define yynerrs         dgnerrs


/* Copy the first part of user declarations.  */
#line 1 "dyngeom_parse.y" /* yacc.c:339  */

  /*#define YYDEBUG 1*/
  #include <stdlib.h>
  #include <string.h>

  #include "mcell_misc.h"
  #include "mcell_objects.h"
  #include "logging.h"
  #include "strfunc.h"
  #include "dyngeom_parse_extras.h"

  typedef void *yyscan_t;
  #include "dyngeom_yacc.h"

  int dglex_init(yyscan_t *ptr_yy_globals) ;
  int dglex_destroy(yyscan_t yyscanner);
  void dgrestart(FILE *infile, yyscan_t scanner);
  int dglex(YYSTYPE *yylval, struct dyngeom_parse_vars *dg_parse, yyscan_t scanner);

  void dgerror(struct dyngeom_parse_vars *dg_parse, yyscan_t scanner, char const *str);

  void object_list_singleton(struct object_list *head, struct object *objp) {
    objp->next = NULL;
    head->obj_tail = head->obj_head = objp;
  }

  void add_object_to_list(struct object_list *head, struct object *objp) {
    objp->next = NULL;
    head->obj_tail = head->obj_tail->next = objp;
  }

  struct dyngeom_parse_vars * create_dg_parse(struct volume *state) {
    struct dyngeom_parse_vars *dg_parse = (struct dyngeom_parse_vars *)malloc(sizeof(struct dyngeom_parse_vars));
    memset(dg_parse, 0, sizeof(struct dyngeom_parse_vars));
    return dg_parse;
  }

  int parse_dg_init(
      struct dyngeom_parse_vars *dg_parse,
      char *dynamic_geometry_filename,
      struct volume *state) {
    dg_parse->include_stack_ptr = 0;
    dg_parse->line_num[0] = 0;
    setup_root_obj_inst(dg_parse, state);
    parse_dg(dg_parse, dynamic_geometry_filename);
    return 0;
  }

  int parse_dg(struct dyngeom_parse_vars *dg_parse,
               char *dynamic_geometry_filename) {
    int cur_stack = dg_parse->include_stack_ptr ++;
    dg_parse->line_num[cur_stack] = 1;
    dg_parse->include_filename[cur_stack] = dynamic_geometry_filename;

    // Open DG MDL
    no_printf("Opening file %s\n", dynamic_geometry_filename);
    FILE *infile=fopen(dynamic_geometry_filename,"r");
    if(!infile)
    {
      no_printf("Couldn't open file for reading\n");
      return 1;
    }

    yyscan_t scanner;
    dglex_init(&scanner);
    dgrestart(infile, scanner);

    // Parse DG MDL
    char const *prev_file;
    prev_file = dg_parse->curr_file;
    dg_parse->curr_file = dynamic_geometry_filename;
    /*extern int dgdebug;*/
    /*dgdebug = 1;*/
    dgparse(dg_parse, scanner);
    dg_parse->curr_file = prev_file;
    -- dg_parse->include_stack_ptr;

    /* Free leftover object names */
    struct name_list *nl;
    while (dg_parse->object_name_list != NULL)
    {
      nl = dg_parse->object_name_list->next;
      free(dg_parse->object_name_list);
      dg_parse->object_name_list = nl;
    }

    /* Clean up! */
    fclose(infile);
    dglex_destroy(scanner);

    return 0;
  }

#define CHECK(a)  do { if ((a) != 0) mcell_error_nodie("Parser fail: %s:%d\n", __FILE__, __LINE__); return 0; } while (0)


#line 169 "dyngeom_yacc.c" /* yacc.c:339  */

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
#ifndef YY_DG_DYNGEOM_YACC_H_INCLUDED
# define YY_DG_DYNGEOM_YACC_H_INCLUDED
/* Debug traces.  */
#ifndef YYDEBUG
# define YYDEBUG 0
#endif
#if YYDEBUG
extern int dgdebug;
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
    INCLUDE_FILE = 275,
    DEFINE_SURFACE_REGIONS = 276,
    STR_VALUE = 277,
    UNARYMINUS = 278
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
#define INCLUDE_FILE 275
#define DEFINE_SURFACE_REGIONS 276
#define STR_VALUE 277
#define UNARYMINUS 278

/* Value type.  */
#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED

union YYSTYPE
{
#line 99 "dyngeom_parse.y" /* yacc.c:355  */

  int tok;
  double dbl;
  char *str;
  long long llival;
  struct sym_entry *sym;
  struct vector3 *vec3;
  struct num_expr_list_head nlist;
  struct object *obj;
  struct object_list obj_list;
  struct region *reg;

#line 268 "dyngeom_yacc.c" /* yacc.c:355  */
};

typedef union YYSTYPE YYSTYPE;
# define YYSTYPE_IS_TRIVIAL 1
# define YYSTYPE_IS_DECLARED 1
#endif



int dgparse (struct dyngeom_parse_vars *dg_parse, yyscan_t scanner);

#endif /* !YY_DG_DYNGEOM_YACC_H_INCLUDED  */

/* Copy the second part of user declarations.  */

#line 284 "dyngeom_yacc.c" /* yacc.c:358  */

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
#define YYFINAL  23
/* YYLAST -- Last index in YYTABLE.  */
#define YYLAST   129

/* YYNTOKENS -- Number of terminals.  */
#define YYNTOKENS  39
/* YYNNTS -- Number of nonterminals.  */
#define YYNNTS  63
/* YYNRULES -- Number of rules.  */
#define YYNRULES  91
/* YYNSTATES -- Number of states.  */
#define YYNSTATES  152

/* YYTRANSLATE[YYX] -- Symbol number corresponding to YYX as returned
   by yylex, with out-of-bounds checking.  */
#define YYUNDEFTOK  2
#define YYMAXUTOK   278

#define YYTRANSLATE(YYX)                                                \
  ((unsigned int) (YYX) <= YYMAXUTOK ? yytranslate[YYX] : YYUNDEFTOK)

/* YYTRANSLATE[TOKEN-NUM] -- Symbol number corresponding to TOKEN-NUM
   as returned by yylex, without out-of-bounds checking.  */
static const yytype_uint8 yytranslate[] =
{
       0,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,    24,     2,
      37,    38,    28,    26,    32,    27,     2,    29,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,    25,     2,
       2,    23,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,    33,     2,    34,    30,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,    35,     2,    36,     2,     2,     2,     2,
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
      15,    16,    17,    18,    19,    20,    21,    22,    31
};

#if YYDEBUG
  /* YYRLINE[YYN] -- Source line where rule number YYN was defined.  */
static const yytype_uint16 yyrline[] =
{
       0,   182,   182,   186,   187,   191,   192,   193,   194,   199,
     202,   205,   209,   212,   213,   217,   218,   221,   222,   227,
     241,   242,   245,   246,   251,   250,   256,   263,   263,   274,
     277,   278,   281,   284,   287,   290,   291,   294,   298,   299,
     300,   305,   308,   305,   316,   316,   320,   323,   324,   328,
     328,   333,   334,   338,   341,   342,   348,   352,   353,   358,
     362,   367,   368,   371,   372,   377,   383,   384,   389,   389,
     395,   402,   406,   407,   408,   413,   435,   438,   441,   442,
     445,   446,   449,   450,   453,   454,   458,   462,   463,   467,
     471,   472
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
  "SCALE", "ROTATE", "INCLUDE_ELEMENTS", "INCLUDE_FILE",
  "DEFINE_SURFACE_REGIONS", "STR_VALUE", "'='", "'&'", "':'", "'+'", "'-'",
  "'*'", "'/'", "'^'", "UNARYMINUS", "','", "'['", "']'", "'{'", "'}'",
  "'('", "')'", "$accept", "mdl_format", "mdl_stmt_list", "mdl_stmt",
  "str_value", "var", "existing_object", "point", "point_or_num",
  "list_range_specs", "range_spec", "meta_object_def", "list_objects",
  "object_ref", "existing_object_ref", "$@1", "new_object_name",
  "instance_def", "$@2", "physical_object_def", "object_def", "new_object",
  "start_object", "end_object", "list_opt_object_cmds", "opt_object_cmd",
  "transformation", "polygon_list_def", "$@3", "@4", "vertex_list_cmd",
  "$@5", "single_vertex", "list_points", "element_connection_cmd", "$@6",
  "list_element_connections", "element_connection",
  "list_opt_polygon_object_cmds", "opt_polygon_object_cmd",
  "element_specifier_list", "element_specifier", "incl_element_list_stmt",
  "list_element_specs", "element_spec", "in_obj_define_surface_regions",
  "list_in_obj_surface_region_defs", "in_obj_surface_region_def", "$@7",
  "new_region", "partition_def", "partition_dimension", "include_stmt",
  "array_value", "array_expr_only", "num_expr", "num_value", "intOrReal",
  "num_expr_only", "existing_num_var", "arith_expr", "str_expr",
  "str_expr_only", YY_NULLPTR
};
#endif

# ifdef YYPRINT
/* YYTOKNUM[NUM] -- (External) token number corresponding to the
   (internal) symbol number NUM (which must be that of a token).  */
static const yytype_uint16 yytoknum[] =
{
       0,   256,   257,   258,   259,   260,   261,   262,   263,   264,
     265,   266,   267,   268,   269,   270,   271,   272,   273,   274,
     275,   276,   277,    61,    38,    58,    43,    45,    42,    47,
      94,   278,    44,    91,    93,   123,   125,    40,    41
};
# endif

#define YYPACT_NINF -49

#define yypact_value_is_default(Yystate) \
  (!!((Yystate) == (-49)))

#define YYTABLE_NINF -33

#define yytable_value_is_error(Yytable_value) \
  0

  /* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
     STATE-NUM.  */
static const yytype_int8 yypact[] =
{
      38,   -49,   -49,   -49,   -49,   -49,    -1,    27,    38,   -49,
      37,   -49,    48,   -49,   -49,   -49,    58,   -49,   -49,    49,
     -49,    63,    54,   -49,   -49,   -49,    42,    45,   -49,   -49,
     -49,   -49,    55,   -49,    42,   -49,    63,     1,   -49,   -49,
      54,    75,    63,   -49,   -49,   -49,    78,   -49,   -49,     4,
       4,     4,   -49,    32,   -49,   -49,   -49,   -49,   -49,   -49,
     -49,   -49,    76,   -49,    -7,    -5,   -49,    69,    46,     1,
     -49,    50,   -49,   -49,    64,    66,    67,   -49,   -49,   -49,
     -49,   -49,    42,     4,   -49,   -49,    45,    51,   -49,    45,
       6,    45,   -49,    77,   -49,   -49,   -13,   -49,    45,    70,
     -49,   -49,   -49,   -49,   -49,   -49,    61,   -49,     4,   -49,
     -49,    26,   -49,   -49,    59,     8,   -49,   -49,     4,    -7,
      62,   -49,   -49,    63,   -49,   -49,   -49,   -49,   -49,    -4,
     -49,    65,   -49,   -49,   -49,    79,    74,    18,   -49,   -49,
      71,   -49,   -49,     4,    39,   -49,    88,     4,   -49,     4,
     -49,   -49
};

  /* YYDEFACT[STATE-NUM] -- Default reduction number in state STATE-NUM.
     Performed when YYTABLE does not specify something else to do.  Zero
     means the default is an error.  */
static const yytype_uint8 yydefact[] =
{
       0,    27,    72,    73,    74,    10,     0,     0,     2,     3,
      26,    30,     0,     8,     7,    29,     0,    31,     6,     0,
       5,     0,     0,     1,     4,    41,     0,     0,    32,    28,
       9,    90,    75,    89,     0,    33,     0,     0,    71,    76,
       0,     0,    35,    20,    22,    23,     0,    83,    82,     0,
       0,     0,    86,     0,    15,    17,    78,    80,    81,    79,
      91,    44,     0,    21,     0,     0,    88,     0,     0,     0,
      77,     0,    49,    42,     0,     0,     0,    34,    19,    36,
      37,    11,     0,     0,    87,    16,     0,     0,    54,     0,
       0,     0,    24,     0,    46,    47,     0,    12,     0,    35,
      38,    13,    39,    84,    14,    85,     0,    35,     0,    45,
      48,     0,    51,    53,     0,     0,    55,    56,     0,     0,
       0,    50,    52,     0,    43,    40,    25,    18,    70,     0,
      66,     0,    65,    67,    68,     0,     0,     0,    57,    59,
       0,    69,    58,     0,     0,    61,    63,     0,    60,     0,
      62,    64
};

  /* YYPGOTO[NTERM-NUM].  */
static const yytype_int8 yypgoto[] =
{
     -49,   -49,   -49,    95,   -49,     0,   -49,   -22,   -49,   -49,
      36,    85,   -49,    68,   -49,   -49,   -49,   -49,   -49,   -49,
      14,    15,   -29,   -12,   -44,   -49,   -49,   -49,   -49,   -49,
     -49,   -49,    12,   -49,   -49,   -49,   -49,    -2,   -49,   -49,
     -49,   -26,   -49,   -49,   -35,   -49,   -49,   -16,   -49,   -49,
     -49,   -49,   -49,   -23,   -49,   -48,   -49,    24,   -49,   -49,
      25,    80,   -49
};

  /* YYDEFGOTO[NTERM-NUM].  */
static const yytype_int16 yydefgoto[] =
{
      -1,     7,     8,     9,    31,    52,    82,    94,   102,    53,
      54,    11,    42,    43,    44,   107,    12,    13,    21,    14,
      15,    16,    36,    78,    64,    79,    80,    17,    34,    88,
      62,    71,    95,    96,    73,    87,   111,   112,    99,   116,
     137,   138,   139,   144,   145,   117,   129,   130,   135,   131,
      18,    19,    20,    97,    39,    55,    56,    57,   104,    58,
      59,    32,    33
};

  /* YYTABLE[YYPACT[STATE-NUM]] -- What to do in state STATE-NUM.  If
     positive, shift that token.  If negative, reduce the rule whose
     number is the opposite.  If YYTABLE_NINF, syntax error.  */
static const yytype_int16 yytable[] =
{
      10,    66,    67,    68,    38,    41,     5,     5,    10,    74,
      75,    76,     5,    47,    48,     5,    47,    48,    47,    48,
      37,    28,    22,   109,    74,    75,    76,    23,    49,    77,
      35,    49,   132,    49,    50,    93,    10,   136,    51,    37,
     -32,    51,    10,    51,   124,     1,     2,     3,     4,     5,
      45,    46,    25,    92,   141,   115,    45,    46,     6,    37,
     120,    26,   121,   119,    69,    81,    70,   100,   101,   106,
     125,   147,    27,   148,     5,   113,    30,    35,    37,    40,
      61,    65,    72,    83,    84,    86,    98,    89,   113,    90,
      91,   114,   108,   118,   123,   146,   127,   140,   136,   146,
     134,   151,   149,    24,   143,    85,    29,   126,   110,   122,
      63,   142,   150,   133,   103,   105,     0,     0,     0,     0,
      60,     0,     0,   128,     0,     0,     0,     0,     0,   128
};

static const yytype_int16 yycheck[] =
{
       0,    49,    50,    51,    27,    34,    11,    11,     8,    16,
      17,    18,    11,    12,    13,    11,    12,    13,    12,    13,
      33,    21,    23,    36,    16,    17,    18,     0,    27,    36,
      35,    27,    36,    27,    33,    83,    36,    19,    37,    33,
       3,    37,    42,    37,    36,     7,     8,     9,    10,    11,
      36,    36,     4,    82,    36,    99,    42,    42,    20,    33,
     108,     3,    36,   107,    32,    65,    34,    89,    90,    91,
     118,    32,    23,    34,    11,    98,    22,    35,    33,    24,
       5,     3,     6,    14,    38,    35,    35,    23,   111,    23,
      23,    21,    15,    32,    35,   143,    34,    23,    19,   147,
      35,   149,    14,     8,    33,    69,    21,   119,    96,   111,
      42,   137,   147,   129,    90,    90,    -1,    -1,    -1,    -1,
      40,    -1,    -1,   123,    -1,    -1,    -1,    -1,    -1,   129
};

  /* YYSTOS[STATE-NUM] -- The (internal number of the) accessing
     symbol of state STATE-NUM.  */
static const yytype_uint8 yystos[] =
{
       0,     7,     8,     9,    10,    11,    20,    40,    41,    42,
      44,    50,    55,    56,    58,    59,    60,    66,    89,    90,
      91,    57,    23,     0,    42,     4,     3,    23,    44,    50,
      22,    43,   100,   101,    67,    35,    61,    33,    92,    93,
      24,    61,    51,    52,    53,    59,    60,    12,    13,    27,
      33,    37,    44,    48,    49,    94,    95,    96,    98,    99,
     100,     5,    69,    52,    63,     3,    94,    94,    94,    32,
      34,    70,     6,    73,    16,    17,    18,    36,    62,    64,
      65,    44,    45,    14,    38,    49,    35,    74,    68,    23,
      23,    23,    61,    94,    46,    71,    72,    92,    35,    77,
      46,    46,    47,    96,    97,    99,    46,    54,    15,    36,
      71,    75,    76,    92,    21,    63,    78,    84,    32,    63,
      94,    36,    76,    35,    36,    94,    62,    34,    44,    85,
      86,    88,    36,    86,    35,    87,    19,    79,    80,    81,
      23,    36,    80,    33,    82,    83,    94,    32,    34,    14,
      83,    94
};

  /* YYR1[YYN] -- Symbol number of symbol that rule YYN derives.  */
static const yytype_uint8 yyr1[] =
{
       0,    39,    40,    41,    41,    42,    42,    42,    42,    43,
      44,    45,    46,    47,    47,    48,    48,    49,    49,    50,
      51,    51,    52,    52,    54,    53,    55,    57,    56,    58,
      59,    59,    60,    61,    62,    63,    63,    64,    65,    65,
      65,    67,    68,    66,    70,    69,    71,    72,    72,    74,
      73,    75,    75,    76,    77,    77,    78,    79,    79,    80,
      81,    82,    82,    83,    83,    84,    85,    85,    87,    86,
      88,    89,    90,    90,    90,    91,    92,    93,    94,    94,
      95,    95,    96,    96,    97,    97,    98,    99,    99,   100,
     101,   101
};

  /* YYR2[YYN] -- Number of symbols on the right hand side of rule YYN.  */
static const yytype_uint8 yyr2[] =
{
       0,     2,     1,     1,     2,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     3,     1,     7,     6,
       1,     2,     1,     1,     0,     7,     1,     0,     3,     1,
       1,     1,     1,     1,     1,     0,     2,     1,     3,     3,
       5,     0,     0,    10,     0,     5,     1,     1,     2,     0,
       5,     1,     2,     1,     0,     2,     1,     1,     2,     1,
       5,     1,     3,     1,     3,     4,     1,     2,     0,     5,
       1,     3,     1,     1,     1,     3,     1,     3,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     3,     2,     1,
       1,     3
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
      yyerror (dg_parse, scanner, YY_("syntax error: cannot back up")); \
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
                  Type, Value, dg_parse, scanner); \
      YYFPRINTF (stderr, "\n");                                           \
    }                                                                     \
} while (0)


/*----------------------------------------.
| Print this symbol's value on YYOUTPUT.  |
`----------------------------------------*/

static void
yy_symbol_value_print (FILE *yyoutput, int yytype, YYSTYPE const * const yyvaluep, struct dyngeom_parse_vars *dg_parse, yyscan_t scanner)
{
  FILE *yyo = yyoutput;
  YYUSE (yyo);
  YYUSE (dg_parse);
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
yy_symbol_print (FILE *yyoutput, int yytype, YYSTYPE const * const yyvaluep, struct dyngeom_parse_vars *dg_parse, yyscan_t scanner)
{
  YYFPRINTF (yyoutput, "%s %s (",
             yytype < YYNTOKENS ? "token" : "nterm", yytname[yytype]);

  yy_symbol_value_print (yyoutput, yytype, yyvaluep, dg_parse, scanner);
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
yy_reduce_print (yytype_int16 *yyssp, YYSTYPE *yyvsp, int yyrule, struct dyngeom_parse_vars *dg_parse, yyscan_t scanner)
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
                                              , dg_parse, scanner);
      YYFPRINTF (stderr, "\n");
    }
}

# define YY_REDUCE_PRINT(Rule)          \
do {                                    \
  if (yydebug)                          \
    yy_reduce_print (yyssp, yyvsp, Rule, dg_parse, scanner); \
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
yydestruct (const char *yymsg, int yytype, YYSTYPE *yyvaluep, struct dyngeom_parse_vars *dg_parse, yyscan_t scanner)
{
  YYUSE (yyvaluep);
  YYUSE (dg_parse);
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
yyparse (struct dyngeom_parse_vars *dg_parse, yyscan_t scanner)
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
      yychar = yylex (&yylval, dg_parse, scanner);
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
        case 11:
#line 205 "dyngeom_parse.y" /* yacc.c:1646  */
    { (yyval.sym) = dg_existing_object(dg_parse, (yyvsp[0].str)); }
#line 1491 "dyngeom_yacc.c" /* yacc.c:1646  */
    break;

  case 12:
#line 209 "dyngeom_parse.y" /* yacc.c:1646  */
    { /*no_printf("point\n");*/ }
#line 1497 "dyngeom_yacc.c" /* yacc.c:1646  */
    break;

  case 13:
#line 212 "dyngeom_parse.y" /* yacc.c:1646  */
    { no_printf("point_or_num\n"); }
#line 1503 "dyngeom_yacc.c" /* yacc.c:1646  */
    break;

  case 14:
#line 213 "dyngeom_parse.y" /* yacc.c:1646  */
    { }
#line 1509 "dyngeom_yacc.c" /* yacc.c:1646  */
    break;

  case 16:
#line 218 "dyngeom_parse.y" /* yacc.c:1646  */
    { }
#line 1515 "dyngeom_yacc.c" /* yacc.c:1646  */
    break;

  case 17:
#line 221 "dyngeom_parse.y" /* yacc.c:1646  */
    { }
#line 1521 "dyngeom_yacc.c" /* yacc.c:1646  */
    break;

  case 18:
#line 222 "dyngeom_parse.y" /* yacc.c:1646  */
    { }
#line 1527 "dyngeom_yacc.c" /* yacc.c:1646  */
    break;

  case 19:
#line 231 "dyngeom_parse.y" /* yacc.c:1646  */
    {
                                                         struct object *the_object = (struct object *) (yyvsp[-5].sym)->value;
                                                         the_object->object_type = META_OBJ;
                                                         add_child_objects(the_object, (yyvsp[-2].obj_list).obj_head, (yyvsp[-2].obj_list).obj_tail);
                                                         (yyval.obj) = the_object;
                                                     }
#line 1538 "dyngeom_yacc.c" /* yacc.c:1646  */
    break;

  case 20:
#line 241 "dyngeom_parse.y" /* yacc.c:1646  */
    { object_list_singleton(& (yyval.obj_list), (yyvsp[0].obj)); }
#line 1544 "dyngeom_yacc.c" /* yacc.c:1646  */
    break;

  case 21:
#line 242 "dyngeom_parse.y" /* yacc.c:1646  */
    { (yyval.obj_list) = (yyvsp[-1].obj_list); add_object_to_list(& (yyval.obj_list), (yyvsp[0].obj)); }
#line 1550 "dyngeom_yacc.c" /* yacc.c:1646  */
    break;

  case 24:
#line 251 "dyngeom_parse.y" /* yacc.c:1646  */
    { dg_deep_copy_object(dg_parse, (struct object *) (yyvsp[-3].sym)->value, (struct object *) (yyvsp[-1].sym)->value); }
#line 1556 "dyngeom_yacc.c" /* yacc.c:1646  */
    break;

  case 25:
#line 253 "dyngeom_parse.y" /* yacc.c:1646  */
    { (yyval.obj) = (struct object *) (yyvsp[-6].sym)->value; }
#line 1562 "dyngeom_yacc.c" /* yacc.c:1646  */
    break;

  case 26:
#line 256 "dyngeom_parse.y" /* yacc.c:1646  */
    { no_printf("new_object_name\n"); }
#line 1568 "dyngeom_yacc.c" /* yacc.c:1646  */
    break;

  case 27:
#line 263 "dyngeom_parse.y" /* yacc.c:1646  */
    { no_printf("INSTANTIATE\n"); dg_parse->current_object = dg_parse->root_instance; }
#line 1574 "dyngeom_yacc.c" /* yacc.c:1646  */
    break;

  case 28:
#line 264 "dyngeom_parse.y" /* yacc.c:1646  */
    {
                                                        no_printf("meta_object_def\n");
                                                        add_child_objects(dg_parse->root_instance, (yyvsp[0].obj), (yyvsp[0].obj));
                                                        dg_parse->current_object = dg_parse->root_object;
                                                     }
#line 1584 "dyngeom_yacc.c" /* yacc.c:1646  */
    break;

  case 29:
#line 274 "dyngeom_parse.y" /* yacc.c:1646  */
    { add_child_objects(dg_parse->root_object, (yyvsp[0].obj), (yyvsp[0].obj)); no_printf("physical_object_def\n"); }
#line 1590 "dyngeom_yacc.c" /* yacc.c:1646  */
    break;

  case 32:
#line 281 "dyngeom_parse.y" /* yacc.c:1646  */
    { no_printf("new_object %s\n", (yyvsp[0].str)); (yyval.sym) = dg_start_object(dg_parse, (yyvsp[0].str)); }
#line 1596 "dyngeom_yacc.c" /* yacc.c:1646  */
    break;

  case 33:
#line 284 "dyngeom_parse.y" /* yacc.c:1646  */
    { no_printf("start_object\n"); }
#line 1602 "dyngeom_yacc.c" /* yacc.c:1646  */
    break;

  case 34:
#line 287 "dyngeom_parse.y" /* yacc.c:1646  */
    { no_printf("end_object\n"); dg_finish_object(dg_parse); }
#line 1608 "dyngeom_yacc.c" /* yacc.c:1646  */
    break;

  case 38:
#line 298 "dyngeom_parse.y" /* yacc.c:1646  */
    { no_printf("TRANSLATE\n"); }
#line 1614 "dyngeom_yacc.c" /* yacc.c:1646  */
    break;

  case 39:
#line 299 "dyngeom_parse.y" /* yacc.c:1646  */
    { no_printf("SCALE\n"); }
#line 1620 "dyngeom_yacc.c" /* yacc.c:1646  */
    break;

  case 40:
#line 300 "dyngeom_parse.y" /* yacc.c:1646  */
    { no_printf("ROTATE\n"); }
#line 1626 "dyngeom_yacc.c" /* yacc.c:1646  */
    break;

  case 41:
#line 305 "dyngeom_parse.y" /* yacc.c:1646  */
    { no_printf("POLYGON_LIST\n"); }
#line 1632 "dyngeom_yacc.c" /* yacc.c:1646  */
    break;

  case 42:
#line 308 "dyngeom_parse.y" /* yacc.c:1646  */
    { (yyval.obj) = dg_new_polygon_list(dg_parse, (yyvsp[-5].str)); }
#line 1638 "dyngeom_yacc.c" /* yacc.c:1646  */
    break;

  case 43:
#line 312 "dyngeom_parse.y" /* yacc.c:1646  */
    { (yyval.obj) = (struct object *) (yyvsp[-3].obj); dg_finish_object(dg_parse); }
#line 1644 "dyngeom_yacc.c" /* yacc.c:1646  */
    break;

  case 44:
#line 316 "dyngeom_parse.y" /* yacc.c:1646  */
    { no_printf("vertex_list_command\n"); }
#line 1650 "dyngeom_yacc.c" /* yacc.c:1646  */
    break;

  case 46:
#line 320 "dyngeom_parse.y" /* yacc.c:1646  */
    { /*no_printf("single_vertex\n");*/ }
#line 1656 "dyngeom_yacc.c" /* yacc.c:1646  */
    break;

  case 47:
#line 323 "dyngeom_parse.y" /* yacc.c:1646  */
    { no_printf("list_points\n"); }
#line 1662 "dyngeom_yacc.c" /* yacc.c:1646  */
    break;

  case 49:
#line 328 "dyngeom_parse.y" /* yacc.c:1646  */
    { no_printf("element_connection_cmd\n"); }
#line 1668 "dyngeom_yacc.c" /* yacc.c:1646  */
    break;

  case 51:
#line 333 "dyngeom_parse.y" /* yacc.c:1646  */
    {}
#line 1674 "dyngeom_yacc.c" /* yacc.c:1646  */
    break;

  case 52:
#line 335 "dyngeom_parse.y" /* yacc.c:1646  */
    {}
#line 1680 "dyngeom_yacc.c" /* yacc.c:1646  */
    break;

  case 53:
#line 338 "dyngeom_parse.y" /* yacc.c:1646  */
    {}
#line 1686 "dyngeom_yacc.c" /* yacc.c:1646  */
    break;

  case 58:
#line 354 "dyngeom_parse.y" /* yacc.c:1646  */
    { no_printf("element_specifier\n"); }
#line 1692 "dyngeom_yacc.c" /* yacc.c:1646  */
    break;

  case 60:
#line 363 "dyngeom_parse.y" /* yacc.c:1646  */
    { no_printf("incl_element_list_stmt\n"); }
#line 1698 "dyngeom_yacc.c" /* yacc.c:1646  */
    break;

  case 63:
#line 371 "dyngeom_parse.y" /* yacc.c:1646  */
    { no_printf("element_spec\n"); }
#line 1704 "dyngeom_yacc.c" /* yacc.c:1646  */
    break;

  case 68:
#line 389 "dyngeom_parse.y" /* yacc.c:1646  */
    { dg_parse->current_region = (yyvsp[-1].reg); }
#line 1710 "dyngeom_yacc.c" /* yacc.c:1646  */
    break;

  case 69:
#line 391 "dyngeom_parse.y" /* yacc.c:1646  */
    { dg_parse->current_region = NULL; }
#line 1716 "dyngeom_yacc.c" /* yacc.c:1646  */
    break;

  case 70:
#line 395 "dyngeom_parse.y" /* yacc.c:1646  */
    { dg_create_region(dg_parse->reg_sym_table, dg_parse->current_object, (yyvsp[0].str)); no_printf("new_region\n"); }
#line 1722 "dyngeom_yacc.c" /* yacc.c:1646  */
    break;

  case 72:
#line 406 "dyngeom_parse.y" /* yacc.c:1646  */
    { (yyval.tok) = X_PARTS; }
#line 1728 "dyngeom_yacc.c" /* yacc.c:1646  */
    break;

  case 73:
#line 407 "dyngeom_parse.y" /* yacc.c:1646  */
    { (yyval.tok) = Y_PARTS; }
#line 1734 "dyngeom_yacc.c" /* yacc.c:1646  */
    break;

  case 74:
#line 408 "dyngeom_parse.y" /* yacc.c:1646  */
    { (yyval.tok) = Z_PARTS; }
#line 1740 "dyngeom_yacc.c" /* yacc.c:1646  */
    break;

  case 75:
#line 413 "dyngeom_parse.y" /* yacc.c:1646  */
    {
                                                          no_printf("include_stmt %s\n", (yyvsp[0].str));
                                                          char *include_path = find_include_file((yyvsp[0].str), dg_parse->curr_file);
                                                          if (include_path == NULL)
                                                          {
                                                            free((yyvsp[0].str));
                                                            return 1;
                                                          }
                                                          if (parse_dg(dg_parse, include_path))
                                                          {
                                                            free(include_path);
                                                            free((yyvsp[0].str));
                                                            return 1;
                                                          }
                                                          free(include_path);
                                                          free((yyvsp[0].str));
                                                      }
#line 1762 "dyngeom_yacc.c" /* yacc.c:1646  */
    break;

  case 76:
#line 435 "dyngeom_parse.y" /* yacc.c:1646  */
    { }
#line 1768 "dyngeom_yacc.c" /* yacc.c:1646  */
    break;

  case 77:
#line 438 "dyngeom_parse.y" /* yacc.c:1646  */
    { (yyval.nlist) = (yyvsp[-1].nlist); }
#line 1774 "dyngeom_yacc.c" /* yacc.c:1646  */
    break;

  case 78:
#line 441 "dyngeom_parse.y" /* yacc.c:1646  */
    { }
#line 1780 "dyngeom_yacc.c" /* yacc.c:1646  */
    break;

  case 80:
#line 445 "dyngeom_parse.y" /* yacc.c:1646  */
    { }
#line 1786 "dyngeom_yacc.c" /* yacc.c:1646  */
    break;

  case 81:
#line 446 "dyngeom_parse.y" /* yacc.c:1646  */
    { (yyval.dbl) = *(double *) (yyvsp[0].sym)->value; }
#line 1792 "dyngeom_yacc.c" /* yacc.c:1646  */
    break;

  case 82:
#line 449 "dyngeom_parse.y" /* yacc.c:1646  */
    { (yyval.dbl) = (yyvsp[0].llival); }
#line 1798 "dyngeom_yacc.c" /* yacc.c:1646  */
    break;

  case 83:
#line 450 "dyngeom_parse.y" /* yacc.c:1646  */
    { }
#line 1804 "dyngeom_yacc.c" /* yacc.c:1646  */
    break;

  case 86:
#line 458 "dyngeom_parse.y" /* yacc.c:1646  */
    { no_printf("existing_num_var\n"); }
#line 1810 "dyngeom_yacc.c" /* yacc.c:1646  */
    break;

  case 87:
#line 462 "dyngeom_parse.y" /* yacc.c:1646  */
    { (yyval.dbl) = (yyvsp[-1].dbl); }
#line 1816 "dyngeom_yacc.c" /* yacc.c:1646  */
    break;

  case 88:
#line 463 "dyngeom_parse.y" /* yacc.c:1646  */
    { (yyval.dbl) = -(yyvsp[0].dbl); }
#line 1822 "dyngeom_yacc.c" /* yacc.c:1646  */
    break;

  case 89:
#line 467 "dyngeom_parse.y" /* yacc.c:1646  */
    { no_printf("str_expr_only\n"); }
#line 1828 "dyngeom_yacc.c" /* yacc.c:1646  */
    break;

  case 90:
#line 471 "dyngeom_parse.y" /* yacc.c:1646  */
    { no_printf("str_value %s\n", (yyval.str)); (yyval.str) = strip_quotes((yyvsp[0].str)); }
#line 1834 "dyngeom_yacc.c" /* yacc.c:1646  */
    break;

  case 91:
#line 472 "dyngeom_parse.y" /* yacc.c:1646  */
    { (yyval.str) = my_strcat((yyvsp[-2].str), (yyvsp[0].str)); }
#line 1840 "dyngeom_yacc.c" /* yacc.c:1646  */
    break;


#line 1844 "dyngeom_yacc.c" /* yacc.c:1646  */
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
      yyerror (dg_parse, scanner, YY_("syntax error"));
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
        yyerror (dg_parse, scanner, yymsgp);
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
                      yytoken, &yylval, dg_parse, scanner);
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
                  yystos[yystate], yyvsp, dg_parse, scanner);
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
  yyerror (dg_parse, scanner, YY_("memory exhausted"));
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
                  yytoken, &yylval, dg_parse, scanner);
    }
  /* Do not reclaim the symbols of the rule whose action triggered
     this YYABORT or YYACCEPT.  */
  YYPOPSTACK (yylen);
  YY_STACK_PRINT (yyss, yyssp);
  while (yyssp != yyss)
    {
      yydestruct ("Cleanup: popping",
                  yystos[*yyssp], yyvsp, dg_parse, scanner);
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
#line 475 "dyngeom_parse.y" /* yacc.c:1906  */


void dgerror(
    struct dyngeom_parse_vars *dg_parse,
    yyscan_t scanner,
    char const *str) {
  mcell_error("%s on line %d in %s\n", str, dg_parse->line_num[dg_parse->include_stack_ptr - 1], dg_parse->curr_file);
}

/*int main(int argc, char *argv[])*/
/*{*/
/*  struct dyngeom_parse_vars *dg_parse = create_dg_parse();*/
/*  parse_dg_init(dg_parse, argv[1]);*/
/*}*/
