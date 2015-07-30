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

  #include <stdio.h>
  #include <stdlib.h>
  #include <string.h>
  #include <math.h>

  #include "sym_table.h"
  #include "mcell_objects.h"
  #include "mcell_structs.h"
  #include "mcell_misc.h"
  #include "mcell_objects.h"

  #define EPS_C 1e-12

  typedef void *yyscan_t;
  void yyerror(char *);
  int yylex(void);
  extern FILE *yyin;

  struct dyngeom_parse_vars *dg_parse;

  struct num_expr_list_head_dg {
    struct num_expr_list *value_head;
    struct num_expr_list *value_tail;
    int value_count;
    int shared;
  };

  struct object_list {
    struct object *obj_head;
    struct object *obj_tail;
  };

  struct dyngeom_parse_vars {
    struct sym_table_head *obj_sym_table;
    struct object *root_object;
    struct object *root_instance;
    struct object *current_object;
    struct name_list *object_name_list;
    struct name_list *object_name_list_end;
  };

  int init_top_level_objs(struct dyngeom_parse_vars *dg_parse);
  void object_list_singleton(struct object_list *head, struct object *objp);
  void add_object_to_list(struct object_list *head, struct object *objp);
  struct vector3 *point_scalar(double val);
  int advance_range_dg(struct num_expr_list_head_dg *list, double tmp_dbl);
  void mcell_free_numeric_list_dg(struct num_expr_list *nel);
  int generate_range(
      struct num_expr_list_head_dg *list,
      double start,
      double end,
      double step);
  struct sym_table *dg_start_object(
      struct dyngeom_parse_vars *dg_parse, char *name);
  void dg_finish_object(struct dyngeom_parse_vars *dg_parse);

  int generate_range(
      struct num_expr_list_head_dg *list,
      double start,
      double end,
      double step) {
    list->value_head = NULL;
    list->value_tail = NULL;
    list->value_count = 0;
    list->shared = 0;

    if (step > 0) {
      for (double tmp_dbl = start;
           tmp_dbl < end || !distinguishable(tmp_dbl, end, EPS_C) ||
               fabs(end - tmp_dbl) <= EPS_C;
           tmp_dbl += step) {
        if (advance_range_dg(list, tmp_dbl))
          return 1;
      }
    } else /* if (step < 0) */
    {
      for (double tmp_dbl = start;
           tmp_dbl > end || !distinguishable(tmp_dbl, end, EPS_C) ||
               fabs(end - tmp_dbl) <= EPS_C;
           tmp_dbl += step) {
        if (advance_range_dg(list, tmp_dbl))
          return 1;
      }
    }
    return 0;
  }

  // This is the same as advance_range in mcell_misc.h, but including that header
  // here causes a number of build problems that are currently difficult to
  // resolve.
  int advance_range_dg(struct num_expr_list_head_dg *list, double tmp_dbl) {
    struct num_expr_list *nel;
    nel = (struct num_expr_list *)malloc(sizeof(struct num_expr_list));
    if (nel == NULL) {
      mcell_free_numeric_list_dg(list->value_head);
      list->value_head = list->value_tail = NULL;
      return 1;
    }
    nel->value = tmp_dbl;
    nel->next = NULL;

    ++list->value_count;
    if (list->value_tail != NULL)
      list->value_tail->next = nel;
    else
      list->value_head = nel;
    list->value_tail = nel;
    return 0;
  }

  void mcell_free_numeric_list_dg(struct num_expr_list *nel) {
    while (nel != NULL) {
      struct num_expr_list *n = nel;
      nel = nel->next;
      free(n);
    }
  }

  struct vector3 *point_scalar(double val) {
    struct vector3 *vec;
    vec = (struct vector3 *)malloc(sizeof(struct vector3));
    if (!vec)
      return NULL;

    vec->x = val;
    vec->y = val;
    vec->z = val;
    return vec;
  }

  void object_list_singleton(struct object_list *head, struct object *objp) {
    objp->next = NULL;
    head->obj_tail = head->obj_head = objp;
  }

  void add_object_to_list(struct object_list *head, struct object *objp) {
    objp->next = NULL;
    head->obj_tail = head->obj_tail->next = objp;
  }

  int init_top_level_objs(struct dyngeom_parse_vars *dg_parse_vars) {
    if ((dg_parse_vars->obj_sym_table = init_symtab(1024)) == NULL) {
      return 1;
    }

    struct sym_table *sym;
    if ((sym = store_sym(
        "WORLD_OBJ", OBJ, dg_parse_vars->obj_sym_table, NULL)) == NULL) {
      return 1;
    }

    dg_parse_vars->root_object = (struct object *)sym->value;
    dg_parse_vars->root_object->object_type = META_OBJ;
    if (!(dg_parse_vars->root_object->last_name = CHECKED_STRDUP_NODIE("", NULL))) {
      return 1;
    }

    if ((sym = store_sym(
        "WORLD_INSTANCE", OBJ, dg_parse_vars->obj_sym_table, NULL)) == NULL) {
      return 1;
    }

    dg_parse_vars->root_instance = (struct object *)sym->value;
    dg_parse_vars->root_instance->object_type = META_OBJ;
    if (!(dg_parse_vars->root_instance->last_name = CHECKED_STRDUP("", NULL))) {
      return 1;
    }

    dg_parse_vars->current_object = dg_parse_vars->root_instance;

    return 0;
  }

  struct sym_table *dg_start_object(
      struct dyngeom_parse_vars *dg_parse_vars,
      char *name) {
    // Create new fully qualified name.
    char *new_name;
    struct object_creation obj_creation;
    obj_creation.object_name_list = dg_parse_vars->object_name_list;
    obj_creation.object_name_list_end = dg_parse_vars->object_name_list_end;
    if ((new_name = push_object_name(&obj_creation, name)) == NULL) {
      free(name);
      return NULL;
    }
    dg_parse_vars->object_name_list = obj_creation.object_name_list;
    dg_parse_vars->object_name_list_end = obj_creation.object_name_list_end;

    // Create the symbol, if it doesn't exist yet.
    struct object *obj_ptr = make_new_object(
        dg_parse_vars->obj_sym_table, new_name, 0);
    if (obj_ptr == NULL) {
      if (name != new_name) {
        free(name);
      }
      free(new_name);
      return NULL;
    }

    struct sym_table *sym_ptr = obj_ptr->sym;
    obj_ptr->last_name = name;

    // Set parent object, make this object "current".
    obj_ptr->parent = dg_parse_vars->current_object;
    dg_parse_vars->current_object = obj_ptr;

    return sym_ptr;
  }

  void dg_finish_object(struct dyngeom_parse_vars *dg_parse_vars) {
    struct object_creation obj_creation;
    obj_creation.object_name_list_end = dg_parse_vars->object_name_list_end;

    pop_object_name(&obj_creation);
    dg_parse_vars->object_name_list_end = obj_creation.object_name_list_end;
    dg_parse_vars->current_object = dg_parse_vars->current_object->parent;
  }

  int parse_dg(char *dynamic_geometry_filename) {
    FILE *fp=fopen(dynamic_geometry_filename,"r");
    if(!fp)
    {
      printf("Couldn't open file for reading\n");
      return 1;
    }
    yyin=fp;
    dg_parse = (struct dyngeom_parse_vars *)malloc(sizeof(struct dyngeom_parse_vars));
    memset(dg_parse, 0, sizeof(struct dyngeom_parse_vars));
    init_top_level_objs(dg_parse);
    yyparse();
    fclose(fp);
    return 0;
  }


#line 303 "dyngeom_yacc.c" /* yacc.c:339  */

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
#line 238 "dyngeom_parse.y" /* yacc.c:355  */

  int tok;
  double dbl;
  char *str;
  long long llival;
  struct sym_table *sym;
  struct vector3 *vec3;
  struct num_expr_list_head_dg nlist;
  struct object *obj;
  struct object_list obj_list;

#line 397 "dyngeom_yacc.c" /* yacc.c:355  */
};
# define YYSTYPE_IS_TRIVIAL 1
# define YYSTYPE_IS_DECLARED 1
#endif


extern YYSTYPE yylval;

int yyparse (void);

#endif /* !YY_YY_DYNGEOM_YACC_H_INCLUDED  */

/* Copy the second part of user declarations.  */

#line 412 "dyngeom_yacc.c" /* yacc.c:358  */

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
#define YYLAST   119

/* YYNTOKENS -- Number of terminals.  */
#define YYNTOKENS  37
/* YYNNTS -- Number of nonterminals.  */
#define YYNNTS  57
/* YYNRULES -- Number of rules.  */
#define YYNRULES  83
/* YYNSTATES -- Number of states.  */
#define YYNSTATES  141

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
       0,   305,   305,   309,   310,   314,   315,   316,   321,   324,
     328,   331,   332,   336,   337,   350,   351,   356,   370,   371,
     374,   375,   379,   385,   392,   392,   404,   407,   408,   411,
     414,   417,   420,   421,   424,   428,   429,   430,   435,   438,
     435,   451,   451,   455,   458,   459,   463,   463,   468,   469,
     473,   476,   477,   483,   487,   488,   493,   497,   502,   503,
     506,   507,   512,   518,   519,   524,   529,   536,   540,   541,
     542,   548,   551,   554,   555,   558,   559,   562,   563,   566,
     567,   571,   575,   576
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
  "list_objects", "object_ref", "existing_object_ref", "new_object_name",
  "instance_def", "$@1", "physical_object_def", "object_def", "new_object",
  "start_object", "end_object", "list_opt_object_cmds", "opt_object_cmd",
  "transformation", "polygon_list_def", "$@2", "$@3", "vertex_list_cmd",
  "$@4", "single_vertex", "list_points", "element_connection_cmd", "$@5",
  "list_element_connections", "element_connection",
  "list_opt_polygon_object_cmds", "opt_polygon_object_cmd",
  "element_specifier_list", "element_specifier", "incl_element_list_stmt",
  "list_element_specs", "element_spec", "in_obj_define_surface_regions",
  "list_in_obj_surface_region_defs", "in_obj_surface_region_def",
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

#define YYPACT_NINF -68

#define yypact_value_is_default(Yystate) \
  (!!((Yystate) == (-68)))

#define YYTABLE_NINF -30

#define yytable_value_is_error(Yytable_value) \
  0

  /* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
     STATE-NUM.  */
static const yytype_int8 yypact[] =
{
      64,   -68,   -68,   -68,   -68,   -68,     6,    64,   -68,    14,
     -68,    27,   -68,   -68,   -68,    38,   -68,   -68,   -10,    42,
     -68,   -68,   -68,     9,    32,   -68,   -68,     9,   -68,    42,
      -3,   -68,   -68,    53,    42,   -68,   -68,   -68,    61,   -68,
     -68,     1,     1,     1,   -68,    29,   -68,   -68,   -68,   -68,
     -68,   -68,   -68,    71,   -68,    28,     7,   -68,    66,    45,
      -3,   -68,    49,   -68,   -68,    62,    63,    65,   -68,   -68,
     -68,   -68,   -68,     9,     1,   -68,   -68,    32,    52,   -68,
      32,     8,    32,   -68,    72,   -68,   -68,    -7,   -68,    32,
      68,   -68,   -68,   -68,   -68,   -68,   -68,    59,    28,     1,
     -68,   -68,    23,   -68,   -68,    57,    31,   -68,   -68,     1,
     -68,    60,   -68,   -68,    42,   -68,   -68,   -68,   -68,     4,
     -68,    58,   -68,   -68,    75,    74,    18,   -68,   -68,    67,
     -68,   -68,     1,    46,   -68,    82,     1,   -68,     1,   -68,
     -68
};

  /* YYDEFACT[STATE-NUM] -- Default reduction number in state STATE-NUM.
     Performed when YYTABLE does not specify something else to do.  Zero
     means the default is an error.  */
static const yytype_uint8 yydefact[] =
{
       0,    24,    68,    69,    70,     8,     0,     2,     3,    23,
      27,     0,     7,     6,    26,     0,    28,     5,     0,     0,
       1,     4,    38,     0,     0,    29,    25,     0,    30,     0,
       0,    67,    71,     0,    32,    18,    20,    21,     0,    78,
      77,     0,     0,     0,    81,     0,    13,    15,    73,    75,
      76,    74,    41,     0,    19,     0,     0,    83,     0,     0,
       0,    72,     0,    46,    39,     0,     0,     0,    31,    17,
      33,    34,     9,     0,     0,    82,    14,     0,     0,    51,
       0,     0,     0,    32,     0,    43,    44,     0,    10,     0,
      32,    35,    11,    36,    79,    12,    80,     0,     0,     0,
      42,    45,     0,    48,    50,     0,     0,    52,    53,     0,
      22,     0,    47,    49,     0,    40,    37,    16,    66,     0,
      63,     0,    62,    64,     0,     0,     0,    54,    56,     0,
      65,    55,     0,     0,    58,    60,     0,    57,     0,    59,
      61
};

  /* YYPGOTO[NTERM-NUM].  */
static const yytype_int8 yypgoto[] =
{
     -68,   -68,   -68,    93,     0,   -68,   -13,   -68,   -68,    41,
      83,   -68,    69,   -68,   -68,   -68,   -68,   -68,    -4,    21,
     -22,    10,   -67,   -68,   -68,   -68,   -68,   -68,   -68,   -68,
      17,   -68,   -68,   -68,   -68,     3,   -68,   -68,   -68,   -20,
     -68,   -68,   -29,   -68,   -68,    -9,   -68,   -68,   -68,   -23,
     -68,   -39,   -68,    30,   -68,   -68,    34
};

  /* YYDEFGOTO[NTERM-NUM].  */
static const yytype_int16 yydefgoto[] =
{
      -1,     6,     7,     8,    44,    73,    85,    93,    45,    46,
      10,    34,    35,    36,    11,    12,    19,    13,    14,    15,
      29,    69,    55,    70,    71,    16,    27,    79,    53,    62,
      86,    87,    64,    78,   102,   103,    90,   107,   126,   127,
     128,   133,   134,   108,   119,   120,   121,    17,    18,    88,
      32,    47,    48,    49,    95,    50,    51
};

  /* YYTABLE[YYPACT[STATE-NUM]] -- What to do in state STATE-NUM.  If
     positive, shift that token.  If negative, reduce the rule whose
     number is the opposite.  If YYTABLE_NINF, syntax error.  */
static const yytype_int16 yytable[] =
{
       9,    31,    57,    58,    59,    33,    20,     9,     5,    39,
      40,    24,     5,    39,    40,     5,    98,   -29,     5,    25,
      39,    40,    41,   106,    30,    37,    41,   100,    42,     9,
      37,    22,    43,    41,     9,    84,    43,   125,   122,    30,
      28,    23,    28,    43,    65,    66,    67,    65,    66,    67,
      38,    83,   130,     5,    30,    38,    72,   112,    52,    60,
     111,    61,    68,    30,    56,   115,   104,    91,    92,    97,
     116,     1,     2,     3,     4,     5,   136,    63,   137,   104,
      74,    75,    77,    80,    81,    89,    82,    99,   105,   109,
     114,   124,   117,   135,   125,   129,   138,   135,   132,   140,
      21,    76,    26,    54,   101,   113,   131,   139,   110,     0,
     123,    94,     0,     0,   118,    96,     0,     0,     0,   118
};

static const yytype_int16 yycheck[] =
{
       0,    24,    41,    42,    43,    27,     0,     7,    11,    12,
      13,    21,    11,    12,    13,    11,    83,     3,    11,    19,
      12,    13,    25,    90,    31,    29,    25,    34,    31,    29,
      34,     4,    35,    25,    34,    74,    35,    19,    34,    31,
      33,     3,    33,    35,    16,    17,    18,    16,    17,    18,
      29,    73,    34,    11,    31,    34,    56,    34,     5,    30,
      99,    32,    34,    31,     3,    34,    89,    80,    81,    82,
     109,     7,     8,     9,    10,    11,    30,     6,    32,   102,
      14,    36,    33,    21,    21,    33,    21,    15,    20,    30,
      33,    33,    32,   132,    19,    21,    14,   136,    31,   138,
       7,    60,    19,    34,    87,   102,   126,   136,    98,    -1,
     119,    81,    -1,    -1,   114,    81,    -1,    -1,    -1,   119
};

  /* YYSTOS[STATE-NUM] -- The (internal number of the) accessing
     symbol of state STATE-NUM.  */
static const yytype_uint8 yystos[] =
{
       0,     7,     8,     9,    10,    11,    38,    39,    40,    41,
      47,    51,    52,    54,    55,    56,    62,    84,    85,    53,
       0,    40,     4,     3,    21,    41,    47,    63,    33,    57,
      31,    86,    87,    57,    48,    49,    50,    55,    56,    12,
      13,    25,    31,    35,    41,    45,    46,    88,    89,    90,
      92,    93,     5,    65,    49,    59,     3,    88,    88,    88,
      30,    32,    66,     6,    69,    16,    17,    18,    34,    58,
      60,    61,    41,    42,    14,    36,    46,    33,    70,    64,
      21,    21,    21,    57,    88,    43,    67,    68,    86,    33,
      73,    43,    43,    44,    90,    91,    93,    43,    59,    15,
      34,    67,    71,    72,    86,    20,    59,    74,    80,    30,
      58,    88,    34,    72,    33,    34,    88,    32,    41,    81,
      82,    83,    34,    82,    33,    19,    75,    76,    77,    21,
      34,    76,    31,    78,    79,    88,    30,    32,    14,    79,
      88
};

  /* YYR1[YYN] -- Symbol number of symbol that rule YYN derives.  */
static const yytype_uint8 yyr1[] =
{
       0,    37,    38,    39,    39,    40,    40,    40,    41,    42,
      43,    44,    44,    45,    45,    46,    46,    47,    48,    48,
      49,    49,    50,    51,    53,    52,    54,    55,    55,    56,
      57,    58,    59,    59,    60,    61,    61,    61,    63,    64,
      62,    66,    65,    67,    68,    68,    70,    69,    71,    71,
      72,    73,    73,    74,    75,    75,    76,    77,    78,    78,
      79,    79,    80,    81,    81,    82,    83,    84,    85,    85,
      85,    86,    87,    88,    88,    89,    89,    90,    90,    91,
      91,    92,    93,    93
};

  /* YYR2[YYN] -- Number of symbols on the right hand side of rule YYN.  */
static const yytype_uint8 yyr2[] =
{
       0,     2,     1,     1,     2,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     3,     1,     7,     6,     1,     2,
       1,     1,     6,     1,     0,     3,     1,     1,     1,     1,
       1,     1,     0,     2,     1,     3,     3,     5,     0,     0,
      10,     0,     5,     1,     1,     2,     0,     5,     1,     2,
       1,     0,     2,     1,     1,     2,     1,     5,     1,     3,
       1,     3,     4,     1,     2,     4,     1,     3,     1,     1,
       1,     1,     3,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     3,     2
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
#line 324 "dyngeom_parse.y" /* yacc.c:1646  */
    {printf("existing_object\n");}
#line 1599 "dyngeom_yacc.c" /* yacc.c:1646  */
    break;

  case 10:
#line 328 "dyngeom_parse.y" /* yacc.c:1646  */
    {printf("point\n");}
#line 1605 "dyngeom_yacc.c" /* yacc.c:1646  */
    break;

  case 11:
#line 331 "dyngeom_parse.y" /* yacc.c:1646  */
    {printf("point_or_num\n");}
#line 1611 "dyngeom_yacc.c" /* yacc.c:1646  */
    break;

  case 12:
#line 332 "dyngeom_parse.y" /* yacc.c:1646  */
    { (yyval.vec3) = point_scalar((yyvsp[0].dbl)); }
#line 1617 "dyngeom_yacc.c" /* yacc.c:1646  */
    break;

  case 14:
#line 337 "dyngeom_parse.y" /* yacc.c:1646  */
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
#line 1633 "dyngeom_yacc.c" /* yacc.c:1646  */
    break;

  case 15:
#line 350 "dyngeom_parse.y" /* yacc.c:1646  */
    { printf("range_spec\n"); }
#line 1639 "dyngeom_yacc.c" /* yacc.c:1646  */
    break;

  case 16:
#line 351 "dyngeom_parse.y" /* yacc.c:1646  */
    { generate_range(&(yyval.nlist), (yyvsp[-5].dbl), (yyvsp[-3].dbl), (yyvsp[-1].dbl)); }
#line 1645 "dyngeom_yacc.c" /* yacc.c:1646  */
    break;

  case 17:
#line 360 "dyngeom_parse.y" /* yacc.c:1646  */
    {
                                                         struct object *the_object = (struct object *) (yyvsp[-5].sym)->value;
                                                         the_object->object_type = META_OBJ;
                                                         add_child_objects(the_object, (yyvsp[-2].obj_list).obj_head, (yyvsp[-2].obj_list).obj_tail);
                                                         (yyval.obj) = the_object;
                                                     }
#line 1656 "dyngeom_yacc.c" /* yacc.c:1646  */
    break;

  case 18:
#line 370 "dyngeom_parse.y" /* yacc.c:1646  */
    { object_list_singleton(& (yyval.obj_list), (yyvsp[0].obj)); }
#line 1662 "dyngeom_yacc.c" /* yacc.c:1646  */
    break;

  case 19:
#line 371 "dyngeom_parse.y" /* yacc.c:1646  */
    { (yyval.obj_list) = (yyvsp[-1].obj_list); add_object_to_list(& (yyval.obj_list), (yyvsp[0].obj)); }
#line 1668 "dyngeom_yacc.c" /* yacc.c:1646  */
    break;

  case 22:
#line 382 "dyngeom_parse.y" /* yacc.c:1646  */
    { (yyval.obj) = (struct object *) (yyvsp[-5].sym)->value; }
#line 1674 "dyngeom_yacc.c" /* yacc.c:1646  */
    break;

  case 23:
#line 385 "dyngeom_parse.y" /* yacc.c:1646  */
    {printf("new_object_name\n");}
#line 1680 "dyngeom_yacc.c" /* yacc.c:1646  */
    break;

  case 24:
#line 392 "dyngeom_parse.y" /* yacc.c:1646  */
    { printf("INSTANTIATE\n"); dg_parse->current_object = dg_parse->root_instance; }
#line 1686 "dyngeom_yacc.c" /* yacc.c:1646  */
    break;

  case 25:
#line 393 "dyngeom_parse.y" /* yacc.c:1646  */
    {
                                                        printf("meta_object_def\n");
                                                        add_child_objects(dg_parse->root_instance, (yyvsp[0].obj), (yyvsp[0].obj));
                                                        dg_parse->current_object = dg_parse->root_object;
                                                      }
#line 1696 "dyngeom_yacc.c" /* yacc.c:1646  */
    break;

  case 26:
#line 404 "dyngeom_parse.y" /* yacc.c:1646  */
    {printf("physical_object_def\n");}
#line 1702 "dyngeom_yacc.c" /* yacc.c:1646  */
    break;

  case 29:
#line 411 "dyngeom_parse.y" /* yacc.c:1646  */
    {printf("new_object\n"); (yyval.sym) = dg_start_object(dg_parse, (yyvsp[0].str));}
#line 1708 "dyngeom_yacc.c" /* yacc.c:1646  */
    break;

  case 30:
#line 414 "dyngeom_parse.y" /* yacc.c:1646  */
    {printf("start_object\n");}
#line 1714 "dyngeom_yacc.c" /* yacc.c:1646  */
    break;

  case 31:
#line 417 "dyngeom_parse.y" /* yacc.c:1646  */
    {printf("end_object\n"); dg_finish_object(dg_parse);}
#line 1720 "dyngeom_yacc.c" /* yacc.c:1646  */
    break;

  case 35:
#line 428 "dyngeom_parse.y" /* yacc.c:1646  */
    {printf("TRANSLATE\n");}
#line 1726 "dyngeom_yacc.c" /* yacc.c:1646  */
    break;

  case 36:
#line 429 "dyngeom_parse.y" /* yacc.c:1646  */
    {printf("SCALE\n");}
#line 1732 "dyngeom_yacc.c" /* yacc.c:1646  */
    break;

  case 37:
#line 430 "dyngeom_parse.y" /* yacc.c:1646  */
    {printf("ROTATE\n");}
#line 1738 "dyngeom_yacc.c" /* yacc.c:1646  */
    break;

  case 38:
#line 435 "dyngeom_parse.y" /* yacc.c:1646  */
    {printf("POLYGON_LIST");}
#line 1744 "dyngeom_yacc.c" /* yacc.c:1646  */
    break;

  case 39:
#line 438 "dyngeom_parse.y" /* yacc.c:1646  */
    {
                                                        //XXX: Need to actually create objects
                                                        /*$<obj>$ = mdl_new_polygon_list(*/
                                                        /*  parse_state, $1, $4.vertex_count, $4.vertex_head,*/
                                                        /*  $5.connection_count, $5.connection_head);*/
                                                      }
#line 1755 "dyngeom_yacc.c" /* yacc.c:1646  */
    break;

  case 40:
#line 447 "dyngeom_parse.y" /* yacc.c:1646  */
    { (yyval.obj) = (struct object *) (yyvsp[-4].obj); }
#line 1761 "dyngeom_yacc.c" /* yacc.c:1646  */
    break;

  case 41:
#line 451 "dyngeom_parse.y" /* yacc.c:1646  */
    {printf("vertex_list_command\n");}
#line 1767 "dyngeom_yacc.c" /* yacc.c:1646  */
    break;

  case 43:
#line 455 "dyngeom_parse.y" /* yacc.c:1646  */
    {printf("single_vertex\n");}
#line 1773 "dyngeom_yacc.c" /* yacc.c:1646  */
    break;

  case 44:
#line 458 "dyngeom_parse.y" /* yacc.c:1646  */
    {printf("list_points\n");}
#line 1779 "dyngeom_yacc.c" /* yacc.c:1646  */
    break;

  case 46:
#line 463 "dyngeom_parse.y" /* yacc.c:1646  */
    {printf("element_connection_cmd\n");}
#line 1785 "dyngeom_yacc.c" /* yacc.c:1646  */
    break;

  case 48:
#line 468 "dyngeom_parse.y" /* yacc.c:1646  */
    {printf("element_connection\n");}
#line 1791 "dyngeom_yacc.c" /* yacc.c:1646  */
    break;

  case 49:
#line 470 "dyngeom_parse.y" /* yacc.c:1646  */
    {printf("element_connection\n");}
#line 1797 "dyngeom_yacc.c" /* yacc.c:1646  */
    break;

  case 50:
#line 473 "dyngeom_parse.y" /* yacc.c:1646  */
    {printf("element_connection\n");}
#line 1803 "dyngeom_yacc.c" /* yacc.c:1646  */
    break;

  case 55:
#line 489 "dyngeom_parse.y" /* yacc.c:1646  */
    {printf("element_specifier\n");}
#line 1809 "dyngeom_yacc.c" /* yacc.c:1646  */
    break;

  case 57:
#line 498 "dyngeom_parse.y" /* yacc.c:1646  */
    {printf("incl_element_list_stmt\n");}
#line 1815 "dyngeom_yacc.c" /* yacc.c:1646  */
    break;

  case 60:
#line 506 "dyngeom_parse.y" /* yacc.c:1646  */
    {printf("element_spec\n");}
#line 1821 "dyngeom_yacc.c" /* yacc.c:1646  */
    break;

  case 66:
#line 529 "dyngeom_parse.y" /* yacc.c:1646  */
    {printf("new_region\n");}
#line 1827 "dyngeom_yacc.c" /* yacc.c:1646  */
    break;

  case 68:
#line 540 "dyngeom_parse.y" /* yacc.c:1646  */
    { (yyval.tok) = X_PARTS; }
#line 1833 "dyngeom_yacc.c" /* yacc.c:1646  */
    break;

  case 69:
#line 541 "dyngeom_parse.y" /* yacc.c:1646  */
    { (yyval.tok) = Y_PARTS; }
#line 1839 "dyngeom_yacc.c" /* yacc.c:1646  */
    break;

  case 70:
#line 542 "dyngeom_parse.y" /* yacc.c:1646  */
    { (yyval.tok) = Z_PARTS; }
#line 1845 "dyngeom_yacc.c" /* yacc.c:1646  */
    break;

  case 71:
#line 548 "dyngeom_parse.y" /* yacc.c:1646  */
    {printf("array_value\n");}
#line 1851 "dyngeom_yacc.c" /* yacc.c:1646  */
    break;

  case 72:
#line 551 "dyngeom_parse.y" /* yacc.c:1646  */
    { printf("array_expr_only\n"); (yyval.nlist) = (yyvsp[-1].nlist); }
#line 1857 "dyngeom_yacc.c" /* yacc.c:1646  */
    break;

  case 73:
#line 554 "dyngeom_parse.y" /* yacc.c:1646  */
    { printf("num_expr\n"); }
#line 1863 "dyngeom_yacc.c" /* yacc.c:1646  */
    break;

  case 75:
#line 558 "dyngeom_parse.y" /* yacc.c:1646  */
    { printf("num_value\n"); }
#line 1869 "dyngeom_yacc.c" /* yacc.c:1646  */
    break;

  case 76:
#line 559 "dyngeom_parse.y" /* yacc.c:1646  */
    { printf("num_value\n"); (yyval.dbl) = *(double *) (yyvsp[0].sym)->value; }
#line 1875 "dyngeom_yacc.c" /* yacc.c:1646  */
    break;

  case 77:
#line 562 "dyngeom_parse.y" /* yacc.c:1646  */
    { printf("LLINTEGER\n"); (yyval.dbl) = (yyvsp[0].llival); }
#line 1881 "dyngeom_yacc.c" /* yacc.c:1646  */
    break;

  case 78:
#line 563 "dyngeom_parse.y" /* yacc.c:1646  */
    { printf("REAL\n"); }
#line 1887 "dyngeom_yacc.c" /* yacc.c:1646  */
    break;

  case 81:
#line 571 "dyngeom_parse.y" /* yacc.c:1646  */
    { printf("existing_num_var\n"); }
#line 1893 "dyngeom_yacc.c" /* yacc.c:1646  */
    break;

  case 82:
#line 575 "dyngeom_parse.y" /* yacc.c:1646  */
    { (yyval.dbl) = (yyvsp[-1].dbl); }
#line 1899 "dyngeom_yacc.c" /* yacc.c:1646  */
    break;

  case 83:
#line 576 "dyngeom_parse.y" /* yacc.c:1646  */
    { (yyval.dbl) = -(yyvsp[0].dbl); }
#line 1905 "dyngeom_yacc.c" /* yacc.c:1646  */
    break;


#line 1909 "dyngeom_yacc.c" /* yacc.c:1646  */
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
#line 580 "dyngeom_parse.y" /* yacc.c:1906  */


void yyerror(char *s) {
  printf("%s\n", s);
}

/*int main(int argc, char *argv[])*/
/*{*/
/*  parse_dg();*/
/*}*/
