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
#line 99 "dyngeom_parse.y" /* yacc.c:1909  */

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

#line 113 "dyngeom_yacc.h" /* yacc.c:1909  */
};

typedef union YYSTYPE YYSTYPE;
# define YYSTYPE_IS_TRIVIAL 1
# define YYSTYPE_IS_DECLARED 1
#endif



int dgparse (struct dyngeom_parse_vars *dg_parse, yyscan_t scanner);

#endif /* !YY_DG_DYNGEOM_YACC_H_INCLUDED  */
