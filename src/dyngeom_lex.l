%{
  #include "config.h"
  #include <stdio.h>
  #include <stdlib.h>
  #include <string.h>
  #include <unistd.h>

  #include "dyngeom_parse_extras.h"
  #include "dyngeom_yacc.h"
  #include "logging.h"

  #define YY_DECL int dglex( YYSTYPE *yylval, struct dyngeom_parse_vars *dg_parse, yyscan_t yyscanner )
  #define YY_NO_UNPUT

  #ifdef __cplusplus
  #define INPUT_FUNC yyinput
  #else
  #define INPUT_FUNC input
  #endif

void *dgalloc(yy_size_t size, yyscan_t yyscanner)
{
  UNUSED(yyscanner);
  return CHECKED_MALLOC(size, "memory buffer for parser");
}

void *dgrealloc(void *ptr, yy_size_t size, yyscan_t yyscanner)
{
  UNUSED(yyscanner);

  ptr = (void *) realloc((char *) ptr, size);
  if (ptr == NULL)
    mcell_allocfailed("Failed to resize memory buffer for parser.");
  return ptr;
}

void dgfree (void * ptr, yyscan_t yyscanner)
{
  UNUSED(yyscanner);

  free((char *) ptr);
}

%}

/* Define state for parsing comments */
%x IN_COMMENT

/* Reentrant lexer allows safer handling of include files from the parser */
%option reentrant

/* Set up function name prefixes and output file name */
%option prefix="dg"

/* Enable condition state stack */
%option stack

/* Disable yywrap feature */
%option noyywrap

/* Eliminate unused symbols to eliminate compiler warnings */
%option noinput
%option nounput
%option noyy_top_state
%option noyyalloc
%option noyyrealloc
%option noyyfree

R ([0-9]+\.)|([0-9]*\.[0-9]+)|([0-9]+[eE][+-]?[0-9]+)|([0-9]+\.[0-9]*[eE][+-]?[0-9]+)|([0-9]*\.[0-9]+[eE][+-]?[0-9]+)
I ([0-9]+)
IDENT [A-Za-z][A-Za-z0-9\_\.]*

%%


"//"[^\n]*"\n"          { dg_parse->line_num[dg_parse->include_stack_ptr - 1] ++; }

<INITIAL>"/*"   {
                  dg_parse->comment_started = dg_parse->line_num[dg_parse->include_stack_ptr - 1];
                  yy_push_state(IN_COMMENT, yyscanner);
                }
<IN_COMMENT>{
[^*/\n]*        { }
"/"[^*\n]       { }
"*"+[^*/\n]*    { }
"/*"            { yy_push_state(YY_START, yyscanner); }
"*"+"/"         { yy_pop_state(yyscanner); }
<<EOF>>         {
                  dg_parse->line_num[dg_parse->include_stack_ptr - 1]--;
                  mcell_error("Unterminated comment started on line %d",
                              dg_parse->comment_started);
                  return 1;
                }
}

[ \t]+			;

"DEFINE_SURFACE_REGIONS"  return(DEFINE_SURFACE_REGIONS);
"ELEMENT_LIST"            return(INCLUDE_ELEMENTS);
"INCLUDE_ELEMENTS"        return(INCLUDE_ELEMENTS);
"OBJECT"          return(OBJECT);
"POLYGON_LIST"    return(POLYGON_LIST);
"VERTEX_LIST"     return(VERTEX_LIST);
"ELEMENT_CONNECTIONS" return(ELEMENT_CONNECTIONS);
"TRANSLATE"       return(TRANSLATE);
"ROTATE"          return(ROTATE);
"SCALE"           return(SCALE);
"INSTANTIATE"     return(INSTANTIATE);
"INCLUDE_FILE"    return(INCLUDE_FILE);
"PARTITION_X"     return(PARTITION_X);
"PARTITION_Y"     return(PARTITION_Y);
"PARTITION_Z"     return(PARTITION_Z);
"TO"              return(TO);
"STEP"            return(STEP);
{I}             {
                  yylval->llival = atoi(yytext);
                  return LLINTEGER;
                }
{R}             {
                  yylval->dbl = atof(yytext);
                  return REAL;
                }
{IDENT}         {
                  yylval->str = strdup(yytext);
                  return VAR;
                }
\"(([^\"]*)|((\\\")*))*\"  {
                             yylval->str = strdup(yytext);
                             return STR_VALUE;
                           }
[\'\,\(\)\/\-\+\=\^\[\]\{\}\|\<\>\*\#\~\@\:\&\;] {return(yytext[0]);}
.			;
<INITIAL,IN_COMMENT>\n+ {dg_parse->line_num[dg_parse->include_stack_ptr - 1] += yyleng; }
<<EOF>>			{ yyterminate(); }
