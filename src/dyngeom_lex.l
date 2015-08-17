%{

  #include "dyngeom_parse_extras.h"
  #include "dyngeom_yacc.h"
  #include <stdlib.h>
  void yyerror(char *);
%}

%option noinput
%option nounput

R ([0-9]+\.)|([0-9]*\.[0-9]+)|([0-9]+[eE][+-]?[0-9]+)|([0-9]+\.[0-9]*[eE][+-]?[0-9]+)|([0-9]*\.[0-9]+[eE][+-]?[0-9]+)
I ([0-9]+)
IDENT [A-Za-z][A-Za-z0-9\_\.]*

%x IN_COMMENT

%%


"//"[^\n]*"\n"        // eat this style of comment (i.e. the double slash);
<INITIAL>{
"/*"                  BEGIN(IN_COMMENT);
}
<IN_COMMENT>{
"*/"                  BEGIN(INITIAL);
[^*\n]+               // eat comment in chunks
"*"                   // eat the lone star
\n                    yylineno++;
}

"OBJECT"          return(OBJECT);
"POLYGON_LIST"    return(POLYGON_LIST);
"VERTEX_LIST"     return(VERTEX_LIST);
"ELEMENT_CONNECTIONS" return(ELEMENT_CONNECTIONS);
"TRANSLATE"       return(TRANSLATE);
"ROTATE"          return(ROTATE);
"SCALE"           return(SCALE);
"INSTANTIATE"     return(INSTANTIATE);
"PARTITION_X"     return(PARTITION_X);
"PARTITION_Y"     return(PARTITION_Y);
"PARTITION_Z"     return(PARTITION_Z);
"TO"              return(TO);
"STEP"            return(STEP);
{I}             {
                  yylval.llival = atoi(yytext);
                  return LLINTEGER;
                }
{R}             {
                  yylval.dbl = atof(yytext);
                  return REAL;
                }
{IDENT}         {
                  yylval.str = strdup(yytext);
                  return VAR;
                }
[ \t]             ;
[\n\'\,\(\)\/\-\+\=\^\[\]\{\}\|\<\>\*\#\~\@\:\&\;] {return(yytext[0]);}

%%

int yywrap(void) {
    return 1;
}
