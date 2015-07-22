%{
  #include "meshparse.h"
  #include "y.tab.h"
  #include <stdlib.h>
  void yyerror(char *);
%}

%option noinput
%option nounput

R ([0-9]+\.)|([0-9]*\.[0-9]+)|([0-9]+[eE][+-]?[0-9]+)|([0-9]+\.[0-9]*[eE][+-]?[0-9]+)|([0-9]*\.[0-9]+[eE][+-]?[0-9]+)
I ([0-9]+)
%%

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
[ \t]             ;
[\'\,\(\)\/\-\+\=\^\[\]\{\}\|\<\>\*\#\~\@\:\&\;] {return(yytext[0]);}
%%

int yywrap(void) {
    return 1;
}
