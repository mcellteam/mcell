%{
  #include <stdio.h>
  #include <string.h>
  #include <stdlib.h>
  #include <unistd.h>
  #include "strfunc.h"
  #include "argparse.h"

  #define YY_DECL int arglex( YYSTYPE *arglval, struct argparse_vars *apvp )

  #define YY_NO_UNPUT

  #ifdef __cplusplus
  #define INPUT_FUNC yyinput
  #else
  #define INPUT_FUNC input
  #endif
%}


%option prefix="arg"
%option outfile="arglex.flex.c"
%option noyywrap

%p 10000
%a 10000 
%o 10000
%e 2000
%n 1000

R ([0-9]+\.)|([0-9]*\.[0-9]+)|([0-9]+[eE][+-]?[0-9]+)|([0-9]+\.[0-9]*[eE][+-]?[0-9]+)|([0-9]*\.[0-9]+[eE][+-]?[0-9]+)
I ([0-9]+)


%%

[ \t]+                  ;
"-help"			{return(HELP_OPT);}
"-info"			{return(INFO_OPT);}
"-logfile"		{return(LOG_FILE_OPT);}
"-logfreq"		{return(LOG_FREQ_OPT);}
"-seed"			{return(SEED_OPT);}
"-iterations"		{return(ITERATIONS_OPT);}
"-checkpoint_infile"	{return(CHECKPOINT_OPT);}
{R}			{sscanf((const char *)argtext,"%lf",&apvp->rval); return(REAL);}
{I}			{sscanf((const char *)argtext,"%d",&apvp->ival); return(INTEGER);}
[[:alnum:][:punct:]]+	 {apvp->cval=my_strdup(argtext); 
			  if(apvp->cval == NULL){
                            printf("Memory allocation error.\n");
                            exit(1);
                          }
			  return(FILE_NAME);}
.			;
<<EOF>>			{return(EOF_TOK);}
