%{
  #include "bngl_parser.hpp"
  #include "parser_utils.h"
    
  #define size_t long // eliminate compiler warning
%}

/* definitions */

/* Enable scanner to generate debug output. disable this for release
 * versions. */
%option debug

/* Eliminate unused symbols to eliminate compiler warnings */
%option nounput
%option noyywrap

/* Set up function name prefixes and output file name */
%option prefix="bngl"

/* TODO: there was a way to store l=file and line information with each token*/ 

DIGIT    [0-9]
ID       [a-zA-Z_][a-zA-Z0-9_]*
R ([0-9]+\.)|([0-9]*\.[0-9]+)|([0-9]+[eE][+-]?[0-9]+)|([0-9]+\.[0-9]*[eE][+-]?[0-9]+)|([0-9]*\.[0-9]+[eE][+-]?[0-9]+)
I ([0-9]+)

%%

    /* comments and whitespaces */

"#"[^\n]*"\n"   { /* line comment */ bngllineno++; }
[ \t]+          ; /* whitespaces are ignored, not sure if this is correct */ 
"\n"            { bngllineno++; }
"\r"            ; /* Windows end line character is ignored */


    /* keywords */
     
"begin"         { return TOK_BEGIN; }
"end"           { return TOK_END; }
"model"         { return TOK_MODEL; }
"parameters"    { return TOK_PARAMETERS; }
"molecule"      { return TOK_MOLECULE; }
"types"         { return TOK_TYPES; }
"reaction"      { return TOK_REACTION; }
"rules"         { return TOK_RULES; }

    /* numbers and identifiers */

{R}             { bngllval.dbl = BNG::convert_to_dbl(bngltext); return TOK_DBL; }
{I}             { bngllval.llong = BNG::convert_dec_to_llong(bngltext); return TOK_LLONG; }
{ID}            { bngllval.str = strdup(yytext); return TOK_ID; }
                        
    /* other characters */
    
"->"            { return TOK_ARROW_RIGHT; }
"<->"           { return TOK_ARROW_BIDIR; }
[\~\,\!\.\(\)\/\+]  { return bngltext[0];}

.               { BNG::errs() << "Unexpected character '" << bngltext[0] << "'.\n"; return bngltext[0]; }
%%


#undef size_t // to disable compiler warning

