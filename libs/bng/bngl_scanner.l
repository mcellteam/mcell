%{
  #include "bngl_parser.hpp"
  #include "parser_utils.h"
    
  #define size_t long // eliminate compiler warning
  
  #define LOC() bngllloc.first_line = bngllineno
  
  extern BNG::ASTContext* g_ctx; // defined in bngl_parser.y
  
%}

/* definitions */

/* Enable scanner to generate debug output. disable this for release
 * versions. */
/* %option debug */

/* Eliminate unused symbols to eliminate compiler warnings */
%option nounput
%option noyywrap

%option yylineno

/* Set up function name prefixes and output file name */
%option prefix="bngl"

/* TODO: there was a way to store l=file and line information with each token*/ 

DIGIT    [0-9]
ID       [a-zA-Z_][a-zA-Z0-9_]*
R ([0-9]+\.)|([0-9]*\.[0-9]+)|([0-9]+[eE][+-]?[0-9]+)|([0-9]+\.[0-9]*[eE][+-]?[0-9]+)|([0-9]*\.[0-9]+[eE][+-]?[0-9]+)
I ([0-9]+)

%%

    /* comments and whitespaces */

"#"[^\n]*"\n"   { /* line comment */ /*bngllineno++;*/ }
[ \t]+          ; /* whitespaces are ignored, not sure if this is correct */ 
"\n"            { /*bngllineno++;*/ }
"\r"            ; /* Windows end line character is ignored */


    /* keywords */
     
"begin"         { LOC(); return TOK_BEGIN; }
"end"           { LOC(); return TOK_END; }
"model"         { LOC(); return TOK_MODEL; }
"parameters"    { LOC(); return TOK_PARAMETERS; }
"molecule"      { LOC(); return TOK_MOLECULE; }
"types"         { LOC(); return TOK_TYPES; }
"reaction"      { LOC(); return TOK_REACTION; }
"rules"         { LOC(); return TOK_RULES; }

    /* numbers and identifiers */

{R}             { LOC(); bngllval.dbl = BNG::convert_to_dbl(bngltext); return TOK_DBL; }
{I}             { LOC(); bngllval.llong = BNG::convert_dec_to_llong(bngltext); return TOK_LLONG; }
{ID}            { LOC(); bngllval.str = g_ctx->insert_to_string_pool(yytext); return TOK_ID; }
                        
    /* other characters */
    
"->"            { LOC(); return TOK_ARROW_RIGHT; }
"<->"           { LOC(); return TOK_ARROW_BIDIR; }
[\~\,\!\.\(\)\/\+]  { LOC(); return bngltext[0];}

.               { LOC(); BNG::errs() << "Unexpected character '" << bngltext[0] << "'.\n"; return bngltext[0]; }
%%


#undef size_t // to disable compiler warning

