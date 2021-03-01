%{
  #include "bngl_parser.hpp"
  #include "bng/parser_utils.h"
    
  #define size_t long // eliminate compiler warning
  
  #define LOC() bngllloc.first_line = bngllineno
  
  extern BNG::ParserContext* g_ctx; // defined in bngl_parser.y
  
%}

/* definitions */

/* Enable scanner to generate debug output. disable this for release
 * versions. */
/* %option debug */

/* Eliminate unused symbols to eliminate compiler warnings */
%option noyywrap

%option yylineno

/* Set up function name prefixes and output file name */
%option prefix="bngl"

ID [a-zA-Z_][a-zA-Z0-9_]*
R ([0-9]+\.)|([0-9]*\.[0-9]+)|([0-9]+[eE][+-]?[0-9]+)|([0-9]+\.[0-9]*[eE][+-]?[0-9]+)|([0-9]*\.[0-9]+[eE][+-]?[0-9]+)
I ([0-9]+)
STR \"([^\\\"]|\\.)*\"

%%

    /* comments and whitespaces */

"#"[^\n]*"\n"   { /* line comment */ LOC(); return TOK_NL; }
[ \t]+          ; /* whitespaces are ignored, not sure if this is correct */ 
"\n"            { LOC(); return TOK_NL; }
"\r"            ; /* Windows end line character is ignored */
"\\\n"	    		; /* Canceling newline preceded by backslash */
"\\\r\n"        ; /* Canceling newline preceded by backslash */
<<EOF>>         {   
                  // we must return newline for EOF only once otherwise 
                  // the parsing won't end  
                  if (g_ctx->get_eof_returned_as_newline()) {
                    return 0;
                  }
                  else {
                    g_ctx->set_eof_returned_as_newline();
                    LOC();
                    return TOK_NL;
                  }
                }
                
                
    /* keywords */
     
"begin"         { LOC(); return TOK_BEGIN; }
"end"           { LOC(); return TOK_END; }
"model"         { LOC(); return TOK_MODEL; }
"parameters"    { LOC(); return TOK_PARAMETERS; }
"molecule"      { LOC(); return TOK_MOLECULE; }
"compartments"  { LOC(); return TOK_COMPARTMENTS; }
"types"         { LOC(); return TOK_TYPES; }
"reaction"      { LOC(); return TOK_REACTION; }
"rules"         { LOC(); return TOK_RULES; }
"seed"          { LOC(); return TOK_SEED; }
"species"       { LOC(); return TOK_SPECIES; }
"observables"   { LOC(); return TOK_OBSERVABLES; }
"actions"       { LOC(); return TOK_ACTIONS; }

"!CPLX" 		{ LOC(); return TOK_SINGLE_CPLX; }

    /* numbers and identifiers */

{R}             { LOC(); bngllval.dbl = BNG::convert_to_dbl(bngltext); return TOK_DBL; } /* must contain decimal point */
{I}             { LOC(); bngllval.llong = BNG::convert_dec_to_llong(bngltext); return TOK_LLONG; }
{ID}            { LOC(); bngllval.str = g_ctx->insert_to_string_pool(yytext); return TOK_ID; }
{STR}           { LOC(); bngllval.str = g_ctx->insert_to_string_pool(yytext); return TOK_STR; } 
                        
    /* other characters */
    
"->"            { LOC(); return TOK_ARROW_RIGHT; }
"<->"           { LOC(); return TOK_ARROW_BIDIR; }
"=>"            { LOC(); return TOK_ARG_ASSIGN; }

[\~\,\!\?\.\(\)\/\+\-\*\^\{\}\"\=\;\@\:]  { LOC(); return bngltext[0];}

.               { LOC(); errs() << "Unexpected character '" << bngltext[0] << "'.\n"; return bngltext[0]; }
%%


#undef size_t // to disable compiler warning

