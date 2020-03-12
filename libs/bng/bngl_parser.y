%{
/* C declarations */

  #include <cstdio>
  #include <cstdarg>
  #include <string>
  
  // Declare stuff from Flex that Bison needs to know about:
  extern int bngllex();
  //extern int bnglparse();
  //extern FILE *bnglin;
 
  void bnglerror(const char *s);
%}

/* bison declarations */

/* Require bison 3.0 or later */
%require "3.0"

/* Add debug output code to generated parser. disable this for release
 * versions. */
%debug

/* write out a header file containing the token defines */
%defines

/* verbose error messages */
%error-verbose

/* set up function name prefixes and output file name */
%define api.prefix  {bngl}

%expect 1

%union {
double dbl;
long long llong;
const char* str; // default, required by flex
}

%token TOK_BEGIN
%token TOK_END
%token TOK_MODEL
%token TOK_PARAMETERS
%token TOK_MOLECULE
%token TOK_TYPES
%token TOK_REACTION
%token TOK_RULES

%token TOK_ID
%token TOK_DBL
%token TOK_LLONG

%token TOK_ARROW_RIGHT
%token TOK_ARROW_BIDIR

%start start

%%
/* grammar rules */

/* TODO: error recovery ?*/

start: 
      sections_list
    | TOK_BEGIN TOK_MODEL sections_list TOK_END TOK_MODEL
    | /* empty - end of file */
;

sections_list:
      sections_list section
    | section

section:
      TOK_BEGIN TOK_PARAMETERS parameters_list_maybe_empty TOK_END TOK_PARAMETERS
    | TOK_BEGIN TOK_MOLECULE TOK_TYPES molecules_list_maybe_empty TOK_END TOK_MOLECULE TOK_TYPES
    | TOK_BEGIN TOK_REACTION TOK_RULES reactions_list_maybe_empty TOK_END TOK_REACTION TOK_RULES
;

/* ---------------- parameters ------------------- */    
parameters_list_maybe_empty:
      parameters_list
    | /* empty */
;

parameters_list:
      parameters_list parameter
    | parameter

parameter:
      TOK_ID expr
;
      
/* ---------------- molecules ------------------- */    
molecules_list_maybe_empty:
      molecules_list
    | /* empty */
;

/* left recursion is preferred */
molecules_list:
      molecules_list molecule
    | molecule
;

/* fully general specification, might contain information on bonds, checked later fir semantic checks */
molecule:
      TOK_ID '(' components_list_maybe_empty ')'
;

components_list_maybe_empty:
      components_list
    | /* empty */
;

/* this rule allows bonds even in molecyle types section, this is checked later in semantic checks */
components_list:
      components_list ',' component
    | component
;

component:
      TOK_ID component_name_with_states_maybe_empty bond_maybe_empty
; 

component_name_with_states_maybe_empty:
      component_states_list
    | /* empty */
    
component_states_list:
      component_states_list component_state
    | component_state
;

component_state:
      '~' TOK_ID
    | '~' TOK_LLONG
;

bond_maybe_empty:
      '!' TOK_LLONG
    | '!' '+'
    | /* empty */
;
    

/* ---------------- reactions ------------------- */    
reactions_list_maybe_empty:
      reactions_list
    | /* empty */
;

reactions_list:
      reactions_list reaction
    | reaction
;

reaction:
      reaction_side reaction_direction reaction_side_maybe_empty rates
;

reaction_side_maybe_empty:
      reaction_side
    | /* empty */
;

reaction_side:
      reaction_side '+' molecule
    | reaction_side '.' molecule
    | molecule
;

reaction_direction:
      TOK_ARROW_RIGHT
    | TOK_ARROW_BIDIR

rates:
      expr ',' expr
    | expr

/* ---------------- expressions --------------------- */
expr:
      TOK_ID
    | TOK_DBL
    | TOK_LLONG 
    

%%

/* additional C code */
extern int bngllineno;


void bnglerror(char const *s) {
  fprintf (stderr, "Line %d: %s\n", (int)bngllineno, s);
}

void bnglerror_fmt(char const *fmt, ...) {
  char buffer[256];
  va_list args;
  va_start (args, fmt);
  vsnprintf (buffer, 256,fmt, args);
  bnglerror(buffer);
  va_end(args);
}

