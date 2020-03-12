
// for top of bngl_parser.hpp
%code requires {	  
#include "ast.h"
}

// for top of bngl_parser.cpp
%code top {	  
#include "ast.h"
}

// for bngl_parser.hpp
%code provides {
	
void create_ast_context();
BNG::ASTContext* get_ast_context();
void delete_ast_context();

}

// for bngl_parser.cpp
%{
  #include <cstdio>
  #include <cstdarg>
  #include <string>
  #include "parser_utils.h"
	
  // Declare stuff from Flex that Bison needs to know about:
  extern int bngllex();
  
  void bnglerror(char const *s);

  // global context used during parsing
  BNG::ASTContext* g_ctx;
%}


%require "3.0"

// add debug output code to generated parser. disable this for release versions. 
%debug

// write out a header file containing the token defines 
%defines

%error-verbose

// set up function name prefixes and output file name 
%define api.prefix {bngl}

// one rr conflict is expected
%expect 1

%union {
  char* str; // default, required by flex
  double dbl;
  long long llong;
  BNG::ASTExprNode* expr;
  BNG::ASTListNode* list;
}

%token TOK_BEGIN
%token TOK_END
%token TOK_MODEL
%token TOK_PARAMETERS
%token TOK_MOLECULE
%token TOK_TYPES
%token TOK_REACTION
%token TOK_RULES

%token <str> TOK_ID
%token <dbl> TOK_DBL
%token <llong> TOK_LLONG

%token TOK_ARROW_RIGHT
%token TOK_ARROW_BIDIR

%type <expr> expr

%start start

%%

// TODO: error recovery 

start: 
      section_list
    | TOK_BEGIN TOK_MODEL section_list TOK_END TOK_MODEL
    | /* empty - end of file */
;

section_list:
      section_list section
    | section

section:
      TOK_BEGIN TOK_PARAMETERS parameter_list_maybe_empty TOK_END TOK_PARAMETERS
    | TOK_BEGIN TOK_MOLECULE TOK_TYPES molecule_list_maybe_empty TOK_END TOK_MOLECULE TOK_TYPES
    | TOK_BEGIN TOK_REACTION TOK_RULES reaction_list_maybe_empty TOK_END TOK_REACTION TOK_RULES
;

// ---------------- parameters ------------------- 
parameter_list_maybe_empty:
      parameter_list
    | /* empty */
;

parameter_list:
      parameter_list parameter
    | parameter

parameter:
      TOK_ID expr {
    	g_ctx->symtab.insert($1, $2);
      }
;
      
// ---------------- molecules -------------------     
molecule_list_maybe_empty:
      molecule_list
    | /* empty */
;

// left recursion is preferred 
molecule_list:
      molecule_list molecule
    | molecule
;

// fully general specification, might contain information on bonds, checked later in semantic checks 
molecule:
      TOK_ID '(' component_list_maybe_empty ')'
;

component_list_maybe_empty:
      component_list
    | /* empty */
;

// this rule allows bonds even in molecule types section, this is checked later in semantic checks 
component_list:
      component_list ',' component
    | component
;

component:
      TOK_ID component_name_with_states_maybe_empty bond_maybe_empty
; 

component_name_with_states_maybe_empty:
      component_state_list
    | /* empty */
    
component_state_list:
      component_state_list component_state
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
    

// ---------------- reactions ------------------- 
reaction_list_maybe_empty:
      reaction_list
    | /* empty */
;

reaction_list:
      reaction_list reaction
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

// ---------------- expressions --------------------- 
// TODO: expressions are just IDs and constants for now
expr:
      TOK_ID { 
        $$ = g_ctx->new_id_node($1); 
        free($1);
        $1 = nullptr;
      } 
    | TOK_DBL { 
        $$ = g_ctx->new_dbl_node($1); 
      } 
    | TOK_LLONG { 
        $$ = g_ctx->new_llong_node($1); 
      } 

%%

void bnglerror(char const *s) {
  BNG::errs() << s << "\n";
}

void create_ast_context() {
  assert(g_ctx == nullptr);
  g_ctx = new BNG::ASTContext();
}

BNG::ASTContext* get_ast_context() {
  return g_ctx;
}

void delete_ast_context() {
  delete g_ctx;
  g_ctx = nullptr;
}
