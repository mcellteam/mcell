
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
	
namespace BNG {
void create_ast_context();
ASTContext* get_ast_context();
void delete_ast_context();
}

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
  
  // used to process TOK_ID tokens
  std::string to_str_and_delete(char*& ptr) {
	std::string res = ptr;
	free(ptr);
	ptr = nullptr;
	return res;
  }
  
%}


%require "3.0"

// add debug output code to generated parser. disable this for release versions. 
//%debug

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
  BNG::ASTExprNode* expr_node;
  BNG::ASTListNode* list_node;
  BNG::ASTStrNode* str_node;
  BNG::ASTComponentNode* component_node;
  BNG::ASTMoleculeNode* molecule_node;
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

%type <expr_node> expr
%type <str_node> bond_maybe_empty
%type <str_node> component_state 
%type <component_node> component
%type <list_node> component_state_list_maybe_empty
%type <list_node> component_state_list
%type <list_node> component_list_maybe_empty
%type <list_node> component_list
%type <molecule_node> molecule
%type <list_node> molecule_list_maybe_empty
%type <list_node> molecule_list

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
    | TOK_BEGIN TOK_MOLECULE TOK_TYPES molecule_list_maybe_empty TOK_END TOK_MOLECULE TOK_TYPES {
    	g_ctx->symtab.insert_molecule_declarations($4, g_ctx);
      }
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
    	g_ctx->symtab.insert($1, $2, g_ctx);
      }
;
      
// ---------------- molecules -------------------     
molecule_list_maybe_empty:
      molecule_list
    | /* empty */ {
       	$$ = g_ctx->new_list_node();
      }
;

// left recursion is preferred 
molecule_list:
      molecule_list molecule {
       	$1->append($2);
        $$ = $1;
      }
    | molecule {
  	    $$ = g_ctx->new_list_node()->append($1);
      }
;

// fully general specification, might contain information on bonds, checked later in semantic checks 
molecule:
      TOK_ID '(' component_list_maybe_empty ')' {
    	$$ = g_ctx->new_molecule_node(
    		to_str_and_delete($1),		
			$3
    	);	
      }
;

component_list_maybe_empty:
      component_list
    | /* empty */ {
      	$$ = g_ctx->new_list_node();
      }
;

component_list:
      component_list ',' component {
      	$1->append($3);
  	    $$ = $1;
      }
    | component {
  	    $$ = g_ctx->new_list_node()->append($1);
      }
;

component:
      TOK_ID component_state_list_maybe_empty bond_maybe_empty {
    	$$ = g_ctx->new_component_node(
    		to_str_and_delete($1),	
			$2,
			$3
    	);
      }
; 

component_state_list_maybe_empty:
      component_state_list
    | /* empty */ {
    	$$ = g_ctx->new_list_node();
      }
;
    
component_state_list:
      component_state_list component_state {
    	$1->append($2);
	    $$ = $1;
      }
    | component_state {
	    $$ = g_ctx->new_list_node()->append($1);
      }
;

component_state:
      '~' TOK_ID {
      	$$ = g_ctx->new_str_node(to_str_and_delete($2));
      }
    | '~' TOK_LLONG {
      	$$ = g_ctx->new_str_node($2);
      }
;

bond_maybe_empty:
      '!' TOK_LLONG {
    	$$ = g_ctx->new_str_node($2);
      }
    | '!' '+' {
    	$$ = g_ctx->new_str_node(BNG::ANY_BOND);
      }
    | /* empty */ {
    	$$ = g_ctx->new_str_node("");
    }
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
        $$ = g_ctx->new_id_node(to_str_and_delete($1)); 
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

namespace BNG {
void create_ast_context() {
  assert(g_ctx == nullptr);
  g_ctx = new BNG::ASTContext();
}

ASTContext* get_ast_context() {
  return g_ctx;
}

void delete_ast_context() {
  delete g_ctx;
  g_ctx = nullptr;
}

}