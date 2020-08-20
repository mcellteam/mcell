// FIXMEs: 
// 1) BNGL (BNG2.pl) requires whitespace to be precisely defined,
//    need to update the parser and tests
// 2) Parentheses are not required for molecule types and molecule instances


// for top of bngl_parser.hpp
%code requires {      
#include "bng/ast.h"
}

// for top of bngl_parser.cpp
%code top {      
#include "bng/ast.h"
}

// for bngl_parser.hpp
%code provides {
    
namespace BNG {
  void create_parser_context();
  ParserContext* get_parser_context();
  void delete_parser_context();
}

}

// for bngl_parser.cpp
%{
  #include <cstdio>
  #include <cstdarg>
  #include <string>
  #include "bng/parser_utils.h"
    
  // Declare stuff from Flex that Bison needs to know about:
  extern int bngllex();
  
  void bnglerror(char const *s);

  // global context used during parsing
  BNG::ParserContext* g_ctx;
%}


%require "3.0"

// add debug output code to generated parser 
// one also needs to set bngldebug to 1
// extern int bngldebug;
// bngldebug = 1;
//%debug

// write out a header file containing the token defines 
%defines

%error-verbose

// set up function name prefixes and output file name 
%define api.prefix {bngl}

// extend yylval (bngllval) with the possibility to store line)
%locations

// One shift-reduce conflict is expected, the reason is that 
// the beginning of the molecule types section is resolved
// only after the first molecule type declaration is parsed as shown here.
// So we need to use GLR parser that is able to do arbitrary lookahead.
//
// begin molecule types
//   A(a)
// end molecule types
// begin reaction rules
//   A(a) +|.
// end reaction rules
//
%glr-parser
%expect 3

%union {
  const char* str;
  double dbl;
  long long llong;
  bool boolean;
  BNG::ASTExprNode* expr_node;
  BNG::ASTListNode* list_node;
  BNG::ASTStrNode* str_node;
  BNG::ASTComponentNode* component_node;
  BNG::ASTMoleculeNode* molecule_node;
}

%token TOK_BEGIN "begin"
%token TOK_END "end"
%token TOK_MODEL "model"
%token TOK_PARAMETERS "parameters"
%token TOK_MOLECULE "molecule"
%token TOK_TYPES "types"
%token TOK_REACTION "reaction"
%token TOK_RULES "rules"
%token TOK_SEED "seed"
%token TOK_SPECIES "species"
%token TOK_OBSERVABLES "observables"

// special token to swithc parser to mode where it parses a single complex instance
%token TOK_SINGLE_CPLX "@CPLX"

%token <str> TOK_ID "identifier"
%token <dbl> TOK_DBL "floating point constant"
%token <llong> TOK_LLONG "integer constant"

%token TOK_ARROW_RIGHT "->"
%token TOK_ARROW_BIDIR "<->"
%token TOK_ARG_ASSIGN  "=>"

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
%type <list_node> rxn_rule_side_maybe_empty
%type <list_node> rxn_rule_side
%type <list_node> rates
%type <boolean> rxn_rule_direction
%type <list_node> cplx_instance

// operator associativities and precendences
// unary minus has really lower precendence than power 
%left '+' '-'
%left '*' '/'
%left UNARYPLUS UNARYMINUS
%left '^'   

%%

// TODO: error recovery 
// this start symbol is for the whole BNGL file
start_bngl: 
    // whole file without begin model .. end model markers
      section_list action_call_list_maybe_empty
      
    // whole file with begin model .. end model markers
    | TOK_BEGIN TOK_MODEL section_list TOK_END TOK_MODEL action_call_list_maybe_empty
    
    // single complex to be parsed, prefixed by a unique string 
    | TOK_SINGLE_CPLX cplx_instance {
    	g_ctx->single_cplx_instance = $2;
    }
    // empty file
    | 
;


section_list:
      section_list section
    | section
;
      
section:
      TOK_BEGIN TOK_PARAMETERS parameter_list_maybe_empty TOK_END TOK_PARAMETERS
    | TOK_BEGIN TOK_MOLECULE TOK_TYPES molecule_list_maybe_empty TOK_END TOK_MOLECULE TOK_TYPES {
        g_ctx->symtab.insert_molecule_declarations($4, g_ctx);
      }
    | TOK_BEGIN TOK_REACTION TOK_RULES rxn_rule_list_maybe_empty TOK_END TOK_REACTION TOK_RULES 
    | TOK_BEGIN TOK_SEED TOK_SPECIES seed_species_list_maybe_empty TOK_END TOK_SEED TOK_SPECIES
    | TOK_BEGIN TOK_SPECIES seed_species_list_maybe_empty TOK_END TOK_SPECIES
    | TOK_BEGIN TOK_OBSERVABLES observables_list_maybe_empty TOK_END TOK_OBSERVABLES
;

// ---------------- parameters ------------------- 
parameter_list_maybe_empty:
      parameter_list
    | /* empty */
;

parameter_list:
      parameter_list parameter
    | parameter
;
      
parameter:
      TOK_ID expr {
        g_ctx->symtab.insert($1, $2, g_ctx);
      }
    | TOK_ID '=' expr {
        g_ctx->symtab.insert($1, $3, g_ctx);
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
        $$ = g_ctx->new_molecule_node($1, $3, @1);    
      }
    | TOK_ID {
        // no components neither parentheses
        $$ = g_ctx->new_molecule_node($1, g_ctx->new_list_node(), @1);    
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
        $$ = g_ctx->new_component_node($1, $2, $3, @1);
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
        $$ = g_ctx->new_str_node($2, @2);
      }
    | '~' TOK_LLONG {
        $$ = g_ctx->new_str_node($2, @2);
      }
;

bond_maybe_empty:
      '!' TOK_LLONG {
        $$ = g_ctx->new_str_node($2, @2);
      }
    | '!' '+' {
        $$ = g_ctx->new_str_node(BNG::BOND_STR_BOUND, @2);
      }
    | '!' '?' {
        $$ = g_ctx->new_str_node(BNG::BOND_STR_ANY, @2);
      }      
    | /* empty */ {
        $$ = g_ctx->new_empty_str_node();
    }
;
    
// ---------------- rxn_rules ------------------- 
rxn_rule_list_maybe_empty:
      rxn_rule_list
    | /* empty */ 
;

rxn_rule_list:
      rxn_rule_list rxn_rule 
    | rxn_rule 
;

rxn_rule:
      rxn_rule_side rxn_rule_direction rxn_rule_side_maybe_empty rates {
         
        BNG::ASTRxnRuleNode* n = g_ctx->new_rxn_rule_node($1, $2, $3, $4);
        g_ctx->add_rxn_rule(n);
      }
;

rxn_rule_side_maybe_empty:
      rxn_rule_side
    | /* empty */ {
        $$ = g_ctx->new_list_node();
      }
;

rxn_rule_side:
      rxn_rule_side '+' molecule {
        $1->append(g_ctx->new_separator_node(BNG::SeparatorType::Plus, @2));
        $1->append($3);
        $$ = $1;
      }
    | rxn_rule_side '.' molecule {
        $1->append(g_ctx->new_separator_node(BNG::SeparatorType::Dot, @2));
        $1->append($3);
        $$ = $1;
      }
    | molecule {
        $$ = g_ctx->new_list_node()->append($1);
      }
;

rxn_rule_direction:
      TOK_ARROW_RIGHT {
        $$ = false;
      }
    | TOK_ARROW_BIDIR {
        $$ = true;
      }
;

rates:
      rates ',' expr {
         $1->append($3);
         $$ = $1; 
      }
    | expr {
        $$ = g_ctx->new_list_node()->append($1);
      }
;

// ---------------- seed species ---------------------

seed_species_list_maybe_empty:
    seed_species_list
    | /* empty */ 
;

seed_species_list:
      seed_species_list seed_species_item 
    | seed_species_item
;

seed_species_item:
    cplx_instance expr {
    
       BNG::ASTSeedSpeciesNode* n = g_ctx->new_seed_species_node($1, $2); 
       g_ctx->add_seed_species(n);
    }
;

// similar to rxn_rule_side only contains one complex instance
cplx_instance:
      cplx_instance '.' molecule {
        $1->append(g_ctx->new_separator_node(BNG::SeparatorType::Dot, @2));
        $1->append($3);
        $$ = $1;
      }
    | molecule {
        $$ = g_ctx->new_list_node()->append($1);
      }
;
      
      
// ---------------- observables ---------------------
// ignored for now
      
observables_list_maybe_empty:
      observables_list
    | /* empty */ 
;

observables_list:
      observables_list observables_item 
    | observables_item
;
      
observables_item:
      TOK_ID TOK_ID cplx_instance_list
;

cplx_instance_list:
      cplx_instance_list ',' cplx_instance
    | cplx_instance
      
// ---------------- action calls ------------------
// ignored

action_call_list_maybe_empty:
      action_call_list
    | /* empty */
;

action_call_list:
      action_call_list action_call
    | action_call
;

action_call:
      TOK_ID '(' action_arg_list_maybe_empty ')'

action_arg_list_maybe_empty:
      '{' action_arg_list '}'
    | /* empty */
;

action_arg_list:
      action_arg_list ',' action_arg
    | action_arg
;

action_arg:
      TOK_ID TOK_ARG_ASSIGN expr_or_str
;
 
expr_or_str:
      expr 
    | '"' TOK_ID '"'  // for now just id, we will need full string parsing in scanner      
;
      
// ---------------- expressions --------------------- 
// TODO: expressions are just IDs and constants for now
expr:
      TOK_ID                        { $$ = g_ctx->new_id_node($1, @1); } 
    | TOK_DBL                       { $$ = g_ctx->new_dbl_node($1, @1); } 
    | TOK_LLONG                     { $$ = g_ctx->new_llong_node($1, @1); }
    | '(' expr ')'                  { $$ = $2; } 
    | expr '+' expr                 { $$ = g_ctx->new_expr_node($1, BNG::ExprType::Add, $3, @2); }
    | expr '-' expr                 { $$ = g_ctx->new_expr_node($1, BNG::ExprType::Sub, $3, @2); }
    | expr '*' expr                 { $$ = g_ctx->new_expr_node($1, BNG::ExprType::Mul, $3, @2); }
    | expr '/' expr                 { $$ = g_ctx->new_expr_node($1, BNG::ExprType::Div, $3, @2); }
    | expr '^' expr                 { $$ = g_ctx->new_expr_node($1, BNG::ExprType::Pow, $3, @2); }
    | '+' expr %prec UNARYPLUS      { $$ = g_ctx->new_expr_node($2, BNG::ExprType::UnaryPlus, nullptr, @1); }
    | '-' expr %prec UNARYMINUS     { $$ = g_ctx->new_expr_node($2, BNG::ExprType::UnaryMinus, nullptr, @1); }
;
    
%%

void bnglerror(char const *s) {
  BNG::errs_loc() << s << "\n";
  g_ctx->inc_error_count();
}

namespace BNG {
void create_parser_context() {
  assert(g_ctx == nullptr);
  g_ctx = new ParserContext();
}

ParserContext* get_parser_context() {
  return g_ctx;
}

void delete_parser_context() {
  delete g_ctx;
  g_ctx = nullptr;
}

}