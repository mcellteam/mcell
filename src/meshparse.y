%{
  #include "meshparse.h"
  void yyerror(char *);
  int yylex(void);
%}

%union {
  int tok;
  double dbl;
  long long llival;
  struct sym_table *sym;
  struct num_expr_list_head nlist;
}

%token       EQUALS
%token       PARTITION_X
%token       PARTITION_Y
%token       PARTITION_Z
%token <str> VAR
%token <dbl> REAL
%token <llival> LLINTEGER
%token       TO
%token       STEP

%type <tok> partition_dimension
%type <nlist> array_value
%type <nlist> range_spec
%type <nlist> list_range_specs
%type <nlist> array_expr_only
%type <dbl> num_expr

%%

mdl_format:
        mdl_stmt_list
;

mdl_stmt_list:
        mdl_stmt
      | mdl_stmt_list mdl_stmt
;

mdl_stmt: 
        partition_def
;

partition_def:
          partition_dimension '=' array_value         {}
;

partition_dimension:
          PARTITION_X                                 { $$ = X_PARTS; }
        | PARTITION_Y                                 { $$ = Y_PARTS; }
        | PARTITION_Z                                 { $$ = Z_PARTS; }
;

array_value: array_expr_only
;

array_expr_only: '[' list_range_specs ']'             { $$ = $2; }
;

list_range_specs:
          range_spec
        | list_range_specs ',' range_spec             {
                                                          if ($1.value_tail)
                                                          {
                                                            $$ = $1;
                                                            $$.value_count += $3.value_count;
                                                            $$.value_tail->next = $3.value_head;
                                                            $$.value_tail = $3.value_tail;
                                                          }
                                                          else
                                                            $$ = $3;
                                                      }
;

range_spec: '[' num_expr TO num_expr STEP num_expr ']' { generate_range(&$$, $2, $4, $6); }
;

num_expr:  LLINTEGER { $$ = $1; }
         | REAL    { $$ = $1; }
;

%%

void yyerror(char *s) {
    fprintf(stderr, "%s\n", s);
}

int main(void) {
    yyparse();
    return 0;
}
