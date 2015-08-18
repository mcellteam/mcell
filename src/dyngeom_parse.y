%{
  #include <stdlib.h>
  #include <string.h>

  #include "mcell_objects.h"
  #include "logging.h"
  #include "dyngeom_parse_extras.h"

  #define EPS_C 1e-12

  typedef void *yyscan_t;
  void yyerror(char *);
  int yylex(void);
  extern FILE *yyin;
  int yyparse(void);

  void object_list_singleton(struct object_list *head, struct object *objp) {
    objp->next = NULL;
    head->obj_tail = head->obj_head = objp;
  }

  void add_object_to_list(struct object_list *head, struct object *objp) {
    objp->next = NULL;
    head->obj_tail = head->obj_tail->next = objp;
  }

  int create_dg_parse() {
    dg_parse = (struct dyngeom_parse_vars *)malloc(sizeof(struct dyngeom_parse_vars));
    memset(dg_parse, 0, sizeof(struct dyngeom_parse_vars));
    init_top_level_objs(dg_parse);
    return 0;
  }

  int parse_dg(char *dynamic_geometry_filename) {
    FILE *fp=fopen(dynamic_geometry_filename,"r");
    if(!fp)
    {
      no_printf("Couldn't open file for reading\n");
      return 1;
    }
    yyin=fp;
    setup_root_obj_inst(dg_parse);
    yyparse();
    fclose(fp);
    return 0;
  }

%}

%union {
  int tok;
  double dbl;
  char *str;
  long long llival;
  struct sym_table *sym;
  struct vector3 *vec3;
  struct num_expr_list_head nlist;
  struct object *obj;
  struct object_list obj_list;
}

%token       OBJECT
%token       POLYGON_LIST
%token       VERTEX_LIST
%token       ELEMENT_CONNECTIONS
%token       INSTANTIATE
%token       PARTITION_X
%token       PARTITION_Y
%token       PARTITION_Z
%token <str> VAR
%token <dbl> REAL
%token <llival> LLINTEGER
%token       TO
%token       STEP
%token       TRANSLATE
%token       SCALE
%token       ROTATE
%token       INCLUDE_ELEMENTS
%token       DEFINE_SURFACE_REGIONS

%type <str> new_object_name
%type <dbl> intOrReal
%type <dbl> num_value
%type <sym> existing_num_var
%type <obj> existing_object_ref
%type <obj> object_def
%type <obj> object_ref
%type <obj_list> list_objects
%type <str> var
%type <obj> polygon_list_def
%type <sym> existing_object
%type <vec3> point
%type <vec3> point_or_num
%type <obj> meta_object_def
%type <sym> new_object
%type <tok> partition_dimension
%type <nlist> array_value
%type <nlist> range_spec
%type <nlist> list_range_specs
%type <nlist> array_expr_only
%type <dbl> num_expr
%type <dbl> num_expr_only
%type <dbl> arith_expr

/* Operator associativities and precendences */
%right '='
%left '&' ':'
%left '+' '-'
%left '*' '/'
%left '^'
%left UNARYMINUS

%%

/* Start */
mdl_format:
        mdl_stmt_list
;

mdl_stmt_list:
        mdl_stmt
      | mdl_stmt_list mdl_stmt
;

mdl_stmt:
        partition_def
      | physical_object_def
      | instance_def
;

/* =================================================================== */
/* Utility definitions */
var: VAR
;

existing_object: var                                 {no_printf("existing_object\n");}

;

point: array_value                                   {no_printf("point\n");}
;

point_or_num: point                                  {no_printf("point_or_num\n");}
            | num_expr_only                          { }
;

list_range_specs:
          range_spec
        | list_range_specs ',' range_spec            { }
;

range_spec: num_expr                                 { no_printf("range_spec\n"); }
        | '[' num_expr TO num_expr STEP num_expr ']' { generate_range(&$$, $2, $4, $6); }
;

/* Object type: Meta-objects */
meta_object_def:
        new_object OBJECT
        start_object
          list_objects
          list_opt_object_cmds
        end_object                                   {
                                                         struct object *the_object = (struct object *) $1->value;
                                                         the_object->object_type = META_OBJ;
                                                         add_child_objects(the_object, $4.obj_head, $4.obj_tail);
                                                         $$ = the_object;
                                                     }

;

list_objects:
        object_ref                                   { object_list_singleton(& $$, $1); }
      | list_objects object_ref                      { $$ = $1; add_object_to_list(& $$, $2); }
;

object_ref: existing_object_ref
          | object_def
;

existing_object_ref:
        new_object OBJECT existing_object
        start_object
          list_opt_object_cmds
        end_object                                   {
                                                         $$ = (struct object *) $1->value;
                                                         $$->object_type = POLY_OBJ;
                                                         dg_create_region(dg_parse->reg_sym_table, $$, "ALL");
                                                     }
;

new_object_name: var                                 {no_printf("new_object_name\n");}
;

/* =================================================================== */
/* Instance definitions */

instance_def:
          INSTANTIATE                                { no_printf("INSTANTIATE\n"); dg_parse->current_object = dg_parse->root_instance; }
          meta_object_def                            {
                                                        no_printf("meta_object_def\n");
                                                        add_child_objects(dg_parse->root_instance, $3, $3);
                                                        dg_parse->current_object = dg_parse->root_object;
                                                     }
;

/* =================================================================== */
/* Object type definitions */

physical_object_def: object_def                      {no_printf("physical_object_def\n");}
;

object_def: meta_object_def
          | polygon_list_def
;

new_object: var                                      {no_printf("new_object\n"); $$ = dg_start_object(dg_parse, $1);}
;

start_object: '{'                                    {no_printf("start_object\n");}
;

end_object: '}'                                      {no_printf("end_object\n"); dg_finish_object(dg_parse);}
;

list_opt_object_cmds:
        | list_opt_object_cmds opt_object_cmd
;

opt_object_cmd: transformation
;

transformation:
          TRANSLATE '=' point                        {no_printf("TRANSLATE\n");}
        | SCALE '=' point_or_num                     {no_printf("SCALE\n");}
        | ROTATE '=' point ',' num_expr              {no_printf("ROTATE\n");}
;

/* Object type: Polygons */
polygon_list_def:
          new_object_name POLYGON_LIST               {no_printf("POLYGON_LIST");}
          start_object
            vertex_list_cmd
            element_connection_cmd                   { $<obj>$->object_type = POLY_OBJ; }
            list_opt_polygon_object_cmds
            list_opt_object_cmds
          '}'
                                                     { $$ = (struct object *) $<obj>6; }
;

vertex_list_cmd:
          VERTEX_LIST {no_printf("vertex_list_command\n");}
          '{' list_points '}'            
;

single_vertex: point                                 {no_printf("single_vertex\n");}
;

list_points: single_vertex                           {no_printf("list_points\n");}
           | list_points single_vertex
;

element_connection_cmd:
          ELEMENT_CONNECTIONS                        {no_printf("element_connection_cmd\n");}
          '{' list_element_connections '}'
;

list_element_connections:
          element_connection                         {no_printf("element_connection\n");}
        | list_element_connections
          element_connection                         {no_printf("element_connection\n");}
;

element_connection: array_value                      {no_printf("element_connection\n");}
;

list_opt_polygon_object_cmds:
        | list_opt_polygon_object_cmds
          opt_polygon_object_cmd
;


opt_polygon_object_cmd:
          in_obj_define_surface_regions
;

element_specifier_list:
          element_specifier
        | element_specifier_list
          element_specifier                          {no_printf("element_specifier\n");}
;

element_specifier:
          incl_element_list_stmt
;

incl_element_list_stmt:
          INCLUDE_ELEMENTS '='
          '[' list_element_specs ']'                 {no_printf("incl_element_list_stmt\n");}
;

list_element_specs:
          element_spec
        | list_element_specs ',' element_spec
;

element_spec: num_expr                               {no_printf("element_spec\n");}
            | num_expr TO num_expr
;


in_obj_define_surface_regions:
          DEFINE_SURFACE_REGIONS '{'
            list_in_obj_surface_region_defs
          '}'
;

list_in_obj_surface_region_defs:
          in_obj_surface_region_def
        | list_in_obj_surface_region_defs
          in_obj_surface_region_def
;

in_obj_surface_region_def:
          new_region '{'
            element_specifier_list
          '}'
;

new_region: var                                      {no_printf("new_region\n");}
;

/* =================================================================== */
/* Partitions */

partition_def:
          partition_dimension '=' array_value
;

partition_dimension:
          PARTITION_X                                { $$ = X_PARTS; }
        | PARTITION_Y                                { $$ = Y_PARTS; }
        | PARTITION_Z                                { $$ = Z_PARTS; }
;

/* =================================================================== */
/* Expressions */

array_value: array_expr_only                         {no_printf("array_value\n");}
;

array_expr_only: '[' list_range_specs ']'            { no_printf("array_expr_only\n"); $$ = $2; }
;

num_expr: num_value                                  { no_printf("num_expr\n"); }
        | arith_expr
;

num_value: intOrReal                                 { no_printf("num_value\n"); }
         | existing_num_var                          { no_printf("num_value\n"); $$ = *(double *) $1->value; }
;

intOrReal: LLINTEGER                                 { no_printf("LLINTEGER\n"); $$ = $1; }
         | REAL                                      { no_printf("REAL\n"); }
;

num_expr_only: intOrReal
             | arith_expr
;


existing_num_var: var                                { no_printf("existing_num_var\n"); }
;

arith_expr:
        '(' num_expr ')'                             { $$ = $2; }
      | '-' num_expr %prec UNARYMINUS                { $$ = -$2; }
;


%%

void yyerror(char *s) {
  mcell_error("%s\n", s);
}

int main(int argc, char *argv[])
{
  create_dg_parse();
  parse_dg(argv[1]);
}
