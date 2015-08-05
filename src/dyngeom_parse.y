%{
  #include <stdio.h>
  #include <stdlib.h>
  #include <string.h>
  #include <math.h>

  #include "sym_table.h"
  #include "mcell_objects.h"
  #include "mcell_structs.h"
  #include "mcell_misc.h"
  #include "mcell_objects.h"
  #include "dyngeom_parse_extras.h"

  #define EPS_C 1e-12

  typedef void *yyscan_t;
  void yyerror(char *);
  int yylex(void);
  extern FILE *yyin;
  int yyparse(void);

  void object_list_singleton(struct object_list *head, struct object *objp);
  void add_object_to_list(struct object_list *head, struct object *objp);
  struct vector3 *point_scalar(double val);

  struct vector3 *point_scalar(double val) {
    struct vector3 *vec;
    vec = (struct vector3 *)malloc(sizeof(struct vector3));
    if (!vec)
      return NULL;

    vec->x = val;
    vec->y = val;
    vec->z = val;
    return vec;
  }

  void object_list_singleton(struct object_list *head, struct object *objp) {
    objp->next = NULL;
    head->obj_tail = head->obj_head = objp;
  }

  void add_object_to_list(struct object_list *head, struct object *objp) {
    objp->next = NULL;
    head->obj_tail = head->obj_tail->next = objp;
  }

  int parse_dg(char *dynamic_geometry_filename) {
    FILE *fp=fopen(dynamic_geometry_filename,"r");
    if(!fp)
    {
      printf("Couldn't open file for reading\n");
      return 1;
    }
    yyin=fp;
    dg_parse = (struct dyngeom_parse_vars *)malloc(sizeof(struct dyngeom_parse_vars));
    memset(dg_parse, 0, sizeof(struct dyngeom_parse_vars));
    init_top_level_objs(dg_parse);
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

existing_object: var                                 {printf("existing_object\n");}

;

point: array_value                                    {printf("point\n");}
;

point_or_num: point                                   {printf("point_or_num\n");}
            | num_expr_only                           { $$ = point_scalar($1); }
;

list_range_specs:
          range_spec
        | list_range_specs ',' range_spec            {
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

range_spec: num_expr                                 { printf("range_spec\n"); }
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
        object_ref                                    { object_list_singleton(& $$, $1); }
      | list_objects object_ref                       { $$ = $1; add_object_to_list(& $$, $2); }
;

object_ref: existing_object_ref
          | object_def
;

existing_object_ref:
        new_object OBJECT existing_object
        start_object
          list_opt_object_cmds
        end_object                                   { $$ = (struct object *) $1->value; $$->object_type = POLY_OBJ; }
;

new_object_name: var                                 {printf("new_object_name\n");}
;

/* =================================================================== */
/* Instance definitions */

instance_def:
          INSTANTIATE                                { printf("INSTANTIATE\n"); dg_parse->current_object = dg_parse->root_instance; }
          meta_object_def                             {
                                                        printf("meta_object_def\n");
                                                        add_child_objects(dg_parse->root_instance, $3, $3);
                                                        dg_parse->current_object = dg_parse->root_object;
                                                      }
;

/* =================================================================== */
/* Object type definitions */

/*physical_object_def: object_def                      {printf("physical_object_def\n"); add_child_objects(dg_parse->root_object, $1, $1);}*/
physical_object_def: object_def                      {printf("physical_object_def\n");}
;

object_def: meta_object_def
          | polygon_list_def
;

new_object: var                                      {printf("new_object\n"); $$ = dg_start_object(dg_parse, $1);}
;

start_object: '{'                                    {printf("start_object\n");}
;

end_object: '}'                                      {printf("end_object\n"); dg_finish_object(dg_parse);}
;

list_opt_object_cmds:
        | list_opt_object_cmds opt_object_cmd
;

opt_object_cmd: transformation
;

transformation:
          TRANSLATE '=' point                        {printf("TRANSLATE\n");}
        | SCALE '=' point_or_num                     {printf("SCALE\n");}
        | ROTATE '=' point ',' num_expr              {printf("ROTATE\n");}
;

/* Object type: Polygons */
polygon_list_def:
          new_object_name POLYGON_LIST               {printf("POLYGON_LIST");}
          start_object
            vertex_list_cmd
            element_connection_cmd                    {
                                                        $<obj>$->object_type = POLY_OBJ;
                                                        //XXX: Need to actually create objects
                                                        /*$<obj>$ = mdl_new_polygon_list(*/
                                                        /*  parse_state, $1, $4.vertex_count, $4.vertex_head,*/
                                                        /*  $5.connection_count, $5.connection_head);*/
                                                      }
            list_opt_polygon_object_cmds
            list_opt_object_cmds
          '}'
                                                     { $$ = (struct object *) $<obj>6; }
;

vertex_list_cmd:
          VERTEX_LIST {printf("vertex_list_command\n");}
          '{' list_points '}'            
;

single_vertex: point                                 {printf("single_vertex\n");}
;

list_points: single_vertex                           {printf("list_points\n");}
           | list_points single_vertex
;

element_connection_cmd:
          ELEMENT_CONNECTIONS                        {printf("element_connection_cmd\n");}
          '{' list_element_connections '}'
;

list_element_connections:
          element_connection                         {printf("element_connection\n");}
        | list_element_connections
          element_connection                         {printf("element_connection\n");}
;

element_connection: array_value                      {printf("element_connection\n");}
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
          element_specifier                          {printf("element_specifier\n");}
;

element_specifier:
          incl_element_list_stmt
;

incl_element_list_stmt:
          INCLUDE_ELEMENTS '='
          '[' list_element_specs ']'                 {printf("incl_element_list_stmt\n");}
;

list_element_specs:
          element_spec
        | list_element_specs ',' element_spec
;

element_spec: num_expr                               {printf("element_spec\n");}
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

new_region: var                                      {printf("new_region\n");}
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

array_value: array_expr_only                         {printf("array_value\n");}
;

array_expr_only: '[' list_range_specs ']'            { printf("array_expr_only\n"); $$ = $2; }
;

num_expr: num_value                                   { printf("num_expr\n"); }
        | arith_expr
;

num_value: intOrReal                                  { printf("num_value\n"); }
         | existing_num_var                           { printf("num_value\n"); $$ = *(double *) $1->value; }
;

intOrReal: LLINTEGER                                  { printf("LLINTEGER\n"); $$ = $1; }
         | REAL                                       { printf("REAL\n"); }
;

num_expr_only: intOrReal
             | arith_expr
;


existing_num_var: var                                 { printf("existing_num_var\n"); }
;

arith_expr:
        '(' num_expr ')'                              { $$ = $2; }
      | '-' num_expr %prec UNARYMINUS                 { $$ = -$2; }
;


%%

void yyerror(char *s) {
  printf("%s\n", s);
}

int main(int argc, char *argv[])
{
  parse_dg(argv[1]);
}
