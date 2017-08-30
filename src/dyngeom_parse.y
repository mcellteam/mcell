%{
  /*#define YYDEBUG 1*/
  #include <stdlib.h>
  #include <string.h>

  #include "mcell_misc.h"
  #include "mcell_objects.h"
  #include "logging.h"
  #include "strfunc.h"
  #include "dyngeom_parse_extras.h"

  typedef void *yyscan_t;
  #include "dyngeom_yacc.h"

  int dglex_init(yyscan_t *ptr_yy_globals) ;
  int dglex_destroy(yyscan_t yyscanner);
  void dgrestart(FILE *infile, yyscan_t scanner);
  int dglex(YYSTYPE *yylval, struct dyngeom_parse_vars *dg_parse, yyscan_t scanner);

  void dgerror(struct dyngeom_parse_vars *dg_parse, yyscan_t scanner, char const *str);

  void object_list_singleton(struct object_list *head, struct object *objp) {
    objp->next = NULL;
    head->obj_tail = head->obj_head = objp;
  }

  void add_object_to_list(struct object_list *head, struct object *objp) {
    objp->next = NULL;
    head->obj_tail = head->obj_tail->next = objp;
  }

  struct dyngeom_parse_vars * create_dg_parse(struct volume *state) {
    struct dyngeom_parse_vars *dg_parse = (struct dyngeom_parse_vars *)malloc(sizeof(struct dyngeom_parse_vars));
    memset(dg_parse, 0, sizeof(struct dyngeom_parse_vars));
    return dg_parse;
  }

  int parse_dg_init(
      struct dyngeom_parse_vars *dg_parse,
      char *dynamic_geometry_filename,
      struct volume *state) {
    dg_parse->include_stack_ptr = 0;
    dg_parse->line_num[0] = 0;
    setup_root_obj_inst(dg_parse, state);
    parse_dg(dg_parse, dynamic_geometry_filename);
    return 0;
  }

  int parse_dg(struct dyngeom_parse_vars *dg_parse,
               char *dynamic_geometry_filename) {
    int cur_stack = dg_parse->include_stack_ptr ++;
    dg_parse->line_num[cur_stack] = 1;
    dg_parse->include_filename[cur_stack] = dynamic_geometry_filename;

    // Open DG MDL
    no_printf("Opening file %s\n", dynamic_geometry_filename);
    FILE *infile=fopen(dynamic_geometry_filename,"r");
    if(!infile)
    {
      no_printf("Couldn't open file for reading\n");
      return 1;
    }

    yyscan_t scanner;
    dglex_init(&scanner);
    dgrestart(infile, scanner);

    // Parse DG MDL
    char const *prev_file;
    prev_file = dg_parse->curr_file;
    dg_parse->curr_file = dynamic_geometry_filename;
    /*extern int dgdebug;*/
    /*dgdebug = 1;*/
    dgparse(dg_parse, scanner);
    dg_parse->curr_file = prev_file;
    -- dg_parse->include_stack_ptr;

    /* Free leftover object names */
    struct name_list *nl;
    while (dg_parse->object_name_list != NULL)
    {
      nl = dg_parse->object_name_list->next;
      free(dg_parse->object_name_list);
      dg_parse->object_name_list = nl;
    }

    /* Clean up! */
    fclose(infile);
    dglex_destroy(scanner);

    return 0;
  }

#define CHECK(a)  do { if ((a) != 0) mcell_error_nodie("Parser fail: %s:%d\n", __FILE__, __LINE__); return 0; } while (0)

%}

/*%debug*/
%union {
  int tok;
  double dbl;
  char *str;
  long long llival;
  struct sym_entry *sym;
  struct vector3 *vec3;
  struct num_expr_list_head nlist;
  struct object *obj;
  struct object_list obj_list;
  struct region *reg;
}

%pure-parser

%lex-param {struct dyngeom_parse_vars *dg_parse}
%lex-param {yyscan_t scanner}
%parse-param {struct dyngeom_parse_vars *dg_parse}
%parse-param {yyscan_t scanner}
%name-prefix "dg"

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
%token       INCLUDE_FILE
%token       DEFINE_SURFACE_REGIONS
%token <str> STR_VALUE

%type <reg> new_region
%type <str> new_object_name
%type <dbl> intOrReal
%type <dbl> num_value
%type <sym> existing_num_var
%type <obj> existing_object_ref
%type <obj> object_def
%type <obj> object_ref
%type <obj_list> list_objects
%type <str> str_value
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
%type <str> str_expr
/*%type <sym> existing_str_var*/
%type <str> str_expr_only

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
        include_stmt
      | partition_def
      | physical_object_def
      | instance_def
;

/* =================================================================== */
/* Utility definitions */
str_value: STR_VALUE
;

var: VAR
;

existing_object: var                                 { $$ = dg_existing_object(dg_parse, $1); }

;

point: array_value                                   { /*no_printf("point\n");*/ }
;

point_or_num: point                                  { no_printf("point_or_num\n"); }
            | num_expr_only                          { }
;

list_range_specs:
          range_spec
        | list_range_specs ',' range_spec            { }
;

range_spec: num_expr                                 { }
        | '[' num_expr TO num_expr STEP num_expr ']' { }
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
        start_object                                 { dg_deep_copy_object(dg_parse, (struct object *) $1->value, (struct object *) $3->value); }
          list_opt_object_cmds
        end_object                                   { $$ = (struct object *) $1->value; }            
;

new_object_name: var                                 { no_printf("new_object_name\n"); }
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

physical_object_def: object_def                      { add_child_objects(dg_parse->root_object, $1, $1); no_printf("physical_object_def\n"); }
;

object_def: meta_object_def
          | polygon_list_def
;

new_object: var                                      { no_printf("new_object %s\n", $1); $$ = dg_start_object(dg_parse, $1); }
;

start_object: '{'                                    { no_printf("start_object\n"); }
;

end_object: '}'                                      { no_printf("end_object\n"); dg_finish_object(dg_parse); }
;

list_opt_object_cmds:
        | list_opt_object_cmds opt_object_cmd
;

opt_object_cmd: transformation
;

transformation:
          TRANSLATE '=' point                        { no_printf("TRANSLATE\n"); }
        | SCALE '=' point_or_num                     { no_printf("SCALE\n"); }
        | ROTATE '=' point ',' num_expr              { no_printf("ROTATE\n"); }
;

/* Object type: Polygons */
polygon_list_def:
          new_object_name POLYGON_LIST               { no_printf("POLYGON_LIST\n"); }
          start_object
            vertex_list_cmd
            element_connection_cmd                   { $<obj>$ = dg_new_polygon_list(dg_parse, $1); }
            list_opt_polygon_object_cmds
            list_opt_object_cmds
          '}'
                                                     { $$ = (struct object *) $<obj>7; dg_finish_object(dg_parse); }
;

vertex_list_cmd:
          VERTEX_LIST                                { no_printf("vertex_list_command\n"); }
          '{' list_points '}'            
;

single_vertex: point                                 { /*no_printf("single_vertex\n");*/ }
;

list_points: single_vertex                           { no_printf("list_points\n"); }
           | list_points single_vertex
;

element_connection_cmd:
          ELEMENT_CONNECTIONS                        { no_printf("element_connection_cmd\n"); }
          '{' list_element_connections '}'
;

list_element_connections:
          element_connection                         {}
        | list_element_connections
          element_connection                         {}
;

element_connection: array_value                      {}
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
          element_specifier                          { no_printf("element_specifier\n"); }
;

element_specifier:
          incl_element_list_stmt
;

incl_element_list_stmt:
          INCLUDE_ELEMENTS '='
          '[' list_element_specs ']'                 { no_printf("incl_element_list_stmt\n"); }
;

list_element_specs:
          element_spec
        | list_element_specs ',' element_spec
;

element_spec: num_expr                               { no_printf("element_spec\n"); }
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
          new_region '{'                             { dg_parse->current_region = $1; }
            element_specifier_list
          '}'                                        { dg_parse->current_region = NULL; }
;

                                                         
new_region: var                                      { dg_create_region(dg_parse->reg_sym_table, dg_parse->current_object, $1); no_printf("new_region\n"); }
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
/* Include files */
include_stmt: INCLUDE_FILE '=' str_expr               {
                                                          no_printf("include_stmt %s\n", $3);
                                                          char *include_path = find_include_file($3, dg_parse->curr_file);
                                                          if (include_path == NULL)
                                                          {
                                                            free($3);
                                                            return 1;
                                                          }
                                                          if (parse_dg(dg_parse, include_path))
                                                          {
                                                            free(include_path);
                                                            free($3);
                                                            return 1;
                                                          }
                                                          free(include_path);
                                                          free($3);
                                                      }
;

/* =================================================================== */
/* Expressions */

array_value: array_expr_only                         { }
;

array_expr_only: '[' list_range_specs ']'            { $$ = $2; }
;

num_expr: num_value                                  { }
        | arith_expr
;

num_value: intOrReal                                 { }
         | existing_num_var                          { $$ = *(double *) $1->value; }
;

intOrReal: LLINTEGER                                 { $$ = $1; }
         | REAL                                      { }
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

str_expr:
        str_expr_only                                { no_printf("str_expr_only\n"); }
;

str_expr_only:
        str_value                                     { no_printf("str_value %s\n", $$); $$ = strip_quotes($1); }
      | str_expr '&' str_expr                         { $$ = my_strcat($1, $3); }
;

%%

void dgerror(
    struct dyngeom_parse_vars *dg_parse,
    yyscan_t scanner,
    char const *str) {
  mcell_error("%s on line %d in %s\n", str, dg_parse->line_num[dg_parse->include_stack_ptr - 1], dg_parse->curr_file);
}

int main(int argc, char *argv[])
{
  /*struct dyngeom_parse_vars *dg_parse = create_dg_parse();*/
  /*parse_dg_init(dg_parse, argv[1]);*/
}
