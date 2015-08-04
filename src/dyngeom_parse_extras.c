#include <stdlib.h>

#include "dyngeom_parse_extras.h"
#include "mcell_objects.h"
#include "mcell_structs.h"
#include "sym_table.h"

int init_top_level_objs(struct dyngeom_parse_vars *dg_parse_vars) {
  if ((dg_parse_vars->obj_sym_table = init_symtab(1024)) == NULL) {
    return 1;
  }

  struct sym_table *sym;
  if ((sym = store_sym(
      "WORLD_OBJ", OBJ, dg_parse_vars->obj_sym_table, NULL)) == NULL) {
    return 1;
  }

  dg_parse_vars->root_object = (struct object *)sym->value;
  dg_parse_vars->root_object->object_type = META_OBJ;
  if (!(dg_parse_vars->root_object->last_name = CHECKED_STRDUP_NODIE("", NULL))) {
    return 1;
  }

  if ((sym = store_sym(
      "WORLD_INSTANCE", OBJ, dg_parse_vars->obj_sym_table, NULL)) == NULL) {
    return 1;
  }

  dg_parse_vars->root_instance = (struct object *)sym->value;
  dg_parse_vars->root_instance->object_type = META_OBJ;
  if (!(dg_parse_vars->root_instance->last_name = CHECKED_STRDUP("", NULL))) {
    return 1;
  }

  dg_parse_vars->current_object = dg_parse_vars->root_instance;

  return 0;
}

struct sym_table *dg_start_object(
    struct dyngeom_parse_vars *dg_parse_vars,
    char *name) {
  // Create new fully qualified name.
  char *new_name;
  struct object_creation obj_creation;
  obj_creation.object_name_list = dg_parse_vars->object_name_list;
  obj_creation.object_name_list_end = dg_parse_vars->object_name_list_end;
  if ((new_name = push_object_name(&obj_creation, name)) == NULL) {
    free(name);
    return NULL;
  }
  dg_parse_vars->object_name_list = obj_creation.object_name_list;
  dg_parse_vars->object_name_list_end = obj_creation.object_name_list_end;

  // Create the symbol, if it doesn't exist yet.
  struct object *obj_ptr = make_new_object(
      dg_parse_vars->obj_sym_table, new_name, 0);
  if (obj_ptr == NULL) {
    if (name != new_name) {
      free(name);
    }
    free(new_name);
    return NULL;
  }

  struct sym_table *sym_ptr = obj_ptr->sym;
  obj_ptr->last_name = name;

  // Set parent object, make this object "current".
  obj_ptr->parent = dg_parse_vars->current_object;
  dg_parse_vars->current_object = obj_ptr;

  return sym_ptr;
}

void dg_finish_object(struct dyngeom_parse_vars *dg_parse_vars) {
  struct object_creation obj_creation;
  obj_creation.object_name_list_end = dg_parse_vars->object_name_list_end;

  pop_object_name(&obj_creation);
  dg_parse_vars->object_name_list_end = obj_creation.object_name_list_end;
  dg_parse_vars->current_object = dg_parse_vars->current_object->parent;
}

