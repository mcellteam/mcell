#include <stdlib.h>

#include "dyngeom_parse_extras.h"
#include "mcell_objects.h"
#include "mcell_structs.h"
#include "sym_table.h"

/***************************************************************************
 init_top_level_objs:
  Create the region and object symbol tables for the DG parser. Create the top
  level object and instance.

 In:  dg_parse_vars: state of dynamic geometry parser
 Out: Returns 1 on error, 0 on success.
***************************************************************************/
int init_top_level_objs(struct dyngeom_parse_vars *dg_parse_vars) {
  // Create the region symbol table
  if ((dg_parse_vars->reg_sym_table = init_symtab(1024)) == NULL) {
    return 1;
  }

  // Create the object symbol table
  if ((dg_parse_vars->obj_sym_table = init_symtab(1024)) == NULL) {
    return 1;
  }

  // Create the top-level/world/root object
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

  // Create the top-level/world/root instance
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

/***************************************************************************
 setup_root_obj_inst:
  This needs to be called before reading every new DG file.

 In:  dg_parse_vars: state of dynamic geometry parser
 Out: Returns 1 on error, 0 on success.
***************************************************************************/
int setup_root_obj_inst(struct dyngeom_parse_vars *dg_parse_vars) {
  struct sym_table *sym;

  sym = retrieve_sym("WORLD_OBJ", dg_parse_vars->obj_sym_table);
  dg_parse_vars->root_object = (struct object *)sym->value;

  sym = retrieve_sym("WORLD_INSTANCE", dg_parse_vars->obj_sym_table);
  dg_parse_vars->root_instance = (struct object *)sym->value;

  dg_parse_vars->current_object = dg_parse_vars->root_instance;
  dg_parse_vars->object_name_list = NULL;
  dg_parse_vars->object_name_list_end = NULL;

  return 0;
}

/***************************************************************************
 dg_start_object:
    Create a new object, adding it to the DG symbol table. The qualified name
    of the object will be built by adding to the object_name_list, and the
    object is made the "current_object" in the mdl parser state. Because of
    these side effects, it is vital to call dg_finish_object at the end of the
    scope of the object created here.

 In:  dg_parse_vars: state of dynamic geometry parser
      name: unqualified object name
 Out: the newly created object
 Note: This is similar to mdl_start_object.
***************************************************************************/
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
      dg_parse_vars->obj_sym_table, new_name, 1);
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

/***************************************************************************
 dg_finish_object:
    "Finishes" a new object, undoing the state changes that occurred when the
    object was "started". This means popping the name off of the object name
    stack, and resetting current_object to its value before this object was
    defined.

 In:  dg_parse_vars: state of dynamic geometry parser
 Out: none
 Note: This is similar to mdl_finish_object.
***************************************************************************/
void dg_finish_object(struct dyngeom_parse_vars *dg_parse_vars) {
  struct object_creation obj_creation;
  obj_creation.object_name_list_end = dg_parse_vars->object_name_list_end;

  pop_object_name(&obj_creation);
  dg_parse_vars->object_name_list_end = obj_creation.object_name_list_end;
  dg_parse_vars->current_object = dg_parse_vars->current_object->parent;
}

/***************************************************************************
 dg_create_region:
    Create a named region on an object.

 In:  reg_sym_table: the symbol table for all the regions
      objp: fully qualified object name
      name: name of the region to define
 Out: The new region.
 Note: This is similar to mdl_create_region.
***************************************************************************/
struct region *dg_create_region(
    struct sym_table_head *reg_sym_table,
    struct object *objp,
    char *name) {
  struct region *rp;
  struct region_list *rlp;
  if ((rp = dg_make_new_region(reg_sym_table, objp->sym->name, name)) == NULL)
    return NULL;
  if ((rlp = CHECKED_MALLOC_STRUCT(struct region_list, "region list")) ==
      NULL) {
    return NULL;
  }
  rp->flags = 0;
  rp->region_last_name = name;
  rp->manifold_flag = IS_MANIFOLD;
  rp->parent = objp;
  rlp->reg = rp;
  rlp->next = objp->regions;
  objp->regions = rlp;
  objp->num_regions++;
  return rp;
}

/***************************************************************************
 dg_make_new_region:
    Create a new region, adding it to the DG symbol table.

    full region names of REG type symbols stored in main symbol table have the
    form:
         metaobj.metaobj.poly,region_last_name

 In:  reg_sym_table: the symbol table for all the regions
      objp: fully qualified object name
      region_last_name: name of the region to define
 Out: The newly created region
 Note: This is similar to mdl_make_new_region.
***************************************************************************/
struct region *dg_make_new_region(
    struct sym_table_head *reg_sym_table,
    char *obj_name,
    char *region_last_name) {
  char *region_name;
  region_name = CHECKED_SPRINTF("%s,%s", obj_name, region_last_name);
  if (region_name == NULL) {
    return NULL;
  }

  struct sym_table *sym_ptr;
  if ((sym_ptr = retrieve_sym(region_name, reg_sym_table)) != NULL) {
    free(region_name);
    return NULL;
  }

  if ((sym_ptr = store_sym(region_name, REG, reg_sym_table, NULL)) == NULL) {
    free(region_name);
    return NULL;
  }

  free(region_name);
  return (struct region *)sym_ptr->value;
}
