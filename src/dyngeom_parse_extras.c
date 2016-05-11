#include <stdlib.h>
#include <string.h>

#include "dyngeom_parse_extras.h"
#include "mcell_structs.h"
#include "sym_table.h"
#include "logging.h"

/***************************************************************************
 setup_root_obj_inst:
  This needs to be called before reading every new DG file.

 In:  dg_parse_vars: state of dynamic geometry parser
 Out: none
***************************************************************************/
void setup_root_obj_inst(
    struct dyngeom_parse_vars *dg_parse_vars,
    struct volume *state) {
  struct sym_entry *sym;

  dg_parse_vars->obj_sym_table = state->obj_sym_table;
  dg_parse_vars->reg_sym_table = state->reg_sym_table;

  sym = retrieve_sym("WORLD_OBJ", state->obj_sym_table);
  dg_parse_vars->root_object = (struct object *)sym->value;

  sym = retrieve_sym("WORLD_INSTANCE", state->obj_sym_table);
  dg_parse_vars->root_instance = (struct object *)sym->value;

  dg_parse_vars->current_object = dg_parse_vars->root_object;
  dg_parse_vars->object_name_list = NULL;
  dg_parse_vars->object_name_list_end = NULL;
}

/***************************************************************************
 dg_start_object:
    Create a new object, adding it to the DG symbol table. The qualified name
    of the object will be built by adding to the object_name_list, and the
    object is made the "current_object" in the mdl parser state. Because of
    these side effects, it is vital to call dg_finish_object at the end of the
    scope of the object created here.

 In:  dg_parse_vars: state of dynamic geometry parser
      obj_name: unqualified object name
 Out: the newly created object or NULL on error
 Note: This is similar to mdl_start_object.
***************************************************************************/
struct sym_entry *dg_start_object(
    struct dyngeom_parse_vars *dg_parse_vars,
    char *obj_name) {
  // Create new fully qualified name.
  char *new_name;
  struct object_creation obj_creation;
  obj_creation.object_name_list = dg_parse_vars->object_name_list;
  obj_creation.object_name_list_end = dg_parse_vars->object_name_list_end;
  if ((new_name = push_object_name(&obj_creation, obj_name)) == NULL) {
    free(obj_name);
    return NULL;
  }
  dg_parse_vars->object_name_list = obj_creation.object_name_list;
  dg_parse_vars->object_name_list_end = obj_creation.object_name_list_end;

  // Create the symbol, if it doesn't exist yet.
  int error_code = 0;
  struct object *obj_ptr = make_new_object(
      dg_parse_vars, dg_parse_vars->obj_sym_table, new_name, &error_code);
  if (obj_ptr == NULL) {
    mcell_error("Object already defined: %s", new_name);
  }

  struct sym_entry *sym_ptr = obj_ptr->sym;
  obj_ptr->last_name = obj_name;

  // Set parent object, make this object "current".
  obj_ptr->parent = dg_parse_vars->current_object;
  dg_parse_vars->current_object = obj_ptr;

  return sym_ptr;
}

/***************************************************************************
 dg_start_object_simple:

 In:  dg_parse_vars: state of dynamic geometry parser
      obj_creation: information about object being created
      name: unqualified object name
 Out: the polygon object
 NOTE: This is similar to start_object.
 XXX: There are too many ways to create objects. This needs to be consolidated
 with dg_start_object.
***************************************************************************/
struct object *dg_start_object_simple(struct dyngeom_parse_vars *dg_parse_vars,
                                      struct object_creation *obj_creation,
                                      char *name) {
  // Create new fully qualified name.
  char *new_name;
  if ((new_name = push_object_name(obj_creation, name)) == NULL) {
    free(name);
    return NULL;
  }

  // Create the symbol, if it doesn't exist yet.
  int error_code = 0;
  struct object *obj_ptr = make_new_object(
      dg_parse_vars,
      dg_parse_vars->obj_sym_table,
      new_name,
      &error_code);
  if (obj_ptr == NULL) {
    mcell_error("Object already defined: %s", new_name);
  }

  obj_ptr->last_name = name;

  // Set parent object, make this object "current".
  obj_ptr->parent = obj_creation->current_object;

  return obj_ptr;
}

/***************************************************************************
 dg_new_polygon_list:

 In:  dg_parse_vars: state of dynamic geometry parser
      obj_ptr:
 Out: the polygon object
 Note: This is similar to mdl_new_polygon_list.
***************************************************************************/
struct object *dg_new_polygon_list(
    struct dyngeom_parse_vars *dg_parse_vars,
    char *obj_name) {
  struct object_creation obj_creation;
  obj_creation.object_name_list = dg_parse_vars->object_name_list;
  obj_creation.object_name_list_end = dg_parse_vars->object_name_list_end;
  obj_creation.current_object = dg_parse_vars->current_object;

  struct object *obj = dg_start_object_simple(
      dg_parse_vars, &obj_creation, obj_name);
  obj->object_type = POLY_OBJ;
  dg_create_region(dg_parse_vars->reg_sym_table, obj, "ALL");

  dg_parse_vars->object_name_list = obj_creation.object_name_list;
  dg_parse_vars->object_name_list_end = obj_creation.object_name_list_end;
  dg_parse_vars->current_object = obj;

  return obj;
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
      reg_name: reg_name of the region to define
 Out: The new region.
 Note: This is similar to mdl_create_region.
***************************************************************************/
struct region *dg_create_region(
    struct sym_table_head *reg_sym_table,
    struct object *objp,
    char *reg_name) {
  struct region *reg;
  struct region_list *reg_list;
  if ((reg = dg_make_new_region(reg_sym_table, objp->sym->name, reg_name)) ==
      NULL) {
    mcell_error("Region already defined: %s", reg_name);
  }
  if ((reg_list = CHECKED_MALLOC_STRUCT(struct region_list, "region list")) ==
      NULL) {
    return NULL;
  }
  reg->flags = 0;
  reg->region_last_name = reg_name;
  reg->manifold_flag = IS_MANIFOLD;
  reg->parent = objp;
  reg_list->reg = reg;
  reg_list->next = objp->regions;
  objp->regions = reg_list;
  objp->num_regions++;
  return reg;
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

  struct sym_entry *sym_ptr;
  if ((sym_ptr = retrieve_sym(region_name, reg_sym_table)) != NULL) {
    free(region_name);
    if (sym_ptr->count == 0) {
      sym_ptr->count = 1;
      return (struct region *)sym_ptr->value;
    }
    else {
      return NULL; 
    }
  }

  if ((sym_ptr = store_sym(region_name, REG, reg_sym_table, NULL)) == NULL) {
    free(region_name);
    return NULL;
  }

  free(region_name);
  return (struct region *)sym_ptr->value;
}

/***************************************************************************
 dg_copy_object_regions:
    Duplicate src_obj's regions on dst_obj.

 In:  dst_obj: destination object
      src_obj: object from which to copy
 Out: 0 on success, 1 on failure
 Note: This is similar to mdl_copy_object_regions.
***************************************************************************/
int dg_copy_object_regions(
    struct dyngeom_parse_vars *dg_parse_vars,
    struct object *dst_obj,
    struct object *src_obj) {
  struct region_list *src_rlp;
  struct region *dst_reg, *src_reg;

  /* Copy each region */
  for (src_rlp = src_obj->regions; src_rlp != NULL; src_rlp = src_rlp->next) {
    src_reg = src_rlp->reg;

    if ((dst_reg = dg_create_region(
        dg_parse_vars->reg_sym_table, dst_obj, src_reg->region_last_name)) == NULL)
      return 1;

    /* Copy over simple region attributes */
    dst_reg->surf_class = src_reg->surf_class;
    dst_reg->flags = src_reg->flags;
    dst_reg->area = src_reg->area;
    dst_reg->bbox = src_reg->bbox;
    dst_reg->manifold_flag = src_reg->manifold_flag;

  }
  return 0;
}

/***************************************************************************
 dg_deep_copy_object:
    Deep copy an object. The destination object should already be added to the
    symbol table, but should be otherwise unpopulated, as no effort is made to
    free any existing data contained in the object.

 In:  dst_obj: object into which to copy
      src_obj: object from which to copy
 Out: The newly created object
 Note: This is similar to mdl_deep_copy_object.
***************************************************************************/
int dg_deep_copy_object(
    struct dyngeom_parse_vars *dg_parse_vars,
    struct object *dst_obj,
    struct object *src_obj) {
  struct object *src_child;

  /* Copy over simple object attributes */
  dst_obj->object_type = src_obj->object_type;

  /* Copy over regions */
  if (dg_copy_object_regions(dg_parse_vars, dst_obj, src_obj))
    return 1;

  struct sym_entry *src_sym = retrieve_sym(src_obj->sym->name, dg_parse_vars->obj_sym_table);
  struct sym_entry *dst_sym = retrieve_sym(dst_obj->sym->name, dg_parse_vars->obj_sym_table);
  src_sym->value = src_obj;
  dst_sym->value = dst_obj;

  switch (dst_obj->object_type) {
  case META_OBJ:
    /* Copy children */
    for (src_child = src_obj->first_child; src_child != NULL;
         src_child = src_child->next) {
      struct object *dst_child;
      char *child_obj_name =
          CHECKED_SPRINTF("%s.%s", dst_obj->sym->name, src_child->last_name);
      if (child_obj_name == NULL)
        return 1;

      /* Create child object */
      int error_code = 0;
      if ((dst_child = make_new_object(dg_parse_vars, dg_parse_vars->obj_sym_table, child_obj_name, &error_code)) == NULL) {
        free(child_obj_name);
        return 1;
      }
      free(child_obj_name);

      /* Copy in last name */
      dst_child->last_name = strdup(src_child->last_name);
      free(src_child->last_name);
      if (dst_child->last_name == NULL)
        return 1;

      /* Recursively copy object and its children */
      if (dg_deep_copy_object(dg_parse_vars, dst_child, src_child))
        return 1;
      dst_child->parent = dst_obj;
      dst_child->next = NULL;
      add_child_objects(dst_obj, dst_child, dst_child);
    }
    break;

  case POLY_OBJ:
  case REL_SITE_OBJ:
  case BOX_OBJ:
  case VOXEL_OBJ:
  default:
    return 0;
  }

  return 0;
}

/***************************************************************************
 dg_existing_object:
    Find an existing object.

 In:  parse_state: parser state
      name: fully qualified object name
 Out: the object
 Note: This is similar to mdl_existing_object.
***************************************************************************/
struct sym_entry *dg_existing_object(struct dyngeom_parse_vars *dg_parse_vars, char *name) {
  // Check to see if it is one of the objects that will be added in
  // the future via a dynamic geometry event.
  return retrieve_sym(name, dg_parse_vars->obj_sym_table);
}

char *find_include_file(char const *path, char const *cur_path) {
  char *candidate = NULL;
  if (path[0] == '/')
    candidate = strdup(path);
  else {
    char *last_slash = strrchr(cur_path, '/');
#ifdef _WIN32
    char *last_bslash = strrchr(cur_path, '\\');
    if (last_bslash > last_slash)
      last_slash = last_bslash;
#endif
    if (last_slash == NULL)
      candidate = strdup(path);
    else
      candidate = CHECKED_SPRINTF("%.*s/%s", (int)(last_slash - cur_path),
                                  cur_path, path);
  }

  return candidate;
}
