/******************************************************************************
 *
 * Copyright (C) 2006-2017 by
 * The Salk Institute for Biological Studies and
 * Pittsburgh Supercomputing Center, Carnegie Mellon University
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301,
 * USA.
 *
******************************************************************************/

#include <assert.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#include "config.h"

#include "logging.h"
#include "sym_table.h"
#include "mcell_species.h"
#include "mcell_release.h"
#include "mcell_init.h"
#include "mcell_objects.h"
#include "dyngeom_parse_extras.h"
#include "mem_util.h"

/* static helper functions */
static int is_region_degenerate(struct region *reg_ptr);

/*************************************************************************
 mcell_create_instance_object:
  Create a new instance object.

 In: state:    the simulation state
               object pointer to store created meta object
 Out: 0 on success; any other integer value is a failure.
      A mesh is created.
*************************************************************************/
MCELL_STATUS
mcell_create_instance_object(MCELL_STATE *state, char *name,
                             struct object **new_obj) {
  // Create the symbol, if it doesn't exist yet.
  int error_code = 0;
  struct object *obj_ptr = make_new_object(
      state->dg_parse, // Need to test that dg_parse actually works here
      state->obj_sym_table,
      name,
      &error_code);
  /*struct object *obj_ptr = make_new_object(state, name, &error_code);*/
  if (obj_ptr == NULL) {
    return MCELL_FAIL;
  }
  obj_ptr->last_name = name;
  obj_ptr->object_type = META_OBJ;

  // instantiate object
  obj_ptr->parent = state->root_instance;
  add_child_objects(state->root_instance, obj_ptr, obj_ptr);

  *new_obj = obj_ptr;

  return MCELL_SUCCESS;
}

MCELL_STATUS mcell_create_periodic_box(
    struct volume *state,
    char *box_name,
    struct vector3 *llf,
    struct vector3 *urb) {
  
  double sf = 100.0; // scaling factor
  struct vector3 llf_sm = { .x=llf->x/sf,
                            .y=llf->y/sf,
                            .z=llf->z/sf,
  };
  struct vector3 urb_sm = { .x=urb->x/sf,
                            .y=urb->y/sf,
                            .z=urb->z/sf,
  };

  struct vertex_list *verts = mcell_add_to_vertex_list(
      urb_sm.x, urb_sm.y, llf_sm.z, NULL);
  verts = mcell_add_to_vertex_list(urb_sm.x, llf_sm.y, llf_sm.z, verts);
  verts = mcell_add_to_vertex_list(llf_sm.x, llf_sm.y, llf_sm.z, verts);
  verts = mcell_add_to_vertex_list(llf_sm.x, urb_sm.y, llf_sm.z, verts);
  verts = mcell_add_to_vertex_list(urb_sm.x, urb_sm.y, urb_sm.z, verts);
  verts = mcell_add_to_vertex_list(urb_sm.x, llf_sm.y, urb_sm.z, verts);
  verts = mcell_add_to_vertex_list(llf_sm.x, llf_sm.y, urb_sm.z, verts);
  verts = mcell_add_to_vertex_list(llf_sm.x, urb_sm.y, urb_sm.z, verts);

  struct element_connection_list *elems =
      mcell_add_to_connection_list(1, 2, 3, NULL);
  elems = mcell_add_to_connection_list(7, 6, 5, elems);
  elems = mcell_add_to_connection_list(0, 4, 5, elems);
  elems = mcell_add_to_connection_list(1, 5, 6, elems);
  elems = mcell_add_to_connection_list(6, 7, 3, elems);
  elems = mcell_add_to_connection_list(0, 3, 7, elems);
  elems = mcell_add_to_connection_list(0, 1, 3, elems);
  elems = mcell_add_to_connection_list(4, 7, 5, elems);
  elems = mcell_add_to_connection_list(1, 0, 5, elems);
  elems = mcell_add_to_connection_list(2, 1, 6, elems);
  elems = mcell_add_to_connection_list(2, 6, 3, elems);
  elems = mcell_add_to_connection_list(4, 0, 7, elems);

  struct subdivided_box *b = CHECKED_MALLOC_STRUCT(
      struct subdivided_box, "subdivided box");
  if (b == NULL)
    return MCELL_FAIL;

  b->nx = b->ny = b->nz = 2;
  if ((b->x = CHECKED_MALLOC_ARRAY(
      double, b->nx, "subdivided box X partitions")) == NULL) {
    free(b);
    return MCELL_FAIL;
  }
  if ((b->y = CHECKED_MALLOC_ARRAY(
      double, b->ny, "subdivided box Y partitions")) == NULL) {
    free(b);
    return MCELL_FAIL;
  }
  if ((b->z = CHECKED_MALLOC_ARRAY(
      double, b->nz, "subdivided box Z partitions")) == NULL) {
    free(b);
    return MCELL_FAIL;
  }

  b->x[0] = llf->x;
  b->x[1] = urb->x;
  b->y[0] = llf->y;
  b->y[1] = urb->y;
  b->z[0] = llf->z;
  b->z[1] = urb->z;
  
  struct object *meta_box = NULL;
  mcell_create_instance_object(state, "PERIODIC_META_BOX", &meta_box);

  struct poly_object polygon = { box_name, verts, 8, elems, 12 };
  struct object *new_mesh = NULL;
  mcell_create_poly_object(state, meta_box, &polygon, &new_mesh);
  new_mesh->periodic_x = true;
  new_mesh->periodic_y = true;
  new_mesh->periodic_z = true;
  new_mesh->object_type = BOX_OBJ;
  state->periodic_box_obj = new_mesh;
  struct polygon_object* p = (struct polygon_object*)(new_mesh->contents);
  p->sb = b;

  return MCELL_SUCCESS;
}

/*************************************************************************
 mcell_create_poly_object:
  Create a new polygon object.

 In: state:    the simulation state
     poly_obj: all the information needed to create the polygon object (name,
               vertices, connections)
 Out: 0 on success; any other integer value is a failure.
      A mesh is created.
*************************************************************************/
MCELL_STATUS
mcell_create_poly_object(MCELL_STATE *state, struct object *parent,
                         struct poly_object *poly_obj,
                         struct object **new_obj) {
  // create qualified object name
  char *qualified_name =
      CHECKED_SPRINTF("%s.%s", parent->sym->name, poly_obj->obj_name);

  // Create the symbol, if it doesn't exist yet.
  int error_code = 0;
  struct object *obj_ptr = make_new_object(
      state->dg_parse, // Need to test that dg_parse actually works here
      state->obj_sym_table,
      qualified_name,
      &error_code);
  /*struct object *obj_ptr = make_new_object(state, qualified_name, &error_code);*/
  if (obj_ptr == NULL) {
    free(qualified_name);
    return MCELL_FAIL;
  }
  obj_ptr->last_name = qualified_name;

  // Create the actual polygon object
  new_polygon_list(state, obj_ptr, poly_obj->num_vert, poly_obj->vertices,
                   poly_obj->num_conn, poly_obj->connections);

  // Do some clean-up.
  remove_gaps_from_regions(obj_ptr);
  if (check_degenerate_polygon_list(obj_ptr)) {
    return MCELL_FAIL;
  }

  // Set the parent of the object to be the root object. Not reciprocal until
  // add_child_objects is called.
  obj_ptr->parent = parent;
  add_child_objects(parent, obj_ptr, obj_ptr);

  *new_obj = obj_ptr;

  return MCELL_SUCCESS;
}

/**************************************************************************
 new_polygon_list:
    Create a new polygon list object.

 In: state: the simulation state
     obj_ptr: contains information about the object (name, etc)
     n_vertices: count of vertices
     vertices: list of vertices
     n_connections: count of walls
     connections: list of walls
 Out: polygon object, or NULL if there was an error
 NOTE: This is similar to mdl_new_polygon_list
**************************************************************************/
struct polygon_object *
new_polygon_list(MCELL_STATE *state, struct object *obj_ptr, int n_vertices,
                 struct vertex_list *vertices, int n_connections,
                 struct element_connection_list *connections) {

  struct polygon_object *poly_obj_ptr =
      allocate_polygon_object("polygon list object");
  if (poly_obj_ptr == NULL) {
    goto failure;
  }

  obj_ptr->object_type = POLY_OBJ;
  obj_ptr->contents = poly_obj_ptr;

  poly_obj_ptr->n_walls = n_connections;
  poly_obj_ptr->n_verts = n_vertices;

  // Allocate and initialize removed sides bitmask
  poly_obj_ptr->side_removed = new_bit_array(poly_obj_ptr->n_walls);
  if (poly_obj_ptr->side_removed == NULL) {
    goto failure;
  }
  set_all_bits(poly_obj_ptr->side_removed, 0);

  // Keep temporarily information about vertices in the form of
  // "parsed_vertices"
  poly_obj_ptr->parsed_vertices = vertices;

  // Copy in vertices and normals
  struct vertex_list *vert_list = poly_obj_ptr->parsed_vertices;
  for (int i = 0; i < poly_obj_ptr->n_verts; i++) {
    // Rescale vertices coordinates
    vert_list->vertex->x *= state->r_length_unit;
    vert_list->vertex->y *= state->r_length_unit;
    vert_list->vertex->z *= state->r_length_unit;
    vert_list = vert_list->next;
  }

  // Allocate wall elements
  struct element_data *elem_data_ptr = NULL;
  if ((elem_data_ptr =
           CHECKED_MALLOC_ARRAY(struct element_data, poly_obj_ptr->n_walls,
                                "polygon list object walls")) == NULL) {
    goto failure;
  }
  poly_obj_ptr->element = elem_data_ptr;

  // Copy in wall elements
  for (int i = 0; i < poly_obj_ptr->n_walls; i++) {
    if (connections->n_verts != 3) {
      // mdlerror(parse_state, "All polygons must have three vertices.");
      goto failure;
    }

    struct element_connection_list *elem_conn_list_temp = connections;
    memcpy(elem_data_ptr[i].vertex_index, connections->indices,
           3 * sizeof(int));
    connections = connections->next;
    free(elem_conn_list_temp->indices);
    free(elem_conn_list_temp);
  }

  // Create object default region on polygon list object:
  struct region *reg_ptr = NULL;
  if ((reg_ptr = mcell_create_region(state, obj_ptr, "ALL")) == NULL) {
    goto failure;
  }
  if ((reg_ptr->element_list_head =
           new_element_list(0, poly_obj_ptr->n_walls - 1)) == NULL) {
    goto failure;
  }

  obj_ptr->n_walls = poly_obj_ptr->n_walls;
  obj_ptr->n_verts = poly_obj_ptr->n_verts;
  if (normalize_elements(reg_ptr, 0)) {
    // mdlerror_fmt(parse_state,
    //             "Error setting up elements in default 'ALL' region in the "
    //             "polygon object '%s'.", sym->name);
    goto failure;
  }

  return poly_obj_ptr;

failure:
  free_connection_list(connections);
  free_vertex_list(vertices);
  if (poly_obj_ptr) {
    if (poly_obj_ptr->element) {
      free(poly_obj_ptr->element);
    }
    if (poly_obj_ptr->side_removed) {
      free_bit_array(poly_obj_ptr->side_removed);
    }
    free(poly_obj_ptr);
  }
  return NULL;
}

/*************************************************************************
 make_new_object:
    Create a new object, adding it to the global symbol table.

 In:  state: system state
      obj_name: fully qualified object name
 Out: the newly created object
*************************************************************************/
struct object *make_new_object(
    struct dyngeom_parse_vars *dg_parse,
    struct sym_table_head *obj_sym_table,
    char *obj_name,
    int *error_code) {

  struct sym_entry *symbol = retrieve_sym(obj_name, obj_sym_table);
  if (symbol != NULL) {
    if (symbol->count == 0) {
      symbol->count = 1;
      return (struct object *)symbol->value;
    }
    else {
      *error_code = 1;
      return NULL; 
    }
  }

  if ((symbol = store_sym(obj_name, OBJ, obj_sym_table, NULL)) == NULL) {
    *error_code = 2;
    return NULL;
  }

  *error_code = 0;
  return (struct object *)symbol->value;
}

/*************************************************************************
 push_object_name:
    Append a name component to the name list. For instance, consider this
    portion of an MDL:

    INSTANTIATE Scene OBJECT
    {
      MetaBox OBJECT {
        MyBox1 OBJECT MyBox{}
        MyBox2 OBJECT MyBox{TRANSLATE = [1, 0, 0]}
      }
    }

    When parsing the line beginning with My_Box1, return Scene.MetaBox.MyBox1


 In:  obj_creation: information about object being created
      name: new name component
 Out: object name stack is updated, returns new qualified object name
*************************************************************************/
char *push_object_name(struct object_creation *obj_creation, char *name) {

  // Initialize object name list
  if (obj_creation->object_name_list == NULL) {
    obj_creation->object_name_list =
        CHECKED_MALLOC_STRUCT(struct name_list, "object name stack");
    if (obj_creation->object_name_list == NULL) {
      return NULL;
    }

    obj_creation->object_name_list->name = NULL;
    obj_creation->object_name_list->prev = NULL;
    obj_creation->object_name_list->next = NULL;
    obj_creation->object_name_list_end = obj_creation->object_name_list;
  }

  /* If the last element is available, just use it.  This typically occurs only
   * for the first item in the list. */
  if (obj_creation->object_name_list_end->name == NULL) {
    obj_creation->object_name_list_end->name = name;
    return obj_creation->object_name_list_end->name;
  }

  // If we've run out of name list components, create a new one
  struct name_list *name_list_ptr;
  if (obj_creation->object_name_list_end->next == NULL) {
    name_list_ptr =
        CHECKED_MALLOC_STRUCT(struct name_list, "object name stack");
    if (name_list_ptr == NULL) {
      return NULL;
    }

    name_list_ptr->next = NULL;
    name_list_ptr->prev = obj_creation->object_name_list_end;
    obj_creation->object_name_list_end->next = name_list_ptr;
  } else {
    name_list_ptr = obj_creation->object_name_list_end->next;
  }

  // Create new name
  name_list_ptr->name =
      CHECKED_SPRINTF("%s.%s", obj_creation->object_name_list_end->name, name);
  if (name_list_ptr->name == NULL) {
    return NULL;
  }

  obj_creation->object_name_list_end = name_list_ptr;

  return name_list_ptr->name;
}

/*************************************************************************
 pop_object_name:
    Remove the trailing name component from the name list. It is expected that
    ownership of the name pointer has passed to someone else, or that the
    pointer has been freed already.

 In:  obj_creation: information about object being created
 Out: object name stack is updated
*************************************************************************/
void pop_object_name(struct object_creation *obj_creation) {
  if (obj_creation->object_name_list_end->name != NULL) {
    free(obj_creation->object_name_list_end->name);
  }
  if (obj_creation->object_name_list_end->prev != NULL) {
    obj_creation->object_name_list_end =
        obj_creation->object_name_list_end->prev;
  } else {
    obj_creation->object_name_list_end->name = NULL;
  }
}

/**************************************************************************
 add_child_objects:
    Adds children to a meta-object, aggregating counts of walls and vertices
    from the children into the specified parent.  The children should already
    have their parent pointers set.  (This must happen earlier so that we can
    resolve and validate region references in certain cases before this
    function is called.)

 In: parent: the parent object
     child_head: pointer to head of child list
     child_tail: pointer to tail of child list
 Out: parent object is updated; child_tail->next pointer is set to NULL
**************************************************************************/
void add_child_objects(struct object *parent, struct object *child_head,
                       struct object *child_tail) {

  if (parent->first_child == NULL) {
    parent->first_child = child_head;
  }
  if (parent->last_child != NULL) {
    parent->last_child->next = child_head;
  }
  parent->last_child = child_tail;
  child_tail->next = NULL;

  while (child_head != NULL) {
    assert(child_head->parent == parent);
    parent->n_walls += child_head->n_walls;
    parent->n_walls_actual += child_head->n_walls_actual;
    parent->n_verts += child_head->n_verts;
    child_head = child_head->next;
  }
}

/*************************************************************************
 common_ancestor:
    Find the nearest common ancestor of two objects

 In: a, b: objects
 Out: their common ancestor in the object tree, or NULL if none exists
*************************************************************************/
struct object *common_ancestor(struct object *a, struct object *b) {
  struct object *pa, *pb;

  for (pa = (a->object_type == META_OBJ) ? a : a->parent; pa != NULL;
       pa = pa->parent) {
    for (pb = (b->object_type == META_OBJ) ? b : b->parent; pb != NULL;
         pb = pb->parent) {
      if (pa == pb) {
        return pa;
      }
    }
  }

  return NULL;
}

/*************************************************************************
 remove_gaps_from_regions:
    Clean up the regions on an object, eliminating any removed walls.

 In: obj_ptr: an object with regions
 Out: Any walls that have been removed from the object are removed from every
      region on that object.
*************************************************************************/
void remove_gaps_from_regions(struct object *obj_ptr) {
  if (obj_ptr->object_type != BOX_OBJ && obj_ptr->object_type != POLY_OBJ) {
    return;
  }

  struct polygon_object *poly_obj_ptr =
      (struct polygon_object *)obj_ptr->contents;
  struct region_list *reg_list;
  for (reg_list = obj_ptr->regions; reg_list != NULL;
       reg_list = reg_list->next) {
    no_printf("Checking region %s\n", reg_list->reg->sym->name);
    if (strcmp(reg_list->reg->region_last_name, "REMOVED") == 0) {
      no_printf("Found a REMOVED region\n");
      reg_list->reg->surf_class = NULL;
      bit_operation(poly_obj_ptr->side_removed, reg_list->reg->membership, '+');
      set_all_bits(reg_list->reg->membership, 0);
    }
  }

  int missing = 0;
  for (int n_side = 0; n_side < poly_obj_ptr->side_removed->nbits; ++n_side) {
    if (get_bit(poly_obj_ptr->side_removed, n_side)) {
      missing++;
    }
  }

  obj_ptr->n_walls_actual = poly_obj_ptr->n_walls - missing;

  for (reg_list = obj_ptr->regions; reg_list != NULL;
       reg_list = reg_list->next) {
    bit_operation(reg_list->reg->membership, poly_obj_ptr->side_removed, '-');
  }
}

/**************************************************************************
 check_degenerate_polygon_list:
    Check a box or polygon list object for degeneracy.

 In: obj_ptr: the object to validate
 Out: 0 if valid, 1 if invalid
**************************************************************************/
int check_degenerate_polygon_list(struct object *obj_ptr) {
  // Check for a degenerate (empty) object and regions.
  struct region_list *reg_list;
  for (reg_list = obj_ptr->regions; reg_list != NULL;
       reg_list = reg_list->next) {
    if (!is_region_degenerate(reg_list->reg)) {
      continue;
    }
    if (strcmp(reg_list->reg->region_last_name, "ALL") == 0) {
      return 1;
    } else if (strcmp(reg_list->reg->region_last_name, "REMOVED") != 0) {
      return 1;
    }
  }
  return 0;
}

/**************************************************************************
 allocate_polygon_object:
    Allocate a polygon object.

 In: desc: object type description
 Out: polygon object on success, NULL on failure
**************************************************************************/
struct polygon_object *allocate_polygon_object(char const *desc) {
  struct polygon_object *poly_obj_ptr;
  if ((poly_obj_ptr = CHECKED_MALLOC_STRUCT(struct polygon_object, desc)) ==
      NULL) {
    return NULL;
  }
  poly_obj_ptr->n_verts = 0;
  poly_obj_ptr->parsed_vertices = NULL;
  poly_obj_ptr->n_walls = 0;
  poly_obj_ptr->element = NULL;
  poly_obj_ptr->sb = NULL;
  poly_obj_ptr->side_removed = NULL;
  poly_obj_ptr->references = 0;
  return poly_obj_ptr;
}

/**************************************************************************
 new_element_list:
    Create a new element list for a region description.

 In: begin: starting side number for this element
     end: ending side number for this element
 Out: element list, or NULL if allocation fails
**************************************************************************/
struct element_list *new_element_list(unsigned int begin, unsigned int end) {

  struct element_list *elem_list =
      CHECKED_MALLOC_STRUCT(struct element_list, "region element");
  if (elem_list == NULL) {
    return NULL;
  }
  elem_list->special = NULL;
  elem_list->next = NULL;
  elem_list->begin = begin;
  elem_list->end = end;
  return elem_list;
}

/**************************************************************************
 free_vertex_list:
    Free a vertex list.

 In: vert_list: vertex to free
 Out: list is freed
**************************************************************************/
void free_vertex_list(struct vertex_list *vert_list) {
  while (vert_list) {
    struct vertex_list *next = vert_list->next;
    free(vert_list->vertex);
    free(vert_list);
    vert_list = next;
  }
}

/**************************************************************************
 free_connection_list:
    Free a connection list.

 In: elem_conn_list: connection list to free
 Out: list is freed
**************************************************************************/
void free_connection_list(struct element_connection_list *elem_conn_list) {
  while (elem_conn_list) {
    struct element_connection_list *next = elem_conn_list->next;
    free(elem_conn_list->indices);
    free(elem_conn_list);
    elem_conn_list = next;
  }
}

/*************************************************************************
 normalize_elements:
    Prepare a region description for use in the simulation by creating a
    membership bitmask on the region object.

 In: reg: a region
     existing: a flag indicating whether the region is being modified or
               created
 Out: returns 1 on failure, 0 on success.  Lists of element specifiers
      are converted into bitmasks that specify whether a wall is in or
      not in that region.  This also handles combinations of regions.
 NOTE: This is similar to mdl_normalize_elements
*************************************************************************/
int normalize_elements(struct region *reg, int existing) {
  struct bit_array *temp = NULL;
  char op;
  unsigned int num_elems;

  if (reg->element_list_head == NULL) {
    return 0;
  }

  struct polygon_object *poly_obj = NULL;
  if (reg->parent->object_type == BOX_OBJ) {
    poly_obj = (struct polygon_object *)reg->parent->contents;
    num_elems = count_cuboid_elements(poly_obj->sb);
  } else {
    num_elems = reg->parent->n_walls;
  }

  struct bit_array *elem_array;
  if (reg->membership == NULL) {
    elem_array = new_bit_array(num_elems);
    if (elem_array == NULL) {
      // mcell_allocfailed("Failed to allocate a region membership bitmask.");
      return 1;
    }
    reg->membership = elem_array;
  } else {
    elem_array = reg->membership;
  }

  if (reg->element_list_head->special == NULL) {
    set_all_bits(elem_array, 0);
  }
  // Special flag for exclusion
  else if ((void *)reg->element_list_head->special ==
           (void *)reg->element_list_head) {
    set_all_bits(elem_array, 1);
  } else {
    if (reg->element_list_head->special->exclude) {
      set_all_bits(elem_array, 1);
    } else {
      set_all_bits(elem_array, 0);
    }
  }

  int i = 0;
  struct element_list *elem_list;
  for (elem_list = reg->element_list_head; elem_list != NULL;
       elem_list = elem_list->next) {
    if (reg->parent->object_type == BOX_OBJ) {
      assert(poly_obj != NULL);
      i = elem_list->begin;
      switch (i) {
      case X_NEG:
        elem_list->begin = 0;
        elem_list->end =
            2 * (poly_obj->sb->ny - 1) * (poly_obj->sb->nz - 1) - 1;
        break;
      case X_POS:
        elem_list->begin = 2 * (poly_obj->sb->ny - 1) * (poly_obj->sb->nz - 1);
        elem_list->end =
            4 * (poly_obj->sb->ny - 1) * (poly_obj->sb->nz - 1) - 1;
        break;
      case Y_NEG:
        elem_list->begin = 4 * (poly_obj->sb->ny - 1) * (poly_obj->sb->nz - 1);
        elem_list->end = elem_list->begin +
                         2 * (poly_obj->sb->nx - 1) * (poly_obj->sb->nz - 1) -
                         1;
        break;
      case Y_POS:
        elem_list->begin = 4 * (poly_obj->sb->ny - 1) * (poly_obj->sb->nz - 1) +
                           2 * (poly_obj->sb->nx - 1) * (poly_obj->sb->nz - 1);
        elem_list->end = elem_list->begin +
                         2 * (poly_obj->sb->nx - 1) * (poly_obj->sb->nz - 1) -
                         1;
        break;
      case Z_NEG:
        elem_list->begin = 4 * (poly_obj->sb->ny - 1) * (poly_obj->sb->nz - 1) +
                           4 * (poly_obj->sb->nx - 1) * (poly_obj->sb->nz - 1);
        elem_list->end = elem_list->begin +
                         2 * (poly_obj->sb->nx - 1) * (poly_obj->sb->ny - 1) -
                         1;
        break;
      case Z_POS:
        elem_list->end = num_elems - 1;
        elem_list->begin = elem_list->end + 1 -
                           2 * (poly_obj->sb->nx - 1) * (poly_obj->sb->ny - 1);
        break;
      case ALL_SIDES:
        elem_list->begin = 0;
        elem_list->end = num_elems - 1;
        break;
      default:
        UNHANDLED_CASE(i);
        /*return 1;*/
      }
    } else if (elem_list->begin >= (u_int)num_elems ||
               elem_list->end >= (u_int)num_elems) {
      // mdlerror_fmt(parse_state,
      //             "Region element specifier refers to sides %u...%u, but
      // polygon has only %u sides.",
      //             elem_list->begin,
      //             elem_list->end,
      //             num_elems);
      return 1;
    }

    if (elem_list->special == NULL) {
      set_bit_range(elem_array, elem_list->begin, elem_list->end, 1);
    } else if ((void *)elem_list->special == (void *)elem_list) {
      set_bit_range(elem_array, elem_list->begin, elem_list->end, 0);
    } else {
      if (elem_list->special->referent != NULL) {
        if (elem_list->special->referent->membership == NULL) {
          if (elem_list->special->referent->element_list_head != NULL) {
            i = normalize_elements(elem_list->special->referent, existing);
            if (i) {
              free_bit_array(temp);
              return i;
            }
          }
        }
        if (elem_list->special->referent->membership != NULL) {
          // What does it mean for the membership array to have length zero?
          if (elem_list->special->referent->membership->nbits == 0) {
            if (elem_list->special->exclude) {
              set_all_bits(elem_array, 0);
            } else {
              set_all_bits(elem_array, 1);
            }
          } else {
            if (elem_list->special->exclude) {
              op = '-';
            } else {
              op = '+';
            }
            bit_operation(elem_array, elem_list->special->referent->membership,
                          op);
          }
        }
      } else {
        int ii;
        if (temp == NULL) {
          temp = new_bit_array(num_elems);
          if (temp == NULL) {
            // mcell_allocfailed("Failed to allocate a region membership
            // bitmask.");
            return 1;
          }
        }
        if (poly_obj == NULL) {
          // mcell_internal_error("Attempt to create a PATCH on a
          // POLYGON_LIST.");
          free_bit_array(temp);
          return 1;
        }
        if (existing) {
          // mcell_internal_error("Attempt to create a PATCH on an already
          // triangulated BOX.");
          free_bit_array(temp);
          return 1;
        }
        if (elem_list->special->exclude) {
          op = '-';
        } else {
          op = '+';
        }

        ii = cuboid_patch_to_bits(poly_obj->sb, &(elem_list->special->corner1),
                                  &(elem_list->special->corner2), temp);
        if (ii) {
          // Something wrong with patch.
          free_bit_array(temp);
          return 1;
        }
        bit_operation(elem_array, temp, op);
      }
    }
  }

  if (temp != NULL) {
    free_bit_array(temp);
  }

  if (existing) {
    bit_operation(
        elem_array,
        ((struct polygon_object *)reg->parent->contents)->side_removed, '-');
  }

  while (reg->element_list_head) {
    struct element_list *next = reg->element_list_head->next;
    if (reg->element_list_head->special) {
      if (reg->element_list_head->special !=
          (struct element_special *)reg->element_list_head) {
        free(reg->element_list_head->special);
      }
    }
    free(reg->element_list_head);
    reg->element_list_head = next;
  }
  return 0;
}

/*************************************************************************
 count_cuboid_elements:
    Trivial utility function that counts # walls in a box

 In: sb: a subdivided box
 Out: the number of walls in the box
*************************************************************************/
int count_cuboid_elements(struct subdivided_box *sb) {
  return 4 * ((sb->nx - 1) * (sb->ny - 1) + (sb->nx - 1) * (sb->nz - 1) +
              (sb->ny - 1) * (sb->nz - 1));
}

/*************************************************************************
 check_patch:

 In: b: a subdivided box (array of subdivision locations along each axis)
     p1: 3D vector that is one corner of the patch
     p2: 3D vector that is the other corner of the patch
     egd: the effector grid density, which limits how fine divisions can be
 Out: returns 0 if the patch is broken, or a bitmask saying which coordinates
      need to be subdivided in order to represent the new patch, if the
      spacings are okay.
 Note: the two corners of the patch must be aligned with a Cartesian plane
      (i.e. it is a planar patch), and must be on the surface of the subdivided
      box.  Furthermore, the coordinates in the first corner must be smaller or
      the same size as the coordinates in the second corner.
*************************************************************************/
int check_patch(struct subdivided_box *b, struct vector3 *p1,
                struct vector3 *p2, double egd) {
  int i = 0;
  int nbits = 0;
  int j;
  const double minspacing = sqrt(2.0 / egd);
  double d;

  if (distinguishable(p1->x, p2->x, EPS_C)) {
    i |= BRANCH_X;
    nbits++;
  }
  if (distinguishable(p1->y, p2->y, EPS_C)) {
    i |= BRANCH_Y;
    nbits++;
  }
  if (distinguishable(p1->z, p2->z, EPS_C)) {
    i |= BRANCH_Z;
    nbits++;
  }

  /* Check that we're a patch on one surface */
  if (nbits != 2)
    return 0;
  if ((i & BRANCH_X) == 0 && (distinguishable(p1->x, b->x[0], EPS_C)) &&
      (distinguishable(p1->x, b->x[b->nx - 1], EPS_C))) {
    return 0;
  }
  if ((i & BRANCH_Y) == 0 && (distinguishable(p1->y, b->y[0], EPS_C)) &&
      (distinguishable(p1->y, b->y[b->ny - 1], EPS_C))) {
    return 0;
  }
  if ((i & BRANCH_Z) == 0 && (distinguishable(p1->z, b->z[0], EPS_C)) &&
      (distinguishable(p1->z, b->z[b->nz - 1], EPS_C))) {
    return 0;
  }

  /* Sanity checks for sizes */
  if ((i & BRANCH_X) != 0 &&
      (p1->x > p2->x || p1->x < b->x[0] || p2->x > b->x[b->nx - 1]))
    return 0;
  if ((i & BRANCH_Y) != 0 &&
      (p1->y > p2->y || p1->y < b->y[0] || p2->y > b->y[b->ny - 1]))
    return 0;
  if ((i & BRANCH_Z) != 0 &&
      (p1->z > p2->z || p1->z < b->z[0] || p2->z > b->z[b->nz - 1]))
    return 0;

  /* Check for sufficient spacing */
  if (i & BRANCH_X) {
    d = p2->x - p1->x;
    if (d > 0 && d < minspacing)
      return 0;
    for (j = 0; j < b->nx; j++) {
      d = fabs(b->x[j] - p1->x);
      if (d > 0 && d < minspacing)
        return 0;
      d = fabs(b->x[j] - p2->x);
      if (d > 0 && d < minspacing)
        return 0;
    }
  }
  if (i & BRANCH_Y) {
    d = p2->y - p1->y;
    if (d > 0 && d < minspacing)
      return 0;
    for (j = 0; j < b->ny; j++) {
      d = fabs(b->y[j] - p1->y);
      if (d > 0 && d < minspacing)
        return 0;
      d = fabs(b->y[j] - p2->y);
      if (d > 0 && d < minspacing)
        return 0;
    }
  }
  if (i & BRANCH_Z) {
    d = p2->z - p1->z;
    if (d > 0 && d < minspacing)
      return 0;
    for (j = 0; j < b->nz; j++) {
      d = fabs(b->z[j] - p1->z);
      if (d > 0 && d < minspacing)
        return 0;
      d = fabs(b->z[j] - p2->z);
      if (d > 0 && d < minspacing)
        return 0;
    }
  }

  return i;
}

/*************************************************************************
 cuboid_patch_to_bits:
    Convert a patch on a cuboid into a bit array representing membership.

 In: subd_box: a subdivided box upon which the patch is located
     v1: the lower-valued corner of the patch
     v2: the other corner
     bit_arr: a bit array to store the results.
 Out: returns 1 on failure, 0 on success.  The surface of the box is considered
      to be tiled with triangles in a particular order, and an array of bits is
      set to be 0 for each triangle that is not in the patch and 1 for each
      triangle that is.  (This is the internal format for regions.)
*************************************************************************/
int cuboid_patch_to_bits(struct subdivided_box *subd_box, struct vector3 *v1,
                         struct vector3 *v2, struct bit_array *bit_arr) {
  int dir_val;
  int patch_bitmask = check_patch(subd_box, v1, v2, GIGANTIC);
  if (!patch_bitmask)
    return 1;
  if ((patch_bitmask & BRANCH_X) == 0) {
    if (!distinguishable(subd_box->x[0], v1->x, EPS_C))
      dir_val = X_NEG;
    else
      dir_val = X_POS;
  } else if ((patch_bitmask & BRANCH_Y) == 0) {
    if (!distinguishable(subd_box->y[0], v1->y, EPS_C))
      dir_val = Y_NEG;
    else
      dir_val = Y_POS;
  } else {
    if (!distinguishable(subd_box->z[0], v1->z, EPS_C))
      dir_val = Z_NEG;
    else
      dir_val = Z_POS;
  }

  int a_lo, a_hi, b_lo, b_hi;
  int line, base;
  switch (dir_val) {
  case X_NEG:
    a_lo = bisect_near(subd_box->y, subd_box->ny, v1->y);
    a_hi = bisect_near(subd_box->y, subd_box->ny, v2->y);
    b_lo = bisect_near(subd_box->z, subd_box->nz, v1->z);
    b_hi = bisect_near(subd_box->z, subd_box->nz, v2->z);
    if (distinguishable(subd_box->y[a_lo], v1->y, EPS_C))
      return 1;
    if (distinguishable(subd_box->y[a_hi], v2->y, EPS_C))
      return 1;
    if (distinguishable(subd_box->z[b_lo], v1->z, EPS_C))
      return 1;
    if (distinguishable(subd_box->z[b_hi], v2->z, EPS_C))
      return 1;
    line = subd_box->ny - 1;
    base = 0;
    break;
  case X_POS:
    a_lo = bisect_near(subd_box->y, subd_box->ny, v1->y);
    a_hi = bisect_near(subd_box->y, subd_box->ny, v2->y);
    b_lo = bisect_near(subd_box->z, subd_box->nz, v1->z);
    b_hi = bisect_near(subd_box->z, subd_box->nz, v2->z);
    if (distinguishable(subd_box->y[a_lo], v1->y, EPS_C))
      return 1;
    if (distinguishable(subd_box->y[a_hi], v2->y, EPS_C))
      return 1;
    if (distinguishable(subd_box->z[b_lo], v1->z, EPS_C))
      return 1;
    if (distinguishable(subd_box->z[b_hi], v2->z, EPS_C))
      return 1;
    line = subd_box->ny - 1;
    base = (subd_box->ny - 1) * (subd_box->nz - 1);
    break;
  case Y_NEG:
    a_lo = bisect_near(subd_box->x, subd_box->nx, v1->x);
    a_hi = bisect_near(subd_box->x, subd_box->nx, v2->x);
    b_lo = bisect_near(subd_box->z, subd_box->nz, v1->z);
    b_hi = bisect_near(subd_box->z, subd_box->nz, v2->z);
    if (distinguishable(subd_box->x[a_lo], v1->x, EPS_C))
      return 1;
    if (distinguishable(subd_box->x[a_hi], v2->x, EPS_C))
      return 1;
    if (distinguishable(subd_box->z[b_lo], v1->z, EPS_C))
      return 1;
    if (distinguishable(subd_box->z[b_hi], v2->z, EPS_C))
      return 1;
    line = subd_box->nx - 1;
    base = 2 * (subd_box->ny - 1) * (subd_box->nz - 1);
    break;
  case Y_POS:
    a_lo = bisect_near(subd_box->x, subd_box->nx, v1->x);
    a_hi = bisect_near(subd_box->x, subd_box->nx, v2->x);
    b_lo = bisect_near(subd_box->z, subd_box->nz, v1->z);
    b_hi = bisect_near(subd_box->z, subd_box->nz, v2->z);
    if (distinguishable(subd_box->x[a_lo], v1->x, EPS_C))
      return 1;
    if (distinguishable(subd_box->x[a_hi], v2->x, EPS_C))
      return 1;
    if (distinguishable(subd_box->z[b_lo], v1->z, EPS_C))
      return 1;
    if (distinguishable(subd_box->z[b_hi], v2->z, EPS_C))
      return 1;
    line = subd_box->nx - 1;
    base = 2 * (subd_box->ny - 1) * (subd_box->nz - 1) +
           (subd_box->nx - 1) * (subd_box->nz - 1);
    break;
  case Z_NEG:
    a_lo = bisect_near(subd_box->x, subd_box->nx, v1->x);
    a_hi = bisect_near(subd_box->x, subd_box->nx, v2->x);
    b_lo = bisect_near(subd_box->y, subd_box->ny, v1->y);
    b_hi = bisect_near(subd_box->y, subd_box->ny, v2->y);
    if (distinguishable(subd_box->x[a_lo], v1->x, EPS_C))
      return 1;
    if (distinguishable(subd_box->x[a_hi], v2->x, EPS_C))
      return 1;
    if (distinguishable(subd_box->y[b_lo], v1->y, EPS_C))
      return 1;
    if (distinguishable(subd_box->y[b_hi], v2->y, EPS_C))
      return 1;
    line = subd_box->nx - 1;
    base = 2 * (subd_box->ny - 1) * (subd_box->nz - 1) +
           2 * (subd_box->nx - 1) * (subd_box->nz - 1);
    break;
  case Z_POS:
    a_lo = bisect_near(subd_box->x, subd_box->nx, v1->x);
    a_hi = bisect_near(subd_box->x, subd_box->nx, v2->x);
    b_lo = bisect_near(subd_box->y, subd_box->ny, v1->y);
    b_hi = bisect_near(subd_box->y, subd_box->ny, v2->y);
    if (distinguishable(subd_box->x[a_lo], v1->x, EPS_C))
      return 1;
    if (distinguishable(subd_box->x[a_hi], v2->x, EPS_C))
      return 1;
    if (distinguishable(subd_box->y[b_lo], v1->y, EPS_C))
      return 1;
    if (distinguishable(subd_box->y[b_hi], v2->y, EPS_C))
      return 1;
    line = subd_box->nx - 1;
    base = 2 * (subd_box->ny - 1) * (subd_box->nz - 1) +
           2 * (subd_box->nx - 1) * (subd_box->nz - 1) +
           (subd_box->nx - 1) * (subd_box->ny - 1);
    break;
  default:
    UNHANDLED_CASE(dir_val);
    /*return 1;*/
  }

  set_all_bits(bit_arr, 0);

  if (a_lo == 0 && a_hi == line) {
    set_bit_range(bit_arr, 2 * (base + line * b_lo + a_lo),
                  2 * (base + line * (b_hi - 1) + (a_hi - 1)) + 1, 1);
  } else {
    for (int i = b_lo; i < b_hi; i++) {
      set_bit_range(bit_arr, 2 * (base + line * i + a_lo),
                    2 * (base + line * i + (a_hi - 1)) + 1, 1);
    }
  }

  return 0;
}

int mcell_check_for_region(char *region_name, struct object *obj_ptr) {
  struct region_list *reg_list;
  for (reg_list = obj_ptr->regions; reg_list != NULL;
       reg_list = reg_list->next) {
    if (strcmp(region_name, reg_list->reg->sym->name) == 0) {
      return 1;
    }
  }
  return 0;
}

/**************************************************************************
 mcell_create_region:
    Create a named region on an object.

 In:  state: the simulation state
      obj_ptr: object upon which to create a region
      name: region name to create
 Out: region object, or NULL if there was an error (region already exists or
      allocation failed)
 NOTE: This is similar to mdl_create_region
**************************************************************************/
struct region *mcell_create_region(MCELL_STATE *state, struct object *obj_ptr,
                                   char *name) {
  struct region *reg_ptr;
  struct region_list *reg_list_ptr;
  no_printf("Creating new region: %s\n", name);
  if ((reg_ptr = make_new_region(state->dg_parse, state, obj_ptr->sym->name, name)) == NULL) {
    return NULL;
  }
  if ((reg_list_ptr =
           CHECKED_MALLOC_STRUCT(struct region_list, "region list")) == NULL) {
    // mdlerror_fmt(parse_state,
    //             "Out of memory while creating object region '%s'",
    //             rp->sym->name);
    return NULL;
  }
  reg_ptr->region_last_name = name;
  reg_ptr->parent = obj_ptr;
  char *region_name = CHECKED_SPRINTF("%s,%s", obj_ptr->sym->name, name);
  if (!mcell_check_for_region(region_name, obj_ptr)) {
    reg_list_ptr->reg = reg_ptr;
    reg_list_ptr->next = obj_ptr->regions;
    obj_ptr->regions = reg_list_ptr;
    obj_ptr->num_regions++;
  }
  else {
    free(reg_list_ptr);
  }
  free(region_name);
  return reg_ptr;
}

/*************************************************************************
 make_new_region:
    Create a new region, adding it to the global symbol table.  The region must
    not be defined yet.  The region is not added to the object's list of
    regions.

    full region names of REG type symbols stored in main symbol table have the
    form:
         metaobj.metaobj.poly,region_last_name

 In:  state: the simulation state
      obj_name: fully qualified object name
      region_last_name: name of the region to define
 Out: The newly created region
 NOTE: This is similar to mdl_make_new_region
*************************************************************************/
struct region *make_new_region(
    struct dyngeom_parse_vars *dg_parse,
    MCELL_STATE *state,
    char *obj_name,
    char *region_last_name) {

  char *region_name;
  region_name = CHECKED_SPRINTF("%s,%s", obj_name, region_last_name);
  if (region_name == NULL) {
    return NULL;
  }

  struct sym_entry *sym_ptr;
  if (((sym_ptr = retrieve_sym(region_name, state->reg_sym_table)) != NULL) && (sym_ptr->count == 0)) {
    free(region_name);
    if (sym_ptr->count == 0) {
      sym_ptr->count = 1;
      return (struct region *)sym_ptr->value;
    }
    else {
      return NULL; 
    }
  }

  if ((sym_ptr = store_sym(region_name, REG, state->reg_sym_table, NULL)) ==
      NULL) {
    free(region_name);
    return NULL;
  }

  free(region_name);
  return (struct region *)sym_ptr->value;
}

/*****************************************************************************
 *
 * mcell_add_to_vertex_list creates a linked list of mesh vertices belonging
 * to a polygon object
 *
 * During the first invocation of this function, NULL should be provided for
 * vertices to initialize a new vertex list. On subsecquent invocations the
 * current vertex_list should be provided as parameter vertices to which the
 * new vertex will be appended.
 *
 *****************************************************************************/
struct vertex_list *mcell_add_to_vertex_list(double x, double y, double z,
                                             struct vertex_list *vertices) {
  struct vertex_list *verts = (struct vertex_list *)CHECKED_MALLOC_STRUCT(
      struct vertex_list, "vertex list");
  if (verts == NULL) {
    return NULL;
  }

  struct vector3 *v =
      (struct vector3 *)CHECKED_MALLOC_STRUCT(struct vector3, "vector");
  if (v == NULL) {
    free(verts);
    return NULL;
  }
  v->x = x;
  v->y = y;
  v->z = z;

  verts->vertex = v;
  verts->next = vertices;

  return verts;
}

/*****************************************************************************
 *
 * mcell_add_to_connection_list creates a linked list of element connections
 * describing a polygon object.
 *
 * During the first invocation of this function, NULL should be provided for
 * elements to initialize a new element connection list. On subsecquent
 * invocations the current element_connection_list should be provided as
 * parameter elements to which the new element connection will be appended.
 *
 *****************************************************************************/
struct element_connection_list *
mcell_add_to_connection_list(int v1, int v2, int v3,
                             struct element_connection_list *elements) {
  struct element_connection_list *elems =
      (struct element_connection_list *)CHECKED_MALLOC_STRUCT(
          struct element_connection_list, "element connection list");
  if (elems == NULL) {
    return NULL;
  }

  int *e = (int *)CHECKED_MALLOC_ARRAY(int, 3, "element connections");
  if (e == NULL) {
    free(elems);
    return NULL;
  }
  e[0] = v1;
  e[1] = v2;
  e[2] = v3;

  elems->n_verts = 3;
  elems->indices = e;
  elems->next = elements;

  return elems;
}

/**************************************************************************
 mcell_set_region_elements:
    Set the elements for a region, normalizing the region if it's on a polygon
    list object.

 In: rgn:  region to receive elements
     elements: elements comprising region
     normalize_now: flag indicating whether to normalize right now
 Out: returns 1 on failure, 0 on success.
 NOTE: Almost the same as mdl_set_region_elements
**************************************************************************/
int mcell_set_region_elements(struct region *rgn, struct element_list *elements,
                              int normalize_now) {
  rgn->element_list_head = elements;
  if (normalize_now)
    return normalize_elements(rgn, 0);
  else
    return 0;
}

/**************************************************************************
 mcell_add_to_region_list:

 In: elements: the list of elements for a region
     region_idx: the index of the region
 Out: the updated list of elements for a region
**************************************************************************/
struct element_list *mcell_add_to_region_list(struct element_list *elements,
                                              u_int region_idx) {

  struct element_list *elem = (struct element_list *)CHECKED_MALLOC_STRUCT(
      struct element_list, "element list");
  if (elem == NULL) {
    return NULL;
  }

  elem->next = elements;
  elem->begin = region_idx;
  elem->end = region_idx;
  elem->special = NULL;

  return elem;
}

/****************************************************************************
 *
 * static helper functions
 *
 ****************************************************************************/

/**************************************************************************
 is_region_degenerate:
    Check a region for degeneracy.

 In: reg_ptr: region to check
 Out: 1 if degenerate, 0 if not
**************************************************************************/
int is_region_degenerate(struct region *reg_ptr) {
  for (int i = 0; i < reg_ptr->membership->nbits; i++) {
    if (get_bit(reg_ptr->membership, i)) {
      return 0;
    }
  }
  return 1;
}

struct sym_entry *
mcell_get_obj_sym(struct object *obj) {
  return obj->sym;
}
struct sym_entry *
mcell_get_reg_sym(struct region *reg) {
  return reg->sym;
}

struct poly_object_list* mcell_add_to_poly_obj_list(
  struct poly_object_list* poly_obj_list,
  char *obj_name,
  struct vertex_list *vertices,
  int num_vert,
  struct element_connection_list *connections,
  int num_conn,
  struct element_list *surf_reg_faces,
  char *reg_name) {

  struct poly_object_list *pobj_list = (struct poly_object_list *)CHECKED_MALLOC_STRUCT(
      struct poly_object_list, "poly obj list");
  if (pobj_list == NULL) {
    return NULL;
  }
  char *object_name = CHECKED_STRDUP(obj_name, "object name");
  char *region_name = CHECKED_STRDUP(reg_name, "object name");

  pobj_list->vertices = vertices;
  pobj_list->num_vert = num_vert;
  pobj_list->connections = connections;
  pobj_list->num_conn = num_conn;
  pobj_list->surf_reg_faces = surf_reg_faces;
  pobj_list->reg_name = region_name;
  pobj_list->obj_name = object_name;
  pobj_list->next = poly_obj_list;

  return pobj_list;
}
