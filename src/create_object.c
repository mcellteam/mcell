/***********************************************************************************
 *                                                                                 *
 * Copyright (C) 2006-2014 by *
 * The Salk Institute for Biological Studies and *
 * Pittsburgh Supercomputing Center, Carnegie Mellon University *
 *                                                                                 *
 * This program is free software; you can redistribute it and/or *
 * modify it under the terms of the GNU General Public License *
 * as published by the Free Software Foundation; either version 2 *
 * of the License, or (at your option) any later version. *
 *                                                                                 *
 * This program is distributed in the hope that it will be useful, *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the *
 * GNU General Public License for more details. *
 *                                                                                 *
 * You should have received a copy of the GNU General Public License *
 * along with this program; if not, write to the Free Software *
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301,
 *USA. *
 *                                                                                 *
 ***********************************************************************************/

#include "create_object.h"
#include "logging.h"
#include "sym_table.h"
#include "mem_util.h"
#include "util.h"
#include "vector.h"

#include <stdlib.h>
#include <assert.h>

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
 transform_translate:
    Apply a translation to the given transformation matrix.

 In:  state: system state
      mat: transformation matrix
      xlat: translation vector
 Out: translation is right-multiplied into the transformation matrix
*************************************************************************/
void transform_translate(MCELL_STATE *state, double (*mat)[4],
                         struct vector3 *xlat) {
  double tm[4][4];
  struct vector3 scaled_xlat = *xlat;
  scaled_xlat.x *= state->r_length_unit;
  scaled_xlat.y *= state->r_length_unit;
  scaled_xlat.z *= state->r_length_unit;
  init_matrix(tm);
  translate_matrix(tm, tm, &scaled_xlat);
  mult_matrix(mat, tm, mat, 4, 4, 4);
}

/*************************************************************************
 transform_scale:
    Apply a scale to the given transformation matrix.

 In:  mat: transformation matrix
      scale: scale vector
 Out: scale is right-multiplied into the transformation matrix
*************************************************************************/
void transform_scale(double (*mat)[4], struct vector3 *scale) {
  double tm[4][4];
  init_matrix(tm);
  scale_matrix(tm, tm, scale);
  mult_matrix(mat, tm, mat, 4, 4, 4);
}

/*************************************************************************
 transform_rotate:
    Apply a rotation to the given transformation matrix.

 In:  mat: transformation matrix
      axis: axis of rotation
      angle: angle of rotation (degrees!)
 Out: 0 on success, 1 on failure; rotation is right-multiplied into the
      transformation matrix
*************************************************************************/
int transform_rotate(double (*mat)[4], struct vector3 *axis, double angle) {
  double tm[4][4];
  if (!distinguishable(vect_length(axis), 0.0, EPS_C)) {
    return 1;
  }
  init_matrix(tm);
  rotate_matrix(tm, tm, axis, angle);
  mult_matrix(mat, tm, mat, 4, 4, 4);
  return 0;
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
