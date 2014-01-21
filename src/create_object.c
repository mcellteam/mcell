/***********************************************************************************
 *                                                                                 *
 * Copyright (C) 2006-2013 by                                                      *
 * The Salk Institute for Biological Studies and                                   *
 * Pittsburgh Supercomputing Center, Carnegie Mellon University                    *
 *                                                                                 *
 * This program is free software; you can redistribute it and/or                   *
 * modify it under the terms of the GNU General Public License                     *
 * as published by the Free Software Foundation; either version 2                  *
 * of the License, or (at your option) any later version.                          *
 *                                                                                 *
 * This program is distributed in the hope that it will be useful,                 *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of                  *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                   *
 * GNU General Public License for more details.                                    *
 *                                                                                 *
 * You should have received a copy of the GNU General Public License               *
 * along with this program; if not, write to the Free Software                     *
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA. *
 *                                                                                 *
 ***********************************************************************************/

#include "create_object.h"
#include "logging.h"
#include "sym_table.h"
#include "mem_util.h"

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
char *
push_object_name(struct object_creation *obj_creation, char *name)
{

  // Initialize object name list 
  if (obj_creation->object_name_list == NULL)
  {
    obj_creation->object_name_list = CHECKED_MALLOC_STRUCT(
      struct name_list, "object name stack");
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
  if (obj_creation->object_name_list_end->name == NULL)
  {
    obj_creation->object_name_list_end->name = name;
    return obj_creation->object_name_list_end->name;
  }

  // If we've run out of name list components, create a new one
  struct name_list *name_list_ptr;
  if (obj_creation->object_name_list_end->next == NULL)
  {
    name_list_ptr = CHECKED_MALLOC_STRUCT(
      struct name_list, "object name stack");
    if (name_list_ptr == NULL) {
      return NULL;
    }

    name_list_ptr->next = NULL;
    name_list_ptr->prev = obj_creation->object_name_list_end;
    obj_creation->object_name_list_end->next = name_list_ptr;
  }
  else {
    name_list_ptr = obj_creation->object_name_list_end->next;
  }

  // Create new name 
  name_list_ptr->name = CHECKED_SPRINTF(
    "%s.%s",
    obj_creation->object_name_list_end->name,
    name);
  if (name_list_ptr->name == NULL) {
    return NULL;
  }

  obj_creation->object_name_list_end = name_list_ptr;

  return name_list_ptr->name;
}



/*************************************************************************
 make_new_object:
    Create a new object, adding it to the global symbol table.

 In:  state: system state
      obj_name: fully qualified object name
 Out: the newly created object
*************************************************************************/
struct object *
make_new_object(MCELL_STATE* state, char *obj_name)
{

  struct sym_table *symbol;
  if ((symbol = store_sym(obj_name, OBJ, state->obj_sym_table, NULL)) == NULL) {
    return NULL;
  }

  return (struct object *) symbol->value;
}



/*************************************************************************
 mdl_pop_object_name:
    Remove the trailing name component fromt the name list.  It is expected
    that ownership of the name pointer has passed to someone else, or that the
    pointer has been freed already.

 In:  obj_creation: information about object being created
 Out: object name stack is updated
*************************************************************************/
void
pop_object_name(struct object_creation *obj_creation)
{
  if (obj_creation->object_name_list_end->name != NULL) {
    free(obj_creation->object_name_list_end->name);
  }
  if (obj_creation->object_name_list_end->prev != NULL) {
    obj_creation->object_name_list_end = obj_creation->object_name_list_end->prev;
  }
  else {
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
void
add_child_objects(struct object *parent,
                  struct object *child_head,
                  struct object *child_tail)
{
  if (parent->first_child == NULL)
    parent->first_child = child_head;
  if (parent->last_child != NULL)
    parent->last_child->next = child_head;
  parent->last_child = child_tail;
  child_tail->next = NULL;
  while (child_head != NULL)
  {
    assert(child_head->parent == parent);
    parent->n_walls        += child_head->n_walls;
    parent->n_walls_actual += child_head->n_walls_actual;
    parent->n_verts        += child_head->n_verts;
    child_head = child_head->next;
  }
}



/*************************************************************************
 start_object:
    Create a new object, adding it to the global symbol table.  the object must
    not be defined yet.  The qualified name of the object will be built by
    adding to the object_name_list, and the object is made the "current_object"
    in the mdl parser state.  Because of these side effects, it is vital to
    call mdl_finish_object at the end of the scope of the object created here.

 In:  obj_creation: information about object being created
      obj_name: unqualified object name
 Out: the newly created object
 NOTE: This is very similar to mdl_start_object, but there is no parse state
*************************************************************************/
struct sym_table *
start_object(MCELL_STATE* state,
             struct object_creation *obj_creation,
             char *name)
{
  // Create new fully qualified name.
  char *new_name;
  if ((new_name = push_object_name(obj_creation, name)) == NULL)
  {
    free(name);
    return NULL;
  }

  // Create the symbol, if it doesn't exist yet.
  struct object *obj_ptr = make_new_object(state, new_name);
  if (obj_ptr == NULL)
  {
    free(name);
    free(new_name);
    return NULL;
  }

  struct sym_table *sym_ptr = obj_ptr->sym;
  obj_ptr->last_name = name;
  no_printf("Creating new object: %s\n", new_name);

  // Set parent object, make this object "current". 
  obj_ptr->parent = obj_creation->current_object;
  obj_creation->current_object = obj_ptr;

  return sym_ptr;
}
