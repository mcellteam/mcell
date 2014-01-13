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

#include "mcell_structs.h"
#include "util.h"
#include "logging.h"
#include "sym_table.h"
#include "mem_util.h"
#include "libmcell.h"
#include "create_geometry.h"

#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <assert.h>



/**************************************************************************
 mdl_finish_polygon_list:
    Finalize the polygon list, cleaning up any state updates that were made
    when we started creating the polygon.

 In: sym_ptr: symbol for the completed polygon
 Out: 1 on failure, 0 on success
**************************************************************************/
int
finish_polygon_list(struct sym_table *sym_ptr)
{
  struct object *obj_ptr = (struct object *) sym_ptr->value;
  remove_gaps_from_regions(obj_ptr);
  no_printf("Polygon list %s defined:\n", sym_ptr->name);
  //no_printf(" n_verts = %d\n", mpvp->current_polygon->n_verts);
  //no_printf(" n_walls = %d\n", mpvp->current_polygon->n_walls);
  if (check_degenerate_polygon_list(obj_ptr)) {
    return 1;
  }
  return 0;
}



/*************************************************************************
 remove_gaps_from_regions:
    Clean up the regions on an object, eliminating any removed walls.

 In: obj_ptr: an object with regions
 Out: Any walls that have been removed from the object are removed from every
      region on that object.
*************************************************************************/
void
remove_gaps_from_regions(struct object *obj_ptr)
{
  if (obj_ptr->object_type != BOX_OBJ && obj_ptr->object_type != POLY_OBJ) {
    return;
  }

  struct polygon_object *poly_obj_ptr= (struct polygon_object*) obj_ptr->contents;
  struct region_list *reg_list;
  for (reg_list = obj_ptr->regions; reg_list != NULL; reg_list = reg_list->next)
  {
    no_printf("Checking region %s\n", reg_list->reg->sym->name);
    if (strcmp(reg_list->reg->region_last_name, "REMOVED") == 0)
    {
      no_printf("Found a REMOVED region\n");
      reg_list->reg->surf_class = NULL;
      bit_operation(poly_obj_ptr->side_removed, reg_list->reg->membership, '+');
      set_all_bits(reg_list->reg->membership, 0);
    }
  }

  int missing=0;
  for (int n_side=0; n_side < poly_obj_ptr->side_removed->nbits; ++n_side)
  {
    if (get_bit(poly_obj_ptr->side_removed, n_side)) {
      missing++;
    }
  }

  obj_ptr->n_walls_actual = poly_obj_ptr->n_walls - missing;

  for (reg_list=obj_ptr->regions; reg_list != NULL; reg_list=reg_list->next) {
    bit_operation(reg_list->reg->membership, poly_obj_ptr->side_removed, '-');
  }

#ifdef DEBUG
  printf("Sides for %s: ", obj_ptr->sym->name);
  for (unsigned int n_side=0; n_side<poly_obj_ptr->side_removed->nbits; ++n_side)
  {
    if (get_bit(poly_obj_ptr->side_removed, n_side)) {
      printf("-");
    }
    else {
      printf("#");
    }
  }
  printf("\n");
  for (reg_list=obj_ptr->regions; reg_list!=NULL; reg_list=reg_list->next)
  {
    printf("Sides for %s: ", reg_list->reg->sym->name);
    for (unsigned int n_side=0; n_side<reg_list->reg->membership->nbits; ++ n_side)
    {
      if (get_bit(reg_list->reg->membership, n_side)) {
        printf("+");
      }
      else {
        printf(".");
      }
    }
    printf("\n");
  }
#endif
}



/**************************************************************************
 check_degenerate_polygon_list:
    Check a box or polygon list object for degeneracy.

 In: obj_ptr: the object to validate
 Out: 0 if valid, 1 if invalid
**************************************************************************/
int
check_degenerate_polygon_list(struct object *obj_ptr)
{
  // Check for a degenerate (empty) object and regions.
  struct region_list *reg_list;
  for (reg_list = obj_ptr->regions; reg_list != NULL; reg_list = reg_list->next)
  {
    if (!is_region_degenerate(reg_list->reg)) {
      continue;
    }
    if (strcmp(reg_list->reg->region_last_name, "ALL") == 0) {
      return 1;
    }
    else if (strcmp(reg_list->reg->region_last_name, "REMOVED") != 0) {
      return 1;
    }
  }
  return 0;
}



/**************************************************************************
 is_region_degenerate:
    Check a region for degeneracy.

 In: reg_ptr: region to check
 Out: 1 if degenerate, 0 if not
**************************************************************************/
int
is_region_degenerate(struct region *reg_ptr)
{ 
  for (int i = 0; i < reg_ptr->membership->nbits; i++)
  {
    if (get_bit(reg_ptr->membership, i)) {
      return 0;
    }
  }
  return 1;
}



/**************************************************************************
 allocate_polygon_object:
    Allocate a polygon object.

 In: desc: object type description
 Out: polygon object on success, NULL on failure
**************************************************************************/
struct polygon_object *
allocate_polygon_object(char const *desc)
{
  struct polygon_object *poly_obj_ptr;
  if ((poly_obj_ptr = CHECKED_MALLOC_STRUCT(struct polygon_object, desc)) == NULL) {
    return NULL;
  }
  poly_obj_ptr->n_verts=0;
  poly_obj_ptr->parsed_vertices=NULL;
  poly_obj_ptr->n_walls=0;
  poly_obj_ptr->element=NULL;
  poly_obj_ptr->sb = NULL;
  poly_obj_ptr->side_removed = NULL;
  return poly_obj_ptr;
}



/**************************************************************************
 new_element_list:
    Create a new element list for a region description.

 In: begin: starting side number for this element
     end: ending side number for this element
 Out: element list, or NULL if allocation fails
**************************************************************************/
struct element_list *
new_element_list(unsigned int begin,
                 unsigned int end)
{

  struct element_list *elem_list = CHECKED_MALLOC_STRUCT(
    struct element_list, "region element");
  if (elem_list == NULL) {
    return NULL;
  }
  elem_list->special=NULL;
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
void
free_vertex_list(struct vertex_list *vert_list)
{
  while (vert_list)
  {
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
void
free_connection_list(struct element_connection_list *elem_conn_list)
{
  while (elem_conn_list)
  {
    struct element_connection_list *next = elem_conn_list->next;
    free(elem_conn_list->indices);
    free(elem_conn_list);
    elem_conn_list = next;
  }
}
