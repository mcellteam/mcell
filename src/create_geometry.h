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

#ifndef CREATE_GEOMETRY_H
#define CREATE_GEOMETRY_H
#include "libmcell.h"

struct element_connection_list
{
  struct element_connection_list *next;
  int n_verts;
  int *indices;
};

// Finalize the polygon list, cleaning up any state updates that were made when
// we started creating the polygon.
int finish_polygon_list(struct sym_table *sym_ptr);
// Clean up the regions on an object, eliminating any removed walls.
void remove_gaps_from_regions(struct object *obj_ptr);
// Check a box or polygon list object for degeneracy.
int check_degenerate_polygon_list(struct object *obj_ptr);
int is_region_degenerate(struct region *reg_ptr);
struct polygon_object *allocate_polygon_object(char const *desc);
struct element_list *new_element_list(unsigned int begin, unsigned int end);
void free_vertex_list(struct vertex_list *vert_list);
void free_connection_list(struct element_connection_list *elem_conn_list);
struct polygon_object * new_polygon_list(
  MCELL_STATE* state,
  struct sym_table *sym,
  int n_vertices,
  struct vertex_list *vertices,
  int n_connections,
  struct element_connection_list *connections);
int normalize_elements(struct region *reg,
                       int existing);
struct region *create_region(MCELL_STATE* state, struct object *objp, char *name);
struct region *make_new_region(MCELL_STATE* state,
                               char *obj_name,
                               char *region_last_name);
int count_cuboid_elements(struct subdivided_box *sb);
int cuboid_patch_to_bits(struct subdivided_box *subd_box,
                         struct vector3 *v1,
                         struct vector3 *v2,
                         struct bit_array *bit_arr);
int check_patch(struct subdivided_box *b,
                struct vector3 *p1,
                struct vector3 *p2,
                double egd);
#endif
