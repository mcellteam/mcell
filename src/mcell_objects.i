/******************************************************************************
 *
 * Copyright (C) 2006-2015 by
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
struct object_creation {
  struct name_list *object_name_list;
  struct name_list *object_name_list_end;
  struct object *current_object;
};

struct poly_object {
  char *obj_name;
  struct vertex_list *vertices;
  int num_vert;
  struct element_connection_list *connections;
  int num_conn;
};

struct poly_object_list {
  char *obj_name;
  struct vertex_list *vertices;
  int num_vert;
  struct element_connection_list *connections;
  int num_conn;
  struct element_list *surf_reg_faces;
  char *reg_name;
  struct poly_object_list *next;
};

%typemap(in) struct object **new_object (struct object *temp) {
  $1 = &temp;
}

%typemap(argout) struct object **new_object {
  %set_output(SWIG_NewPointerObj(SWIG_as_voidptr(*$1), $*1_descriptor, SWIG_POINTER_OWN));
}

/* object creation */
MCELL_STATUS mcell_create_instance_object(MCELL_STATE *state, char *name,
                                          struct object **new_object);

MCELL_STATUS mcell_create_poly_object(MCELL_STATE *state, struct object *parent,
                                      struct poly_object *poly_obj,
                                      struct object **new_object);

struct polygon_object *
new_polygon_list(MCELL_STATE *state, struct object *obj_ptr, int n_vertices,
                 struct vertex_list *vertices, int n_connections,
                 struct element_connection_list *connections);

struct object *make_new_object(
    struct dyngeom_parse_vars *dg_parse,
    struct sym_table_head *obj_sym_table,
    char *obj_name,
    int *error_code);

char *push_object_name(struct object_creation *obj_creation, char *name);

void pop_object_name(struct object_creation *obj_creation);

/* helper functions for creating and deleting vertex and
 * element_connection lists */
struct vertex_list *mcell_add_to_vertex_list(double x, double y, double z,
                                             struct vertex_list *vertices);

void free_vertex_list(struct vertex_list *vert_list);

struct element_connection_list *
mcell_add_to_connection_list(int v1, int v2, int v3,
                             struct element_connection_list *elements);

void free_connection_list(struct element_connection_list *elem_conn_list);

int mcell_set_region_elements(struct region *rgn, struct element_list *elements,
                              int normalize_now);

struct element_list *mcell_add_to_region_list(struct element_list *elements,
                                              unsigned int region_idx);

/* Adds children to a meta-object, aggregating counts of walls and vertices
 * from the children into the specified parent. The children should already
 * have their parent pointers set. */
void add_child_objects(struct object *parent, struct object *child_head,
                       struct object *child_tail);

/* create regions */
struct region *mcell_create_region(MCELL_STATE *state, struct object *objp,
                                   char *name);

struct region *make_new_region(
    struct dyngeom_parse_vars *dg_parse,
    MCELL_STATE *state,
    char *obj_name,
    char *region_last_name);

/* Clean up the regions on an object, eliminating any removed walls. */
void remove_gaps_from_regions(struct object *obj_ptr);

/* lower level helper functions */
int check_degenerate_polygon_list(struct object *obj_ptr);

struct object *common_ancestor(struct object *a, struct object *b);

struct polygon_object *allocate_polygon_object(char const *desc);

struct element_list *new_element_list(unsigned int begin, unsigned int end);

int normalize_elements(struct region *reg, int existing);

int count_cuboid_elements(struct subdivided_box *sb);

int cuboid_patch_to_bits(struct subdivided_box *subd_box, struct vector3 *v1,
                         struct vector3 *v2, struct bit_array *bit_arr);

int check_patch(struct subdivided_box *b, struct vector3 *p1,
                struct vector3 *p2, double egd);

struct sym_entry *mcell_get_obj_sym(struct object *obj);
struct sym_entry *mcell_get_reg_sym(struct region *reg);

struct poly_object_list* mcell_add_to_poly_obj_list(
  struct poly_object_list* poly_obj_list,
  char *obj_name,
  struct vertex_list *vertices,
  int num_vert,
  struct element_connection_list *connections,
  int num_conn,
  struct element_list *surf_reg_faces,
  char *reg_name);
