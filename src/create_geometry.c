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

#ifdef DEBUG
  printf("Sides for %s: ", obj_ptr->sym->name);
  for (unsigned int n_side = 0; n_side < poly_obj_ptr->side_removed->nbits;
       ++n_side) {
    if (get_bit(poly_obj_ptr->side_removed, n_side)) {
      printf("-");
    } else {
      printf("#");
    }
  }
  printf("\n");
  for (reg_list = obj_ptr->regions; reg_list != NULL;
       reg_list = reg_list->next) {
    printf("Sides for %s: ", reg_list->reg->sym->name);
    for (unsigned int n_side = 0; n_side < reg_list->reg->membership->nbits;
         ++n_side) {
      if (get_bit(reg_list->reg->membership, n_side)) {
        printf("+");
      } else {
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

  if (p1->x != p2->x) {
    i |= BRANCH_X;
    nbits++;
  }
  if (p1->y != p2->y) {
    i |= BRANCH_Y;
    nbits++;
  }
  if (p1->z != p2->z) {
    i |= BRANCH_Z;
    nbits++;
  }

  /* Check that we're a patch on one surface */
  if (nbits != 2)
    return 0;
  if ((i & BRANCH_X) == 0 && p1->x != b->x[0] && p1->x != b->x[b->nx - 1])
    return 0;
  if ((i & BRANCH_Y) == 0 && p1->y != b->y[0] && p1->y != b->y[b->ny - 1])
    return 0;
  if ((i & BRANCH_Z) == 0 && p1->z != b->z[0] && p1->z != b->z[b->nz - 1])
    return 0;

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
    if (subd_box->x[0] == v1->x)
      dir_val = X_NEG;
    else
      dir_val = X_POS;
  } else if ((patch_bitmask & BRANCH_Y) == 0) {
    if (subd_box->y[0] == v1->y)
      dir_val = Y_NEG;
    else
      dir_val = Y_POS;
  } else {
    if (subd_box->z[0] == v1->z)
      dir_val = Z_NEG;
    else
      dir_val = Z_POS;
  }

  int a_lo, a_hi, b_lo, b_hi;
  int line, base;
  switch (dir_val) {
  case NODIR:
    return 1;
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
    return 1;
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
        return 1;
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
          return 1;
        }
        if (existing) {
          // mcell_internal_error("Attempt to create a PATCH on an already
          // triangulated BOX.");
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

#ifdef DEBUG
  printf("Normalized membership of %s: ", reg->sym->name);
  for (i = 0; i < reg->membership->nbits; i++) {
    if (get_bit(reg->membership, i)) {
      printf("X");
    } else {
      printf("_");
    }
  }
  printf("\n");
#endif

  return 0;
}

/**************************************************************************
 create_region:
    Create a named region on an object.

 In:  state: the simulation state
      obj_ptr: object upon which to create a region
      name: region name to create
 Out: region object, or NULL if there was an error (region already exists or
      allocation failed)
 NOTE: This is similar to mdl_create_region
**************************************************************************/
struct region *create_region(MCELL_STATE *state, struct object *obj_ptr,
                             char *name) {
  struct region *reg_ptr;
  struct region_list *reg_list_ptr;
  no_printf("Creating new region: %s\n", name);
  if ((reg_ptr = make_new_region(state, obj_ptr->sym->name, name)) == NULL) {
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
  reg_list_ptr->reg = reg_ptr;
  reg_list_ptr->next = obj_ptr->regions;
  obj_ptr->regions = reg_list_ptr;
  obj_ptr->num_regions++;
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
struct region *make_new_region(MCELL_STATE *state, char *obj_name,
                               char *region_last_name) {
  char *region_name;
  region_name = CHECKED_SPRINTF("%s,%s", obj_name, region_last_name);
  if (region_name == NULL) {
    return NULL;
  }

  struct sym_table *sym_ptr;
  if ((sym_ptr = store_sym(region_name, REG, state->reg_sym_table, NULL)) ==
      NULL) {
    free(region_name);
    return NULL;
  }

  free(region_name);
  return (struct region *)sym_ptr->value;
}
