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

#ifndef CREATE_OBJECT_H
#define CREATE_OBJECT_H
#include "libmcell.h"

char *push_object_name(struct object_creation *obj_creation, char *name);

void pop_object_name(struct object_creation *obj_creation);

/* Adds children to a meta-object, aggregating counts of walls and vertices
 * from the children into the specified parent. The children should already
 * have their parent pointers set. */
void add_child_objects(struct object *parent, struct object *child_head,
  struct object *child_tail);

void check_regions(struct object *rootInstance, struct object *child_head);

// Apply a translation to the given transformation matrix.
void transform_translate(MCELL_STATE *state, double (*mat)[4],
                         struct vector3 *xlat);

// Apply a scale to the given transformation matrix.
void transform_scale(double (*mat)[4], struct vector3 *scale);

// Apply a rotation to the given transformation matrix.
int transform_rotate(double (*mat)[4], struct vector3 *axis, double angle);

struct object *common_ancestor(struct object *a, struct object *b);

#endif
