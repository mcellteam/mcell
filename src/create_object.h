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

#ifndef CREATE_OBJECT_H
#define CREATE_OBJECT_H
#include "libmcell.h"

struct object_creation
{
  struct name_list *object_name_list;
  struct name_list *object_name_list_end;
};

char *push_object_name(struct object_creation *obj_creation, char *name);
struct object *make_new_object(MCELL_STATE *state, char *obj_name);
void pop_object_name(struct object_creation *obj_creation);
/* Adds children to a meta-object, aggregating counts of walls and vertices
 * from the children into the specified parent. The children should already
 * have their parent pointers set. */
void add_child_objects(struct object *parent,
                       struct object *child_head,
                       struct object *child_tail);

#endif
