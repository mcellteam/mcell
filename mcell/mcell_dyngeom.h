/******************************************************************************
 *
 * Copyright (C) 2006-2014 by
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
 *****************************************************************************/

#ifndef MCELL_DYNGEOM_H
#define MCELL_DYNGEOM_H

#include "mcell_objects.h"

struct mesh_region_string_buffs {
  struct string_buffer *old_inst_mesh_names;
  struct string_buffer *old_region_names;
};

int mcell_add_dynamic_geometry_file(char *dynamic_geometry_filepath,
                                    struct mdlparse_vars *parse_state);

int mcell_change_geometry(struct volume *state, struct poly_object_list *pobj_list);

#endif
