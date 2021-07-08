/******************************************************************************
 *
 * Copyright (C) 2006-2014 by
 * The Salk Institute for Biological Studies and
 * Pittsburgh Supercomputing Center, Carnegie Mellon University
 *
 * Use of this source code is governed by an MIT-style
 * license that can be found in the LICENSE file or at
 * https://opensource.org/licenses/MIT.
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
