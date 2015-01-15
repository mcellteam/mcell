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
* ****************************************************************************/

#include <stdlib.h>
#include <string.h>

#include "mcell_misc.h"
#include "mcell_dyngeom.h"

int mcell_add_dynamic_geometry_file(char const *dynamic_geometry_filepath,
                                    struct volume *state) {
  char *dynamic_geometry_filename =
      mcell_find_include_file(dynamic_geometry_filepath, state->curr_file);
  state->dynamic_geometry_filename = dynamic_geometry_filename;
  return 0;
}
