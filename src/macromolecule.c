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

#include "config.h"

#include "macromolecule.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>


/*******************************************************************************
 macro_subunit_index:

    Given a macromolecule subunit, find its index within the complex.

    In:  struct abstract_molecule const *subunit - the subunit whose index
                   we'd like to locate
    Out: the subunit's index, or -1 if the molecule is not a subunit.  In
                   normal program operation, -1 should only ever be returned
                   for molecules which represent the complex itself.
********************************************************************************/
int macro_subunit_index(struct abstract_molecule const *subunit) {
  struct abstract_molecule *const *c = subunit->cmplx;
  assert(c != NULL);

  struct complex_species *s = (struct complex_species *)c[0]->properties;
  assert(s->base.flags & IS_COMPLEX);

  for (int i = 0; i < s->num_subunits; ++i)
    if (c[i + 1] == subunit)
      return i;
  return -1;
}
