/******************************************************************************
 *
 * Copyright (C) 2020 by
 * The Salk Institute for Biological Studies
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

#ifndef API_INITIAL_SURFACE_RELEASE_H
#define API_INITIAL_SURFACE_RELEASE_H

#include "generated/gen_initial_surface_release.h"
#include "api/common.h"
#include "api/complex.h"

namespace MCell {
namespace API {

class InitialSurfaceRelease: public GenInitialSurfaceRelease {
public:
  INITIAL_SURFACE_RELEASE_CTOR()

  void check_semantics() const override {
    if (get_num_set(number_to_release, density) != 1) {
      throw ValueError(S("Exactly one of ") + NAME_NUMBER_TO_RELEASE + " and " +
          NAME_DENSITY + " must be set.");
    }

    if (is_set(complex->compartment_name)) {
      throw ValueError(S("Compartment of the complex to be released in ") +
          NAME_CLASS_INITIAL_SURFACE_RELEASE + " must not be set.");
    }
  }

};

} // namespace API
} // namespace MCell

#endif // API_INITIAL_SURFACE_RELEASE_H
