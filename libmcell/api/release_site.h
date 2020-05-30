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

#ifndef API_RELEASE_SITE_H
#define API_RELEASE_SITE_H

#include "generated/gen_release_site.h"
#include "api/common.h"

namespace MCell {
namespace API {

class ReleaseSite: public GenReleaseSite {
public:
  RELEASE_SITE_CTOR()

  void postprocess_in_ctor() override {
    if (is_set(region)) {
      if (shape != Shape::Unset) {
        throw ValueError(S("When ") + NAME_REGION + " is set, "
            "shape must be either unset or set to " + NAME_ENUM_SHAPE + "." + NAME_EV_REGION_EXPR + ".");
      }
      shape = Shape::RegionExpr;
    }
  }

  // actual manual implementation of a semantic check
  void check_semantics() const override {
    GenReleaseSite::check_semantics();
    if (get_num_set(site_diameter, site_radius) > 1) {
      throw ValueError(S("Only either ") + NAME_SITE_DIAMETER + " or " + NAME_SITE_RADIUS + " can be set.");
    }
    if (get_num_set(number_to_release, density) != 1) {
      throw ValueError(S("Exactly one of ") + NAME_NUMBER_TO_RELEASE + " or " + NAME_DENSITY + " must be set.");
    }
  }
};

} // namespace API
} // namespace MCell

#endif // API_RELEASE_SITE_H
