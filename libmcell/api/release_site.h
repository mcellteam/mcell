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
#include "api/molecule_release_info.h"
#include "api/common.h"
#include "api/complex.h"

namespace MCell {
namespace API {

class ReleaseSite: public GenReleaseSite {
public:
  // empty ctor to be used internally
  ReleaseSite() {
    set_all_attributes_as_default_or_unset();
  }

  RELEASE_SITE_CTOR()

  void postprocess_in_ctor() override {
    if (is_set(region) && shape != Shape::REGION_EXPR) {
      if (shape != Shape::UNSET) {
        throw ValueError(S("When ") + NAME_REGION + " is set, "
            "shape must be either unset or set to " + NAME_ENUM_SHAPE + "." + NAME_EV_REGION_EXPR + ".");
      }
      shape = Shape::REGION_EXPR;
    }

    if (is_set(complex) && is_set(complex->compartment_name) && shape != Shape::COMPARTMENT) {
      if (shape != Shape::UNSET) {
        throw ValueError(S("When ") + NAME_COMPARTMENT_NAME + " is set, "
            "shape must be either unset or set to " + NAME_ENUM_SHAPE + "." + NAME_EV_COMPARTMENT + ".");
      }
      shape = Shape::COMPARTMENT;
    }

    if (is_set(molecule_list) && shape != Shape::LIST) {
      if (shape != Shape::UNSET) {
        throw ValueError(S("When ") + NAME_MOLECULE_LIST + " is set, "
            "shape must be either unset or set to " + NAME_ENUM_SHAPE + "." + NAME_EV_LIST + ".");
      }
      shape = Shape::LIST;
    }

    if (is_set(location) && shape != Shape::SPHERICAL) {
      if (shape != Shape::UNSET) {
        throw ValueError(S("When ") + NAME_LOCATION + " is set, "
            "shape must be either unset or set to " + NAME_ENUM_SHAPE + "." + NAME_EV_SPHERICAL + ".");
      }
      shape = Shape::LIST;
    }
  }

  void check_semantics() const override {
    // TODO: add additional checks, e.g. that site_diameter must can be set only for SPHERICAL, etc.
    GenReleaseSite::check_semantics();
    if (release_time < 0) {
      throw ValueError(S("Value of ") + NAME_RELEASE_TIME + " must not be smaller than 0.");
    }

    if (get_num_set(site_diameter, site_radius) > 1) {
      throw ValueError(S("Only either ") + NAME_SITE_DIAMETER + " or " + NAME_SITE_RADIUS + " can be set.");
    }

    if (get_num_set(number_to_release, density, concentration, molecule_list) != 1) {
      throw ValueError(
          S("Exactly one of ") + NAME_NUMBER_TO_RELEASE + ", " + NAME_DENSITY + ", " + NAME_CONCENTRATION +
          " or " + NAME_MOLECULE_LIST + " must be set.");
    }

    if (get_num_set(complex, molecule_list) != 1) {
      throw ValueError(
          S("Exactly one of ") + NAME_COMPLEX + " or " + NAME_MOLECULE_LIST + " must be set.");
    }

    if (shape == Shape::COMPARTMENT &&
        BNG::get_in_or_out_compartment_id(complex->compartment_name) != BNG::COMPARTMENT_ID_INVALID) {
      throw ValueError(
          S(NAME_CLASS_RELEASE_SITE) + " must not use compartment class name " + complex->compartment_name + ".");
    }

    if ( (shape == Shape::COMPARTMENT && get_num_set(region, molecule_list, location) != 0) &&
         (shape != Shape::COMPARTMENT && get_num_set(region, molecule_list, location) != 1)
        ) {
      throw ValueError(
          S("Either compartment or exactly one of ") + NAME_REGION + ", " + NAME_MOLECULE_LIST +
          " or " + NAME_LOCATION + " must be set.");
    }
  }
};

} // namespace API
} // namespace MCell

#endif // API_RELEASE_SITE_H
