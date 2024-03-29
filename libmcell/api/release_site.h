/******************************************************************************
 *
 * Copyright (C) 2020 by
 * The Salk Institute for Biological Studies
 *
 * Use of this source code is governed by an MIT-style
 * license that can be found in the LICENSE file or at
 * https://opensource.org/licenses/MIT.
 *
******************************************************************************/

#ifndef API_RELEASE_SITE_H
#define API_RELEASE_SITE_H

#include "generated/gen_release_site.h"
#include "api/molecule_release_info.h"
#include "api/api_common.h"
#include "api/complex.h"
#include "bng/bng_defines.h"

namespace MCell {
namespace API {

class ReleaseSite: public GenReleaseSite {
public:
  RELEASE_SITE_CTOR()

  void postprocess_in_ctor() override {
    if (is_set(region) && shape != Shape::REGION_EXPR) {
      if (shape != Shape::UNSET) {
        throw ValueError(S("When ") + NAME_REGION + " is set, "
            "shape must be either unset or set to " + NAME_ENUM_SHAPE + "." + NAME_EV_REGION_EXPR + ".");
      }
      shape = Shape::REGION_EXPR;
    }

    if (is_set(complex) && is_set(complex->get_primary_compartment_name()) &&
        shape != Shape::COMPARTMENT && shape != Shape::REGION_EXPR) {
      if (shape != Shape::UNSET) {
        throw ValueError(S("When ") + NAME_COMPARTMENT_NAME + " is set, "
            "shape must be either unset or set to " + NAME_ENUM_SHAPE + "." + NAME_EV_COMPARTMENT + " or " +
            NAME_ENUM_SHAPE + "." + NAME_EV_REGION_EXPR + ".");
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
      shape = Shape::SPHERICAL;
    }

    if (release_probability < 0 || release_probability > 1) {
      throw ValueError(S("Parameter ") + NAME_RELEASE_PROBABILITY + " must be >= 0 and <= 1.");
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

    if (shape != Shape::REGION_EXPR && get_num_set(number_to_release, density, concentration) == 1 &&
        (number_to_release < 0 || density < 0 || concentration < 0)) {
      const char* attribute_name;
      if (is_set(number_to_release)) {
        attribute_name = NAME_NUMBER_TO_RELEASE;
      }
      else if (is_set(density)) {
        attribute_name = NAME_DENSITY;
      }
      else if (is_set(concentration)) {
        attribute_name = NAME_CONCENTRATION;
      }
      else {
        release_assert(false);
      }
      throw ValueError(
          S("Negative release value of ") + attribute_name + " may be set only when " + NAME_REGION + " is set.");
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

    if (shape == Shape::SPHERICAL && !is_set(location)) {
      throw ValueError(
          S("When ") + NAME_SHAPE + " is set to " + NAME_ENUM_SHAPE + "." + NAME_EV_LIST +
          " " + NAME_LOCATION + " must be set.");
    }

    if (is_set(location) && location.size() != 3) {
      throw ValueError(S("Argument ") + NAME_LOCATION + " must be a list containing 3 floats.");
    }
  }
};

} // namespace API
} // namespace MCell

#endif // API_RELEASE_SITE_H
