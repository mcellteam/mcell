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

#ifndef API_SURFACE_PROPERTY_H
#define API_SURFACE_PROPERTY_H

#include "generated/gen_surface_property.h"
#include "api/api_common.h"
#include "bng/bng_defines.h"
#include "api/complex.h"

namespace MCell {
namespace API {

class SurfaceProperty: public GenSurfaceProperty {
public:
  SURFACE_PROPERTY_CTOR()

  void postprocess_in_ctor() override {
    set_all_custom_attributes_to_default();
  }

  void set_all_custom_attributes_to_default() override {
    rxn_rule_id = BNG::RXN_RULE_ID_INVALID;
  }

  void check_semantics() const override {
    GenSurfaceProperty::check_semantics();
    // all checks are done in SurfaceClass::check_semantics that calls 
    // check_semantics_custom
  }
  
  void check_semantics_custom() const {
    GenSurfaceProperty::check_semantics();

    if (type == SurfacePropertyType::UNSET) {
      throw ValueError(S("Attribute '") + NAME_TYPE + "' of " +
          NAME_CLASS_SURFACE_PROPERTY + " objects contained in " + NAME_CLASS_SURFACE_CLASS + " must be set.");
    }
    if (!is_set(affected_complex_pattern) && type != SurfacePropertyType::REACTIVE) {
      throw ValueError(S("Attribute '") + NAME_AFFECTED_COMPLEX_PATTERN + "' of " +
      NAME_CLASS_SURFACE_PROPERTY + " objects contained in " + NAME_CLASS_SURFACE_CLASS + " must be set when " +
      NAME_TYPE + " is different from " + NAME_ENUM_SURFACE_PROPERTY_TYPE + "." + NAME_EV_REACTIVE + ".");
    }
    if (is_set(affected_complex_pattern) && type == SurfacePropertyType::REACTIVE) {
      throw ValueError(S("Attribute '") + NAME_AFFECTED_COMPLEX_PATTERN + "' of " +
      NAME_CLASS_SURFACE_PROPERTY + " objects contained in " + NAME_CLASS_SURFACE_CLASS + " must not be set when " +
      NAME_TYPE + " is " + NAME_ENUM_SURFACE_PROPERTY_TYPE + "." + NAME_EV_REACTIVE + ".");
    }
    if (is_set(affected_complex_pattern) && is_set(affected_complex_pattern->compartment_name)) {
      throw ValueError(S("Attribute '") + NAME_AFFECTED_COMPLEX_PATTERN + "' of " +
          NAME_CLASS_SURFACE_PROPERTY + " must not have a compartment specified.");
    }
    if ((type == SurfacePropertyType::CONCENTRATION_CLAMP || type == SurfacePropertyType::FLUX_CLAMP) &&
        !is_set(concentration)) {
      throw ValueError(S("When ") + NAME_TYPE + " in " + NAME_CLASS_SURFACE_PROPERTY + " is " +
          NAME_ENUM_SURFACE_PROPERTY_TYPE + "." + NAME_EV_CONCENTRATION_CLAMP + " or " +
          NAME_ENUM_SURFACE_PROPERTY_TYPE + "." + NAME_EV_FLUX_CLAMP +
          " then " + NAME_CONCENTRATION + " must be set.");
    }
  }

  // needed when defining a set of SurfaceProperty(s)
  bool operator < (const SurfaceProperty& other) const {
    if (type != other.type) {
      return type < other.type;
    }
    if (concentration != other.concentration) {
      return concentration < other.concentration;
    }

    return affected_complex_pattern->get_canonical_name() <
        other.affected_complex_pattern->get_canonical_name();
  }

  BNG::rxn_rule_id_t rxn_rule_id;
};

} // namespace API
} // namespace MCell

#endif // API_SURFACE_PROPERTY_H
