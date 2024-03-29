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

#ifndef API_SURFACE_CLASS_H
#define API_SURFACE_CLASS_H

#include "generated/gen_surface_class.h"
#include "api/api_common.h"
#include "api/surface_property.h"

namespace MCell {
namespace API {

class SurfaceClass: public GenSurfaceClass {
public:
  SURFACE_CLASS_CTOR()

  void postprocess_in_ctor() override {
    set_all_custom_attributes_to_default();
  }

  void set_all_custom_attributes_to_default() override {
    SurfaceProperty::set_all_custom_attributes_to_default();
    species_id = SPECIES_ID_INVALID;
  }

  void check_semantics() const override {
    GenSurfaceClass::check_semantics(); // does not call further derived classes

    if (properties.empty()) {
      SurfaceProperty::check_semantics_custom();
    }
    else {
      // type of used properties must be set
      for (std::shared_ptr<SurfaceProperty> property: properties) {
        property->check_semantics_custom();
      }
    }
  }

  bool __eq__(const SurfaceClass& other) const override;

  // added methods
  bool is_clamp() const {
    return type == SurfacePropertyType::CONCENTRATION_CLAMP ||
        type == SurfacePropertyType::FLUX_CLAMP;
  }

  // simulation engine mapping
  // this is the species_id created for this surface class, not for the affected species
  species_id_t species_id;
};

} // namespace API
} // namespace MCell

#endif // API_SURFACE_CLASS_H
