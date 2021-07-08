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

#ifndef API_SURFACE_REGION_H
#define API_SURFACE_REGION_H

#include "generated/gen_surface_region.h"
#include "api/api_common.h"
#include "defines.h"

namespace MCell {
namespace API {

class GeometryObject;

class SurfaceRegion: public GenSurfaceRegion {
public:
  SURFACE_REGION_CTOR()

  void postprocess_in_ctor() {
    set_all_custom_attributes_to_default();
  }

  void set_all_custom_attributes_to_default() override {
    Region::set_all_custom_attributes_to_default();

    region_type = RegionType::SURFACE;

    region_id = REGION_ID_INVALID;
    node_type = RegionNodeType::LEAF_SURFACE_REGION;
    parent = nullptr;
  }

  GeometryObject* parent;

  // simulation engine mapping
  region_id_t region_id;
};

} // namespace API
} // namespace MCell

#endif // API_SURFACE_REGION_H
