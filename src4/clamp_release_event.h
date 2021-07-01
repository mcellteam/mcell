/******************************************************************************
 *
 * Copyright (C) 2020 by
 * The Salk Institute for Biological Studies and
 * Pittsburgh Supercomputing Center, Carnegie Mellon University
 *
 * Use of this source code is governed by an MIT-style
 * license that can be found in the LICENSE file or at
 * https://opensource.org/licenses/MIT.
 *
******************************************************************************/

#ifndef SRC4_CLAMP_RELEASE_EVENT_H_
#define SRC4_CLAMP_RELEASE_EVENT_H_

#include <vector>

#include "base_event.h"

namespace MCell {

enum class ClampType {
  INVALID,
  CONCENTRATION,
  FLUX
};

/**
 * Molecules are released at:
 * 1) concentration-clamped surfaces to maintain the desired concentration
 * or
 * 2) flux-clamped surfaces to maintain the desired influx
 *
 * Separate event from releases because it must be run after all releases and
 * viz/count outputs before diffusion.
 */
class ClampReleaseEvent: public BaseEvent {
public:
  ClampReleaseEvent(World* world_) :
    BaseEvent(EVENT_TYPE_INDEX_CLAMP_RELEASE),
    type(ClampType::INVALID),
    species_id(SPECIES_ID_INVALID),
    surf_class_species_id(SPECIES_ID_INVALID),
    concentration(FLT_INVALID),
    orientation(ORIENTATION_NONE),
    scaling_factor(FLT_INVALID),
    world(world_) {
  }

  virtual ~ClampReleaseEvent() {}

  void step() override;

  void dump(const std::string indent) const override;

  void to_data_model(Json::Value& mcell_node) const override;

  void update_cumm_areas_and_scaling();
public:
  ClampType type; // used only when converting to data model
  species_id_t species_id;
  species_id_t surf_class_species_id;
  double concentration;
  orientation_t orientation;
  double scaling_factor;

  std::vector<CummAreaPWallIndexPair> cumm_area_and_pwall_index_pairs;

private:
  World* world;
};

} // namespace mcell


#endif // SRC4_CLAMP_RELEASE_EVENT_H_
