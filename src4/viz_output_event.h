/******************************************************************************
 *
 * Copyright (C) 2019 by
 * The Salk Institute for Biological Studies and
 * Pittsburgh Supercomputing Center, Carnegie Mellon University
 *
 * Use of this source code is governed by an MIT-style
 * license that can be found in the LICENSE file or at
 * https://opensource.org/licenses/MIT.
 *
******************************************************************************/

#ifndef SRC4_VIZ_OUTPUT_EVENT_H_
#define SRC4_VIZ_OUTPUT_EVENT_H_

#include "base_event.h"
#include "mcell_structs_shared.h"

namespace MCell {

class Partition;
class Molecule;

/**
 * Dumps world state either in a textual or cellblender format.
 */
class VizOutputEvent: public BaseEvent {
public:
  VizOutputEvent(World* world_)
    : BaseEvent(EVENT_TYPE_INDEX_VIZ_OUTPUT),
      viz_mode(NO_VIZ_MODE),
      visualize_all_species(false),
      world(world_) {
  }
  virtual ~VizOutputEvent() {}

  void step() override;
  void dump(const std::string ind = "") const override;
  void to_data_model(Json::Value& mcell_node) const override;

  // DiffuseReactEvent must execute only up to this event
  bool is_barrier() const override { return true; }

  bool should_visualize_all_species() const;

  static std::string iterations_to_string(const uint64_t current_iterations, const uint64_t total_iterations);

  viz_mode_t viz_mode;
  std::string file_prefix_name;
  bool visualize_all_species;
  uint_set<species_id_t> species_ids_to_visualize;
  World* world;

private:
  void compute_where_and_norm(
      const Partition& p, const Molecule& m,
      Vec3& where, Vec3& norm
  );

  FILE* create_and_open_output_file_name();
  void output_ascii_molecules();
  void output_cellblender_molecules();
};

} // namespace mcell

#endif // SRC4_VIZ_OUTPUT_EVENT_H_
