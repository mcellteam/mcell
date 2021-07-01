/******************************************************************************
 *
 * Copyright (C) 2021 by
 * The Salk Institute for Biological Studies and
 * Pittsburgh Supercomputing Center, Carnegie Mellon University
 *
 * Use of this source code is governed by an MIT-style
 * license that can be found in the LICENSE file or at
 * https://opensource.org/licenses/MIT.
 *
******************************************************************************/

#ifndef SRC4_BNGL_EXPORTER_H_
#define SRC4_BNGL_EXPORTER_H_

#include <string>

namespace MCell {

namespace API {
enum class BNGSimulationMethod;
}

class World;

class BNGLExporter {
public:
  // returns empty string if everything went well,
  // nonempty string with error message
  std::string export_to_bngl(
      World* world_,
      const std::string& file_name,
      const API::BNGSimulationMethod simulation_method);

private:
  void clear_temporaries();

  std::string set_compartment_volumes_and_areas();

  std::string export_releases_to_bngl_seed_species(
      std::ostream& parameters,
      std::ostream& seed_species) const;

  std::string export_counts_to_bngl_observables(
      std::ostream& observables) const;

  void generate_simulation_action(
      std::ostream& out, const API::BNGSimulationMethod simulation_method) const;

  World* world;
  bool nfsim_export;
};


} // namespace MCell

#endif // SRC4_BNGL_EXPORTER_H_
