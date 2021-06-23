/******************************************************************************
 *
 * Copyright (C) 2021 by
 * The Salk Institute for Biological Studies and
 * Pittsburgh Supercomputing Center, Carnegie Mellon University
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
      const World* world_,
      const std::string& file_name,
      const API::BNGSimulationMethod simulation_method);

private:
  void clear_temporaries();

  std::string export_releases_to_bngl_seed_species(
      std::ostream& parameters,
      std::ostream& seed_species) const;

  std::string export_counts_to_bngl_observables(
      std::ostream& observables) const;

  void generate_simulation_action(
      std::ostream& out, const API::BNGSimulationMethod simulation_method) const;

  const World* world;
  bool nfsim_export;
};


} // namespace MCell

#endif // SRC4_BNGL_EXPORTER_H_
