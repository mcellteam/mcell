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

#ifndef UTILS_BNG_ANALYZER_NFSIM_SPECIES_UNIFIER_H_
#define UTILS_BNG_ANALYZER_NFSIM_SPECIES_UNIFIER_H_

#include <string>
#include <map>

#include "bng/bng.h"

namespace MCell {

class NFSimSpeciesUnifier {
public:
  NFSimSpeciesUnifier()
    : bng_engine(bng_config) {
  }
  bool read_input_file(const std::string& input_file, const bool is_dat);
  bool print_unified_species(const std::string& out_file = "");

  std::map<species_id_t, double> counts_per_unique_species;
  BNG::BNGConfig bng_config;
  BNG::BNGEngine bng_engine;
};

} // namespace MCell

#endif /* UTILS_BNG_ANALYZER_NFSIM_SPECIES_UNIFIER_H_ */
