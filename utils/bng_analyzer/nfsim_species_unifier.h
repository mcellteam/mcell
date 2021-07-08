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

  std::map<BNG::species_id_t, double> counts_per_unique_species;
  BNG::BNGConfig bng_config;
  BNG::BNGEngine bng_engine;
};

} // namespace MCell

#endif /* UTILS_BNG_ANALYZER_NFSIM_SPECIES_UNIFIER_H_ */
