/*
 * nfsim_species_unifier.h
 *
 *  Created on: Aug 19, 2020
 *      Author: ahusar
 */

#ifndef UTILS_BNG_ANALYZER_NFSIM_SPECIES_UNIFIER_H_
#define UTILS_BNG_ANALYZER_NFSIM_SPECIES_UNIFIER_H_

#include <string>

namespace MCell {

class NFSimSpeciesUnifier {
public:
  // TODO
  bool read_species_file(const std::string& input_file) { return true; }
  void print_unified_species() {}
};

}

#endif /* UTILS_BNG_ANALYZER_NFSIM_SPECIES_UNIFIER_H_ */
