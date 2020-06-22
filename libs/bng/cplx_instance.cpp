/*
 * complex_species.cpp
 *
 *  Created on: Jan 9, 2020
 *      Author: ahusar
 */
#include <iostream>
#include <sstream>

#include "bng/ast.h"
#include "bng/bng_engine.h"
#include "bng/cplx_instance.h"
#include "bng/mol_type.h"

using namespace std;

namespace BNG {

void CplxInstance::finalize_flags() {
  assert(!mol_instances.empty() && "There must be at least one molecule type");

  // finalize mol instances first
  for (MolInstance& mp: mol_instances) {
    mp.finalize_flags();
  }

  // volume or surface type
  bool vol_type = true;
  for (MolInstance& mp: mol_instances) {
    mp.finalize_flags();
    if (mp.is_surf()) {
      vol_type = false;
    }
  }
  if (!vol_type) {
    set_flag(SPECIES_CPLX_MOL_FLAG_SURF);
  }

  // CPLX_FLAG_SINGLE_MOL_NO_COMPONENTS
  bool is_simple = true;
  if (mol_instances.size() > 1) {
    is_simple = false;
  }
  if (is_simple) {
    for (MolInstance& mp: mol_instances) {
      if (!mp.component_instances.empty()) {
        is_simple = false;
        break;
      }
    }
  }
  set_flag(SPECIES_CPLX_FLAG_ONE_MOL_NO_COMPONENTS, is_simple);

  set_finalized();
}


bool CplxInstance::matches_complex_pattern_ignore_orientation(const CplxInstance& pattern) const {
  assert(false && "Support for BNG style matching is not implemented yet");
  return false;
}


bool CplxInstance::matches_complex_fully_ignore_orientation(const CplxInstance& pattern) const {
  assert(false && "Support for BNG style matching is not implemented yet");
  return false;
}


std::string CplxInstance::to_str(const BNGData& bng_data, bool in_surf_reaction) const {
  stringstream ss;
  for (size_t i = 0; i < mol_instances.size(); i++) {
    ss << mol_instances[i].to_str(bng_data);

    if (i != mol_instances.size() - 1) {
      ss << ".";
    }
  }
  if (orientation == ORIENTATION_UP) {
    ss << "'";
  }
  else if (orientation == ORIENTATION_DOWN) {
    ss << ",";
  }
  else if (in_surf_reaction && orientation == ORIENTATION_NONE) {
    ss << ";";
  }
  return ss.str();
}

void CplxInstance::dump(const BNGData& bng_data, const bool for_diff, const std::string ind) const {
  if (!for_diff) {
    cout << ind << to_str(bng_data);
  }
  else {
    cout << ind << "orientation: " << orientation << "\n";
    cout << ind << "mol_instances:\n";
    for (size_t i = 0; i < mol_instances.size(); i++) {
      cout << ind << i << ":\n";
      mol_instances[i].dump(bng_data, true, false, ind + "  ");
    }
  }
}

} /* namespace BNG */
