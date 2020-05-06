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

void CplxInstance::finalize() {
  assert(!mol_instances.empty() && "There must be at least one molecule type");

  // finalize mol instances first
  for (MolInstance& mp: mol_instances) {
    mp.finalize();
  }

  // volume or surface type
  bool vol_type = true;
  for (MolInstance& mp: mol_instances) {
    mp.finalize();
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


bool CplxInstance::matches(const CplxInstance& inst, const bool ignore_orientation) const {
  if (is_simple() && inst.is_simple()) {
    // keep it simple for now...
    assert(mol_instances.size() == 1 && inst.mol_instances.size() == 1);
    return mol_instances[0].matches(inst.mol_instances[0], ignore_orientation);
  }
  else {
    assert(false && "Support for BNG style matching is not implemented yet");
    return false;
  }
}


std::string CplxInstance::to_str(const BNGData& bng_data, bool in_reaction) const {
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
  else if (is_surf() && in_reaction && orientation == ORIENTATION_NONE) {
    ss << ";";
  }
  return ss.str();
}

void CplxInstance::dump(const BNGData& bng_data, std::string ind) const {
  cout << ind << to_str(bng_data);
}

} /* namespace BNG */
