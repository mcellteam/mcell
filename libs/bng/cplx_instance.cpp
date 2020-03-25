/*
 * complex_species.cpp
 *
 *  Created on: Jan 9, 2020
 *      Author: ahusar
 */
#include <iostream>

#include "ast.h"
#include "bng_engine.h"
#include "cplx_instance.h"
#include "mol_type.h"

using namespace std;

namespace BNG {

void CplxInstance::finalize() {
  // finalize mol instances first
  for (MolInstance& mp: mol_patterns) {
    mp.finalize();
  }

  // volume or surface type
  bool vol_type = true;
  for (MolInstance& mp: mol_patterns) {
    mp.finalize();
    if (mp.is_surf()) {
      vol_type = false;
    }
  }
  if (vol_type) {
    set_flag(CPLX_MOL_FLAG_VOL);
  }
  else {
    set_flag(CPLX_MOL_FLAG_SURF);
  }

  /*
  // orientation
  bool single_orientation = true;
  orientation_t o;
  for (uint i = 0; i < mol_patterns.size(); i++) {
    if (i == 0) {
      o = mol_patterns[i].orientation;
    }
    else {
      if (mol_patterns[i].orientation != o) {
        single_orientation = false;
      }
    }
  }
  set_flag(CPLX_FLAG_HAS_SINGLE_ORIENTATION, single_orientation);
  if (single_orientation) {
    set_flag(CPLX_FLAG_SINGLE_ORIENTATION_IS_UP, o == ORIENTATION_UP);
  }
  */

  // CPLX_FLAG_SINGLE_MOL_NO_COMPONENTS
  bool is_simple = true;
  if (mol_patterns.size() > 1) {
    is_simple = false;
  }
  if (is_simple) {
    for (MolInstance& mp: mol_patterns) {
      if (!mp.component_instances.empty()) {
        is_simple = false;
        break;
      }
    }
  }
  set_flag(CPLX_FLAG_ONE_MOL_NO_COMPONENTS, is_simple);

  set_finalized();
}


bool CplxInstance::matches(const CplxInstance& inst, const bool ignore_orientation) const {
  if (is_simple() && inst.is_simple()) {
    // keep it simple for now...
    assert(mol_patterns.size() == 1 && inst.mol_patterns.size() == 1);
    return mol_patterns[0].matches(inst.mol_patterns[0], ignore_orientation);
  }
  else {
    assert(false && "TODO");
    return false;
  }
}


void CplxInstance::dump(const BNGData& bng_data) const {
  for (size_t i = 0; i < mol_patterns.size(); i++) {
    mol_patterns[i].dump(bng_data);

    if (i != mol_patterns.size() - 1) {
      cout << ".";
    }
  }
}

} /* namespace BNG */
