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

// this is an empty implementation of some functions from libmcell
// MCell executable (not mcell.so) references but not uses Callbacks class and other
// functions and we do not want to link all the Python libraries it needs

#include "api/reaction_info.h"
#include "api/callbacks.h"
#include "api/mol_wall_hit_info.h"
#include "api/checkpoint_signals.h"

#if PYTHON_VERSION == 39 && !defined(_MSC_VER)

extern "C" {
// including pybind11 with Python 3.9 requires this symbol, however, must not be used
void __attribute__((weak)) _Py_Dealloc(PyObject*) {
  release_assert("must not be called");
}

}
#endif

namespace MCell {
namespace API {

void save_checkpoint_func(double, MCell::API::CheckpointSaveEventContext) {
  release_assert("must not be called");
}

Callbacks::Callbacks(Model*) {
  // empty
}

void Callbacks::do_mol_wall_hit_callbacks(std::shared_ptr<MolWallHitInfo>) {
  release_assert("must not be called");
}

// we also need some implementations for MolWallHitInfo
bool GenMolWallHitInfo::__eq__(const MolWallHitInfo&) const {
  release_assert("must not be called");
  return false;
}

bool GenMolWallHitInfo::eq_nonarray_attributes(const MolWallHitInfo&, const bool) const {
  release_assert("must not be called");
  return false;
}

void Callbacks::do_rxn_callback(std::shared_ptr<ReactionInfo>) {
  release_assert("must not be called");
}

bool GenReactionInfo::__eq__(const ReactionInfo&) const {
  release_assert("must not be called");
  return false;
}

bool GenReactionInfo::eq_nonarray_attributes(const ReactionInfo&, const bool) const {
  release_assert("must not be called");
  return false;
}

} /* namespace API */
} /* namespace MCell */
