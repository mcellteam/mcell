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

// this is an empty implementation of some functions from libmcell
// MCell executable (not mcell.so) references but not uses Callbacks class and other
// functions and we do not want to link all the Python libraries it needs

#include "api/reaction_info.h"
#include "api/callbacks.h"
#include "api/mol_wall_hit_info.h"
#include "api/checkpoint_signals.h"

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
