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

#ifndef API_MOLECULE_RELEASE_INFO_H
#define API_MOLECULE_RELEASE_INFO_H

#include "generated/gen_molecule_release_info.h"
#include "api/common.h"

namespace MCell {
namespace API {

class MoleculeReleaseInfo: public GenMoleculeReleaseInfo {
public:
  MOLECULE_RELEASE_INFO_CTOR()

  void check_semantics() const override {
    if (location.size() != 3) {
      throw ValueError(S("Value of ") + NAME_LOCATION + " must be a triplet of floats.");
    }
  }
};

} // namespace API
} // namespace MCell

#endif // API_MOLECULE_RELEASE_INFO_H
