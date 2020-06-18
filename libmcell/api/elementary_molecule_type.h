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

#ifndef API_ELEMENTARY_MOLECULE_TYPE_H
#define API_ELEMENTARY_MOLECULE_TYPE_H

#include "generated/gen_elementary_molecule_type.h"
#include "api/common.h"
#include "api/elementary_molecule_instance.h"

namespace MCell {
namespace API {

class ElementaryMoleculeType:
    public GenElementaryMoleculeType, public std::enable_shared_from_this<ElementaryMoleculeType> {
public:
  ELEMENTARY_MOLECULE_TYPE_CTOR()

  std::shared_ptr<ElementaryMoleculeInstance> inst(const std::vector<std::shared_ptr<ComponentInstance>> components) override {
    return std::make_shared<ElementaryMoleculeInstance>( shared_from_this(), components);
  }
};

} // namespace API
} // namespace MCell

#endif // API_ELEMENTARY_MOLECULE_TYPE_H
