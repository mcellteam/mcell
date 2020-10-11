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

#ifndef API_VOLUME_COMPARTMENT_H
#define API_VOLUME_COMPARTMENT_H

#include "generated/gen_volume_compartment.h"
#include "api/common.h"
#include "api/region.h"

namespace MCell {
namespace API {

class VolumeCompartment: public GenVolumeCompartment {
public:
  VOLUME_COMPARTMENT_CTOR()

  // simply subtracts all children from this geometry object
  std::shared_ptr<Region> get_volume_compartment_region() override {

    std::shared_ptr<Region> res = std::dynamic_pointer_cast<Region>(geometry_object);

    for (auto child: child_compartments) {
      res = res->__sub__(std::dynamic_pointer_cast<Region>(child->geometry_object));
    }

    return res;
  }
};

} // namespace API
} // namespace MCell

#endif // API_VOLUME_COMPARTMENT_H
