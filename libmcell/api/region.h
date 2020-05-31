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

#ifndef API_REGION_H
#define API_REGION_H

#include "generated/gen_region.h"
#include "api/common.h"

namespace MCell {
namespace API {

class Region: public GenRegion, public std::enable_shared_from_this<Region> {
public:
  REGION_CTOR()

  std::shared_ptr<Region> __add__(std::shared_ptr<Region> other) override {
    return std::make_shared<Region>(RegionNodeType::UNION, shared_from_this(), other);
  }
  std::shared_ptr<Region> __sub__(std::shared_ptr<Region> other) override {
    return std::make_shared<Region>(RegionNodeType::DIFFERENCE, shared_from_this(), other);
  }
  std::shared_ptr<Region> __mul__(std::shared_ptr<Region> other) override {
    return std::make_shared<Region>(RegionNodeType::INTERSECT, shared_from_this(), other);
  }
};

} // namespace API
} // namespace MCell

#endif // API_REGION_H
