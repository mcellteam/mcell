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
#include "api/api_common.h"

namespace MCell {
namespace API {

enum class RegionType {
  UNKNOWN,
  VOLUME,
  SURFACE
};


class Region: public GenRegion, public std::enable_shared_from_this<Region> {
public:
  REGION_CTOR()

  void postprocess_in_ctor() {
    set_all_custom_attributes_to_default();
  }

  void set_all_custom_attributes_to_default() override {
    region_type = RegionType::UNKNOWN;
  }

  std::shared_ptr<Region> __add__(std::shared_ptr<Region> other) override {
    check_has_compatible_type(other);
    std::shared_ptr<Region> res = std::make_shared<Region>(RegionNodeType::UNION, shared_from_this(), other);
    res->region_type = region_type;
    return res;
  }
  std::shared_ptr<Region> __sub__(std::shared_ptr<Region> other) override {
    check_has_compatible_type(other);
    std::shared_ptr<Region> res = std::make_shared<Region>(RegionNodeType::DIFFERENCE, shared_from_this(), other);
    res->region_type = region_type;
    return res;
  }
  std::shared_ptr<Region> __mul__(std::shared_ptr<Region> other) override {
    check_has_compatible_type(other);
    std::shared_ptr<Region> res = std::make_shared<Region>(RegionNodeType::INTERSECT, shared_from_this(), other);
    res->region_type = region_type;
    return res;
  }

  // added
  void check_has_compatible_type(std::shared_ptr<Region> other) {
    if (region_type != other->region_type) {
      throw RuntimeError("When creating regions, one can only combine regions of identical type (volume or surface), error for:\n" +
          to_str() + "and\n" + other->to_str() +
          "-- end of error message related to combining regions of different type --");
    }
  }

  // set to true in GeometryObject, allows to find out whether it is safe to
  // cast this object to a GeometryObject
  RegionType region_type;
};

} // namespace API
} // namespace MCell

#endif // API_REGION_H
