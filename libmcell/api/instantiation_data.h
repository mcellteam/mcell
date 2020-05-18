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

#ifndef API_INSTANTIATION_DATA_H
#define API_INSTANTIATION_DATA_H

#include "../generated/gen_instantiation_data.h"
#include "../api/common.h"

namespace MCell {
namespace API {

class InstantiationData: public GenInstantiationData {
public:

  // from generated template
  void instantiate_release_site(std::shared_ptr<ReleaseSite> s) override {}
  std::shared_ptr<ReleaseSite> find_release_site(const std::string& name) override {return nullptr;}
  void instantiate_geometry_object(std::shared_ptr<GeometryObject> o, const std::string& name = "") override {}
  void find_geometry_object(const std::string& name) override {}

  // added manually
  void dump() const;
};

} // namespace API
} // namespace MCell

#endif // API_INSTANTIATION_DATA_H
