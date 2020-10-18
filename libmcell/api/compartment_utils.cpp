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

#include "compartment_utils.h"

namespace MCell {
namespace API {

bool set_parent_and_children_compartments(
    const std::vector<std::shared_ptr<API::GeometryObject>>& compartment_objects,
    GeometryObjectSetVector& intersecting_objects) {

  intersecting_objects.clear();

  auto obj_cp = std::shared_ptr<API::GeometryObject>(nullptr);
  for (auto obj: compartment_objects) {
    release_assert(obj->name == "CP" || obj->name == "EC"); // limited for now
    if (obj->name == "CP") {
      obj_cp = obj;
    }
  }

  auto obj_ec = std::shared_ptr<API::GeometryObject>(nullptr);
  for (auto obj: compartment_objects) {
    if (obj->name == "EC") {
      obj->child_compartments.insert(obj_cp);
      obj_ec = obj;
    }
  }
  if (is_set(obj_cp) && is_set(obj_ec)) {
    obj_cp->parent_compartment = obj_ec;
  }

  return true;
  // FIXME: this is just temporary
  /*bool ok = compute_containement_mapping(world, counted_objects, contained_in_mapping, intersecting_objects);
  if (!ok) {
    return false;
  }*/
}

} // namespace API
} // namespace MCell
