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

#ifndef API_COMPARTMENT_UTILS_H
#define API_COMPARTMENT_UTILS_H

#include <string>
#include <map>
#include <vector>
#include <memory>

namespace MCell {
namespace API {

class GeometryObject;

void set_parent_and_children_compartments(
    std::vector<std::shared_ptr<API::GeometryObject>>& compartment_objects);

// used also from data model converter
static void get_compartment_names(const std::string& bngl_string, std::vector<std::string>& compartments) {
  size_t i = 0;
  bool in_compartment = false;
  std::string current_name;

  while (i < bngl_string.size()) {
    char c = bngl_string[i];
    if (c == '@') {
      assert(!in_compartment);
      in_compartment = true;
    }
    else if (in_compartment) {
      if ((!isalnum(c) && c != '_')) {
        compartments.push_back(current_name);
        current_name = "";
        in_compartment = false;
      }
      else {
        current_name += c;
      }
    }

    i++;
  }

  if (current_name != "") {
    compartments.push_back(current_name);
  }

}

} // namespace API
} // namespace MCell

#endif // API_COMPARTMENT_UTILS_H
