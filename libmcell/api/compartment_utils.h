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

#ifndef API_COMPARTMENT_UTIL_H
#define API_COMPARTMENT_UTIL_H

#include <string>
#include <map>
#include <vector>

#include "api/geometry_object.h"

namespace MCell {
namespace API {


typedef std::vector<GeometryObjectSet> GeometryObjectSetVector;

void set_parent_and_children_compartments(
    std::vector<std::shared_ptr<API::GeometryObject>>& compartment_objects);

void get_compartment_names(const std::string& bngl_string, std::vector<std::string>& compartments);

} // namespace API
} // namespace MCell

#endif // API_COMPARTMENT_UTIL_H
