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

#ifndef API_GEN_GEOMETRY_UTILS_H
#define API_GEN_GEOMETRY_UTILS_H

#include "../api/common.h"

namespace MCell {
namespace API {

class GeometryObject;

namespace geometry_utils {

std::shared_ptr<GeometryObject> create_box(const std::string& name, const float_t edge_length);

} // namespace geometry_utils

void define_pybinding_geometry_utils(py::module& m);
} // namespace API
} // namespace MCell

#endif // API_GEN_GEOMETRY_UTILS_H
