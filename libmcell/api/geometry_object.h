/******************************************************************************
 *
 * Copyright (C) 2020 by
 * The Salk Institute for Biological Studies and
 * Pittsburgh Supercomputing Center, Carnegie Mellon University
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

#ifndef API_GEOMETRY_OBJECT_H
#define API_GEOMETRY_OBJECT_H

#include "../generated/gen_geometry_object.h"
#include "../api/common.h"

namespace MCell {
namespace API {

class GeometryObject: public GenGeometryObject {
public:
  GEOMETRY_OBJECT_CTOR()

  void check_semantics() const override {
    for (auto& v: vertex_list) {
      if (v.size() != 3) {
        throw ValueError(
            "Each item in the 'vertex_list' argument must be a triplet of floats, error for " +
            vec_nonptr_to_str(v) + ".");
      }
    }

    for (auto& e: element_connections) {
      if (e.size() != 3) {
        throw ValueError(
            "Each item in the 'element_connections' argument must be a triplet of integers, error for " +
            vec_nonptr_to_str(e) + ".");
      }
    }
  }
};

} // namespace API
} // namespace MCell

#endif // API_GEOMETRY_OBJECT_H
