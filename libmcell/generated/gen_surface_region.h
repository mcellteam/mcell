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

#ifndef API_GEN_SURFACE_REGION_H
#define API_GEN_SURFACE_REGION_H

#include "../api/common.h"

namespace MCell {
namespace API {

#define SURFACE_REGION_CTOR() \
    SurfaceRegion( \
        const std::string& name_, \
        const std::vector<int> element_connections_ \
    ) { \
      class_name = "SurfaceRegion"; \
      name = name_; \
      element_connections = element_connections_; \
      postprocess_in_ctor();\
      check_semantics();\
    }

class GenSurfaceRegion: public BaseDataClass {
public:
  void postprocess_in_ctor() override {}
  void check_semantics() const override;
  std::string to_str(const std::string ind="") const override;

  // --- attributes ---
  std::vector<int> element_connections;
  virtual void set_element_connections(const std::vector<int> new_element_connections_) {
    element_connections = new_element_connections_;
  }
  virtual std::vector<int> get_element_connections() const {
    return element_connections;
  }

  // --- methods ---
}; // GenSurfaceRegion

class SurfaceRegion;
py::class_<SurfaceRegion> define_pybinding_SurfaceRegion(py::module& m);
} // namespace API
} // namespace MCell

#endif // API_GEN_SURFACE_REGION_H
