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

#ifndef API_SURFACE_CLASS_H
#define API_SURFACE_CLASS_H

#include "generated/gen_surface_class.h"
#include "api/common.h"
#include "api/surface_property.h"

namespace MCell {
namespace API {

class SurfaceClass: public GenSurfaceClass {
public:
  SURFACE_CLASS_CTOR()

  void check_semantics() const override {
    GenSurfaceClass::check_semantics(); // does not call further derived classes

    if (properties.empty()) {
      GenSurfaceProperty::check_semantics();
    }
    else {
      // TODO: should each type be specified only once or not?
    }
  }
};

} // namespace API
} // namespace MCell

#endif // API_SURFACE_CLASS_H
