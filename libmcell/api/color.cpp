/******************************************************************************
 *
 * Copyright (C) 2021 by
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

#include "api/color.h"

namespace MCell {
namespace API {

void Color::postprocess_in_ctor() {
  if (is_set(red) && is_set(green) && is_set(blue) && is_set(alpha)) {
    check_component_range(red, NAME_RED);
    check_component_range(green, NAME_GREEN);
    check_component_range(blue, NAME_BLUE);
    check_component_range(alpha, NAME_ALPHA);

    components_to_rgba();
  }
  else {
    bool any_component_set = is_set(red) || is_set(green) || is_set(blue);
    if (any_component_set) {
      throw ValueError(S("Either all individual color components must be set or none in initialization of ") +
          NAME_CLASS_COLOR + ".");
    }

    rgba_to_components();
  }
}


void Color::check_component_range(const double value, const char* name) {
  if (value < 0 || value > 1) {
    throw ValueError(S("Value of color component ") + name + " must be within 0..1.");
  }
}


void Color::components_to_rgba() {
  const int MAX = 255;
  rgba =
      (((uint)(red * MAX) & 0xFF) << 24) |
      (((uint)(green * MAX) & 0xFF) << 16) |
      (((uint)(blue * MAX) & 0xFF) << 8) |
      (((uint)(alpha * MAX)) & 0xFF);
}


void Color::rgba_to_components() {
  const double MAX = 255.0;
  red = (((uint)rgba >> 24) & 0xFF) / MAX;
  green = (((uint)rgba >> 16) & 0xFF) / MAX;
  blue = (((uint)rgba >> 8) & 0xFF) / MAX;
  alpha = ((uint)rgba & 0xFF) / MAX;
}

} // namespace API
} // namespace MCell
