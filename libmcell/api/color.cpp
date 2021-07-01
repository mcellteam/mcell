/******************************************************************************
 *
 * Copyright (C) 2021 by
 * The Salk Institute for Biological Studies
 *
 * Use of this source code is governed by an MIT-style
 * license that can be found in the LICENSE file or at
 * https://opensource.org/licenses/MIT.
 *
******************************************************************************/

#include "api/color.h"
#include "geometry.h"

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
  Geometry::rgba_to_components(rgba, red, green, blue, alpha);
}

} // namespace API
} // namespace MCell
