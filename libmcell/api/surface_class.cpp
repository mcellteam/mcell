/******************************************************************************
 *
 * Copyright (C) 2020 by
 * The Salk Institute for Biological Studies
 *
 * Use of this source code is governed by an MIT-style
 * license that can be found in the LICENSE file or at
 * https://opensource.org/licenses/MIT.
 *
******************************************************************************/

#include "api/surface_class.h"
#include "api/surface_property.h"

#include <set>

using namespace std;

namespace MCell {
namespace API {

bool SurfaceClass::__eq__(const SurfaceClass& other) const {
  if (!eq_nonarray_attributes(other)) {
    return false;
  }

  // are surface properties the same?
  std::set<SurfaceProperty> s1;
  for (auto& c: properties) {
    s1.insert(*c);
  }
  std::set<SurfaceProperty> s2;
  for (auto& c: other.properties) {
    s2.insert(*c);
  }
  return s1 == s2;
}


} // namespace API
} // namespace MCell
