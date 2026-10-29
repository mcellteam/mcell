/******************************************************************************
 *
 * Copyright (C) 2026 by
 * The Salk Institute for Biological Studies
 *
 * Use of this source code is governed by an MIT-style
 * license that can be found in the LICENSE file or at
 * https://opensource.org/licenses/MIT.
 *
******************************************************************************/

#ifndef SRC4_PREDICATES_H_
#define SRC4_PREDICATES_H_

namespace MCell {
namespace Predicates {

// Shewchuk's adaptive exact orientation tests (predicates.cpp). The sign of the
// result is always exact: zero only when the points are exactly collinear or
// coplanar. The magnitude is approximate. Points are arrays of 2 or 3 doubles.

// positive if pa, pb and pc occur in counterclockwise order, negative if
// clockwise, zero if collinear; roughly twice the signed area of the triangle
double orient2d(const double* pa, const double* pb, const double* pc);

// positive if pd lies below the plane through pa, pb and pc ("below" means
// that pa, pb and pc appear counterclockwise when viewed from above), negative
// if above, zero if coplanar; roughly six times the signed volume of the
// tetrahedron
double orient3d(const double* pa, const double* pb, const double* pc, const double* pd);

} // namespace Predicates
} // namespace MCell

#endif // SRC4_PREDICATES_H_
