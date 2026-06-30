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

// Windows build needs a fix otherwise this compilation error occurs:
// msys/mingw64/include/c++/10.2.0/cmath:1121:11: error: 'hypot' has not been declared in '::'

#ifndef API_PYBIND_STL_INCLUDE
#define API_PYBIND_STL_INCLUDE

#ifdef _WIN64
// fix for _hypot compilation issue
#define _hypot hypot
#include <cmath>
#endif
#ifdef _MSC_VER
#undef HAVE_UNISTD_H
#undef HAVE_SYS_TIME_H
#endif
// nanobind STL type casters (migrated from pybind11/stl.h). Opt-in per type.
#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/shared_ptr.h>
#include <nanobind/stl/vector.h>
#include <nanobind/stl/map.h>
#include <nanobind/stl/set.h>
#include <nanobind/stl/function.h>
#include <nanobind/stl/bind_vector.h>

#ifndef _WIN64
#undef _hypot
#endif

#endif
