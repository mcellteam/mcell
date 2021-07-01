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
#include "libs/pybind11/include/pybind11/stl.h"

#ifndef _WIN64
#undef _hypot
#endif

#endif
