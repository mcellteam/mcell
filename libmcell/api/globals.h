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

/**
 * This header includes all C++ classes of MCell API.
 * File should be used only from outside of this library to avoid cyclic includes.
 */

#ifndef API_GLOBALS_H
#define API_GLOBALS_H

#include <memory>
#ifdef _MSC_VER
#undef HAVE_UNISTD_H
#undef HAVE_SYS_TIME_H
#endif
#include "pybind11/include/pybind11/pybind11.h" // make sure we won't include the system header
namespace py = pybind11;

namespace MCell {
namespace API {

class Species;

extern std::shared_ptr<Species> AllMolecules;
extern std::shared_ptr<Species> AllVolumeMolecules;
extern std::shared_ptr<Species> AllSurfaceMolecules;

}
}

#endif // API_GLOBALS_H
