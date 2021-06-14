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

/**
 * This header includes all C++ classes of MCell API.
 * File should be used only from outside of this library to avoid cyclic includes.
 */

#ifndef API_GLOBALS_H
#define API_GLOBALS_H

#include <memory>
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
