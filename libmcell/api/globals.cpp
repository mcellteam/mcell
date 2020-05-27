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

#include "api/globals.h"

#include "api/species.h"
#include "generated/gen_constants.h"

namespace MCell {
namespace API {

// we need to have global shared pointers, they are added to Model in its constructor
// having simple global variables would mean that we need to make a copy
// and in model initialization we would be updating the copy (which will not reflect changes to our globals)
std::shared_ptr<Species> AllMolecules(new Species(ALL_MOLECULES.c_str()));
std::shared_ptr<Species> AllVolumeMolecules(new Species(ALL_VOLUME_MOLECULES.c_str()));
std::shared_ptr<Species> AllSurfaceMolecules(new Species(ALL_SURFACE_MOLECULES.c_str()));

} // API
} // MCell


