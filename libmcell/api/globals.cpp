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

#include "api/globals.h"

#include "api/species.h"
#include "generated/gen_constants.h"

namespace MCell {
namespace API {

// we need to have global shared pointers, they are added to Model in its constructor
// having simple global variables would mean that we need to make a copy
// and in model initialization we would be updating the copy (which will not reflect changes to our globals)
std::shared_ptr<Species> AllMolecules(new Species(ALL_MOLECULES, FLT_UNSET, 0));
std::shared_ptr<Species> AllVolumeMolecules(new Species(ALL_VOLUME_MOLECULES, FLT_UNSET, 0));
std::shared_ptr<Species> AllSurfaceMolecules(new Species(ALL_SURFACE_MOLECULES, 0, FLT_UNSET));

} // API
} // MCell


