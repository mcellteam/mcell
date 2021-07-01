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

#include "molecule.h"

#include "world.h"
#include "partition.h"

namespace MCell {
namespace API {

void Molecule::remove() {
  check_initialization();

  Partition& p = world->get_partition(PARTITION_ID_INITIAL);

  if (!p.does_molecule_exist(id)) {
    throw RuntimeError("Molecule with id " + std::to_string(id) + " does not exist anymore.");
  }

  // set that this molecule is defunct
  world->get_partition(PARTITION_ID_INITIAL).get_m(id).set_is_defunct();
}

}
}

