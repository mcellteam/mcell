/******************************************************************************
 *
 * Copyright (C) 2021 by
 * The Salk Institute for Biological Studies and
 * Pittsburgh Supercomputing Center, Carnegie Mellon University
 *
 * Use of this source code is governed by an MIT-style
 * license that can be found in the LICENSE file or at
 * https://opensource.org/licenses/MIT.
 *
******************************************************************************/

#include <iostream>
#include <algorithm>

#include "mol_order_shuffle_event.h"
#include "world.h"
#include "partition.h"

using namespace std;

namespace MCell {

void MolOrderShuffleEvent::dump(const string ind) const {
  cout << ind << "Mol order shuffle event:\n";
  string ind2 = ind + "  ";
  BaseEvent::dump(ind2);
}


void MolOrderShuffleEvent::step() {
  for (Partition& p: world->get_partitions()) {
    p.shuffle_schedulable_molecule_ids();
  }
}

} /* namespace mcell */
