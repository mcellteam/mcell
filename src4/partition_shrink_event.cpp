/******************************************************************************
 *
 * Copyright (C) 2019 by
 * The Salk Institute for Biological Studies and
 * Pittsburgh Supercomputing Center, Carnegie Mellon University
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

#include <iostream>
#include <algorithm>

#include "partition_shrink_event.h"
#include "world.h"
#include "partition.h"

using namespace std;

namespace MCell {

void PartitionShrinkEvent::dump(const string ind) const {
  cout << ind << "Partition shrink event:\n";
  string ind2 = ind + "  ";
  BaseEvent::dump(ind2);
}


void PartitionShrinkEvent::step() {
  for (Partition& p: world->get_partitions()) {
    // we might need more things to do here
    p.shrink_all_volume_molecule_reactants_per_subpart();
  }
}

} /* namespace mcell */
