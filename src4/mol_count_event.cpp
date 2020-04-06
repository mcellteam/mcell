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

#include "mol_count_event.h"
#include "world.h"
#include "partition.h"

using namespace std;

namespace MCell {

void MolCountInfo::dump(const std::string ind) {

  cout << ind << "type: ";
  switch(type) {
    case CountType::Invalid:
      cout << "Invalid";
      break;
    case CountType::World:
      cout << "World";
      break;
    case CountType::EnclosedInObject:
      cout << "EnclosedInObject";
      break;
    default:
      assert(false);
  }
  cout << "\n";

  cout << ind << "buffer_id: \t\t" << buffer_id << " [count_buffer_id_t] \t\t\n";
  cout << ind << "orientation: \t\t" << orientation << " [orientation_t] \t\t\n";
  cout << ind << "species_id: \t\t" << species_id << " [species_id_t] \t\t\n";
  cout << ind << "geometry_object_id: \t\t" << geometry_object_id << " [geometry_object_id_t] \t\t\n";
}


void MolCountEvent::step() {

  // go through all molecules and count them
  PartitionVector& partitions = world->get_partitions();
  std::vector<CountItem> count_items;
  count_items.resize(mol_count_infos.size());

  // initialize new count items
  for (uint i = 0; i < mol_count_infos.size(); i++) {
    count_items[i].time = event_time * world->config.time_unit;
    count_items[i].int_value = 0;
  }

  // for each partition
  for (const Partition& p: partitions) {

    // for each molecule
    for (const Molecule& m: p.get_molecules()) {

      const BNG::Species& species = world->get_all_species().get(m.species_id);
      if (!species.has_count_enclosed_flag()) {
        continue;
      }

      const uint_set<geometry_object_id_t>& enclosing_volumes =
          p.get_enclosing_counted_volumes(m.v.counted_volume_id);

      // for each counting info
      for (uint i = 0; i < mol_count_infos.size(); i++) {
        const MolCountInfo& info = mol_count_infos[i];

        if (info.type == CountType::World) {
          // count the molecule
          count_items[i].inc();
        }
        else if (info.type == CountType::EnclosedInObject) {
          // is the molecule inside of the object that we are checking?
          if (m.v.counted_volume_id == info.geometry_object_id ||
              enclosing_volumes.count(info.geometry_object_id)) {

            count_items[i].inc();
          }
        }
      }

    }
  }

  // check each molecule against what we are checking
  for (uint i = 0; i < mol_count_infos.size(); i++) {
    world->get_count_buffer(mol_count_infos[i].buffer_id).add(count_items[i]);
  }

}


void MolCountEvent::dump(const std::string ind) {
  cout << ind << "Mol count event:\n";
  std::string ind2 = ind + "  ";
  std::string ind4 = ind2 + "  ";
  BaseEvent::dump(ind2);

  cout << ind << " mol_count_infos:\n";
  for(uint i = 0; i < mol_count_infos.size(); i++) {
    cout << ind2 << i << "\n";
    mol_count_infos[i].dump(ind4);
  }

}


}
