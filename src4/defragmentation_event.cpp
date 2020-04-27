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

#include "defragmentation_event.h"
#include "world.h"
#include "partition.h"

using namespace std;

namespace MCell {

void DefragmentationEvent::dump(const string ind) const {
  cout << ind << "Defragmentation event:\n";
  string ind2 = ind + "  ";
  BaseEvent::dump(ind2);
}


void DefragmentationEvent::step() {
  for (Partition& p: world->get_partitions()) {
    vector<Molecule>& volume_molecules = p.get_molecules();
    vector<molecule_index_t>& volume_molecules_id_to_index_mapping = p.get_molecule_id_to_index_mapping();

    vector<Partition::TimeStepMoleculesData>& mols_per_time_step = p.get_molecule_data_per_time_step_array();
    assert(mols_per_time_step.size() == 1 && mols_per_time_step[0].molecule_ids.size() == volume_molecules.size()
        && "For now, volume_molecule_indices_per_time_step[0] must be identical to volume_molecules");
    vector<molecule_id_t>& volume_molecule_ids_per_time_step = mols_per_time_step[0].molecule_ids;

#ifdef DEBUG_DEFRAGMENTATION
    cout << "Defragmentation before sort:\n";
    Molecule::dump_array(volume_molecules);
#endif

    typedef vector<Molecule>::iterator vmit_t;
    vmit_t it_begin = volume_molecules.begin();
    vmit_t it_end = volume_molecules.end();

    // find first defunct molecule
    vmit_t it_first_defunct =  find_if(volume_molecules.begin(), it_end, [](const Molecule & m) -> bool { return m.is_defunct(); });
    vmit_t it_copy_destination = it_first_defunct;
    size_t removed = 0;

    vector<molecule_id_t>::iterator it_indices_begin = volume_molecule_ids_per_time_step.begin();

    while (it_first_defunct != it_end) {

      // then find the next one that is not defunct
      vmit_t it_next_funct = find_if(it_first_defunct, it_end, [](const Molecule & m) -> bool { return !m.is_defunct(); });

      if (it_next_funct == it_end) {
        // are we finished?
        break;
      }

      // then again, find following defunct molecule
      vmit_t it_second_defunct = find_if(it_next_funct, it_end, [](const Molecule & m) -> bool { return m.is_defunct(); });

      // items between it_first_defunct and it_next_funct will be removed
      // for debug mode, we will set their ids in volume_molecules_id_to_index_mapping as invalid
#ifndef NDEBUG
      for (vmit_t it_update_mapping = it_first_defunct; it_update_mapping != it_next_funct; it_update_mapping++) {
        const Molecule& vm = *it_update_mapping;
        assert(vm.is_defunct());
        volume_molecules_id_to_index_mapping[vm.id] = MOLECULE_INDEX_INVALID;
      }
#endif

      // move data: from, to, into position
      std::copy(it_next_funct, it_second_defunct, it_copy_destination);

      // do the same thing to the volume_molecules_per_time_step
      size_t next_funct_index = it_next_funct - it_begin;
      size_t second_defunct_index = it_second_defunct - it_begin;
      size_t copy_destination_index = it_copy_destination - it_begin;
      std::copy(it_indices_begin + next_funct_index, it_indices_begin + second_defunct_index, it_indices_begin + copy_destination_index);

      removed += it_next_funct - it_first_defunct;

      // and also move destination pointer
      it_copy_destination += it_second_defunct - it_next_funct;

       it_first_defunct = it_second_defunct;
    }

    // remove everything after it_copy_destination
    if (removed != 0) {
      volume_molecules.resize(volume_molecules.size() - removed);
      volume_molecule_ids_per_time_step.resize(volume_molecule_ids_per_time_step.size() - removed);
    }

#ifdef DEBUG_DEFRAGMENTATION
    cout << "Defragmentation after defunct removal:\n";
    Molecule::dump_array(volume_molecules);
#endif

    // update mapping
    size_t new_count = volume_molecules.size();
    for (size_t i = 0; i < new_count; i++) {
      const Molecule& vm = volume_molecules[i];
      if (vm.is_defunct()) {
        break;
      }
      // correct index because the molecule could have been moved
      volume_molecules_id_to_index_mapping[vm.id] = i;
    }
  }
}

} /* namespace mcell */
