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
    vector<Molecule>& molecules = p.get_molecules();

    if (molecules.empty()) {
      // simply skip if there are no volume molecules
      continue;
    }

    // first remove all defunct molecules from mols_per_time_step
    vector<Partition::TimeStepMoleculesData>& mols_per_time_step = p.get_molecule_data_per_time_step_array();
    for (Partition::TimeStepMoleculesData& timestep_mol_ids: mols_per_time_step) {
      size_t ids_removed = 0;

      vector<molecule_id_t>& molecule_ids = timestep_mol_ids.molecule_ids;

      typedef vector<molecule_id_t>::iterator id_it_t;
      id_it_t id_it_end = molecule_ids.end();

      // find first defunct molecule
      id_it_t id_it_first_defunct =  find_if(molecule_ids.begin(), id_it_end, [&p](const molecule_id_t id) -> bool { return p.get_m(id).is_defunct(); });
      id_it_t id_it_copy_destination = id_it_first_defunct;

      while (id_it_first_defunct != id_it_end) {
        // then find the next one that is not defunct (might be it_end)
        id_it_t id_it_next_funct = find_if(id_it_first_defunct, id_it_end, [&p](const molecule_id_t id) -> bool { return !p.get_m(id).is_defunct(); });

        // then again, find following defunct molecule
        id_it_t id_it_second_defunct = find_if(id_it_next_funct, id_it_end, [&p](const molecule_id_t id) -> bool { return p.get_m(id).is_defunct(); });

        // move data: from, to, into position
        std::copy(id_it_next_funct, id_it_second_defunct, id_it_copy_destination);

        ids_removed += id_it_next_funct - id_it_first_defunct;

        // and also move destination pointer
        id_it_copy_destination += id_it_second_defunct - id_it_next_funct;

        id_it_first_defunct = id_it_second_defunct;
      }

      // remove everything after it_copy_destination
      if (ids_removed != 0) {
        molecule_ids.resize(molecule_ids.size() - ids_removed);
      }
    }

    // then follow with the molecules array
    vector<molecule_index_t>& molecule_id_to_index_mapping = p.get_molecule_id_to_index_mapping();

#ifdef DEBUG_DEFRAGMENTATION
    cout << "Defragmentation before sort:\n";
    Molecule::dump_array(molecules);
#endif

    size_t removed = 0;

    typedef vector<Molecule>::iterator vmit_t;
    vmit_t it_end = molecules.end();

    // find first defunct molecule
    vmit_t it_first_defunct =  find_if(molecules.begin(), it_end, [](const Molecule & m) -> bool { return m.is_defunct(); });
    vmit_t it_copy_destination = it_first_defunct;

    while (it_first_defunct != it_end) {

      // then find the next one that is not defunct (might be it_end)
      vmit_t it_next_funct = find_if(it_first_defunct, it_end, [](const Molecule & m) -> bool { return !m.is_defunct(); });

      // then again, find following defunct molecule
      vmit_t it_second_defunct = find_if(it_next_funct, it_end, [](const Molecule & m) -> bool { return m.is_defunct(); });

      // items between it_first_defunct and it_next_funct will be removed
      // for debug mode, we will set their ids in volume_molecules_id_to_index_mapping as invalid
#ifndef NDEBUG
      for (vmit_t it_update_mapping = it_first_defunct; it_update_mapping != it_next_funct; it_update_mapping++) {
        const Molecule& vm = *it_update_mapping;
        assert(vm.is_defunct());
        molecule_id_to_index_mapping[vm.id] = MOLECULE_INDEX_INVALID;
      }
#endif

      // move data: from, to, into position
      std::copy(it_next_funct, it_second_defunct, it_copy_destination);

      removed += it_next_funct - it_first_defunct;

      // and also move destination pointer
      it_copy_destination += it_second_defunct - it_next_funct;

      it_first_defunct = it_second_defunct;
    }

    // remove everything after it_copy_destination
    if (removed != 0) {
      molecules.resize(molecules.size() - removed);
    }

#ifdef DEBUG_DEFRAGMENTATION
    cout << "Defragmentation after defunct removal:\n";
    Molecule::dump_array(molecules);
#endif

    // update mapping
    size_t new_count = molecules.size();
    for (size_t i = 0; i < new_count; i++) {
      const Molecule& vm = molecules[i];
      if (vm.is_defunct()) {
        break;
      }
      // correct index because the molecule could have been moved
      molecule_id_to_index_mapping[vm.id] = i;
    }
  }
}

} /* namespace mcell */
