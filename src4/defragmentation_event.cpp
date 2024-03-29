/******************************************************************************
 *
 * Copyright (C) 2019 by
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

    // first cleanup schedulable_molecule_ids by removing all mols that are defunct
    // must be done before the following cleanup
    auto& schedulable_mol_ids = p.get_schedulable_molecule_ids();
    auto it_new_end = remove_if(schedulable_mol_ids.begin(), schedulable_mol_ids.end(),
        [&p](const molecule_id_t id) -> bool { return p.get_m(id).is_defunct(); });
    schedulable_mol_ids.erase(it_new_end, schedulable_mol_ids.end());

    // remove defunct molecules in the molecules array
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
