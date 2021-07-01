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

#include <iostream>
#include <algorithm>

#include "sort_mols_by_subpart_event.h"
#include "world.h"
#include "partition.h"

using namespace std;

namespace MCell {

void SortMolsBySubpartEvent::dump(const string ind) const {
  cout << ind << "Sort mols by subpart event:\n";
  string ind2 = ind + "  ";
  BaseEvent::dump(ind2);
}


struct SubpartComparatorForId
{
  SubpartComparatorForId(const vector<uint>& mempart_indices_)
    : mempart_indices(mempart_indices_) {
  }

  bool operator () (const molecule_id_t id1, const molecule_id_t id2) const {
    assert(id1 < mempart_indices.size());
    assert(id2 < mempart_indices.size());
    return mempart_indices[id1] < mempart_indices[id2];
  }

  const vector<uint>& mempart_indices;
};


struct SubpartComparatorForMol
{
  SubpartComparatorForMol(const vector<uint>& mempart_indices_)
    : mempart_indices(mempart_indices_) {
  }

  bool operator () (const Molecule& m1, const Molecule& m2) const {
    assert(m1.id < mempart_indices.size());
    assert(m2.id < mempart_indices.size());
    return mempart_indices[m1.id] < mempart_indices[m2.id];
  }

  const vector<uint>& mempart_indices;
};


void SortMolsBySubpartEvent::step() {

  for (Partition& p: world->get_partitions()) {
    vector<Molecule>& molecules = p.get_molecules();

    if (molecules.empty()) {
      // simply skip if there are no molecules
      continue;
    }

    // may create too large array when multiple partitions will be used
    // but currently the molecule ID assignment is done in partitions anyway
    vector<uint> mempart_indices(p.get_next_molecule_id_no_increment());

    for (Molecule& m: molecules) {
      const BNG::Species& sp = p.get_species(m.species_id);
      if (m.is_vol() && sp.can_diffuse()) {
        mempart_indices[m.id] = m.v.subpart_index;
      }
      else {
        // set some value that should not collide with subpart indices
        mempart_indices[m.id] = 0x80000000;
      }
    }

    // also sort molecules in the molecules array
    sort(molecules.begin(), molecules.end(), SubpartComparatorForMol(mempart_indices));

    // and update their indices in the molecule_id_to_index_mapping
    std::vector<molecule_index_t>& molecule_id_to_index_mapping = p.get_molecule_id_to_index_mapping();
    for (size_t i = 0; i < molecules.size(); i++) {
      const Molecule& m = molecules[i];
      assert(molecule_id_to_index_mapping.size() < m.id);
      molecule_id_to_index_mapping[m.id] = i;
    }
  }
}

} /* namespace mcell */
