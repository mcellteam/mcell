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
  SubpartComparatorForId(const Partition& p_)
    : p(p_) {
  }

  bool operator () (const molecule_id_t id1, const molecule_id_t id2) const {
    const Molecule& m1 = p.get_m(id1);
    const Molecule& m2 = p.get_m(id2);

    return m1.mempart_index < m2.mempart_index;
  }

  const Partition& p;
};


struct SubpartComparatorForMol
{
  SubpartComparatorForMol(){
  }

  bool operator () (const Molecule& m1, const Molecule& m2) const {
    return m1.mempart_index < m2.mempart_index;
  }
};


void SortMolsBySubpartEvent::step() {
  for (Partition& p: world->get_partitions()) {
    vector<Molecule>& molecules = p.get_molecules();

    if (molecules.empty()) {
      // simply skip if there are no molecules
      continue;
    }

    // update mempart_index
    for (Molecule& m: molecules) {
      const BNG::Species& sp = p.get_all_species().get(m.species_id);
      if (m.is_vol() && sp.can_diffuse()) {
        m.mempart_index = m.v.subpart_index;
      }
      else {
        // set some value that should not collide with subpart indices
        m.mempart_index = 0x80000000;
      }
    }

    // this can be optimized if we would take just volume mols into account
    vector<Partition::TimeStepMoleculesData>& mols_per_time_step = p.get_molecule_data_per_time_step_array();

    for (Partition::TimeStepMoleculesData& time_step_data: mols_per_time_step) {
      vector<molecule_id_t>& molecule_ids = time_step_data.molecule_ids;

      sort(molecule_ids.begin(), molecule_ids.end(), SubpartComparatorForId(p));
    }

    // also sort molecules in the molecules array
    sort(molecules.begin(), molecules.end(), SubpartComparatorForMol());

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
