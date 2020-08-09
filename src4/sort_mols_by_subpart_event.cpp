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


struct SubpartComparator
{
  SubpartComparator(const Partition& p_)
    : p(p_) {
  }

  // Compare 2 Player objects using name
  bool operator () (const molecule_id_t id1, const molecule_id_t id2) const
  {
    const Molecule& m1 = p.get_m(id1);
    const Molecule& m2 = p.get_m(id2);
    const BNG::Species& s1 = p.get_all_species().get(m1.species_id);
    const BNG::Species& s2 = p.get_all_species().get(m2.species_id);

    if (m1.is_vol() && m2.is_vol() && s1.can_diffuse() && s2.can_diffuse()) {
      if (m1.v.subpart_index != m2.v.subpart_index) {
        return m1.v.subpart_index < m2.v.subpart_index;
      }
      else {
        return id1 < id2;
      }
    }
    else {
      return id1 < id2;
    }
  }

  const Partition& p;
};


void SortMolsBySubpartEvent::step() {
  for (Partition& p: world->get_partitions()) {
    vector<Molecule>& volume_molecules = p.get_molecules();

    if (volume_molecules.empty()) {
      // simply skip if there are no volume molecules
      continue;
    }

    // TODO: this can be optimized if we would take just
    // volume mols into account
    vector<Partition::TimeStepMoleculesData>& mols_per_time_step = p.get_molecule_data_per_time_step_array();

    for (Partition::TimeStepMoleculesData& time_step_data: mols_per_time_step) {
      std::vector<molecule_id_t>& molecule_ids = time_step_data.molecule_ids;

      std::sort(molecule_ids.begin(), molecule_ids.end(), SubpartComparator(p));
    }
  }
}

} /* namespace mcell */
