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

// only one can be defined
#define SORT_BY_SUBPART_INDEX
//#define SORT_BY_MEM_CUBE

const int SUBPARTS_PER_MEM_CUBE = 4;

using namespace std;

namespace MCell {

void SortMolsBySubpartEvent::dump(const string ind) const {
  cout << ind << "Sort mols by subpart event:\n";
  string ind2 = ind + "  ";
  BaseEvent::dump(ind2);
}


uint subpart_index_to_mem_cube_index(const Partition& p, const subpart_index_t subpart_index) {
  IVec3 indices;
  p.get_subpart_3d_indices_from_index(subpart_index, indices);

  // "decimate" the subpart indices
  indices = IVec3(indices / IVec3(SUBPARTS_PER_MEM_CUBE));

  assert(indices.x < 256 && indices.y < 256 && indices.z < 256);
  // we just need to encode the memory subpart number in some way
  return indices.x + (indices.y << 8) + (indices.z << 16);
}


inline bool compare_mols(
    const Partition& p,
    const Molecule& m1,
    const Molecule& m2
) {
  const BNG::Species& s1 = p.get_all_species().get(m1.species_id);
  const BNG::Species& s2 = p.get_all_species().get(m2.species_id);

  if (m1.is_vol() && m2.is_vol() && s1.can_diffuse() && s2.can_diffuse()) {
#ifndef SORT_BY_MEM_CUBE
    if (m1.v.subpart_index != m2.v.subpart_index) {
      return m1.v.subpart_index < m2.v.subpart_index;
    }
    else {
      return m1.id < m2.id;
    }
#else
    uint mem_index1 = subpart_index_to_mem_cube_index(p, m1.v.subpart_index);
    uint mem_index2 = subpart_index_to_mem_cube_index(p, m2.v.subpart_index);
    if (mem_index1 != mem_index2) {
      return mem_index1 < mem_index2;
    }
# ifdef SORT_BY_SUBPART_INDEX
    else if (m1.v.subpart_index != m2.v.subpart_index) {
      return m1.v.subpart_index < m2.v.subpart_index;
    }
# endif
    else {
      return m1.id < m2.id;
    }
#endif
  }
  else {
    return m1.id < m2.id;
  }
}


struct SubpartComparatorForId
{
  SubpartComparatorForId(const Partition& p_)
    : p(p_) {
  }

  bool operator () (const molecule_id_t id1, const molecule_id_t id2) const
  {
    const Molecule& m1 = p.get_m(id1);
    const Molecule& m2 = p.get_m(id2);

    return compare_mols(p, m1, m2);
  }

  const Partition& p;
};


struct SubpartComparatorForMol
{
  SubpartComparatorForMol(const Partition& p_)
    : p(p_) {
  }

  bool operator () (const Molecule& m1, const Molecule& m2) const
  {
    return compare_mols(p, m1, m2);
  }

  const Partition& p;
};

void SortMolsBySubpartEvent::step() {
  for (Partition& p: world->get_partitions()) {
    vector<Molecule>& molecules = p.get_molecules();

    if (molecules.empty()) {
      // simply skip if there are no volume molecules
      continue;
    }

    // NOTE: this can be optimized if we would take just volume mols into account
    vector<Partition::TimeStepMoleculesData>& mols_per_time_step = p.get_molecule_data_per_time_step_array();

    for (Partition::TimeStepMoleculesData& time_step_data: mols_per_time_step) {
      vector<molecule_id_t>& molecule_ids = time_step_data.molecule_ids;

      sort(molecule_ids.begin(), molecule_ids.end(), SubpartComparatorForId(p));
    }

    // also sort molecules in the molecules array
    sort(molecules.begin(), molecules.end(), SubpartComparatorForMol(p));

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
