
#include "reactions_info.h"

#include "defines.h"
#include "species.h"
#include "molecule.h"
#include "species.h"

#include <iostream>
#include <sstream>

#include "reaction.h"

using namespace std;

namespace MCell {

void ReactionsInfo::init(const BNG::SpeciesContainer<Species>& all_species) {
  assert(!initialized); // there should not be a reason to reinitialize reactions info

  // create map for fast reaction searches
  for (RxnClass& r: reactions) {
    assert(r.reactants.size() == 1 || r.reactants.size() == 2); // only bimolecular reactions are supported now

    if (r.reactants.size() == 1) {
      // for now we only support only one outcome of a bimolecular reaction
      assert(unimolecular_reactions_map.count(r.reactants[0].species_id) == 0);
      unimolecular_reactions_map[r.reactants[0].species_id] = &r;
    }
    else {
      // check - for now we only support only one outcome of a bimolecular reaction
      if (bimolecular_reactions_map.count(r.reactants[0].species_id) != 0) {
        assert(bimolecular_reactions_map[r.reactants[0].species_id].count(r.reactants[1].species_id) == 0);
      }

      // add mapping for both reactants pairs A + B and also B + A
      // NOTE: if we would sort the species id, we might have just one-directional mapping
      bimolecular_reactions_map[r.reactants[0].species_id][r.reactants[1].species_id] = &r;
      bimolecular_reactions_map[r.reactants[1].species_id][r.reactants[0].species_id] = &r;
    }
  }

  // just to make sure that we have an item for all the species
  const std::vector<Species>& species = all_species.get_species_vector();
  for (const Species& s: species) {
    bimolecular_reactions_map.insert( std::make_pair(s.species_id, SpeciesRxnClassesMap()) );
  }
  assert(bimolecular_reactions_map.size() == all_species.get_count());

  initialized = true;
}

}
