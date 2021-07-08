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

#include "rxn_class_cleanup_event.h"

#include "world.h"

using namespace std;

namespace MCell {


void RxnClassCleanupEvent::dump(const string ind) const {
  cout << ind << "Rxn class cleanup event:\n";
  string ind2 = ind + "  ";
  BaseEvent::dump(ind2);
}


void RxnClassCleanupEvent::step() {

  // for all species
  for (BNG::Species* sp: world->get_all_species().get_species_vector()) {
    release_assert(sp != nullptr);
    if (sp->get_num_instantiations() == 0 && sp->was_instantiated()) {
      // there are no instances, we can efficiently remove all rxn classes for this species

      // remove rxn classes
      world->get_all_rxns().remove_unimol_rxn_classes(sp->id);
      world->get_all_rxns().remove_bimol_rxn_classes(sp->id);

      // clear flag telling that this species was instantiated
      sp->set_was_instantiated(false);

      // and also tell partitions that this species is not known anymore
      // TODO: can this may be done through the was_instantiated flag?
      if (sp->is_vol()) {
        for (Partition& p: world->get_partitions()) {
          p.remove_from_known_vol_species(sp->id);
        }
      }
    }
  }
}


} /* namespace MCell */
