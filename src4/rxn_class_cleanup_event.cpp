/******************************************************************************
 *
 * Copyright (C) 2020 by
 * The Salk Institute for Biological Studies
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
  /*
  // for all species
  for (BNG::Species& sp: world->get_all_species().get_species_vector()) {
    if (sp.get_num_instantiations() == 0) {
      // there are no instances, we can efficiently remove all rxn classes for this species

      // remove rxn classes
      world->get_all_rxns().remove_unimol_rxn_class(sp.id);
      world->get_all_rxns().remove_bimol_rxn_classes(sp.id);

      // and clear flag that this species was instantiated
      sp.set_was_instantiated(false);
    }
  }
  */
}


} /* namespace MCell */
