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
#include <fstream>
#include <sstream>

#include "species_cleanup_event.h"

#include "world.h"

using namespace std;

namespace MCell {


void SpeciesCleanupEvent::dump(const string ind) const {
  cout << ind << "Species cleanup event:\n";
  string ind2 = ind + "  ";
  BaseEvent::dump(ind2);
}


void SpeciesCleanupEvent::remove_unused_reactant_classes() {

  // get a set of used reactant classes
  uint_set<BNG::reactant_class_id_t> used_reactant_classes;
  for (BNG::Species* sp: world->get_all_species().get_species_vector()) {
    release_assert(sp != nullptr);

    if (sp->has_valid_reactant_class_id()) {
      used_reactant_classes.insert(sp->get_reactant_class_id());
    }
  }

  // go through all reactant classes and select those that are not used
  uint_set<BNG::reactant_class_id_t> unused_reactant_classes;
  for (const BNG::ReactantClass* rc: world->get_all_rxns().get_reactant_classes()) {
    if (rc != nullptr && used_reactant_classes.count(rc->id) == 0) {
      unused_reactant_classes.insert(rc->id);
    }
  }

  // remove the unused classes
  for (BNG::reactant_class_id_t id: unused_reactant_classes) {
    // from rxn container including reacting classes
    world->get_all_rxns().remove_reactant_class(id);

    // and from partition's reactants map
    world->get_partition(PARTITION_ID_INITIAL).remove_reactant_class_usage(id);
  }
}


void SpeciesCleanupEvent::step() {

  // remove all rxn classes and all caches
  world->get_all_rxns().reset_caches();

  for (BNG::Species* sp: world->get_all_species().get_species_vector()) {
    release_assert(sp != nullptr);

    // num_instantiations tells us that there are no molecules of these species
    // is_removable - species were created on the fly
    if (sp->get_num_instantiations() == 0 && sp->is_removable()) {

      // tell partitions that this species is not known anymore
      if (sp->is_vol()) {
        for (Partition& p: world->get_partitions()) {
          p.remove_from_known_vol_species(sp->id);
        }
      }

      if (world->config.rxn_and_species_report) {
        stringstream ss;
        ss << sp->id << ": " << sp->to_str() << " - removed\n";
        BNG::append_to_report(world->config.get_species_report_file_name(), ss.str());
      }

      // delete this species
      world->get_all_species().remove(sp->id);

      // and also from caches used by RxnRules
      world->get_all_rxns().remove_species_id_references(sp->id);
    }
  }

  // cleanup the species array, we must remove nullptrs from the species array
  world->get_all_species().defragment();

  // also remove unused reactant classes
  remove_unused_reactant_classes();

}


} /* namespace MCell */
