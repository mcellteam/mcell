
#include "bng/rxn_container.h"
#include "bng/rxn_class.h"

#include <iostream>
#include <sstream>


using namespace std;

namespace BNG {

RxnContainer::~RxnContainer() {
  for (RxnClass* rc: rxn_classes) {
    delete rc;
  }

  for (RxnRule* rxn: rxn_rules) {
    delete rxn;
  }
}


void RxnContainer::reset_caches() {
  for (RxnClass* rc: rxn_classes) {
    delete rc;
  }
  rxn_classes.clear();

  species_processed_for_bimol_rxn_classes.clear();
  species_processed_for_unimol_rxn_classes.clear();
  unimol_rxn_class_map.clear();
  bimol_rxn_class_map.clear();

  // we must also erase all references to rxns classes from rxn rules
  for (RxnRule* rxn: rxn_rules) {
    rxn->reset_rxn_classes_where_used();
  }
}


void RxnContainer::update_all_mols_and_mol_type_compartments() {

  species_id_t all_mols_id = all_species.get_all_molecules_species_id();
  species_id_t all_vol_mols_id = all_species.get_all_volume_molecules_species_id();
  species_id_t all_surf_mols_id = all_species.get_all_surface_molecules_species_id();

  vector<species_id_t> superspecies {
    all_mols_id, all_vol_mols_id, all_surf_mols_id
  };

  // set species flags (we already have all the reactions)
  for (species_id_t species_id: superspecies) {
    BNG::Species& sp = all_species.get(species_id);

    // get reactions, this also creates all reaction classes for the species that we are processing
    // we are getting reactions for all known species because these flags must be initialized
    // before we have any instance of these species
    BNG::ReactantRxnClassesMap* rxn_classes = get_bimol_rxns_for_reactant(Reactant(sp.id, COMPARTMENT_ID_ANY), true);
    if (rxn_classes == nullptr) {
      continue;
    }

    // go through all applicable reactants
    for (auto it: *rxn_classes) {
      const BNG::RxnClass* rxn_class = it.second[COMPARTMENT_ID_ANY];
      assert(rxn_class->is_bimol());

      const BNG::Species& sp2 = all_species.get(rxn_class->get_second_species_id(sp.id));

      // we can use is_vol/is_surf for ALL_VOLUME_MOLECULES and ALL_SURFACE_MOLECULES
      // but not for all_volume molecules because there is no single flag that we can query
      if (sp.is_vol() || sp.id == all_mols_id) {
        if (sp2.is_reactive_surface()) {
          sp.set_flag(BNG::SPECIES_FLAG_CAN_VOLWALL);
          if (sp.id == all_mols_id || sp.id == all_vol_mols_id) {
            all_vol_mols_can_react_with_surface = true;
          }
        }
      }

      if ((sp.is_surf() || sp.id == all_mols_id)) {
        if (sp2.is_reactive_surface()) {
          if (sp.id == all_surf_mols_id || sp.id == all_mols_id) {
            all_surf_mols_can_react_with_surface = true;
          }
        }
      }
    }
  }

  // set flag to all MolTypes that use compartments
  // this does not update already existing complexes,
  for (RxnRule* rxn: rxn_rules) {
    if (rxn->reactants_use_compartments()) {
      for (Cplx& reac: rxn->reactants) {
        if (reac.has_compartment()) {
          // update molecule types
          for (const MolInstance& mi: reac.mol_instances) {
            MolType& mt = bng_data.get_molecule_type(mi.mol_type_id);

            mt.set_flag(SPECIES_CPLX_MOL_FLAG_COMPARTMENT_USED_IN_RXNS);
            // collect compartments
            mt.reactant_compartments.insert(reac.get_compartment_id());
          }
          // and also complexes used as reactants
          reac.set_flag(SPECIES_CPLX_MOL_FLAG_COMPARTMENT_USED_IN_RXNS);
        }
      }
    }
  }

  // update this flag in cplxs of rxn products as well
  // this is run during initialization and no other complexes shoudl exist
  for (RxnRule* rxn: rxn_rules) {
    for (Cplx& reac: rxn->products) {
      reac.update_flag_and_compartments_used_in_rxns();
    }
  }

}


RxnClass* RxnContainer::get_or_create_empty_unimol_rxn_class(const Reactant& reac) {

  ReactantSpeciesIt it_species = unimol_rxn_class_map.find(reac.species_id);

  if (it_species != unimol_rxn_class_map.end()) {
    // map for this species id was already created,
    // was this compartment handled as well?
    ReactantCompartmentIt it_comp = it_species->second.find(reac.compartment_id);
    if (it_comp != it_species->second.end()) {
      return it_comp->second;
    }
  }

  // no created yet
  RxnClass* new_rxn_class = new RxnClass(*this, all_species, bng_config, reac);
  rxn_classes.insert(new_rxn_class);
  unimol_rxn_class_map[reac.species_id][reac.compartment_id] = new_rxn_class;
  return new_rxn_class;
}


void RxnContainer::create_unimol_rxn_classes_for_new_species(const species_id_t species_id) {

  // find all reactions for species id, we don't care about compartment yet
  small_vector<RxnRule*> rxns_for_new_species;
  for (RxnRule* r: rxn_rules) {
    if (r->is_unimol() && r->species_can_be_reactant(species_id, all_species)) {
      rxns_for_new_species.push_back(r);
    }
  }

  if (!rxns_for_new_species.empty()) {

    // for each compartment
    const CompartmentIdSet& reactant_compartments =
        all_species.get(species_id).get_reactant_compartments();
    for (compartment_id_t compartment_id1: reactant_compartments) {
      Reactant reac = Reactant(species_id, compartment_id1);

      // 1) first we need to get to the instance of the reaction class for new_id
      RxnClass* rxn_class = get_or_create_empty_unimol_rxn_class(reac);

      // create reactions classes specific for our species
      for (RxnRule* matching_rxn: rxns_for_new_species) {
        // 2) check if compartment matches
        if (!matching_rxn->reactant_compatment_matches(0, reac.compartment_id)) {
          continue;
        }

        // 3) add the matching_rxn to our rxn class
        rxn_class->add_rxn_rule_no_update(matching_rxn);
      }
      // TODO: initialize rxn class pathways for debugging so that
      // it is clear what can be the products

      if (bng_config.bng_verbosity_level >= 2) {
        cout << "BNG: Created a new unimolecular reaction class:\n";
        rxn_class->dump();
      }
    }
  }
}


RxnClass* RxnContainer::get_or_create_empty_bimol_rxn_class(const Reactant& reac1, const Reactant& reac2) {
  // TODO: simplify, can we use []?
  BimolSpeciesIt it_species1 = bimol_rxn_class_map.find(reac1.species_id);
  if (it_species1 == bimol_rxn_class_map.end()) {
    auto it_species_pair = bimol_rxn_class_map.insert( make_pair(reac1.species_id, BimolCompartmentRxnClassesMap()) );
    it_species1 = it_species_pair.first;
  }
  BimolCompartmentIt it_species1_comp1 = it_species1->second.find(reac1.compartment_id);
  if (it_species1_comp1 == it_species1->second.end()) {
    auto it_comp_pair = it_species1->second.insert( make_pair(reac1.compartment_id, ReactantRxnClassesMap()) );
    it_species1_comp1 = it_comp_pair.first;
  }

  BimolSpeciesIt it_species2 = bimol_rxn_class_map.find(reac2.species_id);
  if (it_species2 == bimol_rxn_class_map.end()) {
    auto it_species_pair = bimol_rxn_class_map.insert( make_pair(reac2.species_id, BimolCompartmentRxnClassesMap()) );
    it_species2 = it_species_pair.first;
  }
  BimolCompartmentIt it_species2_comp2 = it_species2->second.find(reac2.compartment_id);
  if (it_species2_comp2 == it_species2->second.end()) {
    auto it_comp_pair = it_species2->second.insert( make_pair(reac2.compartment_id, ReactantRxnClassesMap()) );
    it_species2_comp2 = it_comp_pair.first;
  }

  ReactantSpeciesIt it_species1_comp1_species2 = it_species1_comp1->second.find(reac2.species_id);

  // do we already have a reaction class for (reac1, reac2)?
  if (it_species1_comp1_species2 != it_species1_comp1->second.end() &&
      it_species1_comp1_species2->second.find(reac2.compartment_id) != it_species1_comp1_species2->second.end()) {
#ifndef NDEBUG
    ReactantSpeciesIt it_species2_comp2_species1 = it_species2_comp2->second.find(reac1.species_id);
    assert(it_species2_comp2_species1 != it_species2_comp2->second.end() && "Map for id1->id2 implies that there must be a map for id2->id1");
    assert(it_species1_comp1_species2->second.count(reac2.compartment_id));
    assert(it_species2_comp2_species1->second.count(reac1.compartment_id));
    assert(it_species1_comp1_species2->second[reac2.compartment_id] == it_species2_comp2_species1->second[reac1.compartment_id] &&
        "Map for id1->id2 must be the same as a map for id2->id1");
#endif
    return it_species1_comp1_species2->second[reac2.compartment_id];
  }
  else {
    // create a new one
    RxnClass* new_rxn_class = new RxnClass(*this, all_species, bng_config, reac1, reac2);
    rxn_classes.insert(new_rxn_class);

    // insert it into maps
    it_species1_comp1->second[reac2.species_id][reac2.compartment_id] = new_rxn_class;
    it_species2_comp2->second[reac1.species_id][reac1.compartment_id] = new_rxn_class;
    return new_rxn_class;
  }
}


// - puts pointers to all corresponding classes to the res_classes_map
// - for bimol rxns, does not reuse already defined rxn class, e.g. when A + B was already created,
//   rxn class for B + A will be created (NOTE: might improve if needed but so far the only issue
//   are reports and printouts
// - called only from get_bimol_rxns_for_reactant
void RxnContainer::create_bimol_rxn_classes_for_new_species(const species_id_t species_id1, const bool for_all_known_species) {

  // find all reactions for species id
  small_vector<RxnRule*> rxns_for_new_species;
  for (RxnRule* r: rxn_rules) {
    if (r->is_bimol() && r->species_can_be_reactant(species_id1, all_species)) {
      rxns_for_new_species.push_back(r);
    }
  }


  // do nothing if this species cannot react
  if (!rxns_for_new_species.empty()) {

    const CompartmentIdSet& reactant_compartments1 =
        all_species.get(species_id1).get_reactant_compartments();
    for (compartment_id_t compartment_id1: reactant_compartments1) {
      Reactant reac1 = Reactant(species_id1, compartment_id1);

      // create reactions classes specific for our species
      // rxn_class->update_rxn_pathways may create new species, therefore we
      // must always read the current all_species contents,
      // on the other hand, polymerizing reactions might cause infinite looping therefore
      // we will limit ourselves to the species that currently exist,
      // the rxn class will be updated once a molecule of the new species (not handled here) will be created
      species_index_t num_species = all_species.get_species_vector().size();
      for (species_index_t i = 0; i < num_species; i++) {
        // reading directly from the species array, not using id
        const Species& species = all_species.get_species_vector()[i];
        species_id_t species_id2 = species.id;

        // TODO: simplify condition - is_species_superclass check does not have to be there
        if (!for_all_known_species &&
            !species.was_instantiated()  &&
            !(reac1.species_id == species_id2) && // we would miss A+A type reactions
            !all_species.is_species_superclass(species_id2)) {
          // we do not care about molecules that do not exist yet, however we must process superclasses
          continue;
        }

        // TODO: can we share RxnClass object for cases when ANY and NONE are the same?
        //       what about cleanup?
        // for each applicable compartment for this species
        const CompartmentIdSet& reactant_compartments2 =
            all_species.get(species_id2).get_reactant_compartments();
        for (compartment_id_t compartment_id2: reactant_compartments2) {

          Reactant reac2(species_id2, compartment_id2);

          // don't we have a rxn class for this pair of special already?
          // (only created in different order, e.g. in A + B and now we are checking for B + A)
          // TODO: move to a function
          // a) second species
          BimolSpeciesIt it_species2 = bimol_rxn_class_map.find(reac2.species_id);
          if (it_species2 != bimol_rxn_class_map.end()) {
            // b) second compartment
            BimolCompartmentIt it_species2_comp2 = it_species2->second.find(reac2.compartment_id);
            if (it_species2_comp2 != it_species2->second.end()) {
              // c) first species
              ReactantSpeciesIt it_species2_comp2_species1 = it_species2_comp2->second.find(reac1.species_id);
              if (it_species2_comp2_species1 != it_species2_comp2->second.end()) {
                // d) first compartment
                ReactantCompartmentIt it_species2_comp2_species1_comp1 =
                    it_species2_comp2_species1->second.find(reac1.compartment_id);
                if (it_species2_comp2_species1_comp1 != it_species2_comp2_species1->second.end()) {
                  // nothing to do here, this rxn class is already correctly mapped in rxn class map
                  continue;
                }
              }
            }
          }

          small_vector<RxnRule*> applicable_rxns;

          for (RxnRule* matching_rxn: rxns_for_new_species) {

            // usually the species must be different but reactions of type A + A are allowed
            if (reac1.species_id == species_id2 && !matching_rxn->species_is_both_bimol_reactants(reac1.species_id, all_species)) {
              continue;
            }

            uint reac1_index, reac2_index;
            bool reactants_match =
                matching_rxn->species_can_be_bimol_reactants(reac1.species_id, species_id2, all_species, &reac1_index, &reac2_index);

            if (reactants_match &&
                matching_rxn->reactant_compatment_matches(reac1_index, reac1.compartment_id) &&
                matching_rxn->reactant_compatment_matches(reac2_index, reac2.compartment_id)) {

              // ok, we have a reaction applicable both to new_id and second_id and compartment matches as well
              // we need to add this rxn to a rxn class for these reactants
              applicable_rxns.push_back(matching_rxn);
            }
          }

          if (!applicable_rxns.empty()) {
            // get to the instance of the reaction class for (new_id, second_id) or (second_id, new_id)
            RxnClass* rxn_class = get_or_create_empty_bimol_rxn_class(reac1, reac2);

            for (RxnRule* rxn: applicable_rxns) {
              rxn_class->add_rxn_rule_no_update(rxn);
            }
            // TODO: init rxn class pathways for debug

            if (bng_config.bng_verbosity_level >= 2) {
              cout <<
                  "BNG: Created or updated a new bimolecular reaction class for species " <<
                  all_species.get(reac1.species_id).name << " (" << reac1.species_id << "):\n";
              rxn_class->dump("  ");
            }
          }
        }
      } // for rxns_for_new_species
    }
  }
}


void RxnContainer::delete_rxn_class(RxnClass* rxn_class) {
  // remove the pointer from the rxn classes container so that they are not deleted
  // when RxnContainer destructor is called
  cout << "\n\n" << std::hex << (void*)rxn_class << "\n";
  rxn_class->dump();

  assert(rxn_classes.count(rxn_class) != 0);
  rxn_classes.erase(rxn_class);

  delete rxn_class;
}


void RxnContainer::remove_unimol_rxn_classes(const species_id_t id) {
  if (species_processed_for_unimol_rxn_classes.count(id) != 0) {
    // forget that we processed this species
    species_processed_for_unimol_rxn_classes.erase(id);

    // and remove the rxn class if there is any
    ReactantSpeciesIt it_species = unimol_rxn_class_map.find(id);
    if (it_species != unimol_rxn_class_map.end()) {

      // remove the rxn classes for all compartments if there are any
      for (ReactantCompartmentPair& pair_species_comp: it_species->second) {
        delete_rxn_class(pair_species_comp.second);
      }

      unimol_rxn_class_map.erase(it_species);
    }
  }
  else {
    assert(unimol_rxn_class_map.count(id) == 0 &&
        "There must be no rxn classes for unprocessed species");
  }
}


void RxnContainer::remove_bimol_rxn_classes(const species_id_t reac1_species_id) {

  // remove the rxn classes and their mappings for the second reactants
  if (species_processed_for_bimol_rxn_classes.count(reac1_species_id) != 0) {
    // forget that we processed this species
    species_processed_for_bimol_rxn_classes.erase(reac1_species_id);

    BimolSpeciesIt it_species1 = bimol_rxn_class_map.find(reac1_species_id);
    if (it_species1 != bimol_rxn_class_map.end()) {
      std::set<species_id_t> reacting_species;

      // for each compartment
      for (BimolCompartmentPair& pair_species1_comp1: it_species1->second) {

        // delete all rxn classes for this species and remember
        // with which species it can react
        for (ReactantSpeciesPair& pair_species1_comp1_species2: pair_species1_comp1.second) {
          reacting_species.insert(pair_species1_comp1_species2.first);

          for (ReactantCompartmentPair& pair_species1_comp1_species2_comp: pair_species1_comp1_species2.second) {
            // delete rxn class only when we did not delete it before, this can happen
            // when there is a rxn in form A@{C,D} + A@{C,D}, deleting A@C + A@D also
            // deletes A@D + A@C
            if (rxn_classes.count(pair_species1_comp1_species2_comp.second) != 0) {
              delete_rxn_class(pair_species1_comp1_species2_comp.second);
            }
          }
        }
      }

      // ok, we removed all A + X mappings and also the rxn classes,
      // we must also remove mappings X + A
      for (species_id_t reac2_species_id: reacting_species) {
        if (reac2_species_id == reac1_species_id) {
          continue;
        }
        BimolSpeciesIt it_species2 = bimol_rxn_class_map.find(reac2_species_id);
        assert(it_species2 != bimol_rxn_class_map.end() &&
            "Reactant must be present in the rxn classes map");

        // for each compartment
        for (BimolCompartmentPair& pair_species2_comp2: it_species2->second) {
          ReactantSpeciesIt it_species2_comp2_species1 = pair_species2_comp2.second.find(reac1_species_id);
          assert(it_species2_comp2_species1 != pair_species2_comp2.second.end() && "Reverse mapping must exist");
          pair_species2_comp2.second.erase(it_species2_comp2_species1);
        }
      }

      // erase the top level entry
      bimol_rxn_class_map.erase(it_species1);
    }
  }
  else {
    assert(bimol_rxn_class_map.count(reac1_species_id) == 0 &&
        "There must be no rxn class maps for unprocessed species");
  }
}


bool RxnContainer::has_bimol_vol_rxns() const {
  for (const BNG::RxnRule* r: rxn_rules) {
    if (r->is_bimol_vol_rxn()) {
      return true;
    }
  }
  return false;
}


void RxnContainer::dump(const bool including_rxn_rules) const {

  for (auto pair_species: unimol_rxn_class_map) {
    for (auto pair_species_comp: pair_species.second) {
      cout <<
          "RxnClass for " <<
          all_species.get(pair_species.first).name << " (" << pair_species.first << "@" << pair_species_comp.first << "):\n";

      const RxnClass* rxn_class = pair_species_comp.second;
      assert(rxn_class != nullptr);
      rxn_class->dump("  ");
      cout << "\n";
    }
  }

  for (auto pair_species1: bimol_rxn_class_map) {
    for (auto pair_species1_comp1: pair_species1.second) {
      for (auto pair_species1_comp1_species2: pair_species1_comp1.second) {
        for (auto pair_species1_comp1_species2_comp2: pair_species1_comp1_species2.second) {
          cout <<
              "RxnClass for " <<
              all_species.get(pair_species1.first).name <<
              " (" << pair_species1.first << "@" << pair_species1_comp1.first << ") + " <<
              all_species.get(pair_species1_comp1_species2.first).name <<
              " (" << pair_species1_comp1_species2.first << "@" << pair_species1_comp1_species2_comp2.first << "):\n";

          const RxnClass* rxn_class = pair_species1_comp1_species2_comp2.second;
          assert(rxn_class != nullptr);
          rxn_class->dump("  ");
          cout << "\n";
        }
      }
    }
  }

  if (including_rxn_rules) {
    for (uint i = 0; i < rxn_rules.size(); i++) {
      const RxnRule* r = rxn_rules[i];
      cout << "RxnRule " << i << ": \n";
      r->dump(true, "  ");
      cout << "\n";
    }
  }

  cout.flush();
}

} // namespace BNG
