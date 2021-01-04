
#include "bng/rxn_container.h"
#include "bng/rxn_class.h"

#include <iostream>
#include <sstream>
#include <bitset>

using namespace std;

namespace BNG {

RxnContainer::~RxnContainer() {
  for (RxnClass* rc: rxn_classes) {
    delete rc;
  }

  for (RxnRule* rxn: rxn_rules) {
    delete rxn;
  }

  for (ReactantClass* rxn: reactant_classes_vector) {
    delete rxn;
  }
}


void RxnContainer::reset_caches() {
  for (RxnClass* rc: rxn_classes) {
    delete rc;
  }
  rxn_classes.clear();

  species_processed_for_bimol_rxn_classes.clear();
  species_processed_for_bimol_rxn_classes.shrink();
  species_processed_for_unimol_rxn_classes.clear();
  species_processed_for_bimol_rxn_classes.shrink();
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
          for (const ElemMol& mi: reac.elem_mols) {
            ElemMolType& mt = bng_data.get_elem_mol_type(mi.elem_mol_type_id);

            mt.set_flag(SPECIES_CPLX_MOL_FLAG_COMPARTMENT_USED_IN_RXNS);
            // collect compartments (ignore in and out)
            mt.reactant_compartments.insert(reac.get_compartment_id(true));
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

  // and also update all existing species
  for (Species* s: all_species.get_species_vector()) {
    s->update_flag_and_compartments_used_in_rxns();
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


void RxnContainer::compute_reacting_classes(const ReactantClass& rc) {
  assert(rc.is_initialized());

  assert(reacting_classes.size() == rc.id);
  reacting_classes.push_back(ReactantClassIdSet());
  ReactantClassIdSet& current_set = reacting_classes.back();

  for (const ReactantClass* reacting_class: reactant_classes_vector) {
    assert(reacting_class->is_initialized());

    // cross check
    // - rule_ids[0] -> reacting_class.rule_ids[1]
    // - rule_ids[1] -> reacting_class.rule_ids[0]
    bool can_react =
        (rc.reaction_id_bitsets[0] & reacting_class->reaction_id_bitsets[1]).any() ||
        (rc.reaction_id_bitsets[1] & reacting_class->reaction_id_bitsets[0]).any();

    if (can_react) {
      // mapping A -> B
      current_set.insert(reacting_class->id);

      // update also the mapping B -> A because this is a new reactant class
      reacting_classes[reacting_class->id].insert(rc.id);
    }
  }
}


// also computes reacting classes if this is a new reactant class
reactant_class_id_t RxnContainer::find_or_add_reactant_class(
    const ReactionIdBitsets& reactions_bitset_per_reactant, const bool target_only) {

  reactant_class_id_t res;
  bool is_new_reactant_class;

  // NOTE: maybe search without creating the ReactantClass object
  ReactantClass* rc = new ReactantClass;
  rc->target_only = target_only;
  rc->reaction_id_bitsets = reactions_bitset_per_reactant;
  const auto it = reactant_classes_set.find(rc);

  if (it != reactant_classes_set.end()) {
    delete rc;
    return (*it)->id;
  }
  else {
    // add and compute reacting classes
    rc->id = next_reactant_class_id;
    next_reactant_class_id++;
    reactant_classes_set.insert(rc);
    reactant_classes_vector.push_back(rc);

    compute_reacting_classes(*rc);
    return rc->id;
  }
}


reactant_class_id_t RxnContainer::compute_reactant_class_for_species(const species_id_t species_id) {

  // prepare reactant class bitsets
  assert(rxn_rules.empty() || rxn_rules.back()->id == rxn_rules.size() - 1);
  ReactionIdBitsets reactions_bitset_per_reactant =
    { boost::dynamic_bitset<>(rxn_rules.size()), boost::dynamic_bitset<>(rxn_rules.size())};
  for (RxnRule* r: rxn_rules) {
    if (r->is_bimol_vol_rxn()) {
      std::vector<uint> indices;
      r->get_reactant_indices(species_id, all_species, indices);
      assert(indices.size() <= 2);

      if (indices.empty()) {
        continue;
      }

      if (indices.size() == 1) {
        // species matches one of reactants
        // set bit on position id
        assert(indices[0] <= 1);
        reactions_bitset_per_reactant[indices[0]][r->id] = 1;
      }
      else {
        // species matches both reactants
        assert(indices[0] + indices[1] == 1);
        reactions_bitset_per_reactant[0][r->id] = 1;
        reactions_bitset_per_reactant[1][r->id] = 1;
      }
    }
  }

  // find or add reactant class based on the computed bitsets, also compute reacting classes
  Species& s = all_species.get(species_id);
  reactant_class_id_t res = find_or_add_reactant_class(reactions_bitset_per_reactant, s.is_target_only());

#if 0
  cout <<
      reactions_bitset_per_reactant[0] << "|" << reactions_bitset_per_reactant[1] <<
      " reactant class id:" << res << ", species id:" << species_id << " (" << rxn_rules.size() << ")\n";
#endif

  return res;
}


// - puts pointers to all corresponding classes to the res_classes_map
// - for bimol rxns, does not reuse already defined rxn class, e.g. when A + B was already created,
//   rxn class for B + A will be created (NOTE: might improve if needed but so far the only issue
//   are reports and printouts
// - called only from get_bimol_rxns_for_reactant
void RxnContainer::create_bimol_rxn_classes_for_new_species(const species_id_t species_id1, const bool for_all_known_species) {

  // find all reactions for species id,
  // also define reactant class
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
        const Species* species = all_species.get_species_vector()[i];
        assert(species != nullptr);
        species_id_t species_id2 = species->id;

        // TODO: simplify condition - is_species_superclass check does not have to be there
        if (!for_all_known_species &&
            !species->was_instantiated()  &&
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

          bool has_rxn_where_both_match_both_patterns = false;

          for (RxnRule* matching_rxn: rxns_for_new_species) {

            // usually the species must be different but reactions of type A + A are allowed
            if (reac1.species_id == species_id2 && !matching_rxn->species_is_both_bimol_reactants(reac1.species_id, all_species)) {
              continue;
            }

            uint reac1_index, reac2_index;
            bool both_match_both_patterns = false;
            bool reactants_match =
                matching_rxn->species_can_be_bimol_reactants(
                    reac1.species_id, species_id2, all_species,
                    &reac1_index, &reac2_index, &both_match_both_patterns);

            if (both_match_both_patterns) {
              has_rxn_where_both_match_both_patterns = true;
            }

            if (reactants_match &&
                matching_rxn->reactant_compatment_matches(reac1_index, reac1.compartment_id) &&
                matching_rxn->reactant_compatment_matches(reac2_index, reac2.compartment_id)) {


              // ok, we have a reaction applicable both to new_id and second_id and compartment matches as well
              // we need to add this rxn to a rxn class for these reactants
              applicable_rxns.push_back(matching_rxn);
            }
          }

          if (!applicable_rxns.empty()) {

            Reactant reacA = reac1;
            Reactant reacB = reac2;

            // if both reactants match both patterns, and compartments are not important,
            // we must sort the reactants so that it does not matter whether the rxn class was first created
            // for reac1 or reac2,
            // otherwise, this might give us different result e.g. based on the frequency of rxn class cleanup
            // NOTE: some cases related to compartments might be missed
            if (reac1.species_id != reac2.species_id &&
                has_rxn_where_both_match_both_patterns &&
                (
                 ((reac1.compartment_id == COMPARTMENT_ID_NONE || reac1.compartment_id == COMPARTMENT_ID_ANY) &&
                  (reac2.compartment_id == COMPARTMENT_ID_NONE || reac2.compartment_id == COMPARTMENT_ID_ANY)) ||
                 reac1.compartment_id == reac2.compartment_id
                )
               ) {
              // order by species name - the one lexicographically smaller will be the first one
              const Species& s1 = all_species.get(reac1.species_id);
              const Species& s2 = all_species.get(reac2.species_id);
              if (s1.name > s2.name) {
                reacA = reac2;
                reacB = reac1;
              }
            }

            // get to the instance of the reaction class for (new_id, second_id) or (second_id, new_id)
            RxnClass* rxn_class = get_or_create_empty_bimol_rxn_class(reacA, reacB);

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
  assert(rxn_classes.count(rxn_class) != 0);
  rxn_classes.erase(rxn_class);

  // destructor also removes existing links from rxn rules to this rxn class
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
  //if (species_processed_for_bimol_rxn_classes.count(reac1_species_id) != 0) {
  if (bimol_rxn_class_map.count(reac1_species_id) != 0) {
    // forget that we processed this species
    // does not necessarily have to be present - we might have not fully processed this species,
    // it might present in the bimol_rxn_class_map only because it was used as a second reactant for
    // another fully processed species
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
          if (it_species2_comp2_species1 != pair_species2_comp2.second.end()) {
            // reverse mapping does not have to exist for all compartments
            pair_species2_comp2.second.erase(it_species2_comp2_species1);
          }
        }
      }

      // erase the top level entry
      bimol_rxn_class_map.erase(it_species1);
    }
  }
}


void RxnContainer::remove_species_id_references(const species_id_t id) {
  // this method is currently used only from SpeciesCleanupEvent
  // we are assuming that there are rxn classes that use this species
  // were removed
  release_assert(species_processed_for_unimol_rxn_classes.count(id) == 0);
  release_assert(species_processed_for_bimol_rxn_classes.count(id) == 0);

  for (RxnRule* rxn: rxn_rules) {
    rxn->remove_species_id_references(id);
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


void RxnContainer::print_periodic_stats() const {
#define ITEM_SIZE(a) "RxnContainer: " << #a << " = " << a.size() << "\n"

  uint64_t applicable_total = 0;
  uint64_t not_applicable_total = 0;
  uint64_t rnx_class_users_total = 0;
  for (const RxnRule* rxn: rxn_rules) {
    applicable_total += rxn->species_applicable_as_any_reactant.size();
    applicable_total += rxn->species_applicable_as_reactant[0].size();
    applicable_total += rxn->species_applicable_as_reactant[1].size();

    not_applicable_total += rxn->species_not_applicable_as_any_reactant.size();
    not_applicable_total += rxn->species_not_applicable_as_reactant[0].size();
    not_applicable_total += rxn->species_not_applicable_as_reactant[1].size();

    rnx_class_users_total += rxn->get_rxn_classed_where_used().size();
  }

  std::cout <<
      ITEM_SIZE(rxn_classes) <<
      ITEM_SIZE(species_processed_for_bimol_rxn_classes) <<
      ITEM_SIZE(species_processed_for_unimol_rxn_classes) <<
      ITEM_SIZE(unimol_rxn_class_map) <<
      ITEM_SIZE(bimol_rxn_class_map) <<
      ITEM_SIZE(bimol_rxn_class_any_orient_compartment_map)
      ITEM_SIZE(rxn_rules) <<
      "RxnContainer: rxn_rules - total species applicable as reactant = " <<
        applicable_total << "\n" <<
      "RxnContainer: rxn_rules - total species not applicable as reactant = " <<
        not_applicable_total << "\n"
      "RxnContainer: rxn_rules - total rxn classes where used = " <<
        rnx_class_users_total << "\n";

#undef ITEM_SIZE
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
