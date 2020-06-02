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

#include "mol_or_rxn_count_event.h"

#include <iostream>

#include "world.h"
#include "partition.h"
#include "datamodel_defines.h"

using namespace std;

namespace MCell {

void MolOrRxnCountTerm::dump(const std::string ind) const {

  cout << ind << "type: ";
  switch(type) {
    case CountType::Invalid:
      cout << "Invalid";
      break;
    case CountType::EnclosedInWorld:
      cout << "World";
      break;
    case CountType::EnclosedInObject:
      cout << "EnclosedInObject";
      break;
    case CountType::PresentOnSurfaceRegion:
      cout << "PresentOnSurfaceRegion";
      break;
    case CountType::RxnCountInWorld:
      cout << "RxnCountInWorld";
      break;
    case CountType::RxnCountInObject:
      cout << "RxnCountInObject";
      break;
    case CountType::RxnCountOnSurfaceRegion:
      cout << "RxnCountOnSurfaceRegion";
      break;
    default:
      assert(false);
  }
  cout << "\n";

  cout << ind << "sign_in_expression: " << sign_in_expression << " [int]\n";
  cout << ind << "orientation: " << orientation << " [orientation_t]\n";
  cout << ind << "species_id: " << species_id << " [species_id_t]\n";
  cout << ind << "rxn_rule_id: " << rxn_rule_id << " [rxn_rule_id_t]\n";
  cout << ind << "geometry_object_id: " << geometry_object_id << " [geometry_object_id_t]\n";
}


std::string MolOrRxnCountTerm::to_data_model_string(const World* world, bool print_positive_sign) const {
  stringstream res;
  assert(sign_in_expression == -1 || sign_in_expression == 1);
  res <<
      ((sign_in_expression == 1) ?
        (print_positive_sign) ? "+" : "":
        "-"
  );

  res << "COUNT[";

  switch(type) {
    case CountType::EnclosedInWorld:
    case CountType::EnclosedInObject:
    case CountType::PresentOnSurfaceRegion:
      res << world->get_all_species().get(species_id).name;
      break;
    case CountType::RxnCountInWorld:
    case CountType::RxnCountInObject:
    case CountType::RxnCountOnSurfaceRegion: {
        string rxn_name = world->get_all_rxns().get_rxn_rule(rxn_rule_id)->name;
        CONVERSION_CHECK(rxn_name != "", "Counted reaction has no name");
        res << rxn_name;
      }
      break;
    default:
      assert(false);
  }

  res << ",";

  switch(type) {
    case CountType::EnclosedInWorld:
    case CountType::RxnCountInWorld:
      res << VALUE_WORLD;
      break;
    case CountType::EnclosedInObject:
    case CountType::RxnCountInObject: {
        string obj_name = world->get_geometry_object(geometry_object_id).name;
        CONVERSION_CHECK(obj_name != "", "Counted object has no name");
        res << obj_name;
      }
      break;
    case CountType::PresentOnSurfaceRegion:
    case CountType::RxnCountOnSurfaceRegion:{
      string reg_name = world->get_region(region_id).name;
      CONVERSION_CHECK(reg_name != "", "Counted region has no name");
      res << DMUtil::get_object_w_region_name(reg_name, false);
    }
    break;

    default:
      assert(false);
  }
  res << "]";
  return res.str();
}


void MolOrRxnCountInfo::dump(const std::string ind) const {

  cout << ind << "buffer_id: " << buffer_id << " [count_buffer_id_t] \t\t\n";
  cout << ind << "terms:\n";
  for (uint i = 0; i < terms.size(); i++) {
    cout << i << ":\n";
    terms[i].dump(ind + "  ");
  }
}

static string basename(const string& s) {
   char sep = '/'; // / is used on Windows as well in the rxn_output file path
   size_t i = s.rfind(sep, s.length());
   if (i != string::npos) {
      return(s.substr(i+1, s.length() - i));
   }
   else {
     return s;
   }
}

static string noext(const string& s) {
   char sep = '.'; // / is used on Windows as well in the rxn_output file path
   size_t i = s.rfind(sep, s.length());
   if (i != string::npos) {
      return(s.substr(0, i));
   }
   else {
     return s;
   }
}

void MolOrRxnCountInfo::to_data_model(const World* world, Json::Value& reaction_output) const {
  DMUtil::json_add_version(reaction_output, JSON_DM_VERSION_1330);

  // MDLString is a general way how to capture the output
  reaction_output[KEY_RXN_OR_MOL] = VALUE_MDLSTRING;
  reaction_output[KEY_DESCRIPTION] = "";
  reaction_output[KEY_PLOTTING_ENABLED] = true;

  // file prefix
  const CountBuffer& buff = world->get_count_buffer(buffer_id);
  string prefix;
  string filename = basename(buff.get_filename());
  prefix = noext(filename);
  reaction_output[KEY_MDL_FILE_PREFIX] = prefix;

  // species or rxn name & location is in mdl_string
  reaction_output[KEY_NAME] = "";
  reaction_output[KEY_MOLECULE_NAME] = "";
  reaction_output[KEY_REACTION_NAME] = "";
  reaction_output[KEY_OBJECT_NAME] = "";
  reaction_output[KEY_REGION_NAME] = "";

  if (terms.size() == 1) {
    string val;
    switch (terms[0].type) {
      case CountType::EnclosedInWorld:
      case CountType::RxnCountInWorld:
        val = VALUE_COUNT_LOCATION_WORLD;
        break;
      case CountType::EnclosedInObject:
      case CountType::RxnCountInObject:
        val = VALUE_COUNT_LOCATION_OBJECT;
        break;
      case CountType::PresentOnSurfaceRegion:
      case CountType::RxnCountOnSurfaceRegion:
        val = VALUE_COUNT_LOCATION_REGION;
        break;
      default:
        assert(false);
    }
    reaction_output[KEY_COUNT_LOCATION] = val;
  }
  else {
    #ifndef NDEBUG
        cerr <<
            "Warning: conversion of count expression with multiple terms to data model's " << KEY_COUNT_LOCATION <<
            " is not fully supported, defaulting to " << VALUE_COUNT_LOCATION_REGION << ".\n";
    #endif
    reaction_output[KEY_COUNT_LOCATION] = VALUE_COUNT_LOCATION_REGION;
  }

  // handled by prefix
  reaction_output[KEY_DATA_FILE_NAME] = "";

  string mdl_string = "";
  for (size_t i = 0; i < terms.size(); i++) {
    mdl_string += terms[i].to_data_model_string(world, i != 0);
  }
  reaction_output[KEY_MDL_STRING] = mdl_string;
}


template<class T>
static uint sum_all_map_items(const T& map) {
  uint res = 0;
  for (const auto it: map) {
    res += it.second;
  }
  return res;
}

static bool wall_is_part_of_region(
    const Partition& p, const wall_index_t wall_index, const region_id_t region_id) {

  assert(wall_index != WALL_INDEX_INVALID);
  assert(region_id != REGION_ID_INVALID);

  // for all the regions of the current molecule
  const Wall& w = p.get_wall(wall_index);
  for (region_index_t ri: w.regions) {
    if (p.get_region(ri).id == region_id) {
      return true;
    }
  }
  return false;
}


// TODO: refactor, too many levels of control logic
void MolOrRxnCountEvent::step() {

  // go through all molecules and count them
  PartitionVector& partitions = world->get_partitions();
  std::vector<CountItem> count_items;
  count_items.resize(mol_count_infos.size());

  // initialize new count items
  for (uint i = 0; i < mol_count_infos.size(); i++) {
    count_items[i].time = event_time * world->config.time_unit;
    count_items[i].int_value = 0;
  }

  species_id_t all_mol_id = world->get_all_species().get_all_molecules_species_id();
  species_id_t all_vol_id = world->get_all_species().get_all_volume_molecules_species_id();
  species_id_t all_surf_id = world->get_all_species().get_all_surface_molecules_species_id();

  // for each partition
  for (Partition& p: partitions) {

    // for each molecule
    // TODO: optimize - we do not need this if we are counting just reactions
    for (const Molecule& m: p.get_molecules()) {

      if (m.is_defunct()) {
        continue;
      }

      const BNG::Species& species = world->get_all_species().get(m.species_id);
      if (!species.has_count_enclosed_flag() && !species.has_count_contents_flag()) {
        continue;
      }

      // volumes that enclose the current molecule
      // might be nullptr if not found or mol is surface mol
      const uint_set<geometry_object_id_t>* enclosing_volumes = nullptr;
      if (m.is_vol()) {
        enclosing_volumes = p.get_enclosing_counted_volumes(m.v.counted_volume_id);
      }

      // for each counting info
      for (uint i = 0; i < mol_count_infos.size(); i++) {
        const MolOrRxnCountInfo& info = mol_count_infos[i];

        for (const MolOrRxnCountTerm& term: info.terms) {
          assert(!term.is_mol_count() || term.species_id != SPECIES_ID_INVALID);

          // does the current molecule match?
          if (term.is_mol_count() &&
              ( term.species_id == m.species_id ||
                term.species_id == all_mol_id ||
               (term.species_id == all_vol_id && m.is_vol()) ||
               (term.species_id == all_surf_id && m.is_surf())
              )
          ) {

            if (term.type == CountType::EnclosedInWorld) {
              // count the molecule
              count_items[i].inc_or_dec(term.sign_in_expression);
            }
            else if (m.is_vol() && term.type == CountType::EnclosedInObject) {
              assert(term.geometry_object_id != GEOMETRY_OBJECT_ID_INVALID);

              // is the molecule inside of the object that we are checking?
              if (m.v.counted_volume_id == term.geometry_object_id ||
                  (enclosing_volumes != nullptr && enclosing_volumes->count(term.geometry_object_id))) {

                count_items[i].inc_or_dec(term.sign_in_expression);
              }
            }
            else if (m.is_surf() && term.type == CountType::PresentOnSurfaceRegion) {
              assert(term.region_id != REGION_ID_INVALID);

              if (wall_is_part_of_region(p, m.s.wall_index, term.region_id)) {
                count_items[i].inc_or_dec(term.sign_in_expression);
              }
            }
          }
        } // for terms
      } // for mol_count_infos
    } // for molecules

    // TODO: optimize - we do not need this if we are counting just molecules
    for (const BNG::RxnRule* rxn: world->get_all_rxns().get_rxn_rules_vector()) {
      if (!rxn->is_counted()) {
        continue;
      }

      // for each counting info
      for (uint i = 0; i < mol_count_infos.size(); i++) {
        const MolOrRxnCountInfo& info = mol_count_infos[i];

        for (const MolOrRxnCountTerm& term: info.terms) {
          assert(!term.is_rxn_count() || term.rxn_rule_id != BNG::RXN_RULE_ID_INVALID);

          // does the current reaction match?
          if (term.is_rxn_count() && term.rxn_rule_id == rxn->id) {

            // get counts from partition
            const CountInGeomObjectMap& counts_in_objects = p.get_rxn_in_volume_count_map(term.rxn_rule_id);
            const CountOnWallMap& counts_on_walls = p.get_rxn_on_surface_count_map(term.rxn_rule_id);

            if (term.type == CountType::RxnCountInWorld) {
              // count all the occurrences

              count_items[i].inc_or_dec(
                  term.sign_in_expression,
                  sum_all_map_items(counts_in_objects)
              );
              count_items[i].inc_or_dec(
                  term.sign_in_expression,
                  sum_all_map_items(counts_on_walls)
              );
            }
            else if (term.type == CountType::RxnCountInObject) {
              assert(term.geometry_object_id != GEOMETRY_OBJECT_ID_INVALID);

              for (const auto it: counts_in_objects) {
                if (it.second != 0) {
                  geometry_object_id_t obj_id_for_this_count = it.first;

                  // volumes that enclose the location where reaction occurred
                  // might be nullptr if not found
                  const uint_set<geometry_object_id_t>* enclosing_volumes =
                      p.get_enclosing_counted_volumes(obj_id_for_this_count);

                  // did the reaction occur in the object that we are checking?
                  if (obj_id_for_this_count == term.geometry_object_id ||
                      (enclosing_volumes != nullptr && enclosing_volumes->count(term.geometry_object_id))
                  ) {
                    count_items[i].inc_or_dec(term.sign_in_expression, it.second);
                  }
                }
              }
            }
            else if (term.type == CountType::RxnCountOnSurfaceRegion) {
              assert(term.region_id != REGION_ID_INVALID);

              for (const auto it: counts_on_walls) {
                if (it.second != 0) {
                  wall_index_t wall_index_for_this_count = it.first;

                  // did the reaction occur on the surface region that we are checking?
                  if (wall_is_part_of_region(p, wall_index_for_this_count, term.region_id)) {
                    count_items[i].inc_or_dec(term.sign_in_expression, it.second);
                  }
                }
              }
            }
          } // for terms
        } // for mol_count_infos
      } // for molecules

    } // for rxn rule

    // we are not reseting counts for another counting event because the result is a
    // sum of all previous iterations

  } // for partition

  // check each molecule against what we are checking
  for (uint i = 0; i < mol_count_infos.size(); i++) {
    world->get_count_buffer(mol_count_infos[i].buffer_id).add(count_items[i]);
  }

}


void MolOrRxnCountEvent::dump(const std::string ind) const {
  cout << ind << "Mol count event:\n";
  std::string ind2 = ind + "  ";
  std::string ind4 = ind2 + "  ";
  BaseEvent::dump(ind2);

  cout << ind << " mol_count_infos:\n";
  for(uint i = 0; i < mol_count_infos.size(); i++) {
    cout << ind2 << i << "\n";
    mol_count_infos[i].dump(ind4);
  }

}


void MolOrRxnCountEvent::to_data_model(Json::Value& mcell_node) const {
  // set global reaction_data_output settings if this is the first
  // counting event that we are converting
  if (mol_count_infos.empty()) {
    return;
  }

  Json::Value& reaction_data_output = mcell_node[KEY_REACTION_DATA_OUTPUT];
  Json::Value& reaction_output_list = reaction_data_output[KEY_REACTION_OUTPUT_LIST];

  for (const MolOrRxnCountInfo& info: mol_count_infos) {
    Json::Value reaction_output;
    info.to_data_model(world, reaction_output);
    reaction_output_list.append(reaction_output);
  }
}

}
