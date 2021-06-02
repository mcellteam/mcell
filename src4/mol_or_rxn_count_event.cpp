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

#include "generated/gen_names.h"

using namespace std;

namespace MCell {

uint MolOrRxnCountTerm::get_num_pattern_matches(const species_id_t species_id) const {
  assert(species_pattern_type != SpeciesPatternType::SpeciesId);

  auto it = species_ids_matching_pattern_w_multiplier_cache.find(species_id);
  if (it == species_ids_matching_pattern_w_multiplier_cache.end()) {
    return 0;
  }
  else if (species_pattern_type == SpeciesPatternType::SpeciesPattern) {
    // each molecule is counted only once
    return 1;
  }
  else {
    assert(species_pattern_type == SpeciesPatternType::MoleculesPattern);
    assert(it->second >= 1);
    return it->second;
  }
}


uint MolOrRxnCountTerm::get_num_molecule_matches(
    const Molecule& m,
    const species_id_t all_mol_id, const species_id_t all_vol_id, const species_id_t all_surf_id) const {
  assert(is_mol_count());
  if (species_pattern_type == SpeciesPatternType::SpeciesId) {
    assert(species_id != SPECIES_ID_INVALID);
    if (
        species_id == m.species_id ||
        species_id == all_mol_id ||
       (species_id == all_vol_id && m.is_vol()) ||
       (species_id == all_surf_id && m.is_surf())
    ) {
      return 1;
    }
    else {
      return 0;
    }
  }
  else {
    assert(!species_molecules_pattern.elem_mols.empty());
    return get_num_pattern_matches(m.species_id);
  }
}


void MolOrRxnCountTerm::dump(const std::string ind) const {

  cout << ind << "type: ";
  switch(type) {
    case CountType::Invalid:
      cout << "Invalid";
      break;
    case CountType::EnclosedInWorld:
      cout << "World";
      break;
    case CountType::EnclosedInVolumeRegion:
      cout << "EnclosedInObject";
      break;
    case CountType::PresentOnSurfaceRegion:
      cout << "PresentOnSurfaceRegion";
      break;
    case CountType::RxnCountInWorld:
      cout << "RxnCountInWorld";
      break;
    case CountType::RxnCountInVolumeRegion:
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

  cout << ind << "type: ";
  switch(species_pattern_type) {
    case SpeciesPatternType::Invalid:
      cout << "Invalid";
      break;
    case SpeciesPatternType::SpeciesId:
      cout << "SpeciesId";
      cout << ind << "species_id: " << species_id << " [species_id_t]\n";
      break;
    case SpeciesPatternType::SpeciesPattern:
      cout << "SpeciesPattern";
      cout << ind << "species_molecules_pattern: " << species_molecules_pattern.to_str() << " [CplxInstance]\n";
      break;
    case SpeciesPatternType::MoleculesPattern:
      cout << "MoleculesPattern";
      cout << ind << "species_molecules_pattern: " << species_molecules_pattern.to_str() << " [CplxInstance]\n";
      break;

    default:
      assert(false);
  }
  cout << "\n";

  cout << ind << "primary_compartment_id: " << primary_compartment_id << " [compartment_id_t]\n";
  cout << ind << "rxn_rule_id: " << rxn_rule_id << " [rxn_rule_id_t]\n";
  if (region_expr.root != nullptr) {
    cout << ind << "region_expr: " << region_expr.root->to_string(nullptr, false) << " [geometry_object_id_t]\n";
  }
  else {
    cout << ind << "region_expr: " << "none\n";
  }
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

  string pattern;
  switch(type) {
    case CountType::EnclosedInWorld:
    case CountType::EnclosedInVolumeRegion:
    case CountType::PresentOnSurfaceRegion:
      if (primary_compartment_id != BNG::COMPARTMENT_ID_NONE) {
        pattern = "@" + world->bng_engine.get_data().get_compartment(primary_compartment_id).name + ":";
      }

      switch (species_pattern_type) {
        case SpeciesPatternType::SpeciesId:
          pattern += world->get_all_species().get(species_id).name;
          break;
        case SpeciesPatternType::SpeciesPattern:
          pattern += species_molecules_pattern.to_str() + MARKER_SPECIES_COMMENT;
          break;
        case SpeciesPatternType::MoleculesPattern:
          pattern += species_molecules_pattern.to_str() + MARKER_MOLECULES_COMMENT;
          break;
        default:
          assert(false);
      }
      break;
    case CountType::RxnCountInWorld:
    case CountType::RxnCountInVolumeRegion:
    case CountType::RxnCountOnSurfaceRegion: {
        string rxn_name = world->get_all_rxns().get(rxn_rule_id)->name;
        CONVERSION_CHECK(rxn_name != "", "Counted reaction has no name");
        pattern = rxn_name;
      }
      break;
    default:
      assert(false);
  }

  string where;
  switch(type) {
    case CountType::EnclosedInWorld:
    case CountType::RxnCountInWorld:
      where = VALUE_WORLD;
      break;
    case CountType::EnclosedInVolumeRegion:
    case CountType::RxnCountInVolumeRegion: {
        if (region_expr.root->op != RegionExprOperator::LEAF_GEOMETRY_OBJECT) {
          errs() << "Data model export with region expressions is not supported yet, use " <<
              API::NAME_CLASS_MODEL << "." << API::NAME_EXPORT_VIZ_DATA_MODEL << " instead, error for " <<
              API::NAME_CLASS_COUNT << " for " <<
              ((type == CountType::EnclosedInVolumeRegion) ?
                  species_molecules_pattern.to_str() :
                  world->get_all_rxns().get(rxn_rule_id)->name) << ".\n";
          exit(1);
        }
        string obj_name = world->get_geometry_object(region_expr.root->geometry_object_id).name;
        CONVERSION_CHECK(obj_name != "", "Counted object has no name");
        where = obj_name;
      }
      break;
    case CountType::PresentOnSurfaceRegion:
    case CountType::RxnCountOnSurfaceRegion:{
      if (region_expr.root->op != RegionExprOperator::LEAF_SURFACE_REGION) {
        errs() << "Data model export with region expressions is not supported yet, use " <<
            API::NAME_CLASS_MODEL << "." << API::NAME_EXPORT_VIZ_DATA_MODEL << " instead, error for " <<
            API::NAME_CLASS_COUNT << " for " <<
            ((type == CountType::PresentOnSurfaceRegion) ?
                species_molecules_pattern.to_str() :
                world->get_all_rxns().get(rxn_rule_id)->name) << ".\n";
        exit(1);
      }
      region_id_t region_id = region_expr.root->region_id;
      // if possible, we should use 2d compartment name because in MCell
      // the encompassing region name is the same as the name of the object,
      // this would lead to the same name for A@CP and A@PM
      bool surf_compartment_used = false;
      const Partition& p = world->get_partition(PARTITION_ID_INITIAL);
      for (const auto& geom_obj: p.get_geometry_objects()) {
        if (geom_obj.encompassing_region_id == region_id &&
            geom_obj.surf_compartment_id != BNG::COMPARTMENT_ID_NONE) {

          const BNG::Compartment& comp = world->bng_engine.get_data().get_compartment(geom_obj.surf_compartment_id);
          assert(!comp.is_3d);
          pattern = "@" + comp.name + ":" + pattern;

          where = VALUE_WORLD;

          surf_compartment_used = true;
          break;
        }
      }

      if (!surf_compartment_used) {
        string reg_name = world->get_region(region_id).name;
        CONVERSION_CHECK(reg_name != "", "Counted region has no name");
        where = DMUtils::get_object_w_region_name(reg_name, false);
      }
    }
    break;

    default:
      assert(false);
  }
  res << pattern << "," << where << "]";
  return res.str();
}


bool MolOrRxnCountItem::counts_mols() const {
  for (auto& t: terms) {
    if (t.is_mol_count()) {
      return true;
    }
  }
  return false;
}


bool MolOrRxnCountItem::counts_rxns() const {
  for (auto& t: terms) {
    if (t.is_rxn_count()) {
      return true;
    }
  }
  return false;
}


bool MolOrRxnCountItem::is_world_mol_count() const {
  // SpeciesId is obsoleted and not handled correctly
  return
      terms.size() == 1 &&
      terms[0].type == CountType::EnclosedInWorld &&
      terms[0].species_pattern_type != SpeciesPatternType::SpeciesId;
}


void MolOrRxnCountItem::dump(const std::string ind) const {

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

void MolOrRxnCountItem::to_data_model(const World* world, Json::Value& reaction_output) const {
  DMUtils::add_version(reaction_output, VER_DM_2018_01_11_1330);

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
      case CountType::EnclosedInVolumeRegion:
      case CountType::RxnCountInVolumeRegion:
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

  if (multiplier != 1) {
    mdl_string = "(" + mdl_string + ")*" + f_to_str(multiplier);
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


static bool counted_volume_matches_region_expr_recursively(
    const CountedVolume& enclosing_volumes,
    const RegionExprNode* node
) {
  assert(node != nullptr);

  if (node->op == RegionExprOperator::LEAF_GEOMETRY_OBJECT) {
    // returns true if the current geometry object is one of the enclosing volumes
    return enclosing_volumes.contained_in_objects.count(node->geometry_object_id) != 0;
  }
  else if (node->has_binary_op()) {
    bool left = counted_volume_matches_region_expr_recursively(enclosing_volumes, node->left);
    bool right = counted_volume_matches_region_expr_recursively(enclosing_volumes, node->right);
    switch (node->op) {
      case RegionExprOperator::DIFFERENCE:
        return left && !right;
      case RegionExprOperator::INTERSECT:
        return left && right;
      case RegionExprOperator::UNION:
        return left || right;
      default:
        assert(false);
        return false;
    }
  }
  else {
    assert(false);
    return false;
  }
}


void MolOrRxnCountEvent::compute_mol_count_item(
    const Partition& p,
    const MolOrRxnCountItem& item,
    const Molecule& m,
    CountItemVector& count_items
) {
  species_id_t all_mol_id = world->get_all_species().get_all_molecules_species_id();
  species_id_t all_vol_id = world->get_all_species().get_all_volume_molecules_species_id();
  species_id_t all_surf_id = world->get_all_species().get_all_surface_molecules_species_id();

  for (const MolOrRxnCountTerm& term: item.terms) {

    // does the current molecule match?
    if (term.is_mol_count()) {

      // num_matches may be > 1 when molecules pattern is used for matching
      uint num_matches = term.get_num_molecule_matches(m, all_mol_id, all_vol_id, all_surf_id);
      if (num_matches == 0) {
        continue;
      }

      if (term.type == CountType::EnclosedInWorld) {
        // count the molecule
        count_items[item.index].inc_or_dec(term.sign_in_expression, num_matches);
      }
      else if (m.is_vol() && term.type == CountType::EnclosedInVolumeRegion) {

        // is the molecule inside of the object/volume region expression that we are checking?
        const CountedVolume& enclosing_volumes = p.get_counted_volume(m.v.counted_volume_index);
        if (counted_volume_matches_region_expr_recursively(enclosing_volumes, term.region_expr.root)) {
          count_items[item.index].inc_or_dec(term.sign_in_expression, num_matches);
        }
      }
      else if (m.is_surf() && term.type == CountType::PresentOnSurfaceRegion) {
        // TODO_COUNTS
        release_assert(term.region_expr.root->op == RegionExprOperator::LEAF_SURFACE_REGION);
        region_id_t region_id = term.region_expr.root->region_id;
        assert(region_id != REGION_ID_INVALID);

        if (wall_is_part_of_region(p, m.s.wall_index, region_id)) {
          count_items[item.index].inc_or_dec(term.sign_in_expression, num_matches);
        }
      }
    }
  } // for terms
}


void MolOrRxnCountEvent::compute_rxn_count_item(
    Partition& p,
    const MolOrRxnCountItem& item,
    const BNG::RxnRule* rxn,
    CountItemVector& count_items
) {
  for (const MolOrRxnCountTerm& term: item.terms) {
    assert(!term.is_rxn_count() || term.rxn_rule_id != BNG::RXN_RULE_ID_INVALID);

    // does the current reaction match?
    if (term.is_rxn_count() && term.rxn_rule_id == rxn->id) {

      // use initial_reactions_count (from previous checkpoint)
      count_items[item.index].value += term.initial_reactions_count;

      // get counts from partition
      const CountInGeomObjectMap& counts_in_objects = p.get_rxn_in_volume_count_map(term.rxn_rule_id);
      const CountOnWallMap& counts_on_walls = p.get_rxn_on_surface_count_map(term.rxn_rule_id);

      if (term.type == CountType::RxnCountInWorld) {
        // count all the occurrences

        count_items[item.index].inc_or_dec(
            term.sign_in_expression,
            sum_all_map_items(counts_in_objects)
        );
        count_items[item.index].inc_or_dec(
            term.sign_in_expression,
            sum_all_map_items(counts_on_walls)
        );
      }
      else if (term.type == CountType::RxnCountInVolumeRegion) {

        // TODO_COUNTS
        release_assert(term.region_expr.root->op == RegionExprOperator::LEAF_GEOMETRY_OBJECT);
        geometry_object_id_t geometry_object_id = term.region_expr.root->geometry_object_id;
        assert(geometry_object_id != GEOMETRY_OBJECT_ID_INVALID);

        for (const auto it: counts_in_objects) {
          if (it.second != 0) {
            counted_volume_index_t counted_volume_index_for_this_count = it.first;

            // volumes that enclose the location where reaction occurred
            // might be nullptr if not found
            const CountedVolume& enclosing_volumes = p.get_counted_volume(counted_volume_index_for_this_count);

            // did the reaction occur in the object that we are checking?
            if (enclosing_volumes.contained_in_objects.count(geometry_object_id) != 0) {
              count_items[item.index].inc_or_dec(term.sign_in_expression, it.second);
            }
          }
        }
      }
      else if (term.type == CountType::RxnCountOnSurfaceRegion) {
        // TODO_COUNTS
        release_assert(term.region_expr.root->op == RegionExprOperator::LEAF_SURFACE_REGION);
        region_id_t region_id = term.region_expr.root->region_id;
        assert(region_id != REGION_ID_INVALID);

        for (const auto it: counts_on_walls) {
          if (it.second != 0) {
            wall_index_t wall_index_for_this_count = it.first;

            // did the reaction occur on the surface region that we are checking?
            if (wall_is_part_of_region(p, wall_index_for_this_count, region_id)) {
              count_items[item.index].inc_or_dec(term.sign_in_expression, it.second);
            }
          }
        }
      }
    }
  } // for terms
}


void MolOrRxnCountEvent::compute_counts(CountItemVector& count_items) {

  // go through all molecules and count them
  PartitionVector& partitions = world->get_partitions();
  count_items.resize(mol_rxn_count_items.size());

  // initialize new count items
  for (uint i = 0; i < mol_rxn_count_items.size(); i++) {
    count_items[i].time = event_time * world->config.time_unit;
    count_items[i].value = 0;
  }

  uint_set<uint> processed_item_indices;

  // to improve performance, first process the species that we are counting in the whole world,
  // but only if there is less species than molecules (very crude heuristics)
  if (world->get_all_species().get_species_vector().size() <
      world->get_partition(PARTITION_ID_INITIAL).get_molecules().size()) {

    for (const BNG::Species* species: world->get_all_species().get_species_vector()) {
      if (species->is_defunct() || species->get_num_instantiations() == 0) {
        continue;
      }

      const CountSpeciesInfo& species_info = get_or_compute_count_species_info(species->id);
      if (species_info.type != CountSpeciesInfoType::Counted) {
        assert(species_info.type == CountSpeciesInfoType::NotCounted);
        continue;
      }

      for (uint index: species_info.world_count_item_indices) {
        const MolOrRxnCountItem& item = mol_rxn_count_items[index];
        assert(item.is_world_mol_count());

        const MolOrRxnCountTerm& term = item.terms[0];
        uint multiplier = term.get_num_pattern_matches(species->id);
        assert(multiplier != 0);

        // use the count of this species as tracked in the BNG library
        count_items[item.index].inc_or_dec(
            term.sign_in_expression, species->get_num_instantiations() * multiplier);

        processed_item_indices.insert(index);
      }
    }
  }

  // for each partition
  for (Partition& p: partitions) {

    // for each molecule (if we are counting them and we did not process all of them already)
    if (count_mols && processed_item_indices.size() != mol_rxn_count_items.size()) {
      for (const Molecule& m: p.get_molecules()) {

        if (m.is_defunct()) {
          continue;
        }

        // check whether we are counting these species at all
        const CountSpeciesInfo& species_info = get_or_compute_count_species_info(m.species_id);
        if (species_info.type != CountSpeciesInfoType::Counted) {
          assert(species_info.type == CountSpeciesInfoType::NotCounted);
          continue;
        }

        // for each counting info
        for (uint i = 0; i < mol_rxn_count_items.size(); i++) {
          // skip already processed count items
          if (processed_item_indices.count(i) != 0) {
            continue;
          }

          compute_mol_count_item(p, mol_rxn_count_items[i], m, count_items);
        }
      } // for molecules
    }

    if (count_rxns) {
      for (const BNG::RxnRule* rxn: world->get_all_rxns().get_rxn_rules_vector()) {
        if (!rxn->is_counted()) {
          continue;
        }

        // for each counting info
        for (uint i = 0; i < mol_rxn_count_items.size(); i++) {
          compute_rxn_count_item(p, mol_rxn_count_items[i], rxn, count_items);
        }
      } // for rxns
    }
  } // for partition

  for (uint i = 0; i < mol_rxn_count_items.size(); i++) {
    // multiply the results by a constant
    // (value is different from 1 when the count expression in MDL had a top level multiplication)
    count_items[i].value *= mol_rxn_count_items[i].multiplier;
  }
}


void MolOrRxnCountEvent::step() {
  CountItemVector count_items;

  compute_counts(count_items);

  // store the counts to buffers
  for (uint i = 0; i < mol_rxn_count_items.size(); i++) {
    world->get_count_buffer(mol_rxn_count_items[i].buffer_id).add(count_items[i]);
  }
}


double MolOrRxnCountEvent::get_single_count_value() {
  CountItemVector count_items;

  compute_counts(count_items);
  assert(count_items.size() == 1);

  return count_items[0].value;
}


void MolOrRxnCountEvent::compute_count_species_info(const species_id_t species_id) {
  CountSpeciesInfo& info = count_species_info[species_id];
  assert(info.type == CountSpeciesInfoType::NotSeen);

  // we saw this species, might be overwritten
  info.type = CountSpeciesInfoType::NotCounted;

  // for each counting info
  for (MolOrRxnCountItem& count_item: mol_rxn_count_items) {

    // and each term
    for (MolOrRxnCountTerm& term: count_item.terms) {
      if (term.is_rxn_count()) {
        // reactions are analyzed in separately in Species::update_rxn_and_custom_flags
        // to determine whether species should have SPECIES_FLAG_NEEDS_COUNTED_VOLUME flag
        continue;
      }
      assert(term.is_mol_count());

      bool matches = false;
      if (term.species_pattern_type == SpeciesPatternType::SpeciesId) {
        if (term.species_id == species_id) {
          matches = true;
        }
      }
      else {
        assert(term.species_pattern_type == SpeciesPatternType::SpeciesPattern ||
            term.species_pattern_type == SpeciesPatternType::MoleculesPattern);

        const BNG::Species& species = world->get_all_species().get(species_id);
        uint num_matches = species.get_pattern_num_matches(term.species_molecules_pattern);

        // also the primary compartment id must match
        if (term.primary_compartment_id != BNG::COMPARTMENT_ID_NONE &&
            term.primary_compartment_id != species.get_primary_compartment_id()) {
          num_matches = 0;
        }

        if (num_matches > 0) {
          // we must also remember that this species id matches the term's pattern
          term.species_ids_matching_pattern_w_multiplier_cache[species_id] = num_matches;
          matches = true;
        }
      }

      if (matches) {
        info.type = CountSpeciesInfoType::Counted;
        if (count_item.is_world_mol_count()) {
          // optimization for faster counting
          info.world_count_item_indices.insert(count_item.index);
        }
        else {
          // there are also other count items besides those in the world_count_item_indices set
          info.all_are_world_mol_counts = false;
        }
      }
    } // for count_item.terms
  } // for mol_count_items
}


inline const CountSpeciesInfo& MolOrRxnCountEvent::get_or_compute_count_species_info(const species_id_t species_id) {
  assert(species_id != SPECIES_ID_INVALID);
  if (species_id >= count_species_info.size()) {
    // extend the array if we did not see these species yet
    count_species_info.resize(species_id + 1);
  }
  CountSpeciesInfo& res = count_species_info[species_id];
  if (res.type == CountSpeciesInfoType::NotSeen) {
    // updates res
    compute_count_species_info(species_id);
  }
  assert(res.type != CountSpeciesInfoType::NotSeen);
  return res;
}


void MolOrRxnCountEvent::dump(const std::string ind) const {
  cout << ind << "MolOrRxnCountEvent:\n";
  std::string ind2 = ind + "  ";
  std::string ind4 = ind2 + "  ";
  BaseEvent::dump(ind2);

  cout << ind << " mol_count_infos:\n";
  for(uint i = 0; i < mol_rxn_count_items.size(); i++) {
    cout << ind2 << i << "\n";
    mol_rxn_count_items[i].dump(ind4);
  }

}


void MolOrRxnCountEvent::to_data_model(Json::Value& mcell_node) const {
  // set global reaction_data_output settings if this is the first
  // counting event that we are converting
  if (mol_rxn_count_items.empty()) {
    return;
  }

  Json::Value& reaction_data_output = mcell_node[KEY_REACTION_DATA_OUTPUT];

  // assuming that there is a single rxn_step value for all
  // (period in second on how often to dump the counted data)
  // at least this is the only supported option in data model right now
  const string& orig_rxn_step = reaction_data_output[KEY_RXN_STEP].asString();
  string new_rxn_step = f_to_str(periodicity_interval * world->config.time_unit);

  if (orig_rxn_step != "" && orig_rxn_step != new_rxn_step) {
    mcell_log(
        "Warning: count events use multiple different time steps, this is not supported "
        "by data model, keeping only value %s (seconds).\n", orig_rxn_step.c_str()
    );
  }
  else {
    reaction_data_output[KEY_RXN_STEP] = new_rxn_step;
  }

  Json::Value& reaction_output_list = reaction_data_output[KEY_REACTION_OUTPUT_LIST];

  for (const MolOrRxnCountItem& info: mol_rxn_count_items) {
    Json::Value reaction_output;
    info.to_data_model(world, reaction_output);
    reaction_output_list.append(reaction_output);
  }
}

}
