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

#ifndef SRC4_MOL_OR_RXN_COUNT_EVENT_H_
#define SRC4_MOL_OR_RXN_COUNT_EVENT_H_

#include "bng/bng.h"

#include "base_event.h"
#include "count_buffer.h"

namespace BNG {
class SpeciesContainer;
}

namespace MCell {

class Partition;
class Molecule;
class World;

enum class SpeciesPatternType {
  Invalid,
  SpeciesId,
  SpeciesPattern,
  MoleculesPattern
};


enum class CountType {
  Invalid,
  EnclosedInWorld,
  EnclosedInVolumeRegion,
  PresentOnSurfaceRegion,
  RxnCountInWorld,
  RxnCountInVolumeRegion,
  RxnCountOnSurfaceRegion,
};


class MolOrRxnCountTerm {
public:
  MolOrRxnCountTerm()
    : type(CountType::Invalid),
      sign_in_expression(0),
      orientation(ORIENTATION_NOT_SET),
      species_pattern_type(SpeciesPatternType::Invalid),
      species_id(SPECIES_ID_INVALID),
      species_molecules_pattern(nullptr),
      rxn_rule_id(BNG::RXN_RULE_ID_INVALID),
      geometry_object_id(GEOMETRY_OBJECT_ID_INVALID),
      region_id(REGION_ID_INVALID)
     {
  }

  bool matches_molecule(
      const Molecule& m,
      const species_id_t all_mol_id, const species_id_t all_vol_id, const species_id_t all_surf_id
  ) const;

  bool is_mol_count() const {
    return type == CountType::EnclosedInWorld || type == CountType::EnclosedInVolumeRegion ||
        type == CountType::PresentOnSurfaceRegion;
  }

  bool is_rxn_count() const {
    return type == CountType::RxnCountInWorld || type == CountType::RxnCountInVolumeRegion ||
        type == CountType::RxnCountOnSurfaceRegion;
  }

  void dump(const std::string ind = "") const;
  std::string to_data_model_string(const World* world, bool print_positive_sign) const;


  CountType type;

  // if sign_in_expression == +1 -> add to the total count
  // if sign_in_expression == -1 -> subtract from the total count
  // 0 - invalid
  int sign_in_expression;

  // TODO: add getters/setters with checks

  // valid when type is EnclosedInWorld, EnclosedInObject or PresentOnSurfaceRegion
  orientation_t orientation;


  SpeciesPatternType species_pattern_type;
  // valid when species_pattern_type is SpeciesId
  species_id_t species_id;
  // valid when species_pattern_type is SpeciesPattern or MoleculesPattern
  BNG::CplxInstance species_molecules_pattern;
  // set in compute_count_species_info based on species_molecules_pattern
  uint_set<species_id_t> species_ids_matching_pattern_cache;

  // valid when type is RxnCountInWorld, RxnCountInObject or RxnOnSurfaceRegion
  BNG::rxn_rule_id_t rxn_rule_id;

  // valid when type is EnclosedInObject or RxnCountInObject
  geometry_object_id_t geometry_object_id;

  // valid when type is PresentOnSurfaceRegion or RxnOnSurfaceRegion
  region_id_t region_id;
};


class MolOrRxnCountItem {
public:
  MolOrRxnCountItem(const count_buffer_id_t buffer_id_)
    : buffer_id(buffer_id_), multiplier(1) {
    assert(buffer_id != COUNT_BUFFER_ID_INVALID);
  }

  void dump(const std::string ind = "") const;
  void to_data_model(const World* world, Json::Value& reaction_output) const;

  // count buffer objects are owned by World
  count_buffer_id_t buffer_id;

  // note: items are shared in MCell3 but so far it seems that
  // we can just count them separately
  std::vector<MolOrRxnCountTerm> terms;

  // value used to multiply the whole result
  float_t multiplier;
};


enum class CountSpeciesInfoType {
  NotSeen,
  Counted,
  NotCounted
};

/**
 * Structure used in caching of information on whether the coutn event counts given species.
 */
struct CountSpeciesInfo {
  CountSpeciesInfo()
    : type(CountSpeciesInfoType::NotSeen),
      multiplicity(UINT_INVALID),
      needs_counted_volume(false) {
  }

  CountSpeciesInfoType type;
  uint multiplicity; // set only when the species are matched with MoleculesPattern
  bool needs_counted_volume;
};

/**
 * Dumps counts of molecules.
 */
class MolOrRxnCountEvent: public BaseEvent {
public:
  MolOrRxnCountEvent(World* world_)
    : BaseEvent(EVENT_TYPE_INDEX_MOL_OR_RXN_COUNT),
      world(world_) {
  }
  virtual ~MolOrRxnCountEvent() {}

  void step() override;

  // DiffuseReactEvent must execute only up to this event
  bool is_barrier() const { return true; }

  void dump(const std::string ind = "") const override;
  void to_data_model(Json::Value& mcell_node) const override;

  void add_mol_count_item(const MolOrRxnCountItem& item) {
    mol_count_items.push_back(item);
  }

  // returns true if this count event count these specific species
  // and whether the species are counted in volume regions so that they
  // need the Molecule's attribute counted_volume to be maintained
  // is cached
  bool species_needs_counted_volume(const species_id_t species_id) {
    return get_or_compute_count_species_info(species_id).needs_counted_volume;
  }

  std::vector<MolOrRxnCountItem> mol_count_items;

  World* world;

private:
  const CountSpeciesInfo& get_or_compute_count_species_info(const species_id_t species_id);
  void compute_count_species_info(const species_id_t species_id);

  // index to this array is species_id
  std::vector<CountSpeciesInfo> count_species_info;
};

} // namespace mcell

#endif // SRC4_MOL_OR_RXN_COUNT_EVENT_H_
