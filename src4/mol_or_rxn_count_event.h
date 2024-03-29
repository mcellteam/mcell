/******************************************************************************
 *
 * Copyright (C) 2019 by
 * The Salk Institute for Biological Studies and
 * Pittsburgh Supercomputing Center, Carnegie Mellon University
 *
 * Use of this source code is governed by an MIT-style
 * license that can be found in the LICENSE file or at
 * https://opensource.org/licenses/MIT.
 *
******************************************************************************/

#ifndef SRC4_MOL_OR_RXN_COUNT_EVENT_H_
#define SRC4_MOL_OR_RXN_COUNT_EVENT_H_

#include "bng/bng.h"

#include "base_event.h"
#include "count_buffer.h"
#include "region_expr.h"

namespace BNG {
class SpeciesContainer;
}

namespace MCell {

class Partition;
class Molecule;
class World;

enum class SpeciesPatternType {
  Invalid,
  SpeciesId, // TODO: remove this variant
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
      primary_compartment_id(BNG::COMPARTMENT_ID_NONE),
      rxn_rule_id(BNG::RXN_RULE_ID_INVALID),
      initial_reactions_count(0)
     {
  }

  uint get_num_pattern_matches(const species_id_t species_id) const;

  uint get_num_molecule_matches(
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

  // valid when type is EnclosedInWorld, EnclosedInObject or PresentOnSurfaceRegion
  orientation_t orientation;

  SpeciesPatternType species_pattern_type;

  // valid when species_pattern_type is SpeciesId, not used in MCell4
  // TODO: remove for MCell4 MDL mode and replace with SpeciesPattern/species_molecules_pattern
  species_id_t species_id;

  // valid when species_pattern_type is SpeciesPattern or MoleculesPattern
  BNG::Cplx species_molecules_pattern;

  // when primary_compartment_id is COMPARTMENT_ID_NONE, it is ignored
  // used only for molecule counting
  BNG::compartment_id_t primary_compartment_id;

  // set in compute_count_species_info based on species_molecules_pattern
  // presence is tested when species_pattern_type is SpeciesPattern
  // the value is used when species_pattern_type is MoleculesPattern
  std::map<species_id_t, uint> species_ids_matching_pattern_w_multiplier_cache;

  // valid when type is RxnCountInWorld, RxnCountInObject or RxnOnSurfaceRegion
  BNG::rxn_rule_id_t rxn_rule_id;

  // valid when type is EnclosedInObject or RxnCountInObject or
  // PresentOnSurfaceRegion or RxnOnSurfaceRegion
  // region_expr.root is nullptr otherwise
  RegionExpr region_expr;

  // holds initial value when resumed from a checkpoint,
  // ignored for molecule counts
  uint64_t initial_reactions_count;
};

typedef small_vector<MolOrRxnCountTerm> MolOrRxnCountTermVector;


class MolOrRxnCountItem {
public:
  MolOrRxnCountItem(
      const count_buffer_id_t buffer_id_, const uint buffer_column_index_)
    : index(UINT_INVALID), buffer_id(buffer_id_), buffer_column_index(buffer_column_index_), multiplier(1) {
    assert(buffer_id != COUNT_BUFFER_ID_INVALID);
    assert(buffer_column_index != UINT_INVALID);
  }

  bool is_world_mol_count() const;

  bool counts_mols() const;
  bool counts_rxns() const;

  void dump(const std::string ind = "") const;
  void to_data_model(const World* world, Json::Value& reaction_output) const;

  // index of this item in MolOrRxnCountEvent's mol_rxn_count_items
  uint index;

  // count buffer objects are owned by World
  count_buffer_id_t buffer_id;

  // index of column in buffer where this value will be stored
  // always 0 when the output format is DAT
  uint buffer_column_index;

  // note: items are shared in MCell3 but so far it seems that
  // we can just count them separately
  MolOrRxnCountTermVector terms;

  // value used to multiply the whole result
  double multiplier;
};

typedef small_vector<MolOrRxnCountItem> MolOrRxnCountItemVector;


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
      all_are_world_mol_counts(true) {
  }

  CountSpeciesInfoType type;

  // when true, all count items that count this species are listed in
  // the world_count_item_indices
  bool all_are_world_mol_counts;

  // indices of count items that count this species in the whole world
  uint_set<uint> world_count_item_indices;
};

/**
 * Dumps counts of molecules.
 */
class MolOrRxnCountEvent: public BaseEvent {
public:
  MolOrRxnCountEvent(World* world_)
    : BaseEvent(EVENT_TYPE_INDEX_MOL_OR_RXN_COUNT),
      world(world_), count_mols(false), count_rxns(false) {
  }
  virtual ~MolOrRxnCountEvent() {
  }

  void step() override;

  // used from MCell4 API
  double get_single_count_value();

  // DiffuseReactEvent must execute only up to this event
  bool is_barrier() const override { return true; }

  void dump(const std::string ind = "") const override;
  void to_data_model(Json::Value& mcell_node) const override;

  void add_mol_count_item(const MolOrRxnCountItem& item) {
    mol_rxn_count_items.push_back(item);
    // set index
    mol_rxn_count_items.back().index = mol_rxn_count_items.size() - 1;

    // set counting flags
    if (item.counts_mols()) {
      count_mols = true;
    }
    if (item.counts_rxns()) {
      count_rxns = true;
    }
  }

  MolOrRxnCountItemVector mol_rxn_count_items;

  World* world;

private:
  const CountSpeciesInfo& get_or_compute_count_species_info(const species_id_t species_id);
  void compute_count_species_info(const species_id_t species_id);

  void compute_mol_count_item(
      const Partition& p,
      const MolOrRxnCountItem& item,
      const Molecule& m,
      CountItemVector& count_items
  );

  void compute_rxn_count_item(
      Partition& p,
      const MolOrRxnCountItem& item,
      const BNG::RxnRule* rxn,
      CountItemVector& count_items
  );

  void compute_counts(CountItemVector& count_items);

  // index to this array is species_id
  std::vector<CountSpeciesInfo> count_species_info;

  // flags to optimize counting
  bool count_mols;
  bool count_rxns;
};

} // namespace mcell

#endif // SRC4_MOL_OR_RXN_COUNT_EVENT_H_
