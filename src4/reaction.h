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


#ifndef SRC4_REACTION_H_
#define SRC4_REACTION_H_

#include "defines.h"

namespace MCell {

class SimulationConfig;
class SpeciesInfo;
class Partition;

class SpeciesWithOrientation {
public:
  SpeciesWithOrientation()
    : species_id(SPECIES_ID_INVALID), orientation(ORIENTATION_NONE), equivalent_product_or_reactant_index(INDEX_INVALID) {
  }
  SpeciesWithOrientation(
      const species_id_t species_id_,
      const orientation_t orientation_)
    : species_id(species_id_), orientation(orientation_), equivalent_product_or_reactant_index(INDEX_INVALID) {
  }

  // ignores equivalent_product_or_reactant_index
  bool operator == (const SpeciesWithOrientation& a) const {
    return species_id == a.species_id && orientation == a.orientation;
  }

  bool operator != (const SpeciesWithOrientation& a) const {
    return !(*this == a);
  }

  bool is_same_tolerate_orientation_none(const species_id_t species_id_, const orientation_t orientation_) const {
    return
        species_id == species_id_ &&
        (orientation == ORIENTATION_NONE || orientation_ == ORIENTATION_NONE || orientation == orientation_);
  }

  bool is_on_both_sides_of_rxn() const {
    return equivalent_product_or_reactant_index != INDEX_INVALID;
  }

  species_id_t species_id;
  orientation_t orientation;

  // if valid (not INDEX_INVALID), this specifies a reactant identical to a product,
  // the value is then the index of Reaction's identical product
  uint equivalent_product_or_reactant_index;

  static void dump_array(const std::vector<SpeciesWithOrientation>& vec, const std::string ind);
};


enum class RxnClassType {
  Standard, // any other reaction than below
  Transparent,
  Reflect,
  AbsorbRegionBorder
};

class RxnClass;

// corresponds to MDL reaction and pathway in MCell3 code
// note sure about the naming Rxn or Reaction?
class Rxn {
public:
  /* Name of MDL reaction */
  std::string name;

  /* Rate constant of this pathway */
  float_t rate_constant;

  /* Cumulative probabilities for (entering) this pathway, based on all pathways of the reaction */
  float_t cum_prob;

  // contains the same reactants as the parent Reaction, but possibly in different order
  // checked for correctness on initialize and it is used as constant anyway
  std::vector<SpeciesWithOrientation> reactants;

  std::vector<SpeciesWithOrientation> products;

  uint get_num_surf_products(const SpeciesInfo& all_species) const;

  uint get_num_players() const {
    return reactants.size() + products.size();
  }

  void initialize(const RxnClass& reaction);
  void dump(const std::string ind) const;
  std::string to_string(const Partition& p) const;

private:
  // asserts in debug mode if the reactants are different
  void debug_check_reactants_against_reaction(const RxnClass& reaction);

  // must be called after initialization
  void update_equivalent_product_indices();

  void move_reused_reactants_to_be_the_first_products();
};

// corresponds to reaction in mcell3 code
class RxnClass {
public:
  RxnClassType type;

  /* Maximum 'p' for region of p-space for all non-cooperative pathways */
  float_t max_fixed_p;

  /* Minimum 'p' for region of p-space which is always in the non-reacting "pathway". (note that
     cooperativity may mean that some values of p less than this still do not produce a reaction) */
  float_t min_noreaction_p;

  // reactants are copied into specific reactions as well
  // because a different order might be needed
  std::vector<SpeciesWithOrientation> reactants;

  std::vector<Rxn> reactions;

public:

  uint get_num_reactions() const {
    return reactions.size();
  }

  const Rxn* get_reaction(reaction_index_t rx_index) const {
    assert(rx_index < reactions.size());
    return &reactions[rx_index];
  }

  void add_and_initialize_reaction(const Rxn& r) {
    reactions.push_back(r);
    reactions.back().initialize(*this);
  }

  // there are no pathways for this type of reactions
  bool is_reflect() const {
    return type == RxnClassType::Reflect;
  }

  bool is_transparent() const {
    return type == RxnClassType::Transparent;
  }

  bool is_absorb() const {
    return type == RxnClassType::AbsorbRegionBorder;
  }

  static void dump_array(const std::vector<RxnClass>& vec);

  void dump(const std::string ind) const;

  std::string to_string(const Partition& p) const;
};


/**
 * Used as a pair molecule id, remaining timestep for molecules newly created in diffusion.
 * Using name action instead of event because events are handled by scheduler and are ordered by time.
 * These actions are simply processes in a queue (FIFO) manner.
 *
 * Used in diffuse_react _event_t and in partition_t.
 */
class DiffuseOrUnimolRxnAction {
public:
  enum class Type {
    DIFFUSE,
    UNIMOL_REACT
  };

  // DIFFUSE action
  DiffuseOrUnimolRxnAction(
      const DiffuseOrUnimolRxnAction::Type type_,
      const molecule_id_t id_,
      const float_t scheduled_time_,
      const WallTileIndexPair& where_created_this_iteration_)
    :
      id(id_),
      scheduled_time(scheduled_time_),
      type(type_),
      unimol_rx(nullptr),
      where_created_this_iteration(where_created_this_iteration_) {

    assert(scheduled_time >= 0.0);
    assert(type == Type::DIFFUSE);
    // position where the molecule was created may be invalid when it was not a result of surface reaction
  }

  // UNIMOL_REACT action
  DiffuseOrUnimolRxnAction(
      const DiffuseOrUnimolRxnAction::Type type_,
      const molecule_id_t id_,
      const float_t scheduled_time_,
      const RxnClass* unimol_rx_)
    :
      id(id_),
      scheduled_time(scheduled_time_),
      type(type_),
      unimol_rx(unimol_rx_) {
    assert(scheduled_time >= 0.0);
    assert(type == Type::UNIMOL_REACT);
    assert(unimol_rx != nullptr);
  }

  // defined because of usage in calendar_t
  const DiffuseOrUnimolRxnAction& operator->() const {
     return *this;
  }

  molecule_id_t id;
  float_t scheduled_time; // this is the scheduled time
  Type type;

  // when type is UNIMOL_REACT
  const RxnClass* unimol_rx;

  // when type is DIFFUSE
  // used to avoid rebinding for surf+vol->surf+vol reactions
  WallTileIndexPair where_created_this_iteration;
};

} // namespace mcell

#endif // SRC4_REACTION_H_
