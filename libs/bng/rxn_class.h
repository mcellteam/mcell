/*
 * RxnClass.h
 *
 *  Created on: Mar 25, 2020
 *      Author: ahusar
 */

#ifndef LIBS_BNG_RXN_CLASS_H_
#define LIBS_BNG_RXN_CLASS_H_

#include <string>
#include <iostream>

#include "bng/bng_defines.h"
#include "bng/rxn_rule.h"

namespace BNG {

class RxnContainer;
class SpeciesContainer;


/**
 * Reaction class contains all applicable reactions for a pair of
 * reactants with their compartment.
 *
 * Reaction classes are created on-the-fly by RxnContainer.
 */
class RxnClass {
public:
  // Standard reaction or special such as Reflect, Transparent or Absorb
  RxnType type;

  // keeping just IDs, with IDs unlike with pointers we are able to check that the species was 'discarded'
  // the rxn class was created specifically for a single or a pair of reactants, therefore
  // the species id is is a specific instance of a complex that matches our pattern
  // and if there are multiple species that match, multiple rxn classes are created
  std::vector<species_id_t> reactant_ids;

  // - information on products for specific reactions, based on all reactions of the class,
  // - the size of this vector is >= than the number of rxn rules in this class because
  //   each rxn rule can have multiple product sets
  // - contains cummulative probabilities
  // - WARNING: do not access directly, only through get_rxn_products_for_pathway (TODO: make friend with RxnRule?),
  //     because the product definitions may not be initialized
  RxnClassPathwayVector pathways;

private:
  // reactions are owned by RxnContainer, order in this vector is important
  std::vector<rxn_rule_id_t> rxn_rule_ids;

  // Maximum probability for region of p-space for all non-cooperative pathways,
  // valid when pathways_and_rates_initialized is true
  float_t max_fixed_p;

  // owned by BNGEngine
  RxnContainer& all_rxns;
  SpeciesContainer& all_species;

  // owned by the simulation engine
  const BNGConfig& bng_config;

  // flags for optimized testing of this rxn class
  bool bimol_vol_rxn_flag;
  bool intermembrane_surf_surf_rxn_flag;

  // flag for initialization of pathways on-demand
  bool pathways_and_rates_initialized;

public:
  RxnClass(
      RxnContainer& all_rxns_, SpeciesContainer& all_species_, const BNGConfig& bng_config_,
      const species_id_t reactant1_id, const species_id_t reactant2_id = SPECIES_ID_INVALID)
    : type(RxnType::Invalid), max_fixed_p(FLT_INVALID),
      all_rxns(all_rxns_), all_species(all_species_), bng_config(bng_config_),
      bimol_vol_rxn_flag(false), intermembrane_surf_surf_rxn_flag(false),
      pathways_and_rates_initialized(false)
    {
    reactant_ids.push_back(reactant1_id);
    if (reactant2_id != SPECIES_ID_INVALID) {
      reactant_ids.push_back(reactant2_id);
    }
  }

  // RxnRules may have pointers to rxn classes that use them, we must
  // remove these links in the destructor
  ~RxnClass();

  uint get_num_reactions() const {
    return rxn_rule_ids.size();
  }

  uint get_num_pathways() const {
    return pathways.size();
  }

  // first query for the max probability or for rxn rate update usually causes
  // rxn pathways to be initialized
  float_t get_max_fixed_p() {
    if (!pathways_and_rates_initialized) {
      init_rxn_pathways_and_rates();
    }
    return max_fixed_p;
  }

  void update_rxn_rates_if_needed(const float_t current_time);

  RxnRule* get_rxn_for_pathway(const rxn_class_pathway_index_t pathway_index);

  const RxnProductsVector& get_rxn_products_for_pathway(const rxn_class_pathway_index_t pathway_index) {
    release_assert(pathways_and_rates_initialized && "Must have been initialized when specific pathways is queried");

    assert((size_t)pathway_index < pathways.size());

    if (pathways[pathway_index].products_are_defined) {
      return pathways[pathway_index].product_species_w_indices;
    }
    else {
      define_rxn_pathway_using_mapping(pathway_index);
      return pathways[pathway_index].product_species_w_indices;
    }
  }

  rxn_class_pathway_index_t get_pathway_index_for_probability(
      const float_t prob, const float_t local_prob_factor);

  orientation_t get_reactant_orientation(uint reactant_index) const;

  // does not do pathways update
  void add_rxn_rule_no_update(RxnRule* r);

  // this function expects that update_rxn_rates_if_needed was called
  // already for the current time
  float_t get_next_time_of_rxn_rate_update() const;

  bool is_standard() const {
    return type == RxnType::Standard;
  }

  // there is exactly one reaction for this type
  bool is_reflect() const {
    return type == RxnType::Reflect;
  }

  bool is_transparent() const {
    return type == RxnType::Transparent;
  }

  bool is_unimol() const {
    return reactant_ids.size() == 1;
  }

  bool is_bimol() const {
    return reactant_ids.size() == 2;
  }

  bool is_bimol_vol_rxn_class() const {
    #ifndef NDEBUG
      debug_check_bimol_vol_rxn_flag();
    #endif
    return bimol_vol_rxn_flag;
  }

  bool is_intermembrane_surf_surf_rxn_class() const {
    return intermembrane_surf_surf_rxn_flag;
  }

  bool is_absorb_region_border() const {
    return type == RxnType::AbsorbRegionBorder;
  }

  bool is_reflect_transparent_or_absorb_region_border() const {
    return
        is_reflect() ||
        is_transparent() ||
        is_absorb_region_border();
  }

  bool is_simple() const;

  bool is_reactant_species_id(species_id_t reactant_id) const {
    if (is_unimol()) {
      return reactant_id == reactant_ids[0];
    }
    else {
      return reactant_id == reactant_ids[0] || reactant_id == reactant_ids[1];
    }
  }

  species_id_t get_second_species_id(species_id_t reactant_id) const {
    assert(is_bimol());
    assert(is_reactant_species_id(reactant_id));
    if (reactant_ids[0] != reactant_id) {
      return reactant_ids[0];
    }
    else {
      return reactant_ids[1];
    }
  }

  // asserts if there is no such reactant
  species_id_t get_reactive_surface_reactant_species_id() const;

  // initializes pathways, called automatically once rates of
  // products of this rxn class are needed, called on-demand
  // because computing even initial pathways without specific product species
  // for complex (and long) reactants may be costly)
  // used also when rates are updated
  void init_rxn_pathways_and_rates(const bool force_update = false);

  void remove_from_rxn_rule_users_list();

  std::string to_str(const std::string ind = "") const;

  static void dump_array(const std::vector<RxnClass>& vec);
  void dump(const std::string ind = "") const;

private:

  void define_rxn_pathway_using_mapping(const rxn_class_pathway_index_t pathway_index);

  void update_variable_rxn_rates(const float_t current_time);

  void debug_check_bimol_vol_rxn_flag() const;

  std::string reactants_to_str() const;

  float_t get_reactant_diffusion(const uint reactant_index) const;
  float_t get_reactant_space_step(const uint reactant_index) const;
  float_t get_reactant_time_step(const uint reactant_index) const;

  float_t compute_pb_factor() const;
};


typedef small_vector<RxnClass*> RxnClassesVector;
typedef std::set<RxnClass*> RxnClassesSet;

} /* namespace BNG */

#endif /* LIBS_BNG_RXN_CLASS_H_ */
