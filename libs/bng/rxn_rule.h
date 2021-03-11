/*
 * rule.h
 *
 *  Created on: Jan 9, 2020
 *      Author: ahusar
 */

#ifndef LIBS_BNG_RXN_RULE_H_
#define LIBS_BNG_RXN_RULE_H_

#include <string>
#include <iostream>

#include "bng/bng_defines.h"
#include "bng/cplx.h"

namespace BNG {

class BNGData;
class SpeciesContainer;
class RxnClass;

/**
 * Used to hold information on a single product and
 * what reaction rule product indices were used to create it.
 * Usually a single product corresponds to one product on the right-hand side of the
 * reaction rule, but in cases such as when a single bond is broken and the complex
 * remains as a whole, one resulting product corresponds to multiple
 * products of the rxn rule.
 */
class ProductSpeciesIdWIndices {
public:
  ProductSpeciesIdWIndices()
    : product_species_id(SPECIES_ID_INVALID) {
  }

  // usual case - one product per rxn rule product
  ProductSpeciesIdWIndices(
      const species_id_t product_species_id_,
      const uint product_index)
    : product_species_id(product_species_id_) {
    rule_product_indices.insert(product_index);
  }

  // general case
  ProductSpeciesIdWIndices(
      const species_id_t product_species_id_,
      const std::set<uint>& rule_product_indices_)
    : product_species_id(product_species_id_),
      rule_product_indices(rule_product_indices_) {
  }

  species_id_t product_species_id;
  // must use container with guaranteed order
  std::set<uint> rule_product_indices;
};

typedef std::vector<ProductSpeciesIdWIndices> RxnProductsVector;

/**
 * Similar as ProductSpeciesWIndices, only uses cplx inst instead of species id.
 */
class ProductSpeciesPtrWIndices {
public:
  ProductSpeciesPtrWIndices(
      Species* product_species_,
      const std::set<uint>& rule_product_indices_)
    : product_species(product_species_),
      rule_product_indices(rule_product_indices_) {
  }

  // does not own this object
  Species* product_species;

  // must use container with guaranteed order
  std::set<uint> rule_product_indices;
};

typedef std::vector<ProductSpeciesPtrWIndices> ProductCplxWIndicesVector;

// - first dimension are individual products that a rxn rule can produce
// - second dimension are individual products
typedef std::vector<std::vector<ProductSpeciesPtrWIndices>> ProductSetsVector;


// TODO: some of these classes rather belong to rxn_class.h
/**
 * Used to hold information on the probability and products for a given reaction rule.
 * One rxn rule can possibly have more different products,
 * e.g. when applying rule A(b~0) -> A(b~1) on reactant A(a~0,b~0).A(a~1,b~0),
 * the result can be either A(a~0,b~1).A(a~1,b~0) or A(a~0,b~0).A(a~1,b~1).
 * So for each of the possible products, one RxnPathway is created in a RxnClass.
 */
class RxnClassPathway {
public:
  RxnClassPathway(
      const rxn_rule_id_t rxn_rule_id_,
      const float_t pathway_prob_,
      const RxnProductsVector& product_species_w_indices_
  ) : rxn_rule_id(rxn_rule_id_),
      pathway_prob(pathway_prob_),
      product_species_w_indices(product_species_w_indices_),
      cum_prob(FLT_INVALID),
      products_are_defined(true) {
  }

  RxnClassPathway(
      const rxn_rule_id_t rxn_rule_id_,
      const float_t pathway_prob_,
      const VertexMapping& rule_mapping_onto_reactants_
  ) : rxn_rule_id(rxn_rule_id_),
      pathway_prob(pathway_prob_),
      rule_mapping_onto_reactants(rule_mapping_onto_reactants_),
      cum_prob(FLT_INVALID),
      products_are_defined(false) {
  }

  // ID of rxn rule from which this pathway was created
  rxn_rule_id_t rxn_rule_id;

  // probability for this specific pathway to be selected when
  // reactants interact (or when an unimol rxn is executed)
  float_t pathway_prob;

  // specific variant of products for this pathway,
  // valid when products_are_defined is true
  RxnProductsVector product_species_w_indices;

  // we did not compute products, only assigned which pattern to reactant
  // mapping is associated with this pathway, needs to be computed once needed
  // valid when products_are_defined is false
  VertexMapping rule_mapping_onto_reactants;

  // cumulative probability - set when used in rxn class
  float_t cum_prob;

  bool products_are_defined;
};

typedef std::vector<RxnClassPathway> RxnClassPathwayVector;


struct CplxIndexPair {
  CplxIndexPair(const uint reactant_index_, const uint product_index_, const bool is_simple_mapping_)
    : reactant_index(reactant_index_), product_index(product_index_), is_simple_mapping(is_simple_mapping_) {
  }

  uint reactant_index;
  uint product_index;
  bool is_simple_mapping; // set to true when this is a mapping between simple complexes
};


struct RxnRateInfo {
  float_t time;
  float_t rate_constant;

  bool operator < (const RxnRateInfo& ri2) const {
    return time < ri2.time;
  }
};

enum class RxnType {
  Invalid,
  // Standard is any other reaction than below,
  // standard rxns are also created from absorptive surfaces and also concentration clamps
  Standard,
  Transparent,
  Reflect,
  AbsorbRegionBorder
};


// BNG reaction rule
// rules are only unidirectional,
// if there is a reversible reaction in BNGL definition,
// two RxnRules are created
class RxnRule: public BaseFlag {
public:
  RxnRule(const BNGData* bng_data_)
    : id(RXN_RULE_ID_INVALID), type(RxnType::Invalid),
      base_rate_constant(FLT_INVALID),
      mol_instances_are_fully_maintained(false),
      next_variable_rate_index(0),
      bng_data(bng_data_)
      {
  }

  // after finalize one should must
  // check_compartments_and_set_orientations (from rxn_compartment_utils)
  // to check that compartments are valid with respect to volume and surface molecules
  // and to set orientations from compartments
  void finalize();

  // NOTE: must be called only after molecule types are fully known, i.e.
  // whether they are surface or volume reactants otherwise MCell
  // orientation won't be printed
  void set_automatic_name(const bool reverse) {
    name = ((reverse)?"rev ":"") + to_str(false, false, false);
  }

  std::string name;
  rxn_rule_id_t id;

  RxnType type;

  // the complex species of reactants are patterns
  CplxVector reactants;
  CplxVector products;

  // base rate constant for this reaction as obtained from the BNGL or other input
  float_t base_rate_constant;

  small_vector<RxnRateInfo> base_variable_rates;

  // set to true if it was possible to do a mapping between reactants and products
  bool mol_instances_are_fully_maintained;

  // caching
  // TODO: we are keeping here all the species from the past, might use some cleanup with species cleanup as well
  uint_set<species_id_t> species_applicable_as_any_reactant; // same as species_applicable_as_reactant[0] for unimol rxns
  uint_set<species_id_t> species_applicable_as_reactant[2]; // [1] is empty for unimol rxns

  uint_set<species_id_t> species_not_applicable_as_any_reactant; // same as species_applicable_as_reactant1 for unimol rxns
  uint_set<species_id_t> species_not_applicable_as_reactant[2]; // [1] is empty for unimol rxns

private:
  // variable reaction rate constants, sorted by time
  // index is initialized to 0
  uint next_variable_rate_index;

  // maintain information on where this reaction was used in order to
  // update all classes if this reaction's rate constant changes due to variable reaction
  // rate MDL feature, not used in
  // WARNING: do not use to figure out whether which species may be affected because this
  // gets periodically cleaned
  std::set<RxnClass*> rxn_classes_where_used;

public:
  // - method used when RxnClass is being created
  // - defines new species when needed and might invalidate Species references
  // TODO: pass reactants as vector
  void define_rxn_pathways_for_specific_reactants(
      SpeciesContainer& all_species,
      const BNGConfig& bng_config,
      const species_id_t reactant_a_species_id,
      const species_id_t reactant_b_species_id,
      const float_t pb_factor,
      RxnClassPathwayVector& pathways
  );


  void define_rxn_pathway_using_mapping(
    SpeciesContainer& all_species,
    const BNGConfig& bng_config,
    const std::vector<species_id_t>& reactant_species,
    RxnClassPathway& pathway
  ) const;

  float_t get_rate_constant() const {
    return base_rate_constant;
  }

  const Cplx& get_cplx_reactant(const uint index) const {
    assert(index <= reactants.size());
    return reactants[index];
  }

  const Cplx& get_cplx_product(const uint index) const {
    assert(index <= products.size());
    return products[index];
  }

  // mcell3 variant of maintaining substances,
  // e.g. for A + B -> A : reactant A is maintained
  bool is_simple_cplx_reactant_on_both_sides_of_rxn_w_identical_compartments(const uint index) const;
  bool is_simple_cplx_product_on_both_sides_of_rxn_w_identical_compartments(const uint index) const;

  bool operator ==(const RxnRule& rr2) {
    // ordering of components in a molecule is not important
    // two component types must have the same id, this is ensured in find_or_add_component_type
    return
        name == rr2.name &&
        reactants == rr2.reactants && products == rr2.products &&
        base_rate_constant == rr2.base_rate_constant;
  }

  // used in semantic check
  bool check_reactants_products_mapping(std::ostream& out);

  void append_reactant(const Cplx& pattern) {
    assert(pattern.get_primary_compartment_id() != COMPARTMENT_ID_INVALID);
    reactants.push_back(pattern);
  }

  void append_product(const Cplx& cplx) {
    assert(cplx.get_primary_compartment_id() != COMPARTMENT_ID_INVALID);
    products.push_back(cplx);
  }

  bool is_unimol() const {
    return reactants.size() == 1;
  }

  bool is_bimol() const {
    return reactants.size() == 2;
  }

  bool is_absorptive_region_rxn() const {
    return is_bimol() && reactants[1].is_reactive_surface() && products.empty();
  }

  bool is_bimol_vol_rxn() const {
    if (is_unimol()) {
      return false;
    }
    else if (is_bimol()) {
      return reactants[0].is_vol() && reactants[1].is_vol();
    }
    else {
      assert(false);
      return false;
    }
  }

  bool is_vol_rxn() const {
    if (is_unimol()) {
      return reactants[0].is_vol();
    }
    else if (is_bimol()) {
      return reactants[0].is_vol() && reactants[1].is_vol();
    }
    else {
      assert(false);
      return false;
    }
  }

  bool is_surf_rxn() const {
    if (is_unimol()) {
      return reactants[0].is_surf();
    }
    else if (is_bimol()) {
      return reactants[0].is_surf() || reactants[1].is_surf();
    }
    else {
      assert(false);
      return false;
    }
  }

  bool is_reactive_surface_rxn() const {
    if (is_unimol()) {
      return false;
    }
    else if (is_bimol()) {
      return reactants[0].is_reactive_surface() || reactants[1].is_reactive_surface();
    }
    else {
      assert(false);
      return false;
    }
  }

  void set_is_intermembrane_surf_rxn() {
    set_flag(RXN_FLAG_INTERMEMBRANE);
  }

  bool is_intermembrane_surf_rxn() const {
    return has_flag(RXN_FLAG_INTERMEMBRANE);
  }

  bool reactants_use_compartments() const {
    if (is_unimol()) {
      return reactants[0].has_compartment();
    }
    else if (is_bimol()) {
      return reactants[0].has_compartment() || reactants[1].has_compartment();
    }
    else {
      assert(false);
      return false;
    }
  }

  void get_bimol_reactant_indices(
      const species_id_t id1, const species_id_t id2,
      const SpeciesContainer& all_species,
      std::vector<std::pair<uint, uint>>& reac_indices);

  void get_reactant_indices(
      const species_id_t species_id, const SpeciesContainer& all_species,
      std::vector<uint>& indices);

  // returns true if species 'id' matches one of the reactants
  // updates local cache
  bool species_can_be_reactant(const species_id_t id, const SpeciesContainer& all_species);

  // returns true if both species can be used as separate reactants for a bimol rxn
  // sets assigned indices if assigned_indexN are not nullptr
  bool species_can_be_bimol_reactants(
      const species_id_t id1, const species_id_t id2, const SpeciesContainer& all_species,
      uint* assigned_index1 = nullptr, uint* assigned_index2 = nullptr, bool* both_match_both_patterns = nullptr
  );

  // returns true if two reactants match each other and species 'id' matches one of the reactants
  bool species_is_both_bimol_reactants(const species_id_t id, const SpeciesContainer& all_species);

  bool get_assigned_cplx_reactant_for_product(
      const uint product_index, const bool only_simple, uint& reactant_index) const;

  void set_is_counted_in_world() {
    set_flag(RXN_FLAG_COUNTED_IN_WORLD);
  }

  void set_is_counted_in_volume_regions() {
    set_flag(RXN_FLAG_COUNTED_IN_VOLUME_REGIONS);
  }

  void set_is_counted_on_surface_regions() {
    set_flag(RXN_FLAG_COUNTED_ON_SURFACE_REGIONS);
  }

  bool is_counted() const {
    return
        has_flag(RXN_FLAG_COUNTED_IN_WORLD) ||
        has_flag(RXN_FLAG_COUNTED_IN_VOLUME_REGIONS) ||
        has_flag(RXN_FLAG_COUNTED_ON_SURFACE_REGIONS);
  }

  bool is_counted_in_volume_regions() const {
    return has_flag(RXN_FLAG_COUNTED_IN_VOLUME_REGIONS);
  }

  bool is_counted_on_surface_regions() const {
    return has_flag(RXN_FLAG_COUNTED_ON_SURFACE_REGIONS);
  }

  bool is_simple() const {
    return has_flag(RXN_FLAG_SIMPLE);
  }

  void add_rxn_class_where_used(RxnClass* rxn_class) {
    rxn_classes_where_used.insert(rxn_class);
  }

  void reset_rxn_classes_where_used() {
    rxn_classes_where_used.clear();
  }

  void remove_rxn_class_where_used(RxnClass* rxn_class) {
    assert(rxn_classes_where_used.count(rxn_class) != 0);
    rxn_classes_where_used.erase(rxn_class);
  }

  void remove_species_id_references(const species_id_t id);

  // returns false when there are no variable rates or we already processed all scheduled times
  bool may_update_rxn_rate() const {
    return next_variable_rate_index < base_variable_rates.size();
  }

  // returns true if rate was updated
  // requester is the rxn class that requested this update
  bool update_variable_rxn_rate(const float_t current_time, const RxnClass* requester);

  float_t get_next_time_of_rxn_rate_update() const {
    if (may_update_rxn_rate()) {
      return base_variable_rates[next_variable_rate_index].time;
    }
    else {
      return TIME_FOREVER;
    }
  }

  // manual update of rxn rate, does nothing if the new rate is the same as existing
  // returns true if the rate was updated
  bool update_rxn_rate(const float_t new_rate);

  std::string to_str(
      const bool with_rate_constant = true, const bool with_name = true, const bool with_id = true) const;
  std::string reactants_to_str() const;
  std::string products_to_str() const;
  void dump(const bool for_diff = false, const std::string ind = "", std::ostream& out = std::cout) const;

private:

  void set_product_compartments(Species* product_species, std::set<uint>& product_indices) const;
  void create_products_from_reactants_graph(
      const BNGData* bng_data,
      Graph& reactants_graph,
      ProductCplxWIndicesVector& created_products
  ) const;

  // BNGL style reaction handling is implemented in this method
  void create_products_for_complex_rxn(
      SpeciesContainer& all_species,
      const BNGConfig& bng_config,
      const std::vector<species_id_t>& reactant_species,
      const float_t pb_factor,
      RxnClassPathwayVector& pathways
  );

  // use other methods that provide caching instead
  void get_reactant_indices_uncached(
        const Cplx& cplx, const SpeciesContainer& all_species,
        std::vector<uint>& indices) const;

  void create_patterns_graph();
  void create_products_graph();

  void move_products_that_are_also_reactants_to_be_the_first_products();

  // checks if it is possible to create a mapping from reactants to products and
  // sets members molecule_instances_are_maintained and mapping
  void compute_reactants_products_mapping();

  bool check_components_mapping(
      const ElemMol& first_mi,
      const ElemMol& second_mi,
      const char* msg,
      std::ostream& out
  );

  bool check_components_states(
      const ElemMol& prod_mi,
      const ElemMol& pat_mi,
      std::ostream& out
  );

  bool may_modify_more_than_one_identical_component() const;
  bool matching_may_produce_multiple_identical_results() const;

  std::string cplx_vector_to_str(const CplxVector& complexes) const;
  void dump_cplx_vector(
      const CplxVector& complexes,
      const std::string ind) const;


  // mutable is needed because of usage in create_products_for_complex_rxn
  // the graphs are not modified, but boost cannot use them as const
  mutable Graph patterns_graph; // graphs based on reactants
  mutable Graph products_graph;
  VertexMapping products_to_patterns_mapping;

  // maps complexes from their pattern to the product, they must use the same compartments
  // used mainly for handling of MCell3 reaction behavior where reactants that stay the same on the 
  // products side are kept unchanged in a reaction
  small_vector<CplxIndexPair> pat_prod_cplx_mapping;

  const BNGData* bng_data; // needed to create results of complex reactions
};

} /* namespace BNG */

#endif /* LIBS_BNG_RXN_RULE_H_ */
