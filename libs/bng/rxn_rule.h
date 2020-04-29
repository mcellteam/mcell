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

#include "bng/cplx_instance.h"

namespace BNG {

class BNGData;
class SpeciesContainer;

struct CplxMolIndex {
  CplxMolIndex()
  : cplx_index(INDEX_INVALID), mol_index(INDEX_INVALID) {
  }

  CplxMolIndex(const uint complex_index_, const uint molecule_index_)
    : cplx_index(complex_index_), mol_index(molecule_index_) {
  }

  uint cplx_index;
  uint mol_index;

  bool operator==(const CplxMolIndex& cmi2) const {
    return cplx_index == cmi2.cplx_index && mol_index == cmi2.mol_index;
  }
};


struct CMIndexPair {
  CMIndexPair(const CplxMolIndex& reactant_cmi_, const CplxMolIndex& product_cmi_)
    : reactant_cmi(reactant_cmi_), product_cmi(product_cmi_) {
  }

  CplxMolIndex reactant_cmi;
  CplxMolIndex product_cmi;

  bool operator==(const CMIndexPair& cmi_pair2) const {
    return reactant_cmi == cmi_pair2.reactant_cmi && product_cmi == cmi_pair2.product_cmi;
  }
};


struct CplxIndexPair {
  CplxIndexPair(const uint reactant_index_, const uint product_index_)
    : reactant_index(reactant_index_), product_index(product_index_) {
  }

  uint reactant_index;
  uint product_index;
};


enum class RxnType {
  Invalid,
  Standard, // any other reaction than below
  Transparent,
  Reflect,
  AbsorbRegionBorder
};


// BNG reaction rule
// rules are only unidirectional,
// if there is a reversible reaction in BNGL definition,
// two
class RxnRule: public BaseFlag {
public:
  std::string name;
  rxn_rule_id_t id;

  RxnType type;

  // the complex species are patterns
  //
  // there is a potential for optimizations, e.g.
  // to make a set of species that match the patterns, but let's keep it
  // for later
  CplxInstanceVector reactants;
  CplxInstanceVector products;

  // set to true if it was possible to do a mapping between reactants and products
  bool mol_instances_are_fully_maintained;

  // matching between molecules of reactants and molecules of products,
  // contains info on what we were able to match, even if
  // mol_instances_are_fully_maintained is false
  small_vector<CMIndexPair> mol_mapping;

  // matching between complexes
  // set only if the complex patterns are identical (MCell-style rxns)
  small_vector<CplxIndexPair> cplx_mapping;

  // caching
  uint_set<species_id_t> species_applicable_as_reactants;
  uint_set<species_id_t> species_not_applicable_as_reactants;

  float_t rate_constant;

private:
  bool finalized;
  uint num_surf_products;

public:
  RxnRule()
    : id(RXN_RULE_ID_INVALID), type(RxnType::Invalid), mol_instances_are_fully_maintained(false), rate_constant(FLT_INVALID),
     finalized(false), num_surf_products(UINT_INVALID) {
  }

  void finalize();

  bool is_finalized() const {
    return finalized;
  }

  const CplxInstance& get_cplx_reactant(const uint index) const {
    assert(index <= reactants.size());
    return reactants[index];
  }

  const CplxInstance& get_cplx_product(const uint index) const {
    assert(index <= products.size());
    return products[index];
  }

  const MolInstance& get_mol_reactant(const CplxMolIndex& cmi) const {
    assert(cmi.mol_index <= reactants[cmi.cplx_index].mol_instances.size());
    return get_cplx_reactant(cmi.cplx_index).mol_instances[cmi.mol_index];
  }

  const MolInstance& get_mol_product(const CplxMolIndex& cmi) const {
    assert(cmi.mol_index <= products[cmi.cplx_index].mol_instances.size());
    return get_cplx_product(cmi.cplx_index).mol_instances[cmi.mol_index];
  }

  // mcell3 variant of maintaining substances,
  // e.g. for A + B -> A : reactant A is maintained
  bool is_cplx_reactant_on_both_sides_of_rxn(const uint index) const;
  bool is_cplx_product_on_both_sides_of_rxn(const uint index) const;

  bool operator ==(const RxnRule& rr2) {
    // ordering of components in a molecule is important
    // two component types must have the same id, this is ensured in find_or_add_component_type
    return
        name == rr2.name &&
        reactants == rr2.reactants && products == rr2.products &&
        rate_constant == rr2.rate_constant;
  }

  void dump(const BNGData& bng_data, const std::string ind = "") const;

  // checks if it is possible to create a mapping from reactants to products and
  // sets members molecule_instances_are_maintained and mapping,
  // might write some error messages to the msgs stream,
  // returns true if errors were encountered
  bool compute_reactants_products_mapping();

  bool compute_reactants_products_mapping_w_error_output(const BNGData& bng_data, std::ostream& out);


  void append_reactant(const CplxInstance& inst) {
    reactants.push_back(inst);
  }

  void append_product(const CplxInstance& inst) {
    products.push_back(inst);
  }

  uint get_num_surf_products() const { // we don't have probably the information that is needed
    assert(finalized);
    return num_surf_products;
  }

  uint get_num_players() const {
    return reactants.size() + products.size();
  }

  bool is_unimol() const {
    return reactants.size() == 1;
  }

  bool is_bimol() const {
    return reactants.size() == 2;
  }

  // returns true if species 'id' matches one of the reactants
  bool species_can_be_reactant(const species_id_t id, const SpeciesContainer& all_species);

  // returns true if two reactants match each other and species 'id' matches one of the reactants
  bool species_is_both_bimol_reactants(const species_id_t id, const SpeciesContainer& all_species);

  bool find_assigned_cplx_reactant_for_product(const uint product_index, uint& reactant_index) const;

  void set_is_counted() {
    set_flag(RXN_FLAG_COUNTED);
  }

  bool is_counted() const {
    return has_flag(RXN_FLAG_COUNTED);
  }

private:

  // returns false if cmi was not found in mapping,
  bool find_assigned_mol_reactant_for_product(const CplxMolIndex& product_cmi, CplxMolIndex& reactant_cmi) const;

  // check if it makes sense to compute mapping at all
  bool has_same_mols_in_reactants_and_products() const;

  // returns false if no fitting product was found
  bool find_most_fitting_unassigned_mol_product(const CplxMolIndex& reactant_cmi, CplxMolIndex& best_product_cmi) const;

  bool compute_mol_reactants_products_mapping(MolInstance& not_matching_mol_inst, CplxMolIndex& not_matching_cmi);

  void compute_cplx_reactants_products_mapping();

  void move_products_that_are_also_reactants_to_be_the_first_products();

  void dump_complex_instance_vector(const BNGData& bng_data, const CplxInstanceVector& complexes) const;

};

} /* namespace BNG */

#endif /* LIBS_BNG_RXN_RULE_H_ */
