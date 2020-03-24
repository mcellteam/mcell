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

#include "bng_defines.h"

#include "cplx_instance.h"

namespace BNG {

class BNGData;

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


// BNG reaction rule
// rules are only unidirectional,
// if there is a reversible reaction in BNGL definition,
// two
class RxnRule {
public:
  std::string name;

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
  // set only if the complex patterns are identical
  small_vector<CplxIndexPair> cplx_mapping;


  float_t rxn_rate;

  const CplxInstance& get_cplx_reactant(const uint index) const {
    assert(index <= reactants.size());
    return reactants[index];
  }

  const CplxInstance& get_cplx_product(const uint index) const {
    assert(index <= products.size());
    return products[index];
  }

  const MolInstance& get_mol_reactant(const CplxMolIndex& cmi) const {
    assert(cmi.mol_index <= reactants[cmi.cplx_index].mol_patterns.size());
    return get_cplx_reactant(cmi.cplx_index).mol_patterns[cmi.mol_index];
  }

  const MolInstance& get_mol_product(const CplxMolIndex& cmi) const {
    assert(cmi.mol_index <= products[cmi.cplx_index].mol_patterns.size());
    return get_cplx_product(cmi.cplx_index).mol_patterns[cmi.mol_index];
  }


  // mcell3 variant of maintaining substances,
  // e.g. for A + B -> A - A is maintained
  // TODO: merge with BNG style of maintaining substances
  bool is_cplx_reactant_on_both_sides_of_rxn(const uint index) const;

  bool is_cplx_product_on_both_sides_of_rxn(const uint index) const;


  bool operator ==(const RxnRule& rr2) {
    // ordering of components in a molecule is important
    // two component types must have the same id, this is ensured in find_or_add_component_type
    return
        name == rr2.name &&
        reactants == rr2.reactants && products == rr2.products &&
        rxn_rate == rr2.rxn_rate;
  }

  void dump(const BNGData& bng_data) const;

  // checks if it is possible to create a mapping from reactants to products and
  // sets members molecule_instances_are_maintained and mapping,
  // might write some error messages to the msgs stream,
  // returns true if errors were encountered
  bool compute_reactants_products_mapping(const BNGData& bng_data, std::ostream& out);

private:

  // returns false if cmi was not found in mapping,
  bool find_assigned_mol_reactant_for_product(const CplxMolIndex& product_cmi, CplxMolIndex& reactant_cmi) const;

  // check if it makes sense to compute mapping at all
  bool has_same_mols_in_reactants_and_products() const;

  // returns false if no fitting product was found
  bool find_most_fitting_unassigned_mol_product(const CplxMolIndex& reactant_cmi, CplxMolIndex& best_product_cmi) const;

  bool compute_mol_reactants_products_mapping(const BNGData& bng_data, std::ostream& out);


  bool find_assigned_cplx_reactant_for_product(const uint product_index, uint& reactant_index) const;
  void compute_cplx_reactants_products_mapping(const BNGData& bng_data);


  void dump_complex_instance_vector(
      const BNGData& bng_data,
      const CplxInstanceVector& complexes) const;
};

} /* namespace BNG */

#endif /* LIBS_BNG_RXN_RULE_H_ */
