/*
 * rule.h
 *
 *  Created on: Jan 9, 2020
 *      Author: ahusar
 */

#ifndef LIBS_BNG_RXN_RULE_H_
#define LIBS_BNG_RXN_RULE_H_

#include <string>

#include "bng_defines.h"

#include "cplx_instance.h"

namespace BNG {

class BNGData;

struct CplxMolIndex {
  CplxMolIndex(const uint complex_index_, const uint molecule_index_)
    : cplx_index(complex_index_), mol_index(molecule_index_) {
  }

  uint cplx_index;
  uint mol_index;
};

struct CMIndexPair {
  CMIndexPair(
      const CplxMolIndex& reactant_cmi_, const CplxMolIndex& product_cmi_)
    : reactant_cmi(reactant_cmi_), product_cmi(product_cmi_) {
  }

  CplxMolIndex reactant_cmi;
  CplxMolIndex product_cmi;
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
  ComplexInstanceVector reactants;
  ComplexInstanceVector products;

  // set to true if it was possible to do a mapping between reactants and products
  bool mol_instances_are_maintained;

  // matching between molecules of reactants and molecules of products,
  // contains information only if molecule_instances_are_maintained is true
  // TODO: private
  small_vector<CMIndexPair> mapping;


  float_t rxn_rate;

  const MolInstance& get_reactant_mol(const CplxMolIndex& cmi) {
    assert(cmi.cplx_index <= reactants.size());
    assert(cmi.mol_index <= reactants[cmi.cplx_index].mol_patterns.size());
    return reactants[cmi.cplx_index].mol_patterns[cmi.mol_index];
  }

  const MolInstance& get_product_mol(const CplxMolIndex& cmi) {
    assert(cmi.cplx_index <= products.size());
    assert(cmi.mol_index <= products[cmi.cplx_index].mol_patterns.size());
    return products[cmi.cplx_index].mol_patterns[cmi.mol_index];
  }

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
  bool compute_reactants_products_mapping(std::stringstream& msgs);

private:
  // returns nullptr if cmi was not found in mapping
  CplxMolIndex* get_assigned_reactant_for_product(CplxMolIndex& product_cmi);

  // check if it makes sense to compute mapping at all
  bool has_same_molecules_in_reactants_and_products();

  void find_most_fitting_unassigned_product(CplxMolIndex& cmi);


  void dump_complex_instance_vector(
      const BNGData& bng_data,
      const ComplexInstanceVector& complexes) const;
};

} /* namespace BNG */

#endif /* LIBS_BNG_RXN_RULE_H_ */
