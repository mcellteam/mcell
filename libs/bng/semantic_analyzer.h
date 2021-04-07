/*
 * semantic_analyzer.h
 *
 *  Created on: Mar 13, 2020
 *      Author: ahusar
 */

#ifndef LIBS_BNG_SEMANTIC_ANALYZER_H_
#define LIBS_BNG_SEMANTIC_ANALYZER_H_

#include <set>
#include <string>
#include <map>

#include "bng/ast.h"
#include "bng/bng_engine.h"

namespace BNG {

class SemanticAnalyzer {
public:
  // modifies context e.g. by resolving expressions
  bool check_and_convert_parsed_file(
      ParserContext* ctx_,
      BNGData* res_bng,
      const std::map<std::string, float_t>& parameter_overrides = std::map<std::string, float_t>()
  );

  bool check_and_convert_single_cplx(
      ParserContext* ctx_, BNGData* res_bng, Cplx& res);

private:
  double evaluate_function_call(ASTExprNode* call_node, const std::vector<double>& arg_values);
  ASTExprNode* evaluate_to_dbl(ASTExprNode* root, std::set<std::string> used_ids={});
  void resolve_rxn_rates();

  void convert_and_evaluate_parameters(const std::map<std::string, float_t>& parameter_overrides);

  state_id_t convert_state_name(const ASTStrNode* s);
  ComponentType convert_component_type(
      const std::string& elem_mol_type_name,
      const ASTComponentNode* c,
      const bool allow_components_to_have_bonds = false
  );
  ElemMolType convert_molecule_type(
      const ASTMolNode* n,
      const bool parsing_single_cplx = false
  );
  void convert_and_store_molecule_types();

  void convert_and_store_compartments();

  void merge_molecule_type_definition(ElemMolType& dstsrc, const ElemMolType& src);
  void collect_molecule_types_molecule_list(
      const ASTListNode* molecule_list,
      std::vector<const ASTMolNode*>& molecule_nodes
  );
  void collect_and_store_implicit_molecule_types();


  ElemMol convert_molecule_pattern(const ASTMolNode* m);
  void convert_cplx(
      const ASTCplxNode* cplx_node,
      Cplx& pattern,
      const bool in_rule_or_observable,
      const bool check_compartments = true
  );
  void convert_rxn_rule_side(
      const ASTListNode* rule_side,
      const bool reactants_side,
      CplxVector& patterns
  );

  void finalize_and_store_rxn_rule(const ASTRxnRuleNode* n, RxnRule& r, const bool forward_direction);
  void convert_and_store_rxn_rules();

  void convert_seed_species();
  void convert_observables();

  void extend_molecule_type_definitions(const ASTCplxNode* cplx_node);

  void define_compartments_used_by_cplx_as_3d_compartments(ASTCplxNode* cplx_node);

  // local copies so that we don't have to pass everything
  // as arguments
  ParserContext* ctx;
  BNGData* bng_data;
};


// auxiliary function also used in Cplx
void insert_compartment_id_to_set_based_on_type(
    const BNGData* bng_data,
    const compartment_id_t cid,
    bool& all_are_none_or_inout,
    bool& has_compartment_none,
    uint_set<compartment_id_t>& vol_compartments,
    uint_set<compartment_id_t>& surf_compartments);

} /* namespace BNG */

#endif /* LIBS_BNG_SEMANTIC_ANALYZER_H_ */
