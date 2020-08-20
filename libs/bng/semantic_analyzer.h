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

#include "bng/ast.h"
#include "bng/bng_engine.h"

namespace BNG {

class SemanticAnalyzer {
public:
  // modifies context e.g. by resolving expressions
  bool check_and_convert_parsed_file(ParserContext* ctx_, BNGData* res_bng);

  bool check_and_convert_single_cplx_instance(
      ParserContext* ctx_, BNGData* res_bng, CplxInstance& res);

private:
  ASTExprNode* evaluate_to_dbl(ASTExprNode* root, std::set<std::string> used_ids={});
  void resolve_rxn_rates();

  void convert_parameters();

  state_id_t convert_state_name(const ASTStrNode* s);
  component_type_id_t convert_component_type(const ASTComponentNode* c);
  MolType convert_molecule_type(
      const ASTMoleculeNode* n,
      const bool allow_same_component_different_state = false
  );
  void convert_and_store_molecule_types();

  void merge_molecule_type_definition(MolType& dstsrc, const MolType& src);
  void collect_molecule_types_molecule_list(
      const ASTListNode* molecule_list,
      std::vector<const ASTMoleculeNode*>& molecule_nodes
  );
  void collect_and_store_implicit_molecule_types();


  MolInstance convert_molecule_pattern(const ASTMoleculeNode* m);
  void convert_complex_pattern(
      const small_vector<const ASTMoleculeNode*>& complex_nodes,
      CplxInstance& pattern
  );
  void convert_cplx_inst_or_rxn_rule_side(
      const ASTListNode* rule_side,
      const bool convert_single_cplx_inst,
      CplxInstanceVector& pattern
  );

  void finalize_and_store_rxn_rule(const ASTRxnRuleNode* n, RxnRule& r, const bool forward_direction);
  void convert_and_store_rxn_rules();

  void convert_seed_species();

  // local copies so that we don't have to pass everything
  // as arguments
  ParserContext* ctx;
  BNGData* bng_data;
};

} /* namespace BNG */

#endif /* LIBS_BNG_SEMANTIC_ANALYZER_H_ */
