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

#include "ast.h"
#include "bng_engine.h"

namespace BNG {

class SemanticAnalyzer {
public:
  // modifies context e.g. by resolving expressions
  bool check_and_convert(ParserContext* ctx_, BNGData* res_bng);

private:
  ASTExprNode* evaluate_to_dbl(ASTExprNode* root, std::set<std::string> used_ids={});
  void resolve_rxn_rates();

  state_id_t convert_state_name(const ASTStrNode* s);
  component_type_id_t convert_component_type(const ASTComponentNode* c);
  MolType convert_molecule_type(const ASTMoleculeNode* n);
  void convert_and_store_molecule_types();

  MolInstance convert_molecule_pattern(const ASTMoleculeNode* m);
  void convert_complex_pattern(const small_vector<const ASTMoleculeNode*>& complex_nodes, CplxInstance& pattern);
  void convert_rxn_rule_side(const ASTListNode* rule_side, CplxInstanceVector& pattern);
  void finalize_and_store_rxn_rule(const ASTRxnRuleNode* n, RxnRule& r, const char* direction_str);
  void convert_and_store_rxn_rules();

  // local copies so that we don't have to pass everything
  // as arguments
  ParserContext* ctx;
  BNGData* bng_data;
};

} /* namespace BNG */

#endif /* LIBS_BNG_SEMANTIC_ANALYZER_H_ */
