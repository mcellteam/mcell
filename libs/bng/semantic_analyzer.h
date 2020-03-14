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
  bool check_and_convert(ASTContext* ctx_, BNGData* res_bng);

private:
  ASTExprNode* evaluate_to_dbl(ASTExprNode* root, std::set<std::string> used_ids={});
  void resolve_rxn_rates();

  state_id_t convert_state_name(const ASTStrNode* s);
  component_type_id_t convert_component_type(const ASTComponentNode* c);
  void convert_molecule_types();

  // map between information from molecule types and AST nodes
  // to be able to determine original source code location

  // local copies so that we don't have to pass everything
  // as arguments
  ASTContext* ctx;
  BNGData* bng_data;
};

} /* namespace BNG */

#endif /* LIBS_BNG_SEMANTIC_ANALYZER_H_ */
