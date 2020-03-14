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


  // local copies so that we don't have to pass everything
  // as arguments
  ASTContext* ctx;
  BNGData* res;
};

} /* namespace BNG */

#endif /* LIBS_BNG_SEMANTIC_ANALYZER_H_ */
