/*
 * semantic_analyzer.h
 *
 *  Created on: Mar 13, 2020
 *      Author: ahusar
 */

#ifndef LIBS_BNG_SEMANTIC_ANALYZER_H_
#define LIBS_BNG_SEMANTIC_ANALYZER_H_

#include "ast.h"
#include "bng_engine.h"

namespace BNG {

class SemanticAnalyzer {
public:
  // modifies context e.g. by resolving expresions
  bool check_and_convert(ASTContext* ctx, BNGData* res);

private:
  void resolve_rxn_rates();

  // local copy so that we son;t have to send everything
  const ASTContext* ctx;
  BNGData* res;
};

} /* namespace BNG */

#endif /* LIBS_BNG_SEMANTIC_ANALYZER_H_ */
