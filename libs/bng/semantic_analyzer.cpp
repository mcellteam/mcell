/*
 * semantic_analyzer.cpp
 *
 *  Created on: Mar 13, 2020
 *      Author: ahusar
 */

#include "semantic_analyzer.h"
#include "parser_utils.h"

namespace BNG {



// compute all reaction rates and replace them with a floating point value (BNGL uses integers only for bond indices)
void SemanticAnalyzer::resolve_rxn_rates() {
  // the only place where expressions are currently used are rates
  // (parameters are evaluated along the way)
  // we are also checking that

  // for each rxn rule
  /*for () {
    // do we have the right number of rates?

    // for each of its rate
    for () {

    }

  }*/

}


// returns true if conversion and semantic checks passed
bool SemanticAnalyzer::check_and_convert(ASTContext* ctx_, BNGData* bng_) {
  assert(ctx_ != nullptr);
  assert(bng_ != nullptr);

  ctx = ctx_;
  res = bng_;

  // first compute all expressions - only rxn rates for now
  resolve_rxn_rates();
  if (ctx->get_error_count() != 0) {
    return false;
  }

  // convert molecule types

  // convert rxn rules
}

} /* namespace BNG */
