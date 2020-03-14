/*
 * semantic_analyzer.cpp
 *
 *  Created on: Mar 13, 2020
 *      Author: ahusar
 */

#include "semantic_analyzer.h"
#include "parser_utils.h"

using namespace std;

namespace BNG {

// returns new node, owned by ctx if any new nodes were created
// called recursively, the set used_ids is copied intentionally every call
ASTExprNode* SemanticAnalyzer::evaluate_to_dbl(ASTExprNode* root, set<string> used_ids) {
  if (root->is_dbl()) {
    // already computed
    return root;
  }

  if (root->is_llong()) {
    return ctx->new_dbl_node(root->get_llong(), root);
  }

  if (root->is_id()) {
    const string& id = root->get_id();
    if (used_ids.count(id) != 0) {
      // TODO: test for this case
      errs(root) <<
          "Cyclic dependence while evaluating an expression, id '" << id << "' was already used.\n";
      ctx->inc_error_count();
      return ctx->new_dbl_node(0, root);
    }

    // find the value in the symbol table
    ASTBaseNode* val = ctx->symtab.get(root->get_id(), root, ctx);
    if (val == nullptr) {
      // error msg was printed by the symbol table
      return ctx->new_dbl_node(0, root);
    }
    if (!val->is_expr()) {
      errs(root) <<
          "Referenced id '" << id << "' cannot be used in an expression.\n";
      ctx->inc_error_count();
      return ctx->new_dbl_node(0, root);
    }

    used_ids.insert(id);
    return evaluate_to_dbl(to_expr(val), used_ids);
  }

  assert(false && "unreachable");
  return nullptr;
}

// compute all reaction rates and replace them with a floating point value (BNGL uses integers only for bond indices)
void SemanticAnalyzer::resolve_rxn_rates() {
  // the only place where expressions are currently used are rates
  // (parameters are evaluated along the way)
  // we are also checking that

  // for each rxn rule
  for (ASTBaseNode* n: ctx->rxn_rules.items) {
    ASTRxnRuleNode* rule = to_rxn_rule(n);

    // do we have the right number of rates?
    if (rule->reversible && rule->rates->size() != 2) {
      errs(rule)
          << "A reversible rule must have exactly 2 rates, "
          << rule->rates->size() << " rates provided.\n";
      ctx->inc_error_count();
    }

    if (!rule->reversible && rule->rates->size() != 1) {
      errs(rule)
          << "A reversible rule must have exactly 2 rates, "
          << rule->rates->size() << " rates provided.\n";
      errs(rule) << "A unidirectional rule must have exactly 1 rate.";
      ctx->inc_error_count();
    }

    // replace each of its rates with a float constant
    for (size_t i = 0; i < rule->rates->items.size(); i++) {
      ASTExprNode* orig_expr = to_expr(rule->rates->items[i]);
      ASTExprNode* new_expr = evaluate_to_dbl(orig_expr);
      // all nodes are owned by context and deleted after parsing has finished
      rule->rates->items[i] = new_expr;
    }
  }
}


// returns true if conversion and semantic checks passed
bool SemanticAnalyzer::check_and_convert(ASTContext* ctx_, BNGData* res_bng) {
  assert(ctx_ != nullptr);
  assert(res_bng != nullptr);

  ctx = ctx_;
  res = res_bng;

  // first compute all expressions - there are only rxn rates for now
  resolve_rxn_rates();
  if (ctx->get_error_count() != 0) {
    return false;
  }

  // convert molecule types

  // convert rxn rules


  return true;
}

} /* namespace BNG */
