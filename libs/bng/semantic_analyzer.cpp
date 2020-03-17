/*
 * semantic_analyzer.cpp
 *
 *  Created on: Mar 13, 2020
 *      Author: ahusar
 */

#include "semantic_analyzer.h"
#include "parser_utils.h"

#include "molecule_type.h"

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
    return evaluate_to_dbl(to_expr_node(val), used_ids);
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
    ASTRxnRuleNode* rule = to_rxn_rule_node(n);

    // do we have the right number of rates?
    if (rule->reversible && rule->rates->size() != 2) {
      errs(rule) <<
          "A reversible rule must have exactly 2 rates, "<<
          rule->rates->size() << " rate(s) provided.\n";
      ctx->inc_error_count();
    }

    if (!rule->reversible && rule->rates->size() != 1) {
      errs(rule) <<
          "A unidirectional rule must have exactly 1 rate, " <<
          rule->rates->size() << " rate(s) provided.\n";
      ctx->inc_error_count();
    }

    // replace each of its rates with a float constant
    for (size_t i = 0; i < rule->rates->items.size(); i++) {
      ASTExprNode* orig_expr = to_expr_node(rule->rates->items[i]);
      ASTExprNode* new_expr = evaluate_to_dbl(orig_expr);
      // all nodes are owned by context and deleted after parsing has finished
      rule->rates->items[i] = new_expr;
    }
  }
}


state_id_t SemanticAnalyzer::convert_state_name(const ASTStrNode* s) {
  assert(s != nullptr);
  return bng_data->find_or_add_state_name(s->str);
}


component_type_id_t SemanticAnalyzer::convert_component_type(const ASTComponentNode* c) {

  ComponentType ct;
  ct.name = c->name;

  // states
  for (const ASTBaseNode* state: c->states->items) {
    state_id_t state_id = convert_state_name(to_str_node(state));
    ct.allowed_state_ids.insert(state_id);
  }

  // bond - ignored, only error is printed
  if (c->bond->str != "") {
    errs(c) <<
        "Definition of a component in the molecule types section must not have a bond, "
        "error for '!" << c->bond->str << "'.\n";
    ctx->inc_error_count();
  }

  return bng_data->find_or_add_component_type(ct);
}


MoleculeType SemanticAnalyzer::convert_molecule_type(const ASTMoleculeNode* n) {
  MoleculeType mt;
  mt.name = n->name;

  for (const ASTBaseNode* c: n->components->items) {
    // component types must not have bonds
    component_type_id_t component_type_id =
        convert_component_type(to_component_node(c));
    mt.component_type_ids.push_back(component_type_id);
  }

  return mt;
}

void SemanticAnalyzer::convert_and_store_molecule_types() {

  // for each molecule (type) from the symbol table
  const ASTSymbolTable::IdToNodeMap& table = ctx->symtab.get_as_map();
  for (const auto& it: table) {
    assert(it.second != nullptr);
    if (it.second->is_molecule()) {
      const ASTMoleculeNode* n = to_molecule_node(it.second);

      MoleculeType mt = convert_molecule_type(n);

      bng_data->find_or_add_molecule_type(mt);
    }
  }
}


component_type_id_t SemanticAnalyzer::convert_component_instance(const ASTComponentNode* c) {

  ComponentType ct;
  ct.name = c->name;

  // states
  for (const ASTBaseNode* state: c->states->items) {
    state_id_t state_id = convert_state_name(to_str_node(state));
    ct.allowed_state_ids.insert(state_id);
  }

  // bond - ignored, only error is printed
  if (c->bond->str != "") {
    errs(c) <<
        "Definition of a component in the molecule types section must not have a bond, "
        "error for '!" << c->bond->str << "'.\n";
    ctx->inc_error_count();
  }

  return bng_data->find_or_add_component_type(ct);
}


ComponentInstance SemanticAnalyzer::convert_component_instance(const ASTComponentNode* n) {

  // check on allowed state?
}


MoleculeTypeInstance SemanticAnalyzer::convert_molecule_type_pattern(const ASTMoleculeNode* m) {

  MoleculeTypeInstance mtp;

  // process and remember ID
  molecule_type_id_t molecule_type_id = bng_data->find_molecule_type_id(m->name);
  if (molecule_type_id == MOLECULE_TYPE_ID_INVALID) {
    errs(m) << "Molecule type with name '" + m->name + "' was not defined.\n";
    ctx->inc_error_count();
    continue;
  }

  mtp.molecule_type_id = molecule_type_id;

  // process component instances for this molecule type pattern/instance

  for (size_t i = 0; i < complex_nodes.size(); i++) {
    //component_instances
    const ASTComponentNode* n = to_component_node(m->components[i]);
  }

  // TODO: semantic check that this is a really an allowed pattern for molecule_type_id
}

// for a pattern it is ok to not to list all components
void SemanticAnalyzer::convert_rxn_complex_pattern(const small_vector<const ASTMoleculeNode*>& complex_nodes, ComplexSpeciesInstance& pattern) {

  for (const ASTMoleculeNode* m: complex_nodes) {
    // molecule ids are based on their name
    pattern.push_back( convert_molecule_type_pattern(m) );
  }

  // TODO: semantic checks on bonds validity

}


// take one side of a reaction rule and create pattern for rule matching
void SemanticAnalyzer::convert_rxn_rule_side(const ASTListNode* rule_side, ComplexSpeciesInstanceVector& patterns) {

  // we need to check each molecule type from each complex
  small_vector<const ASTMoleculeNode*> current_complex_nodes;
  for (size_t i = 0; i < rule_side->items.size(); i++) {
    const ASTBaseNode* n = rule_side->items[i];

    if (n->is_separator()) {
      const ASTSeparatorNode* sep = to_separator(n);

      // separator must separate molecule type patterns,
      // also the last item must not be a separator
      if (current_complex_nodes.empty() || i == rule_side->items.size() - 1) {
        errs(n) <<
            "Invalid use of reaction rule separator '" << sep->to_char() << "'. "
            "It must be used to separate molecule type patterns.\n";
        ctx->inc_error_count();
        continue;
      }

      if (sep->is_plus()) {
        ComplexSpeciesInstance pattern;
        convert_rxn_complex_pattern(current_complex_nodes, pattern);
        patterns.push_back(pattern);

        // start a new complex pattern
        current_complex_nodes.clear();
      }
      // if no need to do anything for '.' -
    }
    else {
      // must be molecule type if not a separator
      const ASTMoleculeNode* m = to_molecule_node(n);
      current_complex_nodes.push_back(m);
    }
  }

  // process final complex
  assert(current_complex_nodes.empty() == rule_side->items.empty() && "Last set can be empty only if the patterns are empty");
  if (!current_complex_nodes.empty()) {
    ComplexSpecies pattern;
    convert_rxn_complex_pattern(current_complex_nodes, pattern);
    patterns.push_back(pattern);
  }
}


void SemanticAnalyzer::convert_and_store_rxn_rules() {

  // for each reaction rule
  for (const ASTBaseNode* n: ctx->rxn_rules.items) {
    const ASTRxnRuleNode* r = to_rxn_rule_node(n);

    ComplexSpeciesInstanceVector reactants;
    convert_rxn_rule_side(r->reactants, reactants);

    ComplexSpeciesInstanceVector products;
    convert_rxn_rule_side(r->products, products);

    RxnRule fwd_rule;
    fwd_rule.name = r->name;
    assert(r->rates->items.size() >= 1);
    fwd_rule.reaction_rate = to_expr_node(r->rates->items[0])->get_dbl();
    fwd_rule.reactants = reactants;
    fwd_rule.products = products;

    if (r->reversible) {
      RxnRule rev_rule;
      rev_rule.name = r->name;

      if (products.empty()) {
        errs(r) <<
            "A reversible rule must have at least one complex pattern " <<
            "on the right side of the reaction rule.\n";
        ctx->inc_error_count();
      }

      assert(r->rates->items.size() == 2);
      fwd_rule.reaction_rate = to_expr_node(r->rates->items[1])->get_dbl();
      fwd_rule.reactants = products;
      fwd_rule.products = reactants;
    }
  }
}



// returns true if conversion and semantic checks passed
bool SemanticAnalyzer::check_and_convert(ASTContext* ctx_, BNGData* res_bng) {
  assert(ctx_ != nullptr);
  assert(res_bng != nullptr);

  ctx = ctx_;
  bng_data = res_bng;

  // first compute all expressions - there are only rxn rates for now
  resolve_rxn_rates();
  if (ctx->get_error_count() != 0) {
    return false;
  }

  // convert molecule types
  // the single purpose of molecule types is to be able to check
  // that reaction rules and also new releases adhere to the molecule type template
  convert_and_store_molecule_types();
  if (ctx->get_error_count() != 0) {
    return false;
  }

  // convert rxn rules
  convert_and_store_rxn_rules();
  if (ctx->get_error_count() != 0) {
    return false;
  }


  return true;
}

} /* namespace BNG */
