/*
 * semantic_analyzer.cpp
 *
 *  Created on: Mar 13, 2020
 *      Author: ahusar
 */

#include <sstream>

#include "bng/semantic_analyzer.h"

#include "bng/bng_data.h"
#include "bng/mol_type.h"
#include "bng/parser_utils.h"


// each semantic check has in comment the name of a test for it

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
      errs_loc(root) <<
          "Cyclic dependence while evaluating an expression, id '" << id << "' was already used.\n"; // test N0012
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
      errs_loc(root) <<
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
      errs_loc(rule) <<
          "A reversible rule must have exactly 2 rates, "<<
          rule->rates->size() << " rate(s) provided.\n"; // test N0015
      ctx->inc_error_count();
    }

    if (!rule->reversible && rule->rates->size() != 1) {
      errs_loc(rule) <<
          "A unidirectional rule must have exactly 1 rate, " <<
          rule->rates->size() << " rate(s) provided.\n"; // test N0014
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
    errs_loc(c) <<
        "Definition of a component in the molecule types section must not have a bond, "
        "error for '!" << c->bond->str << "'.\n"; // test N0100
    ctx->inc_error_count();
  }

  return bng_data->find_or_add_component_type(ct);
}


MolType SemanticAnalyzer::convert_molecule_type(const ASTMoleculeNode* n) {
  MolType mt;
  mt.name = n->name;

  for (const ASTBaseNode* c: n->components->items) {
    // component types must not have bonds
    component_type_id_t component_type_id =
        convert_component_type(to_component_node(c));
    mt.component_type_ids.push_back(component_type_id);
  }

  // when there a components with the same name in on molecule, they must have the same allowed states
  for (size_t i = 0; i < mt.component_type_ids.size(); i++) {

    component_type_id_t ct_i_id = mt.component_type_ids[i];
    const ComponentType& ct_i = bng_data->get_component_type(ct_i_id);

    for (size_t k = i+1; k < mt.component_type_ids.size(); k++) {

      component_type_id_t ct_k_id = mt.component_type_ids[k];
      const ComponentType& ct_k = bng_data->get_component_type(ct_k_id);

      if (ct_i_id != ct_k_id && ct_i.name == ct_k.name) {
        // this is directly a conflict because if the components would have the same name and the same components,
        // their id would be the same
        errs_loc(n) <<
            "Molecule type has 2 components with name '" << ct_i.name << "' but with different states, this is not allowed.\n"; // test N0101
        ctx->inc_error_count();
      }
    }
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

      MolType mt = convert_molecule_type(n);

      bng_data->find_or_add_molecule_type(mt);
    }
  }
}


MolInstance SemanticAnalyzer::convert_molecule_pattern(const ASTMoleculeNode* m) {

  MolInstance mp;
  // there is no support for surface molecules in BNGL yet, so everything must be volume molecule
  mp.set_is_vol();

  // process and remember ID
  mol_type_id_t molecule_type_id = bng_data->find_molecule_type_id(m->name);
  if (molecule_type_id == MOL_TYPE_ID_INVALID) {
    errs_loc(m) << "Molecule type with name '" + m->name + "' was not defined.\n"; // test N0200
    ctx->inc_error_count();
    return mp;
  }

  const MolType& mt = bng_data->get_molecule_type(molecule_type_id);
  mp.mol_type_id = molecule_type_id;

  mp.initialize_components_types(mt);

  uint current_component_index = 0;

  // process component instances for this molecule type pattern/instance
  for (size_t i = 0; i < m->components->items.size(); i++) {
    // component_instances
    const ASTComponentNode* component = to_component_node(m->components->items[i]);

    // which component in the type it is?
    // E.g. when the type is A(a, a, a, b, c, b) then the pattern mathes like this:
    //                       A(a, a,    b, c   ),
    current_component_index = mp.get_corresponding_component_index(*bng_data, mt, component->name, current_component_index);
    if (current_component_index == INDEX_INVALID) {
      errs_loc(m) <<
          "Molecule type '" << mt.name << "' does not declare component '" << component->name <<
          "' or all instances of this component were already used.\n"; // test N0201
      ctx->inc_error_count();
      return mp;
    }

    component_type_id_t component_type_id = mt.component_type_ids[current_component_index];

    // state
    if (component->states->items.size() > 1) {
      errs_loc(component) << "A component might have max. 1 state specified, error for component " << component->name << ".\n"; // test N0202
      ctx->inc_error_count();
      return mp;
    }

    // check and set if state is specified
    if (component->states->items.size() == 1) {
      const string& state_name = to_str_node(component->states->items[0])->str;

      // does this state exist at all?
      state_id_t state_id = bng_data->find_state_id(state_name);
      if (state_id == STATE_ID_INVALID) {
        errs_loc(component) <<
            "Unknown state name '" << state_name << "' for component '" << component->name <<
            "' (this state name was not found for any declared component).\n"; // test N0203
        ctx->inc_error_count();
        return mp;
      }

      // is this state allowed for this component?
      const ComponentType& ct = bng_data->get_component_type(component_type_id);
      if (ct.allowed_state_ids.count(state_id) == 0) {
        errs_loc(component) <<
            "State name '" << state_name << "' was not declared as allowed for component '" << component->name << "'.\n"; // test N0204
        ctx->inc_error_count();
        return mp;
      }

      // finally set the component's state
      mp.component_instances[current_component_index].state_id = state_id;
    }

    // bond
    const string& bond = component->bond->str;
    bond_value_t b = str_to_bond_value(bond);
    if (b == BOND_VALUE_INVALID) {
      // this is already checked by parser, so the a test checks parser message but let's keep this
      // parser test: N0205
      ctx->internal_error(component->bond, "Invalid bond index");
    }

    mp.component_instances[current_component_index].bond_value = b;

    // we also need to remember for further checks that this component was explicitly
    // specified
    mp.component_instances[current_component_index].explicitly_listed_in_pattern = true;

    // need to move component index because we processed this component
    current_component_index++;
  }

  return mp;
}


// for a pattern it is ok to not to list all components
void SemanticAnalyzer::convert_complex_pattern(const small_vector<const ASTMoleculeNode*>& complex_nodes, CplxInstance& pattern) {

  for (const ASTMoleculeNode* m: complex_nodes) {
    // molecule ids are based on their name
    pattern.mol_instances.push_back( convert_molecule_pattern(m) );
  }

  // semantic checks on bonds validity
  map<bond_value_t, vector<uint> > bond_value_to_molecule_index;
  for (uint i = 0; i < pattern.mol_instances.size(); i++) {
    const MolInstance& mi = pattern.mol_instances[i];

    for (const ComponentInstance& compi: mi.component_instances) {
      if (compi.bond_has_numeric_value()) {
        // remember molecule index in this complex for a given bond
        bond_value_to_molecule_index[compi.bond_value].push_back(i);
      }
    }
  }

  for (const auto& it: bond_value_to_molecule_index) {
    // each bond is used exactly twice
    if (it.second.size() != 2) {
      assert(complex_nodes.size() > 0);
      errs_loc(complex_nodes[0]) <<
          "Bond with numerical value '" << it.first << "' must be used exactly twice in a complex pattern of a rule.\n"; // test N0206
      ctx->inc_error_count();
      return;
    }

    // it is used in different molecules of a complex
    if (it.second[0] == it.second[1]) {
      assert(complex_nodes.size() > 0);
      errs_loc(complex_nodes[0]) <<
          "Bond with numerical value '" << it.first << "' must bind different molecules of a complex pattern of a rule.\n"; // test N0208
      ctx->inc_error_count();
      return;
    }
  }

  pattern.finalize();
}


// take one side of a reaction rule and create pattern for rule matching
void SemanticAnalyzer::convert_rxn_rule_side(const ASTListNode* rule_side, CplxInstanceVector& patterns) {

  // we need to check each molecule type from each complex
  small_vector<const ASTMoleculeNode*> current_complex_nodes;
  for (size_t i = 0; i < rule_side->items.size(); i++) {
    const ASTBaseNode* n = rule_side->items[i];

    if (n->is_separator()) {
      const ASTSeparatorNode* sep = to_separator(n);

      // separator must separate molecule type patterns,
      // also the last item must not be a separator
      if (current_complex_nodes.empty() || i == rule_side->items.size() - 1) {
        errs_loc(n) <<
            "Invalid use of reaction rule separator '" << sep->to_char() << "'. "
            "It must be used to separate molecule type patterns.\n"; // should be caught by parser, test N0210
        ctx->inc_error_count();
        continue;
      }

      if (sep->is_plus()) {
        CplxInstance pattern;
        convert_complex_pattern(current_complex_nodes, pattern);
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
    CplxInstance pattern;
    convert_complex_pattern(current_complex_nodes, pattern);
    patterns.push_back(pattern);
  }
}


void SemanticAnalyzer::finalize_and_store_rxn_rule(const ASTRxnRuleNode* n, RxnRule& r, const char* direction_str) {
  // determine mapping from molecule instances on one side to another
  stringstream out;
  bool ok = r.compute_reactants_products_mapping_w_error_output(*bng_data, out);
  if (!ok) {
    errs_loc(n) << out.str() << " (in the " << direction_str << " direction, indices are counted from 0)\n"; // tests N0220, N0230, N0231, N0232
    ctx->inc_error_count();
  }

  bng_data->find_or_add_rxn_rule(r);
}


void SemanticAnalyzer::convert_and_store_rxn_rules() {

  // for each reaction rule
  for (const ASTBaseNode* n: ctx->rxn_rules.items) {
    const ASTRxnRuleNode* r = to_rxn_rule_node(n);

    CplxInstanceVector reactants;
    convert_rxn_rule_side(r->reactants, reactants);

    CplxInstanceVector products;
    convert_rxn_rule_side(r->products, products);

    RxnRule fwd_rule;
    fwd_rule.name = r->name;
    assert(r->rates->items.size() >= 1);
    fwd_rule.rate_constant = to_expr_node(r->rates->items[0])->get_dbl();
    fwd_rule.reactants = reactants;
    fwd_rule.products = products;

    finalize_and_store_rxn_rule(r, fwd_rule, "forward");

    if (r->reversible) {
      RxnRule rev_rule;
      rev_rule.name = r->name;

      if (products.empty()) {
        errs_loc(r) <<
            "A reversible rule must have at least one complex pattern " <<
            "on the right side of the reaction rule.\n"; // caught by parser, test N0211
        ctx->inc_error_count();
      }

      assert(r->rates->items.size() == 2);
      rev_rule.rate_constant = to_expr_node(r->rates->items[1])->get_dbl();
      rev_rule.reactants = products;
      rev_rule.products = reactants;

      finalize_and_store_rxn_rule(r, rev_rule, "reverse");
    }
  }
}


// returns true if conversion and semantic checks passed
bool SemanticAnalyzer::check_and_convert(ParserContext* ctx_, BNGData* res_bng) {
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