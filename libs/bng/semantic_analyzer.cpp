/*
 * semantic_analyzer.cpp
 *
 *  Created on: Mar 13, 2020
 *      Author: ahusar
 */

#include <sstream>
#include <cmath>

#include "bng/semantic_analyzer.h"

#include "bng/bng_data.h"
#include "bng/elem_mol_type.h"
#include "bng/parser_utils.h"
#include "bng/bngl_names.h"

// each semantic check has in comment the name of a test for it

using namespace std;

namespace BNG {

const char* const DIR_FORWARD = "forward";
const char* const DIR_REVERSE = "reverse";


// function names - name and number of arguments, list used by data model to pymcell4 converter
// is in generator_utils.h: mdl_functions_to_py_bngl_map
struct BnglFunctionInfo {
  uint num_arguments;
  double (*eval_1_arg_func_call)(double);
  double (*eval_2_args_func_call)(double, double);
};

static double bngl_max(const double a, const double b) {
  return (a<b)?b:a;
}

static double bngl_min(const double a, const double b) {
  return !(b<a)?a:b;
}

// TODO: check allowed argument ranges e.g. for asin
// TODO: only the intersection of MDL and BNGL functions is supported now, some BNGL functios are missing
static const std::map<std::string, BnglFunctionInfo> bngl_function_infos {
  { "sqrt", {1, sqrt, nullptr} },
  { "exp", {1, exp, nullptr} },
  { "ln", {1, log, nullptr} },
  { "log10", {1, log10, nullptr} },
  { "sin", {1, sin, nullptr} },
  { "cos", {1, cos, nullptr} },
  { "tan", {1, tan, nullptr} },
  { "asin", {1, asin, nullptr} },
  { "acos", {1, acos, nullptr} },
  { "atan", {1, atan, nullptr} },
  { "abs", {1, fabs, nullptr} },
  { "ceil", {1, ceil, nullptr} },
  { "floor", {1, floor, nullptr} },
  { "max", {2, nullptr, bngl_max} },
  { "min", {2, nullptr, bngl_min} }
};

static bool is_thrash_or_null(const string& name) {
  // same check as in nfsim
  return (name == COMPLEX_Null || name == COMPLEX_NULL || name == COMPLEX_null ||
      name == COMPLEX_Trash || name == COMPLEX_TRASH || name == COMPLEX_trash);
}


double SemanticAnalyzer::evaluate_function_call(ASTExprNode* call_node, const std::vector<double>& arg_values) {
  assert(call_node != nullptr);
  const string& name = call_node->get_function_name();
  const auto& func_info_it = bngl_function_infos.find(name);
  if (func_info_it != bngl_function_infos.end()) {
    if (arg_values.size() == func_info_it->second.num_arguments) {
      const BnglFunctionInfo& info = func_info_it->second;
      if (info.num_arguments == 1) {
        return info.eval_1_arg_func_call(arg_values[0]);
      }
      else if (info.num_arguments == 2) {
        return info.eval_2_args_func_call(arg_values[0], arg_values[1]);
      }
      else {
        assert(false);
        return 0;
      }
    }
    else {
      errs_loc(call_node) <<
          "Invalid number of arguments for function '" << name << "', got " << arg_values.size() <<
          " expected " << func_info_it->second.num_arguments << ".\n"; // test TODO
      ctx->inc_error_count();
      return 0;
    }
  }
  else {
    errs_loc(call_node) <<
        "Unknown function '" << name << "' encountered.\n"; // test TODO
    ctx->inc_error_count();
    return 0;
  }
}


// returns new node, owned by ctx if any new nodes were created
// called recursively, the set used_ids is copied intentionally every call
ASTExprNode* SemanticAnalyzer::evaluate_to_dbl(ASTExprNode* root, set<string> used_ids) {
  if (root->is_dbl()) {
    // already computed
    return root;
  }
  else if (root->is_llong()) {
    return ctx->new_dbl_node(root->get_llong(), root);
  }
  else if (root->is_id()) {
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
  else if (root->is_unary_expr()) {
    assert(root->get_left() != nullptr);
    assert(root->get_right() == nullptr);
    double res = evaluate_to_dbl(root->get_left(), used_ids)->get_dbl();
    return ctx->new_dbl_node( (root->get_op() == ExprType::UnaryMinus) ? -res : res, root);
  }
  else if (root->is_binary_expr()) {
    assert(root->get_left() != nullptr);
    assert(root->get_right() != nullptr);
    double res_left = evaluate_to_dbl(root->get_left(), used_ids)->get_dbl();
    double res_right = evaluate_to_dbl(root->get_right(), used_ids)->get_dbl();
    double res = 0;
    switch (root->get_op()) {
      case ExprType::Add:
        res = res_left + res_right;
        break;
      case ExprType::Sub:
        res = res_left - res_right;
        break;
      case ExprType::Mul:
        res = res_left * res_right;
        break;
      case ExprType::Div:
        if (res_right == 0) {
          errs_loc(root) << "Division by zero, left operand evaluated to " << res_left << " and right to 0.\n";
          ctx->inc_error_count();
        }
        else {
          res = res_left / res_right;
        }
        break;
      case ExprType::Pow:
        res = pow(res_left, res_right);
        break;
      default:
        release_assert(false && "Invalid operator");
    }
    return ctx->new_dbl_node(res, root);
  }
  else if (root->is_function_call()) {
    assert(root->get_args() != nullptr);
    // evaluate all arguments
    vector<double> arg_values;
    for (ASTBaseNode* base_arg_node: root->get_args()->items) {
      ASTExprNode* arg_node = to_expr_node(base_arg_node);
      arg_values.push_back(evaluate_to_dbl(arg_node, used_ids)->get_dbl());
    }
    double res = evaluate_function_call(root, arg_values);
    return ctx->new_dbl_node(res, root);
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


// the following conversions do not use the bng_data.parameters map
// if a parameter is in the map parameter_overrides, the supplied value is used instead
void SemanticAnalyzer::convert_and_evaluate_parameters(
    const std::map<std::string, float_t>& parameter_overrides) {

  // first go trough all parameter overrides and either change definition or add
  // this parameter to the symbol table
  ASTSymbolTable::IdToNodeMap& symtab_map = ctx->symtab.get_as_map();
  for (auto it_override: parameter_overrides) {
    // is defined?
    auto it_found_sym = symtab_map.find(it_override.first);
    if (it_found_sym != symtab_map.end()) {
      // and it is a symbol
      if (it_found_sym->second->is_expr()) {
        // override its value
        it_found_sym->second = ctx->new_dbl_node(it_override.second);
      }
      else {
        errs() <<
            "Cannot override symbol " << it_override.first << " that is not a parameter.\n";
        ctx->inc_error_count();
      }
    }
    else {
      // define as a new symbol
      ctx->symtab.insert(
          it_override.first,
          ctx->new_dbl_node(it_override.second),
          ctx
      );
    }
  }

  // every symbol that maps directly into a value is a parameter
  // evaluate them and store
  for (auto it_sym: ctx->symtab.get_as_map()) {
    if (it_sym.second->is_expr()) {
      ASTExprNode* orig_expr = to_expr_node(it_sym.second);
      ASTExprNode* new_expr = evaluate_to_dbl(orig_expr);

      bng_data->add_parameter(it_sym.first, new_expr->get_dbl());
    }
  }
}


state_id_t SemanticAnalyzer::convert_state_name(const ASTStrNode* s) {
  assert(s != nullptr);
  return bng_data->find_or_add_state_name(s->str);
}


ComponentType SemanticAnalyzer::convert_component_type(
    const std::string& elem_mol_type_name,
    const ASTComponentNode* c,
    const bool allow_components_to_have_bonds
) {

  ComponentType ct;
  ct.name = c->name;
  ct.elem_mol_type_name = elem_mol_type_name;

  // states
  for (const ASTBaseNode* state: c->states->items) {
    state_id_t state_id = convert_state_name(to_str_node(state));
    ct.allowed_state_ids.insert(state_id);
  }

  // bond - ignored, only error is printed
  if (!allow_components_to_have_bonds && c->bond->str != "") {
    errs_loc(c) <<
        "Definition of a component in the molecule types section must not have a bond, "
        "error for '!" << c->bond->str << "'.\n"; // checked by parser, original test N0100
    ctx->inc_error_count();
  }

  return ct;
}


ElemMolType SemanticAnalyzer::convert_molecule_type(
    const ASTMolNode* n,
    // when parsing single cplx we don't have the full information from the molecule types section
    // so we allow merging component types and also bonds
    const bool parsing_single_cplx
) {
  ElemMolType mt;
  mt.name = n->name;

  for (const ASTBaseNode* c: n->components->items) {
    const ASTComponentNode* comp = to_component_node(c);

    // did we already define a component for this molecule with this name?
    component_type_id_t ct_id = bng_data->find_component_type_id(mt, comp->name);
    if (ct_id == COMPONENT_TYPE_ID_INVALID) {
      // new component type
      ComponentType ct = convert_component_type(mt.name, to_component_node(c), parsing_single_cplx);
      component_type_id_t new_ct_id = bng_data->find_or_add_component_type(ct, parsing_single_cplx);
      mt.component_type_ids.push_back(new_ct_id);
    }
    else {
      // check or add states used in the new component against
      // what we defined before
      ComponentType& existing_ct = bng_data->get_component_type(ct_id);
      ComponentType new_ct = convert_component_type(mt.name, comp, parsing_single_cplx);
      assert(existing_ct.name == new_ct.name);
      if (!parsing_single_cplx) {
        if (existing_ct.allowed_state_ids != new_ct.allowed_state_ids) {
          errs_loc(n) <<
              "Molecule type has 2 components with name '" << existing_ct.name <<
              "' but with different states, this is not allowed.\n"; // test N0101
          ctx->inc_error_count();
        }
      }
      else {
        // while parsing a single cplx, we didn't get the molecule types definitions
        // so we do not know how the component definitions look like, let's merge the allowed states
        existing_ct.allowed_state_ids.insert(new_ct.allowed_state_ids.begin(), new_ct.allowed_state_ids.end());
      }

      // append it to the molecule type's components
      mt.component_type_ids.push_back(ct_id);
    }
  }

  return mt;
}


void SemanticAnalyzer::convert_and_store_molecule_types() {

  // for each molecule (type) from the symbol table
  const ASTSymbolTable::IdToNodeMap& table = ctx->symtab.get_as_map();
  for (const auto& it: table) {
    assert(it.second != nullptr);
    if (it.second->is_mol()) {
      const ASTMolNode* n = to_molecule_node(it.second);

      if (n->has_compartment()) {
        errs_loc(n) <<
            "Compartments cannot be used in 'molecule types' section, error for molecule '" << n->name << "'.\n"; // test N0302
        ctx->inc_error_count();
        continue;
      }
      ElemMolType mt = convert_molecule_type(n);

      bng_data->find_or_add_elem_mol_type(mt);
    }
  }
}


void SemanticAnalyzer::convert_and_store_compartments() {
  // first define all compartments without their parents and children
  for (size_t i = 0; i < ctx->compartments.items.size(); i++) {
    const ASTCompartmentNode* n = to_compartment_node(ctx->compartments.items[i]);
    Compartment c;

    // name
    if (bng_data->find_compartment_id(n->name) != COMPARTMENT_ID_INVALID) {
      errs_loc(n) <<
          "Compartment '" << n->name << "' was already defined.\n"; // test N0302
      ctx->inc_error_count();
      continue;
    }

    if (n->name == DEFAULT_COMPARTMENT_NAME) {
      errs_loc(n) <<
          "Compartment name '" << n->name << "' is reserved, the specification will be ignored.\n"; // test TODO
      continue;
    }

    c.name = n->name;

    // dimensions
    if (n->dimensions != 2 && n->dimensions != 3) {
      errs_loc(n) <<
          "Compartment '" << n->name << "' has invalid dimension of value " << n->dimensions <<
          ", the only values allowed are 2 or 3.\n"; // test N0300
      ctx->inc_error_count();
      continue;
    }
    c.is_3d = n->dimensions == 3;

    // volume
    ASTExprNode* evaluated_volume = evaluate_to_dbl(n->volume);
    float_t volume = evaluated_volume->get_dbl();
    if (volume < 0) {
      errs_loc(n) <<
          "Compartment '" << n->name << "' has negative volume " << volume << ".\n"; // test N0301
      ctx->inc_error_count();
      continue;
    }
    c.set_volume(volume);

    bng_data->add_compartment(c);
  }

  if (ctx->get_error_count()) {
    return;
  }

  // now define their parents and children
  for (size_t i = 0; i < ctx->compartments.items.size(); i++) {
    const ASTCompartmentNode* n = to_compartment_node(ctx->compartments.items[i]);

    Compartment* c = bng_data->find_compartment(n->name);
    if (c == nullptr) {
      // ignoring default_compartment
      assert(n->name == DEFAULT_COMPARTMENT_NAME);
      continue;
    }
    if (n->parent_name != "") {

      compartment_id_t parent_compartment_id = bng_data->find_compartment_id(n->parent_name);
      if (parent_compartment_id == COMPARTMENT_ID_INVALID) {
        errs_loc(n) <<
            "Compartment's '" << n->name << "' parent '" << n->parent_name << "' was not defined.\n"; // test N0303
        ctx->inc_error_count();
        continue;
      }
      // set parent
      c->parent_compartment_id = parent_compartment_id;

      Compartment& parent = bng_data->get_compartment(parent_compartment_id);

      // check
      if (c->is_3d == parent.is_3d) {
        errs_loc(n) <<
            "Parent compartment " + parent.name +
            " must be of different dimensionality than its child '" + c->name + "'."; // test TODO
        ctx->inc_error_count();
        continue;
      }

      // set 'c' as a child of its parent
      parent.children_compartments.insert(c->id);
    }
  }

  if (ctx->get_error_count()) {
    return;
  }

  for (size_t i = 0; i < ctx->compartments.items.size(); i++) {
    const ASTCompartmentNode* n = to_compartment_node(ctx->compartments.items[i]);

    // check that the lowest level children is 3d
    const Compartment* c = bng_data->find_compartment(n->name);
    if (c == nullptr) {
      // ignoring default_compartment
      assert(n->name == DEFAULT_COMPARTMENT_NAME);
      continue;
    }
    if (!c->has_children() && !c->is_3d) {
      errs() <<
          "Compartment without sub-compartments '" + c->name + "' must be a 3D compartment.\n"; // test TODO
      ctx->inc_error_count();
    }
  }
}


void SemanticAnalyzer::collect_molecule_types_molecule_list(
    const ASTListNode* cplx_list,
    vector<const ASTMolNode*>& molecule_nodes
) {
  for (size_t i = 0; i < cplx_list->items.size(); i++) {
    const ASTCplxNode* cplx = to_cplx_node(cplx_list->items[i]);

    molecule_nodes.insert(molecule_nodes.end(), cplx->mols.begin(), cplx->mols.end());
  }
}


void SemanticAnalyzer::collect_and_store_implicit_molecule_types() {
  // go through reaction rules and seed species to define molecule types

  // first collect all molecule nodes
  // for each rxn rule
  vector<const ASTMolNode*> found_mol_nodes;
  for (const ASTBaseNode* n: ctx->rxn_rules.items) {
    const ASTRxnRuleNode* r = to_rxn_rule_node(n);
    collect_molecule_types_molecule_list(r->reactants, found_mol_nodes);
    collect_molecule_types_molecule_list(r->products, found_mol_nodes);
  }

  for (const ASTBaseNode* n: ctx->seed_species.items) {
    const ASTSeedSpeciesNode* ss = to_seed_species_node(n);
    found_mol_nodes.insert(found_mol_nodes.end(),
        ss->cplx->mols.begin(), ss->cplx->mols.end());
  }

  // sort by name and skip those that are already known
  map<string, vector<const ASTMolNode*>> mol_uses_with_same_name;
  for (const ASTMolNode* n: found_mol_nodes) {
    elem_mol_type_id_t mt_id = bng_data->find_elem_mol_type_id(n->name);
    if (mt_id == MOL_TYPE_ID_INVALID) {
      mol_uses_with_same_name[n->name].push_back(n);
    }
  }

  // and merge into a single definition
  // for each different name
  for (auto same_name_it: mol_uses_with_same_name) {

    // create a map of used states per each component
    map<string, set<string>> component_state_names;
    // also count the maximum number of components
    map<string, uint> max_component_count_per_all_mts;
    for (const ASTMolNode* mn: same_name_it.second) {

      map<string, uint> max_component_count_per_single_mt;
      for (const ASTBaseNode* cnb: mn->components->items) {
        const ASTComponentNode* cn = to_component_node(cnb);

        // count occurrence
        if (max_component_count_per_single_mt.count(cn->name) == 0) {
          max_component_count_per_single_mt[cn->name] = 1;
        }
        else {
          max_component_count_per_single_mt[cn->name]++;
        }

        // collect state names
        for (const ASTBaseNode* snb: cn->states->items) {
          const ASTStrNode* sn = to_str_node(snb);
          component_state_names[cn->name].insert(sn->str);
        }
      }

      // update max count of these components
      for (auto single_max: max_component_count_per_single_mt) {
        const string& comp_name = single_max.first;
        if (max_component_count_per_all_mts.count(comp_name) == 0) {
          // count for this component was not set
          max_component_count_per_all_mts[comp_name] = single_max.second;
        }
        else {
          // overwrite smaller
          if (max_component_count_per_all_mts[comp_name] < single_max.second) {
            max_component_count_per_all_mts[comp_name] = single_max.second;
          }
        }
      }
    }

    // we finally know how the components will look like, lets create it
    // component_state_names - component name and allowed states
    // max_component_count_per_all_mts - how many components per molecule type are there
    ElemMolType new_mt;
    new_mt.name = same_name_it.first;

    if (is_thrash_or_null(new_mt.name)) {
      continue;
    }

    for (auto comp_info_it: max_component_count_per_all_mts) {
      ComponentType new_ct;
      new_ct.name = comp_info_it.first;
      new_ct.elem_mol_type_name = new_mt.name;

      for (const string& s: component_state_names[new_ct.name]) {
        state_id_t s_id = bng_data->find_or_add_state_name(s);
        new_ct.allowed_state_ids.insert(s_id);
      }

      component_type_id_t ct_id = bng_data->find_or_add_component_type(new_ct);
      for (uint i = 0; i < comp_info_it.second; i++) {
        new_mt.component_type_ids.push_back(ct_id);
      }
    }

    bng_data->find_or_add_elem_mol_type(new_mt);
  }
}


ElemMol SemanticAnalyzer::convert_molecule_pattern(const ASTMolNode* m) {

  ElemMol mi;
  // there is no support for surface molecules in BNGL yet, so everything must be volume molecule
  mi.set_is_vol();

  // process and remember ID
  elem_mol_type_id_t molecule_type_id = bng_data->find_elem_mol_type_id(m->name);
  if (molecule_type_id == MOL_TYPE_ID_INVALID) {
    errs_loc(m) << "Molecule type with name '" + m->name + "' was not defined.\n"; // test N0200
    ctx->inc_error_count();
    return mi;
  }

  const ElemMolType& mt = bng_data->get_elem_mol_type(molecule_type_id);
  mi.elem_mol_type_id = molecule_type_id;

  // process compartment
  compartment_id_t cid = COMPARTMENT_ID_NONE;
  if (m->compartment != nullptr) {
    string compartment_name = m->compartment->str;
    if (compartment_name != "") {
      cid = bng_data->find_compartment_id(compartment_name);
      if (cid == COMPARTMENT_ID_INVALID) {
        errs_loc(m) <<
            "Compartment '" << compartment_name << "' was not defined.\n"; // test XXX
        ctx->inc_error_count();
        return mi;
      }
    }
  }
  mi.compartment_id = cid;
  
  // make a multiset of components type ids so that we can check that
  // out molecule instance does not use wrong or too many components
  multiset<component_type_id_t> remaining_component_ids;
  for (component_type_id_t component_type_id: mt.component_type_ids) {
    remaining_component_ids.insert(component_type_id);
  }


  uint current_component_index = 0;

  // process component instances for this molecule type pattern/instance
  for (size_t i = 0; i < m->components->items.size(); i++) {
    // component_instances
    const ASTComponentNode* component = to_component_node(m->components->items[i]);

    // search in
    component_type_id_t component_type_id =
        bng_data->find_component_type_id(mt, component->name);

    if (component_type_id == COMPONENT_TYPE_ID_INVALID) {
      errs_loc(m) <<
          "Molecule type '" << mt.name << "' does not declare component '" << component->name << "'.\n"; // test N0201
      ctx->inc_error_count();
      return mi;
    }

    // didn't we use the component too many times?
    if (remaining_component_ids.count(component_type_id) == 0) {
      errs_loc(m) <<
          "Molecule type's '" << mt.name << "' component '" << component->name << "' is used too many times.\n"; // test N0102
      ctx->inc_error_count();
      return mi;
    }
    remaining_component_ids.erase(remaining_component_ids.find(component_type_id));

    // add this component to our instance
    mi.components.push_back(Component(component_type_id));

    // state
    if (component->states->items.size() > 1) {
      // checked by parser, original test N0202
      errs_loc(component) <<
          "A component might have max. 1 state specified, error for component " << component->name << ".\n";
      ctx->inc_error_count();
      return mi;
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
        return mi;
      }

      // is this state allowed for this component?
      const ComponentType& ct = bng_data->get_component_type(component_type_id);
      if (ct.allowed_state_ids.count(state_id) == 0) {
        errs_loc(component) <<
            "State name '" << state_name << "' was not declared as allowed for component '" << component->name << "'.\n"; // test N0204
        ctx->inc_error_count();
        return mi;
      }

      // finally set the component's state
      mi.components[current_component_index].state_id = state_id;
    }

    // bond
    const string& bond = component->bond->str;
    bond_value_t b = str_to_bond_value(bond);
    if (b == BOND_VALUE_INVALID) {
      // this is already checked by parser, so the a test checks parser message but let's keep this
      // parser test: N0205
      ctx->internal_error(component->bond, "Invalid bond index");
    }

    mi.components.back().bond_value = b;

    // need to move component index because we processed this component
    current_component_index++;
  }

  return mi;
}


void insert_compartment_id_to_set_based_on_type(
    const BNGData* bng_data,
    const compartment_id_t cid,
    bool& all_are_none_or_inout,
    bool& has_compartment_none,
    uint_set<compartment_id_t>& vol_compartments,
    uint_set<compartment_id_t>& surf_compartments) {

  if (cid == COMPARTMENT_ID_NONE) {
    has_compartment_none = true;
  }
  else if (cid == COMPARTMENT_ID_IN ||
      cid == COMPARTMENT_ID_OUT) {
    // continue
  }
  else if (bng_data->get_compartment(cid).is_3d) {
    vol_compartments.insert(cid);
    all_are_none_or_inout = false;
  }
  else {
    surf_compartments.insert(cid);
    all_are_none_or_inout = false;
  }
}


// for a pattern it is ok to not to list all components
void SemanticAnalyzer::convert_cplx(
    const ASTCplxNode* cplx_node,
    Cplx& bng_cplx,
    const bool in_rule_or_observable,
    const bool check_compartments
) {
  for (const ASTMolNode* m: cplx_node->mols) {

    // molecule ids are based on their name
    bng_cplx.elem_mols.push_back( convert_molecule_pattern(m) );
    if (ctx->get_error_count() != 0) {
      return;
    }
  }

  // semantic checks on bonds validity
  map<bond_value_t, vector<uint> > bond_value_to_molecule_index;
  for (uint i = 0; i < bng_cplx.elem_mols.size(); i++) {
    const ElemMol& mi = bng_cplx.elem_mols[i];

    for (const Component& compi: mi.components) {
      if (compi.bond_has_numeric_value()) {
        // remember molecule index in this complex for a given bond
        bond_value_to_molecule_index[compi.bond_value].push_back(i);
      }
    }
  }

  for (const auto& it: bond_value_to_molecule_index) {
    // each bond is used exactly twice
    if (it.second.size() != 2) {
      assert(cplx_node->size() > 0);
      errs_loc(cplx_node->mols[0]) <<
          "Bond with numerical value '" << it.first << "' must be used exactly twice in a complex pattern of a rule.\n"; // test N0206
      ctx->inc_error_count();
      return;
    }
  }


  // global compartment
  compartment_id_t global_compartment_id = COMPARTMENT_ID_NONE;
  if (cplx_node->compartment != nullptr) {
    string compartment_name = cplx_node->compartment->str;
    if (compartment_name != "") {
      global_compartment_id = bng_data->find_compartment_id(compartment_name);
      if (global_compartment_id == COMPARTMENT_ID_INVALID) {
        errs_loc(cplx_node) <<
            "Compartment '" << compartment_name << "' was not defined.\n"; // tests N0305, N0306
        ctx->inc_error_count();
        return;
      }
    }
  }

  // apply global compartment to elem mols that were not specified
  if (global_compartment_id != COMPARTMENT_ID_NONE) {
    for (auto& em: bng_cplx.elem_mols) {
      if (em.compartment_id == COMPARTMENT_ID_NONE) {
        em.compartment_id = global_compartment_id;
      }
    }
  }

  if (check_compartments) {
    // check that compartments are used consistently
    // we do not know yet whether elementary molecules are surface or not, but
    // compartments were already defined in case we are parsing whole BNGL file
    uint_set<compartment_id_t> vol_compartments;
    uint_set<compartment_id_t> surf_compartments;
    bool all_are_none_or_inout = true;
    bool has_compartment_none = false;
    for (const auto& em: bng_cplx.elem_mols) {
      insert_compartment_id_to_set_based_on_type(
          bng_data, em.compartment_id,
          all_are_none_or_inout, has_compartment_none, vol_compartments, surf_compartments);
    }

    uint_set<compartment_id_t> all_vol_surf_compartment_ids;
    all_vol_surf_compartment_ids.insert(vol_compartments.begin(), vol_compartments.end());
    all_vol_surf_compartment_ids.insert(surf_compartments.begin(), surf_compartments.end());

    if (!all_are_none_or_inout && surf_compartments.empty() && vol_compartments.size() > 1) {
      errs_loc(cplx_node->mols[0]) <<
          "The maximum number of compartments that a volume complex may use is 1, error for '" << bng_cplx.to_str() << "'.\n"; // test N307
      ctx->inc_error_count();
      return;
    }
    if (!all_are_none_or_inout && surf_compartments.size() > 1) {
      errs_loc(cplx_node->mols[0]) <<
          "The maximum number of surface compartments that a surface complex may use is 1, error for '" << bng_cplx.to_str() << "'.\n"; // test XXX
      ctx->inc_error_count();
      return;
    }
    if (!all_are_none_or_inout && surf_compartments.size() == 1 && vol_compartments.size() > 2) {
      errs_loc(cplx_node->mols[0]) <<
          "The maximum number of volume compartments that a surface complex may use is 2, error for '" << bng_cplx.to_str() << "'.\n"; // test XXX
      ctx->inc_error_count();
      return;
    }


    bng_cplx.finalize_cplx();
    if (!bng_cplx.is_connected() && !in_rule_or_observable) {
      errs_loc(cplx_node->mols[0]) <<
          "All complexes that are not patterns must be currently fully connected, error for '" << bng_cplx.to_str() << "'.\n"; // test XXX
      ctx->inc_error_count();
      return;
    }
  }
}


// take one side of a reaction rule and create pattern for rule matching
void SemanticAnalyzer::convert_rxn_rule_side(
    const ASTListNode* rule_side,
    const bool reactants_side,
    CplxVector& patterns) {

  // we need to check each molecule type from each complex
  std::vector<const ASTMolNode*> current_complex_nodes;
  for (size_t i = 0; i < rule_side->items.size(); i++) {
    const ASTCplxNode* cplx = to_cplx_node(rule_side->items[i]);

    if (cplx->size() == 1) {
      const ASTMolNode* m = to_molecule_node(cplx->mols[0]);
      if (is_thrash_or_null(m->name)) {
        if (reactants_side) {
          errs_loc(m) <<
              "Null/Trash product cannot be used on the reactants side of a reeaction rule.\n"; // test N0620
          ctx->inc_error_count();
          return;
        }
        else {
          // ok, ignore
          continue;
        }
      }
    }

    Cplx pattern(bng_data);
    convert_cplx(cplx, pattern, true);
    if (ctx->get_error_count() > 0) {
      return;
    }
    pattern.finalize_cplx();
    patterns.push_back(pattern);
  }
}


void SemanticAnalyzer::finalize_and_store_rxn_rule(
    const ASTRxnRuleNode* n, RxnRule& r, const bool forward_direction) {

  // check in/out
  if (r.reactants.size() == 2) {
    if (r.reactants[0].has_compartment_class_in_out() && r.reactants[1].has_compartment_class_in_out()) {
      errs_loc(n) << "Maximum one reactant may use @" << COMPARTMENT_NAME_IN << " or @" << COMPARTMENT_NAME_OUT <<
          " compartment class. Reported for reaction in the " <<
          (forward_direction ? DIR_FORWARD : DIR_REVERSE) << " direction.\n"; // TEST 600
      ctx->inc_error_count();
      return;
    }
  }

  r.finalize();

  // determine mapping from molecule instances on one side to another
  stringstream out;
  bool ok = r.check_reactants_products_mapping(out);
  if (!ok) {
    // tests N0220, N0230, N0231, N0232
    errs_loc(n) << out.str() <<
        " Reported for reaction in the " << (forward_direction ? DIR_FORWARD : DIR_REVERSE) << " direction.\n";
    ctx->inc_error_count();
    return;
  }

  bng_data->find_or_add_rxn_rule(r);
}


void SemanticAnalyzer::convert_and_store_rxn_rules() {

  // for each reaction rule
  for (const ASTBaseNode* n: ctx->rxn_rules.items) {
    const ASTRxnRuleNode* r = to_rxn_rule_node(n);

    CplxVector reactants;
    convert_rxn_rule_side(r->reactants, true, reactants);

    CplxVector products;
    convert_rxn_rule_side(r->products, false, products);

    if (ctx->get_error_count() > 0) {
      return;
    }

    RxnRule fwd_rule(bng_data);
    fwd_rule.type = RxnType::Standard;
    fwd_rule.name = r->name;
    assert(r->rates->items.size() >= 1);
    fwd_rule.base_rate_constant = to_expr_node(r->rates->items[0])->get_dbl();
    fwd_rule.reactants = reactants;
    fwd_rule.products = products;

    finalize_and_store_rxn_rule(r, fwd_rule, true);

    if (r->reversible) {
      RxnRule rev_rule(bng_data);
      rev_rule.type = RxnType::Standard;
      rev_rule.name = r->name;

      if (products.empty()) {
        errs_loc(r) <<
            "A reversible rule must have at least one complex pattern " <<
            "on the right side of the reaction rule.\n"; // caught by parser, test N0211
        ctx->inc_error_count();
      }

      assert(r->rates->items.size() == 2);
      rev_rule.base_rate_constant = to_expr_node(r->rates->items[1])->get_dbl();
      rev_rule.reactants = products;
      rev_rule.products = reactants;

      finalize_and_store_rxn_rule(r, rev_rule, false);
    }
  }
}


void SemanticAnalyzer::convert_seed_species() {
  for (const ASTBaseNode* n: ctx->seed_species.items) {
    const ASTSeedSpeciesNode* ss_node = to_seed_species_node(n);

    SeedSpecies ss(bng_data);

    convert_cplx(ss_node->cplx, ss.cplx, false);
    if (ctx->get_error_count() != 0) {
      return;
    }

    if (ss.cplx.has_compartment_class_in_out()) {
      errs_loc(ss_node->cplx) <<
          "It is not allowed to use compartment @" << compartment_id_to_str(ss.cplx.get_primary_compartment_id()) <<
          " in the seed species section.\n"; // test N0601
      ctx->inc_error_count();
      return;
    }

    uint_set<compartment_id_t> used_compartments;
    ss.cplx.get_used_compartments(used_compartments);
    if (used_compartments.count(COMPARTMENT_ID_NONE) && used_compartments.size() > 1) {
      errs_loc(n) <<
          "In the seed species section either all elementary molecules must have their compartment specified or none, error for '" <<
          ss.cplx.to_str() << "'.\n"; // test N0621
      ctx->inc_error_count();
      return;
    }


    ASTExprNode* orig_expr = to_expr_node(ss_node->count);
    ASTExprNode* new_expr = evaluate_to_dbl(orig_expr);
    ss.count = new_expr->get_dbl();

    bng_data->add_seed_species(ss);
  }
}


void SemanticAnalyzer::convert_observables() {
  for (const ASTBaseNode* n: ctx->observables.items) {
    const ASTObservableNode* o_node = to_observable_node(n);

    Observable o;

    if (o_node->type == OBSERVABLE_MOLECULES) {
      o.type = ObservableType::Molecules;
    }
    else if (o_node->type == OBSERVABLE_SPECIES) {
      o.type = ObservableType::Species;
    }
    else {
      errs_loc(o_node) <<
          "Invalid observable type '" << o_node->type << "', the allowed values are" <<
          OBSERVABLE_MOLECULES << " or " << OBSERVABLE_SPECIES << ".\n"; // TODO test
      ctx->inc_error_count();
      continue;
    }

    o.name = o_node->name;

    for (const ASTBaseNode* base_pat: o_node->cplx_patterns->items) {
      const ASTCplxNode* cplx_pat = to_cplx_node(base_pat);
      Cplx cplx(bng_data);
      convert_cplx(cplx_pat, cplx, true);
      if (ctx->get_error_count() != 0) {
        return;
      }

      if (cplx.has_compartment_class_in_out()) {
        errs_loc(cplx_pat) <<
            "It is not allowed to use compartment @" << compartment_id_to_str(cplx.get_primary_compartment_id()) <<
            " in the observables section.\n"; // test N0602
        ctx->inc_error_count();
        return;
      }

      o.patterns.push_back(cplx);
    }

    bng_data->add_observable(o);
  }
}


// returns true if conversion and semantic checks passed
bool SemanticAnalyzer::check_and_convert_parsed_file(
    ParserContext* ctx_,
    BNGData* res_bng,
    const std::map<std::string, float_t>& parameter_overrides) {

  assert(ctx_ != nullptr);
  assert(res_bng != nullptr);

  ctx = ctx_;
  bng_data = res_bng;

  // the following conversions do not use the bng_data.parameters map
  convert_and_evaluate_parameters(parameter_overrides);
  if (ctx->get_error_count() != 0) {
    return false;
  }

  // first compute all reaction rates
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

  convert_and_store_compartments();
  if (ctx->get_error_count() != 0) {
    return false;
  }

  // molecule types do not have to be defined,
  // in this case we will define them based on what we found in the
  // seed species, reactions (and observables - not supported yet)
  collect_and_store_implicit_molecule_types();
  if (ctx->get_error_count() != 0) {
    return false;
  }

  // convert rxn rules
  convert_and_store_rxn_rules();
  if (ctx->get_error_count() != 0) {
    return false;
  }

  convert_seed_species();
  if (ctx->get_error_count() != 0) {
    return false;
  }

  convert_observables();
  if (ctx->get_error_count() != 0) {
    return false;
  }

  return true;
}


// analyze AST for a single complex and extend molecule type definitions by what it contains
void SemanticAnalyzer::extend_molecule_type_definitions(const ASTCplxNode* cplx_node) {

  for (const ASTMolNode* mol: cplx_node->mols) {
    // was this molecule type already defined?
    elem_mol_type_id_t mt_id = bng_data->find_elem_mol_type_id(mol->name);
    if (mt_id == MOL_TYPE_ID_INVALID) {
      // no, simply add as a new one
      ElemMolType mt = convert_molecule_type(mol, true);
      bng_data->find_or_add_elem_mol_type(mt);
    }
    else {
      // yes - extend existing components definitions
      ElemMolType& mt = bng_data->get_elem_mol_type(mt_id);

      // for each component of the new cplx
      for (const ASTBaseNode* c: mol->components->items) {

        const ASTComponentNode* comp_node = to_component_node(c);

        // find whether a component with the same name is already present
        // and extend its allowed states
        ComponentType* ct = nullptr;
        for (component_type_id_t existing_ct_id: mt.component_type_ids) {
          // we must limit ourselves to the components already defined for this molecule
          ComponentType* ct_existing = &bng_data->get_component_type(existing_ct_id);
          if (comp_node->name == ct_existing->name) {
            ct = ct_existing;
            break;
          }
        }
        if (ct != nullptr) {
          // found existing component, extend its allowed states
          for (const ASTBaseNode* state: comp_node->states->items) {
            state_id_t state_id = convert_state_name(to_str_node(state));
            ct->allowed_state_ids.insert(state_id);
          }
        }
        else {
          // there is no such component, this may happen in cases such as:
          // CaMKII(l!1,Y286~0,cam!2).CaM(C~0,N~0,camkii!2).CaMKII(r!1,Y286~P)
          // we must add the component to the molecule type
          ComponentType ct = convert_component_type(mt.name, to_component_node(c), true);
          component_type_id_t new_ct_id = bng_data->find_or_add_component_type(ct);
          mt.component_type_ids.push_back(new_ct_id);
        }
      }
    }
  }
}


void SemanticAnalyzer::define_compartments_used_by_cplx_as_3d_compartments(
    ASTCplxNode* cplx_node) {

  // to be used only for single cplx mode
  set<string> compartment_names;
  if (cplx_node->compartment != nullptr && cplx_node->compartment->str != "") {
    compartment_names.insert(cplx_node->compartment->str);
  }
  for (ASTMolNode* mol: cplx_node->mols) {
    if (mol->compartment != nullptr && mol->compartment->str != "") {
      compartment_names.insert(mol->compartment->str);
    }
  }

  for (const string& n: compartment_names) {
    compartment_id_t in_out_id = get_in_or_out_compartment_id(n);
    if (in_out_id == COMPARTMENT_ID_INVALID) {
      // is a standard compartment, add it
      Compartment c;
      c.name = n;
      c.is_3d = true;
      bng_data->add_compartment(c);
    }
  }
}


// returns true if conversion and semantic checks passed,
// resulting complex is stored into res
// called only from parse_single_cplx_string
bool SemanticAnalyzer::check_and_convert_single_cplx(
    ParserContext* ctx_, BNGData* res_bng, Cplx& res) {
  assert(ctx_ != nullptr);
  assert(res_bng != nullptr);

  ctx = ctx_;
  bng_data = res_bng;

  // define molecule types first
  extend_molecule_type_definitions(ctx->single_cplx);
  if (ctx->get_error_count() != 0) {
    return false;
  }

  define_compartments_used_by_cplx_as_3d_compartments(ctx->single_cplx);

  convert_cplx(ctx->single_cplx, res, true, false);
  if (ctx->get_error_count() != 0) {
    return false;
  }

  return true;
}

} /* namespace BNG */
