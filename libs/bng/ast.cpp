#include <cstdlib>
#include <cerrno>
#include <iostream>
#include <string>

#include "bng/ast.h"
#include "bng/parser_utils.h"
#include "bng/bngl_names.h"

using namespace std;

namespace BNG {

static const std::string IND2 = "  ";
static const std::string IND4 = "    ";

// ------------------------------- ASTBaseNode ------------------------
void ASTBaseNode::dump(const std::string ind) const {
  if (has_loc) {
    // not printing file names because they are the same for now
    cout << ind << "line: " << line << "\n";
  }
}


// ------------------------------- ASTExprNode ------------------------
void ASTExprNode::dump(const std::string ind) const {
  cout << ind;
  switch (expr_type) {
    case ExprType::Id:
      cout << "Id: " << id;
      break;
    case ExprType::Dbl:
      cout << "Dbl: " << dbl;
      break;
    case ExprType::Llong:
      cout << "Llong: " << llong;
      break;
    default:
      cout << "expression";
  }

  cout << "\n";
  ASTBaseNode::dump(ind);
}


// ------------------------------- ASTStrNode ------------------------
void ASTStrNode::dump(const std::string ind) const {
  cout << ind << "str: '" << str << "'";
  ASTBaseNode::dump(ind);
}


// ------------------------------- ASTListNode ------------------------
void ASTListNode::dump(const std::string ind) const {
  if (items.empty()) {
    cout << ind << "(empty)\n";
  }
  else {
    for (size_t i = 0; i < items.size(); i++) {
      assert(items[i] != nullptr);
      cout << ind << i << ": \n";
      items[i]->dump(ind + IND2);
    }
    cout << "\n";
  }
  ASTBaseNode::dump(ind);
}


// ------------------------------- ASTComponentNode ------------------------
void ASTComponentNode::dump(const std::string ind) const {
  cout << ind << "component: name='" << name << "'\n";
  cout << ind << "  states:\n";
  assert(states != nullptr);
  states->dump(ind + IND4);
  assert(bond != nullptr);
  cout << ind << "  bond:\n";
  bond->dump(ind + IND4);
  cout << "\n";
  ASTBaseNode::dump(ind);
}


// ------------------------------- ASTMoleculeNode ------------------------
void ASTMolNode::dump(const std::string ind) const {
  cout << ind << "molecule: name='" << name << "'\n";
  cout << ind << "  components:\n";
  assert(components != nullptr);
  components->dump(ind + IND4);
  if (compartment != nullptr) {
    cout << ind << "  compartment:\n";
    compartment->dump(ind + IND4);
  }
  ASTBaseNode::dump(ind);
}


// ------------------------------- ASTMoleculeNode ------------------------
void ASTCompartmentNode::dump(const std::string ind) const {
  cout << ind <<
      "compartment: name='" << name <<
      "', dimensions " << dimensions <<
      ", volume (expression, not evaluated)"
      ", parent: '" << parent_name << "'\n";

  ASTBaseNode::dump(ind);
}


// ------------------------------- ASTCplxNode ------------------------
void ASTCplxNode::dump(const std::string ind) const {
  if (compartment != nullptr) {
    cout << ind << "  compartment:\n";
    compartment->dump(ind + IND4);
  }
  if (mols.empty()) {
    cout << ind << "(empty)\n";
  }
  else {
    for (size_t i = 0; i < mols.size(); i++) {
      assert(mols[i] != nullptr);
      cout << ind << i << ": \n";
      mols[i]->dump(ind + IND2);
    }
    cout << "\n";
  }
  ASTBaseNode::dump(ind);
}


// ------------------------------- ASTRxnRuleNode ------------------------
void ASTRxnRuleNode::dump(const std::string ind) const {
  cout << ind << "reaction rule: name='" << name << "', reversible: " << (reversible?"true":"false") << "\n";
  cout << ind << "  reactants:\n";
  assert(reactants != nullptr);
  reactants->dump(ind + IND4);
  cout << ind << "  products:\n";
  assert(products != nullptr);
  products->dump(ind + IND4);
  cout << ind << "  rates:\n";
  assert(rates != nullptr);
  rates->dump(ind + IND4);
  ASTBaseNode::dump(ind);
}


// ------------------------------- ASTSeedSpeciesNode ------------------------
void ASTSeedSpeciesNode::dump(const std::string ind) const {
  cout << ind << "seed species item:\n";
  cout << ind << "  cplx:\n";
  assert(cplx != nullptr);
  cplx->dump(ind + IND4);
  cout << ind << "  count:\n";
  assert(count != nullptr);
  count->dump(ind + IND4);
  ASTBaseNode::dump(ind);
}


// ------------------------------- ASTObservableNode ------------------------
void ASTObservableNode::dump(const std::string ind) const {
  cout << ind << "observable item:\n";
  cout << ind << "  type:" << type << "\n";
  cout << ind << "  name:" << name << "\n";
  cout << ind << "  cplx_patterns:\n";
  for (ASTBaseNode* n: cplx_patterns->items) {
    const ASTCplxNode* cplx = to_cplx_node(n);
    cplx->dump(ind + IND4);
  }
  ASTBaseNode::dump(ind);
}


// ------------------------------- ASTSymbolTable ------------------------
void ASTSymbolTable::insert(const std::string id, ASTBaseNode* node, ParserContext* ctx) {
  if (table.count(id) != 0) {
    errs_loc() << "Symbol '" << id << "' was already defined.\n";
    ctx->inc_error_count();
  }

  // special handling of MCELL_REDEFINE_ symbols
  if (id.find(MCELL_REDEFINE_PREFIX) == 0) {
    // if symbol starts with the redefine prefix, overwrite it
    string id_to_redefine = id.substr(strlen(MCELL_REDEFINE_PREFIX));
    table[id_to_redefine] = node;
  }
  else {
    table[id] = node;
  }
}


ASTBaseNode* ASTSymbolTable::get(const std::string& id, ASTBaseNode* loc, ParserContext* ctx) const {
  auto it = table.find(id);
  if (it == table.end()) {
    errs_loc() << "Symbol '" << id << "' is not defined.\n";
    ctx->inc_error_count();
    return nullptr;
  }
  else {
    return it->second;
  }
}


void ASTSymbolTable::insert_molecule_declarations(const ASTListNode* molecule_node_list, ParserContext* ctx) {
  for (ASTBaseNode* n: molecule_node_list->items) {
    assert(n->is_mol());
    ASTMolNode* mn = to_molecule_node(n);
    insert(mn->name, mn, ctx);
  }
}


void ASTSymbolTable::dump() {
  cout << "ASTSymbolTable:\n";
  for (const auto& item: table) {
    cout << IND2 << item.first << " = \n";
    item.second->dump(IND4);
  }
}


// ------------------------------- ASTContext ----------------------------
ParserContext::~ParserContext() {

  for (ASTBaseNode* n: all_nodes) {
    delete n;
  }
  all_nodes.clear();
}


ASTExprNode* ParserContext::new_id_node(const std::string& id, const BNGLLTYPE& loc) {
  ASTExprNode* n = new ASTExprNode;
  n->set_id(id);
  n->set_loc(current_file, loc);
  remember_node(n);
  return n;
}


ASTExprNode* ParserContext::new_dbl_node(const double val, const BNGLLTYPE& loc) {
  ASTExprNode* n = new ASTExprNode;
  n->set_dbl(val);
  n->set_loc(current_file, loc);
  remember_node(n);
  return n;
}


ASTExprNode* ParserContext::new_dbl_node(const double val, const ASTBaseNode* loc) {
  ASTExprNode* n = new ASTExprNode;
  n->set_dbl(val);
  n->set_loc(loc->file, loc->line);
  remember_node(n);
  return n;
}


ASTExprNode* ParserContext::new_dbl_node(const double val) {
  ASTExprNode* n = new ASTExprNode;
  n->set_dbl(val);
  n->has_loc = false;
  remember_node(n);
  return n;
}

ASTExprNode* ParserContext::new_llong_node(const long long val, const BNGLLTYPE& loc) {
  ASTExprNode* n = new ASTExprNode;
  n->set_llong(val);
  n->set_loc(current_file, loc);
  remember_node(n);
  return n;
}


ASTExprNode* ParserContext::new_expr_node(
    ASTExprNode* left, const ExprType op, ASTExprNode* right, const BNGLLTYPE& loc) {
  assert((op == ExprType::UnaryPlus || op == ExprType::UnaryMinus) == (right == nullptr));
  ASTExprNode* n = new ASTExprNode;
  n->set_left(left);
  n->set_type(op); // operator
  n->set_right(right);
  n->set_loc(current_file, loc);
  remember_node(n);
  return n;
}


ASTExprNode* ParserContext::new_expr_node(
    const std::string& function_name,
    ASTListNode* arguments,
    const BNGLLTYPE& loc) {

  assert(function_name != "");
  ASTExprNode* n = new ASTExprNode;
  n->set_type(ExprType::FunctionCall);
  n->set_function_name(function_name);
  n->set_function_arguments(arguments);
  n->set_loc(current_file, loc);
  remember_node(n);
  return n;
}


ASTStrNode* ParserContext::new_empty_str_node() {
  ASTStrNode* n = new ASTStrNode;
  n->str = "";
  remember_node(n);
  return n;
}


ASTStrNode* ParserContext::new_str_node(const std::string str, const BNGLLTYPE& loc) {
  ASTStrNode* n = new ASTStrNode;
  n->str = str;
  n->set_loc(current_file, loc);
  remember_node(n);
  return n;
}


ASTStrNode* ParserContext::new_str_node(const long long val_to_str, const BNGLLTYPE& loc) {
  ASTStrNode* n = new ASTStrNode;
  n->str = to_string(val_to_str);
  n->set_loc(current_file, loc);
  remember_node(n);
  return n;
}


ASTListNode* ParserContext::new_list_node() {
  ASTListNode* n = new ASTListNode();
  remember_node(n);
  return n;
}


ASTComponentNode* ParserContext::new_component_node(
    const std::string& name,
    ASTListNode* state_list,
    ASTStrNode* bond,
    const BNGLLTYPE& loc
) {
  ASTComponentNode* n = new ASTComponentNode();
  n->name = name;
  n->states = state_list;
  n->bond = bond;
  n->set_loc(current_file, loc);
  remember_node(n);
  return n;
}


ASTMolNode* ParserContext::new_molecule_node(
    const std::string& name,
    ASTListNode* component_list,
    ASTStrNode* compartment,
    const BNGLLTYPE& loc
) {
  ASTMolNode* n = new ASTMolNode();
  n->name = name;
  n->components = component_list;
  n->compartment = compartment;
  n->set_loc(current_file, loc);
  remember_node(n);
  return n;
}


ASTCompartmentNode* ParserContext::new_compartment_node(
    const std::string& name,
    const int dimensions,
    ASTExprNode* volume,
    const std::string parent_name,
    const BNGLLTYPE& loc
) {
  ASTCompartmentNode* n = new ASTCompartmentNode();
  n->name = name;
  n->dimensions = dimensions;
  n->volume = volume;
  n->parent_name = parent_name;
  n->set_loc(current_file, loc);
  remember_node(n);
  return n;
}


ASTCplxNode* ParserContext::new_cplx_node(ASTMolNode* first_mol) {
  ASTCplxNode* n = new ASTCplxNode(first_mol);
  n->set_loc(first_mol);

  remember_node(n);
  return n;
}


ASTRxnRuleNode* ParserContext::new_rxn_rule_node(
    ASTStrNode* name,
    ASTListNode* reactants,
    const bool reversible,
    ASTListNode* products,
    ASTListNode* rates
) {
  ASTRxnRuleNode* n = new ASTRxnRuleNode();
  n->name = name->str;
  n->reactants = reactants;
  n->reversible = reversible;
  n->products = products;
  n->rates = rates;

  // use the first reactant as the location
  assert(reactants->items.size() >= 1);
  ASTCplxNode* cplx = to_cplx_node(reactants->items[0]);
  n->set_loc(cplx->mols[0]);

  remember_node(n);
  return n;
}


ASTSeedSpeciesNode* ParserContext::new_seed_species_node(
    ASTCplxNode* cplx,
    ASTExprNode* count
) {
  ASTSeedSpeciesNode* n = new ASTSeedSpeciesNode();
  n->cplx = cplx;
  n->count = count;

  // use the first molecule of the complex as the location
  n->set_loc(cplx);

  remember_node(n);
  return n;
}


ASTObservableNode* ParserContext::new_observable_node(
    const std::string& type,
    const std::string& name,
    ASTListNode* cplx_patterns,
    const BNGLLTYPE& loc
) {
  ASTObservableNode* n = new ASTObservableNode();
  n->type = type;
  n->name = name;
  n->cplx_patterns = cplx_patterns;
  n->set_loc(current_file, loc);

  remember_node(n);
  return n;
}


void ParserContext::print_error_report() {
  if (errors != 0) {
    cerr << "Compilation failed, there were " << errors << " errors.\n";
  }
}


void ParserContext::internal_error(const ASTBaseNode* loc, const std::string msg) {
  errs_loc(loc) << "INTERNAL: " << msg;
  exit(2);
}


void ParserContext::dump() {
  cout << "-- ASTContext dump --\n";
  symtab.dump();
  cout << "compartments:\n";
  compartments.dump(IND2);
  cout << "reaction rules:\n";
  rxn_rules.dump(IND2);
  cout << "seed species:\n";
  seed_species.dump(IND2);
  cout << "observables:\n";
  observables.dump(IND2);
}


bond_value_t str_to_bond_value(const std::string& s) {
  if (s == "") {
    return BOND_VALUE_UNBOUND;
  }
  else if (s == BOND_STR_ANY) {
    return BOND_VALUE_ANY; // !?
  }
  else if (s == BOND_STR_BOUND) {
    return BOND_VALUE_BOUND; // !+
  }
  else {
    // try to convert

    char* end;
    long long int res;

    errno = 0; // note: errno is thread-local
    res = strtoll(s.c_str(), &end, 10);

    // conversion error
    if (errno != 0 || *end != '\0') {
      return BOND_VALUE_INVALID;
    }

    // range check
    if (res < 0 || res >= BOND_VALUE_UNBOUND) {
      return BOND_VALUE_INVALID;
    }

     // ok
    return res;
  }
}

} // namespace BNG
