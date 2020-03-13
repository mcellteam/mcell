#include <iostream>

#include "ast.h"
#include "parser_utils.h"

using namespace std;

namespace BNG {

static const std::string IND2 = "  ";
static const std::string IND4 = "    ";

// ------------------------------- ASTBaseNode ------------------------
void ASTBaseNode::dump(const std::string ind) {
  // empty for now, should be used to print file and line information
}

// ------------------------------- ASTExprNode ------------------------
void ASTExprNode::dump(const std::string ind) {
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
      assert(false);
  }

  cout << "\n";
}


// ------------------------------- ASTStrNode ------------------------
void ASTStrNode::dump(const std::string ind) {
  cout << ind << "str: '" << str << "'";
}



// ------------------------------- ASTStrNode ------------------------
void ASTSeparatorNode::dump(const std::string ind) {
  string s;
  switch (separator_type) {
    case SeparatorType::Dot:
      s = ".";
      break;
    case SeparatorType::Plus:
      s = "+";
      break;
    default:
      assert(false);
      s = "error";
  }
  cout << ind << "separator: '" << s << "'\n";
}


// ------------------------------- ASTListNode ------------------------
void ASTListNode::dump(const std::string ind) {
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
}

// ------------------------------- ASTComponentNode ------------------------
void ASTComponentNode::dump(const std::string ind) {
  cout << ind << "component: name='" << name << "'\n";
  cout << ind << "  states:\n";
  assert(states != nullptr);
  states->dump(ind + IND4);
  assert(bond != nullptr);
  cout << ind << "  bond:\n";
  bond->dump(ind + IND4);
  cout << "\n";
}

// ------------------------------- ASTMoleculeNode ------------------------
void ASTMoleculeNode::dump(const std::string ind) {
  cout << ind << "molecule: name='" << name << "'\n";
  cout << ind << "  components:\n";
  assert(components != nullptr);
  components->dump(ind + IND4);
}

// ------------------------------- ASTMoleculeNode ------------------------
void ASTRxnRuleNode::dump(const std::string ind) {
  cout << ind << "reaction rule: name='" << name << "', reversible: " << (reversible?"true":"false") << "\n";
  cout << ind << "  reactants:\n";
  assert(reactants != nullptr);
  reactants->dump(ind + IND4);
  cout << ind << "  products:\n";
  assert(reactants != nullptr);
  reactants->dump(ind + IND4);
  cout << ind << "  rates:\n";
  assert(rates != nullptr);
  rates->dump(ind + IND4);
}

// ------------------------------- ASTSymbolTable ------------------------
void ASTSymbolTable::insert(const std::string id, ASTBaseNode* node, ASTContext* ctx) {
  if (table.count(id) != 0) {
    errs() << "Symbol '" << id << "' was already defined.\n";
    ctx->inc_error_count();
  }

  table[id] = node;
}

void ASTSymbolTable::insert_molecule_declarations(const ASTListNode* molecule_node_list, ASTContext* ctx) {
  for (ASTBaseNode* n: molecule_node_list->items) {
    assert(n->node_type == NodeType::Molecule);
    ASTMoleculeNode* mn = dynamic_cast<ASTMoleculeNode*>(n);
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
ASTContext::~ASTContext() {

  for (ASTBaseNode* n: all_nodes) {
    delete n;
  }
  all_nodes.clear();
}

ASTExprNode* ASTContext::new_id_node(const std::string& id) {
  ASTExprNode* n = new ASTExprNode;
  n->set_id(id);
  remember_node(n);
  return n;
}

ASTExprNode* ASTContext::new_dbl_node(const double val) {
  ASTExprNode* n = new ASTExprNode;
  n->set_dbl(val);
  remember_node(n);
  return n;
}

ASTExprNode* ASTContext::new_llong_node(const long long val) {
  ASTExprNode* n = new ASTExprNode;
  n->set_llong(val);
  remember_node(n);
  return n;
}

ASTStrNode* ASTContext::new_str_node(const std::string str) {
  ASTStrNode* n = new ASTStrNode;
  n->str = str;
  remember_node(n);
  return n;
}

ASTStrNode* ASTContext::new_str_node(const long long val_to_str) {
  ASTStrNode* n = new ASTStrNode;
  n->str = to_string(val_to_str);
  remember_node(n);
  return n;
}

ASTListNode* ASTContext::new_list_node() {
  ASTListNode* n = new ASTListNode();
  remember_node(n);
  return n;
}

ASTSeparatorNode* ASTContext::new_separator_node(const SeparatorType type) {
  ASTSeparatorNode* n = new ASTSeparatorNode();
  n->separator_type = type;
  remember_node(n);
  return n;
}

ASTComponentNode* ASTContext::new_component_node(
    const std::string& name,
    ASTListNode* state_list,
    ASTStrNode* bond
) {
  ASTComponentNode* n = new ASTComponentNode();
  n->name = name;
  n->states = state_list;
  n->bond = bond;
  remember_node(n);
  return n;
}

ASTMoleculeNode* ASTContext::new_molecule_node(
    const std::string& name,
    ASTListNode* component_list
) {
  ASTMoleculeNode* n = new ASTMoleculeNode();
  n->name = name;
  n->components = component_list;
  remember_node(n);
  return n;
}

ASTRxnRuleNode* ASTContext::new_reaction_rule_node(
    ASTListNode* products,
    const bool reversible,
    ASTListNode* reactants,
    ASTListNode* rates
) {
  ASTRxnRuleNode* n = new ASTRxnRuleNode();
  n->products = products;
  n->reversible = reversible;
  n->reactants = reactants;
  n->rates = rates;
  remember_node(n);
  return n;
}

void ASTContext::print_error_report() {
  if (errors != 0) {
    cerr << "Compilation failed, there were " << errors << " errors.";
  }
}

void ASTContext::dump() {
  cout << "-- ASTContext dump --\n";
  symtab.dump();
  cout << "reaction rules:\n";
  reaction_rules.dump(IND2);
}

} // namespace BNG
