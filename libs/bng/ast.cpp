#include <iostream>

#include "ast.h"
#include "parser_utils.h"

using namespace std;

namespace BNG {

static const std::string IND2 = "  ";

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

// ------------------------------- ASTListNode ------------------------
void ASTListNode::dump(const std::string ind) {
  if (items.empty()) {
    cout << ind << "(empty)\n";
  }
  else {
    for (size_t i = 0; i < items.size(); i++) {
      assert(items[i] != nullptr);
      cout << ind << i << ": ";
      items[i]->dump(ind + IND2);
    }
    cout << "\n";
  }
}

// ------------------------------- ASTComponentNode ------------------------
void ASTComponentNode::dump(const std::string ind) {
  cout << "component: name='" << name << "', states:\n";
  assert(state_list != nullptr);
  state_list->dump(ind + IND2);
  assert(bond != nullptr);
  cout << ind << " " << "bond: ";
  bond->dump(ind);
  cout << "\n";
}

// ------------------------------- ASTStrNode ------------------------
void ASTStrNode::dump(const std::string ind) {
  cout << "str: '" << str << "'";
}

// ------------------------------- ASTMoleculeNode ------------------------
void ASTMoleculeNode::dump(const std::string ind) {
  cout << ind << "molecule: name='" << name << "', components:\n";
  assert(component_list != nullptr);
  component_list->dump(ind + IND2);
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
    item.second->dump(IND2 + IND2);
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

ASTComponentNode* ASTContext::new_component_node(
    const std::string& name,
    ASTListNode* state_list,
    ASTStrNode* bond
) {
  ASTComponentNode* n = new ASTComponentNode();
  n->name = name;
  n->state_list = state_list;
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
  n->component_list = component_list;
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
}

} // namespace BNG
