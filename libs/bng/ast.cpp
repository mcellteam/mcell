#include <iostream>

#include "ast.h"
#include "parser_utils.h"

using namespace std;

namespace BNG {

// ------------------------------- ASTBaseNode ------------------------
void ASTBaseNode::dump(const std::string ind) {
  // empty for now
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

// ------------------------------- ASTBaseNode ------------------------
void ASTListNode::dump(const std::string ind) {
  // empty
}

// ------------------------------- ASTSymbolTable ------------------------
void ASTSymbolTable::insert(const std::string id, ASTBaseNode* node) {
  if (table.count(id) != 0) {
    errs() << "Symbol '" << id << "' was already defined.\n";
  }

  table[id] = node;
}

void ASTSymbolTable::dump() {
  cout << "ASTSymbolTable:\n";
  for (const auto& item: table) {
    cout << "  " << item.first << " = \n";
    item.second->dump("    ");
  }
}

// ------------------------------- ASTContext ----------------------------
ASTContext::~ASTContext() {

  for (ASTBaseNode* n: all_nodes) {
    delete n;
  }
  all_nodes.clear();
}

ASTExprNode* ASTContext::new_id_node(const std::string id) {
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

void ASTContext::dump() {
  cout << "-- ASTContext dump --\n";
  symtab.dump();
}

} // namespace BNG
