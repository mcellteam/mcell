#ifndef LIBS_BNG_AST_H_
#define LIBS_BNG_AST_H_

/**
 * This file defines abstract syntax tree for BNGL.
 * Compilation must be done in two steps because
 * BNGL does not require definition before use
 * (at least in expressions).
 */

#include <string>

#include "bng_defines.h"

namespace BNG {

class ASTContext;

enum class NodeType {
  Invalid,
  Expr,
  List,
  Molecule,
  Component,
  Str
};

enum class ExprType {
  Invalid,
  Id,
  Dbl,
  Llong
};

static const std::string ANY_BOND = "+";

class ASTBaseNode {
public:
  ASTBaseNode()
    : node_type(NodeType::Invalid), line(0), file(nullptr) {
  }
  virtual ~ASTBaseNode() {}
  virtual void dump(const std::string ind);

  NodeType node_type;

  // TODO
  int line;
  const char* file; // owned by ASTContext
};


class ASTExprNode: public ASTBaseNode {
public:
  ASTExprNode()
    : expr_type(ExprType::Invalid),
      left(nullptr), right(nullptr), dbl(0), llong(0) {
    node_type = NodeType::Expr;
  }
  void dump(const std::string ind) override;

  ExprType expr_type;

  void set_id(const std::string id_) {
    expr_type = ExprType::Id;
    id = id_;
  }

  void set_dbl(const double dbl_) {
    expr_type = ExprType::Dbl;
    dbl = dbl_;
  }

  void set_llong(const long long llong_) {
    expr_type = ExprType::Llong;
    llong = llong_;
  }

  // getters that check the type
private:
  ASTExprNode* left;
  ASTExprNode* right;
  std::string id;
  double dbl;
  long long llong;
};


// general node to store strings
class ASTStrNode: public ASTBaseNode {
public:
  ASTStrNode() {
    node_type = NodeType::Str;
  }
  void dump(const std::string ind) override;

  std::string str;
};


class ASTListNode: public ASTBaseNode {
public:
  ASTListNode() {}
  void dump(const std::string ind) override;

  ASTListNode* append(ASTBaseNode* n) {
    assert(n != nullptr);
    items.push_back(n);
    return this;
  }

  std::vector<ASTBaseNode*> items;
};


class ASTComponentNode: public ASTBaseNode {
public:
  ASTComponentNode()
    : state_list(nullptr), bond(nullptr) {
    node_type = NodeType::Component;
  }
  void dump(const std::string ind) override;

  std::string name;
  ASTListNode* state_list;
  ASTStrNode* bond;
};


// specific instance of a molecule with states
// when defined in molecule types section, it is inserted
// into the symbol table
// when used in reaction rule, it is referenced from the
// reaction itself
class ASTMoleculeNode: public ASTBaseNode {
public:
  ASTMoleculeNode()
    : component_list(nullptr) {
    node_type = NodeType::Molecule;
  }
  void dump(const std::string ind) override;

  std::string name;
  ASTListNode* component_list;
};


class ASTSymbolTable {
public:
  void insert(const std::string id, ASTBaseNode* node, ASTContext* ctx);
  void insert_molecule_declarations(const ASTListNode* molecule_node_list, ASTContext* ctx);
  ASTBaseNode* get(const std::string& id) const;

  void dump();
private:
  // we do not hold constant pointers because one for instance might
  // want to evaluate an expression
  std::map<std::string, ASTBaseNode*> table;
};


// contains and owns all created nodes
// also contains symbol table
class ASTContext {
public:
  ASTContext()
    : errors(0) {
  }

  // frees all owned nodes
  ~ ASTContext();

  ASTExprNode* new_id_node(const std::string& id);
  ASTExprNode* new_dbl_node(const double val);
  ASTExprNode* new_llong_node(const long long val);

  ASTStrNode* new_str_node(const std::string str);
  ASTStrNode* new_str_node(const long long val_to_str);

  ASTListNode* new_list_node();

  ASTComponentNode* new_component_node(
      const std::string& name,
      ASTListNode* state_list,
      ASTStrNode* bond
  );

  ASTMoleculeNode* new_molecule_node(
      const std::string& name,
      ASTListNode* component_list
  );

  void inc_error_count() {
    errors++;
  }

  void print_error_report();

  // contains parameters from parameters sections
  // molecules from molecule types sections
  // and rules from reaction rules sections
  ASTSymbolTable symtab;

  void dump();

private:
  int errors;

  void remember_node(ASTBaseNode* n) {
    assert(all_nodes.count(n) == 0 && "Cannot add one node twice to be deleted");
    all_nodes.insert(n);
  }

  // set that contains all nodes, used to free up memory
  std::set<ASTBaseNode*> all_nodes;

  // input file names
  std::vector<std::string> file_names;
};

} // namespace BNG


#endif // LIBS_BNG_AST_H_
