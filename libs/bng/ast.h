#ifndef LIBS_BNG_AST_H_
#define LIBS_BNG_AST_H_

#include <string>

#include "bng_defines.h"

namespace BNG {

enum class NodeType {
  Invalid,
  Expr,
  List
};

enum class ExprType {
  Invalid,
  Id,
  Dbl,
  Llong
};


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

  void dump(const std::string ind) override;

  // getters that check the type
private:
  ASTExprNode* left;
  ASTExprNode* right;
  std::string id;
  double dbl;
  long long llong;
};


class ASTListNode: public ASTBaseNode {
public:
  ASTListNode() {}
  void dump(const std::string ind) override;

  std::vector<ASTBaseNode*> items;
};


class ASTSymbolTable {
public:
  void insert(const std::string id, ASTBaseNode* node);
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
  ASTContext() {
  }

  // frees all owned nodes
  ~ ASTContext();

  ASTExprNode* new_id_node(const std::string id);
  ASTExprNode* new_dbl_node(const double val);
  ASTExprNode* new_llong_node(const long long val);

  ASTSymbolTable symtab;

  void dump();

private:
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
