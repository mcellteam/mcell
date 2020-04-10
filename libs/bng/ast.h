#ifndef LIBS_BNG_AST_H_
#define LIBS_BNG_AST_H_

/**
 * This file defines abstract syntax tree for BNGL.
 * Compilation must be done in two steps because
 * BNGL does not require definition before use
 * (at least in expressions).
 */

#include <string>

#include "bng/bng_defines.h"

// The declaration below should come from #include "bngl_parser.hpp"
// but this include causes many conflicts
#if ! defined BNGLLTYPE && ! defined BNGLLTYPE_IS_DECLARED
struct BNGLLTYPE
{
  int first_line;
  int first_column;
  int last_line;
  int last_column;
};
# define BNGLLTYPE_IS_DECLARED 1
#endif

namespace BNG {

class ParserContext;
class ASTExprNode;
class ASTRxnRuleNode;


enum class NodeType {
  Invalid,
  Expr,
  List,
  Molecule,
  Component,
  RxnRule,
  Str,
  Separator
};


enum class ExprType {
  Invalid,
  Id,
  Dbl,
  Llong
};


enum class SeparatorType {
  Invalid,
  Plus,
  Dot
};


static const std::string BOND_STR_ANY = "+";


class ASTBaseNode {
public:
  ASTBaseNode()
    : node_type(NodeType::Invalid), has_loc(false), line(0), file(nullptr) {
  }
  virtual ~ASTBaseNode() {}
  virtual void dump(const std::string ind);

  void set_loc(const char* file_, const BNGLLTYPE& loc) {
    file = file_;
    line = loc.first_line;
    has_loc = true;
  }

  void set_loc(const char* file_, const int line_) {
    file = file_;
    line = line_;
    has_loc = true;
  }

  bool is_expr() const {
    return node_type == NodeType::Expr;
  }

  bool is_str() const {
    return node_type == NodeType::Str;
  }

  bool is_molecule() const {
    return node_type == NodeType::Molecule;
  }

  bool is_component() const {
    return node_type == NodeType::Component;
  }

  bool is_rxn_rule() const {
    return node_type == NodeType::RxnRule;
  }

  bool is_separator() const {
    return node_type == NodeType::Separator;
  }

  NodeType node_type;
  bool has_loc;
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

  bool is_id() const {
    return expr_type == ExprType::Id;
  }

  bool is_dbl() const {
    return expr_type == ExprType::Dbl;
  }

  bool is_llong() const {
    return expr_type == ExprType::Llong;
  }

  const std::string& get_id() const {
    assert(is_id());
    return id;
  }

  double get_dbl() const {
    assert(is_dbl());
    return dbl;
  }

  double get_llong() const {
    assert(is_llong());
    return llong;
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


// separators represent:
//  '+' - adding a reactant or a product to a side of a reaction)
//  '.' - creating a complex with another molecule
class ASTSeparatorNode: public ASTBaseNode {
public:
  ASTSeparatorNode()
    : separator_type(SeparatorType::Invalid) {
    node_type = NodeType::Separator;
  }

  bool is_dot() const {
    return separator_type == SeparatorType::Dot;
  }

  bool is_plus() const {
    return separator_type == SeparatorType::Plus;
  }

  void dump(const std::string ind) override;

  char to_char() const {
    assert(separator_type != SeparatorType::Invalid);
    return (separator_type == SeparatorType::Dot) ? '.' : '+';
  }

  SeparatorType separator_type;
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

  size_t size() {
    return items.size();
  }

  std::vector<ASTBaseNode*> items;
};


class ASTComponentNode: public ASTBaseNode {
public:
  ASTComponentNode()
    : states(nullptr), bond(nullptr) {
    node_type = NodeType::Component;
  }
  void dump(const std::string ind) override;

  std::string name;
  ASTListNode* states;
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
    : components(nullptr) {
    node_type = NodeType::Molecule;
  }
  void dump(const std::string ind) override;

  std::string name;
  ASTListNode* components;
};


class ASTRxnRuleNode: public ASTBaseNode {
public:
  ASTRxnRuleNode()
    : reactants(nullptr), reversible(false), products(nullptr), rates(nullptr) {
    node_type = NodeType::RxnRule;
  }
  void dump(const std::string ind) override;

  std::string generate_name() const;

  std::string name; // generated automatically for now
  ASTListNode* reactants;
  bool reversible;
  ASTListNode* products;
  ASTListNode* rates;
};


class ASTSymbolTable {
public:
  typedef std::map<std::string, ASTBaseNode*> IdToNodeMap;

  void insert(const std::string id, ASTBaseNode* node, ParserContext* ctx);
  void insert_molecule_declarations(const ASTListNode* molecule_node_list, ParserContext* ctx);

  // if symbol does was not defined, returns null and prints out error message
  ASTBaseNode* get(const std::string& id, ASTBaseNode* loc, ParserContext* ctx) const;

  const IdToNodeMap& get_as_map() const {
    return table;
  }

  void dump();
private:
  // we do not hold constant pointers because one for instance might
  // want to evaluate an expression
  IdToNodeMap table;
};


// contains and owns all created nodes
// also contains symbol table
class ParserContext {
public:
  ParserContext()
    : errors(0), current_file(nullptr) {
  }

  // frees all owned nodes
  ~ ParserContext();

  // ------------- AST node manipulation -----------------

  ASTExprNode* new_id_node(const std::string& id, const BNGLLTYPE& loc);
  ASTExprNode* new_dbl_node(const double val, const BNGLLTYPE& loc);
  ASTExprNode* new_dbl_node(const double val, const ASTBaseNode* loc);
  ASTExprNode* new_llong_node(const long long val, const BNGLLTYPE& loc);

  ASTStrNode* new_empty_str_node();
  ASTStrNode* new_str_node(const std::string str, const BNGLLTYPE& loc);
  ASTStrNode* new_str_node(const long long val_to_str, const BNGLLTYPE& loc);

  ASTListNode* new_list_node();

  ASTSeparatorNode* new_separator_node(const SeparatorType type, const BNGLLTYPE& loc);

  ASTComponentNode* new_component_node(
      const std::string& name,
      ASTListNode* states,
      ASTStrNode* bond,
      const BNGLLTYPE& loc
  );

  ASTMoleculeNode* new_molecule_node(
      const std::string& name,
      ASTListNode* components,
      const BNGLLTYPE& loc
  );

  ASTRxnRuleNode* new_rxn_rule_node(
      ASTListNode* reactants,
      const bool reversible,
      ASTListNode* products,
      ASTListNode* rates
  );

  void add_rxn_rule(ASTRxnRuleNode* n) {
    assert(n != nullptr);
    rxn_rules.append(n);
  }

  // contains parameters from parameters sections
  // molecules from molecule types sections
  // and rules from reaction rules sections
  ASTSymbolTable symtab;

  ASTListNode rxn_rules;

  // ------------- other parsing utilities -----------------

  void inc_error_count() {
    errors++;
  }

  int get_error_count() const {
    return errors;
  }

  void internal_error(const ASTBaseNode* loc, const std::string msg);
  void print_error_report();

  void set_current_file_name(const char* file_name) {
    // insert the file name into a set that owns all the name strings
    auto pair_it_bool = file_names.insert(file_name);
    current_file = pair_it_bool.first->c_str();
  }

  const char* get_current_file_name() const {
    return current_file;
  }

  void dump();

  // context owns all strings parsed as IDs
  const char* insert_to_string_pool(const char* s) {
    auto it = parser_string_pool.insert(s);
    return it.first->c_str();
  }

private:
  int errors;

  void remember_node(ASTBaseNode* n) {
    assert(all_nodes.count(n) == 0 && "Cannot add one node twice to be deleted");
    all_nodes.insert(n);
  }

  // set that contains all nodes, used to free up memory
  std::set<ASTBaseNode*> all_nodes;

  // all input file names are owned by this set and only pointers are passed around
  const char* current_file;
  std::set<std::string> file_names;


  // container for all parsed strings,
  // with GLR parser we must not delete strings returned as the attribute of a token,
  // so we must maintain them somewhere, should not have much impact on performance
  std::set<std::string> parser_string_pool;
};


// ------------- dynamic cast utilities -----------------

// these cannot be methods of ASTBaseNode because the target data types are undefined
// when ASTBaseNode is being defined

static inline ASTExprNode* to_expr_node(ASTBaseNode* n) {
  assert(n != nullptr);
  assert(n->is_expr());
  return dynamic_cast<ASTExprNode*>(n);
}

static inline const ASTStrNode* to_str_node(const ASTBaseNode* n) {
  assert(n != nullptr);
  assert(n->is_str());
  return dynamic_cast<const ASTStrNode*>(n);
}

static inline ASTMoleculeNode* to_molecule_node(ASTBaseNode* n) {
  assert(n != nullptr);
  assert(n->is_molecule());
  return dynamic_cast<ASTMoleculeNode*>(n);
}

static inline const ASTMoleculeNode* to_molecule_node(const ASTBaseNode* n) {
  assert(n != nullptr);
  assert(n->is_molecule());
  return dynamic_cast<const ASTMoleculeNode*>(n);
}

static inline const ASTComponentNode* to_component_node(const ASTBaseNode* n) {
  assert(n != nullptr);
  assert(n->is_component());
  return dynamic_cast<const ASTComponentNode*>(n);
}

static inline ASTRxnRuleNode* to_rxn_rule_node(ASTBaseNode* n) {
  assert(n != nullptr);
  assert(n->is_rxn_rule());
  return dynamic_cast<ASTRxnRuleNode*>(n);
}

static inline const ASTRxnRuleNode* to_rxn_rule_node(const ASTBaseNode* n) {
  assert(n != nullptr);
  assert(n->is_rxn_rule());
  return dynamic_cast<const ASTRxnRuleNode*>(n);
}

static inline const ASTSeparatorNode* to_separator(const ASTBaseNode* n) {
  assert(n != nullptr);
  assert(n->is_separator());
  return dynamic_cast<const ASTSeparatorNode*>(n);
}

// returns bond index, BOND_VALUE_NO_BOND, BOND_VALUE_ANY or BOND_VALUE_INVALID
bond_value_t str_to_bond_value(const std::string& s);

} // namespace BNG

#endif // LIBS_BNG_AST_H_
