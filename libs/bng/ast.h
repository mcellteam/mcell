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
class ASTListNode;


enum class NodeType {
  Invalid,
  Expr,
  List,
  Component,
  Mol,
  Compartment,
  Cplx,
  RxnRule,
  SeedSpecies,
  Observable,
  Str
};


enum class ExprType {
  Invalid,
  Id,
  Dbl,
  Llong,
  // operators
  UnaryPlus,
  UnaryMinus,

  Add,
  Sub,
  Mul,
  Div,
  Pow,

  FunctionCall
};


static const std::string BOND_STR_ANY = "?";
static const std::string BOND_STR_BOUND = "+";


class ASTBaseNode {
public:
  ASTBaseNode()
    : node_type(NodeType::Invalid), has_loc(false), line(0), file(nullptr) {
  }
  virtual ~ASTBaseNode() {}
  virtual void dump(const std::string ind) const;

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

  void set_loc(const ASTBaseNode* node) {
    assert(node->has_loc);
    line = node->line;
    file = node->file;
    has_loc = true;
  }

  bool is_expr() const {
    return node_type == NodeType::Expr;
  }

  bool is_str() const {
    return node_type == NodeType::Str;
  }

  bool is_list() const {
    return node_type == NodeType::List;
  }

  bool is_component() const {
    return node_type == NodeType::Component;
  }

  bool is_mol() const {
    return node_type == NodeType::Mol;
  }

  bool is_compartment() const {
    return node_type == NodeType::Compartment;
  }

  bool is_cplx() const {
    return node_type == NodeType::Cplx;
  }

  bool is_rxn_rule() const {
    return node_type == NodeType::RxnRule;
  }

  bool is_seed_species() const {
    return node_type == NodeType::SeedSpecies;
  }

  bool is_observable() const {
    return node_type == NodeType::Observable;
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
      left(nullptr), right(nullptr), args(nullptr), dbl(0), llong(0) {
    node_type = NodeType::Expr;
  }
  void dump(const std::string ind) const override;

  ExprType expr_type;

  void set_id(const std::string id_) {
    expr_type = ExprType::Id;
    id = id_;
  }

  void set_function_name(const std::string function_name_) {
    expr_type = ExprType::FunctionCall;
    function_name = function_name_;
  }

  void set_dbl(const double dbl_) {
    expr_type = ExprType::Dbl;
    dbl = dbl_;
  }

  void set_llong(const long long llong_) {
    expr_type = ExprType::Llong;
    llong = llong_;
  }

  void set_left(ASTExprNode* left_) {
    left = left_;
  }

  void set_right(ASTExprNode* right_) {
    right = right_;
  }

  void set_function_arguments(ASTListNode* args_) {
    args = args_;
  }

  void set_type(const ExprType op) {
    expr_type = op;
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

  bool is_unary_expr() const {
    return
        expr_type == ExprType::UnaryPlus ||
        expr_type == ExprType::UnaryMinus;
  }

  bool is_binary_expr() const {
    return
        expr_type == ExprType::Add ||
        expr_type == ExprType::Sub ||
        expr_type == ExprType::Mul ||
        expr_type == ExprType::Div ||
        expr_type == ExprType::Pow;
  }

  bool is_function_call() const {
    return expr_type == ExprType::FunctionCall;
  }

  const std::string& get_id() const {
    assert(is_id());
    return id;
  }

  const std::string& get_function_name() const {
    assert(is_function_call());
    return function_name;
  }

  double get_dbl() const {
    assert(is_dbl());
    return dbl;
  }

  double get_llong() const {
    assert(is_llong());
    return llong;
  }

  ASTExprNode* get_left() const {
    assert(!is_function_call());
    return left;
  }

  ASTExprNode* get_right() const {
    assert(!is_function_call());
    return right;
  }

  ExprType get_op() const {
    assert(is_unary_expr() || is_binary_expr());
    return expr_type;
  }

  ASTListNode* get_args() const {
    assert(is_function_call());
    return args;
  }

private:
  ASTExprNode* left;
  ASTExprNode* right;
  ASTListNode* args; // list of ASTExprNode*
  std::string id;
  std::string function_name;
  double dbl;
  long long llong;
};


// general node to store strings
class ASTStrNode: public ASTBaseNode {
public:
  ASTStrNode() {
    node_type = NodeType::Str;
  }
  void dump(const std::string ind) const override;

  std::string str;
};


class ASTListNode: public ASTBaseNode {
public:
  ASTListNode() {
    node_type = NodeType::List;
  }
  void dump(const std::string ind) const override;

  ASTListNode* append(ASTBaseNode* n) {
    assert(n != nullptr);
    items.push_back(n);
    return this;
  }

  size_t size() const {
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
  void dump(const std::string ind) const override;

  std::string name;
  ASTListNode* states;
  ASTStrNode* bond;
};


// specific instance of a molecule with states
// when defined in molecule types section, it is inserted
// into the symbol table
// when used in reaction rule, it is referenced from the
// reaction itself
// TODO: rename to ASTMolNode
class ASTMolNode: public ASTBaseNode {
public:
  ASTMolNode()
    : components(nullptr), compartment(nullptr) {
    node_type = NodeType::Mol;
  }
  void dump(const std::string ind) const override;

  bool has_compartment() const {
    return compartment != nullptr;
  }

  std::string name;
  ASTListNode* components;
  ASTStrNode* compartment; // nullptr when compartment is not set
};


// declaration of a compartment, it is referenced by its name in complexes
class ASTCompartmentNode: public ASTBaseNode {
public:
  ASTCompartmentNode()
    : dimensions(0), volume(nullptr)  {
    node_type = NodeType::Compartment;
  }
  void dump(const std::string ind) const override;

  std::string name;
  int dimensions; // expression that must evaluate to a constant 2 or 3
  ASTExprNode* volume;
  std::string parent_name; // parent name may be "" - this specifies a top level compartment
};


class ASTCplxNode: public ASTBaseNode {
public:
  ASTCplxNode(ASTMolNode* first_mol)
    : compartment(nullptr) {
    node_type = NodeType::Cplx;
    mols.push_back(first_mol);
  }
  void dump(const std::string ind) const override;

  ASTCplxNode* append(ASTMolNode* n) {
    assert(n != nullptr);
    mols.push_back(n);
    return this;
  }

  size_t size() const {
    return mols.size();
  }

  std::vector<ASTMolNode*> mols;

  // list node is also used for complexes and we need to put a complexes
  // compartment somewhere
  ASTStrNode* compartment;
};


class ASTRxnRuleNode: public ASTBaseNode {
public:
  ASTRxnRuleNode()
    : reactants(nullptr), reversible(false), products(nullptr), rates(nullptr) {
    node_type = NodeType::RxnRule;
  }
  void dump(const std::string ind) const override;

  std::string generate_name() const;

  std::string name; // set by user or generated automatically when empty
  ASTListNode* reactants;
  bool reversible;
  ASTListNode* products;
  ASTListNode* rates;
};


class ASTSeedSpeciesNode: public ASTBaseNode {
public:
  ASTSeedSpeciesNode()
    : cplx(nullptr), count(nullptr) {
    node_type = NodeType::SeedSpecies;
  }
  void dump(const std::string ind) const override;

  ASTCplxNode* cplx; // cplx to be released
  ASTExprNode* count;
};


class ASTObservableNode: public ASTBaseNode {
public:
  ASTObservableNode()
    : cplx_patterns(nullptr) {
    node_type = NodeType::Observable;
  }
  void dump(const std::string ind) const override;

  std::string type;
  std::string name;
  // list of complexes that are represented by ASTListNode as well
  ASTListNode* cplx_patterns;
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

  IdToNodeMap& get_as_map() {
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
    : single_cplx(nullptr), errors(0), eof_returned_as_newline(false),
      current_file(nullptr) {
  }

  // frees all owned nodes
  ~ ParserContext();

  // ------------- AST node manipulation -----------------

  ASTExprNode* new_id_node(const std::string& id, const BNGLLTYPE& loc);
  ASTExprNode* new_dbl_node(const double val, const BNGLLTYPE& loc);
  ASTExprNode* new_dbl_node(const double val, const ASTBaseNode* loc);
  ASTExprNode* new_dbl_node(const double val);
  ASTExprNode* new_llong_node(const long long val, const BNGLLTYPE& loc);

  ASTExprNode* new_expr_node(
      ASTExprNode* left,
      const ExprType op,
      ASTExprNode* right,
      const BNGLLTYPE& loc
  );

  ASTExprNode* new_expr_node(
      const std::string& function_name,
      ASTListNode* arguments,
      const BNGLLTYPE& loc
  );

  ASTStrNode* new_empty_str_node();
  ASTStrNode* new_str_node(const std::string str, const BNGLLTYPE& loc);
  ASTStrNode* new_str_node(const long long val_to_str, const BNGLLTYPE& loc);

  ASTListNode* new_list_node();

  ASTComponentNode* new_component_node(
      const std::string& name,
      ASTListNode* states,
      ASTStrNode* bond,
      const BNGLLTYPE& loc
  );

  ASTMolNode* new_molecule_node(
      const std::string& name,
      ASTListNode* components,
      ASTStrNode* compartment,
      const BNGLLTYPE& loc
  );

  ASTCompartmentNode* new_compartment_node(
      const std::string& name,
      const int dimensions,
      ASTExprNode* volume,
      const std::string parent_name, // may be ""
      const BNGLLTYPE& loc
  );

  // location is given by the first added molecule node
  ASTCplxNode* new_cplx_node(ASTMolNode* first_mol);

  ASTRxnRuleNode* new_rxn_rule_node(
      ASTStrNode* name,
      ASTListNode* reactants,
      const bool reversible,
      ASTListNode* products,
      ASTListNode* rates
  );

  ASTSeedSpeciesNode* new_seed_species_node(
      ASTCplxNode* cplx,
      ASTExprNode* count
  );

  ASTObservableNode* new_observable_node(
      const std::string& type,
      const std::string& name,
      ASTListNode* cplx_patterns,
      const BNGLLTYPE& loc // location may be used to report wrong observable type
  );

  void add_compartment(ASTCompartmentNode* n) {
    assert(n != nullptr);
    compartments.append(n);
  }

  void add_rxn_rule(ASTRxnRuleNode* n) {
    assert(n != nullptr);
    rxn_rules.append(n);
  }

  void add_seed_species(ASTSeedSpeciesNode* n) {
    assert(n != nullptr);
    seed_species.append(n);
  }

  void add_observable(ASTObservableNode* n) {
    assert(n != nullptr);
    observables.append(n);
  }

  // contains parameters from parameters sections
  // molecules from molecule types sections
  // and rules from reaction rules sections
  ASTSymbolTable symtab;

  ASTListNode compartments;

  ASTListNode rxn_rules;

  ASTListNode seed_species;

  ASTListNode observables;

  // single_cplx is set for mode when the parse parses a
  // single complex string with parse_single_cplx_string(),
  // no need to delete it because deletion is handled by vector 'all_nodes'
  ASTCplxNode* single_cplx;

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

  // handling of a missing newline at the end of the file
  void set_eof_returned_as_newline() {
    eof_returned_as_newline = true;
  }

  bool get_eof_returned_as_newline() {
    return eof_returned_as_newline;
  }

private:
  int errors;
  bool eof_returned_as_newline;

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

static inline const ASTListNode* to_list_node(const ASTBaseNode* n) {
  assert(n != nullptr);
  assert(n->is_list());
  return dynamic_cast<const ASTListNode*>(n);
}

static inline ASTMolNode* to_molecule_node(ASTBaseNode* n) {
  assert(n != nullptr);
  assert(n->is_mol());
  return dynamic_cast<ASTMolNode*>(n);
}

static inline const ASTComponentNode* to_component_node(const ASTBaseNode* n) {
  assert(n != nullptr);
  assert(n->is_component());
  return dynamic_cast<const ASTComponentNode*>(n);
}

static inline const ASTMolNode* to_mol_node(const ASTBaseNode* n) {
  assert(n != nullptr);
  assert(n->is_mol());
  return dynamic_cast<const ASTMolNode*>(n);
}

static inline const ASTCompartmentNode* to_compartment_node(const ASTBaseNode* n) {
  assert(n != nullptr);
  assert(n->is_compartment());
  return dynamic_cast<const ASTCompartmentNode*>(n);
}

static inline const ASTCplxNode* to_cplx_node(const ASTBaseNode* n) {
  assert(n != nullptr);
  assert(n->is_cplx());
  return dynamic_cast<const ASTCplxNode*>(n);
}

static inline ASTCplxNode* to_cplx_node(ASTBaseNode* n) {
  assert(n != nullptr);
  assert(n->is_cplx());
  return dynamic_cast<ASTCplxNode*>(n);
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

static inline const ASTSeedSpeciesNode* to_seed_species_node(const ASTBaseNode* n) {
  assert(n != nullptr);
  assert(n->is_seed_species());
  return dynamic_cast<const ASTSeedSpeciesNode*>(n);
}


static inline const ASTObservableNode* to_observable_node(const ASTBaseNode* n) {
  assert(n != nullptr);
  assert(n->is_observable());
  return dynamic_cast<const ASTObservableNode*>(n);
}


// returns bond index, BOND_VALUE_BOUND, BOND_VALUE_ANY or BOND_VALUE_UNBOUND
bond_value_t str_to_bond_value(const std::string& s);

} // namespace BNG

#endif // LIBS_BNG_AST_H_
