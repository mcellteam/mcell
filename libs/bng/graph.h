/*
 * graph.h
 *
 *  Created on: Jun 22, 2020
 *      Author: ahusar
 */

#ifndef LIBS_BNG_GRAPH_H_
#define LIBS_BNG_GRAPH_H_

#include <iostream>

#include "bng/mol_type.h"
#include "bng/mol_instance.h"

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/iteration_macros.hpp>

namespace BNG {

class BNGData;

struct Node {
  Node()
    : is_mol(true), mol(nullptr), component(nullptr), used_in_rxn_product(true),
      product_index(INDEX_INVALID), ordering_index(INDEX_INVALID), modified_ordering_index(ordering_index) {
  }

  Node(MolInstance* mol_)
    : is_mol(true), mol(mol_), component(nullptr), used_in_rxn_product(true),
      product_index(INDEX_INVALID), ordering_index(INDEX_INVALID), modified_ordering_index(INDEX_INVALID) {
  }

  Node(ComponentInstance* component_)
    : is_mol(false), mol(nullptr), component(component_), used_in_rxn_product(true),
      product_index(INDEX_INVALID), ordering_index(INDEX_INVALID), modified_ordering_index(INDEX_INVALID) {
  }

  bool compare(const Node& n2) const {
    const Node& n1 = *this;

    if (n1.modified_ordering_index != n2.modified_ordering_index) {
      // extra ordering needed to get unique products,
      // modified_ordering_index is set to INDEX_INVALID by default and
      // set toa different value only when needed
      return false;
    }
    else if (n1.is_mol != n2.is_mol) {
      return false;
    }
    else if (n1.is_mol) {
      assert(n2.is_mol);
      assert(n1.mol != nullptr);
      assert(n2.mol != nullptr);

      // molecule
      return n1.mol->mol_type_id == n2.mol->mol_type_id;
    }
    else {
      assert(!n1.is_mol && !n2.is_mol);
      assert(n1.component != nullptr);
      assert(n2.component != nullptr);

      // component
      if (n1.component->component_type_id != n2.component->component_type_id) {
        // must be the same
        return false;
      }

      // state
      if (n1.component->state_id != STATE_ID_DONT_CARE && n2.component->state_id != STATE_ID_DONT_CARE) {
        // must be the same or don't care for one of the compared nodes
        if (n1.component->state_id != n2.component->state_id) {
          return false;
        }
      }

      // bond
      // comparing !?
      if (n1.component->bond_value == BOND_VALUE_ANY ||
        n2.component->bond_value == BOND_VALUE_ANY) {
        return true;
      }

      // comparing !+
      if (n1.component->bond_value == BOND_VALUE_BOUND) {
        // it is ok when the second node has !+ as well
        return n2.component->bond_value != BOND_VALUE_UNBOUND;
      }
      else if (n1.component->bond_value == BOND_VALUE_UNBOUND) {
        // no bond means that there must be no bond on the other side either
        return n2.component->bond_value == BOND_VALUE_UNBOUND;
      }
      if (n2.component->bond_value == BOND_VALUE_BOUND) {
        return n1.component->bond_value != BOND_VALUE_UNBOUND;
      }
      else if (n2.component->bond_value == BOND_VALUE_UNBOUND) {
        return n1.component->bond_value == BOND_VALUE_UNBOUND;
      }

      // we do not care about actual bond values because to what is the component connected is
      // represented by the graph
      return true;
    }
  }

  std::string to_str(const BNGData* bng_data = nullptr) const;

  bool is_mol;
  MolInstance* mol;
  ComponentInstance* component;

  // for reaction handling, default is true
  bool used_in_rxn_product;
  // also for reactions - specifies index of the reaction product
  // needed to match patterns onto resulting products
  uint product_index;

  // auxiliary ordering markers, used when determining whether two product sets are identical
  // ordering_index is set to all components and when a reactant changes by 
  // a rule application, its value is copied to the modified_ordering_index
  uint ordering_index;
  uint modified_ordering_index;
};

std::ostream & operator<<(std::ostream &out, const Node& n);

// molecule instance/pattern is a node in the complex graph, edge is a bond
typedef boost::property<boost::vertex_name_t, Node, boost::property< boost::vertex_index_t, int > > MtVertexProperty;

// NOTE: we can experiment which underlying types are more efficient (vecS/listS/setS...)
// we are not using any edge property
typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, MtVertexProperty> Graph;

typedef boost::property_map<Graph, boost::vertex_name_t >::type VertexNameMap;

// FIXME: this has to be a map, not a vector of pairs..
typedef std::map<Graph::vertex_descriptor, Graph::vertex_descriptor> VertexMapping;
typedef std::vector<VertexMapping> VertexMappingVector;

typedef Graph::vertex_descriptor vertex_descriptor_t;
typedef Graph::edge_descriptor edge_descriptor_t;

typedef boost::graph_traits<Graph>::vertex_iterator vertex_iter_t; // TODO: use this typedef wherever needed

// finds all subgraph isomorphism mappings of pattern graph on cplx graph
void get_subgraph_isomorphism_mappings(
    Graph& pattern,
    Graph& cplx,
    const bool only_first_match,
    VertexMappingVector& res
);

void dump_graph(const Graph& g_const, const BNGData* bng_data = nullptr, const std::string ind = "");
void dump_graph_mapping(const VertexMapping& mapping);

} // namespace BNG

#endif /* LIBS_BNG_GRAPH_H_ */
