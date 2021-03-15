/*
 * graph.h
 *
 *  Created on: Jun 22, 2020
 *      Author: ahusar
 */

#ifndef LIBS_BNG_GRAPH_H_
#define LIBS_BNG_GRAPH_H_

#include <iostream>

#include "bng/elem_mol_type.h"
#include "bng/elem_mol.h"

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/iteration_macros.hpp>
#define BOOST_POOL_NO_MT
#include <boost/pool/pool_alloc.hpp>
#undef BOOST_POOL_NO_MT

namespace boost {
  // NOTE: not sure if we need to put this really into the boost namespace
  struct pool_listS { };

  // examples of allocators usage are here: https://theboostcpplibraries.com/boost.pool
  // TODO: figure out whether we need to use purge_memory
  //
  // note from https://www.boost.org/doc/libs/1_65_0/libs/pool/doc/html/boost/fast_pool_allocator.html:
  // The underlying singleton_pool used by the this allocator constructs a pool instance
  // that is never freed. This means that memory allocated by the allocator can be still
  // used after main() has completed, but may mean that some memory checking programs will
  // complain about leaks.

  template <class ValueType>
  struct container_gen<pool_listS, ValueType> {
    typedef std::list<ValueType, boost::fast_pool_allocator<ValueType>> type;
  };

  template <>
  struct parallel_edge_traits<pool_listS> {
    typedef allow_parallel_edge_tag type; 
  };

}


namespace BNG {

class BNGData;

struct Node {
  Node()
    : is_mol(true), mol(nullptr), component(nullptr), used_in_rxn_product(true),
      product_index(INDEX_INVALID), reactant_pattern_index(INDEX_INVALID),
      ordering_index(INDEX_INVALID), modified_ordering_index(ordering_index) {
  }

  Node(ElemMol* mol_)
    : is_mol(true), mol(mol_), component(nullptr), used_in_rxn_product(true),
      product_index(INDEX_INVALID), reactant_pattern_index(INDEX_INVALID),
      ordering_index(INDEX_INVALID), modified_ordering_index(INDEX_INVALID) {
  }

  Node(Component* component_)
    : is_mol(false), mol(nullptr), component(component_), used_in_rxn_product(true),
      product_index(INDEX_INVALID), reactant_pattern_index(INDEX_INVALID),
      ordering_index(INDEX_INVALID), modified_ordering_index(INDEX_INVALID) {
  }

  // returns true if n2 matches pattern represented by 'this'
  bool compare(const Node& n2) const {
    const Node& n1 = *this;

    if (n1.reactant_pattern_index != INDEX_INVALID &&
        n2.reactant_pattern_index != INDEX_INVALID &&
        n1.reactant_pattern_index != n2.reactant_pattern_index) {
      // in bimol rxns, we must check that the pattern matched the correct reactant
      // checked only when both are set
      return false;
    }
    else if (n1.modified_ordering_index != n2.modified_ordering_index) {
      // extra ordering needed to get unique products,
      // modified_ordering_index is set to INDEX_INVALID by default and
      // set to a different value only when needed
      return false;
    }
    else if (n1.is_mol != n2.is_mol) {
      return false;
    }
    else if (n1.is_mol) {
      // molecule
      assert(n2.is_mol);
      assert(n1.mol != nullptr);
      assert(n2.mol != nullptr);

      // molecule type
      if (n1.mol->elem_mol_type_id != n2.mol->elem_mol_type_id) {
        return false;
      }

      // does compartment match?
      if (n1.mol->compartment_id != COMPARTMENT_ID_NONE &&
          !is_in_out_compartment_id(n1.mol->compartment_id) &&
          n1.mol->compartment_id != n2.mol->compartment_id) {
        return false;
      }

      return true;
    }
    else {
      // component
      assert(!n1.is_mol && !n2.is_mol);
      assert(n1.component != nullptr);
      assert(n2.component != nullptr);

      // component type
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
  ElemMol* mol;
  Component* component;

  // for reaction handling, default is true
  bool used_in_rxn_product;
  // also for reactions - specifies index of the reaction product
  // needed to match patterns onto resulting products
  uint product_index;

  // used when creating product sets - in bimol rxns one reactant pattern must match
  // exactly one reactant, not both
  uint reactant_pattern_index;

  // auxiliary ordering markers, used when determining whether two product sets are identical
  // ordering_index is set to all components and when a reactant changes by 
  // a rule application, its value is copied to the modified_ordering_index
  uint ordering_index;
  uint modified_ordering_index;
};

std::ostream & operator<<(std::ostream &out, const Node& n);

// molecule instance/pattern is a node in the complex graph, edge is a bond
typedef boost::property<boost::vertex_name_t, Node, boost::property< boost::vertex_index_t, int > > MtVertexProperty;

// one can choose different underlying types (vecS/listS/setS...) but
// vecS seems to be the most efficient
// last 7th unlisted argument EdgeList must be listS (default) otherwise it is not possible to remove edges
typedef boost::adjacency_list<
    boost::vecS, boost::vecS, boost::undirectedS, MtVertexProperty,
    boost::no_property, boost::no_property,
    boost::pool_listS
    > Graph;

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
