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

struct Node {
  Node()
    : is_mol(true), mol(nullptr), component(nullptr) {
  }

  Node(const MolInstance* mol_)
    : is_mol(true), mol(mol_), component(nullptr) {
  }

  Node(const ComponentInstance* component_)
    : is_mol(false), mol(nullptr), component(component_) {
  }

  bool compare(const Node& n2) const {
    const Node& n1 = *this;

    if (n1.is_mol != n2.is_mol) {
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
      // comparing !+
      if (n1.component->bond_value == BOND_VALUE_ANY) {
        // it is ok when the second node has !+ as well
        return n2.component->bond_value != BOND_VALUE_NO_BOND;
      }
      if (n2.component->bond_value == BOND_VALUE_ANY) {
        return n1.component->bond_value != BOND_VALUE_NO_BOND;
      }

      // we do not care about actual bond values because to what is the component connected is
      // represented by the graph
      return true;
    }
  }

  bool is_mol;
  const MolInstance* mol;
  const ComponentInstance* component;
};

static std::ostream & operator<<(std::ostream &out, const Node& n);

// molecule instance/pattern is a node in the complex graph, edge is a bond
typedef boost::property<boost::vertex_name_t, Node, boost::property< boost::vertex_index_t, int > > MtVertexProperty;

// using a vecS graphs => the index maps are implicit.
// we are not using any edge property
typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::bidirectionalS, MtVertexProperty> Graph;

typedef boost::property_map<Graph, boost::vertex_name_t >::type VertexNameMap;

typedef std::vector<std::pair<Graph::vertex_descriptor, Graph::vertex_descriptor>> MappingVector;
typedef std::vector<MappingVector> MultipleMappingsVector;

// finds all subgraph isomorphism mappings of pattern graph on cplx graph
void get_subgraph_isomorphism_mappings(Graph& pattern, Graph& cplx, const bool only_first_match, MultipleMappingsVector& res);

} // namespace BNG

#endif /* LIBS_BNG_GRAPH_H_ */
